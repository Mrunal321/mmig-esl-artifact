/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file mmig_optimizer.hpp
  \brief mMIG optimization pipeline wrapper
*/

#pragma once

#include "cleanup.hpp"
#include "mmig_algebraic_rewriting.hpp"
#include "mmig_balancing.hpp"
#include "mmig_cec_guard.hpp"
#include "mmig_cut_rewriting.hpp"
#include "mmig_exact_rewriting.hpp"
#include "mmig_cone_polarity_flip.hpp"
#include "mmig_inv_optimization.hpp"
#include "mmig_inv_propagation.hpp"
#include "mmig_minority_seeding.hpp"
#include "mmig_refactoring.hpp"
#include "mmig_resubstitution.hpp"

#include "../traits.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/depth_view.hpp"

#include <algorithm>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace mockturtle
{

enum class mmig_stage_mode : uint8_t
{
  pre,
  post,
  both
};

enum class mmig_advanced_flow : uint8_t
{
  round_robin,
  dac19_default,
  dac19_area,
  compress2rs,
  paper2019,
  epfl,
  legacy
};

struct mmig_optimizer_params
{
  uint32_t max_iterations{ 16u };
  bool enable_minority_seeding{ true };
  mmig_seed_mode seeding_mode{ mmig_seed_mode::both };
  uint32_t seeding_max_candidates{ 24u };
  uint32_t seeding_max_rounds{ 1u };
  bool seeding_incremental_reranking{ false };
  bool enable_sr5{ false };
  bool auto_sr5_with_cec{ true };  /* P7: auto-enable SR5 when CEC is active */
  mmig_stage_mode stage_mode{ mmig_stage_mode::both };
  bool run_cec{ true };
  bool enable_advanced{ false };
  uint32_t advanced_rounds{ 1u };
  mmig_advanced_flow advanced_flow{ mmig_advanced_flow::round_robin };
  bool advanced_allow_zero_gain{ false };
  bool advanced_run_resub{ true };
  bool advanced_run_cut{ true };
  bool advanced_run_refactor{ true };
  bool advanced_run_exact{ true };
  bool advanced_run_balance{ true };
  bool advanced_tuned_policy{ true };
  bool advanced_minority_aware_cost{ true };
  bool advanced_reward_minority_count{ false };
  uint32_t advanced_resub_candidate_cec_budget{ 0u };
  uint32_t advanced_exact_gate_limit{ 0u };
  uint32_t advanced_exact_min_depth{ 0u };
  uint32_t advanced_exact_min_inverted_edges{ 0u };
  bool advanced_exact_require_minority{ false };
  bool advanced_rollback_on_objective_regression{ true };
  bool advanced_stop_on_stagnation{ true };
  uint32_t advanced_max_stagnation_rounds{ 1u };
  bool enable_dont_cares{ false };           /* P2: don't-care-aware resub/refactor/exact */
  uint32_t dont_care_gate_threshold{ 10000u }; /* P2: only for circuits below this size */
  bool enable_interleaved_seeding{ false };  /* P3: re-seed minorities between advanced rounds */
  uint32_t interleaved_seeding_interval{ 2u }; /* P3: re-seed every N rounds */
  bool enable_dual_inv_propagation{ true };  /* P1: dual-inversion propagation */
  bool rewrite_rank_candidates{ true };
  bool rewrite_normalize_complemented_inner{ true };
  bool resub_rank_candidates{ true };
  uint32_t resub_candidate_pool_size{ 24u };
  uint32_t resub_full_score_limit{ 6u };
  bool verbose{ false };
};

struct mmig_optimizer_stats
{
  stopwatch<>::duration time_total{ 0 };
  mmig_cec_guard_stats cec{};
  mmig_minority_seeding_stats seeding_pre{};
  mmig_algebraic_rewriting_stats rewrite_pre{};
  mmig_algebraic_rewriting_stats rewrite_post{};
  mmig_inv_propagation_stats inv_propagation_pre{};
  mmig_inv_propagation_stats inv_propagation_post{};
  mmig_inv_optimization_stats inv_optimization_pre{};
  mmig_inv_optimization_stats inv_optimization_post{};
  mmig_algebraic_rewriting_stats advanced_mighty{};
  mmig_resubstitution_stats advanced_resub{};
  mmig_cut_rewriting_stats advanced_cut{};
  mmig_refactoring_stats advanced_refactor{};
  mmig_exact_rewriting_stats advanced_exact{};
  mmig_balancing_stats advanced_balance{};
  uint32_t advanced_rounds_executed{ 0u };
  uint32_t advanced_rounds_improved{ 0u };
  uint32_t advanced_rounds_rolled_back{ 0u };
  uint32_t advanced_stagnation_stops{ 0u };
  uint32_t advanced_depth_before{ 0u };
  uint32_t advanced_depth_after{ 0u };
  uint32_t advanced_gates_before{ 0u };
  uint32_t advanced_gates_after{ 0u };
  uint32_t advanced_majority_before{ 0u };
  uint32_t advanced_majority_after{ 0u };
  uint32_t advanced_minority_before{ 0u };
  uint32_t advanced_minority_after{ 0u };
  uint32_t advanced_inverted_edges_before{ 0u };
  uint32_t advanced_inverted_edges_after{ 0u };
  uint32_t advanced_effective_inverted_edges_before{ 0u };
  uint32_t advanced_effective_inverted_edges_after{ 0u };
  uint64_t advanced_effective_area_score_before{ 0u };
  uint64_t advanced_effective_area_score_after{ 0u };
};

namespace detail
{

template<class Ntk, class Fn>
void run_mmig_stage( Ntk& ntk,
                     std::string const& stage_name,
                     mmig_optimizer_params const& ps,
                     mmig_cec_guard_stats& cec_stats,
                     Fn&& fn )
{
  mmig_cec_guard_params guard_ps;
  guard_ps.run_cec = ps.run_cec;
  guard_ps.verbose = ps.verbose;
  run_with_mmig_cec_guard( ntk, stage_name, std::forward<Fn>( fn ), guard_ps, &cec_stats );
}

inline void accumulate_resubstitution_stats( resubstitution_stats& dst, resubstitution_stats const& src )
{
  dst.time_total += src.time_total;
  dst.time_divs += src.time_divs;
  dst.time_resub += src.time_resub;
  dst.time_callback += src.time_callback;
  dst.num_total_divisors += src.num_total_divisors;
  dst.estimated_gain += src.estimated_gain;
  if ( dst.initial_size == 0u )
  {
    dst.initial_size = src.initial_size;
  }
}

inline void accumulate_mmig_resubstitution_stats( mmig_resubstitution_stats& dst, mmig_resubstitution_stats const& src )
{
  dst.time_total += src.time_total;
  accumulate_resubstitution_stats( dst.resub, src.resub );
  dst.num_candidate_checks += src.num_candidate_checks;
  dst.num_candidate_accepted += src.num_candidate_accepted;
  dst.num_candidate_rejected += src.num_candidate_rejected;
  dst.num_candidate_inconclusive += src.num_candidate_inconclusive;
  dst.num_candidate_skipped_non_mixed += src.num_candidate_skipped_non_mixed;
  dst.num_candidate_skipped_budget += src.num_candidate_skipped_budget;
  dst.time_candidate_cec += src.time_candidate_cec;
}

inline void accumulate_cut_rewriting_stats( cut_rewriting_stats& dst, cut_rewriting_stats const& src )
{
  dst.time_total += src.time_total;
  dst.time_cuts += src.time_cuts;
  dst.time_rewriting += src.time_rewriting;
  dst.time_mis += src.time_mis;
}

inline void accumulate_rewrite_stats( rewrite_stats& dst, rewrite_stats const& src )
{
  dst.time_total += src.time_total;
  dst.estimated_gain += src.estimated_gain;
  dst.candidates += src.candidates;
}

inline void accumulate_mmig_cut_rewriting_stats( mmig_cut_rewriting_stats& dst, mmig_cut_rewriting_stats const& src )
{
  dst.time_total += src.time_total;
  accumulate_cut_rewriting_stats( dst.cut, src.cut );
  dst.num_depth_rollbacks += src.num_depth_rollbacks;
  dst.depth_before = src.depth_before;
  dst.depth_after = src.depth_after;
}

inline void accumulate_balancing_stats( balancing_stats& dst, balancing_stats const& src )
{
  dst.time_total += src.time_total;
  dst.cut_enumeration_st.time_total += src.cut_enumeration_st.time_total;
  dst.cut_enumeration_st.time_truth_table += src.cut_enumeration_st.time_truth_table;
}

inline void accumulate_mmig_balancing_stats( mmig_balancing_stats& dst, mmig_balancing_stats const& src )
{
  dst.time_total += src.time_total;
  accumulate_balancing_stats( dst.balancing, src.balancing );
  dst.num_depth_rollbacks += src.num_depth_rollbacks;
  dst.depth_before = src.depth_before;
  dst.depth_after = src.depth_after;
}

inline void accumulate_refactoring_stats( refactoring_stats& dst, refactoring_stats const& src )
{
  dst.time_total += src.time_total;
  dst.time_mffc += src.time_mffc;
  dst.time_refactoring += src.time_refactoring;
  dst.time_simulation += src.time_simulation;
}

inline void accumulate_mmig_refactoring_stats( mmig_refactoring_stats& dst, mmig_refactoring_stats const& src )
{
  dst.time_total += src.time_total;
  accumulate_refactoring_stats( dst.refac, src.refac );
  dst.num_depth_rollbacks += src.num_depth_rollbacks;
  dst.depth_before = src.depth_before;
  dst.depth_after = src.depth_after;
}

inline void accumulate_mmig_exact_rewriting_stats( mmig_exact_rewriting_stats& dst, mmig_exact_rewriting_stats const& src )
{
  dst.time_total += src.time_total;
  accumulate_rewrite_stats( dst.rewrite, src.rewrite );
  dst.num_depth_rollbacks += src.num_depth_rollbacks;
  dst.num_gate_skips += src.num_gate_skips;
  dst.depth_before = src.depth_before;
  dst.depth_after = src.depth_after;
  dst.observed_gates = src.observed_gates;
  dst.observed_minority = src.observed_minority;
  dst.observed_inverted_edges = src.observed_inverted_edges;
}

inline void accumulate_mmig_algebraic_rewriting_stats( mmig_algebraic_rewriting_stats& dst, mmig_algebraic_rewriting_stats const& src )
{
  dst.time_total += src.time_total;
  dst.num_iterations += src.num_iterations;
  dst.num_candidates += src.num_candidates;
  dst.num_evaluated += src.num_evaluated;
  dst.num_applied += src.num_applied;
  dst.num_rejected += src.num_rejected;
  dst.num_depth_improvements += src.num_depth_improvements;
  dst.num_area_improvements += src.num_area_improvements;
  dst.num_inverter_improvements += src.num_inverter_improvements;
  for ( uint32_t i = 0u; i < dst.per_rule_applied.size(); ++i )
  {
    dst.per_rule_applied[i] += src.per_rule_applied[i];
    dst.per_rule_rejected[i] += src.per_rule_rejected[i];
  }
}

enum class mmig_advanced_stage : uint8_t
{
  mighty_area,
  resub,
  cut,
  cut_zero,
  refactor,
  refactor_zero,
  exact_rewriting,
  balancing
};

struct mmig_advanced_step
{
  mmig_advanced_stage stage;
  std::optional<uint32_t> rs_cut_override{};
  std::optional<uint32_t> rs_insert_override{};
};

template<class Ntk>
struct mmig_objective_cost
{
  uint32_t depth{ 0u };
  uint32_t gates{ 0u };
  uint32_t majority_gates{ 0u };
  uint32_t minority_gates{ 0u };
  uint32_t inverted_edges{ 0u };
  uint32_t effective_inverted_edges{ 0u };
  uint64_t effective_area_score{ 0u };
};

template<class Ntk>
mmig_objective_cost<Ntk> compute_mmig_objective_cost( Ntk const& ntk )
{
  mmig_objective_cost<Ntk> cost{};
  depth_view<Ntk> dv{ ntk };
  cost.depth = dv.depth();
  cost.gates = ntk.num_gates();
  ntk.foreach_gate( [&]( auto const& n ) {
    if ( ntk.is_min( n ) )
    {
      ++cost.minority_gates;
    }
    else
    {
      ++cost.majority_gates;
    }
  } );

  uint64_t inverted = 0u;
  ntk.foreach_gate( [&]( auto const& n ) {
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( ntk.is_complemented( f ) )
      {
        ++inverted;
      }
    } );
  } );
  ntk.foreach_po( [&]( auto const& s, auto ) {
    if ( ntk.is_complemented( s ) )
    {
      ++inverted;
    }
  } );
  cost.inverted_edges = static_cast<uint32_t>( std::min<uint64_t>( inverted, std::numeric_limits<uint32_t>::max() ) );
  auto const minority_credit = std::min( cost.minority_gates, cost.inverted_edges );
  cost.effective_inverted_edges = cost.inverted_edges - minority_credit;

  /* P5: Adaptive cost function weight based on inversion density.
   * Inversion-heavy circuits (crypto) benefit from heavier inverter penalty,
   * while gate-dominated circuits (arithmetic) use the standard weight. */
  uint64_t inv_weight = 1u;
  if ( cost.gates > 0u )
  {
    auto const inv_ratio = static_cast<double>( cost.effective_inverted_edges ) / cost.gates;
    inv_weight = inv_ratio > 0.8 ? 4u : ( inv_ratio > 0.4 ? 2u : 1u );
  }
  cost.effective_area_score = static_cast<uint64_t>( cost.gates ) * 1024u + cost.effective_inverted_edges * inv_weight;
  return cost;
}

template<class Ntk>
bool is_mmig_objective_better( mmig_objective_cost<Ntk> const& lhs,
                               mmig_objective_cost<Ntk> const& rhs,
                               bool depth_mode,
                               bool minority_aware,
                               bool reward_minority_count )
{
  if ( !minority_aware )
  {
    if ( depth_mode )
    {
      if ( lhs.depth != rhs.depth )
      {
        return lhs.depth < rhs.depth;
      }
      if ( lhs.gates != rhs.gates )
      {
        return lhs.gates < rhs.gates;
      }
      if ( lhs.inverted_edges != rhs.inverted_edges )
      {
        return lhs.inverted_edges < rhs.inverted_edges;
      }
      return false;
    }

    if ( lhs.gates != rhs.gates )
    {
      return lhs.gates < rhs.gates;
    }
    if ( lhs.depth != rhs.depth )
    {
      return lhs.depth < rhs.depth;
    }
    if ( lhs.inverted_edges != rhs.inverted_edges )
    {
      return lhs.inverted_edges < rhs.inverted_edges;
    }
    return false;
  }

  if ( depth_mode )
  {
    if ( lhs.depth != rhs.depth )
    {
      return lhs.depth < rhs.depth;
    }
    if ( lhs.gates != rhs.gates )
    {
      return lhs.gates < rhs.gates;
    }
    if ( lhs.effective_inverted_edges != rhs.effective_inverted_edges )
    {
      return lhs.effective_inverted_edges < rhs.effective_inverted_edges;
    }
    if ( lhs.inverted_edges != rhs.inverted_edges )
    {
      return lhs.inverted_edges < rhs.inverted_edges;
    }
    if ( reward_minority_count && lhs.minority_gates != rhs.minority_gates )
    {
      return lhs.minority_gates > rhs.minority_gates;
    }
    return false;
  }

  if ( lhs.effective_area_score != rhs.effective_area_score )
  {
    return lhs.effective_area_score < rhs.effective_area_score;
  }
  if ( lhs.depth != rhs.depth )
  {
    return lhs.depth < rhs.depth;
  }
  if ( lhs.gates != rhs.gates )
  {
    return lhs.gates < rhs.gates;
  }
  if ( lhs.effective_inverted_edges != rhs.effective_inverted_edges )
  {
    return lhs.effective_inverted_edges < rhs.effective_inverted_edges;
  }
  if ( lhs.inverted_edges != rhs.inverted_edges )
  {
    return lhs.inverted_edges < rhs.inverted_edges;
  }
  if ( reward_minority_count && lhs.minority_gates != rhs.minority_gates )
  {
    return lhs.minority_gates > rhs.minority_gates;
  }
  return false;
}

template<class Ntk>
uint32_t auto_resub_candidate_budget( Ntk const& ntk, bool run_cec )
{
  if ( !run_cec )
  {
    return 0u;
  }
  auto const gates = ntk.num_gates();
  if ( gates <= 1000u )
  {
    return 128u;
  }
  if ( gates <= 5000u )
  {
    return 64u;
  }
  if ( gates <= 12000u )
  {
    return 24u;
  }
  return 8u;
}

template<class Ntk>
uint32_t auto_cut_limit( Ntk const& ntk, bool depth_mode )
{
  auto const gates = ntk.num_gates();
  if ( gates > 20000u )
  {
    return depth_mode ? 6u : 8u;
  }
  if ( gates > 8000u )
  {
    return depth_mode ? 8u : 9u;
  }
  if ( gates > 2500u )
  {
    return 10u;
  }
  return depth_mode ? 12u : 14u;
}

template<class Ntk>
uint32_t auto_refactor_max_pis( Ntk const& ntk, bool depth_mode )
{
  auto const gates = ntk.num_gates();
  if ( gates > 20000u )
  {
    return 4u;
  }
  if ( gates > 7000u )
  {
    return 5u;
  }
  return depth_mode ? 6u : 5u;
}

template<class Ntk>
uint32_t auto_exact_gate_limit( Ntk const& ntk, bool run_cec, bool depth_mode )
{
  (void)ntk;
  if ( run_cec )
  {
    return depth_mode ? 2500u : 3500u;
  }
  return depth_mode ? 5000u : 7000u;
}

template<class Ntk>
uint32_t auto_exact_min_inverted_edges( Ntk const& ntk )
{
  auto const gates = ntk.num_gates();
  if ( gates == 0u )
  {
    return 0u;
  }
  return std::max( 16u, gates / 10u );
}

inline std::vector<mmig_advanced_step> build_mmig_advanced_flow( mmig_optimizer_params const& ps )
{
  std::vector<mmig_advanced_step> steps;
  auto add_stage = [&]( mmig_advanced_stage stage ) {
    steps.push_back( mmig_advanced_step{ stage, std::nullopt, std::nullopt } );
  };
  auto add_resub = [&]( uint32_t cut, uint32_t inserts ) {
    steps.push_back( mmig_advanced_step{ mmig_advanced_stage::resub, cut, inserts } );
  };
  auto add_exact = [&]() {
    steps.push_back( mmig_advanced_step{ mmig_advanced_stage::exact_rewriting, std::nullopt, std::nullopt } );
  };
  auto add_balance = [&]() {
    steps.push_back( mmig_advanced_step{ mmig_advanced_stage::balancing, std::nullopt, std::nullopt } );
  };
  auto add_resub_with_depth_cap = [&]( uint32_t cut, uint32_t inserts ) {
    const auto capped_inserts = ( cut > 8u && inserts > 1u ) ? 1u : inserts;
    add_resub( cut, capped_inserts );
  };

  switch ( ps.advanced_flow )
  {
  case mmig_advanced_flow::round_robin:
    if ( ps.advanced_run_resub )
    {
      add_stage( mmig_advanced_stage::resub );
    }
    if ( ps.advanced_run_cut )
    {
      add_stage( mmig_advanced_stage::cut );
    }
    if ( ps.advanced_run_refactor )
    {
      add_stage( mmig_advanced_stage::refactor );
    }
    if ( ps.advanced_run_exact )
    {
      add_exact();
    }
    if ( ps.advanced_run_balance )
    {
      add_balance();
    }
    break;
  case mmig_advanced_flow::compress2rs:
    add_stage( mmig_advanced_stage::mighty_area );
    add_resub_with_depth_cap( 6u, 1u );
    add_stage( mmig_advanced_stage::cut );
    add_resub_with_depth_cap( 6u, 2u );
    add_stage( mmig_advanced_stage::refactor );
    add_resub_with_depth_cap( 8u, 1u );
    add_stage( mmig_advanced_stage::mighty_area );
    add_resub_with_depth_cap( 8u, 2u );
    add_stage( mmig_advanced_stage::cut );
    add_resub_with_depth_cap( 10u, 1u );
    add_stage( mmig_advanced_stage::cut_zero );
    add_resub_with_depth_cap( 10u, 2u );
    add_stage( mmig_advanced_stage::mighty_area );
    add_resub_with_depth_cap( 12u, 1u );
    add_stage( mmig_advanced_stage::refactor_zero );
    add_resub_with_depth_cap( 12u, 2u );
    add_stage( mmig_advanced_stage::cut_zero );
    add_stage( mmig_advanced_stage::mighty_area );
    if ( ps.advanced_run_exact )
    {
      add_exact();
    }
    if ( ps.advanced_run_balance )
    {
      add_balance();
    }
    break;
  case mmig_advanced_flow::paper2019:
    /* Exact sequence from Riener et al., "Scalable Generic Logic Synthesis:
       One Approach to Rule Them All", Section 3.1:
       bz; rs -c 6; rw; rs -c 6 -d 2; rf; rs -c 8; bz;
       rs -c 8 -d 2; rw; rs -c 10; rwz; rs -c 10 -d 2;
       bz; rs -c 12; rfz; rs -c 12 -d 2; rwz; bz. */
    add_stage( mmig_advanced_stage::mighty_area );
    add_resub( 6u, 1u );
    add_stage( mmig_advanced_stage::cut );
    add_resub( 6u, 2u );
    add_stage( mmig_advanced_stage::refactor );
    add_resub( 8u, 1u );
    add_stage( mmig_advanced_stage::mighty_area );
    add_resub( 8u, 2u );
    add_stage( mmig_advanced_stage::cut );
    add_resub( 10u, 1u );
    add_stage( mmig_advanced_stage::cut_zero );
    add_resub( 10u, 2u );
    add_stage( mmig_advanced_stage::mighty_area );
    add_resub( 12u, 1u );
    add_stage( mmig_advanced_stage::refactor_zero );
    add_resub( 12u, 2u );
    add_stage( mmig_advanced_stage::cut_zero );
    add_stage( mmig_advanced_stage::mighty_area );
    break;
  case mmig_advanced_flow::epfl:
    add_stage( mmig_advanced_stage::mighty_area );
    add_resub_with_depth_cap( 6u, 1u );
    add_stage( mmig_advanced_stage::cut );
    add_resub_with_depth_cap( 6u, 2u );
    add_stage( mmig_advanced_stage::refactor );
    add_resub_with_depth_cap( 8u, 1u );
    add_stage( mmig_advanced_stage::mighty_area );
    add_stage( mmig_advanced_stage::cut_zero );
    add_resub_with_depth_cap( 8u, 2u );
    add_stage( mmig_advanced_stage::refactor_zero );
    add_resub_with_depth_cap( 10u, 1u );
    add_stage( mmig_advanced_stage::cut );
    add_stage( mmig_advanced_stage::mighty_area );
    if ( ps.advanced_run_exact )
    {
      add_exact();
    }
    if ( ps.advanced_run_balance )
    {
      add_balance();
    }
    add_stage( mmig_advanced_stage::mighty_area );
    break;
  case mmig_advanced_flow::dac19_default:
  case mmig_advanced_flow::dac19_area:
    add_resub( 8u, 2u );
    add_stage( mmig_advanced_stage::cut );
    add_resub( 8u, 2u );
    add_stage( mmig_advanced_stage::refactor );
    add_stage( mmig_advanced_stage::cut );
    if ( ps.advanced_allow_zero_gain )
    {
      add_stage( mmig_advanced_stage::cut_zero );
    }
    add_resub( 8u, 2u );
    add_stage( ps.advanced_allow_zero_gain ? mmig_advanced_stage::refactor_zero : mmig_advanced_stage::refactor );
    add_stage( mmig_advanced_stage::cut );
    if ( ps.advanced_run_exact )
    {
      add_exact();
    }
    if ( ps.advanced_run_balance )
    {
      add_balance();
    }
    break;
  case mmig_advanced_flow::legacy:
    add_resub( 8u, 2u );
    add_stage( mmig_advanced_stage::cut );
    add_resub( 8u, 2u );
    add_stage( mmig_advanced_stage::refactor );
    add_stage( mmig_advanced_stage::cut );
    if ( ps.advanced_allow_zero_gain )
    {
      add_stage( mmig_advanced_stage::cut_zero );
    }
    add_resub( 8u, 2u );
    add_stage( mmig_advanced_stage::refactor );
    if ( ps.advanced_allow_zero_gain )
    {
      add_stage( mmig_advanced_stage::refactor_zero );
    }
    add_stage( mmig_advanced_stage::cut );
    if ( ps.advanced_run_exact )
    {
      add_exact();
    }
    if ( ps.advanced_run_balance )
    {
      add_balance();
    }
    break;
  }

  return steps;
}

template<class Ntk>
void run_mmig_advanced_rounds( Ntk& ntk, mmig_optimizer_params const& ps, mmig_optimizer_stats& st, bool depth_mode )
{
  if ( !ps.enable_advanced )
  {
    return;
  }

  auto const steps = build_mmig_advanced_flow( ps );
  if ( steps.empty() )
  {
    return;
  }

  auto const prefix = depth_mode ? "mmig-depth:" : "mmig-area:";
  auto const rounds = std::max( 1u, ps.advanced_rounds );
  auto const stagnation_limit = ps.advanced_stop_on_stagnation ? ps.advanced_max_stagnation_rounds : 0u;
  auto const epfl_profile = ps.advanced_flow == mmig_advanced_flow::epfl;
  uint32_t stagnation_rounds = 0u;
  auto const overall_before = compute_mmig_objective_cost( ntk );
  st.advanced_depth_before = overall_before.depth;
  st.advanced_gates_before = overall_before.gates;
  st.advanced_majority_before = overall_before.majority_gates;
  st.advanced_minority_before = overall_before.minority_gates;
  st.advanced_inverted_edges_before = overall_before.inverted_edges;
  st.advanced_effective_inverted_edges_before = overall_before.effective_inverted_edges;
  st.advanced_effective_area_score_before = overall_before.effective_area_score;

  for ( uint32_t round = 1u; round <= rounds; ++round )
  {
    ++st.advanced_rounds_executed;
    auto const before = compute_mmig_objective_cost( ntk );
    auto backup = ntk.clone();

    uint32_t step_index = 0u;
    for ( auto const& step : steps )
    {
      ++step_index;
      switch ( step.stage )
      {
      case mmig_advanced_stage::mighty_area:
        run_mmig_stage( ntk, prefix + std::string( "adv-mighty:r" ) + std::to_string( round ) + ":s" + std::to_string( step_index ), ps, st.cec, [&]( auto& candidate ) {
          mmig_algebraic_rewriting_params rw_ps;
          rw_ps.max_iterations = ps.max_iterations;
          rw_ps.allow_area_increase = false;
          rw_ps.allow_non_improving = false;
          rw_ps.enable_sr5 = ps.enable_sr5 || ( ps.auto_sr5_with_cec && ps.run_cec ); /* P7 */
          rw_ps.rank_candidates_globally = ps.rewrite_rank_candidates;
          rw_ps.normalize_complemented_inner = ps.rewrite_normalize_complemented_inner;
          rw_ps.max_candidates_per_iteration = epfl_profile ? 512u : rw_ps.max_candidates_per_iteration;
          rw_ps.use_mffc_cost = epfl_profile && !depth_mode;
          rw_ps.verbose = ps.verbose;

          mmig_algebraic_rewriting_stats round_st{};
          mmig_algebraic_rewriting( candidate, rw_ps, &round_st );
          accumulate_mmig_algebraic_rewriting_stats( st.advanced_mighty, round_st );
          candidate = cleanup_dangling( candidate );
        } );
        break;
      case mmig_advanced_stage::resub:
        if ( !ps.advanced_run_resub )
        {
          break;
        }
        run_mmig_stage( ntk, prefix + std::string( "adv-resub:r" ) + std::to_string( round ) + ":s" + std::to_string( step_index ), ps, st.cec, [&]( auto& candidate ) {
          auto const tuned = ps.advanced_tuned_policy;
          auto const gates = candidate.num_gates();
          mmig_resubstitution_params rsp;
          auto const default_max_pis = tuned ? ( gates > 15000u ? 6u : 8u ) : 8u;
          rsp.resub_ps.max_pis = step.rs_cut_override.value_or( default_max_pis );
          rsp.resub_ps.max_divisors = tuned ? ( gates > 15000u ? 96u : ( gates > 5000u ? 128u : 160u ) ) : 128u;
          if ( epfl_profile )
          {
            rsp.resub_ps.max_divisors = gates > 15000u ? 128u : ( gates > 5000u ? 192u : 256u );
          }
          auto const default_max_inserts = tuned ? ( gates > 10000u ? 1u : 2u ) : 2u;
          rsp.resub_ps.max_inserts = step.rs_insert_override.value_or( default_max_inserts );
          rsp.resub_ps.preserve_depth = depth_mode;
          rsp.resub_ps.use_dont_cares = ps.enable_dont_cares && ( candidate.num_gates() <= ps.dont_care_gate_threshold ); /* P2 */
          rsp.resub_ps.window_size = tuned ? ( depth_mode ? ( gates > 10000u ? 6u : 8u ) : ( gates > 10000u ? 8u : 10u ) ) : ( depth_mode ? 8u : 10u );
          if ( epfl_profile )
          {
            rsp.resub_ps.window_size = depth_mode ? ( gates > 15000u ? 8u : 10u ) : ( gates > 15000u ? 10u : 12u );
          }
          rsp.resub_ps.rank_mig_candidates = ps.resub_rank_candidates;
          rsp.resub_ps.mig_candidate_pool_size = epfl_profile ? std::max( ps.resub_candidate_pool_size, 40u ) : ps.resub_candidate_pool_size;
          rsp.resub_ps.mig_full_score_limit = epfl_profile ? std::max( ps.resub_full_score_limit, 10u ) : ps.resub_full_score_limit;
          rsp.enable_candidate_cec = ps.run_cec;
          rsp.verify_mixed_risk_only = true;
          if ( rsp.enable_candidate_cec )
          {
            rsp.max_candidate_cec_checks = ps.advanced_resub_candidate_cec_budget > 0u ? ps.advanced_resub_candidate_cec_budget : auto_resub_candidate_budget( candidate, ps.run_cec );
            if ( epfl_profile )
            {
              auto const epfl_budget = gates > 15000u ? 16u : ( gates > 5000u ? 32u : 96u );
              rsp.max_candidate_cec_checks = std::max( rsp.max_candidate_cec_checks, epfl_budget );
            }
            rsp.skip_candidate_when_cec_budget_exhausted = true;
          }
          rsp.verbose = ps.verbose;

          mmig_resubstitution_stats round_st{};
          mmig_resubstitution( candidate, rsp, &round_st );
          accumulate_mmig_resubstitution_stats( st.advanced_resub, round_st );
          candidate = cleanup_dangling( candidate );
        } );
        break;
      case mmig_advanced_stage::cut:
      case mmig_advanced_stage::cut_zero:
        if ( !ps.advanced_run_cut )
        {
          break;
        }
        run_mmig_stage( ntk, prefix + std::string( step.stage == mmig_advanced_stage::cut_zero ? "adv-cutz:r" : "adv-cut:r" ) + std::to_string( round ) + ":s" + std::to_string( step_index ), ps, st.cec, [&]( auto& candidate ) {
          auto const tuned = ps.advanced_tuned_policy;
          mmig_cut_rewriting_params crp;
          crp.cut_ps.cut_enumeration_ps.cut_size = 4u;
          crp.cut_ps.cut_enumeration_ps.cut_limit = tuned ? auto_cut_limit( candidate, depth_mode ) : 10u;
          if ( epfl_profile )
          {
            crp.cut_ps.cut_enumeration_ps.cut_limit = std::min( 18u, crp.cut_ps.cut_enumeration_ps.cut_limit + 4u );
          }
          crp.cut_ps.allow_zero_gain = ( step.stage == mmig_advanced_stage::cut_zero );
          crp.cut_ps.preserve_depth = depth_mode;
          crp.depth_rollback_on_regression = depth_mode;
          crp.verbose = ps.verbose;

          mmig_cut_rewriting_stats round_st{};
          mmig_cut_rewriting( candidate, crp, &round_st );
          accumulate_mmig_cut_rewriting_stats( st.advanced_cut, round_st );
          candidate = cleanup_dangling( candidate );
        } );
        break;
      case mmig_advanced_stage::refactor:
      case mmig_advanced_stage::refactor_zero:
        if ( !ps.advanced_run_refactor )
        {
          break;
        }
        run_mmig_stage( ntk, prefix + std::string( step.stage == mmig_advanced_stage::refactor_zero ? "adv-refactorz:r" : "adv-refactor:r" ) + std::to_string( round ) + ":s" + std::to_string( step_index ), ps, st.cec, [&]( auto& candidate ) {
          auto const tuned = ps.advanced_tuned_policy;
          auto const gates = candidate.num_gates();
          mmig_refactoring_params rfp;
          rfp.refac_ps.max_pis = tuned ? auto_refactor_max_pis( candidate, depth_mode ) : 6u;
          rfp.refac_ps.use_reconvergence_cut = !tuned || gates <= 20000u;
          rfp.refac_ps.allow_zero_gain = ( step.stage == mmig_advanced_stage::refactor_zero );
          rfp.refac_ps.use_dont_cares = ps.enable_dont_cares && ( candidate.num_gates() <= ps.dont_care_gate_threshold ); /* P2 */
          rfp.preserve_depth = depth_mode;
          rfp.depth_rollback_on_regression = depth_mode;
          rfp.verbose = ps.verbose;

          mmig_refactoring_stats round_st{};
          mmig_refactoring( candidate, rfp, &round_st );
          accumulate_mmig_refactoring_stats( st.advanced_refactor, round_st );
          candidate = cleanup_dangling( candidate );
        } );
        break;
      case mmig_advanced_stage::exact_rewriting:
        if ( !ps.advanced_run_exact )
        {
          break;
        }
        run_mmig_stage( ntk, prefix + std::string( "adv-exact:r" ) + std::to_string( round ) + ":s" + std::to_string( step_index ), ps, st.cec, [&]( auto& candidate ) {
          auto const tuned = ps.advanced_tuned_policy;
          mmig_exact_rewriting_params erp;
          erp.rw_ps.cut_enumeration_ps.cut_size = 4u;
          erp.rw_ps.cut_enumeration_ps.cut_limit = 10u;
          erp.rw_ps.allow_zero_gain = ps.advanced_allow_zero_gain;
          erp.rw_ps.use_dont_cares = ps.enable_dont_cares && ( candidate.num_gates() <= ps.dont_care_gate_threshold ); /* P2 */
          erp.preserve_depth = depth_mode;
          erp.depth_rollback_on_regression = depth_mode;
          erp.gate_enable = tuned;
          erp.gate_max_gates = ps.advanced_exact_gate_limit > 0u ? ps.advanced_exact_gate_limit : auto_exact_gate_limit( candidate, ps.run_cec, depth_mode );
          if ( epfl_profile && ps.advanced_exact_gate_limit == 0u )
          {
            erp.gate_max_gates = std::max( erp.gate_max_gates, depth_mode ? 5000u : 7000u );
          }
          erp.gate_min_depth = ps.advanced_exact_min_depth > 0u ? ps.advanced_exact_min_depth : ( depth_mode ? 8u : 0u );
          erp.gate_min_inverted_edges = ps.advanced_exact_min_inverted_edges > 0u ? ps.advanced_exact_min_inverted_edges : auto_exact_min_inverted_edges( candidate );
          if ( epfl_profile && ps.advanced_exact_min_inverted_edges == 0u )
          {
            erp.gate_min_inverted_edges = std::max( 8u, candidate.num_gates() / 20u );
          }
          erp.gate_require_minority = ps.advanced_exact_require_minority;
          erp.verbose = ps.verbose;

          mmig_exact_rewriting_stats round_st{};
          mmig_exact_rewriting( candidate, erp, &round_st );
          accumulate_mmig_exact_rewriting_stats( st.advanced_exact, round_st );
          candidate = cleanup_dangling( candidate );
        } );
        break;
      case mmig_advanced_stage::balancing:
        if ( !ps.advanced_run_balance )
        {
          break;
        }
        run_mmig_stage( ntk, prefix + std::string( "adv-balance:r" ) + std::to_string( round ) + ":s" + std::to_string( step_index ), ps, st.cec, [&]( auto& candidate ) {
          mmig_balancing_params mbp;
          mbp.bal_ps.cut_enumeration_ps.cut_size = 4u;
          mbp.bal_ps.cut_enumeration_ps.cut_limit = 16u;
          mbp.bal_ps.only_on_critical_path = depth_mode;
          mbp.preserve_depth = depth_mode;
          mbp.depth_rollback_on_regression = depth_mode;
          mbp.verbose = ps.verbose;

          mmig_balancing_stats round_st{};
          mmig_balancing( candidate, mbp, &round_st );
          accumulate_mmig_balancing_stats( st.advanced_balance, round_st );
          candidate = cleanup_dangling( candidate );
        } );
        break;
      }
    }

    /* P3: Inter-round minority re-seeding */
    auto const reseed_interval = ps.interleaved_seeding_interval;
    if ( ps.enable_interleaved_seeding &&
         reseed_interval > 0u &&
         round > 1u &&
         ( round % reseed_interval == 0u ) )
    {
      run_mmig_stage( ntk, prefix + std::string( "adv-reseed:r" ) + std::to_string( round ), ps, st.cec, [&]( auto& candidate ) {
        mmig_minority_seeding_params seed_ps;
        seed_ps.mode = mmig_seed_mode::level1; /* fast L1-only */
        seed_ps.max_candidates = 12u;
        seed_ps.max_rounds = 1u;
        seed_ps.allow_area_increase = false;
        seed_ps.verify_each_candidate = ps.run_cec;
        mmig_minority_seeding( candidate, seed_ps, nullptr );
        candidate = cleanup_dangling( candidate );
      } );
    }

    auto after = compute_mmig_objective_cost( ntk );
    auto const improved = is_mmig_objective_better( after, before, depth_mode, ps.advanced_minority_aware_cost, ps.advanced_reward_minority_count );
    auto const regressed = is_mmig_objective_better( before, after, depth_mode, ps.advanced_minority_aware_cost, ps.advanced_reward_minority_count );

    if ( regressed && ps.advanced_rollback_on_objective_regression )
    {
      ntk = std::move( backup );
      after = before;
      ++st.advanced_rounds_rolled_back;
    }

    if ( improved )
    {
      ++st.advanced_rounds_improved;
      stagnation_rounds = 0u;
    }
    else
    {
      ++stagnation_rounds;
      if ( stagnation_limit > 0u && stagnation_rounds >= stagnation_limit )
      {
        ++st.advanced_stagnation_stops;
        break;
      }
    }
  }

  auto const overall_after = compute_mmig_objective_cost( ntk );
  st.advanced_depth_after = overall_after.depth;
  st.advanced_gates_after = overall_after.gates;
  st.advanced_majority_after = overall_after.majority_gates;
  st.advanced_minority_after = overall_after.minority_gates;
  st.advanced_inverted_edges_after = overall_after.inverted_edges;
  st.advanced_effective_inverted_edges_after = overall_after.effective_inverted_edges;
  st.advanced_effective_area_score_after = overall_after.effective_area_score;
}

template<class Ntk>
void run_mmig_depth_pipeline( Ntk& ntk, mmig_optimizer_params const& ps, mmig_optimizer_stats& st )
{
  auto run_pre_stages = [&]() {
    if ( ps.enable_minority_seeding )
    {
      run_mmig_stage( ntk, "mmig-depth:seed-pre", ps, st.cec, [&]( auto& candidate ) {
        mmig_minority_seeding_params seed_ps;
        seed_ps.mode = ps.seeding_mode;
        seed_ps.max_candidates = ps.seeding_max_candidates;
        seed_ps.max_rounds = ps.seeding_max_rounds;
        seed_ps.enable_incremental_reranking = ps.seeding_incremental_reranking;
        seed_ps.allow_area_increase = true;
        seed_ps.allow_non_improving = true;
        seed_ps.verify_each_candidate = ps.run_cec;
        seed_ps.prefer_l1_when_verifying = ps.run_cec;
        seed_ps.min_l1_for_verification = 1;
        seed_ps.verbose = ps.verbose;
        mmig_minority_seeding( candidate, seed_ps, &st.seeding_pre );
        candidate = cleanup_dangling( candidate );
      } );
    }

    run_mmig_stage( ntk, "mmig-depth:rewrite-pre", ps, st.cec, [&]( auto& candidate ) {
      mmig_algebraic_rewriting_params rw_ps;
      rw_ps.max_iterations = ps.max_iterations;
      rw_ps.allow_area_increase = true;
      rw_ps.allow_non_improving = false;
      rw_ps.enable_sr5 = ps.enable_sr5 || ( ps.auto_sr5_with_cec && ps.run_cec ); /* P7 */
      rw_ps.rank_candidates_globally = ps.rewrite_rank_candidates;
      rw_ps.normalize_complemented_inner = ps.rewrite_normalize_complemented_inner;
      rw_ps.max_candidates_per_iteration = ps.advanced_flow == mmig_advanced_flow::epfl ? 512u : rw_ps.max_candidates_per_iteration;
      rw_ps.use_mffc_cost = false;
      rw_ps.verbose = ps.verbose;
      mmig_algebraic_rewriting( candidate, rw_ps, &st.rewrite_pre );
      candidate = cleanup_dangling( candidate );
    } );

    run_mmig_stage( ntk, "mmig-depth:inv-prop-pre", ps, st.cec, [&]( auto& candidate ) {
      mmig_inv_propagation_params inv_ps;
      inv_ps.enable_dual_inversion = ps.enable_dual_inv_propagation; /* P1 */
      mmig_inv_propagation( candidate, inv_ps, &st.inv_propagation_pre );
      candidate = cleanup_dangling( candidate );
    } );

    run_mmig_stage( ntk, "mmig-depth:inv-opt-pre", ps, st.cec, [&]( auto& candidate ) {
      mmig_inv_optimization( candidate, &st.inv_optimization_pre );
      candidate = cleanup_dangling( candidate );
    } );
  };

  auto run_post_stages = [&]() {
    run_mmig_stage( ntk, "mmig-depth:rewrite-post", ps, st.cec, [&]( auto& candidate ) {
      mmig_algebraic_rewriting_params rw_ps;
      rw_ps.max_iterations = ps.max_iterations;
      rw_ps.allow_area_increase = false;
      rw_ps.allow_non_improving = false;
      rw_ps.enable_sr5 = ps.enable_sr5 || ( ps.auto_sr5_with_cec && ps.run_cec ); /* P7 */
      rw_ps.rank_candidates_globally = ps.rewrite_rank_candidates;
      rw_ps.normalize_complemented_inner = ps.rewrite_normalize_complemented_inner;
      rw_ps.max_candidates_per_iteration = ps.advanced_flow == mmig_advanced_flow::epfl ? 512u : rw_ps.max_candidates_per_iteration;
      rw_ps.use_mffc_cost = false;
      rw_ps.verbose = ps.verbose;
      mmig_algebraic_rewriting( candidate, rw_ps, &st.rewrite_post );
      candidate = cleanup_dangling( candidate );
    } );

    run_mmig_stage( ntk, "mmig-depth:cone-flip-post", ps, st.cec, [&]( auto& candidate ) {
      mmig_cone_polarity_flip_params cpf_ps;
      mmig_cone_polarity_flip( candidate, cpf_ps );
      candidate = cleanup_dangling( candidate );
    } );

    run_mmig_stage( ntk, "mmig-depth:inv-prop-post", ps, st.cec, [&]( auto& candidate ) {
      mmig_inv_propagation_params inv_ps;
      inv_ps.enable_dual_inversion = ps.enable_dual_inv_propagation; /* P1 */
      mmig_inv_propagation( candidate, inv_ps, &st.inv_propagation_post );
      candidate = cleanup_dangling( candidate );
    } );

    run_mmig_stage( ntk, "mmig-depth:inv-opt-post", ps, st.cec, [&]( auto& candidate ) {
      mmig_inv_optimization( candidate, &st.inv_optimization_post );
      candidate = cleanup_dangling( candidate );
    } );
  };

  if ( ps.stage_mode == mmig_stage_mode::pre || ps.stage_mode == mmig_stage_mode::both )
  {
    run_pre_stages();
  }
  if ( ps.stage_mode == mmig_stage_mode::post || ps.stage_mode == mmig_stage_mode::both )
  {
    run_post_stages();
  }

  run_mmig_advanced_rounds( ntk, ps, st, true );
}

template<class Ntk>
void run_mmig_area_pipeline( Ntk& ntk, mmig_optimizer_params const& ps, mmig_optimizer_stats& st )
{
  auto run_pre_stages = [&]() {
    if ( ps.enable_minority_seeding )
    {
      run_mmig_stage( ntk, "mmig-area:seed-pre", ps, st.cec, [&]( auto& candidate ) {
        mmig_minority_seeding_params seed_ps;
        seed_ps.mode = ps.seeding_mode;
        seed_ps.max_candidates = ps.seeding_max_candidates;
        seed_ps.max_rounds = ps.seeding_max_rounds;
        seed_ps.enable_incremental_reranking = ps.seeding_incremental_reranking;
        seed_ps.allow_area_increase = true;
        seed_ps.allow_non_improving = true;
        seed_ps.verify_each_candidate = ps.run_cec;
        seed_ps.prefer_l1_when_verifying = ps.run_cec;
        seed_ps.min_l1_for_verification = 1;
        seed_ps.verbose = ps.verbose;
        mmig_minority_seeding( candidate, seed_ps, &st.seeding_pre );
        candidate = cleanup_dangling( candidate );
      } );
    }

    run_mmig_stage( ntk, "mmig-area:rewrite-pre", ps, st.cec, [&]( auto& candidate ) {
      mmig_algebraic_rewriting_params rw_ps;
      rw_ps.max_iterations = ps.max_iterations;
      rw_ps.allow_area_increase = false;
      rw_ps.allow_non_improving = false;
      rw_ps.enable_sr5 = ps.enable_sr5 || ( ps.auto_sr5_with_cec && ps.run_cec ); /* P7 */
      rw_ps.rank_candidates_globally = ps.rewrite_rank_candidates;
      rw_ps.normalize_complemented_inner = ps.rewrite_normalize_complemented_inner;
      rw_ps.max_candidates_per_iteration = ps.advanced_flow == mmig_advanced_flow::epfl ? 512u : rw_ps.max_candidates_per_iteration;
      rw_ps.use_mffc_cost = ps.advanced_flow == mmig_advanced_flow::epfl;
      rw_ps.verbose = ps.verbose;
      mmig_algebraic_rewriting( candidate, rw_ps, &st.rewrite_pre );
      candidate = cleanup_dangling( candidate );
    } );

    run_mmig_stage( ntk, "mmig-area:inv-prop-pre", ps, st.cec, [&]( auto& candidate ) {
      mmig_inv_propagation_params inv_ps;
      inv_ps.enable_dual_inversion = ps.enable_dual_inv_propagation; /* P1 */
      mmig_inv_propagation( candidate, inv_ps, &st.inv_propagation_pre );
      candidate = cleanup_dangling( candidate );
    } );

    run_mmig_stage( ntk, "mmig-area:inv-opt-pre", ps, st.cec, [&]( auto& candidate ) {
      mmig_inv_optimization( candidate, &st.inv_optimization_pre );
      candidate = cleanup_dangling( candidate );
    } );
  };

  auto run_post_stages = [&]() {
    run_mmig_stage( ntk, "mmig-area:rewrite-post", ps, st.cec, [&]( auto& candidate ) {
      mmig_algebraic_rewriting_params rw_ps;
      rw_ps.max_iterations = ps.max_iterations;
      rw_ps.allow_area_increase = false;
      rw_ps.allow_non_improving = false;
      rw_ps.enable_sr5 = ps.enable_sr5 || ( ps.auto_sr5_with_cec && ps.run_cec ); /* P7 */
      rw_ps.rank_candidates_globally = ps.rewrite_rank_candidates;
      rw_ps.normalize_complemented_inner = ps.rewrite_normalize_complemented_inner;
      rw_ps.max_candidates_per_iteration = ps.advanced_flow == mmig_advanced_flow::epfl ? 512u : rw_ps.max_candidates_per_iteration;
      rw_ps.use_mffc_cost = ps.advanced_flow == mmig_advanced_flow::epfl;
      rw_ps.verbose = ps.verbose;
      mmig_algebraic_rewriting( candidate, rw_ps, &st.rewrite_post );
      candidate = cleanup_dangling( candidate );
    } );

    run_mmig_stage( ntk, "mmig-area:cone-flip-post", ps, st.cec, [&]( auto& candidate ) {
      mmig_cone_polarity_flip_params cpf_ps;
      mmig_cone_polarity_flip( candidate, cpf_ps );
      candidate = cleanup_dangling( candidate );
    } );

    run_mmig_stage( ntk, "mmig-area:inv-prop-post", ps, st.cec, [&]( auto& candidate ) {
      mmig_inv_propagation_params inv_ps;
      inv_ps.enable_dual_inversion = ps.enable_dual_inv_propagation; /* P1 */
      mmig_inv_propagation( candidate, inv_ps, &st.inv_propagation_post );
      candidate = cleanup_dangling( candidate );
    } );

    run_mmig_stage( ntk, "mmig-area:inv-opt-post", ps, st.cec, [&]( auto& candidate ) {
      mmig_inv_optimization( candidate, &st.inv_optimization_post );
      candidate = cleanup_dangling( candidate );
    } );
  };

  if ( ps.stage_mode == mmig_stage_mode::pre || ps.stage_mode == mmig_stage_mode::both )
  {
    run_pre_stages();
  }
  if ( ps.stage_mode == mmig_stage_mode::post || ps.stage_mode == mmig_stage_mode::both )
  {
    run_post_stages();
  }

  run_mmig_advanced_rounds( ntk, ps, st, false );
}

} // namespace detail

template<class Ntk>
void mmig_depth_optimization( Ntk& ntk, mmig_optimizer_params const& ps = {}, mmig_optimizer_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( std::is_same_v<typename Ntk::base_type, mig_network>, "Ntk is not MIG-based" );
  static_assert( has_is_min_v<Ntk>, "Ntk does not implement is_min" );
  static_assert( has_create_min_v<Ntk>, "Ntk does not implement create_min" );

  mmig_optimizer_stats st{};
  {
    stopwatch t( st.time_total );
    detail::run_mmig_depth_pipeline( ntk, ps, st );
  }

  if ( pst != nullptr )
  {
    *pst = st;
  }
}

template<class Ntk>
void mmig_area_optimization( Ntk& ntk, mmig_optimizer_params const& ps = {}, mmig_optimizer_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( std::is_same_v<typename Ntk::base_type, mig_network>, "Ntk is not MIG-based" );
  static_assert( has_is_min_v<Ntk>, "Ntk does not implement is_min" );
  static_assert( has_create_min_v<Ntk>, "Ntk does not implement create_min" );

  mmig_optimizer_stats st{};
  {
    stopwatch t( st.time_total );
    detail::run_mmig_area_pipeline( ntk, ps, st );
  }

  if ( pst != nullptr )
  {
    *pst = st;
  }
}

} // namespace mockturtle
