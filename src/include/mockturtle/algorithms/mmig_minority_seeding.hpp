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
  \file mmig_minority_seeding.hpp
  \brief mMIG minority seeding from majority candidates (L1/L2/L3)
*/

#pragma once

#include "cleanup.hpp"
#include "equivalence_checking.hpp"
#include "miter.hpp"

#include "../networks/mig.hpp"
#include "../traits.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/depth_view.hpp"
#include "../views/fanout_view.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <optional>
#include <type_traits>
#include <vector>

namespace mockturtle
{

enum class mmig_seed_mode : uint8_t
{
  level1,
  level2,
  level3,
  both
};

struct mmig_minority_seeding_params
{
  mmig_seed_mode mode{ mmig_seed_mode::both };
  uint32_t max_candidates{ 24u };
  uint32_t max_rounds{ 1u };
  bool enable_input_complement_form{ false };
  bool verify_each_candidate{ false };
  bool prefer_l1_when_verifying{ true };
  int32_t min_l1_for_verification{ 1 };
  bool allow_area_increase{ true };
  bool allow_non_improving{ true };
  /* After each accepted seed within a round, recompute L1/L2/L3 scores for
   * all remaining candidates and re-sort.  Fixes the stale-score problem at
   * the cost of O(candidates) re-scoring per applied seed. */
  bool enable_incremental_reranking{ false };
  bool verbose{ false };
};

struct mmig_minority_seeding_stats
{
  stopwatch<>::duration time_total{ 0 };
  uint32_t num_rounds{ 0 };
  uint32_t num_scored{ 0 };
  uint32_t num_evaluated{ 0 };
  uint32_t num_selected{ 0 };
  uint32_t num_applied{ 0 };
  uint32_t num_equivalence_checks{ 0 };
  uint32_t num_equivalence_rejected{ 0 };
  uint32_t num_filtered_low_l1{ 0 };
  int32_t total_depth_gain{ 0 };
  int32_t total_area_gain{ 0 };
  int32_t total_inv_gain{ 0 };
  uint32_t num_form_output_complement{ 0 };
  uint32_t num_form_input_complement{ 0 };
};

namespace detail
{

enum class mmig_seed_form : uint8_t
{
  output_complement,
  input_complement
};

template<class Ntk>
struct mmig_seed_cost
{
  int32_t depth{ 0 };
  int32_t gates{ 0 };
  int32_t inverted_edges{ 0 };
};

template<class Ntk>
struct mmig_seed_candidate
{
  node<Ntk> n{};
  int32_t l1_score{ 0 };
  int32_t l2_score{ 0 };
  int32_t l3_score{ 0 };
};

template<class Ntk>
struct mmig_seed_eval
{
  node<Ntk> n{};
  mmig_seed_form form{ mmig_seed_form::output_complement };
  mmig_seed_cost<Ntk> after{};
  int32_t depth_gain{ 0 };
  int32_t area_gain{ 0 };
  int32_t inv_gain{ 0 };
};

template<class Ntk>
class mmig_minority_seeding_impl
{
public:
  mmig_minority_seeding_impl( Ntk& ntk, mmig_minority_seeding_params const& ps, mmig_minority_seeding_stats& st )
      : ntk( ntk ),
        ps( ps ),
        st( st )
  {
  }

  void run()
  {
    stopwatch t( st.time_total );

    for ( uint32_t round = 0u; round < ps.max_rounds; ++round )
    {
      ++st.num_rounds;

      auto candidates = rank_candidates();
      if ( candidates.empty() )
        break;

      if ( candidates.size() > ps.max_candidates )
        candidates.resize( ps.max_candidates );
      st.num_selected += static_cast<uint32_t>( candidates.size() );

      auto before = network_cost( ntk );
      bool applied_in_round = false;

      for ( std::size_t i = 0u; i < candidates.size(); ++i )
      {
        auto& cand = candidates[i];
        if ( ntk.is_dead( cand.n ) || !ntk.is_maj( cand.n ) )
          continue;
        if ( ps.verify_each_candidate && cand.l1_score < ps.min_l1_for_verification )
        {
          ++st.num_filtered_low_l1;
          continue;
        }

        auto eval = evaluate_node( cand.n, before );
        if ( !eval.has_value() )
          continue;
        if ( !accept( *eval ) )
          continue;
        if ( ps.verify_each_candidate && !verify_candidate( eval->n, eval->form ) )
        {
          ++st.num_equivalence_rejected;
          continue;
        }

        apply_eval( *eval );
        ++st.num_applied;
        st.total_depth_gain += eval->depth_gain;
        st.total_area_gain += eval->area_gain;
        st.total_inv_gain += eval->inv_gain;
        if ( eval->form == mmig_seed_form::output_complement )
          ++st.num_form_output_complement;
        else
          ++st.num_form_input_complement;

        applied_in_round = true;
        before = network_cost( ntk );

        /* Incremental reranking: recompute L1/L2/L3 for remaining candidates
         * so that downstream scores reflect the already-applied seeds.  This
         * fixes the stale-score problem: without it, the 24th candidate's
         * scores were computed before any seeds were applied. */
        if ( ps.enable_incremental_reranking && i + 1u < candidates.size() )
        {
          fanout_view<Ntk> fntk{ ntk };
          for ( std::size_t j = i + 1u; j < candidates.size(); ++j )
          {
            if ( ntk.is_dead( candidates[j].n ) )
              continue;
            candidates[j].l1_score = gain_level1( fntk, candidates[j].n );
            candidates[j].l2_score = gain_level2( fntk, candidates[j].n, candidates[j].l1_score );
            candidates[j].l3_score = gain_level3( fntk, candidates[j].n, candidates[j].l2_score );
          }
          std::stable_sort(
              candidates.begin() + static_cast<std::ptrdiff_t>( i ) + 1,
              candidates.end(),
              [&]( auto const& a, auto const& b ) { return cmp_candidates( a, b ); } );
        }
      }

      if ( !applied_in_round )
        break;

      ntk = cleanup_dangling( ntk );
    }
  }

private:
  std::vector<mmig_seed_candidate<Ntk>> rank_candidates()
  {
    fanout_view<Ntk> fntk{ ntk };
    std::vector<mmig_seed_candidate<Ntk>> ranked;
    ranked.reserve( ntk.num_gates() );

    ntk.foreach_gate( [&]( auto const& n ) {
      if ( !ntk.is_maj( n ) || ntk.is_dead( n ) )
      {
        return;
      }

      mmig_seed_candidate<Ntk> cand{};
      cand.n = n;
      cand.l1_score = gain_level1( fntk, n );
      cand.l2_score = gain_level2( fntk, n, cand.l1_score );
      cand.l3_score = gain_level3( fntk, n, cand.l2_score );
      ranked.push_back( cand );
    } );

    st.num_scored += static_cast<uint32_t>( ranked.size() );
    sort_candidates( ranked );
    return ranked;
  }

  /* Stable-sort a candidate list using the configured mode comparator. */
  void sort_candidates( std::vector<mmig_seed_candidate<Ntk>>& v ) const
  {
    std::stable_sort( v.begin(), v.end(), [&]( auto const& a, auto const& b ) {
      return cmp_candidates( a, b );
    } );
  }

  bool cmp_candidates( mmig_seed_candidate<Ntk> const& a, mmig_seed_candidate<Ntk> const& b ) const
  {
    if ( ps.verify_each_candidate && ps.prefer_l1_when_verifying )
    {
      if ( a.l1_score != b.l1_score )
        return a.l1_score > b.l1_score;
      if ( a.l2_score != b.l2_score )
        return a.l2_score > b.l2_score;
      return a.n < b.n;
    }

    auto const pa = primary_key( a );
    auto const pb = primary_key( b );
    if ( pa != pb )
      return pa > pb;
    auto const sa = secondary_key( a );
    auto const sb = secondary_key( b );
    if ( sa != sb )
      return sa > sb;
    return a.n < b.n;
  }

  int32_t primary_key( mmig_seed_candidate<Ntk> const& c ) const
  {
    switch ( ps.mode )
    {
    case mmig_seed_mode::level1:
      return c.l1_score;
    case mmig_seed_mode::level2:
      return c.l2_score;
    case mmig_seed_mode::level3:
      return c.l3_score;
    case mmig_seed_mode::both:
    default:
      /* Keep original l1/l2 behavior — l3 is opt-in via level3 mode */
      return std::max( c.l1_score, c.l2_score );
    }
  }

  int32_t secondary_key( mmig_seed_candidate<Ntk> const& c ) const
  {
    switch ( ps.mode )
    {
    case mmig_seed_mode::level1:
      return c.l2_score;
    case mmig_seed_mode::level2:
      return c.l1_score;
    case mmig_seed_mode::level3:
      return c.l2_score;
    case mmig_seed_mode::both:
    default:
      return c.l1_score + c.l2_score;
    }
  }

  int32_t gain_level1( fanout_view<Ntk> const& fntk, node<Ntk> const& n ) const
  {
    int32_t gain = 0;

    fntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( fntk.is_constant( fntk.get_node( f ) ) )
      {
        return;
      }
      gain += fntk.is_complemented( f ) ? 1 : -1;
    } );

    fntk.foreach_fanout( n, [&]( auto const& parent ) {
      gain += is_complemented_parent( fntk, parent, n ) ? 1 : -1;
    } );

    fntk.foreach_po( [&]( auto const& s ) {
      if ( fntk.get_node( s ) == n )
      {
        gain += fntk.is_complemented( s ) ? 1 : -1;
      }
    } );

    return gain;
  }

  int32_t gain_level2( fanout_view<Ntk> const& fntk, node<Ntk> const& n, int32_t l1_gain ) const
  {
    int32_t gain = l1_gain;
    fntk.foreach_fanout( n, [&]( auto const& parent ) {
      int32_t sub = gain_level1( fntk, parent );
      if ( is_complemented_parent( fntk, parent, n ) )
      {
        sub -= 2;
      }
      else
      {
        sub += 2;
      }
      if ( sub > 0 )
      {
        gain += sub;
      }
    } );
    return gain;
  }

  /* L3 gain: for every grandparent of n (parent's fanout), compute L1 gain at
   * grandparent adjusted by the polarity of the parent→grandparent edge.  The
   * ±2 adjustment mirrors gain_level2: flipping n propagates a polarity change
   * through the parent, which changes whether the grandparent sees a complement
   * or non-complement on that incoming edge.  Only positive contributions are
   * accumulated (same conservative filter as L2). */
  int32_t gain_level3( fanout_view<Ntk> const& fntk, node<Ntk> const& n, int32_t l2_gain ) const
  {
    int32_t gain = l2_gain;
    fntk.foreach_fanout( n, [&]( auto const& parent ) {
      fntk.foreach_fanout( parent, [&]( auto const& grandparent ) {
        int32_t sub = gain_level1( fntk, grandparent );
        if ( is_complemented_parent( fntk, grandparent, parent ) )
          sub -= 2;
        else
          sub += 2;
        if ( sub > 0 )
          gain += sub;
      } );
    } );
    return gain;
  }

  bool is_complemented_parent( fanout_view<Ntk> const& fntk, node<Ntk> const& parent, node<Ntk> const& child ) const
  {
    bool found = false;
    bool complement = false;
    fntk.foreach_fanin( parent, [&]( auto const& f ) {
      if ( fntk.get_node( f ) == child )
      {
        found = true;
        complement = fntk.is_complemented( f );
      }
    } );
    return found && complement;
  }

  std::optional<mmig_seed_eval<Ntk>> evaluate_node( node<Ntk> const& n, mmig_seed_cost<Ntk> const& before )
  {
    auto eval_output = evaluate_form( n, mmig_seed_form::output_complement, before );
    std::optional<mmig_seed_eval<Ntk>> eval_input = std::nullopt;
    if ( ps.enable_input_complement_form )
    {
      eval_input = evaluate_form( n, mmig_seed_form::input_complement, before );
    }

    if ( !eval_output.has_value() && !eval_input.has_value() )
    {
      return std::nullopt;
    }
    if ( !eval_output.has_value() )
    {
      return eval_input;
    }
    if ( !eval_input.has_value() )
    {
      return eval_output;
    }
    return better_eval( *eval_input, *eval_output ) ? eval_input : eval_output;
  }

  std::optional<mmig_seed_eval<Ntk>> evaluate_form( node<Ntk> const& n, mmig_seed_form form, mmig_seed_cost<Ntk> const& before )
  {
    ++st.num_evaluated;
    auto trial = ntk.clone();
    if ( trial.is_dead( n ) || !trial.is_maj( n ) )
    {
      return std::nullopt;
    }

    auto replacement = make_replacement( trial, n, form );
    trial.substitute_node( n, replacement );
    trial = cleanup_dangling( trial );

    auto const after = network_cost( trial );

    mmig_seed_eval<Ntk> eval{};
    eval.n = n;
    eval.form = form;
    eval.after = after;
    eval.depth_gain = before.depth - after.depth;
    eval.area_gain = before.gates - after.gates;
    eval.inv_gain = before.inverted_edges - after.inverted_edges;
    return eval;
  }

  signal<Ntk> make_replacement( Ntk& net, node<Ntk> const& n, mmig_seed_form form ) const
  {
    std::array<signal<Ntk>, 3u> fanins{};
    net.foreach_fanin( n, [&]( auto const& f, auto i ) { fanins[i] = f; } );

    if ( form == mmig_seed_form::output_complement )
    {
      return !net.create_min( fanins[0u], fanins[1u], fanins[2u] );
    }
    return net.create_min( !fanins[0u], !fanins[1u], !fanins[2u] );
  }

  bool better_eval( mmig_seed_eval<Ntk> const& lhs, mmig_seed_eval<Ntk> const& rhs ) const
  {
    if ( lhs.depth_gain != rhs.depth_gain )
    {
      return lhs.depth_gain > rhs.depth_gain;
    }
    if ( lhs.area_gain != rhs.area_gain )
    {
      return lhs.area_gain > rhs.area_gain;
    }
    if ( lhs.inv_gain != rhs.inv_gain )
    {
      return lhs.inv_gain > rhs.inv_gain;
    }
    if ( lhs.form != rhs.form )
    {
      return lhs.form == mmig_seed_form::output_complement;
    }
    return lhs.n < rhs.n;
  }

  bool accept( mmig_seed_eval<Ntk> const& eval ) const
  {
    if ( ps.allow_non_improving )
    {
      return true;
    }

    if ( eval.depth_gain < 0 )
    {
      return false;
    }
    if ( !ps.allow_area_increase && eval.depth_gain > 0 && eval.area_gain < 0 )
    {
      return false;
    }
    if ( eval.depth_gain == 0 && eval.area_gain < 0 )
    {
      return false;
    }
    if ( eval.depth_gain == 0 && eval.area_gain == 0 && eval.inv_gain < 0 )
    {
      return false;
    }
    return true;
  }

  void apply_eval( mmig_seed_eval<Ntk> const& eval )
  {
    if ( ntk.is_dead( eval.n ) || !ntk.is_maj( eval.n ) )
    {
      return;
    }
    auto replacement = make_replacement( ntk, eval.n, eval.form );
    ntk.substitute_node( eval.n, replacement );
    if ( ps.verbose )
    {
      std::cout << "[mmig-seed] node=" << eval.n
                << " form=" << ( eval.form == mmig_seed_form::output_complement ? "out-comp" : "in-comp" )
                << " d=" << eval.depth_gain
                << " a=" << eval.area_gain
                << " i=" << eval.inv_gain << "\n";
    }
  }

  bool verify_candidate( node<Ntk> const& n, mmig_seed_form form )
  {
    ++st.num_equivalence_checks;
    auto trial = ntk.clone();
    if ( trial.is_dead( n ) || !trial.is_maj( n ) )
    {
      return false;
    }
    auto replacement = make_replacement( trial, n, form );
    trial.substitute_node( n, replacement );

    auto const maybe_miter = miter<Ntk>( ntk, trial );
    if ( !maybe_miter.has_value() )
    {
      return false;
    }
    auto const eq = equivalence_checking( *maybe_miter );
    return eq.has_value() && *eq;
  }

  mmig_seed_cost<Ntk> network_cost( Ntk const& net ) const
  {
    depth_view depth_net{ net };
    mmig_seed_cost<Ntk> cost{};
    cost.depth = static_cast<int32_t>( depth_net.depth() );
    cost.gates = static_cast<int32_t>( net.num_gates() );
    cost.inverted_edges = count_inverted_edges( net );
    return cost;
  }

  int32_t count_inverted_edges( Ntk const& net ) const
  {
    int32_t count = 0;
    net.foreach_gate( [&]( auto const& n ) {
      net.foreach_fanin( n, [&]( auto const& f ) {
        if ( net.is_complemented( f ) )
        {
          ++count;
        }
      } );
    } );
    net.foreach_po( [&]( auto const& s ) {
      if ( net.is_complemented( s ) )
      {
        ++count;
      }
    } );
    return count;
  }

private:
  Ntk& ntk;
  mmig_minority_seeding_params const& ps;
  mmig_minority_seeding_stats& st;
};

} // namespace detail

template<class Ntk>
void mmig_minority_seeding( Ntk& ntk,
                            mmig_minority_seeding_params const& ps = {},
                            mmig_minority_seeding_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( std::is_same_v<typename Ntk::base_type, mig_network>, "Ntk is not MIG-based" );
  static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement foreach_gate" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement foreach_fanin" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement foreach_po" );
  static_assert( has_is_maj_v<Ntk>, "Ntk does not implement is_maj" );
  static_assert( has_is_min_v<Ntk>, "Ntk does not implement is_min" );
  static_assert( has_create_min_v<Ntk>, "Ntk does not implement create_min" );
  static_assert( has_substitute_node_v<Ntk>, "Ntk does not implement substitute_node" );
  static_assert( has_is_dead_v<Ntk>, "Ntk does not implement is_dead" );
  static_assert( has_num_gates_v<Ntk>, "Ntk does not implement num_gates" );
  static_assert( has_clone_v<Ntk>, "Ntk does not implement clone" );

  mmig_minority_seeding_stats st{};
  detail::mmig_minority_seeding_impl<Ntk> impl( ntk, ps, st );
  impl.run();
  if ( pst != nullptr )
  {
    *pst = st;
  }
}

} // namespace mockturtle
