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
  \file mmig_algebraic_rewriting.hpp
  \brief mMIG algebraic rewriting

  \author Mockturtle-mMIG contributors
*/

#pragma once

#include "cleanup.hpp"

#include "../networks/mig.hpp"
#include "../traits.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/depth_view.hpp"
#include "../views/topo_view.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <optional>
#include <type_traits>
#include <vector>

namespace mockturtle
{

enum class mmig_rule_id : uint8_t
{
  associativity_1 = 0,
  associativity_2,
  distributivity_1,
  distributivity_2,
  sr1,
  sr2,
  sr3,
  sr4,
  sr5,
  sr6,
  snr1,
  snr2,
  snr3,
  snr4,
  snr5,
  snr6,
  relevance_1,
  relevance_2,
  num_rules
};

struct mmig_algebraic_rewriting_params
{
  uint32_t max_iterations{ 16u };
  uint32_t max_candidates_per_iteration{ 256u };
  uint32_t max_stagnation{ 3u };
  bool allow_area_increase{ false };
  bool allow_non_improving{ false };
  bool enable_sr5{ false };
  bool rank_candidates_globally{ true };
  bool normalize_complemented_inner{ true };
  bool use_mffc_cost{ false };
  bool verbose{ false };
};

struct mmig_algebraic_rewriting_stats
{
  stopwatch<>::duration time_total{ 0 };
  uint32_t num_iterations{ 0 };
  uint32_t num_candidates{ 0 };
  uint32_t num_evaluated{ 0 };
  uint32_t num_applied{ 0 };
  uint32_t num_rejected{ 0 };
  uint32_t num_depth_improvements{ 0 };
  uint32_t num_area_improvements{ 0 };
  uint32_t num_inverter_improvements{ 0 };
  std::array<uint32_t, static_cast<uint32_t>( mmig_rule_id::num_rules )> per_rule_applied{};
  std::array<uint32_t, static_cast<uint32_t>( mmig_rule_id::num_rules )> per_rule_rejected{};
};

namespace detail
{

enum class mmig_gate_kind : uint8_t
{
  maj,
  min
};

enum class mmig_rule_class : uint8_t
{
  none,
  sr,
  snr
};

enum class mmig_var_id : uint8_t
{
  x = 0u,
  y,
  z,
  u,
  v,
  w,
  num_vars
};

struct mmig_var_pattern
{
  mmig_var_id id;
  bool complemented{ false };
};

struct mmig_rule_pattern
{
  mmig_gate_kind outer_kind;
  mmig_gate_kind inner_kind;
  std::array<mmig_var_pattern, 2u> outer_vars;
  std::array<mmig_var_pattern, 3u> inner_vars;
};

template<class Ntk>
struct mmig_bindings
{
  std::array<signal<Ntk>, static_cast<uint32_t>( mmig_var_id::num_vars )> values{};
  std::array<uint8_t, static_cast<uint32_t>( mmig_var_id::num_vars )> assigned{};

  void clear()
  {
    assigned.fill( 0u );
  }

  bool has( mmig_var_id id ) const
  {
    return assigned[static_cast<uint32_t>( id )] != 0u;
  }

  signal<Ntk> variable( mmig_var_id id ) const
  {
    return values[static_cast<uint32_t>( id )];
  }

  signal<Ntk> literal( mmig_var_pattern const& pat ) const
  {
    auto const base = values[static_cast<uint32_t>( pat.id )];
    return pat.complemented ? !base : base;
  }

  bool bind( mmig_var_pattern const& pat, signal<Ntk> const& actual )
  {
    auto const idx = static_cast<uint32_t>( pat.id );
    auto const expected_base = pat.complemented ? !actual : actual;
    if ( assigned[idx] == 0u )
    {
      values[idx] = expected_base;
      assigned[idx] = 1u;
      return true;
    }

    auto const expected = pat.complemented ? !values[idx] : values[idx];
    return expected == actual;
  }
};

template<class Ntk>
struct mmig_rewrite_candidate
{
  node<Ntk> root{};
  signal<Ntk> replacement{};
  mmig_rule_id rule{ mmig_rule_id::associativity_1 };
  mmig_rule_class rule_class{ mmig_rule_class::none };
  mmig_bindings<Ntk> bindings{};
};

template<class Ntk>
struct mmig_scored_candidate
{
  mmig_rewrite_candidate<Ntk> candidate;
  int32_t depth_gain{ std::numeric_limits<int32_t>::min() };
  int32_t area_gain{ std::numeric_limits<int32_t>::min() };
  int32_t inv_gain{ std::numeric_limits<int32_t>::min() };
  int32_t sr_snr_preference{ 0 };
};

template<class Ntk>
struct mmig_cheap_candidate
{
  mmig_rewrite_candidate<Ntk> candidate;
  int32_t level_gain{ std::numeric_limits<int32_t>::min() };
  int32_t estimated_level{ std::numeric_limits<int32_t>::max() };
  int32_t new_nodes{ std::numeric_limits<int32_t>::max() };
  int32_t mffc_size{ 0 };
  int32_t mffc_gain{ std::numeric_limits<int32_t>::min() };
  int32_t inv_burden{ std::numeric_limits<int32_t>::max() };
  int32_t fanout_penalty{ std::numeric_limits<int32_t>::max() };
  int32_t hash_bonus{ std::numeric_limits<int32_t>::min() };
  int32_t sr_snr_preference{ 0 };
};

struct mmig_network_cost
{
  int32_t depth{ 0 };
  int32_t gates{ 0 };
  int32_t inverted_edges{ 0 };
};

constexpr std::array<std::array<uint32_t, 2u>, 2u> k_permutations2{ {
    { { 0u, 1u } },
    { { 1u, 0u } },
} };

constexpr std::array<std::array<uint32_t, 3u>, 6u> k_permutations3{ {
    { { 0u, 1u, 2u } },
    { { 0u, 2u, 1u } },
    { { 1u, 0u, 2u } },
    { { 1u, 2u, 0u } },
    { { 2u, 0u, 1u } },
    { { 2u, 1u, 0u } },
} };

constexpr mmig_var_pattern x( bool complemented = false )
{
  return { mmig_var_id::x, complemented };
}

constexpr mmig_var_pattern y( bool complemented = false )
{
  return { mmig_var_id::y, complemented };
}

constexpr mmig_var_pattern z( bool complemented = false )
{
  return { mmig_var_id::z, complemented };
}

constexpr mmig_var_pattern u( bool complemented = false )
{
  return { mmig_var_id::u, complemented };
}

constexpr mmig_var_pattern v( bool complemented = false )
{
  return { mmig_var_id::v, complemented };
}

inline char const* mmig_rule_name( mmig_rule_id id )
{
  switch ( id )
  {
  case mmig_rule_id::associativity_1:
    return "A1";
  case mmig_rule_id::associativity_2:
    return "A2";
  case mmig_rule_id::distributivity_1:
    return "D1";
  case mmig_rule_id::distributivity_2:
    return "D2";
  case mmig_rule_id::sr1:
    return "SR1";
  case mmig_rule_id::sr2:
    return "SR2";
  case mmig_rule_id::sr3:
    return "SR3";
  case mmig_rule_id::sr4:
    return "SR4";
  case mmig_rule_id::sr5:
    return "SR5";
  case mmig_rule_id::sr6:
    return "SR6";
  case mmig_rule_id::snr1:
    return "SNR1";
  case mmig_rule_id::snr2:
    return "SNR2";
  case mmig_rule_id::snr3:
    return "SNR3";
  case mmig_rule_id::snr4:
    return "SNR4";
  case mmig_rule_id::snr5:
    return "SNR5";
  case mmig_rule_id::snr6:
    return "SNR6";
  case mmig_rule_id::relevance_1:
    return "R1";
  case mmig_rule_id::relevance_2:
    return "R2";
  default:
    return "unknown";
  }
}

template<class Ntk>
class mmig_algebraic_rewriting_impl
{
public:
  mmig_algebraic_rewriting_impl( Ntk& ntk, mmig_algebraic_rewriting_params const& ps, mmig_algebraic_rewriting_stats& st )
      : ntk( ntk ),
        ps( ps ),
        st( st )
  {
  }

  void run()
  {
    stopwatch t( st.time_total );
    uint32_t stagnation_counter = 0u;

    for ( uint32_t iteration = 0u; iteration < ps.max_iterations; ++iteration )
    {
      ++st.num_iterations;

      auto const before = network_cost( ntk );
      auto best = find_best_candidate( before );
      if ( !best.has_value() )
      {
        break;
      }

      ntk.substitute_node( best->candidate.root, best->candidate.replacement );
      ntk = cleanup_dangling( ntk );

      auto const after = network_cost( ntk );
      ++st.num_applied;
      ++st.per_rule_applied[static_cast<uint32_t>( best->candidate.rule )];
      if ( best->depth_gain > 0 )
      {
        ++st.num_depth_improvements;
      }
      if ( best->area_gain > 0 )
      {
        ++st.num_area_improvements;
      }
      if ( best->inv_gain > 0 )
      {
        ++st.num_inverter_improvements;
      }

      if ( ps.verbose )
      {
        std::cout << "[mmig-rw] iter=" << ( iteration + 1u )
                  << " rule=" << mmig_rule_name( best->candidate.rule )
                  << " depth_gain=" << best->depth_gain
                  << " area_gain=" << best->area_gain
                  << " inv_gain=" << best->inv_gain << "\n";
      }

      auto const progressed = ( after.depth < before.depth ) ||
                              ( after.depth == before.depth && after.gates < before.gates ) ||
                              ( after.depth == before.depth && after.gates == before.gates && after.inverted_edges < before.inverted_edges );

      if ( progressed )
      {
        stagnation_counter = 0u;
      }
      else
      {
        ++stagnation_counter;
        if ( ps.max_stagnation > 0u && stagnation_counter >= ps.max_stagnation )
        {
          break;
        }
      }
    }
  }

private:
  template<class Fn>
  std::optional<mmig_rewrite_candidate<Ntk>> match_rule( node<Ntk> const& n,
                                                         mmig_rule_id rule,
                                                         mmig_rule_class rule_class,
                                                         mmig_rule_pattern const& lhs,
                                                         Fn&& build_replacement ) const
  {
    if ( !is_kind( n, lhs.outer_kind ) )
    {
      return std::nullopt;
    }

    auto const original_size = ntk.size();
    auto const outer_fanins = fanins_of( n );
    std::optional<mmig_cheap_candidate<Ntk>> best;
    for ( uint32_t inner_pos = 0u; inner_pos < 3u; ++inner_pos )
    {
      auto const inner_signal = outer_fanins[inner_pos];
      auto const inner_node = ntk.get_node( inner_signal );
      if ( ntk.is_constant( inner_node ) || ntk.is_pi( inner_node ) )
      {
        continue;
      }

      auto inner_kind = ntk.is_maj( inner_node ) ? mmig_gate_kind::maj : mmig_gate_kind::min;
      auto inner_fanins = fanins_of( inner_node );
      if ( ntk.is_complemented( inner_signal ) )
      {
        if ( !ps.normalize_complemented_inner )
        {
          continue;
        }
        inner_kind = ( inner_kind == mmig_gate_kind::maj ) ? mmig_gate_kind::min : mmig_gate_kind::maj;
        for ( auto& f : inner_fanins )
        {
          f = !f;
        }
      }

      if ( inner_kind != lhs.inner_kind )
      {
        continue;
      }

      std::array<uint32_t, 2u> outer_positions{};
      uint32_t p = 0u;
      for ( uint32_t i = 0u; i < 3u; ++i )
      {
        if ( i != inner_pos )
        {
          outer_positions[p++] = i;
        }
      }

      for ( auto const& outer_perm : k_permutations2 )
      {
        mmig_bindings<Ntk> prefix_bindings{};
        prefix_bindings.clear();
        if ( !prefix_bindings.bind( lhs.outer_vars[0], outer_fanins[outer_positions[outer_perm[0]]] ) )
        {
          continue;
        }
        if ( !prefix_bindings.bind( lhs.outer_vars[1], outer_fanins[outer_positions[outer_perm[1]]] ) )
        {
          continue;
        }

        for ( auto const& inner_perm : k_permutations3 )
        {
          auto bindings = prefix_bindings;
          if ( !bindings.bind( lhs.inner_vars[0], inner_fanins[inner_perm[0]] ) )
          {
            continue;
          }
          if ( !bindings.bind( lhs.inner_vars[1], inner_fanins[inner_perm[1]] ) )
          {
            continue;
          }
          if ( !bindings.bind( lhs.inner_vars[2], inner_fanins[inner_perm[2]] ) )
          {
            continue;
          }

          auto replacement = std::invoke( build_replacement, bindings );
          if ( replacement == ntk.make_signal( n ) )
          {
            continue;
          }

          mmig_rewrite_candidate<Ntk> cand{};
          cand.root = n;
          cand.replacement = replacement;
          cand.rule = rule;
          cand.rule_class = rule_class;
          cand.bindings = bindings;

          auto scored = score_candidate_cheap( cand, original_size );
          if ( !best.has_value() || better_cheap_candidate( scored, *best ) )
          {
            best = scored;
          }
        }
      }
    }

    if ( best.has_value() )
    {
      return best->candidate;
    }
    return std::nullopt;
  }

  std::vector<mmig_rewrite_candidate<Ntk>> collect_candidates() const
  {
    std::vector<mmig_cheap_candidate<Ntk>> ranked;
    ranked.reserve( ps.max_candidates_per_iteration );

    topo_view topo{ ntk };
    topo.foreach_node( [&]( auto const& n ) {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      {
        return true;
      }
      if ( !ntk.is_maj( n ) && !ntk.is_min( n ) )
      {
        return true;
      }

      auto add = [&]( auto&& maybe_candidate ) {
        if ( !maybe_candidate.has_value() )
        {
          return;
        }

        if ( !ps.rank_candidates_globally )
        {
          if ( ranked.size() < ps.max_candidates_per_iteration )
          {
            ranked.push_back( score_candidate_cheap( *maybe_candidate ) );
          }
          return;
        }

        auto scored = score_candidate_cheap( *maybe_candidate );
        if ( ranked.size() < ps.max_candidates_per_iteration )
        {
          ranked.push_back( scored );
          return;
        }

        auto worst = std::min_element( ranked.begin(), ranked.end(), [&]( auto const& a, auto const& b ) {
          return better_cheap_candidate( b, a );
        } );
        if ( worst != ranked.end() && better_cheap_candidate( scored, *worst ) )
        {
          *worst = scored;
        }
      };

      add( try_associativity_1( n ) );
      add( try_associativity_2( n ) );
      add( try_distributivity_1( n ) );
      add( try_distributivity_2( n ) );
      add( try_sr1( n ) );
      add( try_sr2( n ) );
      add( try_sr3( n ) );
      add( try_sr4( n ) );
      if ( ps.enable_sr5 )
      {
        add( try_sr5( n ) );
      }
      add( try_sr6( n ) );
      add( try_snr1( n ) );
      add( try_snr2( n ) );
      add( try_snr3( n ) );
      add( try_snr4( n ) );
      add( try_snr5( n ) );
      add( try_snr6( n ) );
      add( try_relevance_1( n ) );
      add( try_relevance_2( n ) );

      return !ps.rank_candidates_globally || ranked.size() < ps.max_candidates_per_iteration;
    } );

    std::sort( ranked.begin(), ranked.end(), [&]( auto const& a, auto const& b ) {
      return better_cheap_candidate( a, b );
    } );

    std::vector<mmig_rewrite_candidate<Ntk>> cands;
    cands.reserve( ranked.size() );
    for ( auto const& cand : ranked )
    {
      cands.push_back( cand.candidate );
    }
    return cands;
  }

  std::optional<mmig_scored_candidate<Ntk>> find_best_candidate( mmig_network_cost const& before ) const
  {
    auto candidates = collect_candidates();
    st.num_candidates += static_cast<uint32_t>( candidates.size() );
    if ( candidates.empty() )
    {
      return std::nullopt;
    }

    std::optional<mmig_scored_candidate<Ntk>> best;
    for ( auto const& cand : candidates )
    {
      ++st.num_evaluated;
      auto scored = evaluate_candidate( cand, before );
      if ( !scored.has_value() )
      {
        ++st.num_rejected;
        ++st.per_rule_rejected[static_cast<uint32_t>( cand.rule )];
        continue;
      }

      if ( !best.has_value() || better_candidate( *scored, *best ) )
      {
        best = scored;
      }
    }

    return best;
  }

  std::optional<mmig_scored_candidate<Ntk>> evaluate_candidate( mmig_rewrite_candidate<Ntk> const& cand,
                                                                mmig_network_cost const& before ) const
  {
    auto candidate_ntk = ntk.clone();
    candidate_ntk.substitute_node( cand.root, cand.replacement );
    candidate_ntk = cleanup_dangling( candidate_ntk );

    auto const after = network_cost( candidate_ntk );

    mmig_scored_candidate<Ntk> scored{};
    scored.candidate = cand;
    scored.depth_gain = before.depth - after.depth;
    scored.area_gain = before.gates - after.gates;
    scored.inv_gain = before.inverted_edges - after.inverted_edges;
    scored.sr_snr_preference = sr_snr_preference( cand );

    if ( !ps.allow_non_improving )
    {
      if ( scored.depth_gain < 0 )
      {
        return std::nullopt;
      }
      if ( scored.depth_gain > 0 && !ps.allow_area_increase && scored.area_gain < 0 )
      {
        return std::nullopt;
      }
      if ( scored.depth_gain == 0 && scored.area_gain < 0 )
      {
        return std::nullopt;
      }
      if ( scored.depth_gain == 0 && scored.area_gain == 0 && scored.inv_gain <= 0 )
      {
        return std::nullopt;
      }
    }

    return scored;
  }

  bool better_candidate( mmig_scored_candidate<Ntk> const& lhs, mmig_scored_candidate<Ntk> const& rhs ) const
  {
    if ( ps.allow_area_increase )
    {
      if ( lhs.depth_gain != rhs.depth_gain )
      {
        return lhs.depth_gain > rhs.depth_gain;
      }
      if ( lhs.area_gain != rhs.area_gain )
      {
        return lhs.area_gain > rhs.area_gain;
      }
    }
    else
    {
      if ( lhs.area_gain != rhs.area_gain )
      {
        return lhs.area_gain > rhs.area_gain;
      }
      if ( lhs.depth_gain != rhs.depth_gain )
      {
        return lhs.depth_gain > rhs.depth_gain;
      }
    }
    if ( lhs.inv_gain != rhs.inv_gain )
    {
      return lhs.inv_gain > rhs.inv_gain;
    }
    if ( lhs.sr_snr_preference != rhs.sr_snr_preference )
    {
      return lhs.sr_snr_preference > rhs.sr_snr_preference;
    }
    if ( lhs.candidate.root != rhs.candidate.root )
    {
      return lhs.candidate.root < rhs.candidate.root;
    }
    if ( lhs.candidate.rule != rhs.candidate.rule )
    {
      return static_cast<uint32_t>( lhs.candidate.rule ) < static_cast<uint32_t>( rhs.candidate.rule );
    }
    return lhs.candidate.replacement < rhs.candidate.replacement;
  }

  int32_t sr_snr_preference( mmig_rewrite_candidate<Ntk> const& cand ) const
  {
    if ( cand.rule_class == mmig_rule_class::none )
    {
      return 0;
    }
    if ( !cand.bindings.has( mmig_var_id::x ) || !cand.bindings.has( mmig_var_id::y ) ||
         !cand.bindings.has( mmig_var_id::z ) || !cand.bindings.has( mmig_var_id::u ) )
    {
      return 0;
    }

    auto const tx = signal_level( cand.bindings.variable( mmig_var_id::x ) );
    auto const ty = signal_level( cand.bindings.variable( mmig_var_id::y ) );
    auto const tz = signal_level( cand.bindings.variable( mmig_var_id::z ) );
    auto const tu = signal_level( cand.bindings.variable( mmig_var_id::u ) );

    if ( tu > std::max( { tx, ty, tz } ) )
    {
      return cand.rule_class == mmig_rule_class::sr ? 1 : -1;
    }
    if ( tz > std::max( { tx, ty, tu } ) )
    {
      return cand.rule_class == mmig_rule_class::snr ? 1 : -1;
    }
    return 0;
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_associativity_1( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::min, mmig_gate_kind::maj,
        { { x(), z() } },
        { { x(), u(), y() } } };
    return match_rule( n, mmig_rule_id::associativity_1, mmig_rule_class::none, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::min, b.literal( x() ), b.literal( u() ),
                          create_gate( mmig_gate_kind::maj, b.literal( x() ), b.literal( y() ), b.literal( z() ) ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_associativity_2( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::maj, mmig_gate_kind::min,
        { { x(), z() } },
        { { u(), x(), y() } } };
    return match_rule( n, mmig_rule_id::associativity_2, mmig_rule_class::none, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::maj, b.literal( u( true ) ), b.literal( z() ),
                          create_gate( mmig_gate_kind::min, b.literal( x( true ) ), b.literal( y() ), b.literal( z( true ) ) ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_distributivity_1( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::min, mmig_gate_kind::maj,
        { { x(), y() } },
        { { u(), v(), z() } } };
    return match_rule( n, mmig_rule_id::distributivity_1, mmig_rule_class::none, lhs, [this]( auto const& b ) {
      auto const left = create_gate( mmig_gate_kind::maj, b.literal( x() ), b.literal( y() ), b.literal( u() ) );
      auto const right = create_gate( mmig_gate_kind::maj, b.literal( x() ), b.literal( y() ), b.literal( v() ) );
      return create_gate( mmig_gate_kind::min, left, right, b.literal( z() ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_distributivity_2( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::maj, mmig_gate_kind::min,
        { { x(), y() } },
        { { u(), v(), z() } } };
    return match_rule( n, mmig_rule_id::distributivity_2, mmig_rule_class::none, lhs, [this]( auto const& b ) {
      auto const left = create_gate( mmig_gate_kind::min, b.literal( x( true ) ), b.literal( y( true ) ), b.literal( u() ) );
      auto const right = create_gate( mmig_gate_kind::min, b.literal( x( true ) ), b.literal( y( true ) ), b.literal( v() ) );
      return create_gate( mmig_gate_kind::maj, left, right, b.literal( z( true ) ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_sr1( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::min, mmig_gate_kind::maj,
        { { x(), u() } },
        { { u( true ), y(), z() } } };
    return match_rule( n, mmig_rule_id::sr1, mmig_rule_class::sr, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::min, b.literal( x() ), b.literal( z() ),
                          create_gate( mmig_gate_kind::maj, b.literal( x() ), b.literal( y() ), b.literal( z() ) ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_sr2( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::min, mmig_gate_kind::maj,
        { { x(), u() } },
        { { u(), y(), z( true ) } } };
    return match_rule( n, mmig_rule_id::sr2, mmig_rule_class::sr, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::min, b.literal( x() ), b.literal( z() ),
                          create_gate( mmig_gate_kind::maj, b.literal( x() ), b.literal( u() ), b.literal( y( true ) ) ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_sr3( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::min, mmig_gate_kind::maj,
        { { x(), u( true ) } },
        { { u(), y( true ), z() } } };
    return match_rule( n, mmig_rule_id::sr3, mmig_rule_class::sr, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::min, b.literal( x( true ) ), b.literal( z() ),
                          create_gate( mmig_gate_kind::maj, b.literal( x( true ) ), b.literal( u( true ) ), b.literal( y( true ) ) ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_sr4( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::maj, mmig_gate_kind::min,
        { { x(), u() } },
        { { u(), y(), z( true ) } } };
    return match_rule( n, mmig_rule_id::sr4, mmig_rule_class::sr, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::maj, b.literal( x() ), b.literal( z() ),
                          create_gate( mmig_gate_kind::min, b.literal( u( true ) ), b.literal( y() ), b.literal( x( true ) ) ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_sr5( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::maj, mmig_gate_kind::min,
        { { x( true ), u( true ) } },
        { { u( true ), y(), z() } } };
    return match_rule( n, mmig_rule_id::sr5, mmig_rule_class::sr, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::maj, b.literal( x( true ) ), b.literal( u( true ) ),
                          create_gate( mmig_gate_kind::min, b.literal( x() ), b.literal( y() ), b.literal( z() ) ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_sr6( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::maj, mmig_gate_kind::min,
        { { x(), u( true ) } },
        { { u( true ), y(), z() } } };
    return match_rule( n, mmig_rule_id::sr6, mmig_rule_class::sr, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::maj, b.literal( x() ), b.literal( z( true ) ),
                          create_gate( mmig_gate_kind::min, b.literal( x( true ) ), b.literal( u() ), b.literal( y() ) ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_snr1( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::min, mmig_gate_kind::maj,
        { { x(), u() } },
        { { u(), y(), z( true ) } } };
    return match_rule( n, mmig_rule_id::snr1, mmig_rule_class::snr, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::min, b.literal( z( true ) ), b.literal( u() ),
                          create_gate( mmig_gate_kind::maj, b.literal( u() ), b.literal( x() ), b.literal( y() ) ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_snr2( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::min, mmig_gate_kind::maj,
        { { x(), u( true ) } },
        { { u( true ), y(), z() } } };
    return match_rule( n, mmig_rule_id::snr2, mmig_rule_class::snr, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::min, b.literal( z() ), b.literal( u( true ) ),
                          create_gate( mmig_gate_kind::maj, b.literal( u( true ) ), b.literal( x() ), b.literal( y() ) ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_snr3( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::min, mmig_gate_kind::maj,
        { { x( true ), u() } },
        { { u(), y(), z() } } };
    return match_rule( n, mmig_rule_id::snr3, mmig_rule_class::snr, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::min, b.literal( z() ), b.literal( u() ),
                          create_gate( mmig_gate_kind::maj, b.literal( u() ), b.literal( x( true ) ), b.literal( y() ) ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_snr4( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::maj, mmig_gate_kind::min,
        { { x( true ), u( true ) } },
        { { u(), y(), z() } } };
    return match_rule( n, mmig_rule_id::snr4, mmig_rule_class::snr, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::maj, b.literal( z( true ) ), b.literal( u( true ) ),
                          create_gate( mmig_gate_kind::min, b.literal( u() ), b.literal( x() ), b.literal( y() ) ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_snr5( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::maj, mmig_gate_kind::min,
        { { x( true ), u( true ) } },
        { { u(), y( true ), z() } } };
    return match_rule( n, mmig_rule_id::snr5, mmig_rule_class::snr, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::maj, b.literal( y() ), b.literal( u( true ) ),
                          create_gate( mmig_gate_kind::min, b.literal( x() ), b.literal( u() ), b.literal( z() ) ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_snr6( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::maj, mmig_gate_kind::min,
        { { x(), u() } },
        { { u( true ), y(), z() } } };
    return match_rule( n, mmig_rule_id::snr6, mmig_rule_class::snr, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::maj, b.literal( z( true ) ), b.literal( u() ),
                          create_gate( mmig_gate_kind::min, b.literal( u( true ) ), b.literal( y() ), b.literal( x( true ) ) ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_relevance_1( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::min, mmig_gate_kind::maj,
        { { x(), y() } },
        { { x(), y(), z() } } };
    return match_rule( n, mmig_rule_id::relevance_1, mmig_rule_class::none, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::min, b.literal( x() ), b.literal( y() ), b.literal( z() ) );
    } );
  }

  std::optional<mmig_rewrite_candidate<Ntk>> try_relevance_2( node<Ntk> const& n ) const
  {
    mmig_rule_pattern lhs{
        mmig_gate_kind::maj, mmig_gate_kind::min,
        { { x(), y() } },
        { { x(), y(), z() } } };
    return match_rule( n, mmig_rule_id::relevance_2, mmig_rule_class::none, lhs, [this]( auto const& b ) {
      return create_gate( mmig_gate_kind::maj, b.literal( x() ), b.literal( y() ), b.literal( z( true ) ) );
    } );
  }

  bool is_kind( node<Ntk> const& n, mmig_gate_kind kind ) const
  {
    return ( kind == mmig_gate_kind::maj ) ? ntk.is_maj( n ) : ntk.is_min( n );
  }

  std::array<signal<Ntk>, 3u> fanins_of( node<Ntk> const& n ) const
  {
    std::array<signal<Ntk>, 3u> fs{};
    ntk.foreach_fanin( n, [&]( auto const& f, auto i ) { fs[i] = f; } );
    return fs;
  }

  signal<Ntk> create_gate( mmig_gate_kind kind, signal<Ntk> const& a, signal<Ntk> const& b, signal<Ntk> const& c ) const
  {
    return ( kind == mmig_gate_kind::maj ) ? ntk.create_maj( a, b, c ) : ntk.create_min( a, b, c );
  }

  int32_t node_level( node<Ntk> const& n ) const
  {
    if constexpr ( has_level_v<Ntk> )
    {
      return static_cast<int32_t>( ntk.level( n ) );
    }
    else
    {
      depth_view depth_ntk{ ntk };
      return static_cast<int32_t>( depth_ntk.level( n ) );
    }
  }

  int32_t signal_level( signal<Ntk> const& s ) const
  {
    auto const n = ntk.get_node( s );
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
    {
      return 0;
    }
    return node_level( n );
  }

  int32_t candidate_new_nodes( signal<Ntk> const& s, uint32_t original_size ) const
  {
    auto const n = ntk.get_node( s );
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) || n < original_size )
    {
      return 0;
    }

    int32_t total = 1;
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      total += candidate_new_nodes( f, original_size );
    } );
    return total;
  }

  int32_t candidate_hash_bonus( signal<Ntk> const& s, uint32_t original_size ) const
  {
    auto const n = ntk.get_node( s );
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
    {
      return 0;
    }
    if ( n < original_size )
    {
      return 1;
    }

    int32_t total = 0;
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      total += candidate_hash_bonus( f, original_size );
    } );
    return total;
  }

  int32_t candidate_inv_burden( signal<Ntk> const& s, uint32_t original_size ) const
  {
    auto const n = ntk.get_node( s );
    int32_t total = ntk.is_complemented( s ) ? 1 : 0;
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) || n < original_size )
    {
      return total;
    }

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      total += candidate_inv_burden( f, original_size );
    } );
    return total;
  }

  int32_t candidate_estimated_level( signal<Ntk> const& s, uint32_t original_size ) const
  {
    auto const n = ntk.get_node( s );
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
    {
      return 0;
    }
    if ( n < original_size )
    {
      return node_level( n );
    }

    int32_t level = 0;
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      level = std::max( level, candidate_estimated_level( f, original_size ) );
    } );
    return level + 1;
  }

  int32_t candidate_fanout_penalty( mmig_bindings<Ntk> const& bindings ) const
  {
    std::array<node<Ntk>, static_cast<uint32_t>( mmig_var_id::num_vars )> seen{};
    uint32_t seen_size = 0u;
    int32_t total = 0;
    for ( uint32_t i = 0u; i < static_cast<uint32_t>( mmig_var_id::num_vars ); ++i )
    {
      if ( bindings.assigned[i] == 0u )
      {
        continue;
      }
      auto const n = ntk.get_node( bindings.values[i] );
      if ( ntk.is_constant( n ) )
      {
        continue;
      }
      if ( std::find( seen.begin(), seen.begin() + seen_size, n ) != seen.begin() + seen_size )
      {
        continue;
      }
      seen[seen_size++] = n;
      total += static_cast<int32_t>( ntk.fanout_size( n ) );
    }
    return total;
  }

  int32_t estimate_mffc_size_rec( node<Ntk> const& n, bool root, std::vector<node<Ntk>>& seen ) const
  {
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
    {
      return 0;
    }
    if ( std::find( seen.begin(), seen.end(), n ) != seen.end() )
    {
      return 0;
    }
    if ( !root && ntk.fanout_size( n ) != 1u )
    {
      return 0;
    }

    seen.push_back( n );
    int32_t total = 1;
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      total += estimate_mffc_size_rec( ntk.get_node( f ), false, seen );
    } );
    return total;
  }

  int32_t estimate_mffc_size( node<Ntk> const& n ) const
  {
    if ( !ps.use_mffc_cost )
    {
      return 1;
    }
    std::vector<node<Ntk>> seen;
    seen.reserve( 16u );
    return estimate_mffc_size_rec( n, true, seen );
  }

  mmig_cheap_candidate<Ntk> score_candidate_cheap( mmig_rewrite_candidate<Ntk> const& cand, uint32_t original_size ) const
  {
    mmig_cheap_candidate<Ntk> scored{};
    scored.candidate = cand;
    scored.estimated_level = candidate_estimated_level( cand.replacement, original_size );
    scored.level_gain = node_level( cand.root ) - scored.estimated_level;
    scored.new_nodes = candidate_new_nodes( cand.replacement, original_size );
    scored.mffc_size = estimate_mffc_size( cand.root );
    scored.mffc_gain = scored.mffc_size - scored.new_nodes;
    scored.inv_burden = candidate_inv_burden( cand.replacement, original_size );
    scored.fanout_penalty = candidate_fanout_penalty( cand.bindings );
    scored.hash_bonus = candidate_hash_bonus( cand.replacement, original_size );
    scored.sr_snr_preference = sr_snr_preference( cand );
    return scored;
  }

  mmig_cheap_candidate<Ntk> score_candidate_cheap( mmig_rewrite_candidate<Ntk> const& cand ) const
  {
    return score_candidate_cheap( cand, ntk.size() );
  }

  bool better_cheap_candidate( mmig_cheap_candidate<Ntk> const& lhs, mmig_cheap_candidate<Ntk> const& rhs ) const
  {
    if ( ps.allow_area_increase )
    {
      if ( lhs.level_gain != rhs.level_gain )
      {
        return lhs.level_gain > rhs.level_gain;
      }
      if ( ps.use_mffc_cost && lhs.mffc_gain != rhs.mffc_gain )
      {
        return lhs.mffc_gain > rhs.mffc_gain;
      }
      if ( lhs.hash_bonus != rhs.hash_bonus )
      {
        return lhs.hash_bonus > rhs.hash_bonus;
      }
    }
    else
    {
      if ( ps.use_mffc_cost && lhs.mffc_gain != rhs.mffc_gain )
      {
        return lhs.mffc_gain > rhs.mffc_gain;
      }
      if ( lhs.new_nodes != rhs.new_nodes )
      {
        return lhs.new_nodes < rhs.new_nodes;
      }
      if ( lhs.hash_bonus != rhs.hash_bonus )
      {
        return lhs.hash_bonus > rhs.hash_bonus;
      }
    }

    if ( lhs.inv_burden != rhs.inv_burden )
    {
      return lhs.inv_burden < rhs.inv_burden;
    }
    if ( lhs.sr_snr_preference != rhs.sr_snr_preference )
    {
      return lhs.sr_snr_preference > rhs.sr_snr_preference;
    }
    if ( lhs.estimated_level != rhs.estimated_level )
    {
      return lhs.estimated_level < rhs.estimated_level;
    }
    if ( lhs.fanout_penalty != rhs.fanout_penalty )
    {
      return lhs.fanout_penalty < rhs.fanout_penalty;
    }
    if ( lhs.candidate.root != rhs.candidate.root )
    {
      return lhs.candidate.root < rhs.candidate.root;
    }
    if ( lhs.candidate.rule != rhs.candidate.rule )
    {
      return static_cast<uint32_t>( lhs.candidate.rule ) < static_cast<uint32_t>( rhs.candidate.rule );
    }
    return lhs.candidate.replacement < rhs.candidate.replacement;
  }

  int32_t count_inverted_edges( Ntk const& net ) const
  {
    int32_t inv_edges = 0;
    net.foreach_gate( [&]( auto const& n ) {
      net.foreach_fanin( n, [&]( auto const& f ) {
        if ( net.is_complemented( f ) )
        {
          ++inv_edges;
        }
      } );
    } );
    net.foreach_po( [&]( auto const& s ) {
      if ( net.is_complemented( s ) )
      {
        ++inv_edges;
      }
    } );
    return inv_edges;
  }

  mmig_network_cost network_cost( Ntk const& net ) const
  {
    depth_view depth_net{ net };
    mmig_network_cost cost{};
    cost.depth = static_cast<int32_t>( depth_net.depth() );
    cost.gates = static_cast<int32_t>( net.num_gates() );
    cost.inverted_edges = count_inverted_edges( net );
    return cost;
  }

private:
  Ntk& ntk;
  mmig_algebraic_rewriting_params const& ps;
  mmig_algebraic_rewriting_stats& st;
};

} // namespace detail

template<class Ntk>
void mmig_algebraic_rewriting( Ntk& ntk,
                               mmig_algebraic_rewriting_params const& ps = {},
                               mmig_algebraic_rewriting_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( std::is_same_v<typename Ntk::base_type, mig_network>, "Ntk is not MIG-based" );
  static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement foreach_gate" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement foreach_fanin" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement foreach_po" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement fanout_size" );
  static_assert( has_is_maj_v<Ntk>, "Ntk does not implement is_maj" );
  static_assert( has_is_min_v<Ntk>, "Ntk does not implement is_min" );
  static_assert( has_create_maj_v<Ntk>, "Ntk does not implement create_maj" );
  static_assert( has_create_min_v<Ntk>, "Ntk does not implement create_min" );
  static_assert( has_substitute_node_v<Ntk>, "Ntk does not implement substitute_node" );
  static_assert( has_num_gates_v<Ntk>, "Ntk does not implement num_gates" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement size" );
  static_assert( has_make_signal_v<Ntk>, "Ntk does not implement make_signal" );
  static_assert( has_clone_v<Ntk>, "Ntk does not implement clone" );

  mmig_algebraic_rewriting_stats st{};
  detail::mmig_algebraic_rewriting_impl<Ntk> impl( ntk, ps, st );
  impl.run();
  if ( pst != nullptr )
  {
    *pst = st;
  }
}

} // namespace mockturtle
