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
  \file mmig_cone_polarity_flip.hpp
  \brief Cone-level MAJ<->MIN polarity flip pass.

  An mMIG-only optimisation: pick a connected cone of gates, flip every MAJ
  to a MIN (and vice versa) keeping the same fanin signals, and rebalance
  inverter polarity at the cone boundary. A pure MIG flow cannot perform
  this move because it has no MIN gate to flip to.

  Identity used:
      MIN(a, b, c) = !MAJ(a, b, c)
  Flipping a node's gate type therefore negates its output. To preserve
  function, every edge crossing OUT of the flipped cone must have its
  complement bit toggled. Edges entering the cone from outside are
  unchanged. Edges fully inside the cone are unchanged (both endpoints
  flip their semantics, double negation cancels).

  The "boundary inverter delta" determines whether the flip is profitable:
  for each external edge (cone -> outside), toggling its complement bit
  either ADDS or REMOVES an inverter. Sum the deltas; flip iff < 0.
*/

#pragma once

#include <cstdint>
#include <unordered_set>
#include <vector>

#include "../traits.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/fanout_view.hpp"

namespace mockturtle
{

struct mmig_cone_polarity_flip_params
{
  /* Maximum cone size (in gates) considered for a single flip attempt.
   * Larger cones explore more space but cost more per attempt.            */
  uint32_t max_cone_size{ 8u };

  /* Minimum required inverter reduction to accept a flip. Values > 0 add
   * a margin so we don't churn on neutral flips.                          */
  int32_t min_gain{ 1 };

  /* Iterate until no more profitable flips found, or this many sweeps.   */
  uint32_t max_sweeps{ 4u };

  bool verbose{ false };
};

struct mmig_cone_polarity_flip_stats
{
  uint32_t num_cones_examined{ 0u };
  uint32_t num_cones_flipped{ 0u };
  int32_t total_inverter_gain{ 0 };
  stopwatch<>::duration time_total{ 0 };
};

namespace detail
{

template<class Ntk>
class mmig_cone_polarity_flip_impl
{
public:
  mmig_cone_polarity_flip_impl( Ntk& ntk,
                                mmig_cone_polarity_flip_params const& ps,
                                mmig_cone_polarity_flip_stats& st )
      : _ntk( ntk ), _ps( ps ), _st( st )
  {
  }

  void run()
  {
    stopwatch t( _st.time_total );

    for ( uint32_t sweep = 0u; sweep < _ps.max_sweeps; ++sweep )
    {
      bool changed = false;
      _ntk.foreach_gate( [&]( auto const& seed ) {
        std::unordered_set<node<Ntk>> cone;
        if ( !grow_cone( seed, cone ) )
        {
          return true;
        }
        ++_st.num_cones_examined;

        int32_t const delta = boundary_inverter_delta( cone );
        if ( delta <= -_ps.min_gain )
        {
          flip_cone( cone );
          ++_st.num_cones_flipped;
          _st.total_inverter_gain += -delta;
          changed = true;
        }
        return true;
      } );

      if ( !changed )
      {
        break;
      }
    }
  }

private:
  /* Grow a connected cone starting at `seed` using a deterministic BFS
   * over fanin edges, capped at max_cone_size. We only include gates
   * (not PIs / constants) so the boundary is always well-defined.        */
  bool grow_cone( node<Ntk> const& seed, std::unordered_set<node<Ntk>>& cone )
  {
    if ( !_ntk.is_constant( seed ) && !_ntk.is_pi( seed ) )
    {
      cone.insert( seed );
    }
    else
    {
      return false;
    }

    std::vector<node<Ntk>> frontier{ seed };
    while ( !frontier.empty() && cone.size() < _ps.max_cone_size )
    {
      node<Ntk> const cur = frontier.back();
      frontier.pop_back();

      _ntk.foreach_fanin( cur, [&]( auto const& f ) {
        auto const child = _ntk.get_node( f );
        if ( _ntk.is_constant( child ) || _ntk.is_pi( child ) )
        {
          return;
        }
        if ( cone.size() >= _ps.max_cone_size )
        {
          return;
        }
        if ( cone.insert( child ).second )
        {
          frontier.push_back( child );
        }
      } );
    }

    return cone.size() >= 2u;
  }

  /* For each edge (cone-node -> outside-consumer), flipping the cone
   * toggles the consumer's complement bit on that edge. A complemented
   * edge would become non-complemented (-1 inverter) and vice versa
   * (+1 inverter). Sum across all such edges, INCLUDING POs.            */
  int32_t boundary_inverter_delta( std::unordered_set<node<Ntk>> const& cone )
  {
    int32_t delta = 0;

    for ( auto const& n : cone )
    {
      _ntk.foreach_fanout( n, [&]( auto const& parent ) {
        if ( cone.count( parent ) )
        {
          return;
        }
        _ntk.foreach_fanin( parent, [&]( auto const& f ) {
          if ( _ntk.get_node( f ) == n )
          {
            delta += _ntk.is_complemented( f ) ? -1 : 1;
          }
        } );
      } );
    }

    _ntk.foreach_po( [&]( auto const& s, auto ) {
      if ( cone.count( _ntk.get_node( s ) ) )
      {
        delta += _ntk.is_complemented( s ) ? -1 : 1;
      }
    } );

    return delta;
  }

  /* Flip every node in the cone: replace MAJ with MIN (and vice versa)
   * keeping the same fanin signals. Toggle the complement bit on every
   * edge that exits the cone. Internal edges flip on both ends and
   * therefore stay semantically identical.                              */
  void flip_cone( std::unordered_set<node<Ntk>> const& cone )
  {
    if constexpr ( !has_is_min_v<Ntk> || !has_create_min_v<Ntk> )
    {
      return;
    }

    for ( auto const& n : cone )
    {
      signal<Ntk> fanins[3]{};
      _ntk.foreach_fanin( n, [&]( auto const& f, auto i ) { fanins[i] = f; } );

      signal<Ntk> replacement{};
      if ( _ntk.is_min( n ) )
      {
        replacement = !_ntk.create_maj( fanins[0], fanins[1], fanins[2] );
      }
      else
      {
        replacement = !_ntk.create_min( fanins[0], fanins[1], fanins[2] );
      }

      _ntk.substitute_node( n, replacement );
      _ntk.replace_in_outputs( n, replacement );
    }
  }

  fanout_view<Ntk> _ntk;
  mmig_cone_polarity_flip_params const& _ps;
  mmig_cone_polarity_flip_stats& _st;
};

} // namespace detail

template<class Ntk>
void mmig_cone_polarity_flip( Ntk& ntk,
                              mmig_cone_polarity_flip_params const& ps = {},
                              mmig_cone_polarity_flip_stats* pst = nullptr )
{
  static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement foreach_gate" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement foreach_fanin" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement foreach_po" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement get_node" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement is_complemented" );
  static_assert( has_create_maj_v<Ntk>, "Ntk does not implement create_maj" );
  static_assert( has_substitute_node_v<Ntk>, "Ntk does not implement substitute_node" );
  static_assert( has_replace_in_outputs_v<Ntk>, "Ntk does not implement replace_in_outputs" );

  mmig_cone_polarity_flip_stats st;
  detail::mmig_cone_polarity_flip_impl<Ntk> impl( ntk, ps, st );
  impl.run();
  if ( pst != nullptr )
  {
    *pst = st;
  }
}

} // namespace mockturtle
