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
  \file mmig_exact_rewriting.hpp
  \brief mMIG-oriented exact-library rewriting wrapper
*/

#pragma once

#include "cleanup.hpp"
#include "node_resynthesis/mig_npn.hpp"
#include "rewrite.hpp"

#include "../traits.hpp"
#include "../utils/stopwatch.hpp"
#include "../utils/tech_library.hpp"
#include "../views/depth_view.hpp"

#include <cstdint>
#include <limits>
#include <type_traits>

namespace mockturtle
{

struct mmig_exact_rewriting_params
{
  mmig_exact_rewriting_params()
  {
    rw_ps.cut_enumeration_ps.cut_size = 4u;
    rw_ps.cut_enumeration_ps.cut_limit = 10u;
    rw_ps.cut_enumeration_ps.minimize_truth_table = true;
    rw_ps.preserve_depth = true;
    rw_ps.allow_multiple_structures = true;
    rw_ps.allow_zero_gain = false;
    rw_ps.use_dont_cares = false;
    rw_ps.window_size = 8u;
    rw_ps.verbose = false;
  }

  rewrite_params rw_ps{};
  bool preserve_depth{ true };
  bool depth_rollback_on_regression{ true };
  bool gate_enable{ false };
  uint32_t gate_max_gates{ 0u };
  uint32_t gate_min_depth{ 0u };
  uint32_t gate_min_inverted_edges{ 0u };
  bool gate_require_minority{ false };
  bool verbose{ false };
};

struct mmig_exact_rewriting_stats
{
  stopwatch<>::duration time_total{ 0 };
  rewrite_stats rewrite{};
  uint32_t num_depth_rollbacks{ 0 };
  uint32_t num_gate_skips{ 0 };
  uint32_t depth_before{ 0 };
  uint32_t depth_after{ 0 };
  uint32_t observed_gates{ 0u };
  uint32_t observed_minority{ 0u };
  uint32_t observed_inverted_edges{ 0u };
};

namespace detail
{

inline exact_library<mig_network>& mmig_exact_library_instance()
{
  static auto lib = []() {
    mig_npn_resynthesis npn_resyn;
    return exact_library<mig_network>( npn_resyn );
  }();
  return lib;
}

template<class Ntk>
uint32_t network_depth_exact( Ntk const& ntk )
{
  depth_view<Ntk> dv{ ntk };
  return dv.depth();
}

template<class Ntk>
uint32_t count_minority_exact( Ntk const& ntk )
{
  uint32_t count = 0u;
  ntk.foreach_gate( [&]( auto const& n ) {
    if ( ntk.is_min( n ) )
    {
      ++count;
    }
  } );
  return count;
}

template<class Ntk>
uint32_t count_inverted_edges_exact( Ntk const& ntk )
{
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
  return static_cast<uint32_t>( std::min<uint64_t>( inverted, std::numeric_limits<uint32_t>::max() ) );
}

} // namespace detail

template<class Ntk>
void mmig_exact_rewriting( Ntk& ntk,
                           mmig_exact_rewriting_params const& ps = {},
                           mmig_exact_rewriting_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( std::is_same_v<typename Ntk::base_type, mig_network>, "Network type is not MIG-based" );
  static_assert( has_is_min_v<Ntk>, "Ntk does not implement is_min" );
  static_assert( has_clone_v<Ntk>, "Ntk does not implement clone" );

  mmig_exact_rewriting_stats st{};
  {
    stopwatch t( st.time_total );
    st.depth_before = detail::network_depth_exact( ntk );
    st.observed_gates = ntk.num_gates();
    st.observed_minority = detail::count_minority_exact( ntk );
    st.observed_inverted_edges = detail::count_inverted_edges_exact( ntk );

    if ( ps.gate_enable )
    {
      bool blocked = false;
      if ( ps.gate_max_gates > 0u && st.observed_gates > ps.gate_max_gates )
      {
        blocked = true;
      }
      if ( ps.gate_min_depth > 0u && st.depth_before < ps.gate_min_depth )
      {
        blocked = true;
      }
      if ( ps.gate_min_inverted_edges > 0u && st.observed_inverted_edges < ps.gate_min_inverted_edges )
      {
        blocked = true;
      }
      if ( ps.gate_require_minority && st.observed_minority == 0u )
      {
        blocked = true;
      }

      if ( blocked )
      {
        ++st.num_gate_skips;
        st.depth_after = st.depth_before;
        if ( pst != nullptr )
        {
          *pst = st;
        }
        return;
      }
    }

    auto local_ps = ps.rw_ps;
    local_ps.preserve_depth = ps.preserve_depth;
    if ( ps.verbose )
    {
      local_ps.verbose = true;
    }

    auto backup = ntk.clone();

    auto& lib = detail::mmig_exact_library_instance();
    rewrite( ntk, lib, local_ps, &st.rewrite );
    ntk = cleanup_dangling( ntk );
    st.depth_after = detail::network_depth_exact( ntk );

    if ( ps.preserve_depth && ps.depth_rollback_on_regression && st.depth_after > st.depth_before )
    {
      ++st.num_depth_rollbacks;
      ntk = std::move( backup );
      st.depth_after = st.depth_before;
    }
  }

  if ( pst != nullptr )
  {
    *pst = st;
  }
}

} // namespace mockturtle
