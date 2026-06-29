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
  \file mmig_cut_rewriting.hpp
  \brief mMIG-oriented cut rewriting wrapper
*/

#pragma once

#include "cleanup.hpp"
#include "cut_rewriting.hpp"
#include "node_resynthesis/mig_npn.hpp"

#include "../traits.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/depth_view.hpp"

#include <cstdint>
#include <type_traits>

namespace mockturtle
{

struct mmig_cut_rewriting_params
{
  mmig_cut_rewriting_params()
  {
    cut_ps.cut_enumeration_ps.cut_size = 4u;
    cut_ps.cut_enumeration_ps.cut_limit = 10u;
    cut_ps.allow_zero_gain = false;
    cut_ps.preserve_depth = true;
    cut_ps.verbose = false;
  }

  cut_rewriting_params cut_ps{};
  bool depth_rollback_on_regression{ true };
  bool verbose{ false };
};

struct mmig_cut_rewriting_stats
{
  stopwatch<>::duration time_total{ 0 };
  cut_rewriting_stats cut{};
  uint32_t num_depth_rollbacks{ 0 };
  uint32_t depth_before{ 0 };
  uint32_t depth_after{ 0 };
};

namespace detail
{

template<class Ntk>
uint32_t network_depth( Ntk const& ntk )
{
  depth_view<Ntk> dv{ ntk };
  return dv.depth();
}

} // namespace detail

template<class Ntk>
void mmig_cut_rewriting( Ntk& ntk,
                         mmig_cut_rewriting_params const& ps = {},
                         mmig_cut_rewriting_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( std::is_same_v<typename Ntk::base_type, mig_network>, "Network type is not MIG-based" );
  static_assert( has_is_min_v<Ntk>, "Ntk does not implement is_min" );
  static_assert( has_clone_v<Ntk>, "Ntk does not implement clone" );

  mmig_cut_rewriting_stats st{};
  {
    stopwatch t( st.time_total );
    st.depth_before = detail::network_depth( ntk );

    auto local_ps = ps.cut_ps;
    if ( ps.verbose )
    {
      local_ps.verbose = true;
    }

    auto backup = ntk.clone();
    mig_npn_resynthesis resyn;
    auto rewritten = cut_rewriting<Ntk>( ntk, resyn, local_ps, &st.cut );
    rewritten = cleanup_dangling( rewritten );
    st.depth_after = detail::network_depth( rewritten );

    if ( local_ps.preserve_depth && ps.depth_rollback_on_regression && st.depth_after > st.depth_before )
    {
      ++st.num_depth_rollbacks;
      ntk = std::move( backup );
      st.depth_after = st.depth_before;
    }
    else
    {
      ntk = std::move( rewritten );
    }
  }

  if ( pst != nullptr )
  {
    *pst = st;
  }
}

} // namespace mockturtle
