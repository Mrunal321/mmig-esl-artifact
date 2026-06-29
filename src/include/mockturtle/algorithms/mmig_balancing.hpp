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
  \file mmig_balancing.hpp
  \brief mMIG-oriented balancing wrapper
*/

#pragma once

#include "balancing.hpp"
#include "cleanup.hpp"

#include "../networks/mig.hpp"
#include "../traits.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/depth_view.hpp"

#include <cstdint>
#include <type_traits>

namespace mockturtle
{

struct mmig_balancing_params
{
  mmig_balancing_params()
  {
    bal_ps.cut_enumeration_ps.cut_size = 4u;
    bal_ps.cut_enumeration_ps.cut_limit = 16u;
    bal_ps.only_on_critical_path = false;
    bal_ps.progress = false;
    bal_ps.verbose = false;
  }

  balancing_params bal_ps{};
  bool preserve_depth{ true };
  bool depth_rollback_on_regression{ true };
  bool verbose{ false };
};

struct mmig_balancing_stats
{
  stopwatch<>::duration time_total{ 0 };
  balancing_stats balancing{};
  uint32_t num_depth_rollbacks{ 0 };
  uint32_t depth_before{ 0 };
  uint32_t depth_after{ 0 };
};

namespace detail
{

template<class Ntk>
uint32_t network_depth_balance( Ntk const& ntk )
{
  depth_view<Ntk> dv{ ntk };
  return dv.depth();
}

} // namespace detail

template<class Ntk>
void mmig_balancing( Ntk& ntk,
                     mmig_balancing_params const& ps = {},
                     mmig_balancing_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( std::is_same_v<typename Ntk::base_type, mig_network>, "Network type is not MIG-based" );
  static_assert( has_is_min_v<Ntk>, "Ntk does not implement is_min" );
  static_assert( has_clone_v<Ntk>, "Ntk does not implement clone" );

  mmig_balancing_stats st{};
  {
    stopwatch t( st.time_total );
    st.depth_before = detail::network_depth_balance( ntk );

    auto local_ps = ps.bal_ps;
    if ( ps.verbose )
    {
      local_ps.verbose = true;
    }

    auto backup = ntk.clone();

    sop_rebalancing<Ntk> balance_fn;
    rebalancing_function_t<Ntk> fn = balance_fn;
    auto balanced = balancing( ntk, fn, local_ps, &st.balancing );
    balanced = cleanup_dangling( balanced );
    st.depth_after = detail::network_depth_balance( balanced );

    if ( ps.preserve_depth && ps.depth_rollback_on_regression && st.depth_after > st.depth_before )
    {
      ++st.num_depth_rollbacks;
      ntk = std::move( backup );
      st.depth_after = st.depth_before;
    }
    else
    {
      ntk = std::move( balanced );
    }
  }

  if ( pst != nullptr )
  {
    *pst = st;
  }
}

} // namespace mockturtle
