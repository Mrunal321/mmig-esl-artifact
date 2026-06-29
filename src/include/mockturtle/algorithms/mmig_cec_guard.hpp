#pragma once

#include <cstdint>
#include <functional>
#include <iostream>
#include <optional>
#include <string>
#include <utility>

#include "cleanup.hpp"
#include "equivalence_checking.hpp"
#include "miter.hpp"
#include "../traits.hpp"
#include "../utils/stopwatch.hpp"

namespace mockturtle
{

struct mmig_cec_guard_params
{
  bool run_cec = true;
  bool verbose = false;
};

struct mmig_cec_guard_stats
{
  uint32_t accepted = 0;
  uint32_t rejected = 0;
  uint32_t inconclusive = 0;
  stopwatch<>::duration time_cec{ 0 };
};

template<class Ntk, class Fn>
bool run_with_mmig_cec_guard( Ntk& ntk,
                              std::string const& stage_name,
                              Fn&& fn,
                              mmig_cec_guard_params const& ps = {},
                              mmig_cec_guard_stats* pst = nullptr )
{
  static_assert( has_clone_v<Ntk>, "Ntk does not implement clone" );

  mmig_cec_guard_stats st{};
  if ( pst != nullptr )
  {
    st = *pst;
  }

  auto candidate = ntk.clone();
  std::invoke( std::forward<Fn>( fn ), candidate );
  candidate = cleanup_dangling( candidate );

  if ( !ps.run_cec )
  {
    ntk = std::move( candidate );
    ++st.accepted;
    if ( ps.verbose )
    {
      std::cout << "[mmig-cec] " << stage_name << ": accepted (CEC disabled)\n";
    }
    if ( pst != nullptr )
    {
      *pst = st;
    }
    return true;
  }

  const auto maybe_miter = miter<Ntk>( ntk, candidate );
  if ( !maybe_miter.has_value() )
  {
    ++st.rejected;
    if ( ps.verbose )
    {
      std::cout << "[mmig-cec] " << stage_name << ": rejected (miter construction failed)\n";
    }
    if ( pst != nullptr )
    {
      *pst = st;
    }
    return false;
  }

  std::optional<bool> equivalent;
  {
    stopwatch t( st.time_cec );
    equivalent = equivalence_checking( *maybe_miter );
  }

  if ( !equivalent.has_value() )
  {
    ++st.inconclusive;
    ++st.rejected;
    if ( ps.verbose )
    {
      std::cout << "[mmig-cec] " << stage_name << ": rejected (inconclusive CEC)\n";
    }
    if ( pst != nullptr )
    {
      *pst = st;
    }
    return false;
  }

  if ( *equivalent )
  {
    ntk = std::move( candidate );
    ++st.accepted;
    if ( ps.verbose )
    {
      std::cout << "[mmig-cec] " << stage_name << ": accepted\n";
    }
    if ( pst != nullptr )
    {
      *pst = st;
    }
    return true;
  }

  ++st.rejected;
  if ( ps.verbose )
  {
    std::cout << "[mmig-cec] " << stage_name << ": rejected (not equivalent)\n";
  }
  if ( pst != nullptr )
  {
    *pst = st;
  }
  return false;
}

} // namespace mockturtle
