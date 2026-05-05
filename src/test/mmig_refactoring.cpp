#include <catch.hpp>

#include <mockturtle/algorithms/equivalence_checking.hpp>
#include <mockturtle/algorithms/miter.hpp>
#include <mockturtle/algorithms/mmig_refactoring.hpp>
#include <mockturtle/networks/mig.hpp>

using namespace mockturtle;

namespace
{

void check_equivalent( mig_network const& a, mig_network const& b )
{
  auto const m = miter<mig_network>( a, b );
  REQUIRE( m.has_value() );
  auto const eq = equivalence_checking( *m );
  REQUIRE( eq.has_value() );
  CHECK( *eq );
}

} // namespace

TEST_CASE( "mmig refactoring preserves equivalence on mixed MAJ/MIN network", "[mmig_refactoring]" )
{
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();

  auto const n1 = ntk.create_maj( a, b, c );
  auto const n2 = ntk.create_min( n1, d, !c );
  auto const n3 = ntk.create_maj( n1, n2, e );
  auto const n4 = ntk.create_min( n3, b, d );
  ntk.create_po( n3 );
  ntk.create_po( n4 );

  auto const original = ntk.clone();

  mmig_refactoring_params ps;
  ps.refac_ps.max_pis = 6u;
  ps.refac_ps.use_reconvergence_cut = true;
  ps.refac_ps.allow_zero_gain = false;
  ps.refac_ps.use_dont_cares = false;
  ps.preserve_depth = true;
  ps.depth_rollback_on_regression = true;

  mmig_refactoring_stats st{};
  mmig_refactoring( ntk, ps, &st );

  CHECK( st.depth_after <= st.depth_before );
  check_equivalent( original, ntk );
}

