#include <catch.hpp>

#include <mockturtle/algorithms/equivalence_checking.hpp>
#include <mockturtle/algorithms/miter.hpp>
#include <mockturtle/algorithms/mmig_balancing.hpp>
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

TEST_CASE( "mmig balancing preserves equivalence on mixed MAJ/MIN network", "[mmig_balancing]" )
{
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();

  auto const n1 = ntk.create_maj( a, b, c );
  auto const n2 = ntk.create_min( n1, d, !c );
  auto const n3 = ntk.create_maj( n2, e, b );
  auto const n4 = ntk.create_min( n3, !a, d );
  auto const n5 = ntk.create_maj( n4, n2, a );
  ntk.create_po( n3 );
  ntk.create_po( n4 );
  ntk.create_po( n5 );

  auto const original = ntk.clone();

  mmig_balancing_params ps;
  ps.bal_ps.cut_enumeration_ps.cut_size = 4u;
  ps.bal_ps.cut_enumeration_ps.cut_limit = 16u;
  ps.bal_ps.only_on_critical_path = true;
  ps.preserve_depth = true;
  ps.depth_rollback_on_regression = true;

  mmig_balancing_stats st{};
  mmig_balancing( ntk, ps, &st );

  CHECK( st.depth_after <= st.depth_before );
  check_equivalent( original, ntk );
}

