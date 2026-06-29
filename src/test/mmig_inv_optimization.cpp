#include <catch.hpp>

#include <mockturtle/algorithms/mmig_inv_optimization.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/networks/mig.hpp>

using namespace mockturtle;

namespace
{

int32_t count_inverted_edges( mig_network const& ntk )
{
  int32_t count = 0;
  ntk.foreach_gate( [&]( auto const& n ) {
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( ntk.is_complemented( f ) )
      {
        ++count;
      }
    } );
  } );
  ntk.foreach_po( [&]( auto const& s ) {
    if ( ntk.is_complemented( s ) )
    {
      ++count;
    }
  } );
  return count;
}

void check_equivalent( mig_network const& a, mig_network const& b )
{
  auto const ta = simulate<kitty::static_truth_table<5u>>( a );
  auto const tb = simulate<kitty::static_truth_table<5u>>( b );
  REQUIRE( ta.size() == tb.size() );
  for ( uint32_t i = 0u; i < ta.size(); ++i )
  {
    CHECK( ta[i] == tb[i] );
  }
}

} // namespace

TEST_CASE( "mmig inv optimization is equivalence-safe on mixed MAJ/MIN", "[mmig_inv_optimization]" )
{
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();

  auto const n = ntk.create_maj( a, b, c );
  auto const p1 = ntk.create_maj( !n, d, e );
  auto const p2 = ntk.create_maj( !n, a, d );
  auto const p3 = ntk.create_min( !n, b, e );
  auto const p4 = ntk.create_min( !n, c, d );
  auto const p5 = ntk.create_maj( !n, e, b );

  ntk.create_po( !n );
  ntk.create_po( p1 );
  ntk.create_po( p2 );
  ntk.create_po( p3 );
  ntk.create_po( p4 );
  ntk.create_po( p5 );

  auto const original = ntk.clone();
  auto const inv_before = count_inverted_edges( ntk );

  mmig_inv_optimization_stats st{};
  mmig_inv_optimization( ntk, &st );

  auto const inv_after = count_inverted_edges( ntk );
  CHECK( inv_after <= inv_before );
  check_equivalent( original, ntk );
}
