#include <catch.hpp>

#include <mockturtle/algorithms/equivalence_checking.hpp>
#include <mockturtle/algorithms/miter.hpp>
#include <mockturtle/algorithms/mmig_minority_seeding.hpp>
#include <mockturtle/networks/mig.hpp>

using namespace mockturtle;

namespace
{

uint32_t count_min_nodes( mig_network const& ntk )
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

void check_equivalent( mig_network const& a, mig_network const& b )
{
  auto const m = miter<mig_network>( a, b );
  REQUIRE( m.has_value() );
  auto const eq = equivalence_checking( *m );
  REQUIRE( eq.has_value() );
  CHECK( *eq );
}

} // namespace

TEST_CASE( "mmig minority seeding introduces minority nodes and preserves equivalence", "[mmig_minority_seeding]" )
{
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const n1 = ntk.create_maj( a, b, c );
  auto const n2 = ntk.create_maj( !a, c, d );
  auto const n3 = ntk.create_maj( n1, n2, !b );
  ntk.create_po( n3 );

  auto const original = ntk;
  CHECK( count_min_nodes( ntk ) == 0u );

  mmig_minority_seeding_params ps;
  ps.mode = mmig_seed_mode::both;
  ps.max_candidates = 4u;
  ps.max_rounds = 2u;
  ps.allow_non_improving = true;
  ps.allow_area_increase = true;

  mmig_minority_seeding_stats st{};
  mmig_minority_seeding( ntk, ps, &st );

  CHECK( st.num_applied > 0u );
  CHECK( count_min_nodes( ntk ) > 0u );
  check_equivalent( original, ntk );
}

TEST_CASE( "mmig minority seeding is deterministic for fixed params", "[mmig_minority_seeding]" )
{
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();

  auto const n1 = ntk.create_maj( a, b, c );
  auto const n2 = ntk.create_maj( b, c, d );
  auto const n3 = ntk.create_maj( n1, n2, !a );
  auto const n4 = ntk.create_maj( n1, !n2, d );
  ntk.create_po( n3 );
  ntk.create_po( n4 );

  auto ntk1 = ntk.clone();
  auto ntk2 = ntk.clone();

  mmig_minority_seeding_params ps;
  ps.mode = mmig_seed_mode::level2;
  ps.max_candidates = 3u;
  ps.max_rounds = 1u;
  ps.allow_non_improving = true;

  mmig_minority_seeding_stats st1{}, st2{};
  mmig_minority_seeding( ntk1, ps, &st1 );
  mmig_minority_seeding( ntk2, ps, &st2 );

  CHECK( ntk1.num_gates() == ntk2.num_gates() );
  CHECK( count_min_nodes( ntk1 ) == count_min_nodes( ntk2 ) );
  CHECK( st1.num_applied == st2.num_applied );
  check_equivalent( ntk1, ntk2 );
}

TEST_CASE( "mmig minority seeding filters low-L1 candidates during candidate CEC", "[mmig_minority_seeding]" )
{
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();
  auto const f = ntk.create_pi();
  auto const g = ntk.create_pi();
  auto const h = ntk.create_pi();

  auto const n1 = ntk.create_maj( a, b, c );
  auto const n2 = ntk.create_maj( d, e, f );
  auto const n3 = ntk.create_maj( n1, n2, g );
  auto const n4 = ntk.create_maj( n1, n3, h );
  auto const n5 = ntk.create_maj( n2, n4, a );
  auto const n6 = ntk.create_maj( n3, n5, b );
  ntk.create_po( n4 );
  ntk.create_po( n6 );
  ntk.create_po( ntk.create_maj( n5, n6, c ) );

  auto const original = ntk.clone();
  auto ntk_loose = original.clone();
  auto ntk_filtered = original.clone();

  mmig_minority_seeding_params ps_loose;
  ps_loose.mode = mmig_seed_mode::both;
  ps_loose.max_candidates = 24u;
  ps_loose.max_rounds = 1u;
  ps_loose.verify_each_candidate = true;
  ps_loose.prefer_l1_when_verifying = false;
  ps_loose.min_l1_for_verification = -1000;
  ps_loose.allow_non_improving = true;
  ps_loose.allow_area_increase = true;

  mmig_minority_seeding_stats st_loose{};
  mmig_minority_seeding( ntk_loose, ps_loose, &st_loose );

  auto ps_filtered = ps_loose;
  ps_filtered.prefer_l1_when_verifying = true;
  ps_filtered.min_l1_for_verification = 1;

  mmig_minority_seeding_stats st_filtered{};
  mmig_minority_seeding( ntk_filtered, ps_filtered, &st_filtered );

  CHECK( st_filtered.num_filtered_low_l1 > 0u );
  CHECK( st_filtered.num_equivalence_checks <= st_loose.num_equivalence_checks );
  check_equivalent( original, ntk_loose );
  check_equivalent( original, ntk_filtered );
}
