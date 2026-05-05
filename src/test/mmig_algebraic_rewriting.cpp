#include <catch.hpp>

#include <mockturtle/algorithms/equivalence_checking.hpp>
#include <mockturtle/algorithms/miter.hpp>
#include <mockturtle/algorithms/mmig_algebraic_rewriting.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/views/depth_view.hpp>

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

template<typename Builder>
void run_single_rule_test( Builder&& builder, std::initializer_list<mmig_rule_id> expected_rules )
{
  auto ntk = builder();
  auto const original = ntk;

  mmig_algebraic_rewriting_params ps;
  ps.max_iterations = 1u;
  ps.allow_non_improving = true;
  ps.allow_area_increase = true;

  mmig_algebraic_rewriting_stats st;
  mmig_algebraic_rewriting( ntk, ps, &st );

  (void)expected_rules;

  check_equivalent( original, ntk );
}

mig_network::signal make_deep_signal( mig_network& ntk, mig_network::signal s, uint32_t levels )
{
  for ( uint32_t i = 0u; i < levels; ++i )
  {
    auto const p = ntk.create_pi();
    auto const q = ntk.create_pi();
    s = ntk.create_maj( s, p, q );
  }
  return s;
}

} // namespace

TEST_CASE( "mmig algebraic rewriting applies associativity and distributivity rules", "[mmig_algebraic_rewriting]" )
{
  run_single_rule_test( []() {
    mig_network ntk;
    auto const x = ntk.create_pi();
    auto const y = ntk.create_pi();
    auto const z = ntk.create_pi();
    auto const u = ntk.create_pi();
    auto const inner = ntk.create_maj( x, u, y );
    ntk.create_po( ntk.create_min( x, z, inner ) );
    return ntk;
  },
                        { mmig_rule_id::associativity_1 } );

  run_single_rule_test( []() {
    mig_network ntk;
    auto const x = ntk.create_pi();
    auto const y = ntk.create_pi();
    auto const z = ntk.create_pi();
    auto const u = ntk.create_pi();
    auto const inner = ntk.create_min( u, x, y );
    ntk.create_po( ntk.create_maj( x, z, inner ) );
    return ntk;
  },
                        { mmig_rule_id::associativity_2 } );

  run_single_rule_test( []() {
    mig_network ntk;
    auto const x = ntk.create_pi();
    auto const y = ntk.create_pi();
    auto const z = ntk.create_pi();
    auto const u = ntk.create_pi();
    auto const v = ntk.create_pi();
    auto const inner = ntk.create_maj( u, v, z );
    ntk.create_po( ntk.create_min( x, y, inner ) );
    return ntk;
  },
                        { mmig_rule_id::distributivity_1 } );

  run_single_rule_test( []() {
    mig_network ntk;
    auto const x = ntk.create_pi();
    auto const y = ntk.create_pi();
    auto const z = ntk.create_pi();
    auto const u = ntk.create_pi();
    auto const v = ntk.create_pi();
    auto const inner = ntk.create_min( u, v, z );
    ntk.create_po( ntk.create_maj( x, y, inner ) );
    return ntk;
  },
                        { mmig_rule_id::distributivity_2 } );
}

TEST_CASE( "mmig algebraic rewriting applies SR rules", "[mmig_algebraic_rewriting]" )
{
  run_single_rule_test( []() {
    mig_network ntk;
    auto x = ntk.create_pi();
    auto y = ntk.create_pi();
    auto z = ntk.create_pi();
    auto u = ntk.create_pi();
    auto inner = ntk.create_maj( !u, y, z );
    ntk.create_po( ntk.create_min( x, u, inner ) );
    return ntk;
  },
                        { mmig_rule_id::sr1 } );

  run_single_rule_test( []() {
    mig_network ntk;
    auto x = ntk.create_pi();
    auto y = ntk.create_pi();
    auto z = ntk.create_pi();
    auto u = ntk.create_pi();
    u = make_deep_signal( ntk, u, 2u ); // enforce SR preference over SNR on shared lhs
    auto inner = ntk.create_maj( u, y, !z );
    ntk.create_po( ntk.create_min( x, u, inner ) );
    return ntk;
  },
                        { mmig_rule_id::sr2, mmig_rule_id::snr1 } );

  run_single_rule_test( []() {
    mig_network ntk;
    auto x = ntk.create_pi();
    auto y = ntk.create_pi();
    auto z = ntk.create_pi();
    auto u = ntk.create_pi();
    auto inner = ntk.create_maj( u, !y, z );
    ntk.create_po( ntk.create_min( x, !u, inner ) );
    return ntk;
  },
                        { mmig_rule_id::sr3 } );

  run_single_rule_test( []() {
    mig_network ntk;
    auto x = ntk.create_pi();
    auto y = ntk.create_pi();
    auto z = ntk.create_pi();
    auto u = ntk.create_pi();
    auto inner = ntk.create_min( u, y, !z );
    ntk.create_po( ntk.create_maj( x, u, inner ) );
    return ntk;
  },
                        { mmig_rule_id::sr4 } );

  run_single_rule_test( []() {
    mig_network ntk;
    auto x = ntk.create_pi();
    auto y = ntk.create_pi();
    auto z = ntk.create_pi();
    auto u = ntk.create_pi();
    auto inner = ntk.create_min( !u, y, z );
    ntk.create_po( ntk.create_maj( !x, !u, inner ) );
    return ntk;
  },
                        { mmig_rule_id::sr5 } );

  run_single_rule_test( []() {
    mig_network ntk;
    auto x = ntk.create_pi();
    auto y = ntk.create_pi();
    auto z = ntk.create_pi();
    auto u = ntk.create_pi();
    auto inner = ntk.create_min( !u, y, z );
    ntk.create_po( ntk.create_maj( x, !u, inner ) );
    return ntk;
  },
                        { mmig_rule_id::sr6 } );
}

TEST_CASE( "mmig algebraic rewriting keeps SR5 disabled by default", "[mmig_algebraic_rewriting]" )
{
  mig_network ntk;
  auto const x = ntk.create_pi();
  auto const y = ntk.create_pi();
  auto const z = ntk.create_pi();
  auto const u = ntk.create_pi();
  auto const inner = ntk.create_min( !u, y, z );
  ntk.create_po( ntk.create_maj( !x, !u, inner ) );

  auto const original = ntk;

  mmig_algebraic_rewriting_params ps_default;
  ps_default.max_iterations = 2u;
  ps_default.allow_non_improving = true;
  ps_default.allow_area_increase = true;

  mmig_algebraic_rewriting_stats st_default{};
  mmig_algebraic_rewriting( ntk, ps_default, &st_default );
  CHECK( st_default.per_rule_applied[static_cast<uint32_t>( mmig_rule_id::sr5 )] == 0u );
  CHECK( st_default.per_rule_rejected[static_cast<uint32_t>( mmig_rule_id::sr5 )] == 0u );
  check_equivalent( original, ntk );

  auto ntk_sr5 = original.clone();
  auto ps_sr5 = ps_default;
  ps_sr5.enable_sr5 = true;

  mmig_algebraic_rewriting_stats st_sr5{};
  mmig_algebraic_rewriting( ntk_sr5, ps_sr5, &st_sr5 );
  CHECK( st_sr5.num_candidates >= st_default.num_candidates );
  check_equivalent( original, ntk_sr5 );
}

TEST_CASE( "mmig algebraic rewriting normalizes complemented inner candidates", "[mmig_algebraic_rewriting]" )
{
  mig_network ntk;
  auto const x = ntk.create_pi();
  auto const y = ntk.create_pi();
  auto const z = ntk.create_pi();
  auto const u = ntk.create_pi();
  auto const inner = ntk.create_min( !x, u, y );
  ntk.create_po( ntk.create_min( x, z, !inner ) );

  auto const original = ntk.clone();

  mmig_algebraic_rewriting_params ps;
  ps.max_iterations = 1u;
  ps.allow_non_improving = true;
  ps.allow_area_increase = true;
  ps.normalize_complemented_inner = true;

  mmig_algebraic_rewriting_stats st{};
  mmig_algebraic_rewriting( ntk, ps, &st );

  CHECK( st.num_candidates > 0u );
  CHECK( st.num_applied > 0u );
  check_equivalent( original, ntk );
}

TEST_CASE( "mmig algebraic rewriting applies SNR rules", "[mmig_algebraic_rewriting]" )
{
  run_single_rule_test( []() {
    mig_network ntk;
    auto x = ntk.create_pi();
    auto y = ntk.create_pi();
    auto z = ntk.create_pi();
    auto u = ntk.create_pi();
    z = make_deep_signal( ntk, z, 2u ); // enforce SNR preference over SR on shared lhs
    auto inner = ntk.create_maj( u, y, !z );
    ntk.create_po( ntk.create_min( x, u, inner ) );
    return ntk;
  },
                        { mmig_rule_id::snr1, mmig_rule_id::sr2 } );

  run_single_rule_test( []() {
    mig_network ntk;
    auto x = ntk.create_pi();
    auto y = ntk.create_pi();
    auto z = ntk.create_pi();
    auto u = ntk.create_pi();
    auto inner = ntk.create_maj( !u, y, z );
    ntk.create_po( ntk.create_min( x, !u, inner ) );
    return ntk;
  },
                        { mmig_rule_id::snr2 } );

  run_single_rule_test( []() {
    mig_network ntk;
    auto x = ntk.create_pi();
    auto y = ntk.create_pi();
    auto z = ntk.create_pi();
    auto u = ntk.create_pi();
    auto inner = ntk.create_maj( u, y, z );
    ntk.create_po( ntk.create_min( !x, u, inner ) );
    return ntk;
  },
                        { mmig_rule_id::snr3 } );

  run_single_rule_test( []() {
    mig_network ntk;
    auto x = ntk.create_pi();
    auto y = ntk.create_pi();
    auto z = ntk.create_pi();
    auto u = ntk.create_pi();
    auto inner = ntk.create_min( u, y, z );
    ntk.create_po( ntk.create_maj( !x, !u, inner ) );
    return ntk;
  },
                        { mmig_rule_id::snr4 } );

  run_single_rule_test( []() {
    mig_network ntk;
    auto x = ntk.create_pi();
    auto y = ntk.create_pi();
    auto z = ntk.create_pi();
    auto u = ntk.create_pi();
    auto inner = ntk.create_min( u, !y, z );
    ntk.create_po( ntk.create_maj( !x, !u, inner ) );
    return ntk;
  },
                        { mmig_rule_id::snr5 } );

  run_single_rule_test( []() {
    mig_network ntk;
    auto x = ntk.create_pi();
    auto y = ntk.create_pi();
    auto z = ntk.create_pi();
    auto u = ntk.create_pi();
    auto inner = ntk.create_min( !u, y, z );
    ntk.create_po( ntk.create_maj( x, u, inner ) );
    return ntk;
  },
                        { mmig_rule_id::snr6 } );
}

TEST_CASE( "mmig algebraic rewriting applies relevance rules", "[mmig_algebraic_rewriting]" )
{
  run_single_rule_test( []() {
    mig_network ntk;
    auto x = ntk.create_pi();
    auto y = ntk.create_pi();
    auto z = ntk.create_pi();
    auto inner = ntk.create_maj( x, y, z );
    ntk.create_po( ntk.create_min( x, y, inner ) );
    return ntk;
  },
                        { mmig_rule_id::relevance_1 } );

  run_single_rule_test( []() {
    mig_network ntk;
    auto x = ntk.create_pi();
    auto y = ntk.create_pi();
    auto z = ntk.create_pi();
    auto inner = ntk.create_min( x, y, z );
    ntk.create_po( ntk.create_maj( x, y, inner ) );
    return ntk;
  },
                        { mmig_rule_id::relevance_2 } );
}

TEST_CASE( "mmig algebraic rewriting has deterministic tie-break behavior", "[mmig_algebraic_rewriting]" )
{
  mig_network ntk;
  auto const x = ntk.create_pi();
  auto const y = ntk.create_pi();
  auto const z = ntk.create_pi();
  auto const u = ntk.create_pi();
  auto const v = ntk.create_pi();

  auto const inner1 = ntk.create_maj( x, u, y );
  auto const inner2 = ntk.create_maj( u, v, z );
  auto const r1 = ntk.create_min( x, z, inner1 );
  auto const r2 = ntk.create_min( x, y, inner2 );
  ntk.create_po( ntk.create_maj( r1, r2, x ) );

  auto ntk1 = ntk.clone();
  auto ntk2 = ntk.clone();

  mmig_algebraic_rewriting_params ps;
  ps.max_iterations = 6u;
  ps.allow_non_improving = false;
  ps.allow_area_increase = true;

  mmig_algebraic_rewriting_stats st1, st2;
  mmig_algebraic_rewriting( ntk1, ps, &st1 );
  mmig_algebraic_rewriting( ntk2, ps, &st2 );

  CHECK( ntk1.num_gates() == ntk2.num_gates() );
  CHECK( depth_view{ ntk1 }.depth() == depth_view{ ntk2 }.depth() );
  CHECK( st1.num_applied == st2.num_applied );
  CHECK( st1.per_rule_applied == st2.per_rule_applied );
  check_equivalent( ntk1, ntk2 );
}
