// bench_nf2.cc — L3 for the QCD_Nf2/no_mesons flow set.
//
// Same shape as bench.cc; separate file only because the flow class, the interpolator set and the
// config format differ. Nf2's generated headers take DiFfRG::JSONValue, which is a deprecated alias
// for ConfigTree (config_tree.hh:243), so they build against the current install unchanged — the
// earlier assumption that Nf2 was unbuildable was wrong.
//
// This is the code that matters most for the interpolator-sharing result: Nf2 does far more spline
// lookups per quadrature point than YangMills (ZAqbq4 alone does 119 over 18 distinct arguments),
// so the win here was previously only modelled from SASS counts. This measures it.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <DiFfRG/DiFfRG.hh>
#include <DiFfRG/common/config_tree.hh>

#include "flows/flows.hh"

using namespace DiFfRG;

namespace
{
  constexpr uint p_grid_size = 64;

  const char *env_or(const char *k, const char *d)
  {
    const char *v = std::getenv(k);
    return (v && *v) ? v : d;
  }

  struct Stat {
    double med, lo, hi;
  };

  Stat summarize(std::vector<double> &v)
  {
    std::sort(v.begin(), v.end());
    return {v[v.size() / 2], v.front(), v.back()};
  }
} // namespace

int main(int argc, char **argv)
{
  const std::string cfg_path = env_or("NTB_CONFIG", "parameter.json");
  const int reps = std::atoi(env_or("NTB_REPS", "11"));
  const int warmup = std::atoi(env_or("NTB_WARMUP", "3"));
  const std::string tag = env_or("NTB_TAG", "baseline");

  Init init(argc, argv, cfg_path);
  const ConfigurationHelper cfg_helper = init.get_configuration_helper();
  const ConfigTree &config = cfg_helper.get_config();

  using Coordinates1D = LogarithmicCoordinates1D<double>;
  const Coordinates1D coords(p_grid_size, config.get_double("/discretization/p_grid_min"),
                             config.get_double("/discretization/p_grid_max"),
                             config.get_double("/discretization/p_grid_bias"));

  QCDFlows flows(config);
  const double k = config.get_double("/physical/Lambda");
  const double etaQ = 0.1;
  flows.set_k(k);

  using Interp = SplineInterpolator1D<double, Coordinates1D, GPU_memory>;
  Interp ZA3(coords), ZAcbc(coords), ZA4(coords), ZAqbq1(coords), ZAqbq4(coords), ZAqbq7(coords),
      dtZc(coords), Zc(coords), dtZA(coords), ZA(coords), zq(coords), ZAAqbq1(coords),
      ZAAqbq2(coords);

  std::vector<double> buf(p_grid_size);
  auto fill_with = [&](Interp &interp, double a, double b) {
    for (uint i = 0; i < p_grid_size; ++i) {
      const double p = coords.forward(i);
      buf[i] = a + b / (1. + p * p);
    }
    interp.update(buf.data());
  };
  fill_with(ZA3, 1.2, 0.8);
  fill_with(ZAcbc, 1.1, 0.5);
  fill_with(ZA4, 1.5, 1.0);
  fill_with(ZAqbq1, 1.0, 0.6);
  fill_with(ZAqbq4, 0.9, 0.4);
  fill_with(ZAqbq7, 0.8, 0.3);
  fill_with(Zc, 1.0, 0.7);
  fill_with(ZA, 1.0, 3.0);
  fill_with(zq, 1.0, 0.2);
  fill_with(ZAAqbq1, 0.5, 0.2);
  fill_with(ZAAqbq2, 0.4, 0.15);
  fill_with(dtZA, 0.05, 0.30);
  fill_with(dtZc, 0.02, 0.15);

  const auto args = device::tie(k, etaQ, ZA3, ZAcbc, ZA4, ZAqbq1, ZAqbq4, ZAqbq7, dtZc, Zc, dtZA, ZA,
                                zq, ZAAqbq1, ZAAqbq2);

  std::vector<double> dest(p_grid_size, 0.);

  auto time_flow = [&](const char *name, auto &integrator) {
    std::vector<double> ms;
    ms.reserve(reps);
    for (int r = 0; r < warmup + reps; ++r) {
      Kokkos::fence();
      const auto t0 = std::chrono::steady_clock::now();
      integrator.map(dest.data(), coords, args);
      Kokkos::fence();
      const auto t1 = std::chrono::steady_clock::now();
      if (r >= warmup) ms.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    const Stat s = summarize(ms);
    double sum = 0., amax = 0.;
    for (double d : dest) {
      sum += d;
      amax = std::max(amax, std::abs(d));
    }
    std::printf("%-8s %-10s %9.4f %9.4f %9.4f   %+.15e %.6e\n", tag.c_str(), name, s.med, s.lo, s.hi,
                sum, amax);
    std::fflush(stdout);
  };

  // Fixed wall-clock burn to get past the clock ramp — see bench.cc for why this must not be an
  // "until it stops improving" heuristic.
  {
    const double burn_ms = std::atof(env_or("NTB_BURN_MS", "2500"));
    const auto t0 = std::chrono::steady_clock::now();
    while (std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() <
           burn_ms) {
      flows.ZA4.map(dest.data(), coords, args);
      Kokkos::fence();
    }
  }

  std::printf("%-8s %-10s %9s %9s %9s   %-22s %s\n", "tag", "flow", "med_ms", "min_ms", "max_ms",
              "checksum", "max_abs");
  time_flow("Zc", flows.Zc);
  time_flow("ZA", flows.ZA);
  time_flow("Zq", flows.Zq);
  time_flow("ZAcbc", flows.ZAcbc);
  time_flow("ZA3", flows.ZA3);
  time_flow("ZAqbq1", flows.ZAqbq1);
  time_flow("ZAqbq7", flows.ZAqbq7);
  time_flow("ZAqbq4", flows.ZAqbq4);
  time_flow("ZAAqbq1", flows.ZAAqbq1);
  time_flow("ZAAqbq2", flows.ZAAqbq2);
  time_flow("ZA4", flows.ZA4);

  return 0;
}
