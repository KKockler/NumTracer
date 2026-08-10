// bench.cc — L3: per-flow map() timing for the vendored YangMills flow set.
//
// This drives the REAL DiFfRG integrator (QuadratureIntegrator via Integrator_p2_*) over the REAL
// production grid orders read from the vendored parameter.toml, so a number here transfers to the
// production solve. What it deliberately does NOT do is run the ODE: the point is to isolate the
// cost of one map() per flow, which the full solve buries under 240 timesteps.
//
// The interpolators are filled with a smooth analytic profile rather than physical dressings. That
// is legitimate for timing (the spline lookup cost is data-independent) and it makes the harness
// self-contained and deterministic. It is NOT legitimate for physics, and nothing here claims to be.
//
// Reported: median of NTB_REPS timed map() calls with a warmup, plus the min/max spread. On a
// 115 W laptop part the spread is the honest error bar — a delta smaller than it is not a result,
// which is why L2 (SASS) is the grade for small deltas.

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
  const std::string cfg_path = env_or("NTB_CONFIG", "parameter.toml");
  const int reps = std::atoi(env_or("NTB_REPS", "11"));
  const int warmup = std::atoi(env_or("NTB_WARMUP", "3"));
  const std::string tag = env_or("NTB_TAG", "baseline");

  auto note = [](const char *s) { std::fprintf(stderr, "[ntb] %s\n", s); std::fflush(stderr); };

  // Init owns the deal.II/MPI/Kokkos bring-up; QuadratureProvider refuses to construct without it.
  note("Init");
  Init init(argc, argv, cfg_path);
  note("config");
  // get_configuration_helper() returns BY VALUE — the helper has to be kept alive, otherwise the
  // ConfigTree reference dangles the moment the full-expression ends.
  const ConfigurationHelper cfg_helper = init.get_configuration_helper();
  const ConfigTree &config = cfg_helper.get_config();

  using Coordinates1D = LogarithmicCoordinates1D<double>;
  const Coordinates1D coords(p_grid_size, config.get_double("/discretization/p_grid_min"),
                             config.get_double("/discretization/p_grid_max"),
                             config.get_double("/discretization/p_grid_bias"));

  note("coords ok, constructing flows");
  YangMillsFlows flows(config);
  note("flows constructed");
  const double k = config.get_double("/physical/Lambda");
  flows.set_k(k);

  // --- fixtures ------------------------------------------------------------------------------
  // A smooth, strictly positive profile on the production grid. Deterministic, no RNG, so two
  // harness runs are comparable and a variant's output can be diffed against the baseline's.
  SplineInterpolator1D<double, Coordinates1D, GPU_memory> dtZc(coords), dtZA(coords), ZA(coords),
      Zc(coords), ZA4(coords), ZAcbc(coords), ZA3(coords);

  std::vector<double> buf(p_grid_size);
  auto fill_with = [&](auto &interp, double a, double b) {
    for (uint i = 0; i < p_grid_size; ++i) {
      const double p = coords.forward(i);
      buf[i] = a + b / (1. + p * p);
    }
    interp.update(buf.data());
  };
  fill_with(ZA3, 1.2, 0.8);
  fill_with(ZAcbc, 1.1, 0.5);
  fill_with(ZA4, 1.5, 1.0);
  fill_with(ZA, 1.0, 3.0);
  fill_with(Zc, 1.0, 0.7);
  fill_with(dtZA, 0.05, 0.30);
  fill_with(dtZc, 0.02, 0.15);

  note("fixtures done");
  const auto args = device::tie(k, ZA3, ZAcbc, ZA4, dtZc, Zc, dtZA, ZA);

  std::vector<double> dest(p_grid_size, 0.);

  // --- timing --------------------------------------------------------------------------------
  auto time_flow = [&](const char *name, auto &integrator) {
    std::vector<double> ms;
    ms.reserve(reps);
    for (int r = 0; r < warmup + reps; ++r) {
      Kokkos::fence();
      const auto t0 = std::chrono::steady_clock::now();
      integrator.map(dest.data(), coords, args);
      Kokkos::fence(); // without this we time the launch, not the execution
      const auto t1 = std::chrono::steady_clock::now();
      if (r >= warmup)
        ms.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    const Stat s = summarize(ms);
    // checksum makes the correctness gate possible: a variant must reproduce it
    double sum = 0., amax = 0.;
    for (double d : dest) {
      sum += d;
      amax = std::max(amax, std::abs(d));
    }
    std::printf("%-8s %-10s %9.4f %9.4f %9.4f   %+.15e %.6e\n", tag.c_str(), name, s.med, s.lo,
                s.hi, sum, amax);
    std::fflush(stdout);
  };

  // Global clock ramp. This laptop idles at ~210 MHz and needs order a SECOND of sustained load to
  // reach its 2445 MHz boost, so whichever measurement runs first comes out 20-40% slow. This has
  // to be a fixed wall-clock burn, not an "until it stops improving" heuristic: a single noisy
  // iteration satisfies such a heuristic immediately and it exits while the clock is still ramping.
  // Getting this wrong manufactures speedups out of nothing -- an unmatched-tile "1.3-1.5x" and a
  // LaunchBounds "1.14x" both evaporated once this burn was long enough.
  {
    const double burn_ms = std::atof(env_or("NTB_BURN_MS", "2500"));
    const auto t0 = std::chrono::steady_clock::now();
    double last = 0.;
    while (std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() <
           burn_ms) {
      const auto a = std::chrono::steady_clock::now();
      flows.ZA4.map(dest.data(), coords, args);
      Kokkos::fence();
      last = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count();
    }
    std::fprintf(stderr, "[ntb] clocks settled after %.0f ms burn, last ZA4 = %.3f ms\n", burn_ms,
                 last);
  }

  std::printf("%-8s %-10s %9s %9s %9s   %-22s %s\n", "tag", "flow", "med_ms", "min_ms", "max_ms",
              "checksum", "max_abs");
  note("starting timing");
  time_flow("ZA", flows.ZA);
  time_flow("Zc", flows.Zc);
  time_flow("ZA3", flows.ZA3);
  time_flow("ZAcbc", flows.ZAcbc);
  time_flow("ZA4", flows.ZA4);

  return 0;
}
