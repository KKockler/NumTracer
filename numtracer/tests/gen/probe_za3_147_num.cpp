#define KOKKOS_INLINE_FUNCTION inline
#define KOKKOS_FORCEINLINE_FUNCTION inline
#define KOKKOS_FUNCTION
#define __host__
#define __device__
#include <complex>
#include <cmath>
#include <random>
#include <cstdio>
#include <cstring>
#include "ZA3_147_num_kernels.hh"
template<int N> static inline double powr(double x){ double r=1.0; int n=N<0?-N:N; for(int i=0;i<n;++i) r*=x; return N<0?1.0/r:r; }
using std::pow; using std::sqrt; using std::sin; using std::cos; using std::tan; using std::exp; using std::log; using std::fma; using std::fabs;
static inline std::complex<double> fma(const std::complex<double>&a,const std::complex<double>&b,const std::complex<double>&c){return a*b+c;}
template<class T> using complex = std::complex<T>;
static inline double ntStub(double seed, double x){ double h = std::sin(seed*0.1031 + x*0.3127 + 1.7)*43758.5453; return 0.4 + 0.5*(h - std::floor(h)); }
static inline auto probe_full(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k)
{
  const double cosl1p1 = cos1;
  const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
  const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
  double fenv[(DiFfRG::za3_147_num::nenv) > 0 ? (DiFfRG::za3_147_num::nenv) : 1];
  const double dr_0 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_1 = ntStub(38926., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_2 = ntStub(14883., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_3 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_4 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_5 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_6 = ntStub(38926., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_7 = ntStub(14883., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_8 = ntStub(30415., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_9 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_10 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_11 = ntStub(38926., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_12 = ntStub(14883., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_13 = ntStub(30415., 1. * l1);
  const double dr_14 = -powr<-1>(l1) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1)) * ntStub(85453., 1. * k) - ntStub(85453., 1. * l1);
  const double dr_15 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_16 = ntStub(38926., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_17 = ntStub(14883., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_18 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_19 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_20 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_21 = ntStub(38926., 0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_22 = ntStub(14883., 0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_23 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_24 = ntStub(38926., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_25 = ntStub(14883., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_26 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_27 = ntStub(38926., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_28 = ntStub(14883., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_29 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_30 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  DiFfRG::za3_147_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20, dr_21, dr_22, dr_23, dr_24, dr_25, dr_26, dr_27, dr_28, dr_29, dr_30);
  const auto _interp21 = ntStub(62323., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp22 = ntStub(29441., 1. * l1);
  const auto _interp23 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp24 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp25 = ntStub(13492., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp26 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp27 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp28 = ntStub(29441., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp29 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
  const auto _interp31 = ntStub(62323., 0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp32 = ntStub(29441., 1. * sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp33 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _interp35 = ntStub(62323., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp36 = ntStub(29441., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp37 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
  const auto _interp43 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(fma(_interp23, _interp24, fma(_interp22, powr<2>(l1), 0.)));
  const auto _den4 = powr<-1>(fma(_interp23, _interp29, fma(_interp28, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
  const auto _den5 = powr<-1>(fma(_interp23, _interp33, fma(_interp32, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _den11 = powr<-1>(fma(_interp23, _interp37, fma(_interp36, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
  const auto _interp20 = ntStub(13076., 0.5 * sqrt(2. * powr<2>(l1) - 2. * cosl1p1 * l1 * p + 3. * powr<2>(p)));
  const auto _interp30 = ntStub(13076., 0.5 * sqrt(2. * powr<2>(l1) - 2. * cosl1p2 * l1 * p + 3. * powr<2>(p)));
  const auto _interp34 = ntStub(13076., 0.5 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp44 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp45 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp46 = ntStub(27736., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (2. * cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp47 = ntStub(20621., 1. * l1);
  const auto _interp48 = ntStub(20621., 1. * k);
  const auto _interp49 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp50 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp51 = ntStub(20621., 1.02 * k);
  const auto _interp52 = ntStub(95661., 1. * k);
  const auto _interp53 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp54 = ntStub(27736., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp55 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp56 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp57 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 1.5 * powr<2>(p)));
  const auto _interp58 = ntStub(20621., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp1 = ntStub(27191., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp2 = ntStub(85453., 1. * k);
  const auto _interp3 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp4 = ntStub(67246.00000000001, 1. * k);
  const auto _interp5 = ntStub(85453., 1.02 * k);
  const auto _interp6 = ntStub(30415., 1. * l1);
  const auto _interp7 = ntStub(85453., 1. * l1);
  const auto _interp8 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp9 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp10 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp11 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp12 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp13 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp14 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp15 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp16 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp17 = ntStub(30415., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp18 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp19 = ntStub(85453., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp38 = ntStub(62323., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (2. * cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp39 = ntStub(62323., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp40 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp41 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 1.5 * powr<2>(p)));
  const auto _interp42 = ntStub(29441., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _den13 = powr<-1>(powr<2>(_interp8) + powr<2>(_interp2 * _interp9 + _interp10 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
  const auto _den14 = powr<-1>(powr<2>(_interp14) + powr<2>(_interp15 * _interp2 + _interp16 * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
  const auto _den15 = powr<-1>(powr<2>(_interp17) + powr<2>(_interp18 * _interp2 + _interp19 * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
  const auto _den16 = powr<-1>(powr<2>(_interp11) + powr<2>(_interp13) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + _interp12 * _interp2 * (_interp12 * _interp2 + 2. * _interp13 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))));
  // clang-format off
  using _T = decltype(_den1 + _den11 + _den13 + _den14 + _den15 + _den16 + _den2 + _den4 + _den5 + DiFfRG::za3_147_num::tr0(fenv) + DiFfRG::za3_147_num::tr10(fenv) + DiFfRG::za3_147_num::tr11(fenv) + DiFfRG::za3_147_num::tr12(fenv) + DiFfRG::za3_147_num::tr13(fenv) + DiFfRG::za3_147_num::tr14(fenv) + DiFfRG::za3_147_num::tr15(fenv) + DiFfRG::za3_147_num::tr16(fenv) + DiFfRG::za3_147_num::tr17(fenv) + DiFfRG::za3_147_num::tr1(fenv) + DiFfRG::za3_147_num::tr2(fenv) + DiFfRG::za3_147_num::tr3(fenv) + DiFfRG::za3_147_num::tr4(fenv) + DiFfRG::za3_147_num::tr5(fenv) + DiFfRG::za3_147_num::tr6(fenv) + DiFfRG::za3_147_num::tr7(fenv) + DiFfRG::za3_147_num::tr8(fenv) + DiFfRG::za3_147_num::tr9(fenv) + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp6 + _interp7 + _interp8 + _interp9 + cosl1p1 + cosl1p2 + k + l1 + p);
  // clang-format on
  _T _acc{};
  { // subkernel 1
    const auto _den3 = powr<-2>(_interp24 * _interp48 + _interp47 * powr<2>(l1));
    const auto _den7 = -powr<-1>(_interp29 * _interp48 + _interp49 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _den8 = -powr<-1>(_interp33 * _interp48 + _interp55 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
    const auto _den9 = -powr<-1>(_interp43 * _interp48 + _interp58 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
    const auto _den12 = -powr<-1>(_interp37 * _interp48 + _interp50 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
    const auto _cse1 = powr<-2>(p);
    const auto _cse2 = -_interp23;
    const auto _cse3 = _cse2 + _interp27;
    const auto _cse4 = 50. * _cse3 * _den1 * powr<6>(k);
    const auto _cse5 = _cse4 + _interp26;
    const auto _cse6 = _cse5 * _interp24;
    const auto _cse7 = _interp23 * _interp25;
    const auto _cse8 = _cse6 + _cse7;
    const auto _cse9 = -_interp48;
    const auto _cse10 = _cse9 + _interp51;
    const auto _cse11 = 50. * _cse10;
    const auto _cse12 = _cse11 + _interp52;
    const auto _cse13 = _cse12 * _interp24;
    const auto _cse14 = _interp25 * _interp48;
    const auto _cse15 = _cse13 + _cse14;
    // clang-format off
    _acc += 0.001683501683501684 * fma(-1., _cse1 * _cse8 * _den2 * _den4 * DiFfRG::za3_147_num::tr3(fenv) * _interp20 * _interp21, fma(-1., _cse1 * _cse8 * _den2 * _den5 * DiFfRG::za3_147_num::tr4(fenv) * _interp30 * _interp31, fma(-1., _cse1 * _cse8 * _den11 * _den2 * DiFfRG::za3_147_num::tr5(fenv) * _interp34 * _interp35, fma(-1., _cse1 * _cse15 * _den12 * _den3 * _den7 * DiFfRG::za3_147_num::tr8(fenv) * _interp44 * _interp45 * _interp46, fma(_cse1, _cse15 * _den12 * _den3 * _den7 * DiFfRG::za3_147_num::tr9(fenv) * _interp44 * _interp45 * _interp46, fma(-1., _cse1 * _cse15 * _den12 * _den3 * _den8 * DiFfRG::za3_147_num::tr10(fenv) * _interp45 * _interp53 * _interp54, fma(_cse1, _cse15 * _den12 * _den3 * _den8 * DiFfRG::za3_147_num::tr11(fenv) * _interp45 * _interp53 * _interp54, fma(_cse1, _cse15 * _den3 * _den7 * _den9 * DiFfRG::za3_147_num::tr6(fenv) * _interp44 * _interp56 * _interp57, fma(-1., _cse1 * _cse15 * _den3 * _den7 * _den9 * DiFfRG::za3_147_num::tr7(fenv) * _interp44 * _interp56 * _interp57, 0.)))))))));
    // clang-format on
  }
  { // subkernel 2
    const auto _den6 = powr<-1>(_interp23 * _interp43 + _interp42 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
    const auto _den10 = powr<-2>(powr<2>(_interp6) + powr<2>(_interp2 * _interp3 + _interp7 * l1));
    const auto _cse1 = powr<-1>(l1);
    const auto _cse2 = powr<-2>(p);
    const auto _cse3 = -_cse1 * _interp1 * _interp2;
    const auto _cse4 = -_interp2;
    const auto _cse5 = _cse4 + _interp5;
    const auto _cse6 = 50. * _cse5;
    const auto _cse7 = _cse6 + _interp4;
    const auto _cse8 = -_cse1 * _cse7 * _interp3;
    const auto _cse9 = _cse3 + _cse8;
    const auto _cse10 = -_interp23;
    const auto _cse11 = _cse10 + _interp27;
    const auto _cse12 = 50. * _cse11 * _den1 * powr<6>(k);
    const auto _cse13 = _cse12 + _interp26;
    // clang-format off
    _acc += 0.001683501683501684 * fma(-2., _cse2 * _cse9 * _den10 * _den13 * _den15 * DiFfRG::za3_147_num::tr12(fenv), fma(2., _cse2 * _cse9 * _den10 * _den13 * _den15 * DiFfRG::za3_147_num::tr13(fenv), fma(2., _cse2 * _cse9 * _den10 * _den13 * _den16 * DiFfRG::za3_147_num::tr14(fenv), fma(-2., _cse2 * _cse9 * _den10 * _den13 * _den16 * DiFfRG::za3_147_num::tr15(fenv), fma(2., _cse2 * _cse9 * _den10 * _den14 * _den16 * DiFfRG::za3_147_num::tr16(fenv), fma(-2., _cse2 * _cse9 * _den10 * _den14 * _den16 * DiFfRG::za3_147_num::tr17(fenv), fma(-1., _cse13 * _cse2 * _den11 * _den2 * _den4 * DiFfRG::za3_147_num::tr2(fenv) * _interp21 * _interp24 * _interp35 * _interp38, fma(-1., _cse2 * _den11 * _den2 * _den4 * DiFfRG::za3_147_num::tr2(fenv) * _interp21 * _interp23 * _interp25 * _interp35 * _interp38, fma(-1., _cse13 * _cse2 * _den11 * _den2 * _den5 * DiFfRG::za3_147_num::tr1(fenv) * _interp24 * _interp31 * _interp35 * _interp39, fma(-1., _cse2 * _den11 * _den2 * _den5 * DiFfRG::za3_147_num::tr1(fenv) * _interp23 * _interp25 * _interp31 * _interp35 * _interp39, fma(-1., _cse13 * _cse2 * _den2 * _den4 * _den6 * DiFfRG::za3_147_num::tr0(fenv) * _interp21 * _interp24 * _interp40 * _interp41, fma(-1., _cse2 * _den2 * _den4 * _den6 * DiFfRG::za3_147_num::tr0(fenv) * _interp21 * _interp23 * _interp25 * _interp40 * _interp41, 0.))))))))))));
    // clang-format on
  }
  return _acc;
}
static inline auto probe_proj(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k)
{
  const double cosl1p1 = cos1;
  const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
  const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
  double fenv[(DiFfRG::za3_147_num::nenv) > 0 ? (DiFfRG::za3_147_num::nenv) : 1];
  const double dr_0 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_1 = ntStub(38926., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_2 = ntStub(14883., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_3 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_4 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_5 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_6 = ntStub(38926., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_7 = ntStub(14883., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_8 = ntStub(30415., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_9 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_10 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_11 = ntStub(38926., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_12 = ntStub(14883., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_13 = ntStub(30415., 1. * l1);
  const double dr_14 = -powr<-1>(l1) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1)) * ntStub(85453., 1. * k) - ntStub(85453., 1. * l1);
  const double dr_15 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_16 = ntStub(38926., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_17 = ntStub(14883., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_18 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_19 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_20 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_21 = ntStub(38926., 0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_22 = ntStub(14883., 0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_23 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_24 = ntStub(38926., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_25 = ntStub(14883., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_26 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_27 = ntStub(38926., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_28 = ntStub(14883., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
  const double dr_29 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const double dr_30 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  DiFfRG::za3_147_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20, dr_21, dr_22, dr_23, dr_24, dr_25, dr_26, dr_27, dr_28, dr_29, dr_30);
  const auto _interp21 = ntStub(62323., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp22 = ntStub(29441., 1. * l1);
  const auto _interp23 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp24 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp25 = ntStub(13492., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp26 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp27 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp28 = ntStub(29441., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp29 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
  const auto _interp31 = ntStub(62323., 0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp32 = ntStub(29441., 1. * sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp33 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _interp35 = ntStub(62323., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp36 = ntStub(29441., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
  const auto _interp37 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
  const auto _interp43 = ntStub(89189., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(fma(_interp23, _interp24, fma(_interp22, powr<2>(l1), 0.)));
  const auto _den4 = powr<-1>(fma(_interp23, _interp29, fma(_interp28, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
  const auto _den5 = powr<-1>(fma(_interp23, _interp33, fma(_interp32, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
  const auto _den11 = powr<-1>(fma(_interp23, _interp37, fma(_interp36, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
  const auto _interp20 = ntStub(13076., 0.5 * sqrt(2. * powr<2>(l1) - 2. * cosl1p1 * l1 * p + 3. * powr<2>(p)));
  const auto _interp30 = ntStub(13076., 0.5 * sqrt(2. * powr<2>(l1) - 2. * cosl1p2 * l1 * p + 3. * powr<2>(p)));
  const auto _interp34 = ntStub(13076., 0.5 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp44 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp45 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp46 = ntStub(27736., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (2. * cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp47 = ntStub(20621., 1. * l1);
  const auto _interp48 = ntStub(20621., 1. * k);
  const auto _interp49 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp50 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp51 = ntStub(20621., 1.02 * k);
  const auto _interp52 = ntStub(95661., 1. * k);
  const auto _interp53 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp54 = ntStub(27736., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp55 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp56 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp57 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 1.5 * powr<2>(p)));
  const auto _interp58 = ntStub(20621., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp1 = ntStub(27191., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp2 = ntStub(85453., 1. * k);
  const auto _interp3 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp4 = ntStub(67246.00000000001, 1. * k);
  const auto _interp5 = ntStub(85453., 1.02 * k);
  const auto _interp6 = ntStub(30415., 1. * l1);
  const auto _interp7 = ntStub(85453., 1. * l1);
  const auto _interp8 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp9 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp10 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp11 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp12 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp13 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
  const auto _interp14 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp15 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp16 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp17 = ntStub(30415., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp18 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp19 = ntStub(85453., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp38 = ntStub(62323., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (2. * cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp39 = ntStub(62323., 0.5773502691896258 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 3. * powr<2>(p)));
  const auto _interp40 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp41 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 1.5 * powr<2>(p)));
  const auto _interp42 = ntStub(29441., 1. * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _den13 = powr<-1>(powr<2>(_interp8) + powr<2>(_interp2 * _interp9 + _interp10 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
  const auto _den14 = powr<-1>(powr<2>(_interp14) + powr<2>(_interp15 * _interp2 + _interp16 * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
  const auto _den15 = powr<-1>(powr<2>(_interp17) + powr<2>(_interp18 * _interp2 + _interp19 * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
  const auto _den16 = powr<-1>(powr<2>(_interp11) + powr<2>(_interp13) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + _interp12 * _interp2 * (_interp12 * _interp2 + 2. * _interp13 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))));
  // clang-format off
  using _T = decltype(_den1 + _den11 + _den13 + _den14 + _den15 + _den16 + _den2 + _den4 + _den5 + DiFfRG::za3_147_num::tr0(fenv) + DiFfRG::za3_147_num::tr10(fenv) + DiFfRG::za3_147_num::tr11(fenv) + DiFfRG::za3_147_num::tr12(fenv) + DiFfRG::za3_147_num::tr13(fenv) + DiFfRG::za3_147_num::tr14(fenv) + DiFfRG::za3_147_num::tr15(fenv) + DiFfRG::za3_147_num::tr16(fenv) + DiFfRG::za3_147_num::tr17(fenv) + DiFfRG::za3_147_num::tr1(fenv) + DiFfRG::za3_147_num::tr2(fenv) + DiFfRG::za3_147_num::tr3(fenv) + DiFfRG::za3_147_num::tr4(fenv) + DiFfRG::za3_147_num::tr5(fenv) + DiFfRG::za3_147_num::tr6(fenv) + DiFfRG::za3_147_num::tr7(fenv) + DiFfRG::za3_147_num::tr8(fenv) + DiFfRG::za3_147_num::tr9(fenv) + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp6 + _interp7 + _interp8 + _interp9 + cosl1p1 + cosl1p2 + k + l1 + p);
  // clang-format on
  _T _acc{};
  { // subkernel 1
    const auto _den3 = powr<-2>(_interp24 * _interp48 + _interp47 * powr<2>(l1));
    const auto _den7 = -powr<-1>(_interp29 * _interp48 + _interp49 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _den8 = -powr<-1>(_interp33 * _interp48 + _interp55 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
    const auto _den9 = -powr<-1>(_interp43 * _interp48 + _interp58 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
    const auto _den12 = -powr<-1>(_interp37 * _interp48 + _interp50 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
    const auto _cse1 = powr<-2>(p);
    const auto _cse2 = -_interp23;
    const auto _cse3 = _cse2 + _interp27;
    const auto _cse4 = 50. * _cse3 * _den1 * powr<6>(k);
    const auto _cse5 = _cse4 + _interp26;
    const auto _cse6 = _cse5 * _interp24;
    const auto _cse7 = _interp23 * _interp25;
    const auto _cse8 = _cse6 + _cse7;
    const auto _cse9 = -_interp48;
    const auto _cse10 = _cse9 + _interp51;
    const auto _cse11 = 50. * _cse10;
    const auto _cse12 = _cse11 + _interp52;
    const auto _cse13 = _cse12 * _interp24;
    const auto _cse14 = _interp25 * _interp48;
    const auto _cse15 = _cse13 + _cse14;
    // clang-format off
    _acc += 0.001683501683501684 * fma(-1., _cse1 * _cse8 * _den2 * _den4 * DiFfRG::za3_147_num::tr3(fenv) * _interp20 * _interp21, fma(-1., _cse1 * _cse8 * _den2 * _den5 * DiFfRG::za3_147_num::tr4(fenv) * _interp30 * _interp31, fma(-1., _cse1 * _cse8 * _den11 * _den2 * DiFfRG::za3_147_num::tr5(fenv) * _interp34 * _interp35, fma(-1., _cse1 * _cse15 * _den12 * _den3 * _den7 * DiFfRG::za3_147_num::tr8(fenv) * _interp44 * _interp45 * _interp46, fma(_cse1, _cse15 * _den12 * _den3 * _den7 * DiFfRG::za3_147_num::tr9(fenv) * _interp44 * _interp45 * _interp46, fma(-1., _cse1 * _cse15 * _den12 * _den3 * _den8 * DiFfRG::za3_147_num::tr10(fenv) * _interp45 * _interp53 * _interp54, fma(_cse1, _cse15 * _den12 * _den3 * _den8 * DiFfRG::za3_147_num::tr11(fenv) * _interp45 * _interp53 * _interp54, fma(_cse1, _cse15 * _den3 * _den7 * _den9 * DiFfRG::za3_147_num::tr6(fenv) * _interp44 * _interp56 * _interp57, fma(-1., _cse1 * _cse15 * _den3 * _den7 * _den9 * DiFfRG::za3_147_num::tr7(fenv) * _interp44 * _interp56 * _interp57, 0.)))))))));
    // clang-format on
  }
  { // subkernel 2
    const auto _den6 = powr<-1>(_interp23 * _interp43 + _interp42 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
    const auto _den10 = powr<-2>(powr<2>(_interp6) + powr<2>(_interp2 * _interp3 + _interp7 * l1));
    const auto _cse1 = powr<-1>(l1);
    const auto _cse2 = powr<-2>(p);
    const auto _cse3 = -_cse1 * _interp1 * _interp2;
    const auto _cse4 = -_interp2;
    const auto _cse5 = _cse4 + _interp5;
    const auto _cse6 = 50. * _cse5;
    const auto _cse7 = _cse6 + _interp4;
    const auto _cse8 = -_cse1 * _cse7 * _interp3;
    const auto _cse9 = _cse3 + _cse8;
    const auto _cse10 = -_interp23;
    const auto _cse11 = _cse10 + _interp27;
    const auto _cse12 = 50. * _cse11 * _den1 * powr<6>(k);
    const auto _cse13 = _cse12 + _interp26;
    // clang-format off
    _acc += 0.001683501683501684 * fma(-2., _cse2 * _cse9 * _den10 * _den13 * _den15 * DiFfRG::za3_147_num::tr12(fenv), fma(2., _cse2 * _cse9 * _den10 * _den13 * _den15 * DiFfRG::za3_147_num::tr13(fenv), fma(2., _cse2 * _cse9 * _den10 * _den13 * _den16 * DiFfRG::za3_147_num::tr14(fenv), fma(-2., _cse2 * _cse9 * _den10 * _den13 * _den16 * DiFfRG::za3_147_num::tr15(fenv), fma(2., _cse2 * _cse9 * _den10 * _den14 * _den16 * DiFfRG::za3_147_num::tr16(fenv), fma(-2., _cse2 * _cse9 * _den10 * _den14 * _den16 * DiFfRG::za3_147_num::tr17(fenv), fma(-1., _cse13 * _cse2 * _den11 * _den2 * _den4 * DiFfRG::za3_147_num::tr2(fenv) * _interp21 * _interp24 * _interp35 * _interp38, fma(-1., _cse2 * _den11 * _den2 * _den4 * DiFfRG::za3_147_num::tr2(fenv) * _interp21 * _interp23 * _interp25 * _interp35 * _interp38, fma(-1., _cse13 * _cse2 * _den11 * _den2 * _den5 * DiFfRG::za3_147_num::tr1(fenv) * _interp24 * _interp31 * _interp35 * _interp39, fma(-1., _cse2 * _den11 * _den2 * _den5 * DiFfRG::za3_147_num::tr1(fenv) * _interp23 * _interp25 * _interp31 * _interp35 * _interp39, fma(-1., _cse13 * _cse2 * _den2 * _den4 * _den6 * DiFfRG::za3_147_num::tr0(fenv) * _interp21 * _interp24 * _interp40 * _interp41, fma(-1., _cse2 * _den2 * _den4 * _den6 * DiFfRG::za3_147_num::tr0(fenv) * _interp21 * _interp23 * _interp25 * _interp40 * _interp41, 0.))))))))))));
    // clang-format on
  }
  return _acc;
}
int main(int argc, char** argv){
  const char* outf=nullptr; const char* macro=nullptr;
  for(int i=1;i<argc;++i){ if(!std::strcmp(argv[i],"-o") && i+1<argc) outf=argv[++i];
                           else if(!std::strcmp(argv[i],"-m") && i+1<argc) macro=argv[++i]; }
  std::mt19937_64 rng(12345); std::uniform_real_distribution<double> U(0.25,3.0),Uc(-0.9,0.9),Uph(0.1,6.2);
  double mim=0,mdiff=0,mre=0,mrim=0,mrdiff=0; long ok=0;
  for(int n=0;n<4000;++n){ double l1 = U(rng); double cos1 = Uc(rng); double cos2 = Uc(rng); double p = U(rng); double k = U(rng);
    std::complex<double> f = probe_full(l1, cos1, cos2, p, k); std::complex<double> pj = probe_proj(l1, cos1, cos2, p, k);
    double im=std::imag(f), re=std::real(f), df=std::abs(f-pj);
    if(std::isfinite(im)&&std::isfinite(re)&&std::isfinite(df)){ mim=std::max(mim,std::fabs(im)); mdiff=std::max(mdiff,df); mre=std::max(mre,std::fabs(re));
      mrim=std::max(mrim, std::fabs(im)/(std::fabs(re)+1.0)); mrdiff=std::max(mrdiff, df/(std::abs(f)+1.0)); ++ok; } }
  if(ok < 1){ std::fprintf(stderr, "[probe] no usable points\n"); return 2; }
  const int verdict = (mrim > 1.e-9) ? 0 : ((mrdiff <= 1.e-9) ? 2 : 1);
  std::printf("%.10e %.10e %.10e %.10e %.10e %ld %d\n", mim, mdiff, mre, mrim, mrdiff, ok, verdict);
  if(outf && macro){ std::FILE* f = std::fopen(outf, "w");
    if(!f){ std::fprintf(stderr, "[probe] cannot write %s\n", outf); return 3; }
    std::fprintf(f, "// GENERATED by the numtrace step — do not edit.\n");
    std::fprintf(f, "// 2 = Pure (imaginary coefficients dropped), 1 = RePart (re/im split), 0 = complex.\n");
    std::fprintf(f, "#pragma once\n#define %s %d\n", macro, verdict);
    std::fclose(f); }
  return 0; }
