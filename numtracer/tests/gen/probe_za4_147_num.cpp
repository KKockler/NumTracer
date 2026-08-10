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
#include "ZA4_147_num_kernels.hh"
template<int N> static inline double powr(double x){ double r=1.0; int n=N<0?-N:N; for(int i=0;i<n;++i) r*=x; return N<0?1.0/r:r; }
using std::pow; using std::sqrt; using std::sin; using std::cos; using std::tan; using std::exp; using std::log; using std::fma; using std::fabs;
static inline std::complex<double> fma(const std::complex<double>&a,const std::complex<double>&b,const std::complex<double>&c){return a*b+c;}
template<class T> using complex = std::complex<T>;
static inline double ntStub(double seed, double x){ double h = std::sin(seed*0.1031 + x*0.3127 + 1.7)*43758.5453; return 0.4 + 0.5*(h - std::floor(h)); }
static inline auto probe_full(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k)
{
  const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
  const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
  const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
  const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
  double fenv[(DiFfRG::za4_147_num::nenv) > 0 ? (DiFfRG::za4_147_num::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_1 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_2 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_3 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_4 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_5 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_6 = ntStub(30415., 1. * l1);
  const double dr_7 = -powr<-1>(l1) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1)) * ntStub(85453., 1. * k) - ntStub(85453., 1. * l1);
  DiFfRG::za4_147_num::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7);
  const auto _interp1 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp2 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
  const auto _interp3 = ntStub(12206., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
  const auto _interp4 = ntStub(12206., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
  const auto _interp5 = ntStub(27191., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp6 = ntStub(85453., 1. * k);
  const auto _interp7 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp8 = ntStub(67246.00000000001, 1. * k);
  const auto _interp9 = ntStub(85453., 1.02 * k);
  const auto _interp10 = ntStub(30415., 1. * l1);
  const auto _interp11 = ntStub(85453., 1. * l1);
  const auto _interp12 = ntStub(30415., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp13 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
  const auto _interp14 = ntStub(85453., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp15 = ntStub(30415., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
  const auto _interp16 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p), 0.)));
  const auto _interp17 = ntStub(85453., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
  const auto _interp18 = ntStub(30415., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
  const auto _interp19 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
  const auto _interp20 = ntStub(85453., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
  const auto _interp21 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp22 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
  const auto _interp23 = ntStub(14883., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
  const auto _interp24 = ntStub(14883., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
  const auto _interp25 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp26 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
  const auto _interp27 = ntStub(38926., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
  const auto _interp28 = ntStub(38926., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
  const auto _den7 = powr<-2>(fma(powr<-2>(l1) * powr<2>(_interp6 * _interp7 + _interp11 * l1), powr<2>(l1), powr<2>(_interp10)));
  const auto _den13 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), powr<2>(_interp14 + _interp13 * _interp6 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)))), powr<2>(_interp12)));
  const auto _den14 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p), powr<2>(_interp20 + _interp19 * _interp6 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)))), powr<2>(_interp18)));
  const auto _den15 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p), powr<2>(_interp17 + _interp16 * _interp6 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)))), powr<2>(_interp15)));
  const auto _interp29 = ntStub(13076., 0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
  const auto _interp30 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp31 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp32 = ntStub(29441., 1. * l1);
  const auto _interp33 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp34 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp35 = ntStub(13492., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp36 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp37 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp38 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp39 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp40 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp41 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp42 = ntStub(13076., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp43 = ntStub(62323., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp44 = ntStub(29441., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp45 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p1 - 3.23606797749979 * cosl1p2) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp46 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp47 = ntStub(62323., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp48 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp49 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp50 = ntStub(13076., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp51 = ntStub(29441., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
  const auto _interp52 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p2 - 3.23606797749979 * cosl1p3) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp53 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp54 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp55 = ntStub(27736., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp56 = ntStub(27736., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp57 = ntStub(20621., 1. * l1);
  const auto _interp58 = ntStub(20621., 1. * k);
  const auto _interp59 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp60 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp61 = ntStub(20621., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp62 = ntStub(20621., 1.02 * k);
  const auto _interp63 = ntStub(95661., 1. * k);
  const auto _den8 = 3. * powr<-1>(3. * _interp33 * _interp45 + _interp44 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _den9 = 3. * powr<-1>(3. * _interp33 * _interp52 + _interp51 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
  const auto _den11 = -3. * powr<-1>(3. * _interp45 * _interp58 + _interp61 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  // clang-format off
  using _T = decltype(complex<double>(0.,0.002551020408163265) + _den11 + _den13 + _den14 + _den15 + _den7 + _den8 + _den9 + DiFfRG::za4_147_num::tr0(fenv) + DiFfRG::za4_147_num::tr10(fenv) + DiFfRG::za4_147_num::tr11(fenv) + DiFfRG::za4_147_num::tr12(fenv) + DiFfRG::za4_147_num::tr13(fenv) + DiFfRG::za4_147_num::tr14(fenv) + DiFfRG::za4_147_num::tr15(fenv) + DiFfRG::za4_147_num::tr16(fenv) + DiFfRG::za4_147_num::tr17(fenv) + DiFfRG::za4_147_num::tr18(fenv) + DiFfRG::za4_147_num::tr19(fenv) + DiFfRG::za4_147_num::tr1(fenv) + DiFfRG::za4_147_num::tr20(fenv) + DiFfRG::za4_147_num::tr21(fenv) + DiFfRG::za4_147_num::tr22(fenv) + DiFfRG::za4_147_num::tr23(fenv) + DiFfRG::za4_147_num::tr24(fenv) + DiFfRG::za4_147_num::tr25(fenv) + DiFfRG::za4_147_num::tr26(fenv) + DiFfRG::za4_147_num::tr27(fenv) + DiFfRG::za4_147_num::tr28(fenv) + DiFfRG::za4_147_num::tr29(fenv) + DiFfRG::za4_147_num::tr2(fenv) + DiFfRG::za4_147_num::tr30(fenv) + DiFfRG::za4_147_num::tr31(fenv) + DiFfRG::za4_147_num::tr32(fenv) + DiFfRG::za4_147_num::tr33(fenv) + DiFfRG::za4_147_num::tr34(fenv) + DiFfRG::za4_147_num::tr35(fenv) + DiFfRG::za4_147_num::tr36(fenv) + DiFfRG::za4_147_num::tr37(fenv) + DiFfRG::za4_147_num::tr38(fenv) + DiFfRG::za4_147_num::tr39(fenv) + DiFfRG::za4_147_num::tr3(fenv) + DiFfRG::za4_147_num::tr40(fenv) + DiFfRG::za4_147_num::tr41(fenv) + DiFfRG::za4_147_num::tr42(fenv) + DiFfRG::za4_147_num::tr43(fenv) + DiFfRG::za4_147_num::tr44(fenv) + DiFfRG::za4_147_num::tr45(fenv) + DiFfRG::za4_147_num::tr46(fenv) + DiFfRG::za4_147_num::tr47(fenv) + DiFfRG::za4_147_num::tr48(fenv) + DiFfRG::za4_147_num::tr49(fenv) + DiFfRG::za4_147_num::tr4(fenv) + DiFfRG::za4_147_num::tr50(fenv) + DiFfRG::za4_147_num::tr51(fenv) + DiFfRG::za4_147_num::tr52(fenv) + DiFfRG::za4_147_num::tr53(fenv) + DiFfRG::za4_147_num::tr54(fenv) + DiFfRG::za4_147_num::tr55(fenv) + DiFfRG::za4_147_num::tr56(fenv) + DiFfRG::za4_147_num::tr57(fenv) + DiFfRG::za4_147_num::tr58(fenv) + DiFfRG::za4_147_num::tr59(fenv) + DiFfRG::za4_147_num::tr5(fenv) + DiFfRG::za4_147_num::tr60(fenv) + DiFfRG::za4_147_num::tr61(fenv) + DiFfRG::za4_147_num::tr62(fenv) + DiFfRG::za4_147_num::tr63(fenv) + DiFfRG::za4_147_num::tr64(fenv) + DiFfRG::za4_147_num::tr65(fenv) + DiFfRG::za4_147_num::tr66(fenv) + DiFfRG::za4_147_num::tr67(fenv) + DiFfRG::za4_147_num::tr68(fenv) + DiFfRG::za4_147_num::tr69(fenv) + DiFfRG::za4_147_num::tr6(fenv) + DiFfRG::za4_147_num::tr70(fenv) + DiFfRG::za4_147_num::tr71(fenv) + DiFfRG::za4_147_num::tr72(fenv) + DiFfRG::za4_147_num::tr73(fenv) + DiFfRG::za4_147_num::tr74(fenv) + DiFfRG::za4_147_num::tr75(fenv) + DiFfRG::za4_147_num::tr76(fenv) + DiFfRG::za4_147_num::tr77(fenv) + DiFfRG::za4_147_num::tr78(fenv) + DiFfRG::za4_147_num::tr79(fenv) + DiFfRG::za4_147_num::tr7(fenv) + DiFfRG::za4_147_num::tr80(fenv) + DiFfRG::za4_147_num::tr81(fenv) + DiFfRG::za4_147_num::tr82(fenv) + DiFfRG::za4_147_num::tr83(fenv) + DiFfRG::za4_147_num::tr84(fenv) + DiFfRG::za4_147_num::tr85(fenv) + DiFfRG::za4_147_num::tr8(fenv) + DiFfRG::za4_147_num::tr9(fenv) + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp7 + _interp8 + _interp9 + cosl1p1 + cosl1p2 + cosl1p3 + k + l1 + p);
  // clang-format on
  _T _acc{};
  { // subkernel 1
    const auto _den1 = powr<-1>(1. + powr<6>(k));
    const auto _den2 = powr<-2>(_interp33 * _interp34 + _interp32 * powr<2>(l1));
    const auto _den3 = powr<-2>(_interp34 * _interp58 + _interp57 * powr<2>(l1));
    const auto _den4 = powr<-1>(_interp33 * _interp49 + _interp48 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _den5 = powr<-1>(_interp33 * _interp39 + _interp38 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
    const auto _den6 = -powr<-1>(_interp49 * _interp58 + _interp59 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _den10 = powr<-1>(_interp33 * _interp41 + _interp40 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
    const auto _den12 = -powr<-1>(_interp41 * _interp58 + _interp60 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
    // clang-format off
    _acc += fma(0.01020408163265306, _den11 * _den12 * _den3 * _den6 * DiFfRG::za4_147_num::tr4(fenv) * _interp53 * _interp54 * _interp55 * _interp56 * (_interp35 * _interp58 + _interp34 * (50. * (-_interp58 + _interp62) + _interp63)), fma(0.002551020408163265, _den10 * _den2 * _den5 * DiFfRG::za4_147_num::tr2(fenv) * _interp29 * _interp30 * _interp31 * (_interp33 * _interp35 + _interp34 * (_interp36 + 50. * _den1 * (-_interp33 + _interp37) * powr<6>(k))), fma(0.00510204081632653, _den10 * _den2 * _den8 * DiFfRG::za4_147_num::tr1(fenv) * _interp31 * _interp42 * _interp43 * (_interp33 * _interp35 + _interp34 * (_interp36 + 50. * _den1 * (-_interp33 + _interp37) * powr<6>(k))), fma(0.00510204081632653, _den10 * _den2 * _den4 * _den8 * DiFfRG::za4_147_num::tr0(fenv) * _interp31 * _interp43 * _interp46 * _interp47 * (_interp33 * _interp35 + _interp34 * (_interp36 + 50. * _den1 * (-_interp33 + _interp37) * powr<6>(k))), fma(0.001275510204081632, _den2 * _den9 * DiFfRG::za4_147_num::tr3(fenv) * powr<2>(_interp50) * (_interp33 * _interp35 + _interp34 * (_interp36 + 50. * _den1 * (-_interp33 + _interp37) * powr<6>(k))), fma(complex<double>(0.,-0.00510204081632653), _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr53(fenv) * _interp1 * _interp2 * _interp23 * _interp24 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), fma(complex<double>(0.,-0.002551020408163265), _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr57(fenv) * _interp2 * _interp21 * _interp23 * _interp24 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), fma(complex<double>(0.,-0.002551020408163265), _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr83(fenv) * _interp1 * _interp22 * _interp23 * _interp24 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), fma(complex<double>(0.,-0.001275510204081632), _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr85(fenv) * _interp21 * _interp22 * _interp23 * _interp24 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), fma(-0.00510204081632653, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr54(fenv) * _interp2 * _interp23 * _interp24 * _interp25 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), fma(-0.002551020408163265, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr84(fenv) * _interp22 * _interp23 * _interp24 * _interp25 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), fma(0.00510204081632653, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr55(fenv) * _interp1 * _interp23 * _interp24 * _interp26 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), fma(0.002551020408163265, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr58(fenv) * _interp21 * _interp23 * _interp24 * _interp26 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), 0.)))))))))))));
    // clang-format on
  }
  { // subkernel 2
    const auto _cse1_k2 = powr<-1>(l1);
    const auto _cse2_k2 = -_cse1_k2 * _interp5 * _interp6;
    const auto _cse3_k2 = -_interp6;
    const auto _cse4_k2 = _cse3_k2 + _interp9;
    const auto _cse5_k2 = 50. * _cse4_k2;
    const auto _cse6_k2 = _cse5_k2 + _interp8;
    const auto _cse7_k2 = -_cse1_k2 * _cse6_k2 * _interp7;
    const auto _cse8_k2 = _cse2_k2 + _cse7_k2;
    // clang-format off
    _acc += 0.002551020408163265 * fma(complex<double>(0.,-2.), _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr56(fenv) * _interp23 * _interp24 * _interp25 * _interp26, fma(-4., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr43(fenv) * _interp1 * _interp2 * _interp24 * _interp27, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr50(fenv) * _interp2 * _interp21 * _interp24 * _interp27, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr79(fenv) * _interp1 * _interp22 * _interp24 * _interp27, fma(-1., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr82(fenv) * _interp21 * _interp22 * _interp24 * _interp27, fma(complex<double>(0.,4.), _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr44(fenv) * _interp2 * _interp24 * _interp25 * _interp27, fma(complex<double>(0.,2.), _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr80(fenv) * _interp22 * _interp24 * _interp25 * _interp27, fma(complex<double>(0.,-4.), _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr47(fenv) * _interp1 * _interp24 * _interp26 * _interp27, fma(complex<double>(0.,-2.), _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr52(fenv) * _interp21 * _interp24 * _interp26 * _interp27, fma(-4., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr48(fenv) * _interp24 * _interp25 * _interp26 * _interp27, fma(-4., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr31(fenv) * _interp1 * _interp2 * _interp23 * _interp28, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr38(fenv) * _interp2 * _interp21 * _interp23 * _interp28, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr73(fenv) * _interp1 * _interp22 * _interp23 * _interp28, 0.)))))))))))));
    // clang-format on
  }
  { // subkernel 3
    const auto _cse1_k3 = powr<-1>(l1);
    const auto _cse2_k3 = -_cse1_k3 * _interp5 * _interp6;
    const auto _cse3_k3 = -_interp6;
    const auto _cse4_k3 = _cse3_k3 + _interp9;
    const auto _cse5_k3 = 50. * _cse4_k3;
    const auto _cse6_k3 = _cse5_k3 + _interp8;
    const auto _cse7_k3 = -_cse1_k3 * _cse6_k3 * _interp7;
    const auto _cse8_k3 = _cse2_k3 + _cse7_k3;
    // clang-format off
    _acc += complex<double>(0.,0.002551020408163265) * fma(complex<double>(0.,1.), _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr76(fenv) * _interp21 * _interp22 * _interp23 * _interp28, fma(4., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr32(fenv) * _interp2 * _interp23 * _interp25 * _interp28, fma(2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr74(fenv) * _interp22 * _interp23 * _interp25 * _interp28, fma(-4., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr35(fenv) * _interp1 * _interp23 * _interp26 * _interp28, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr40(fenv) * _interp21 * _interp23 * _interp26 * _interp28, fma(complex<double>(0.,4.), _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr36(fenv) * _interp23 * _interp25 * _interp26 * _interp28, fma(8., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr11(fenv) * _interp1 * _interp2 * _interp27 * _interp28, fma(4., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr24(fenv) * _interp2 * _interp21 * _interp27 * _interp28, fma(4., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr65(fenv) * _interp1 * _interp22 * _interp27 * _interp28, fma(2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr70(fenv) * _interp21 * _interp22 * _interp27 * _interp28, fma(complex<double>(0.,-8.), _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr12(fenv) * _interp2 * _interp25 * _interp27 * _interp28, fma(complex<double>(0.,-4.), _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr66(fenv) * _interp22 * _interp25 * _interp27 * _interp28, fma(complex<double>(0.,8.), _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr19(fenv) * _interp1 * _interp26 * _interp27 * _interp28, 0.)))))))))))));
    // clang-format on
  }
  { // subkernel 4
    const auto _cse1_k4 = powr<-1>(l1);
    const auto _cse2_k4 = -_cse1_k4 * _interp5 * _interp6;
    const auto _cse3_k4 = -_interp6;
    const auto _cse4_k4 = _cse3_k4 + _interp9;
    const auto _cse5_k4 = 50. * _cse4_k4;
    const auto _cse6_k4 = _cse5_k4 + _interp8;
    const auto _cse7_k4 = -_cse1_k4 * _cse6_k4 * _interp7;
    const auto _cse8_k4 = _cse2_k4 + _cse7_k4;
    // clang-format off
    _acc += complex<double>(0.,0.002551020408163265) * fma(complex<double>(0.,4.), _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr28(fenv) * _interp21 * _interp26 * _interp27 * _interp28, fma(8., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr20(fenv) * _interp25 * _interp26 * _interp27 * _interp28, fma(-4., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr41(fenv) * _interp1 * _interp2 * _interp24 * _interp3, fma(-2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr49(fenv) * _interp2 * _interp21 * _interp24 * _interp3, fma(-2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr77(fenv) * _interp1 * _interp22 * _interp24 * _interp3, fma(-1., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr81(fenv) * _interp21 * _interp22 * _interp24 * _interp3, fma(complex<double>(0.,4.), _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr42(fenv) * _interp2 * _interp24 * _interp25 * _interp3, fma(complex<double>(0.,2.), _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr78(fenv) * _interp22 * _interp24 * _interp25 * _interp3, fma(complex<double>(0.,-4.), _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr45(fenv) * _interp1 * _interp24 * _interp26 * _interp3, fma(complex<double>(0.,-2.), _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr51(fenv) * _interp21 * _interp24 * _interp26 * _interp3, fma(-4., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr46(fenv) * _interp24 * _interp25 * _interp26 * _interp3, fma(complex<double>(0.,8.), _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr9(fenv) * _interp1 * _interp2 * _interp28 * _interp3, fma(complex<double>(0.,4.), _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr23(fenv) * _interp2 * _interp21 * _interp28 * _interp3, 0.)))))))))))));
    // clang-format on
  }
  { // subkernel 5
    const auto _cse1_k5 = powr<-1>(l1);
    const auto _cse2_k5 = -_cse1_k5 * _interp5 * _interp6;
    const auto _cse3_k5 = -_interp6;
    const auto _cse4_k5 = _cse3_k5 + _interp9;
    const auto _cse5_k5 = 50. * _cse4_k5;
    const auto _cse6_k5 = _cse5_k5 + _interp8;
    const auto _cse7_k5 = -_cse1_k5 * _cse6_k5 * _interp7;
    const auto _cse8_k5 = _cse2_k5 + _cse7_k5;
    // clang-format off
    _acc += complex<double>(0.,0.002551020408163265) * fma(complex<double>(0.,4.), _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr63(fenv) * _interp1 * _interp22 * _interp28 * _interp3, fma(complex<double>(0.,2.), _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr69(fenv) * _interp21 * _interp22 * _interp28 * _interp3, fma(8., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr10(fenv) * _interp2 * _interp25 * _interp28 * _interp3, fma(4., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr64(fenv) * _interp22 * _interp25 * _interp28 * _interp3, fma(-8., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr17(fenv) * _interp1 * _interp26 * _interp28 * _interp3, fma(-4., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr27(fenv) * _interp21 * _interp26 * _interp28 * _interp3, fma(complex<double>(0.,8.), _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr18(fenv) * _interp25 * _interp26 * _interp28 * _interp3, fma(-4., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr29(fenv) * _interp1 * _interp2 * _interp23 * _interp4, fma(-2., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr37(fenv) * _interp2 * _interp21 * _interp23 * _interp4, fma(-2., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr71(fenv) * _interp1 * _interp22 * _interp23 * _interp4, fma(-1., _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr75(fenv) * _interp21 * _interp22 * _interp23 * _interp4, fma(complex<double>(0.,4.), _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr30(fenv) * _interp2 * _interp23 * _interp25 * _interp4, fma(complex<double>(0.,2.), _cse8_k5 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr72(fenv) * _interp22 * _interp23 * _interp25 * _interp4, 0.)))))))))))));
    // clang-format on
  }
  { // subkernel 6
    const auto _cse1_k6 = powr<-1>(l1);
    const auto _cse2_k6 = -_cse1_k6 * _interp5 * _interp6;
    const auto _cse3_k6 = -_interp6;
    const auto _cse4_k6 = _cse3_k6 + _interp9;
    const auto _cse5_k6 = 50. * _cse4_k6;
    const auto _cse6_k6 = _cse5_k6 + _interp8;
    const auto _cse7_k6 = -_cse1_k6 * _cse6_k6 * _interp7;
    const auto _cse8_k6 = _cse2_k6 + _cse7_k6;
    // clang-format off
    _acc += complex<double>(0.,-0.00510204081632653) * fma(complex<double>(0.,2.), _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr33(fenv) * _interp1 * _interp23 * _interp26 * _interp4, fma(complex<double>(0.,1.), _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr39(fenv) * _interp21 * _interp23 * _interp26 * _interp4, fma(2., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr34(fenv) * _interp23 * _interp25 * _interp26 * _interp4, fma(complex<double>(0.,-4.), _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr7(fenv) * _interp1 * _interp2 * _interp27 * _interp4, fma(complex<double>(0.,-2.), _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr22(fenv) * _interp2 * _interp21 * _interp27 * _interp4, fma(complex<double>(0.,-2.), _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr61(fenv) * _interp1 * _interp22 * _interp27 * _interp4, fma(complex<double>(0.,-1.), _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr68(fenv) * _interp21 * _interp22 * _interp27 * _interp4, fma(-4., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr8(fenv) * _interp2 * _interp25 * _interp27 * _interp4, fma(-2., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr62(fenv) * _interp22 * _interp25 * _interp27 * _interp4, fma(4., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr15(fenv) * _interp1 * _interp26 * _interp27 * _interp4, fma(2., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr26(fenv) * _interp21 * _interp26 * _interp27 * _interp4, fma(complex<double>(0.,-4.), _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr16(fenv) * _interp25 * _interp26 * _interp27 * _interp4, fma(4., _cse8_k6 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr5(fenv) * _interp1 * _interp2 * _interp3 * _interp4, 0.)))))))))))));
    // clang-format on
  }
  { // subkernel 7
    const auto _cse1_k7 = powr<-1>(l1);
    const auto _cse2_k7 = -_cse1_k7 * _interp5 * _interp6;
    const auto _cse3_k7 = -_interp6;
    const auto _cse4_k7 = _cse3_k7 + _interp9;
    const auto _cse5_k7 = 50. * _cse4_k7;
    const auto _cse6_k7 = _cse5_k7 + _interp8;
    const auto _cse7_k7 = -_cse1_k7 * _cse6_k7 * _interp7;
    const auto _cse8_k7 = _cse2_k7 + _cse7_k7;
    // clang-format off
    _acc += 0.00510204081632653 * fma(complex<double>(0.,-2.), _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr21(fenv) * _interp2 * _interp21 * _interp3 * _interp4, fma(complex<double>(0.,-2.), _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr59(fenv) * _interp1 * _interp22 * _interp3 * _interp4, fma(complex<double>(0.,-1.), _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr67(fenv) * _interp21 * _interp22 * _interp3 * _interp4, fma(-4., _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr6(fenv) * _interp2 * _interp25 * _interp3 * _interp4, fma(-2., _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr60(fenv) * _interp22 * _interp25 * _interp3 * _interp4, fma(4., _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr13(fenv) * _interp1 * _interp26 * _interp3 * _interp4, fma(2., _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr25(fenv) * _interp21 * _interp26 * _interp3 * _interp4, fma(complex<double>(0.,-4.), _cse8_k7 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr14(fenv) * _interp25 * _interp26 * _interp3 * _interp4, 0.))))))));
    // clang-format on
  }
  return _acc;
}
static inline auto probe_proj(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k)
{
  const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
  const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
  const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
  const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
  double fenv[(DiFfRG::za4_147_num::nenv) > 0 ? (DiFfRG::za4_147_num::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_1 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
  const double dr_2 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_3 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
  const double dr_4 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_5 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const double dr_6 = ntStub(30415., 1. * l1);
  const double dr_7 = -powr<-1>(l1) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1)) * ntStub(85453., 1. * k) - ntStub(85453., 1. * l1);
  DiFfRG::za4_147_num::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7);
  const auto _interp1 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
  const auto _interp2 = ntStub(12206., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
  const auto _interp3 = ntStub(12206., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
  const auto _interp4 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp5 = ntStub(27191., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp6 = ntStub(85453., 1. * k);
  const auto _interp7 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1), 0.)));
  const auto _interp8 = ntStub(67246.00000000001, 1. * k);
  const auto _interp9 = ntStub(85453., 1.02 * k);
  const auto _interp10 = ntStub(30415., 1. * l1);
  const auto _interp11 = ntStub(85453., 1. * l1);
  const auto _interp12 = ntStub(30415., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp13 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
  const auto _interp14 = ntStub(85453., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp15 = ntStub(30415., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
  const auto _interp16 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p), 0.)));
  const auto _interp17 = ntStub(85453., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
  const auto _interp18 = ntStub(30415., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
  const auto _interp19 = ntStub(80376., fma(1., powr<2>(k), fma(1.618033988749895, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
  const auto _interp20 = ntStub(85453., 1. * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
  const auto _interp21 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
  const auto _interp22 = ntStub(14883., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
  const auto _interp23 = ntStub(14883., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
  const auto _interp24 = ntStub(12206., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp25 = ntStub(38926., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
  const auto _interp26 = ntStub(14883., 0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
  const auto _interp27 = ntStub(38926., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
  const auto _interp28 = ntStub(38926., 0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
  const auto _den7 = powr<-2>(fma(powr<-2>(l1) * powr<2>(_interp6 * _interp7 + _interp11 * l1), powr<2>(l1), powr<2>(_interp10)));
  const auto _den13 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), powr<2>(_interp14 + _interp13 * _interp6 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)))), powr<2>(_interp12)));
  const auto _den14 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p), powr<2>(_interp20 + _interp19 * _interp6 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)))), powr<2>(_interp18)));
  const auto _den15 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p), powr<2>(_interp17 + _interp16 * _interp6 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)))), powr<2>(_interp15)));
  const auto _interp29 = ntStub(13076., 0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
  const auto _interp30 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp31 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp32 = ntStub(29441., 1. * l1);
  const auto _interp33 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp34 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp35 = ntStub(13492., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp36 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp37 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp38 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp39 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
  const auto _interp40 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp41 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp42 = ntStub(13076., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp43 = ntStub(62323., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp44 = ntStub(29441., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp45 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p1 - 3.23606797749979 * cosl1p2) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp46 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp47 = ntStub(62323., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp48 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp49 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp50 = ntStub(13076., 0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp51 = ntStub(29441., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
  const auto _interp52 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1) + (-3.23606797749979 * cosl1p2 - 3.23606797749979 * cosl1p3) * l1 * p + 2.157378651666527 * powr<2>(p));
  const auto _interp53 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp54 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp55 = ntStub(27736., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
  const auto _interp56 = ntStub(27736., 0.4714045207910317 * sqrt(3. * powr<2>(l1) - 3. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 5. * powr<2>(p)));
  const auto _interp57 = ntStub(20621., 1. * l1);
  const auto _interp58 = ntStub(20621., 1. * k);
  const auto _interp59 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
  const auto _interp60 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
  const auto _interp61 = ntStub(20621., 0.5773502691896257 * sqrt(3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _interp62 = ntStub(20621., 1.02 * k);
  const auto _interp63 = ntStub(95661., 1. * k);
  const auto _den8 = 3. * powr<-1>(3. * _interp33 * _interp45 + _interp44 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  const auto _den9 = 3. * powr<-1>(3. * _interp33 * _interp52 + _interp51 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
  const auto _den11 = -3. * powr<-1>(3. * _interp45 * _interp58 + _interp61 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
  // clang-format off
  using _T = decltype(_den11 + _den13 + _den14 + _den15 + _den7 + _den8 + _den9 + DiFfRG::za4_147_num::tr0(fenv) + DiFfRG::za4_147_num::tr12(fenv) + DiFfRG::za4_147_num::tr13(fenv) + DiFfRG::za4_147_num::tr16(fenv) + DiFfRG::za4_147_num::tr18(fenv) + DiFfRG::za4_147_num::tr19(fenv) + DiFfRG::za4_147_num::tr1(fenv) + DiFfRG::za4_147_num::tr22(fenv) + DiFfRG::za4_147_num::tr23(fenv) + DiFfRG::za4_147_num::tr25(fenv) + DiFfRG::za4_147_num::tr28(fenv) + DiFfRG::za4_147_num::tr2(fenv) + DiFfRG::za4_147_num::tr30(fenv) + DiFfRG::za4_147_num::tr31(fenv) + DiFfRG::za4_147_num::tr33(fenv) + DiFfRG::za4_147_num::tr36(fenv) + DiFfRG::za4_147_num::tr38(fenv) + DiFfRG::za4_147_num::tr39(fenv) + DiFfRG::za4_147_num::tr3(fenv) + DiFfRG::za4_147_num::tr42(fenv) + DiFfRG::za4_147_num::tr43(fenv) + DiFfRG::za4_147_num::tr45(fenv) + DiFfRG::za4_147_num::tr48(fenv) + DiFfRG::za4_147_num::tr4(fenv) + DiFfRG::za4_147_num::tr50(fenv) + DiFfRG::za4_147_num::tr51(fenv) + DiFfRG::za4_147_num::tr54(fenv) + DiFfRG::za4_147_num::tr55(fenv) + DiFfRG::za4_147_num::tr58(fenv) + DiFfRG::za4_147_num::tr60(fenv) + DiFfRG::za4_147_num::tr61(fenv) + DiFfRG::za4_147_num::tr63(fenv) + DiFfRG::za4_147_num::tr66(fenv) + DiFfRG::za4_147_num::tr68(fenv) + DiFfRG::za4_147_num::tr69(fenv) + DiFfRG::za4_147_num::tr6(fenv) + DiFfRG::za4_147_num::tr72(fenv) + DiFfRG::za4_147_num::tr73(fenv) + DiFfRG::za4_147_num::tr76(fenv) + DiFfRG::za4_147_num::tr78(fenv) + DiFfRG::za4_147_num::tr79(fenv) + DiFfRG::za4_147_num::tr7(fenv) + DiFfRG::za4_147_num::tr82(fenv) + DiFfRG::za4_147_num::tr84(fenv) + DiFfRG::za4_147_num::tr9(fenv) + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp7 + _interp8 + _interp9 + cosl1p1 + cosl1p2 + cosl1p3 + k + l1 + p);
  // clang-format on
  _T _acc{};
  { // subkernel 1
    const auto _den1 = powr<-1>(1. + powr<6>(k));
    const auto _den2 = powr<-2>(_interp33 * _interp34 + _interp32 * powr<2>(l1));
    const auto _den3 = powr<-2>(_interp34 * _interp58 + _interp57 * powr<2>(l1));
    const auto _den4 = powr<-1>(_interp33 * _interp49 + _interp48 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _den5 = powr<-1>(_interp33 * _interp39 + _interp38 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
    const auto _den6 = -powr<-1>(_interp49 * _interp58 + _interp59 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
    const auto _den10 = powr<-1>(_interp33 * _interp41 + _interp40 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
    const auto _den12 = -powr<-1>(_interp41 * _interp58 + _interp60 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
    // clang-format off
    _acc += fma(0.01020408163265306, _den11 * _den12 * _den3 * _den6 * DiFfRG::za4_147_num::tr4(fenv) * _interp53 * _interp54 * _interp55 * _interp56 * (_interp35 * _interp58 + _interp34 * (50. * (-_interp58 + _interp62) + _interp63)), fma(0.002551020408163265, _den10 * _den2 * _den5 * DiFfRG::za4_147_num::tr2(fenv) * _interp29 * _interp30 * _interp31 * (_interp33 * _interp35 + _interp34 * (_interp36 + 50. * _den1 * (-_interp33 + _interp37) * powr<6>(k))), fma(0.00510204081632653, _den10 * _den2 * _den8 * DiFfRG::za4_147_num::tr1(fenv) * _interp31 * _interp42 * _interp43 * (_interp33 * _interp35 + _interp34 * (_interp36 + 50. * _den1 * (-_interp33 + _interp37) * powr<6>(k))), fma(0.00510204081632653, _den10 * _den2 * _den4 * _den8 * DiFfRG::za4_147_num::tr0(fenv) * _interp31 * _interp43 * _interp46 * _interp47 * (_interp33 * _interp35 + _interp34 * (_interp36 + 50. * _den1 * (-_interp33 + _interp37) * powr<6>(k))), fma(0.001275510204081632, _den2 * _den9 * DiFfRG::za4_147_num::tr3(fenv) * powr<2>(_interp50) * (_interp33 * _interp35 + _interp34 * (_interp36 + 50. * _den1 * (-_interp33 + _interp37) * powr<6>(k))), fma(0.01020408163265306, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr45(fenv) * _interp2 * _interp23 * _interp24 * _interp25 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), fma(0.00510204081632653, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr55(fenv) * _interp22 * _interp23 * _interp24 * _interp25 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), fma(0.00510204081632653, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr51(fenv) * _interp2 * _interp23 * _interp25 * _interp26 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), fma(0.002551020408163265, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr58(fenv) * _interp22 * _interp23 * _interp25 * _interp26 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), fma(-0.01020408163265306, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr43(fenv) * _interp1 * _interp23 * _interp24 * _interp27 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), fma(-0.00510204081632653, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr79(fenv) * _interp21 * _interp23 * _interp24 * _interp27 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), fma(-0.00510204081632653, _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr50(fenv) * _interp1 * _interp23 * _interp26 * _interp27 * (-_interp5 * _interp6 * powr<-1>(l1) - _interp7 * (_interp8 + 50. * (-_interp6 + _interp9)) * powr<-1>(l1)), 0.))))))))))));
    // clang-format on
  }
  { // subkernel 2
    const auto _cse1_k2 = powr<-1>(l1);
    const auto _cse2_k2 = -_cse1_k2 * _interp5 * _interp6;
    const auto _cse3_k2 = -_interp6;
    const auto _cse4_k2 = _cse3_k2 + _interp9;
    const auto _cse5_k2 = 50. * _cse4_k2;
    const auto _cse6_k2 = _cse5_k2 + _interp8;
    const auto _cse7_k2 = -_cse1_k2 * _cse6_k2 * _interp7;
    const auto _cse8_k2 = _cse2_k2 + _cse7_k2;
    // clang-format off
    _acc += 0.002551020408163265 * fma(-1., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr82(fenv) * _interp21 * _interp23 * _interp26 * _interp27, fma(-8., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr9(fenv) * _interp1 * _interp2 * _interp24 * _interp28, fma(-4., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr63(fenv) * _interp2 * _interp21 * _interp24 * _interp28, fma(-4., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr31(fenv) * _interp1 * _interp22 * _interp24 * _interp28, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr73(fenv) * _interp21 * _interp22 * _interp24 * _interp28, fma(-4., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr23(fenv) * _interp1 * _interp2 * _interp26 * _interp28, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr69(fenv) * _interp2 * _interp21 * _interp26 * _interp28, fma(-2., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr38(fenv) * _interp1 * _interp22 * _interp26 * _interp28, fma(-1., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr76(fenv) * _interp21 * _interp22 * _interp26 * _interp28, fma(-8., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr19(fenv) * _interp24 * _interp25 * _interp27 * _interp28, fma(-4., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr28(fenv) * _interp25 * _interp26 * _interp27 * _interp28, fma(8., _cse8_k2 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr13(fenv) * _interp2 * _interp24 * _interp25 * _interp3, 0.))))))))))));
    // clang-format on
  }
  { // subkernel 3
    const auto _cse1_k3 = powr<-1>(l1);
    const auto _cse2_k3 = -_cse1_k3 * _interp5 * _interp6;
    const auto _cse3_k3 = -_interp6;
    const auto _cse4_k3 = _cse3_k3 + _interp9;
    const auto _cse5_k3 = 50. * _cse4_k3;
    const auto _cse6_k3 = _cse5_k3 + _interp8;
    const auto _cse7_k3 = -_cse1_k3 * _cse6_k3 * _interp7;
    const auto _cse8_k3 = _cse2_k3 + _cse7_k3;
    // clang-format off
    _acc += 0.002551020408163265 * fma(4., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr33(fenv) * _interp22 * _interp24 * _interp25 * _interp3, fma(4., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr25(fenv) * _interp2 * _interp25 * _interp26 * _interp3, fma(2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr39(fenv) * _interp22 * _interp25 * _interp26 * _interp3, fma(-8., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr7(fenv) * _interp1 * _interp24 * _interp27 * _interp3, fma(-4., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr61(fenv) * _interp21 * _interp24 * _interp27 * _interp3, fma(-4., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr22(fenv) * _interp1 * _interp26 * _interp27 * _interp3, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr68(fenv) * _interp21 * _interp26 * _interp27 * _interp3, fma(-4., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr42(fenv) * _interp1 * _interp2 * _interp23 * _interp4, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr78(fenv) * _interp2 * _interp21 * _interp23 * _interp4, fma(-2., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr54(fenv) * _interp1 * _interp22 * _interp23 * _interp4, fma(-1., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr84(fenv) * _interp21 * _interp22 * _interp23 * _interp4, fma(-4., _cse8_k3 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr48(fenv) * _interp23 * _interp25 * _interp27 * _interp4, 0.))))))))))));
    // clang-format on
  }
  { // subkernel 4
    const auto _cse1_k4 = powr<-1>(l1);
    const auto _cse2_k4 = -_cse1_k4 * _interp5 * _interp6;
    const auto _cse3_k4 = -_interp6;
    const auto _cse4_k4 = _cse3_k4 + _interp9;
    const auto _cse5_k4 = 50. * _cse4_k4;
    const auto _cse6_k4 = _cse5_k4 + _interp8;
    const auto _cse7_k4 = -_cse1_k4 * _cse6_k4 * _interp7;
    const auto _cse8_k4 = _cse2_k4 + _cse7_k4;
    // clang-format off
    _acc += 0.00510204081632653 * fma(-4., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr18(fenv) * _interp2 * _interp25 * _interp28 * _interp4, fma(-2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr36(fenv) * _interp22 * _interp25 * _interp28 * _interp4, fma(4., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr12(fenv) * _interp1 * _interp27 * _interp28 * _interp4, fma(2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr66(fenv) * _interp21 * _interp27 * _interp28 * _interp4, fma(-4., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr6(fenv) * _interp1 * _interp2 * _interp3 * _interp4, fma(-2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr60(fenv) * _interp2 * _interp21 * _interp3 * _interp4, fma(-2., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr30(fenv) * _interp1 * _interp22 * _interp3 * _interp4, fma(-1., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr72(fenv) * _interp21 * _interp22 * _interp3 * _interp4, fma(-4., _cse8_k4 * _den13 * _den14 * _den15 * _den7 * DiFfRG::za4_147_num::tr16(fenv) * _interp25 * _interp27 * _interp3 * _interp4, 0.)))))))));
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
  for(int n=0;n<4000;++n){ double l1 = U(rng); double cos1 = Uc(rng); double cos2 = Uc(rng); double phi = Uph(rng); double p = U(rng); double k = U(rng);
    std::complex<double> f = probe_full(l1, cos1, cos2, phi, p, k); std::complex<double> pj = probe_proj(l1, cos1, cos2, phi, p, k);
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
