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
#include "ZA_num_kernels.hh"
template<int N> static inline double powr(double x){ double r=1.0; int n=N<0?-N:N; for(int i=0;i<n;++i) r*=x; return N<0?1.0/r:r; }
using std::pow; using std::sqrt; using std::sin; using std::cos; using std::tan; using std::exp; using std::log; using std::fma; using std::fabs;
static inline std::complex<double> fma(const std::complex<double>&a,const std::complex<double>&b,const std::complex<double>&c){return a*b+c;}
template<class T> using complex = std::complex<T>;
static inline double ntStub(double seed, double x){ double h = std::sin(seed*0.1031 + x*0.3127 + 1.7)*43758.5453; return 0.4 + 0.5*(h - std::floor(h)); }
static inline auto probe_full(const double& l1, const double& cos1, const double& p, const double& k)
{
  double fenv[(DiFfRG::za_num::nenv) > 0 ? (DiFfRG::za_num::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * l1);
  const double dr_1 = -sqrt(powr<-1>(powr<2>(l1))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1)) * ntStub(85453., 1. * k) - ntStub(85453., 1. * l1);
  const double dr_2 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const double dr_3 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  DiFfRG::za_num::fill(fenv, l1, cos1, p, dr_0, dr_1, dr_2, dr_3);
  const auto _interp1 = ntStub(12206., 0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
  const auto _interp2 = ntStub(27191., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp3 = ntStub(85453., 1. * k);
  const auto _interp4 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp5 = ntStub(67246.00000000001, 1. * k);
  const auto _interp6 = ntStub(85453., 1.02 * k);
  const auto _interp7 = ntStub(30415., 1. * l1);
  const auto _interp8 = ntStub(85453., 1. * l1);
  const auto _interp9 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp10 = ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp11 = ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp12 = ntStub(13076., 0.7071067811865475 * sqrt(powr<2>(l1) + powr<2>(p)));
  const auto _interp13 = ntStub(29441., 1. * l1);
  const auto _interp14 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp15 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp16 = ntStub(13492., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp17 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp18 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp19 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
  const auto _interp20 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp21 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp22 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
  const auto _interp23 = ntStub(20621., 1. * l1);
  const auto _interp24 = ntStub(20621., 1. * k);
  const auto _interp25 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp26 = ntStub(20621., 1.02 * k);
  const auto _interp27 = ntStub(95661., 1. * k);
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(_interp14 * _interp15 + _interp13 * powr<2>(l1));
  const auto _den3 = powr<-2>(_interp15 * _interp24 + _interp23 * powr<2>(l1));
  const auto _den4 = powr<-1>(_interp14 * _interp21 + _interp20 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _den5 = -powr<-1>(_interp21 * _interp24 + _interp25 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _den6 = powr<-2>(powr<2>(_interp7) + powr<2>(_interp3 * _interp4 + _interp8 * l1));
  const auto _den7 = powr<-1>(powr<2>(_interp9) + powr<2>(_interp10 * _interp3 + _interp11 * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
  const auto _cse1 = powr<-1>(l1);
  const auto _cse2 = powr<-2>(p);
  const auto _cse3 = -_cse1 * _interp2 * _interp3;
  const auto _cse4 = -_interp3;
  const auto _cse5 = _cse4 + _interp6;
  const auto _cse6 = 50. * _cse5;
  const auto _cse7 = _cse6 + _interp5;
  const auto _cse8 = -_cse1 * _cse7 * _interp4;
  const auto _cse9 = _cse3 + _cse8;
  const auto _cse10 = _interp14 * _interp16;
  const auto _cse11 = -_interp14;
  const auto _cse12 = _cse11 + _interp18;
  const auto _cse13 = 50. * _cse12 * _den1 * powr<6>(k);
  const auto _cse14 = _cse13 + _interp17;
  const auto _cse15 = _cse14 * _interp15;
  const auto _cse16 = _cse10 + _cse15;
  const auto _cse17 = _interp16 * _interp24;
  const auto _cse18 = -_interp24;
  const auto _cse19 = _cse18 + _interp26;
  const auto _cse20 = 50. * _cse19;
  const auto _cse21 = _cse20 + _interp27;
  const auto _cse22 = _cse21 * _interp15;
  const auto _cse23 = _cse17 + _cse22;
  return 0.02083333333333333 * fma(complex<double>(0.,4.), _cse2 * _cse9 * _den6 * _den7 * DiFfRG::za_num::tr4(fenv) * powr<2>(_interp1), fma(complex<double>(0.,-4.), _cse2 * _cse9 * _den6 * _den7 * DiFfRG::za_num::tr5(fenv) * powr<2>(_interp1), fma(-1., _cse16 * _cse2 * _den2 * DiFfRG::za_num::tr1(fenv) * _interp12, fma(-2., _cse16 * _cse2 * _den2 * _den4 * DiFfRG::za_num::tr0(fenv) * powr<2>(_interp19), fma(-2., _cse2 * _cse23 * _den3 * _den5 * DiFfRG::za_num::tr2(fenv) * powr<2>(_interp22), fma(2., _cse2 * _cse23 * _den3 * _den5 * DiFfRG::za_num::tr3(fenv) * powr<2>(_interp22), 0.))))));
}
static inline auto probe_proj(const double& l1, const double& cos1, const double& p, const double& k)
{
  double fenv[(DiFfRG::za_num::nenv) > 0 ? (DiFfRG::za_num::nenv) : 1];
  const double dr_0 = ntStub(30415., 1. * l1);
  const double dr_1 = -sqrt(powr<-1>(powr<2>(l1))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1)) * ntStub(85453., 1. * k) - ntStub(85453., 1. * l1);
  const double dr_2 = ntStub(30415., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const double dr_3 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))) * ntStub(80376., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))) * ntStub(85453., 1. * k) - ntStub(85453., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  DiFfRG::za_num::fill(fenv, l1, cos1, p, dr_0, dr_1, dr_2, dr_3);
  const auto _interp1 = ntStub(13076., 0.7071067811865475 * sqrt(powr<2>(l1) + powr<2>(p)));
  const auto _interp2 = ntStub(29441., 1. * l1);
  const auto _interp3 = ntStub(29441., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp4 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp5 = ntStub(13492., 1. * powr<2>(k) + 1.618033988749895 * powr<2>(l1));
  const auto _interp6 = ntStub(18969., 1. * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp7 = ntStub(29441., 1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
  const auto _interp8 = ntStub(62323., 0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
  const auto _interp9 = ntStub(29441., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp10 = ntStub(89189., 1. * powr<2>(k) + 1.618033988749895 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp11 = ntStub(27736., 0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
  const auto _interp12 = ntStub(20621., 1. * l1);
  const auto _interp13 = ntStub(20621., 1. * k);
  const auto _interp14 = ntStub(20621., 1. * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _interp15 = ntStub(20621., 1.02 * k);
  const auto _interp16 = ntStub(95661., 1. * k);
  const auto _den1 = powr<-1>(1. + powr<6>(k));
  const auto _den2 = powr<-2>(_interp3 * _interp4 + _interp2 * powr<2>(l1));
  const auto _den3 = powr<-2>(_interp13 * _interp4 + _interp12 * powr<2>(l1));
  const auto _den4 = powr<-1>(_interp10 * _interp3 + _interp9 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _den5 = -powr<-1>(_interp10 * _interp13 + _interp14 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
  const auto _cse1 = powr<-2>(p);
  const auto _cse2 = -_interp3;
  const auto _cse3 = _cse2 + _interp7;
  const auto _cse4 = 50. * _cse3 * _den1 * powr<6>(k);
  const auto _cse5 = _cse4 + _interp6;
  const auto _cse6 = _cse5 * _interp4;
  const auto _cse7 = _interp3 * _interp5;
  const auto _cse8 = _cse6 + _cse7;
  const auto _cse9 = -_interp13;
  const auto _cse10 = _cse9 + _interp15;
  const auto _cse11 = 50. * _cse10;
  const auto _cse12 = _cse11 + _interp16;
  const auto _cse13 = _cse12 * _interp4;
  const auto _cse14 = _interp13 * _interp5;
  const auto _cse15 = _cse13 + _cse14;
  return 0.02083333333333333 * fma(-1., _cse1 * _cse8 * _den2 * DiFfRG::za_num::tr1(fenv) * _interp1, fma(-2., _cse1 * _cse15 * _den3 * _den5 * DiFfRG::za_num::tr2(fenv) * powr<2>(_interp11), fma(2., _cse1 * _cse15 * _den3 * _den5 * DiFfRG::za_num::tr3(fenv) * powr<2>(_interp11), fma(-2., _cse1 * _cse8 * _den2 * _den4 * DiFfRG::za_num::tr0(fenv) * powr<2>(_interp8), 0.))));
}
int main(int argc, char** argv){
  const char* outf=nullptr; const char* macro=nullptr;
  for(int i=1;i<argc;++i){ if(!std::strcmp(argv[i],"-o") && i+1<argc) outf=argv[++i];
                           else if(!std::strcmp(argv[i],"-m") && i+1<argc) macro=argv[++i]; }
  std::mt19937_64 rng(12345); std::uniform_real_distribution<double> U(0.25,3.0),Uc(-0.9,0.9),Uph(0.1,6.2);
  double mim=0,mdiff=0,mre=0,mrim=0,mrdiff=0; long ok=0;
  for(int n=0;n<4000;++n){ double l1 = U(rng); double cos1 = Uc(rng); double p = U(rng); double k = U(rng);
    std::complex<double> f = probe_full(l1, cos1, p, k); std::complex<double> pj = probe_proj(l1, cos1, p, k);
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
