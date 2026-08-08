#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Ortho12_kernels.hh"

namespace numtracer_kernels
{
  class Ortho12_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p1m, const double& p2m, const double& cosP)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::ortho12::nenv) > 0 ? (numtracer_kernels::ortho12::nenv) : 1];
      numtracer_kernels::ortho12::fill(fenv, l1, cos1, cos2, p1m, p2m, cosP);
      const auto _den1 = powr<-2>(fma(powr<2>(cosP), powr<2>(p2m), fma(-1., powr<2>(p1m + cosP * p2m), 0.)));
      const auto _den2 = powr<-1>(fma(powr<2>(cosP), powr<2>(p2m), fma(2., -(-1. + powr<2>(cosP)) * powr<2>(p2m), fma(2., -p2m * (cosP * p1m + p2m), powr<2>(p1m + cosP * p2m)))));
      const auto _den3 = powr<-2>(fma(-1., powr<2>(p2m) * (powr<2>(p1m) + 2. * cosP * p1m * p2m + powr<2>(p2m)), powr<2>(cosP * p1m * p2m + powr<2>(p2m))));
      const auto _den4 = powr<-1>(fma(-1., powr<2>(p2m) * (powr<2>(p1m) + 2. * cosP * p1m * p2m + powr<2>(p2m)), powr<2>(cosP * p1m * p2m + powr<2>(p2m))));
      const auto _den5 = powr<-1>(fma(-1., powr<2>(cosP * p1m * p2m + powr<2>(p2m)), fma(powr<2>(p2m), powr<2>(p1m) + 2. * cosP * p1m * p2m + powr<2>(p2m), 0.)));
      // clang-format off
      using _T = decltype(_den1 + _den2 + _den3 + _den4 + _den5 + numtracer_kernels::ortho12::tr0(fenv) + numtracer_kernels::ortho12::tr100(fenv) + numtracer_kernels::ortho12::tr101(fenv) + numtracer_kernels::ortho12::tr102(fenv) + numtracer_kernels::ortho12::tr103(fenv) + numtracer_kernels::ortho12::tr104(fenv) + numtracer_kernels::ortho12::tr105(fenv) + numtracer_kernels::ortho12::tr106(fenv) + numtracer_kernels::ortho12::tr107(fenv) + numtracer_kernels::ortho12::tr108(fenv) + numtracer_kernels::ortho12::tr109(fenv) + numtracer_kernels::ortho12::tr10(fenv) + numtracer_kernels::ortho12::tr110(fenv) + numtracer_kernels::ortho12::tr111(fenv) + numtracer_kernels::ortho12::tr112(fenv) + numtracer_kernels::ortho12::tr113(fenv) + numtracer_kernels::ortho12::tr114(fenv) + numtracer_kernels::ortho12::tr115(fenv) + numtracer_kernels::ortho12::tr116(fenv) + numtracer_kernels::ortho12::tr117(fenv) + numtracer_kernels::ortho12::tr118(fenv) + numtracer_kernels::ortho12::tr119(fenv) + numtracer_kernels::ortho12::tr11(fenv) + numtracer_kernels::ortho12::tr120(fenv) + numtracer_kernels::ortho12::tr121(fenv) + numtracer_kernels::ortho12::tr122(fenv) + numtracer_kernels::ortho12::tr123(fenv) + numtracer_kernels::ortho12::tr124(fenv) + numtracer_kernels::ortho12::tr125(fenv) + numtracer_kernels::ortho12::tr126(fenv) + numtracer_kernels::ortho12::tr127(fenv) + numtracer_kernels::ortho12::tr128(fenv) + numtracer_kernels::ortho12::tr129(fenv) + numtracer_kernels::ortho12::tr12(fenv) + numtracer_kernels::ortho12::tr130(fenv) + numtracer_kernels::ortho12::tr131(fenv) + numtracer_kernels::ortho12::tr132(fenv) + numtracer_kernels::ortho12::tr133(fenv) + numtracer_kernels::ortho12::tr134(fenv) + numtracer_kernels::ortho12::tr135(fenv) + numtracer_kernels::ortho12::tr136(fenv) + numtracer_kernels::ortho12::tr137(fenv) + numtracer_kernels::ortho12::tr138(fenv) + numtracer_kernels::ortho12::tr139(fenv) + numtracer_kernels::ortho12::tr13(fenv) + numtracer_kernels::ortho12::tr140(fenv) + numtracer_kernels::ortho12::tr141(fenv) + numtracer_kernels::ortho12::tr142(fenv) + numtracer_kernels::ortho12::tr143(fenv) + numtracer_kernels::ortho12::tr144(fenv) + numtracer_kernels::ortho12::tr145(fenv) + numtracer_kernels::ortho12::tr146(fenv) + numtracer_kernels::ortho12::tr147(fenv) + numtracer_kernels::ortho12::tr148(fenv) + numtracer_kernels::ortho12::tr149(fenv) + numtracer_kernels::ortho12::tr14(fenv) + numtracer_kernels::ortho12::tr150(fenv) + numtracer_kernels::ortho12::tr151(fenv) + numtracer_kernels::ortho12::tr152(fenv) + numtracer_kernels::ortho12::tr153(fenv) + numtracer_kernels::ortho12::tr154(fenv) + numtracer_kernels::ortho12::tr155(fenv) + numtracer_kernels::ortho12::tr156(fenv) + numtracer_kernels::ortho12::tr157(fenv) + numtracer_kernels::ortho12::tr158(fenv) + numtracer_kernels::ortho12::tr159(fenv) + numtracer_kernels::ortho12::tr15(fenv) + numtracer_kernels::ortho12::tr160(fenv) + numtracer_kernels::ortho12::tr161(fenv) + numtracer_kernels::ortho12::tr162(fenv) + numtracer_kernels::ortho12::tr163(fenv) + numtracer_kernels::ortho12::tr164(fenv) + numtracer_kernels::ortho12::tr165(fenv) + numtracer_kernels::ortho12::tr166(fenv) + numtracer_kernels::ortho12::tr167(fenv) + numtracer_kernels::ortho12::tr168(fenv) + numtracer_kernels::ortho12::tr169(fenv) + numtracer_kernels::ortho12::tr16(fenv) + numtracer_kernels::ortho12::tr170(fenv) + numtracer_kernels::ortho12::tr171(fenv) + numtracer_kernels::ortho12::tr172(fenv) + numtracer_kernels::ortho12::tr173(fenv) + numtracer_kernels::ortho12::tr174(fenv) + numtracer_kernels::ortho12::tr175(fenv) + numtracer_kernels::ortho12::tr176(fenv) + numtracer_kernels::ortho12::tr177(fenv) + numtracer_kernels::ortho12::tr178(fenv) + numtracer_kernels::ortho12::tr179(fenv) + numtracer_kernels::ortho12::tr17(fenv) + numtracer_kernels::ortho12::tr180(fenv) + numtracer_kernels::ortho12::tr181(fenv) + numtracer_kernels::ortho12::tr182(fenv) + numtracer_kernels::ortho12::tr183(fenv) + numtracer_kernels::ortho12::tr184(fenv) + numtracer_kernels::ortho12::tr185(fenv) + numtracer_kernels::ortho12::tr186(fenv) + numtracer_kernels::ortho12::tr187(fenv) + numtracer_kernels::ortho12::tr188(fenv) + numtracer_kernels::ortho12::tr189(fenv) + numtracer_kernels::ortho12::tr18(fenv) + numtracer_kernels::ortho12::tr190(fenv) + numtracer_kernels::ortho12::tr191(fenv) + numtracer_kernels::ortho12::tr192(fenv) + numtracer_kernels::ortho12::tr193(fenv) + numtracer_kernels::ortho12::tr194(fenv) + numtracer_kernels::ortho12::tr195(fenv) + numtracer_kernels::ortho12::tr196(fenv) + numtracer_kernels::ortho12::tr197(fenv) + numtracer_kernels::ortho12::tr198(fenv) + numtracer_kernels::ortho12::tr199(fenv) + numtracer_kernels::ortho12::tr19(fenv) + numtracer_kernels::ortho12::tr1(fenv) + numtracer_kernels::ortho12::tr200(fenv) + numtracer_kernels::ortho12::tr201(fenv) + numtracer_kernels::ortho12::tr202(fenv) + numtracer_kernels::ortho12::tr203(fenv) + numtracer_kernels::ortho12::tr204(fenv) + numtracer_kernels::ortho12::tr205(fenv) + numtracer_kernels::ortho12::tr206(fenv) + numtracer_kernels::ortho12::tr207(fenv) + numtracer_kernels::ortho12::tr208(fenv) + numtracer_kernels::ortho12::tr209(fenv) + numtracer_kernels::ortho12::tr20(fenv) + numtracer_kernels::ortho12::tr210(fenv) + numtracer_kernels::ortho12::tr211(fenv) + numtracer_kernels::ortho12::tr212(fenv) + numtracer_kernels::ortho12::tr213(fenv) + numtracer_kernels::ortho12::tr214(fenv) + numtracer_kernels::ortho12::tr215(fenv) + numtracer_kernels::ortho12::tr216(fenv) + numtracer_kernels::ortho12::tr217(fenv) + numtracer_kernels::ortho12::tr218(fenv) + numtracer_kernels::ortho12::tr219(fenv) + numtracer_kernels::ortho12::tr21(fenv) + numtracer_kernels::ortho12::tr220(fenv) + numtracer_kernels::ortho12::tr221(fenv) + numtracer_kernels::ortho12::tr222(fenv) + numtracer_kernels::ortho12::tr223(fenv) + numtracer_kernels::ortho12::tr224(fenv) + numtracer_kernels::ortho12::tr225(fenv) + numtracer_kernels::ortho12::tr226(fenv) + numtracer_kernels::ortho12::tr227(fenv) + numtracer_kernels::ortho12::tr228(fenv) + numtracer_kernels::ortho12::tr229(fenv) + numtracer_kernels::ortho12::tr22(fenv) + numtracer_kernels::ortho12::tr230(fenv) + numtracer_kernels::ortho12::tr231(fenv) + numtracer_kernels::ortho12::tr232(fenv) + numtracer_kernels::ortho12::tr233(fenv) + numtracer_kernels::ortho12::tr234(fenv) + numtracer_kernels::ortho12::tr235(fenv) + numtracer_kernels::ortho12::tr236(fenv) + numtracer_kernels::ortho12::tr237(fenv) + numtracer_kernels::ortho12::tr238(fenv) + numtracer_kernels::ortho12::tr239(fenv) + numtracer_kernels::ortho12::tr23(fenv) + numtracer_kernels::ortho12::tr240(fenv) + numtracer_kernels::ortho12::tr241(fenv) + numtracer_kernels::ortho12::tr242(fenv) + numtracer_kernels::ortho12::tr243(fenv) + numtracer_kernels::ortho12::tr244(fenv) + numtracer_kernels::ortho12::tr245(fenv) + numtracer_kernels::ortho12::tr246(fenv) + numtracer_kernels::ortho12::tr247(fenv) + numtracer_kernels::ortho12::tr248(fenv) + numtracer_kernels::ortho12::tr249(fenv) + numtracer_kernels::ortho12::tr24(fenv) + numtracer_kernels::ortho12::tr250(fenv) + numtracer_kernels::ortho12::tr251(fenv) + numtracer_kernels::ortho12::tr252(fenv) + numtracer_kernels::ortho12::tr253(fenv) + numtracer_kernels::ortho12::tr254(fenv) + numtracer_kernels::ortho12::tr255(fenv) + numtracer_kernels::ortho12::tr256(fenv) + numtracer_kernels::ortho12::tr257(fenv) + numtracer_kernels::ortho12::tr258(fenv) + numtracer_kernels::ortho12::tr259(fenv) + numtracer_kernels::ortho12::tr25(fenv) + numtracer_kernels::ortho12::tr260(fenv) + numtracer_kernels::ortho12::tr261(fenv) + numtracer_kernels::ortho12::tr262(fenv) + numtracer_kernels::ortho12::tr26(fenv) + numtracer_kernels::ortho12::tr27(fenv) + numtracer_kernels::ortho12::tr28(fenv) + numtracer_kernels::ortho12::tr29(fenv) + numtracer_kernels::ortho12::tr2(fenv) + numtracer_kernels::ortho12::tr30(fenv) + numtracer_kernels::ortho12::tr31(fenv) + numtracer_kernels::ortho12::tr32(fenv) + numtracer_kernels::ortho12::tr33(fenv) + numtracer_kernels::ortho12::tr34(fenv) + numtracer_kernels::ortho12::tr35(fenv) + numtracer_kernels::ortho12::tr36(fenv) + numtracer_kernels::ortho12::tr37(fenv) + numtracer_kernels::ortho12::tr38(fenv) + numtracer_kernels::ortho12::tr39(fenv) + numtracer_kernels::ortho12::tr3(fenv) + numtracer_kernels::ortho12::tr40(fenv) + numtracer_kernels::ortho12::tr41(fenv) + numtracer_kernels::ortho12::tr42(fenv) + numtracer_kernels::ortho12::tr43(fenv) + numtracer_kernels::ortho12::tr44(fenv) + numtracer_kernels::ortho12::tr45(fenv) + numtracer_kernels::ortho12::tr46(fenv) + numtracer_kernels::ortho12::tr47(fenv) + numtracer_kernels::ortho12::tr48(fenv) + numtracer_kernels::ortho12::tr49(fenv) + numtracer_kernels::ortho12::tr4(fenv) + numtracer_kernels::ortho12::tr50(fenv) + numtracer_kernels::ortho12::tr51(fenv) + numtracer_kernels::ortho12::tr52(fenv) + numtracer_kernels::ortho12::tr53(fenv) + numtracer_kernels::ortho12::tr54(fenv) + numtracer_kernels::ortho12::tr55(fenv) + numtracer_kernels::ortho12::tr56(fenv) + numtracer_kernels::ortho12::tr57(fenv) + numtracer_kernels::ortho12::tr58(fenv) + numtracer_kernels::ortho12::tr59(fenv) + numtracer_kernels::ortho12::tr5(fenv) + numtracer_kernels::ortho12::tr60(fenv) + numtracer_kernels::ortho12::tr61(fenv) + numtracer_kernels::ortho12::tr62(fenv) + numtracer_kernels::ortho12::tr63(fenv) + numtracer_kernels::ortho12::tr64(fenv) + numtracer_kernels::ortho12::tr65(fenv) + numtracer_kernels::ortho12::tr66(fenv) + numtracer_kernels::ortho12::tr67(fenv) + numtracer_kernels::ortho12::tr68(fenv) + numtracer_kernels::ortho12::tr69(fenv) + numtracer_kernels::ortho12::tr6(fenv) + numtracer_kernels::ortho12::tr70(fenv) + numtracer_kernels::ortho12::tr71(fenv) + numtracer_kernels::ortho12::tr72(fenv) + numtracer_kernels::ortho12::tr73(fenv) + numtracer_kernels::ortho12::tr74(fenv) + numtracer_kernels::ortho12::tr75(fenv) + numtracer_kernels::ortho12::tr76(fenv) + numtracer_kernels::ortho12::tr77(fenv) + numtracer_kernels::ortho12::tr78(fenv) + numtracer_kernels::ortho12::tr79(fenv) + numtracer_kernels::ortho12::tr7(fenv) + numtracer_kernels::ortho12::tr80(fenv) + numtracer_kernels::ortho12::tr81(fenv) + numtracer_kernels::ortho12::tr82(fenv) + numtracer_kernels::ortho12::tr83(fenv) + numtracer_kernels::ortho12::tr84(fenv) + numtracer_kernels::ortho12::tr85(fenv) + numtracer_kernels::ortho12::tr86(fenv) + numtracer_kernels::ortho12::tr87(fenv) + numtracer_kernels::ortho12::tr88(fenv) + numtracer_kernels::ortho12::tr89(fenv) + numtracer_kernels::ortho12::tr8(fenv) + numtracer_kernels::ortho12::tr90(fenv) + numtracer_kernels::ortho12::tr91(fenv) + numtracer_kernels::ortho12::tr92(fenv) + numtracer_kernels::ortho12::tr93(fenv) + numtracer_kernels::ortho12::tr94(fenv) + numtracer_kernels::ortho12::tr95(fenv) + numtracer_kernels::ortho12::tr96(fenv) + numtracer_kernels::ortho12::tr97(fenv) + numtracer_kernels::ortho12::tr98(fenv) + numtracer_kernels::ortho12::tr99(fenv) + numtracer_kernels::ortho12::tr9(fenv) + cosP + p1m + p2m);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        _acc += 0.002604166666666666 * fma(831., _den4 * numtracer_kernels::ortho12::tr105(fenv), fma(-2202., _den4 * numtracer_kernels::ortho12::tr116(fenv), fma(-2334., _den4 * numtracer_kernels::ortho12::tr117(fenv), fma(3176., _den2 * numtracer_kernels::ortho12::tr118(fenv), fma(2382., _den4 * numtracer_kernels::ortho12::tr119(fenv), fma(-1203., _den4 * numtracer_kernels::ortho12::tr120(fenv), fma(-2514., _den4 * numtracer_kernels::ortho12::tr121(fenv), fma(2586., _den4 * numtracer_kernels::ortho12::tr122(fenv), fma(2658., _den4 * numtracer_kernels::ortho12::tr123(fenv), fma(1268., numtracer_kernels::ortho12::tr12(fenv), 0.))))))))));
      }
      { // subkernel 2
        _acc += 0.001302083333333333 * fma(32., numtracer_kernels::ortho12::tr14(fenv), fma(136., numtracer_kernels::ortho12::tr15(fenv), fma(978., _den4 * numtracer_kernels::ortho12::tr151(fenv) * numtracer_kernels::ortho12::tr223(fenv) * numtracer_kernels::ortho12::tr224(fenv), fma(579., _den4 * numtracer_kernels::ortho12::tr153(fenv) * numtracer_kernels::ortho12::tr227(fenv) * numtracer_kernels::ortho12::tr228(fenv), fma(1194., _den4 * numtracer_kernels::ortho12::tr155(fenv) * numtracer_kernels::ortho12::tr231(fenv) * numtracer_kernels::ortho12::tr232(fenv), fma(669., _den4 * numtracer_kernels::ortho12::tr157(fenv) * numtracer_kernels::ortho12::tr235(fenv) * numtracer_kernels::ortho12::tr236(fenv), fma(-1362., _den4 * numtracer_kernels::ortho12::tr30(fenv), fma(2008., _den2 * numtracer_kernels::ortho12::tr31(fenv), fma(-770.9999999999999, _den4 * numtracer_kernels::ortho12::tr33(fenv), fma(-1614., _den4 * numtracer_kernels::ortho12::tr34(fenv), 0.))))))))));
      }
      { // subkernel 3
        _acc += 0.002604166666666666 * fma(2694., _den4 * numtracer_kernels::ortho12::tr36(fenv), fma(-2742., _den4 * numtracer_kernels::ortho12::tr37(fenv), fma(1437., _den4 * numtracer_kernels::ortho12::tr39(fenv), fma(88., numtracer_kernels::ortho12::tr40(fenv), fma(2946., _den4 * numtracer_kernels::ortho12::tr41(fenv), fma(-104., numtracer_kernels::ortho12::tr42(fenv), fma(1509., _den4 * numtracer_kernels::ortho12::tr43(fenv), fma(-4456., numtracer_kernels::ortho12::tr46(fenv), fma(-248., numtracer_kernels::ortho12::tr61(fenv), fma(5608., numtracer_kernels::ortho12::tr6(fenv), 0.))))))))));
      }
      { // subkernel 4
        const auto _cse1 = powr<2>(cosP) * powr<2>(p2m);
        const auto _cse2 = -p1m;
        const auto _cse3 = -cosP * p2m;
        const auto _cse4 = _cse2 + _cse3;
        const auto _cse5 = -powr<2>(_cse4);
        const auto _cse6 = _cse1 + _cse5;
        _acc += 0.002604166666666666 * fma(-831., _cse6 * _den4 * numtracer_kernels::ortho12::tr104(fenv), fma(-951., _cse6 * _den4 * numtracer_kernels::ortho12::tr109(fenv), fma(807., _cse6 * _den4 * numtracer_kernels::ortho12::tr10(fenv), fma(-141., _cse6 * _den4 * numtracer_kernels::ortho12::tr19(fenv), fma(159., _cse6 * _den4 * numtracer_kernels::ortho12::tr21(fenv), fma(184., numtracer_kernels::ortho12::tr7(fenv), fma(501., _den4 * numtracer_kernels::ortho12::tr84(fenv), fma(-519., _den4 * numtracer_kernels::ortho12::tr86(fenv), fma(-723., _den4 * numtracer_kernels::ortho12::tr96(fenv), fma(752.9999999999999, _den4 * numtracer_kernels::ortho12::tr99(fenv), 0.))))))))));
      }
      { // subkernel 5
        const auto _cse1 = powr<2>(cosP) * powr<2>(p2m);
        const auto _cse2 = -p1m;
        const auto _cse3 = -cosP * p2m;
        const auto _cse4 = _cse2 + _cse3;
        const auto _cse5 = -powr<2>(_cse4);
        const auto _cse6 = _cse1 + _cse5;
        // clang-format off
        _acc += 0.00390625 * fma(-326., _cse6 * _den4 * numtracer_kernels::ortho12::tr150(fenv) * numtracer_kernels::ortho12::tr221(fenv) * numtracer_kernels::ortho12::tr222(fenv), fma(-193., _cse6 * _den4 * numtracer_kernels::ortho12::tr152(fenv) * numtracer_kernels::ortho12::tr225(fenv) * numtracer_kernels::ortho12::tr226(fenv), fma(-398., _cse6 * _den4 * numtracer_kernels::ortho12::tr154(fenv) * numtracer_kernels::ortho12::tr229(fenv) * numtracer_kernels::ortho12::tr230(fenv), fma(-223., _cse6 * _den4 * numtracer_kernels::ortho12::tr156(fenv) * numtracer_kernels::ortho12::tr233(fenv) * numtracer_kernels::ortho12::tr234(fenv), fma(206., _cse6 * _den4 * numtracer_kernels::ortho12::tr23(fenv), fma(-214., _cse6 * _den4 * numtracer_kernels::ortho12::tr25(fenv), fma(334., _cse6 * _den4 * numtracer_kernels::ortho12::tr26(fenv), fma(-346., _cse6 * _den4 * numtracer_kernels::ortho12::tr28(fenv), fma(257., _cse6 * _den4 * numtracer_kernels::ortho12::tr32(fenv), fma(-86., _cse6 * _den4 * numtracer_kernels::ortho12::tr64(fenv), 0.))))))))));
        // clang-format on
      }
      { // subkernel 6
        const auto _cse1 = powr<2>(cosP) * powr<2>(p2m);
        const auto _cse2 = -p1m;
        const auto _cse3 = -cosP * p2m;
        const auto _cse4 = _cse2 + _cse3;
        const auto _cse5 = -powr<2>(_cse4);
        const auto _cse6 = _cse1 + _cse5;
        _acc += 0.00390625 * fma(-71., _cse6 * _den4 * numtracer_kernels::ortho12::tr66(fenv), fma(-158., _cse6 * _den4 * numtracer_kernels::ortho12::tr68(fenv), fma(-89., _cse6 * _den4 * numtracer_kernels::ortho12::tr70(fenv), fma(202., _cse6 * _den4 * numtracer_kernels::ortho12::tr72(fenv), fma(131., _cse6 * _den4 * numtracer_kernels::ortho12::tr76(fenv), fma(278., _cse6 * _den4 * numtracer_kernels::ortho12::tr78(fenv), fma(151., _cse6 * _den4 * numtracer_kernels::ortho12::tr80(fenv), fma(326., _cse6 * _den4 * numtracer_kernels::ortho12::tr81(fenv), fma(-334., _cse6 * _den4 * numtracer_kernels::ortho12::tr83(fenv), fma(346., _cse6 * _den4 * numtracer_kernels::ortho12::tr85(fenv), 0.))))))))));
      }
      { // subkernel 7
        const auto _cse1 = -powr<2>(cosP);
        const auto _cse2 = 1. + _cse1;
        const auto _cse3 = -p1m;
        const auto _cse4 = -cosP * p2m;
        const auto _cse5 = _cse3 + _cse4;
        const auto _cse6 = powr<2>(cosP) * powr<2>(p2m);
        const auto _cse7 = 2. * _cse2 * powr<2>(p2m);
        const auto _cse8 = -_cse2 * powr<2>(p2m);
        const auto _cse9 = _cse5 * cosP * p2m;
        const auto _cse10 = _cse8 + _cse9;
        const auto _cse11 = -2. * _cse10;
        const auto _cse12 = _cse11 + powr<2>(_cse5) + _cse6 + _cse7;
        const auto _cse13 = -powr<2>(_cse5);
        const auto _cse14 = _cse13 + _cse6;
        _acc += 0.001302083333333333 * fma(-1002., _cse12 * _den4 * numtracer_kernels::ortho12::tr27(fenv), fma(1038., _cse12 * _den4 * numtracer_kernels::ortho12::tr29(fenv), fma(8912., _cse12 * numtracer_kernels::ortho12::tr4(fenv), fma(579., _cse14 * _den4 * numtracer_kernels::ortho12::tr87(fenv), fma(1194., _cse14 * _den4 * numtracer_kernels::ortho12::tr89(fenv), fma(1362., _cse14 * _den4 * numtracer_kernels::ortho12::tr8(fenv), fma(669., _cse14 * _den4 * numtracer_kernels::ortho12::tr91(fenv), fma(1446., _cse14 * _den4 * numtracer_kernels::ortho12::tr95(fenv), fma(-1506., _cse14 * _den4 * numtracer_kernels::ortho12::tr98(fenv), fma(-2008., _cse14 * _den2 * numtracer_kernels::ortho12::tr9(fenv), 0.))))))))));
      }
      { // subkernel 8
        const auto _cse1 = -powr<2>(cosP);
        const auto _cse2 = 1. + _cse1;
        const auto _cse3 = -p1m;
        const auto _cse4 = -cosP * p2m;
        const auto _cse5 = _cse3 + _cse4;
        const auto _cse6 = powr<2>(cosP) * powr<2>(p2m);
        const auto _cse7 = 2. * _cse2 * powr<2>(p2m);
        const auto _cse8 = -_cse2 * powr<2>(p2m);
        const auto _cse9 = _cse5 * cosP * p2m;
        const auto _cse10 = _cse8 + _cse9;
        const auto _cse11 = 2. * _cse10;
        const auto _cse12 = _cse11 + powr<2>(_cse5) + _cse6 + _cse7;
        const auto _cse13 = -2. * _cse10;
        const auto _cse14 = _cse13 + powr<2>(_cse5) + _cse6 + _cse7;
        _acc += 0.00390625 * fma(-2908., _cse12 * _den4 * numtracer_kernels::ortho12::tr0(fenv), fma(-626., _cse12 * _den4 * numtracer_kernels::ortho12::tr108(fenv), fma(634., _cse12 * _den4 * numtracer_kernels::ortho12::tr110(fenv), fma(706., _cse12 * _den4 * numtracer_kernels::ortho12::tr115(fenv), fma(-566., _cse12 * _den4 * numtracer_kernels::ortho12::tr11(fenv), fma(1347., _cse12 * _den3 * numtracer_kernels::ortho12::tr124(fenv), fma(-326., _cse14 * _den4 * numtracer_kernels::ortho12::tr82(fenv), fma(-193., _cse14 * _den4 * numtracer_kernels::ortho12::tr88(fenv), fma(-398., _cse14 * _den4 * numtracer_kernels::ortho12::tr90(fenv), fma(-223., _cse14 * _den4 * numtracer_kernels::ortho12::tr92(fenv), 0.))))))))));
      }
      { // subkernel 9
        const auto _cse1 = -powr<2>(cosP);
        const auto _cse2 = 1. + _cse1;
        const auto _cse3 = -p1m;
        const auto _cse4 = -cosP * p2m;
        const auto _cse5 = _cse3 + _cse4;
        const auto _cse6 = powr<2>(cosP) * powr<2>(p2m);
        const auto _cse7 = 2. * _cse2 * powr<2>(p2m);
        const auto _cse8 = -_cse2 * powr<2>(p2m);
        const auto _cse9 = _cse5 * cosP * p2m;
        const auto _cse10 = _cse8 + _cse9;
        const auto _cse11 = 2. * _cse10;
        const auto _cse12 = _cse11 + powr<2>(_cse5) + _cse6 + _cse7;
        // clang-format off
        _acc += 0.001953125 * fma(-1388., _cse12 * _den4 * numtracer_kernels::ortho12::tr13(fenv), fma(1437., _cse12 * _den3 * numtracer_kernels::ortho12::tr125(fenv) * numtracer_kernels::ortho12::tr171(fenv) * numtracer_kernels::ortho12::tr172(fenv), fma(2946., _cse12 * _den3 * numtracer_kernels::ortho12::tr126(fenv) * numtracer_kernels::ortho12::tr173(fenv) * numtracer_kernels::ortho12::tr174(fenv), fma(1509., _cse12 * _den3 * numtracer_kernels::ortho12::tr127(fenv) * numtracer_kernels::ortho12::tr175(fenv) * numtracer_kernels::ortho12::tr176(fenv), fma(4376., _cse12 * _den4 * numtracer_kernels::ortho12::tr128(fenv) * numtracer_kernels::ortho12::tr177(fenv) * numtracer_kernels::ortho12::tr178(fenv), fma(-2836., _cse12 * _den4 * numtracer_kernels::ortho12::tr17(fenv), fma(-4456., _cse12 * _den4 * numtracer_kernels::ortho12::tr129(fenv) * numtracer_kernels::ortho12::tr179(fenv) * numtracer_kernels::ortho12::tr180(fenv), fma(-4696., _cse12 * _den4 * numtracer_kernels::ortho12::tr130(fenv) * numtracer_kernels::ortho12::tr181(fenv) * numtracer_kernels::ortho12::tr182(fenv), fma(-19232., _cse12 * _den1 * numtracer_kernels::ortho12::tr131(fenv) * numtracer_kernels::ortho12::tr183(fenv) * numtracer_kernels::ortho12::tr184(fenv), fma(-10256., _cse12 * _den1 * numtracer_kernels::ortho12::tr132(fenv) * numtracer_kernels::ortho12::tr185(fenv) * numtracer_kernels::ortho12::tr186(fenv), 0.))))))))));
        // clang-format on
      }
      { // subkernel 10
        const auto _cse1 = -powr<2>(cosP);
        const auto _cse2 = 1. + _cse1;
        const auto _cse3 = -p1m;
        const auto _cse4 = -cosP * p2m;
        const auto _cse5 = _cse3 + _cse4;
        const auto _cse6 = powr<2>(cosP) * powr<2>(p2m);
        const auto _cse7 = 2. * _cse2 * powr<2>(p2m);
        const auto _cse8 = -_cse2 * powr<2>(p2m);
        const auto _cse9 = _cse5 * cosP * p2m;
        const auto _cse10 = _cse8 + _cse9;
        const auto _cse11 = 2. * _cse10;
        const auto _cse12 = _cse11 + powr<2>(_cse5) + _cse6 + _cse7;
        // clang-format off
        _acc += 0.00390625 * fma(-10352., _cse12 * _den1 * numtracer_kernels::ortho12::tr133(fenv) * numtracer_kernels::ortho12::tr187(fenv) * numtracer_kernels::ortho12::tr188(fenv), fma(94., _cse12 * _den4 * numtracer_kernels::ortho12::tr18(fenv), fma(-5272., _cse12 * _den1 * numtracer_kernels::ortho12::tr134(fenv) * numtracer_kernels::ortho12::tr189(fenv) * numtracer_kernels::ortho12::tr190(fenv), fma(6056., _cse12 * _den1 * _den5 * numtracer_kernels::ortho12::tr138(fenv) * numtracer_kernels::ortho12::tr197(fenv) * numtracer_kernels::ortho12::tr198(fenv), fma(3236., _cse12 * _den1 * _den5 * numtracer_kernels::ortho12::tr139(fenv) * numtracer_kernels::ortho12::tr199(fenv) * numtracer_kernels::ortho12::tr200(fenv), fma(6568., _cse12 * _den1 * _den5 * numtracer_kernels::ortho12::tr140(fenv) * numtracer_kernels::ortho12::tr201(fenv) * numtracer_kernels::ortho12::tr202(fenv), fma(3308., _cse12 * _den1 * _den5 * numtracer_kernels::ortho12::tr141(fenv) * numtracer_kernels::ortho12::tr203(fenv) * numtracer_kernels::ortho12::tr204(fenv), fma(86., _cse12 * _den4 * numtracer_kernels::ortho12::tr142(fenv) * numtracer_kernels::ortho12::tr205(fenv) * numtracer_kernels::ortho12::tr206(fenv), fma(71., _cse12 * _den4 * numtracer_kernels::ortho12::tr143(fenv) * numtracer_kernels::ortho12::tr207(fenv) * numtracer_kernels::ortho12::tr208(fenv), fma(-106., _cse12 * _den4 * numtracer_kernels::ortho12::tr20(fenv), 0.))))))))));
        // clang-format on
      }
      { // subkernel 11
        const auto _cse1 = -powr<2>(cosP);
        const auto _cse2 = 1. + _cse1;
        const auto _cse3 = -p1m;
        const auto _cse4 = -cosP * p2m;
        const auto _cse5 = _cse3 + _cse4;
        const auto _cse6 = powr<2>(cosP) * powr<2>(p2m);
        const auto _cse7 = 2. * _cse2 * powr<2>(p2m);
        const auto _cse8 = -_cse2 * powr<2>(p2m);
        const auto _cse9 = _cse5 * cosP * p2m;
        const auto _cse10 = _cse8 + _cse9;
        const auto _cse11 = 2. * _cse10;
        const auto _cse12 = _cse11 + powr<2>(_cse5) + _cse6 + _cse7;
        // clang-format off
        _acc += 0.00390625 * fma(158., _cse12 * _den4 * numtracer_kernels::ortho12::tr144(fenv) * numtracer_kernels::ortho12::tr209(fenv) * numtracer_kernels::ortho12::tr210(fenv), fma(89., _cse12 * _den4 * numtracer_kernels::ortho12::tr145(fenv) * numtracer_kernels::ortho12::tr211(fenv) * numtracer_kernels::ortho12::tr212(fenv), fma(202., _cse12 * _den4 * numtracer_kernels::ortho12::tr146(fenv) * numtracer_kernels::ortho12::tr213(fenv) * numtracer_kernels::ortho12::tr214(fenv), fma(131., _cse12 * _den4 * numtracer_kernels::ortho12::tr147(fenv) * numtracer_kernels::ortho12::tr215(fenv) * numtracer_kernels::ortho12::tr216(fenv), fma(278., _cse12 * _den4 * numtracer_kernels::ortho12::tr148(fenv) * numtracer_kernels::ortho12::tr217(fenv) * numtracer_kernels::ortho12::tr218(fenv), fma(151., _cse12 * _den4 * numtracer_kernels::ortho12::tr149(fenv) * numtracer_kernels::ortho12::tr219(fenv) * numtracer_kernels::ortho12::tr220(fenv), fma(-206., _cse12 * _den4 * numtracer_kernels::ortho12::tr22(fenv), fma(214., _cse12 * _den4 * numtracer_kernels::ortho12::tr24(fenv), fma(1329., _cse12 * _den3 * numtracer_kernels::ortho12::tr170(fenv) * numtracer_kernels::ortho12::tr261(fenv) * numtracer_kernels::ortho12::tr262(fenv), fma(-331., _cse12 * _den4 * numtracer_kernels::ortho12::tr35(fenv), 0.))))))))));
        // clang-format on
      }
      { // subkernel 12
        const auto _cse1 = -powr<2>(cosP);
        const auto _cse2 = 1. + _cse1;
        const auto _cse3 = -p1m;
        const auto _cse4 = -cosP * p2m;
        const auto _cse5 = _cse3 + _cse4;
        const auto _cse6 = powr<2>(cosP) * powr<2>(p2m);
        const auto _cse7 = 2. * _cse2 * powr<2>(p2m);
        const auto _cse8 = -_cse2 * powr<2>(p2m);
        const auto _cse9 = _cse5 * cosP * p2m;
        const auto _cse10 = _cse8 + _cse9;
        const auto _cse11 = 2. * _cse10;
        const auto _cse12 = _cse11 + powr<2>(_cse5) + _cse6 + _cse7;
        _acc += 0.001302083333333333 * fma(-4113., _cse12 * _den3 * numtracer_kernels::ortho12::tr38(fenv), fma(6108., _cse12 * _den4 * numtracer_kernels::ortho12::tr44(fenv), fma(3378., _cse12 * _den4 * numtracer_kernels::ortho12::tr48(fenv), fma(6852.000000000001, _cse12 * _den4 * numtracer_kernels::ortho12::tr49(fenv), fma(-29136., _cse12 * _den1 * numtracer_kernels::ortho12::tr51(fenv), fma(29424., _cse12 * _den1 * numtracer_kernels::ortho12::tr52(fenv), fma(-8292., _cse12 * _den4 * numtracer_kernels::ortho12::tr54(fenv), fma(-11216., _cse12 * _den1 * numtracer_kernels::ortho12::tr55(fenv), fma(8412., _cse12 * _den4 * numtracer_kernels::ortho12::tr56(fenv), fma(8868., _cse12 * _den4 * numtracer_kernels::ortho12::tr59(fenv), 0.))))))))));
      }
      { // subkernel 13
        const auto _cse1 = -powr<2>(cosP);
        const auto _cse2 = 1. + _cse1;
        const auto _cse3 = -p1m;
        const auto _cse4 = -cosP * p2m;
        const auto _cse5 = _cse3 + _cse4;
        const auto _cse6 = powr<2>(cosP) * powr<2>(p2m);
        const auto _cse7 = 2. * _cse2 * powr<2>(p2m);
        const auto _cse8 = -_cse2 * powr<2>(p2m);
        const auto _cse9 = _cse5 * cosP * p2m;
        const auto _cse10 = _cse8 + _cse9;
        const auto _cse11 = 2. * _cse10;
        const auto _cse12 = _cse11 + powr<2>(_cse5) + _cse6 + _cse7;
        _acc += 0.00390625 * fma(-2644., _cse12 * _den4 * numtracer_kernels::ortho12::tr5(fenv), fma(6087.999999999999, _cse12 * _den1 * _den5 * numtracer_kernels::ortho12::tr60(fenv), fma(-6152., _cse12 * _den1 * _den5 * numtracer_kernels::ortho12::tr62(fenv), fma(86., _cse12 * _den4 * numtracer_kernels::ortho12::tr63(fenv), fma(71., _cse12 * _den4 * numtracer_kernels::ortho12::tr65(fenv), fma(158., _cse12 * _den4 * numtracer_kernels::ortho12::tr67(fenv), fma(89., _cse12 * _den4 * numtracer_kernels::ortho12::tr69(fenv), fma(-202., _cse12 * _den4 * numtracer_kernels::ortho12::tr71(fenv), fma(206., _cse12 * _den4 * numtracer_kernels::ortho12::tr73(fenv), fma(-214., _cse12 * _den4 * numtracer_kernels::ortho12::tr74(fenv), 0.))))))))));
      }
      { // subkernel 14
        const auto _cse1 = -powr<2>(cosP);
        const auto _cse2 = 1. + _cse1;
        const auto _cse3 = -p1m;
        const auto _cse4 = -cosP * p2m;
        const auto _cse5 = _cse3 + _cse4;
        const auto _cse6 = powr<2>(cosP) * powr<2>(p2m);
        const auto _cse7 = 2. * _cse2 * powr<2>(p2m);
        const auto _cse8 = -_cse2 * powr<2>(p2m);
        const auto _cse9 = _cse5 * cosP * p2m;
        const auto _cse10 = _cse8 + _cse9;
        const auto _cse11 = 2. * _cse10;
        const auto _cse12 = _cse11 + powr<2>(_cse5) + _cse6 + _cse7;
        const auto _cse13 = -powr<2>(_cse5);
        const auto _cse14 = _cse13 + _cse6;
        _acc += 0.0009765625 * fma(770.9999999999999, _cse12 * _cse14 * _den3 * numtracer_kernels::ortho12::tr101(fenv), fma(1614., _cse12 * _cse14 * _den3 * numtracer_kernels::ortho12::tr103(fenv), fma(1698., _cse12 * _cse14 * _den3 * numtracer_kernels::ortho12::tr106(fenv), fma(993., _cse12 * _cse14 * _den3 * numtracer_kernels::ortho12::tr111(fenv), fma(2082., _cse12 * _cse14 * _den3 * numtracer_kernels::ortho12::tr113(fenv), fma(1446., _cse12 * _cse14 * _den3 * numtracer_kernels::ortho12::tr159(fenv) * numtracer_kernels::ortho12::tr239(fenv) * numtracer_kernels::ortho12::tr240(fenv), fma(-1506., _cse12 * _cse14 * _den3 * numtracer_kernels::ortho12::tr161(fenv) * numtracer_kernels::ortho12::tr243(fenv) * numtracer_kernels::ortho12::tr244(fenv), fma(-524., _cse12 * _den4 * numtracer_kernels::ortho12::tr75(fenv), fma(-1112., _cse12 * _den4 * numtracer_kernels::ortho12::tr77(fenv), fma(-604., _cse12 * _den4 * numtracer_kernels::ortho12::tr79(fenv), 0.))))))))));
      }
      { // subkernel 15
        const auto _cse1 = powr<2>(cosP) * powr<2>(p2m);
        const auto _cse2 = -p1m;
        const auto _cse3 = -cosP * p2m;
        const auto _cse4 = _cse2 + _cse3;
        const auto _cse5 = -powr<2>(cosP);
        const auto _cse6 = 1. + _cse5;
        const auto _cse7 = -powr<2>(_cse4);
        const auto _cse8 = _cse1 + _cse7;
        const auto _cse9 = 2. * _cse6 * powr<2>(p2m);
        const auto _cse10 = -_cse6 * powr<2>(p2m);
        const auto _cse11 = _cse4 * cosP * p2m;
        const auto _cse12 = _cse10 + _cse11;
        const auto _cse13 = 2. * _cse12;
        const auto _cse14 = _cse1 + _cse13 + powr<2>(_cse4) + _cse9;
        const auto _cse15 = -2. * _cse12;
        const auto _cse16 = _cse1 + _cse15 + powr<2>(_cse4) + _cse9;
        // clang-format off
        _acc += 0.001953125 * fma(-2252., _cse14 * _cse16 * _den4 * numtracer_kernels::ortho12::tr16(fenv), fma(-4568., _cse14 * _cse16 * _den4 * numtracer_kernels::ortho12::tr1(fenv), fma(-831., _cse14 * _cse8 * _den3 * numtracer_kernels::ortho12::tr163(fenv) * numtracer_kernels::ortho12::tr247(fenv) * numtracer_kernels::ortho12::tr248(fenv), fma(939., _cse14 * _cse8 * _den3 * numtracer_kernels::ortho12::tr164(fenv) * numtracer_kernels::ortho12::tr249(fenv) * numtracer_kernels::ortho12::tr250(fenv), fma(-951., _cse14 * _cse8 * _den3 * numtracer_kernels::ortho12::tr166(fenv) * numtracer_kernels::ortho12::tr253(fenv) * numtracer_kernels::ortho12::tr254(fenv), fma(-1059., _cse14 * _cse8 * _den3 * numtracer_kernels::ortho12::tr168(fenv) * numtracer_kernels::ortho12::tr257(fenv) * numtracer_kernels::ortho12::tr258(fenv), fma(-4072., _cse14 * _cse16 * _den4 * numtracer_kernels::ortho12::tr3(fenv), fma(-4376., _cse14 * _cse16 * _den4 * numtracer_kernels::ortho12::tr45(fenv), fma(4456., _cse14 * _cse16 * _den4 * numtracer_kernels::ortho12::tr47(fenv), fma(681., _cse14 * _cse8 * _den3 * numtracer_kernels::ortho12::tr94(fenv), 0.))))))))));
        // clang-format on
      }
      { // subkernel 16
        const auto _cse1 = -powr<2>(cosP);
        const auto _cse2 = 1. + _cse1;
        const auto _cse3 = -p1m;
        const auto _cse4 = -cosP * p2m;
        const auto _cse5 = _cse3 + _cse4;
        const auto _cse6 = powr<2>(cosP) * powr<2>(p2m);
        const auto _cse7 = 2. * _cse2 * powr<2>(p2m);
        const auto _cse8 = -_cse2 * powr<2>(p2m);
        const auto _cse9 = _cse5 * cosP * p2m;
        const auto _cse10 = _cse8 + _cse9;
        const auto _cse11 = 2. * _cse10;
        const auto _cse12 = _cse11 + powr<2>(_cse5) + _cse6 + _cse7;
        // clang-format off
        _acc += 0.0009765625 * fma(-1698., powr<2>(_cse12) * _den3 * numtracer_kernels::ortho12::tr107(fenv), fma(-993., powr<2>(_cse12) * _den3 * numtracer_kernels::ortho12::tr112(fenv), fma(-2082., powr<2>(_cse12) * _den3 * numtracer_kernels::ortho12::tr114(fenv), fma(11056., powr<2>(_cse12) * _den1 * _den4 * numtracer_kernels::ortho12::tr135(fenv) * numtracer_kernels::ortho12::tr191(fenv) * numtracer_kernels::ortho12::tr192(fenv), fma(-11216., powr<2>(_cse12) * _den1 * _den4 * numtracer_kernels::ortho12::tr136(fenv) * numtracer_kernels::ortho12::tr193(fenv) * numtracer_kernels::ortho12::tr194(fenv), fma(-11824., powr<2>(_cse12) * _den1 * _den4 * numtracer_kernels::ortho12::tr137(fenv) * numtracer_kernels::ortho12::tr195(fenv) * numtracer_kernels::ortho12::tr196(fenv), fma(-1878., powr<2>(_cse12) * _den3 * numtracer_kernels::ortho12::tr165(fenv) * numtracer_kernels::ortho12::tr251(fenv) * numtracer_kernels::ortho12::tr252(fenv), fma(1902., powr<2>(_cse12) * _den3 * numtracer_kernels::ortho12::tr167(fenv) * numtracer_kernels::ortho12::tr255(fenv) * numtracer_kernels::ortho12::tr256(fenv), fma(2118., powr<2>(_cse12) * _den3 * numtracer_kernels::ortho12::tr169(fenv) * numtracer_kernels::ortho12::tr259(fenv) * numtracer_kernels::ortho12::tr260(fenv), fma(-18784., _cse10 * _cse12 * _den4 * numtracer_kernels::ortho12::tr50(fenv), fma(9392., _cse12 * powr<2>(_cse5) * _den4 * numtracer_kernels::ortho12::tr50(fenv), fma(9392., _cse12 * _cse6 * _den4 * numtracer_kernels::ortho12::tr50(fenv), fma(9392., _cse12 * _cse7 * _den4 * numtracer_kernels::ortho12::tr50(fenv), 0.)))))))))))));
        // clang-format on
      }
      { // subkernel 17
        const auto _cse1 = -powr<2>(cosP);
        const auto _cse2 = 1. + _cse1;
        const auto _cse3 = -p1m;
        const auto _cse4 = -cosP * p2m;
        const auto _cse5 = _cse3 + _cse4;
        const auto _cse6 = powr<2>(cosP) * powr<2>(p2m);
        const auto _cse7 = 2. * _cse2 * powr<2>(p2m);
        const auto _cse8 = -_cse2 * powr<2>(p2m);
        const auto _cse9 = _cse5 * cosP * p2m;
        const auto _cse10 = _cse8 + _cse9;
        const auto _cse11 = 2. * _cse10;
        const auto _cse12 = _cse11 + powr<2>(_cse5) + _cse6 + _cse7;
        const auto _cse13 = _cse2 * powr<2>(p2m);
        const auto _cse14 = _cse13 + _cse6;
        const auto _cse15 = _cse13 + powr<2>(_cse5);
        const auto _cse16 = 3. * powr<2>(_cse14);
        const auto _cse17 = -4. * powr<2>(_cse10);
        const auto _cse18 = -2. * _cse14 * _cse15;
        const auto _cse19 = 3. * powr<2>(_cse15);
        const auto _cse20 = _cse16 + _cse17 + _cse18 + _cse19;
        // clang-format off
        _acc += 0.0001085069444444444 * fma(-2313., _cse20 * _den3 * numtracer_kernels::ortho12::tr100(fenv), fma(-4842., _cse20 * _den3 * numtracer_kernels::ortho12::tr102(fenv), fma(-4338., _cse20 * _den3 * numtracer_kernels::ortho12::tr158(fenv) * numtracer_kernels::ortho12::tr237(fenv) * numtracer_kernels::ortho12::tr238(fenv), fma(4518., _cse20 * _den3 * numtracer_kernels::ortho12::tr160(fenv) * numtracer_kernels::ortho12::tr241(fenv) * numtracer_kernels::ortho12::tr242(fenv), fma(4986., _cse20 * _den3 * numtracer_kernels::ortho12::tr162(fenv) * numtracer_kernels::ortho12::tr245(fenv) * numtracer_kernels::ortho12::tr246(fenv), fma(-3328., powr<2>(_cse10) * _den2 * numtracer_kernels::ortho12::tr2(fenv), fma(3328., _cse14 * _cse15 * _den2 * numtracer_kernels::ortho12::tr2(fenv), fma(95184., powr<2>(_cse12) * _den1 * _den4 * numtracer_kernels::ortho12::tr53(fenv), fma(51048., powr<2>(_cse12) * _den1 * _den4 * numtracer_kernels::ortho12::tr57(fenv), fma(104688., powr<2>(_cse12) * _den1 * _den4 * numtracer_kernels::ortho12::tr58(fenv), fma(-4086., _cse20 * _den3 * numtracer_kernels::ortho12::tr93(fenv), 0.)))))))))));
        // clang-format on
      }
      { // subkernel 18
        const auto _cse1 = -powr<2>(cosP);
        const auto _cse2 = 1. + _cse1;
        const auto _cse3 = powr<2>(cosP) * powr<2>(p2m);
        const auto _cse4 = _cse2 * powr<2>(p2m);
        const auto _cse5 = _cse3 + _cse4;
        const auto _cse6 = -p1m;
        const auto _cse7 = -cosP * p2m;
        const auto _cse8 = _cse6 + _cse7;
        const auto _cse9 = _cse4 + powr<2>(_cse8);
        _acc += -0.6536458333333334 * fma(-3., powr<2>(_cse5) * _den2 * _den4 * numtracer_kernels::ortho12::tr97(fenv), fma(2., _cse5 * _cse9 * _den2 * _den4 * numtracer_kernels::ortho12::tr97(fenv), fma(-3., powr<2>(_cse9) * _den2 * _den4 * numtracer_kernels::ortho12::tr97(fenv), fma(4., powr<2>(_cse8) * _den2 * _den4 * numtracer_kernels::ortho12::tr97(fenv) * powr<2>(cosP) * powr<2>(p2m), fma(-8., _cse2 * _cse8 * _den2 * _den4 * numtracer_kernels::ortho12::tr97(fenv) * cosP * powr<3>(p2m), fma(4., powr<2>(_cse2) * _den2 * _den4 * numtracer_kernels::ortho12::tr97(fenv) * powr<4>(p2m), 0.))))));
      }
      return _acc;
    }

    static inline auto constant()
    {
      return 0.;
    }
  };
}
using numtracer_kernels::Ortho12_kernel;