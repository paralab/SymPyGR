
b_rhs3[pp] = kograd(0, beta0[pp]) + kograd(1, beta0[pp]) + kograd(2, beta0[pp]);
//--
b_rhs5[pp] = kograd(0, beta2[pp]) + kograd(1, beta2[pp]) + kograd(2, beta2[pp]);
//--
b_rhs4[pp] = kograd(0, beta1[pp]) + kograd(1, beta1[pp]) + kograd(2, beta1[pp]);
// Dendro: reduced ops: 15
// Dendro: }}} 
// Dendro vectorized code: {{{
//--
  double v0 = *(kograd_0_beta0+pp );
  double v1 = *(kograd_1_beta0+pp );
  double v2 = *(kograd_2_beta0+pp );
  double v3 = dadd(v2, v1);
  double v4 = dadd(v3, v0);
  b_rhs3[pp] = v4;
//--
  v0 = *(kograd_0_beta2+pp );
  v1 = *(kograd_1_beta2+pp );
  v2 = *(kograd_2_beta2+pp );
  v3 = dadd(v2, v1);
  v4 = dadd(v3, v0);
  b_rhs5[pp] = v4;
//--
  v0 = *(kograd_0_beta1+pp );
  v1 = *(kograd_1_beta1+pp );
  v2 = *(kograd_2_beta1+pp );
  v3 = dadd(v2, v1);
  v4 = dadd(v3, v0);
  b_rhs4[pp] = v4;
// Dendro vectorized code: }}}
// Dendro: {{{
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}}
// Dendro vectorized code: {{{
// Dendro vectorized code: }}}
// Dendro: {{{
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}}
// Dendro vectorized code: {{{
// Dendro vectorized code: }}}
// Dendro: {{{
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}}
// Dendro vectorized code: {{{
// Dendro vectorized code: }}}
// Dendro: {{{
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}}
// Dendro vectorized code: {{{
// Dendro vectorized code: }}}
// Dendro: {{{
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}}
// Dendro vectorized code: {{{
// Dendro vectorized code: }}}
// Dendro: {{{
// Dendro: original ops: 1263
// Dendro: printing temp variables
double DENDRO_0 = -DENDRO_872;
double DENDRO_1 = -DENDRO_858;
double DENDRO_8 = DENDRO_853 + DENDRO_854;
double DENDRO_10 = DENDRO_941*(DENDRO_969 - DENDRO_970) + DENDRO_968;

// Dendro: printing variables
//--
B_rhs2[pp] = -B2[pp]*eta - DENDRO_978*lambda[3] + DENDRO_980 + lambda[2]*(beta0[pp]*agrad(0, B2[pp]) + beta1[pp]*agrad(1, B2[pp]) + beta2[pp]*agrad(2, B2[pp]));
//--
At_rhs12[pp] = At1[pp]*DENDRO_15 + At2[pp]*DENDRO_12 + At3[pp]*DENDRO_16 - At4[pp]*DENDRO_18 + At5[pp]*DENDRO_13 + DENDRO_4*DENDRO_870 + DENDRO_56*(DENDRO_108*(DENDRO_124*DENDRO_46*(DENDRO_826 + DENDRO_847) + DENDRO_125*(DENDRO_689 - DENDRO_712) + DENDRO_125*(DENDRO_185*DENDRO_821 - DENDRO_397*DENDRO_821 + DENDRO_456*DENDRO_857) - DENDRO_126*DENDRO_308 + DENDRO_138*(DENDRO_197*DENDRO_219 + DENDRO_831) + DENDRO_138*(DENDRO_203 + DENDRO_833 + DENDRO_860) + DENDRO_149*(DENDRO_254*DENDRO_827 + DENDRO_728) - DENDRO_149*(DENDRO_262*DENDRO_855 - DENDRO_873) + DENDRO_149*(DENDRO_370*DENDRO_819 + DENDRO_873) - DENDRO_149*(DENDRO_143*DENDRO_874 + DENDRO_163*DENDRO_857 + DENDRO_729) - DENDRO_149*(DENDRO_152*DENDRO_874 + DENDRO_684*DENDRO_821 - DENDRO_706) + DENDRO_155*(DENDRO_1 + DENDRO_227*DENDRO_821 + DENDRO_71*DENDRO_857) - DENDRO_155*(DENDRO_205*DENDRO_370 + DENDRO_849 - DENDRO_871) + DENDRO_155*(DENDRO_456*DENDRO_855 + DENDRO_848 + DENDRO_871) - DENDRO_161*(DENDRO_0 + DENDRO_205*DENDRO_358) - DENDRO_161*(DENDRO_0 + DENDRO_684*DENDRO_819) - DENDRO_161*(DENDRO_584*DENDRO_827 - DENDRO_722) - DENDRO_161*(DENDRO_150*DENDRO_173 - DENDRO_724 + 0.25*DENDRO_852) - DENDRO_161*(DENDRO_158*DENDRO_821 - DENDRO_704 + DENDRO_869) + DENDRO_161*(-DENDRO_456*DENDRO_827 + DENDRO_572 + DENDRO_674) - DENDRO_229*(DENDRO_244*DENDRO_525 + DENDRO_275*DENDRO_522 + DENDRO_294*DENDRO_527 + DENDRO_652) - DENDRO_251*DENDRO_288 - DENDRO_260*DENDRO_270 + DENDRO_357*(DENDRO_358*DENDRO_819 - DENDRO_710) - DENDRO_357*(-DENDRO_150*DENDRO_821 + DENDRO_699 + DENDRO_730) + DENDRO_357*(DENDRO_262*DENDRO_827 + DENDRO_263*DENDRO_827 - DENDRO_451) - DENDRO_362*(DENDRO_240*DENDRO_827 - DENDRO_433) - DENDRO_640*DENDRO_837 + DENDRO_658 + DENDRO_659 + DENDRO_661 + DENDRO_662 + DENDRO_663 + DENDRO_664 + DENDRO_665 + DENDRO_666 + DENDRO_667 + DENDRO_668 + DENDRO_669 + DENDRO_670 + DENDRO_671 + DENDRO_714 + DENDRO_719 + DENDRO_720 + DENDRO_726) - 12*DENDRO_637 - 6.0*DENDRO_638*(DENDRO_152*DENDRO_50 - DENDRO_238*DENDRO_35 + DENDRO_240*DENDRO_46 - DENDRO_343*DENDRO_60*gt4[pp]) + DENDRO_640*DENDRO_804 + DENDRO_817*(DENDRO_60*(DENDRO_517 + DENDRO_640*DENDRO_65) + DENDRO_644) + DENDRO_846*(DENDRO_645 - DENDRO_646 - DENDRO_647 + DENDRO_650)) + DENDRO_6*DENDRO_870 - alpha[pp]*(-At4[pp]*K[pp] + DENDRO_34*DENDRO_844 + DENDRO_841*DENDRO_862 + DENDRO_845*DENDRO_863) + beta0[pp]*agrad(0, At4[pp]) + beta1[pp]*agrad(1, At4[pp]) + beta2[pp]*agrad(2, At4[pp]);
//--
At_rhs00[pp] = At0[pp]*DENDRO_3 + DENDRO_11*DENDRO_26 - DENDRO_24*DENDRO_4 - DENDRO_24*DENDRO_6 + DENDRO_25*DENDRO_9 + DENDRO_56*(-DENDRO_108*(DENDRO_139*(DENDRO_209 + DENDRO_210) + DENDRO_139*(DENDRO_197*DENDRO_76 + DENDRO_215) + DENDRO_149*(DENDRO_213 + 0.5*DENDRO_214) + DENDRO_149*(DENDRO_216 - DENDRO_218) + DENDRO_149*(DENDRO_222 + DENDRO_225) + DENDRO_149*(DENDRO_173*DENDRO_70 + DENDRO_203) + DENDRO_149*(DENDRO_179*DENDRO_77 + DENDRO_204) - DENDRO_154*(DENDRO_207 + DENDRO_208) - DENDRO_154*(DENDRO_205*DENDRO_76 + DENDRO_206) - DENDRO_155*(1.0*DENDRO_208 + DENDRO_211) - DENDRO_155*(DENDRO_205*DENDRO_227 - DENDRO_228*DENDRO_93) + DENDRO_161*(1.0*DENDRO_210 + DENDRO_212) + DENDRO_161*(1.0*DENDRO_152*DENDRO_93 + DENDRO_197*DENDRO_227) - DENDRO_168 - DENDRO_169*DENDRO_174 - DENDRO_175*DENDRO_180 - DENDRO_182*(DENDRO_150 + DENDRO_166) - DENDRO_183*DENDRO_184*DENDRO_50 - DENDRO_193*(DENDRO_185 + DENDRO_188) - DENDRO_197*DENDRO_200 - DENDRO_201*DENDRO_93 + DENDRO_229*(DENDRO_231 + DENDRO_232*DENDRO_233 + DENDRO_234*DENDRO_80 + DENDRO_235*DENDRO_236) + DENDRO_270*DENDRO_69 + DENDRO_288*DENDRO_76 + DENDRO_308*DENDRO_72 + DENDRO_336*gt0[pp]) - 12*DENDRO_57 + DENDRO_805*gt0[pp] + 12*DENDRO_81 + DENDRO_83*(DENDRO_89 - DENDRO_93) + 12*DENDRO_94*(DENDRO_100 + DENDRO_107 - DENDRO_96 - DENDRO_98)) + alpha[pp]*(At0[pp]*K[pp] - DENDRO_34*(DENDRO_36 + DENDRO_38 - DENDRO_41) - DENDRO_42*(-At0[pp]*DENDRO_50 + DENDRO_43 + DENDRO_48) + DENDRO_51*(-DENDRO_52 + DENDRO_53 + DENDRO_55)) + beta0[pp]*agrad(0, At0[pp]) + beta1[pp]*agrad(1, At0[pp]) + beta2[pp]*agrad(2, At0[pp]);
//--
At_rhs22[pp] = -At5[pp]*DENDRO_18 + At5[pp]*DENDRO_22 - At5[pp]*DENDRO_5 + DENDRO_15*DENDRO_26 + DENDRO_16*DENDRO_861 - DENDRO_56*(DENDRO_108*(-DENDRO_134*DENDRO_205*DENDRO_374 + DENDRO_134*DENDRO_308 + DENDRO_149*(DENDRO_398 - 1.0*DENDRO_400) + DENDRO_149*(0.25*DENDRO_875 + 1.0*DENDRO_876) - DENDRO_154*(DENDRO_877 + DENDRO_878) - DENDRO_154*(DENDRO_134*DENDRO_173 + DENDRO_76*DENDRO_857) - DENDRO_155*(DENDRO_394 - 1.0*DENDRO_396) - DENDRO_155*(0.25*DENDRO_877 + 1.0*DENDRO_878) + DENDRO_155*(-DENDRO_402*DENDRO_857 + DENDRO_856) + DENDRO_161*(DENDRO_202*DENDRO_821 + 0.5*DENDRO_859) + DENDRO_161*(DENDRO_205*DENDRO_250 + DENDRO_871) + DENDRO_161*(DENDRO_219*DENDRO_819 + DENDRO_850) + DENDRO_161*(DENDRO_78*DENDRO_868 + 0.25*DENDRO_825) - DENDRO_173*DENDRO_377 - DENDRO_175*DENDRO_240*DENDRO_819 + DENDRO_229*(DENDRO_232*DENDRO_281 + DENDRO_234*DENDRO_255 + DENDRO_235*DENDRO_301 + DENDRO_403) + DENDRO_240*DENDRO_270 + DENDRO_247*DENDRO_288 + DENDRO_336*gt5[pp] + DENDRO_362*(DENDRO_875 + DENDRO_876) + DENDRO_362*(DENDRO_134*DENDRO_821 + DENDRO_152*DENDRO_857) - DENDRO_372 - DENDRO_376*DENDRO_54*DENDRO_855 - DENDRO_378*DENDRO_821 - DENDRO_379*(DENDRO_188 + DENDRO_408) - DENDRO_380*(DENDRO_253 + DENDRO_405)) + 12*DENDRO_337 - 12*DENDRO_339 - DENDRO_805*gt5[pp] + 12*DENDRO_82*(DENDRO_345 - DENDRO_346 + DENDRO_347 - DENDRO_352) - DENDRO_864*(DENDRO_344 - DENDRO_857)) - alpha[pp]*(At5[pp]*DENDRO_328*DENDRO_845 - At5[pp]*K[pp] + DENDRO_51*DENDRO_844 + DENDRO_841*DENDRO_863) + beta0[pp]*agrad(0, At5[pp]) + beta1[pp]*agrad(1, At5[pp]) + beta2[pp]*agrad(2, At5[pp]);
//--
At_rhs02[pp] = At0[pp]*DENDRO_15 + At1[pp]*DENDRO_16 - At2[pp]*DENDRO_5 + At4[pp]*DENDRO_9 + At5[pp]*DENDRO_11 + DENDRO_2*DENDRO_838 + DENDRO_56*(DENDRO_108*(DENDRO_124*DENDRO_40*(DENDRO_823 + DENDRO_852) + DENDRO_125*(DENDRO_587 - DENDRO_610) + DENDRO_125*(DENDRO_173*DENDRO_185 + DENDRO_77*DENDRO_857 - DENDRO_856) - DENDRO_129*DENDRO_270 + DENDRO_138*(DENDRO_202*DENDRO_205 + 1.0*DENDRO_206) + DENDRO_138*(DENDRO_184*DENDRO_77 + DENDRO_184*DENDRO_78 + DENDRO_211) - DENDRO_149*(DENDRO_1 + DENDRO_70*DENDRO_857 + 0.25*DENDRO_859) - DENDRO_149*(-DENDRO_605 + DENDRO_606 + DENDRO_632) + DENDRO_149*(DENDRO_397*DENDRO_819 + DENDRO_849 - DENDRO_850) - DENDRO_149*(DENDRO_578*DENDRO_855 + DENDRO_848 + DENDRO_850) + DENDRO_154*(DENDRO_134*DENDRO_184 + DENDRO_174) + DENDRO_155*(-DENDRO_205*DENDRO_397 + DENDRO_8) + DENDRO_155*(DENDRO_77*DENDRO_855 + DENDRO_8) + DENDRO_155*(DENDRO_173*DENDRO_227 - DENDRO_184*DENDRO_220 + DENDRO_73*DENDRO_857) - DENDRO_161*(0.25*DENDRO_214 + DENDRO_851) - DENDRO_161*(DENDRO_829 + DENDRO_851) - DENDRO_161*(DENDRO_834 + DENDRO_860) - DENDRO_161*(DENDRO_184*DENDRO_584 + DENDRO_832) - DENDRO_186*DENDRO_288 - DENDRO_229*(DENDRO_136*DENDRO_525 + DENDRO_192*DENDRO_522 + DENDRO_289*DENDRO_527 + DENDRO_521) - DENDRO_308*DENDRO_74 + DENDRO_357*(DENDRO_197*DENDRO_250 + 0.25*DENDRO_820) - DENDRO_510*DENDRO_837 + DENDRO_533 + DENDRO_535 + DENDRO_537 + DENDRO_538 + DENDRO_540 + DENDRO_542 + DENDRO_544 + DENDRO_545 + DENDRO_546 + DENDRO_547 + DENDRO_548 + DENDRO_549 + DENDRO_550 + DENDRO_614 + DENDRO_621 + DENDRO_626 + DENDRO_628 + DENDRO_629 + DENDRO_634 + DENDRO_635 + DENDRO_636 - DENDRO_818*(DENDRO_825 + DENDRO_847)) - 12*DENDRO_503 + 6.0*DENDRO_504 + DENDRO_510*DENDRO_804 - DENDRO_816*(DENDRO_514 + DENDRO_515 - DENDRO_516 - DENDRO_519) - DENDRO_846*(DENDRO_506 + DENDRO_507 - DENDRO_508 - DENDRO_512)) + DENDRO_6*DENDRO_838 - alpha[pp]*(-At2[pp]*K[pp] + DENDRO_34*DENDRO_841 + DENDRO_42*DENDRO_844 + DENDRO_51*DENDRO_845) + beta0[pp]*agrad(0, At2[pp]) + beta1[pp]*agrad(1, At2[pp]) + beta2[pp]*agrad(2, At2[pp]);
//--
B_rhs1[pp] = -B1[pp]*eta + DENDRO_10 - DENDRO_16*DENDRO_958 + DENDRO_266*DENDRO_946 - DENDRO_4*DENDRO_957 - DENDRO_9*DENDRO_959 - DENDRO_914*DENDRO_940 + DENDRO_914*DENDRO_973 - DENDRO_924*DENDRO_948 + DENDRO_924*DENDRO_974 - DENDRO_941*(DENDRO_914*DENDRO_943 + DENDRO_976) - DENDRO_941*(DENDRO_924*DENDRO_952 + DENDRO_975) + DENDRO_957*DENDRO_960 + lambda[2]*(beta0[pp]*agrad(0, B1[pp]) + beta1[pp]*agrad(1, B1[pp]) + beta2[pp]*agrad(2, B1[pp])) - lambda[3]*(DENDRO_962 + DENDRO_963 + DENDRO_964);
//--
At_rhs01[pp] = At0[pp]*DENDRO_12 - At1[pp]*DENDRO_7 + At2[pp]*DENDRO_13 + At3[pp]*DENDRO_9 + At4[pp]*DENDRO_11 + DENDRO_2*DENDRO_806 + DENDRO_4*DENDRO_806 + DENDRO_56*(DENDRO_108*(DENDRO_124*DENDRO_54*(DENDRO_824 + DENDRO_825 + DENDRO_826) + DENDRO_125*(-DENDRO_381 + DENDRO_593 + DENDRO_691) - DENDRO_131*DENDRO_288 + DENDRO_138*(DENDRO_184*DENDRO_70 + DENDRO_184*DENDRO_71 + DENDRO_212) + DENDRO_138*(DENDRO_197*DENDRO_202 + 0.5*DENDRO_215 + DENDRO_456*DENDRO_93) - DENDRO_139*(DENDRO_140*DENDRO_184 + DENDRO_180) - DENDRO_149*(DENDRO_578*DENDRO_827 - DENDRO_777) - DENDRO_149*(-DENDRO_694 + DENDRO_795 + DENDRO_796) + DENDRO_149*(DENDRO_173*DENDRO_465 - DENDRO_179*DENDRO_684 + DENDRO_724) + DENDRO_155*(DENDRO_833 + DENDRO_834) + DENDRO_155*(DENDRO_184*DENDRO_456 + DENDRO_832) + DENDRO_155*(DENDRO_225 + DENDRO_798 + DENDRO_799) + DENDRO_155*(DENDRO_829 + DENDRO_830 + DENDRO_831) - DENDRO_161*(DENDRO_70*DENDRO_827 - DENDRO_774) - DENDRO_161*(DENDRO_133*DENDRO_836 + DENDRO_263*DENDRO_93 - DENDRO_793) - DENDRO_161*(DENDRO_143*DENDRO_836 + DENDRO_152*DENDRO_836 + DENDRO_262*DENDRO_93) + DENDRO_161*(-DENDRO_158*DENDRO_179 + DENDRO_167*DENDRO_184 + DENDRO_791) - DENDRO_164*DENDRO_270 - DENDRO_229*(DENDRO_146*DENDRO_525 + DENDRO_221*DENDRO_522 + DENDRO_296*DENDRO_527 + DENDRO_744) - DENDRO_308*DENDRO_67 - DENDRO_357*(DENDRO_766 - DENDRO_767) - DENDRO_357*(-DENDRO_150*DENDRO_179 + DENDRO_782 + DENDRO_828) - DENDRO_357*(-DENDRO_197*DENDRO_358 + DENDRO_786 + DENDRO_835) - DENDRO_733*DENDRO_837 + DENDRO_750 + DENDRO_751 + DENDRO_752 + DENDRO_753 + DENDRO_754 + DENDRO_755 + DENDRO_756 + DENDRO_757 + DENDRO_758 + DENDRO_759 + DENDRO_760 + DENDRO_761 + DENDRO_762 + DENDRO_764 + DENDRO_765 + DENDRO_769 + DENDRO_771 + DENDRO_775 + DENDRO_779 - DENDRO_818*(DENDRO_140*DENDRO_173 + DENDRO_823) - DENDRO_818*(DENDRO_197*DENDRO_240 + DENDRO_205*DENDRO_238 + DENDRO_820)) - 12*DENDRO_731 + 6.0*DENDRO_732*(DENDRO_140*DENDRO_37 - DENDRO_143*DENDRO_54 - DENDRO_46*DENDRO_69 + DENDRO_60*DENDRO_87*gt1[pp]) + DENDRO_733*DENDRO_804 + DENDRO_816*(DENDRO_738 - DENDRO_739 - DENDRO_740 + DENDRO_742) + DENDRO_817*(DENDRO_60*(DENDRO_509 + DENDRO_65*DENDRO_733) + DENDRO_737)) - alpha[pp]*(-At1[pp]*K[pp] + DENDRO_34*DENDRO_808 + DENDRO_42*DENDRO_811 + DENDRO_51*DENDRO_815) + beta0[pp]*agrad(0, At1[pp]) + beta1[pp]*agrad(1, At1[pp]) + beta2[pp]*agrad(2, At1[pp]);
//--
At_rhs11[pp] = -At3[pp]*DENDRO_18 + At3[pp]*DENDRO_19 - At3[pp]*DENDRO_7 + DENDRO_12*DENDRO_25 + DENDRO_13*DENDRO_861 + DENDRO_56*(DENDRO_108*(DENDRO_139*(DENDRO_447 - DENDRO_866) - DENDRO_139*(DENDRO_140*DENDRO_179 - DENDRO_450) - DENDRO_139*(DENDRO_197*DENDRO_238 - DENDRO_459) - DENDRO_140*DENDRO_308 + DENDRO_149*(DENDRO_451 - 1.0*DENDRO_865) - DENDRO_149*(DENDRO_466 - 1.0*DENDRO_468) - DENDRO_149*(DENDRO_469*DENDRO_819 + DENDRO_471) + DENDRO_155*(DENDRO_462 - DENDRO_464) + DENDRO_155*(DENDRO_156*DENDRO_821 + DENDRO_179*DENDRO_456) + DENDRO_155*(DENDRO_456*DENDRO_819 + DENDRO_867) + DENDRO_155*(DENDRO_71*DENDRO_868 + DENDRO_869) + DENDRO_161*(DENDRO_452 - 1.0*DENDRO_866) + DENDRO_161*(DENDRO_473 + DENDRO_828) + DENDRO_161*(DENDRO_475 + DENDRO_835) + DENDRO_179*DENDRO_440 + DENDRO_197*DENDRO_443 - DENDRO_229*(DENDRO_232*DENDRO_284 + DENDRO_234*DENDRO_266 + DENDRO_235*DENDRO_304 + DENDRO_476) - DENDRO_238*DENDRO_288 - DENDRO_257*DENDRO_270 - DENDRO_336*gt3[pp] + DENDRO_362*(DENDRO_445 - DENDRO_865) - DENDRO_362*(DENDRO_140*DENDRO_821 - DENDRO_458) - DENDRO_362*(DENDRO_238*DENDRO_819 - DENDRO_449) + DENDRO_40*DENDRO_438*DENDRO_827 + DENDRO_437 + DENDRO_439*(DENDRO_253 + DENDRO_369) + DENDRO_442*DENDRO_821 + DENDRO_444*(DENDRO_166 + DENDRO_477)) - 12*DENDRO_419 + 12*DENDRO_58*(DENDRO_429 + DENDRO_60*(DENDRO_424 + DENDRO_425*DENDRO_65)) + DENDRO_805*gt3[pp] + DENDRO_83*(DENDRO_284 - DENDRO_422*(-DENDRO_84 + DENDRO_85 + DENDRO_86)) + DENDRO_864*(DENDRO_304 - DENDRO_422*(-DENDRO_103 + DENDRO_341 + DENDRO_342))) - alpha[pp]*(-At3[pp]*K[pp] + DENDRO_34*DENDRO_811 + DENDRO_808*DENDRO_862 + DENDRO_815*DENDRO_863) + beta0[pp]*agrad(0, At3[pp]) + beta1[pp]*agrad(1, At3[pp]) + beta2[pp]*agrad(2, At3[pp]);
//--
K_rhs[pp] = -DENDRO_879*chi[pp]*(-DENDRO_337 + DENDRO_638*(DENDRO_338*DENDRO_882 + DENDRO_413) + DENDRO_82*(DENDRO_33*DENDRO_404 + DENDRO_33*DENDRO_407 + DENDRO_33*DENDRO_410 - DENDRO_60*(DENDRO_348 - DENDRO_351)) + DENDRO_880*(DENDRO_338*DENDRO_881 + DENDRO_412)) - DENDRO_883*chi[pp]*(-DENDRO_419 + DENDRO_58*(DENDRO_33*DENDRO_480 + DENDRO_33*DENDRO_481 + DENDRO_428 - DENDRO_60*(DENDRO_423 - DENDRO_427)) + DENDRO_638*(DENDRO_421*DENDRO_882 + DENDRO_482) + DENDRO_732*(DENDRO_421*DENDRO_884 + DENDRO_479)) - DENDRO_885*chi[pp]*(-DENDRO_57 + DENDRO_732*(DENDRO_499 + DENDRO_59*DENDRO_884) + DENDRO_880*(DENDRO_500 + DENDRO_59*DENDRO_881) + DENDRO_94*(DENDRO_100 + DENDRO_33*DENDRO_501 + DENDRO_33*DENDRO_502 - DENDRO_60*(DENDRO_101 - DENDRO_106))) + DENDRO_886*chi[pp]*(-DENDRO_637 + DENDRO_887*DENDRO_94*(DENDRO_656 + DENDRO_882*gt4[pp]) + DENDRO_888*(DENDRO_33*DENDRO_653 + DENDRO_33*DENDRO_654 - DENDRO_60*(DENDRO_61 - DENDRO_649) + DENDRO_645) + DENDRO_889*(DENDRO_33*DENDRO_655 - DENDRO_60*(DENDRO_62 - DENDRO_641) + DENDRO_642 + DENDRO_643)) - DENDRO_890*chi[pp]*(-DENDRO_503 + DENDRO_58*DENDRO_887*(DENDRO_526 + DENDRO_881*gt2[pp]) + DENDRO_888*(DENDRO_33*DENDRO_523 + DENDRO_33*DENDRO_524 + DENDRO_508 + DENDRO_60*(DENDRO_511 - DENDRO_63)) + DENDRO_891*(DENDRO_33*DENDRO_528 + DENDRO_33*DENDRO_529 + DENDRO_516 + DENDRO_60*(DENDRO_518 - DENDRO_62))) + DENDRO_892*chi[pp]*(-DENDRO_731 + DENDRO_82*DENDRO_887*(DENDRO_745 + DENDRO_884*gt1[pp]) + DENDRO_889*(DENDRO_33*DENDRO_746 - DENDRO_60*(DENDRO_63 - DENDRO_734) + DENDRO_735 + DENDRO_736) + DENDRO_891*(DENDRO_33*DENDRO_747 + DENDRO_33*DENDRO_748 - DENDRO_60*(DENDRO_61 - DENDRO_741) + DENDRO_738)) + (1.0L/3.0L)*alpha[pp]*(At0[pp]*DENDRO_893*DENDRO_898 + At1[pp]*DENDRO_903*DENDRO_914 + At2[pp]*DENDRO_903*DENDRO_904 + At3[pp]*DENDRO_893*DENDRO_901 + At4[pp]*DENDRO_903*DENDRO_924 + At5[pp]*DENDRO_893*DENDRO_902 + pow(K[pp], 2)) + beta0[pp]*agrad(0, K[pp]) + beta1[pp]*agrad(1, K[pp]) + beta2[pp]*agrad(2, K[pp]);
//--
B_rhs0[pp] = -B0[pp]*eta - DENDRO_925*lambda[3] + DENDRO_961 + lambda[2]*(beta0[pp]*agrad(0, B0[pp]) + beta1[pp]*agrad(1, B0[pp]) + beta2[pp]*agrad(2, B0[pp]));
//--
Gt_rhs1[pp] = DENDRO_10 + DENDRO_124*DENDRO_16*DENDRO_287 + DENDRO_124*DENDRO_307*DENDRO_9 + DENDRO_4*DENDRO_977 - DENDRO_827*DENDRO_946 + DENDRO_940*DENDRO_972 + DENDRO_941*(DENDRO_943*DENDRO_972 - DENDRO_976) + DENDRO_941*(DENDRO_952*DENDRO_971 - DENDRO_975) + DENDRO_948*DENDRO_971 - DENDRO_960*DENDRO_977 - DENDRO_971*DENDRO_974 - DENDRO_972*DENDRO_973;
//--
Gt_rhs0[pp] = DENDRO_961;
//--
Gt_rhs2[pp] = DENDRO_980;
// Dendro: reduced ops: 1255
// Dendro: }}} 
// Dendro vectorized code: {{{
  double v0 = DENDRO_872;
  double v1 = dmul(v0, negone);
  double DENDRO_0 = v1;
  v0 = DENDRO_858;
  v1 = dmul(v0, negone);
  double DENDRO_1 = v1;
  v0 = DENDRO_853;
  v1 = DENDRO_854;
  double v2 = dadd(v1, v0);
  double DENDRO_8 = v2;
  v0 = DENDRO_968;
  v1 = DENDRO_941;
  v2 = DENDRO_969;
  double v3 = DENDRO_970;
  double v4 = dmul(v3, negone);
  double v5 = dadd(v4, v2);
  double v6 = dmul(v5, v1);
  double v7 = dadd(v6, v0);
  double DENDRO_10 = v7;
//--
  v0 = DENDRO_980;
  v1 = lambda[2];
  v2 = beta0[pp];
  v3 = *(agrad_0_B2+pp );
  v4 = dmul(v3, v2);
  v5 = beta1[pp];
  v6 = *(agrad_1_B2+pp );
  v7 = dmul(v6, v5);
  double v8 = beta2[pp];
  double v9 = *(agrad_2_B2+pp );
  double v10 = dmul(v9, v8);
  double v11 = dadd(v10, v7);
  double v12 = dadd(v11, v4);
  double v13 = dmul(v12, v1);
  double v14 = B2[pp];
  double v15 = eta;
  double v16 = dmul(v15, v14);
  double v17 = dmul(v16, negone);
  double v18 = DENDRO_978;
  double v19 = lambda[3];
  double v20 = dmul(v19, v18);
  double v21 = dmul(v20, negone);
  double v22 = dadd(v21, v17);
  double v23 = dadd(v22, v13);
  double v24 = dadd(v23, v0);
  B_rhs2[pp] = v24;
//--
  v0 = At1[pp];
  v1 = DENDRO_15;
  v2 = dmul(v1, v0);
  v3 = At2[pp];
  v4 = DENDRO_12;
  v5 = dmul(v4, v3);
  v6 = At3[pp];
  v7 = DENDRO_16;
  v8 = dmul(v7, v6);
  v9 = At5[pp];
  v10 = DENDRO_13;
  v11 = dmul(v10, v9);
  v12 = DENDRO_4;
  v13 = DENDRO_870;
  v14 = dmul(v13, v12);
  v15 = DENDRO_56;
  v16 = DENDRO_108;
  v17 = DENDRO_658;
  v18 = DENDRO_659;
  v19 = DENDRO_661;
  v20 = DENDRO_662;
  v21 = DENDRO_663;
  v22 = DENDRO_664;
  v23 = DENDRO_665;
  v24 = DENDRO_666;
  double v25 = DENDRO_667;
  double v26 = DENDRO_668;
  double v27 = DENDRO_669;
  double v28 = DENDRO_670;
  double v29 = DENDRO_671;
  double v30 = DENDRO_714;
  double v31 = DENDRO_719;
  double v32 = DENDRO_720;
  double v33 = DENDRO_726;
  double v34 = DENDRO_125;
  double v35 = DENDRO_689;
  double v36 = DENDRO_712;
  double v37 = dmul(v36, negone);
  double v38 = dadd(v37, v35);
  double v39 = dmul(v38, v34);
  double v40 = DENDRO_125;
  double v41 = DENDRO_185;
  double v42 = DENDRO_821;
  double v43 = dmul(v42, v41);
  double v44 = DENDRO_456;
  double v45 = DENDRO_857;
  double v46 = dmul(v45, v44);
  double v47 = DENDRO_397;
  double v48 = DENDRO_821;
  double v49 = dmul(v48, v47);
  double v50 = dmul(v49, negone);
  double v51 = dadd(v50, v46);
  double v52 = dadd(v51, v43);
  double v53 = dmul(v52, v40);
  double v54 = DENDRO_138;
  double v55 = DENDRO_831;
  double v56 = DENDRO_197;
  double v57 = DENDRO_219;
  double v58 = dmul(v57, v56);
  double v59 = dadd(v58, v55);
  double v60 = dmul(v59, v54);
  double v61 = DENDRO_138;
  double v62 = DENDRO_203;
  double v63 = DENDRO_833;
  double v64 = DENDRO_860;
  double v65 = dadd(v64, v63);
  double v66 = dadd(v65, v62);
  double v67 = dmul(v66, v61);
  double v68 = DENDRO_149;
  double v69 = DENDRO_728;
  double v70 = DENDRO_254;
  double v71 = DENDRO_827;
  double v72 = dmul(v71, v70);
  double v73 = dadd(v72, v69);
  double v74 = dmul(v73, v68);
  double v75 = DENDRO_149;
  double v76 = DENDRO_873;
  double v77 = DENDRO_370;
  double v78 = DENDRO_819;
  double v79 = dmul(v78, v77);
  double v80 = dadd(v79, v76);
  double v81 = dmul(v80, v75);
  double v82 = DENDRO_155;
  double v83 = DENDRO_1;
  double v84 = DENDRO_227;
  double v85 = DENDRO_821;
  double v86 = dmul(v85, v84);
  double v87 = DENDRO_71;
  double v88 = DENDRO_857;
  double v89 = dmul(v88, v87);
  double v90 = dadd(v89, v86);
  double v91 = dadd(v90, v83);
  double v92 = dmul(v91, v82);
  double v93 = DENDRO_155;
  double v94 = DENDRO_848;
  double v95 = DENDRO_871;
  double v96 = DENDRO_456;
  double v97 = DENDRO_855;
  double v98 = dmul(v97, v96);
  double v99 = dadd(v98, v95);
  double v100 = dadd(v99, v94);
  double v101 = dmul(v100, v93);
  double v102 = DENDRO_161;
  double v103 = DENDRO_572;
  double v104 = DENDRO_674;
  double v105 = DENDRO_456;
  double v106 = DENDRO_827;
  double v107 = dmul(v106, v105);
  double v108 = dmul(v107, negone);
  double v109 = dadd(v108, v104);
  double v110 = dadd(v109, v103);
  double v111 = dmul(v110, v102);
  double v112 = DENDRO_357;
  double v113 = DENDRO_710;
  double v114 = dmul(v113, negone);
  double v115 = DENDRO_358;
  double v116 = DENDRO_819;
  double v117 = dmul(v116, v115);
  double v118 = dadd(v117, v114);
  double v119 = dmul(v118, v112);
  double v120 = DENDRO_357;
  double v121 = DENDRO_451;
  double v122 = dmul(v121, negone);
  double v123 = DENDRO_262;
  double v124 = DENDRO_827;
  double v125 = dmul(v124, v123);
  double v126 = DENDRO_263;
  double v127 = DENDRO_827;
  double v128 = dmul(v127, v126);
  double v129 = dadd(v128, v125);
  double v130 = dadd(v129, v122);
  double v131 = dmul(v130, v120);
  double v132 = DENDRO_126;
  double v133 = DENDRO_308;
  double v134 = dmul(v133, v132);
  double v135 = dmul(v134, negone);
  double v136 = DENDRO_149;
  double v137 = DENDRO_873;
  double v138 = dmul(v137, negone);
  double v139 = DENDRO_262;
  double v140 = DENDRO_855;
  double v141 = dmul(v140, v139);
  double v142 = dadd(v141, v138);
  double v143 = dmul(v142, v136);
  double v144 = dmul(v143, negone);
  double v145 = DENDRO_149;
  double v146 = DENDRO_729;
  double v147 = DENDRO_143;
  double v148 = DENDRO_874;
  double v149 = dmul(v148, v147);
  double v150 = DENDRO_163;
  double v151 = DENDRO_857;
  double v152 = dmul(v151, v150);
  double v153 = dadd(v152, v149);
  double v154 = dadd(v153, v146);
  double v155 = dmul(v154, v145);
  double v156 = dmul(v155, negone);
  double v157 = DENDRO_149;
  double v158 = DENDRO_706;
  double v159 = dmul(v158, negone);
  double v160 = DENDRO_152;
  double v161 = DENDRO_874;
  double v162 = dmul(v161, v160);
  double v163 = DENDRO_684;
  double v164 = DENDRO_821;
  double v165 = dmul(v164, v163);
  double v166 = dadd(v165, v162);
  double v167 = dadd(v166, v159);
  double v168 = dmul(v167, v157);
  double v169 = dmul(v168, negone);
  double v170 = DENDRO_155;
  double v171 = DENDRO_849;
  double v172 = DENDRO_871;
  double v173 = dmul(v172, negone);
  double v174 = DENDRO_205;
  double v175 = DENDRO_370;
  double v176 = dmul(v175, v174);
  double v177 = dadd(v176, v173);
  double v178 = dadd(v177, v171);
  double v179 = dmul(v178, v170);
  double v180 = dmul(v179, negone);
  double v181 = DENDRO_161;
  double v182 = DENDRO_0;
  double v183 = DENDRO_205;
  double v184 = DENDRO_358;
  double v185 = dmul(v184, v183);
  double v186 = dadd(v185, v182);
  double v187 = dmul(v186, v181);
  double v188 = dmul(v187, negone);
  double v189 = DENDRO_161;
  double v190 = DENDRO_0;
  double v191 = DENDRO_684;
  double v192 = DENDRO_819;
  double v193 = dmul(v192, v191);
  double v194 = dadd(v193, v190);
  double v195 = dmul(v194, v189);
  double v196 = dmul(v195, negone);
  double v197 = DENDRO_161;
  double v198 = DENDRO_722;
  double v199 = dmul(v198, negone);
  double v200 = DENDRO_584;
  double v201 = DENDRO_827;
  double v202 = dmul(v201, v200);
  double v203 = dadd(v202, v199);
  double v204 = dmul(v203, v197);
  double v205 = dmul(v204, negone);
  double v206 = DENDRO_161;
  double v207 = DENDRO_869;
  double v208 = DENDRO_704;
  double v209 = dmul(v208, negone);
  double v210 = DENDRO_158;
  double v211 = DENDRO_821;
  double v212 = dmul(v211, v210);
  double v213 = dadd(v212, v209);
  double v214 = dadd(v213, v207);
  double v215 = dmul(v214, v206);
  double v216 = dmul(v215, negone);
  double v217 = DENDRO_161;
  double v218 = DENDRO_724;
  double v219 = dmul(v218, negone);
  double v220 = 0.250000000000000;
  double v221 = DENDRO_852;
  double v222 = dmul(v221, v220);
  double v223 = DENDRO_150;
  double v224 = DENDRO_173;
  double v225 = dmul(v224, v223);
  double v226 = dadd(v225, v222);
  double v227 = dadd(v226, v219);
  double v228 = dmul(v227, v217);
  double v229 = dmul(v228, negone);
  double v230 = DENDRO_229;
  double v231 = DENDRO_652;
  double v232 = DENDRO_244;
  double v233 = DENDRO_525;
  double v234 = dmul(v233, v232);
  double v235 = DENDRO_275;
  double v236 = DENDRO_522;
  double v237 = dmul(v236, v235);
  double v238 = DENDRO_294;
  double v239 = DENDRO_527;
  double v240 = dmul(v239, v238);
  double v241 = dadd(v240, v237);
  double v242 = dadd(v241, v234);
  double v243 = dadd(v242, v231);
  double v244 = dmul(v243, v230);
  double v245 = dmul(v244, negone);
  double v246 = DENDRO_251;
  double v247 = DENDRO_288;
  double v248 = dmul(v247, v246);
  double v249 = dmul(v248, negone);
  double v250 = DENDRO_260;
  double v251 = DENDRO_270;
  double v252 = dmul(v251, v250);
  double v253 = dmul(v252, negone);
  double v254 = DENDRO_357;
  double v255 = DENDRO_699;
  double v256 = DENDRO_730;
  double v257 = DENDRO_150;
  double v258 = DENDRO_821;
  double v259 = dmul(v258, v257);
  double v260 = dmul(v259, negone);
  double v261 = dadd(v260, v256);
  double v262 = dadd(v261, v255);
  double v263 = dmul(v262, v254);
  double v264 = dmul(v263, negone);
  double v265 = DENDRO_362;
  double v266 = DENDRO_433;
  double v267 = dmul(v266, negone);
  double v268 = DENDRO_240;
  double v269 = DENDRO_827;
  double v270 = dmul(v269, v268);
  double v271 = dadd(v270, v267);
  double v272 = dmul(v271, v265);
  double v273 = dmul(v272, negone);
  double v274 = DENDRO_640;
  double v275 = DENDRO_837;
  double v276 = dmul(v275, v274);
  double v277 = dmul(v276, negone);
  double v278 = DENDRO_124;
  double v279 = DENDRO_46;
  double v280 = DENDRO_826;
  double v281 = DENDRO_847;
  double v282 = dadd(v281, v280);
  double v283 = dmul(v282, v279);
  double v284 = dmul(v283, v278);
  double v285 = dadd(v284, v277);
  double v286 = dadd(v285, v273);
  double v287 = dadd(v286, v264);
  double v288 = dadd(v287, v253);
  double v289 = dadd(v288, v249);
  double v290 = dadd(v289, v245);
  double v291 = dadd(v290, v229);
  double v292 = dadd(v291, v216);
  double v293 = dadd(v292, v205);
  double v294 = dadd(v293, v196);
  double v295 = dadd(v294, v188);
  double v296 = dadd(v295, v180);
  double v297 = dadd(v296, v169);
  double v298 = dadd(v297, v156);
  double v299 = dadd(v298, v144);
  double v300 = dadd(v299, v135);
  double v301 = dadd(v300, v131);
  double v302 = dadd(v301, v119);
  double v303 = dadd(v302, v111);
  double v304 = dadd(v303, v101);
  double v305 = dadd(v304, v92);
  double v306 = dadd(v305, v81);
  double v307 = dadd(v306, v74);
  double v308 = dadd(v307, v67);
  double v309 = dadd(v308, v60);
  double v310 = dadd(v309, v53);
  double v311 = dadd(v310, v39);
  double v312 = dadd(v311, v33);
  double v313 = dadd(v312, v32);
  double v314 = dadd(v313, v31);
  double v315 = dadd(v314, v30);
  double v316 = dadd(v315, v29);
  double v317 = dadd(v316, v28);
  double v318 = dadd(v317, v27);
  double v319 = dadd(v318, v26);
  double v320 = dadd(v319, v25);
  double v321 = dadd(v320, v24);
  double v322 = dadd(v321, v23);
  double v323 = dadd(v322, v22);
  double v324 = dadd(v323, v21);
  double v325 = dadd(v324, v20);
  double v326 = dadd(v325, v19);
  double v327 = dadd(v326, v18);
  double v328 = dadd(v327, v17);
  double v329 = dmul(v328, v16);
  double v330 = DENDRO_640;
  double v331 = DENDRO_804;
  double v332 = dmul(v331, v330);
  double v333 = DENDRO_817;
  double v334 = DENDRO_644;
  double v335 = DENDRO_60;
  double v336 = DENDRO_517;
  double v337 = DENDRO_640;
  double v338 = DENDRO_65;
  double v339 = dmul(v338, v337);
  double v340 = dadd(v339, v336);
  double v341 = dmul(v340, v335);
  double v342 = dadd(v341, v334);
  double v343 = dmul(v342, v333);
  double v344 = DENDRO_846;
  double v345 = DENDRO_645;
  double v346 = DENDRO_650;
  double v347 = DENDRO_646;
  double v348 = dmul(v347, negone);
  double v349 = DENDRO_647;
  double v350 = dmul(v349, negone);
  double v351 = dadd(v350, v348);
  double v352 = dadd(v351, v346);
  double v353 = dadd(v352, v345);
  double v354 = dmul(v353, v344);
  double v355 = 12.0;
  double v356 = DENDRO_637;
  double v357 = dmul(v356, v355);
  double v358 = dmul(v357, negone);
  double v359 = 6.00000000000000;
  double v360 = DENDRO_638;
  double v361 = DENDRO_152;
  double v362 = DENDRO_50;
  double v363 = dmul(v362, v361);
  double v364 = DENDRO_240;
  double v365 = DENDRO_46;
  double v366 = dmul(v365, v364);
  double v367 = DENDRO_238;
  double v368 = DENDRO_35;
  double v369 = dmul(v368, v367);
  double v370 = dmul(v369, negone);
  double v371 = DENDRO_343;
  double v372 = DENDRO_60;
  double v373 = gt4[pp];
  double v374 = dmul(v373, v372);
  double v375 = dmul(v374, v371);
  double v376 = dmul(v375, negone);
  double v377 = dadd(v376, v370);
  double v378 = dadd(v377, v366);
  double v379 = dadd(v378, v363);
  double v380 = dmul(v379, v360);
  double v381 = dmul(v380, v359);
  double v382 = dmul(v381, negone);
  double v383 = dadd(v382, v358);
  double v384 = dadd(v383, v354);
  double v385 = dadd(v384, v343);
  double v386 = dadd(v385, v332);
  double v387 = dadd(v386, v329);
  double v388 = dmul(v387, v15);
  double v389 = DENDRO_6;
  double v390 = DENDRO_870;
  double v391 = dmul(v390, v389);
  double v392 = beta0[pp];
  double v393 = *(agrad_0_At4+pp );
  double v394 = dmul(v393, v392);
  double v395 = beta1[pp];
  double v396 = *(agrad_1_At4+pp );
  double v397 = dmul(v396, v395);
  double v398 = beta2[pp];
  double v399 = *(agrad_2_At4+pp );
  double v400 = dmul(v399, v398);
  double v401 = At4[pp];
  double v402 = DENDRO_18;
  double v403 = dmul(v402, v401);
  double v404 = dmul(v403, negone);
  double v405 = alpha[pp];
  double v406 = DENDRO_34;
  double v407 = DENDRO_844;
  double v408 = dmul(v407, v406);
  double v409 = DENDRO_841;
  double v410 = DENDRO_862;
  double v411 = dmul(v410, v409);
  double v412 = DENDRO_845;
  double v413 = DENDRO_863;
  double v414 = dmul(v413, v412);
  double v415 = At4[pp];
  double v416 = K[pp];
  double v417 = dmul(v416, v415);
  double v418 = dmul(v417, negone);
  double v419 = dadd(v418, v414);
  double v420 = dadd(v419, v411);
  double v421 = dadd(v420, v408);
  double v422 = dmul(v421, v405);
  double v423 = dmul(v422, negone);
  double v424 = dadd(v423, v404);
  double v425 = dadd(v424, v400);
  double v426 = dadd(v425, v397);
  double v427 = dadd(v426, v394);
  double v428 = dadd(v427, v391);
  double v429 = dadd(v428, v388);
  double v430 = dadd(v429, v14);
  double v431 = dadd(v430, v11);
  double v432 = dadd(v431, v8);
  double v433 = dadd(v432, v5);
  double v434 = dadd(v433, v2);
  At_rhs12[pp] = v434;
//--
  v0 = At0[pp];
  v1 = DENDRO_3;
  v2 = dmul(v1, v0);
  v3 = DENDRO_11;
  v4 = DENDRO_26;
  v5 = dmul(v4, v3);
  v6 = DENDRO_25;
  v7 = DENDRO_9;
  v8 = dmul(v7, v6);
  v9 = DENDRO_56;
  v10 = 12.0;
  v11 = DENDRO_81;
  v12 = dmul(v11, v10);
  v13 = DENDRO_805;
  v14 = gt0[pp];
  v15 = dmul(v14, v13);
  v16 = DENDRO_83;
  v17 = DENDRO_89;
  v18 = DENDRO_93;
  v19 = dmul(v18, negone);
  v20 = dadd(v19, v17);
  v21 = dmul(v20, v16);
  v22 = 12.0;
  v23 = DENDRO_57;
  v24 = dmul(v23, v22);
  v25 = dmul(v24, negone);
  v26 = DENDRO_108;
  v27 = DENDRO_168;
  v28 = dmul(v27, negone);
  v29 = DENDRO_139;
  v30 = DENDRO_209;
  v31 = DENDRO_210;
  v32 = dadd(v31, v30);
  v33 = dmul(v32, v29);
  v34 = DENDRO_139;
  v35 = DENDRO_215;
  v36 = DENDRO_197;
  v37 = DENDRO_76;
  v38 = dmul(v37, v36);
  v39 = dadd(v38, v35);
  v40 = dmul(v39, v34);
  v41 = DENDRO_149;
  v42 = DENDRO_203;
  v43 = DENDRO_173;
  v44 = DENDRO_70;
  v45 = dmul(v44, v43);
  v46 = dadd(v45, v42);
  v47 = dmul(v46, v41);
  v48 = DENDRO_149;
  v49 = DENDRO_204;
  v50 = DENDRO_179;
  v51 = DENDRO_77;
  v52 = dmul(v51, v50);
  v53 = dadd(v52, v49);
  v54 = dmul(v53, v48);
  v55 = DENDRO_149;
  v56 = DENDRO_213;
  v57 = 0.500000000000000;
  v58 = DENDRO_214;
  v59 = dmul(v58, v57);
  v60 = dadd(v59, v56);
  v61 = dmul(v60, v55);
  v62 = DENDRO_149;
  v63 = DENDRO_216;
  v64 = DENDRO_218;
  v65 = dmul(v64, negone);
  v66 = dadd(v65, v63);
  v67 = dmul(v66, v62);
  v68 = DENDRO_149;
  v69 = DENDRO_222;
  v70 = DENDRO_225;
  v71 = dadd(v70, v69);
  v72 = dmul(v71, v68);
  v73 = DENDRO_161;
  v74 = DENDRO_212;
  v75 = 1.00000000000000;
  v76 = DENDRO_210;
  v77 = dmul(v76, v75);
  v78 = dadd(v77, v74);
  v79 = dmul(v78, v73);
  v80 = DENDRO_161;
  v81 = DENDRO_197;
  v82 = DENDRO_227;
  v83 = dmul(v82, v81);
  v84 = 1.00000000000000;
  v85 = DENDRO_152;
  v86 = DENDRO_93;
  v87 = dmul(v86, v85);
  v88 = dmul(v87, v84);
  v89 = dadd(v88, v83);
  v90 = dmul(v89, v80);
  v91 = DENDRO_229;
  v92 = DENDRO_231;
  v93 = DENDRO_232;
  v94 = DENDRO_233;
  v95 = dmul(v94, v93);
  v96 = DENDRO_234;
  v97 = DENDRO_80;
  v98 = dmul(v97, v96);
  v99 = DENDRO_235;
  v100 = DENDRO_236;
  v101 = dmul(v100, v99);
  v102 = dadd(v101, v98);
  v103 = dadd(v102, v95);
  v104 = dadd(v103, v92);
  v105 = dmul(v104, v91);
  v106 = DENDRO_270;
  v107 = DENDRO_69;
  v108 = dmul(v107, v106);
  v109 = DENDRO_288;
  v110 = DENDRO_76;
  v111 = dmul(v110, v109);
  v112 = DENDRO_308;
  v113 = DENDRO_72;
  v114 = dmul(v113, v112);
  v115 = DENDRO_336;
  v116 = gt0[pp];
  v117 = dmul(v116, v115);
  v118 = DENDRO_154;
  v119 = DENDRO_206;
  v120 = DENDRO_205;
  v121 = DENDRO_76;
  v122 = dmul(v121, v120);
  v123 = dadd(v122, v119);
  v124 = dmul(v123, v118);
  v125 = dmul(v124, negone);
  v126 = DENDRO_154;
  v127 = DENDRO_207;
  v128 = DENDRO_208;
  v129 = dadd(v128, v127);
  v130 = dmul(v129, v126);
  v131 = dmul(v130, negone);
  v132 = DENDRO_155;
  v133 = DENDRO_211;
  v134 = 1.00000000000000;
  v135 = DENDRO_208;
  v136 = dmul(v135, v134);
  v137 = dadd(v136, v133);
  v138 = dmul(v137, v132);
  v139 = dmul(v138, negone);
  v140 = DENDRO_155;
  v141 = DENDRO_205;
  v142 = DENDRO_227;
  v143 = dmul(v142, v141);
  v144 = DENDRO_228;
  v145 = DENDRO_93;
  v146 = dmul(v145, v144);
  v147 = dmul(v146, negone);
  v148 = dadd(v147, v143);
  v149 = dmul(v148, v140);
  v150 = dmul(v149, negone);
  v151 = DENDRO_169;
  v152 = DENDRO_174;
  v153 = dmul(v152, v151);
  v154 = dmul(v153, negone);
  v155 = DENDRO_175;
  v156 = DENDRO_180;
  v157 = dmul(v156, v155);
  v158 = dmul(v157, negone);
  v159 = DENDRO_182;
  v160 = DENDRO_150;
  v161 = DENDRO_166;
  v162 = dadd(v161, v160);
  v163 = dmul(v162, v159);
  v164 = dmul(v163, negone);
  v165 = DENDRO_193;
  v166 = DENDRO_185;
  v167 = DENDRO_188;
  v168 = dadd(v167, v166);
  v169 = dmul(v168, v165);
  v170 = dmul(v169, negone);
  v171 = DENDRO_197;
  v172 = DENDRO_200;
  v173 = dmul(v172, v171);
  v174 = dmul(v173, negone);
  v175 = DENDRO_201;
  v176 = DENDRO_93;
  v177 = dmul(v176, v175);
  v178 = dmul(v177, negone);
  v179 = DENDRO_183;
  v180 = DENDRO_184;
  v181 = DENDRO_50;
  v182 = dmul(v181, v180);
  v183 = dmul(v182, v179);
  v184 = dmul(v183, negone);
  v185 = dadd(v184, v178);
  v186 = dadd(v185, v174);
  v187 = dadd(v186, v170);
  v188 = dadd(v187, v164);
  v189 = dadd(v188, v158);
  v190 = dadd(v189, v154);
  v191 = dadd(v190, v150);
  v192 = dadd(v191, v139);
  v193 = dadd(v192, v131);
  v194 = dadd(v193, v125);
  v195 = dadd(v194, v117);
  v196 = dadd(v195, v114);
  v197 = dadd(v196, v111);
  v198 = dadd(v197, v108);
  v199 = dadd(v198, v105);
  v200 = dadd(v199, v90);
  v201 = dadd(v200, v79);
  v202 = dadd(v201, v72);
  v203 = dadd(v202, v67);
  v204 = dadd(v203, v61);
  v205 = dadd(v204, v54);
  v206 = dadd(v205, v47);
  v207 = dadd(v206, v40);
  v208 = dadd(v207, v33);
  v209 = dadd(v208, v28);
  v210 = dmul(v209, v26);
  v211 = dmul(v210, negone);
  v212 = 12.0;
  v213 = DENDRO_94;
  v214 = DENDRO_100;
  v215 = DENDRO_107;
  v216 = DENDRO_96;
  v217 = dmul(v216, negone);
  v218 = DENDRO_98;
  v219 = dmul(v218, negone);
  v220 = dadd(v219, v217);
  v221 = dadd(v220, v215);
  v222 = dadd(v221, v214);
  v223 = dmul(v222, v213);
  v224 = dmul(v223, v212);
  v225 = dadd(v224, v211);
  v226 = dadd(v225, v25);
  v227 = dadd(v226, v21);
  v228 = dadd(v227, v15);
  v229 = dadd(v228, v12);
  v230 = dmul(v229, v9);
  v231 = alpha[pp];
  v232 = At0[pp];
  v233 = K[pp];
  v234 = dmul(v233, v232);
  v235 = DENDRO_51;
  v236 = DENDRO_53;
  v237 = DENDRO_55;
  v238 = DENDRO_52;
  v239 = dmul(v238, negone);
  v240 = dadd(v239, v237);
  v241 = dadd(v240, v236);
  v242 = dmul(v241, v235);
  v243 = DENDRO_34;
  v244 = DENDRO_36;
  v245 = DENDRO_38;
  v246 = DENDRO_41;
  v247 = dmul(v246, negone);
  v248 = dadd(v247, v245);
  v249 = dadd(v248, v244);
  v250 = dmul(v249, v243);
  v251 = dmul(v250, negone);
  v252 = DENDRO_42;
  v253 = DENDRO_43;
  v254 = DENDRO_48;
  v255 = At0[pp];
  v256 = DENDRO_50;
  v257 = dmul(v256, v255);
  v258 = dmul(v257, negone);
  v259 = dadd(v258, v254);
  v260 = dadd(v259, v253);
  v261 = dmul(v260, v252);
  v262 = dmul(v261, negone);
  v263 = dadd(v262, v251);
  v264 = dadd(v263, v242);
  v265 = dadd(v264, v234);
  v266 = dmul(v265, v231);
  v267 = beta0[pp];
  v268 = *(agrad_0_At0+pp );
  v269 = dmul(v268, v267);
  v270 = beta1[pp];
  v271 = *(agrad_1_At0+pp );
  v272 = dmul(v271, v270);
  v273 = beta2[pp];
  v274 = *(agrad_2_At0+pp );
  v275 = dmul(v274, v273);
  v276 = DENDRO_24;
  v277 = DENDRO_4;
  v278 = dmul(v277, v276);
  v279 = dmul(v278, negone);
  v280 = DENDRO_24;
  v281 = DENDRO_6;
  v282 = dmul(v281, v280);
  v283 = dmul(v282, negone);
  v284 = dadd(v283, v279);
  v285 = dadd(v284, v275);
  v286 = dadd(v285, v272);
  v287 = dadd(v286, v269);
  v288 = dadd(v287, v266);
  v289 = dadd(v288, v230);
  v290 = dadd(v289, v8);
  v291 = dadd(v290, v5);
  v292 = dadd(v291, v2);
  At_rhs00[pp] = v292;
//--
  v0 = At5[pp];
  v1 = DENDRO_22;
  v2 = dmul(v1, v0);
  v3 = DENDRO_15;
  v4 = DENDRO_26;
  v5 = dmul(v4, v3);
  v6 = DENDRO_16;
  v7 = DENDRO_861;
  v8 = dmul(v7, v6);
  v9 = beta0[pp];
  v10 = *(agrad_0_At5+pp );
  v11 = dmul(v10, v9);
  v12 = beta1[pp];
  v13 = *(agrad_1_At5+pp );
  v14 = dmul(v13, v12);
  v15 = beta2[pp];
  v16 = *(agrad_2_At5+pp );
  v17 = dmul(v16, v15);
  v18 = At5[pp];
  v19 = DENDRO_18;
  v20 = dmul(v19, v18);
  v21 = dmul(v20, negone);
  v22 = At5[pp];
  v23 = DENDRO_5;
  v24 = dmul(v23, v22);
  v25 = dmul(v24, negone);
  v26 = DENDRO_56;
  v27 = 12.0;
  v28 = DENDRO_337;
  v29 = dmul(v28, v27);
  v30 = DENDRO_108;
  v31 = DENDRO_372;
  v32 = dmul(v31, negone);
  v33 = DENDRO_134;
  v34 = DENDRO_308;
  v35 = dmul(v34, v33);
  v36 = DENDRO_149;
  v37 = DENDRO_398;
  v38 = 1.00000000000000;
  v39 = DENDRO_400;
  v40 = dmul(v39, v38);
  v41 = dmul(v40, negone);
  v42 = dadd(v41, v37);
  v43 = dmul(v42, v36);
  v44 = DENDRO_149;
  v45 = 0.250000000000000;
  v46 = DENDRO_875;
  v47 = dmul(v46, v45);
  v48 = 1.00000000000000;
  v49 = DENDRO_876;
  v50 = dmul(v49, v48);
  v51 = dadd(v50, v47);
  v52 = dmul(v51, v44);
  v53 = DENDRO_155;
  v54 = DENDRO_856;
  v55 = DENDRO_402;
  v56 = DENDRO_857;
  v57 = dmul(v56, v55);
  v58 = dmul(v57, negone);
  v59 = dadd(v58, v54);
  v60 = dmul(v59, v53);
  v61 = DENDRO_161;
  v62 = DENDRO_850;
  v63 = DENDRO_219;
  v64 = DENDRO_819;
  v65 = dmul(v64, v63);
  v66 = dadd(v65, v62);
  v67 = dmul(v66, v61);
  v68 = DENDRO_161;
  v69 = DENDRO_871;
  v70 = DENDRO_205;
  v71 = DENDRO_250;
  v72 = dmul(v71, v70);
  v73 = dadd(v72, v69);
  v74 = dmul(v73, v68);
  v75 = DENDRO_161;
  v76 = 0.250000000000000;
  v77 = DENDRO_825;
  v78 = dmul(v77, v76);
  v79 = DENDRO_78;
  v80 = DENDRO_868;
  v81 = dmul(v80, v79);
  v82 = dadd(v81, v78);
  v83 = dmul(v82, v75);
  v84 = DENDRO_161;
  v85 = 0.500000000000000;
  v86 = DENDRO_859;
  v87 = dmul(v86, v85);
  v88 = DENDRO_202;
  v89 = DENDRO_821;
  v90 = dmul(v89, v88);
  v91 = dadd(v90, v87);
  v92 = dmul(v91, v84);
  v93 = DENDRO_229;
  v94 = DENDRO_403;
  v95 = DENDRO_232;
  v96 = DENDRO_281;
  v97 = dmul(v96, v95);
  v98 = DENDRO_234;
  v99 = DENDRO_255;
  v100 = dmul(v99, v98);
  v101 = DENDRO_235;
  v102 = DENDRO_301;
  v103 = dmul(v102, v101);
  v104 = dadd(v103, v100);
  v105 = dadd(v104, v97);
  v106 = dadd(v105, v94);
  v107 = dmul(v106, v93);
  v108 = DENDRO_240;
  v109 = DENDRO_270;
  v110 = dmul(v109, v108);
  v111 = DENDRO_247;
  v112 = DENDRO_288;
  v113 = dmul(v112, v111);
  v114 = DENDRO_336;
  v115 = gt5[pp];
  v116 = dmul(v115, v114);
  v117 = DENDRO_362;
  v118 = DENDRO_875;
  v119 = DENDRO_876;
  v120 = dadd(v119, v118);
  v121 = dmul(v120, v117);
  v122 = DENDRO_362;
  v123 = DENDRO_134;
  v124 = DENDRO_821;
  v125 = dmul(v124, v123);
  v126 = DENDRO_152;
  v127 = DENDRO_857;
  v128 = dmul(v127, v126);
  v129 = dadd(v128, v125);
  v130 = dmul(v129, v122);
  v131 = DENDRO_154;
  v132 = DENDRO_877;
  v133 = DENDRO_878;
  v134 = dadd(v133, v132);
  v135 = dmul(v134, v131);
  v136 = dmul(v135, negone);
  v137 = DENDRO_154;
  v138 = DENDRO_134;
  v139 = DENDRO_173;
  v140 = dmul(v139, v138);
  v141 = DENDRO_76;
  v142 = DENDRO_857;
  v143 = dmul(v142, v141);
  v144 = dadd(v143, v140);
  v145 = dmul(v144, v137);
  v146 = dmul(v145, negone);
  v147 = DENDRO_155;
  v148 = DENDRO_394;
  v149 = 1.00000000000000;
  v150 = DENDRO_396;
  v151 = dmul(v150, v149);
  v152 = dmul(v151, negone);
  v153 = dadd(v152, v148);
  v154 = dmul(v153, v147);
  v155 = dmul(v154, negone);
  v156 = DENDRO_155;
  v157 = 0.250000000000000;
  v158 = DENDRO_877;
  v159 = dmul(v158, v157);
  v160 = 1.00000000000000;
  v161 = DENDRO_878;
  v162 = dmul(v161, v160);
  v163 = dadd(v162, v159);
  v164 = dmul(v163, v156);
  v165 = dmul(v164, negone);
  v166 = DENDRO_173;
  v167 = DENDRO_377;
  v168 = dmul(v167, v166);
  v169 = dmul(v168, negone);
  v170 = DENDRO_378;
  v171 = DENDRO_821;
  v172 = dmul(v171, v170);
  v173 = dmul(v172, negone);
  v174 = DENDRO_379;
  v175 = DENDRO_188;
  v176 = DENDRO_408;
  v177 = dadd(v176, v175);
  v178 = dmul(v177, v174);
  v179 = dmul(v178, negone);
  v180 = DENDRO_380;
  v181 = DENDRO_253;
  v182 = DENDRO_405;
  v183 = dadd(v182, v181);
  v184 = dmul(v183, v180);
  v185 = dmul(v184, negone);
  v186 = DENDRO_134;
  v187 = DENDRO_205;
  v188 = DENDRO_374;
  v189 = dmul(v188, v187);
  v190 = dmul(v189, v186);
  v191 = dmul(v190, negone);
  v192 = DENDRO_175;
  v193 = DENDRO_240;
  v194 = DENDRO_819;
  v195 = dmul(v194, v193);
  v196 = dmul(v195, v192);
  v197 = dmul(v196, negone);
  v198 = DENDRO_376;
  v199 = DENDRO_54;
  v200 = DENDRO_855;
  v201 = dmul(v200, v199);
  v202 = dmul(v201, v198);
  v203 = dmul(v202, negone);
  v204 = dadd(v203, v197);
  v205 = dadd(v204, v191);
  v206 = dadd(v205, v185);
  v207 = dadd(v206, v179);
  v208 = dadd(v207, v173);
  v209 = dadd(v208, v169);
  v210 = dadd(v209, v165);
  v211 = dadd(v210, v155);
  v212 = dadd(v211, v146);
  v213 = dadd(v212, v136);
  v214 = dadd(v213, v130);
  v215 = dadd(v214, v121);
  v216 = dadd(v215, v116);
  v217 = dadd(v216, v113);
  v218 = dadd(v217, v110);
  v219 = dadd(v218, v107);
  v220 = dadd(v219, v92);
  v221 = dadd(v220, v83);
  v222 = dadd(v221, v74);
  v223 = dadd(v222, v67);
  v224 = dadd(v223, v60);
  v225 = dadd(v224, v52);
  v226 = dadd(v225, v43);
  v227 = dadd(v226, v35);
  v228 = dadd(v227, v32);
  v229 = dmul(v228, v30);
  v230 = 12.0;
  v231 = DENDRO_339;
  v232 = dmul(v231, v230);
  v233 = dmul(v232, negone);
  v234 = DENDRO_805;
  v235 = gt5[pp];
  v236 = dmul(v235, v234);
  v237 = dmul(v236, negone);
  v238 = DENDRO_864;
  v239 = DENDRO_344;
  v240 = DENDRO_857;
  v241 = dmul(v240, negone);
  v242 = dadd(v241, v239);
  v243 = dmul(v242, v238);
  v244 = dmul(v243, negone);
  v245 = 12.0;
  v246 = DENDRO_82;
  v247 = DENDRO_345;
  v248 = DENDRO_347;
  v249 = DENDRO_346;
  v250 = dmul(v249, negone);
  v251 = DENDRO_352;
  v252 = dmul(v251, negone);
  v253 = dadd(v252, v250);
  v254 = dadd(v253, v248);
  v255 = dadd(v254, v247);
  v256 = dmul(v255, v246);
  v257 = dmul(v256, v245);
  v258 = dadd(v257, v244);
  v259 = dadd(v258, v237);
  v260 = dadd(v259, v233);
  v261 = dadd(v260, v229);
  v262 = dadd(v261, v29);
  v263 = dmul(v262, v26);
  v264 = dmul(v263, negone);
  v265 = alpha[pp];
  v266 = DENDRO_51;
  v267 = DENDRO_844;
  v268 = dmul(v267, v266);
  v269 = DENDRO_841;
  v270 = DENDRO_863;
  v271 = dmul(v270, v269);
  v272 = At5[pp];
  v273 = K[pp];
  v274 = dmul(v273, v272);
  v275 = dmul(v274, negone);
  v276 = At5[pp];
  v277 = DENDRO_328;
  v278 = DENDRO_845;
  v279 = dmul(v278, v277);
  v280 = dmul(v279, v276);
  v281 = dadd(v280, v275);
  v282 = dadd(v281, v271);
  v283 = dadd(v282, v268);
  v284 = dmul(v283, v265);
  v285 = dmul(v284, negone);
  v286 = dadd(v285, v264);
  v287 = dadd(v286, v25);
  v288 = dadd(v287, v21);
  v289 = dadd(v288, v17);
  v290 = dadd(v289, v14);
  v291 = dadd(v290, v11);
  v292 = dadd(v291, v8);
  v293 = dadd(v292, v5);
  v294 = dadd(v293, v2);
  At_rhs22[pp] = v294;
//--
  v0 = At0[pp];
  v1 = DENDRO_15;
  v2 = dmul(v1, v0);
  v3 = At1[pp];
  v4 = DENDRO_16;
  v5 = dmul(v4, v3);
  v6 = At4[pp];
  v7 = DENDRO_9;
  v8 = dmul(v7, v6);
  v9 = At5[pp];
  v10 = DENDRO_11;
  v11 = dmul(v10, v9);
  v12 = DENDRO_2;
  v13 = DENDRO_838;
  v14 = dmul(v13, v12);
  v15 = DENDRO_56;
  v16 = 6.00000000000000;
  v17 = DENDRO_504;
  v18 = dmul(v17, v16);
  v19 = DENDRO_108;
  v20 = DENDRO_533;
  v21 = DENDRO_535;
  v22 = DENDRO_537;
  v23 = DENDRO_538;
  v24 = DENDRO_540;
  v25 = DENDRO_542;
  v26 = DENDRO_544;
  v27 = DENDRO_545;
  v28 = DENDRO_546;
  v29 = DENDRO_547;
  v30 = DENDRO_548;
  v31 = DENDRO_549;
  v32 = DENDRO_550;
  v33 = DENDRO_614;
  v34 = DENDRO_621;
  v35 = DENDRO_626;
  v36 = DENDRO_628;
  v37 = DENDRO_629;
  v38 = DENDRO_634;
  v39 = DENDRO_635;
  v40 = DENDRO_636;
  v41 = DENDRO_125;
  v42 = DENDRO_587;
  v43 = DENDRO_610;
  v44 = dmul(v43, negone);
  v45 = dadd(v44, v42);
  v46 = dmul(v45, v41);
  v47 = DENDRO_125;
  v48 = DENDRO_856;
  v49 = dmul(v48, negone);
  v50 = DENDRO_173;
  v51 = DENDRO_185;
  v52 = dmul(v51, v50);
  v53 = DENDRO_77;
  v54 = DENDRO_857;
  v55 = dmul(v54, v53);
  v56 = dadd(v55, v52);
  v57 = dadd(v56, v49);
  v58 = dmul(v57, v47);
  v59 = DENDRO_138;
  v60 = 1.00000000000000;
  v61 = DENDRO_206;
  v62 = dmul(v61, v60);
  v63 = DENDRO_202;
  v64 = DENDRO_205;
  v65 = dmul(v64, v63);
  v66 = dadd(v65, v62);
  v67 = dmul(v66, v59);
  v68 = DENDRO_138;
  v69 = DENDRO_211;
  v70 = DENDRO_184;
  v71 = DENDRO_77;
  v72 = dmul(v71, v70);
  v73 = DENDRO_184;
  v74 = DENDRO_78;
  v75 = dmul(v74, v73);
  v76 = dadd(v75, v72);
  v77 = dadd(v76, v69);
  v78 = dmul(v77, v68);
  v79 = DENDRO_149;
  v80 = DENDRO_849;
  v81 = DENDRO_850;
  v82 = dmul(v81, negone);
  v83 = DENDRO_397;
  v84 = DENDRO_819;
  v85 = dmul(v84, v83);
  v86 = dadd(v85, v82);
  v87 = dadd(v86, v80);
  v88 = dmul(v87, v79);
  v89 = DENDRO_154;
  v90 = DENDRO_174;
  v91 = DENDRO_134;
  v92 = DENDRO_184;
  v93 = dmul(v92, v91);
  v94 = dadd(v93, v90);
  v95 = dmul(v94, v89);
  v96 = DENDRO_155;
  v97 = DENDRO_8;
  v98 = DENDRO_77;
  v99 = DENDRO_855;
  v100 = dmul(v99, v98);
  v101 = dadd(v100, v97);
  v102 = dmul(v101, v96);
  v103 = DENDRO_155;
  v104 = DENDRO_8;
  v105 = DENDRO_205;
  v106 = DENDRO_397;
  v107 = dmul(v106, v105);
  v108 = dmul(v107, negone);
  v109 = dadd(v108, v104);
  v110 = dmul(v109, v103);
  v111 = DENDRO_155;
  v112 = DENDRO_173;
  v113 = DENDRO_227;
  v114 = dmul(v113, v112);
  v115 = DENDRO_73;
  v116 = DENDRO_857;
  v117 = dmul(v116, v115);
  v118 = DENDRO_184;
  v119 = DENDRO_220;
  v120 = dmul(v119, v118);
  v121 = dmul(v120, negone);
  v122 = dadd(v121, v117);
  v123 = dadd(v122, v114);
  v124 = dmul(v123, v111);
  v125 = DENDRO_357;
  v126 = 0.250000000000000;
  v127 = DENDRO_820;
  v128 = dmul(v127, v126);
  v129 = DENDRO_197;
  v130 = DENDRO_250;
  v131 = dmul(v130, v129);
  v132 = dadd(v131, v128);
  v133 = dmul(v132, v125);
  v134 = DENDRO_129;
  v135 = DENDRO_270;
  v136 = dmul(v135, v134);
  v137 = dmul(v136, negone);
  v138 = DENDRO_149;
  v139 = DENDRO_1;
  v140 = 0.250000000000000;
  v141 = DENDRO_859;
  v142 = dmul(v141, v140);
  v143 = DENDRO_70;
  v144 = DENDRO_857;
  v145 = dmul(v144, v143);
  v146 = dadd(v145, v142);
  v147 = dadd(v146, v139);
  v148 = dmul(v147, v138);
  v149 = dmul(v148, negone);
  v150 = DENDRO_149;
  v151 = DENDRO_606;
  v152 = DENDRO_632;
  v153 = DENDRO_605;
  v154 = dmul(v153, negone);
  v155 = dadd(v154, v152);
  v156 = dadd(v155, v151);
  v157 = dmul(v156, v150);
  v158 = dmul(v157, negone);
  v159 = DENDRO_149;
  v160 = DENDRO_848;
  v161 = DENDRO_850;
  v162 = DENDRO_578;
  v163 = DENDRO_855;
  v164 = dmul(v163, v162);
  v165 = dadd(v164, v161);
  v166 = dadd(v165, v160);
  v167 = dmul(v166, v159);
  v168 = dmul(v167, negone);
  v169 = DENDRO_161;
  v170 = DENDRO_829;
  v171 = DENDRO_851;
  v172 = dadd(v171, v170);
  v173 = dmul(v172, v169);
  v174 = dmul(v173, negone);
  v175 = DENDRO_161;
  v176 = DENDRO_832;
  v177 = DENDRO_184;
  v178 = DENDRO_584;
  v179 = dmul(v178, v177);
  v180 = dadd(v179, v176);
  v181 = dmul(v180, v175);
  v182 = dmul(v181, negone);
  v183 = DENDRO_161;
  v184 = DENDRO_834;
  v185 = DENDRO_860;
  v186 = dadd(v185, v184);
  v187 = dmul(v186, v183);
  v188 = dmul(v187, negone);
  v189 = DENDRO_161;
  v190 = DENDRO_851;
  v191 = 0.250000000000000;
  v192 = DENDRO_214;
  v193 = dmul(v192, v191);
  v194 = dadd(v193, v190);
  v195 = dmul(v194, v189);
  v196 = dmul(v195, negone);
  v197 = DENDRO_186;
  v198 = DENDRO_288;
  v199 = dmul(v198, v197);
  v200 = dmul(v199, negone);
  v201 = DENDRO_229;
  v202 = DENDRO_521;
  v203 = DENDRO_136;
  v204 = DENDRO_525;
  v205 = dmul(v204, v203);
  v206 = DENDRO_192;
  v207 = DENDRO_522;
  v208 = dmul(v207, v206);
  v209 = DENDRO_289;
  v210 = DENDRO_527;
  v211 = dmul(v210, v209);
  v212 = dadd(v211, v208);
  v213 = dadd(v212, v205);
  v214 = dadd(v213, v202);
  v215 = dmul(v214, v201);
  v216 = dmul(v215, negone);
  v217 = DENDRO_308;
  v218 = DENDRO_74;
  v219 = dmul(v218, v217);
  v220 = dmul(v219, negone);
  v221 = DENDRO_510;
  v222 = DENDRO_837;
  v223 = dmul(v222, v221);
  v224 = dmul(v223, negone);
  v225 = DENDRO_818;
  v226 = DENDRO_825;
  v227 = DENDRO_847;
  v228 = dadd(v227, v226);
  v229 = dmul(v228, v225);
  v230 = dmul(v229, negone);
  v231 = DENDRO_124;
  v232 = DENDRO_40;
  v233 = DENDRO_823;
  v234 = DENDRO_852;
  v235 = dadd(v234, v233);
  v236 = dmul(v235, v232);
  v237 = dmul(v236, v231);
  v238 = dadd(v237, v230);
  v239 = dadd(v238, v224);
  v240 = dadd(v239, v220);
  v241 = dadd(v240, v216);
  v242 = dadd(v241, v200);
  v243 = dadd(v242, v196);
  v244 = dadd(v243, v188);
  v245 = dadd(v244, v182);
  v246 = dadd(v245, v174);
  v247 = dadd(v246, v168);
  v248 = dadd(v247, v158);
  v249 = dadd(v248, v149);
  v250 = dadd(v249, v137);
  v251 = dadd(v250, v133);
  v252 = dadd(v251, v124);
  v253 = dadd(v252, v110);
  v254 = dadd(v253, v102);
  v255 = dadd(v254, v95);
  v256 = dadd(v255, v88);
  v257 = dadd(v256, v78);
  v258 = dadd(v257, v67);
  v259 = dadd(v258, v58);
  v260 = dadd(v259, v46);
  v261 = dadd(v260, v40);
  v262 = dadd(v261, v39);
  v263 = dadd(v262, v38);
  v264 = dadd(v263, v37);
  v265 = dadd(v264, v36);
  v266 = dadd(v265, v35);
  v267 = dadd(v266, v34);
  v268 = dadd(v267, v33);
  v269 = dadd(v268, v32);
  v270 = dadd(v269, v31);
  v271 = dadd(v270, v30);
  v272 = dadd(v271, v29);
  v273 = dadd(v272, v28);
  v274 = dadd(v273, v27);
  v275 = dadd(v274, v26);
  v276 = dadd(v275, v25);
  v277 = dadd(v276, v24);
  v278 = dadd(v277, v23);
  v279 = dadd(v278, v22);
  v280 = dadd(v279, v21);
  v281 = dadd(v280, v20);
  v282 = dmul(v281, v19);
  v283 = DENDRO_510;
  v284 = DENDRO_804;
  v285 = dmul(v284, v283);
  v286 = 12.0;
  v287 = DENDRO_503;
  v288 = dmul(v287, v286);
  v289 = dmul(v288, negone);
  v290 = DENDRO_816;
  v291 = DENDRO_514;
  v292 = DENDRO_515;
  v293 = DENDRO_516;
  v294 = dmul(v293, negone);
  v295 = DENDRO_519;
  v296 = dmul(v295, negone);
  v297 = dadd(v296, v294);
  v298 = dadd(v297, v292);
  v299 = dadd(v298, v291);
  v300 = dmul(v299, v290);
  v301 = dmul(v300, negone);
  v302 = DENDRO_846;
  v303 = DENDRO_506;
  v304 = DENDRO_507;
  v305 = DENDRO_508;
  v306 = dmul(v305, negone);
  v307 = DENDRO_512;
  v308 = dmul(v307, negone);
  v309 = dadd(v308, v306);
  v310 = dadd(v309, v304);
  v311 = dadd(v310, v303);
  v312 = dmul(v311, v302);
  v313 = dmul(v312, negone);
  v314 = dadd(v313, v301);
  v315 = dadd(v314, v289);
  v316 = dadd(v315, v285);
  v317 = dadd(v316, v282);
  v318 = dadd(v317, v18);
  v319 = dmul(v318, v15);
  v320 = DENDRO_6;
  v321 = DENDRO_838;
  v322 = dmul(v321, v320);
  v323 = beta0[pp];
  v324 = *(agrad_0_At2+pp );
  v325 = dmul(v324, v323);
  v326 = beta1[pp];
  v327 = *(agrad_1_At2+pp );
  v328 = dmul(v327, v326);
  v329 = beta2[pp];
  v330 = *(agrad_2_At2+pp );
  v331 = dmul(v330, v329);
  v332 = At2[pp];
  v333 = DENDRO_5;
  v334 = dmul(v333, v332);
  v335 = dmul(v334, negone);
  v336 = alpha[pp];
  v337 = DENDRO_34;
  v338 = DENDRO_841;
  v339 = dmul(v338, v337);
  v340 = DENDRO_42;
  v341 = DENDRO_844;
  v342 = dmul(v341, v340);
  v343 = DENDRO_51;
  v344 = DENDRO_845;
  v345 = dmul(v344, v343);
  v346 = At2[pp];
  v347 = K[pp];
  v348 = dmul(v347, v346);
  v349 = dmul(v348, negone);
  v350 = dadd(v349, v345);
  v351 = dadd(v350, v342);
  v352 = dadd(v351, v339);
  v353 = dmul(v352, v336);
  v354 = dmul(v353, negone);
  v355 = dadd(v354, v335);
  v356 = dadd(v355, v331);
  v357 = dadd(v356, v328);
  v358 = dadd(v357, v325);
  v359 = dadd(v358, v322);
  v360 = dadd(v359, v319);
  v361 = dadd(v360, v14);
  v362 = dadd(v361, v11);
  v363 = dadd(v362, v8);
  v364 = dadd(v363, v5);
  v365 = dadd(v364, v2);
  At_rhs02[pp] = v365;
//--
  v0 = DENDRO_10;
  v1 = DENDRO_266;
  v2 = DENDRO_946;
  v3 = dmul(v2, v1);
  v4 = DENDRO_914;
  v5 = DENDRO_973;
  v6 = dmul(v5, v4);
  v7 = DENDRO_924;
  v8 = DENDRO_974;
  v9 = dmul(v8, v7);
  v10 = DENDRO_957;
  v11 = DENDRO_960;
  v12 = dmul(v11, v10);
  v13 = lambda[2];
  v14 = beta0[pp];
  v15 = *(agrad_0_B1+pp );
  v16 = dmul(v15, v14);
  v17 = beta1[pp];
  v18 = *(agrad_1_B1+pp );
  v19 = dmul(v18, v17);
  v20 = beta2[pp];
  v21 = *(agrad_2_B1+pp );
  v22 = dmul(v21, v20);
  v23 = dadd(v22, v19);
  v24 = dadd(v23, v16);
  v25 = dmul(v24, v13);
  v26 = B1[pp];
  v27 = eta;
  v28 = dmul(v27, v26);
  v29 = dmul(v28, negone);
  v30 = DENDRO_16;
  v31 = DENDRO_958;
  v32 = dmul(v31, v30);
  v33 = dmul(v32, negone);
  v34 = DENDRO_4;
  v35 = DENDRO_957;
  v36 = dmul(v35, v34);
  v37 = dmul(v36, negone);
  v38 = DENDRO_9;
  v39 = DENDRO_959;
  v40 = dmul(v39, v38);
  v41 = dmul(v40, negone);
  v42 = DENDRO_914;
  v43 = DENDRO_940;
  v44 = dmul(v43, v42);
  v45 = dmul(v44, negone);
  v46 = DENDRO_924;
  v47 = DENDRO_948;
  v48 = dmul(v47, v46);
  v49 = dmul(v48, negone);
  v50 = DENDRO_941;
  v51 = DENDRO_975;
  v52 = DENDRO_924;
  v53 = DENDRO_952;
  v54 = dmul(v53, v52);
  v55 = dadd(v54, v51);
  v56 = dmul(v55, v50);
  v57 = dmul(v56, negone);
  v58 = DENDRO_941;
  v59 = DENDRO_976;
  v60 = DENDRO_914;
  v61 = DENDRO_943;
  v62 = dmul(v61, v60);
  v63 = dadd(v62, v59);
  v64 = dmul(v63, v58);
  v65 = dmul(v64, negone);
  v66 = lambda[3];
  v67 = DENDRO_962;
  v68 = DENDRO_963;
  v69 = DENDRO_964;
  v70 = dadd(v69, v68);
  v71 = dadd(v70, v67);
  v72 = dmul(v71, v66);
  v73 = dmul(v72, negone);
  v74 = dadd(v73, v65);
  v75 = dadd(v74, v57);
  v76 = dadd(v75, v49);
  v77 = dadd(v76, v45);
  v78 = dadd(v77, v41);
  v79 = dadd(v78, v37);
  v80 = dadd(v79, v33);
  v81 = dadd(v80, v29);
  v82 = dadd(v81, v25);
  v83 = dadd(v82, v12);
  v84 = dadd(v83, v9);
  v85 = dadd(v84, v6);
  v86 = dadd(v85, v3);
  v87 = dadd(v86, v0);
  B_rhs1[pp] = v87;
//--
  v0 = At0[pp];
  v1 = DENDRO_12;
  v2 = dmul(v1, v0);
  v3 = At2[pp];
  v4 = DENDRO_13;
  v5 = dmul(v4, v3);
  v6 = At3[pp];
  v7 = DENDRO_9;
  v8 = dmul(v7, v6);
  v9 = At4[pp];
  v10 = DENDRO_11;
  v11 = dmul(v10, v9);
  v12 = DENDRO_2;
  v13 = DENDRO_806;
  v14 = dmul(v13, v12);
  v15 = DENDRO_4;
  v16 = DENDRO_806;
  v17 = dmul(v16, v15);
  v18 = DENDRO_56;
  v19 = DENDRO_108;
  v20 = DENDRO_750;
  v21 = DENDRO_751;
  v22 = DENDRO_752;
  v23 = DENDRO_753;
  v24 = DENDRO_754;
  v25 = DENDRO_755;
  v26 = DENDRO_756;
  v27 = DENDRO_757;
  v28 = DENDRO_758;
  v29 = DENDRO_759;
  v30 = DENDRO_760;
  v31 = DENDRO_761;
  v32 = DENDRO_762;
  v33 = DENDRO_764;
  v34 = DENDRO_765;
  v35 = DENDRO_769;
  v36 = DENDRO_771;
  v37 = DENDRO_775;
  v38 = DENDRO_779;
  v39 = DENDRO_125;
  v40 = DENDRO_593;
  v41 = DENDRO_691;
  v42 = DENDRO_381;
  v43 = dmul(v42, negone);
  v44 = dadd(v43, v41);
  v45 = dadd(v44, v40);
  v46 = dmul(v45, v39);
  v47 = DENDRO_138;
  v48 = DENDRO_212;
  v49 = DENDRO_184;
  v50 = DENDRO_70;
  v51 = dmul(v50, v49);
  v52 = DENDRO_184;
  v53 = DENDRO_71;
  v54 = dmul(v53, v52);
  v55 = dadd(v54, v51);
  v56 = dadd(v55, v48);
  v57 = dmul(v56, v47);
  v58 = DENDRO_138;
  v59 = 0.500000000000000;
  v60 = DENDRO_215;
  v61 = dmul(v60, v59);
  v62 = DENDRO_197;
  v63 = DENDRO_202;
  v64 = dmul(v63, v62);
  v65 = DENDRO_456;
  v66 = DENDRO_93;
  v67 = dmul(v66, v65);
  v68 = dadd(v67, v64);
  v69 = dadd(v68, v61);
  v70 = dmul(v69, v58);
  v71 = DENDRO_149;
  v72 = DENDRO_724;
  v73 = DENDRO_173;
  v74 = DENDRO_465;
  v75 = dmul(v74, v73);
  v76 = DENDRO_179;
  v77 = DENDRO_684;
  v78 = dmul(v77, v76);
  v79 = dmul(v78, negone);
  v80 = dadd(v79, v75);
  v81 = dadd(v80, v72);
  v82 = dmul(v81, v71);
  v83 = DENDRO_155;
  v84 = DENDRO_832;
  v85 = DENDRO_184;
  v86 = DENDRO_456;
  v87 = dmul(v86, v85);
  v88 = dadd(v87, v84);
  v89 = dmul(v88, v83);
  v90 = DENDRO_155;
  v91 = DENDRO_833;
  v92 = DENDRO_834;
  v93 = dadd(v92, v91);
  v94 = dmul(v93, v90);
  v95 = DENDRO_155;
  v96 = DENDRO_225;
  v97 = DENDRO_798;
  v98 = DENDRO_799;
  v99 = dadd(v98, v97);
  v100 = dadd(v99, v96);
  v101 = dmul(v100, v95);
  v102 = DENDRO_155;
  v103 = DENDRO_829;
  v104 = DENDRO_830;
  v105 = DENDRO_831;
  v106 = dadd(v105, v104);
  v107 = dadd(v106, v103);
  v108 = dmul(v107, v102);
  v109 = DENDRO_161;
  v110 = DENDRO_791;
  v111 = DENDRO_167;
  v112 = DENDRO_184;
  v113 = dmul(v112, v111);
  v114 = DENDRO_158;
  v115 = DENDRO_179;
  v116 = dmul(v115, v114);
  v117 = dmul(v116, negone);
  v118 = dadd(v117, v113);
  v119 = dadd(v118, v110);
  v120 = dmul(v119, v109);
  v121 = DENDRO_131;
  v122 = DENDRO_288;
  v123 = dmul(v122, v121);
  v124 = dmul(v123, negone);
  v125 = DENDRO_139;
  v126 = DENDRO_180;
  v127 = DENDRO_140;
  v128 = DENDRO_184;
  v129 = dmul(v128, v127);
  v130 = dadd(v129, v126);
  v131 = dmul(v130, v125);
  v132 = dmul(v131, negone);
  v133 = DENDRO_149;
  v134 = DENDRO_777;
  v135 = dmul(v134, negone);
  v136 = DENDRO_578;
  v137 = DENDRO_827;
  v138 = dmul(v137, v136);
  v139 = dadd(v138, v135);
  v140 = dmul(v139, v133);
  v141 = dmul(v140, negone);
  v142 = DENDRO_149;
  v143 = DENDRO_795;
  v144 = DENDRO_796;
  v145 = DENDRO_694;
  v146 = dmul(v145, negone);
  v147 = dadd(v146, v144);
  v148 = dadd(v147, v143);
  v149 = dmul(v148, v142);
  v150 = dmul(v149, negone);
  v151 = DENDRO_161;
  v152 = DENDRO_774;
  v153 = dmul(v152, negone);
  v154 = DENDRO_70;
  v155 = DENDRO_827;
  v156 = dmul(v155, v154);
  v157 = dadd(v156, v153);
  v158 = dmul(v157, v151);
  v159 = dmul(v158, negone);
  v160 = DENDRO_161;
  v161 = DENDRO_793;
  v162 = dmul(v161, negone);
  v163 = DENDRO_133;
  v164 = DENDRO_836;
  v165 = dmul(v164, v163);
  v166 = DENDRO_263;
  v167 = DENDRO_93;
  v168 = dmul(v167, v166);
  v169 = dadd(v168, v165);
  v170 = dadd(v169, v162);
  v171 = dmul(v170, v160);
  v172 = dmul(v171, negone);
  v173 = DENDRO_161;
  v174 = DENDRO_143;
  v175 = DENDRO_836;
  v176 = dmul(v175, v174);
  v177 = DENDRO_152;
  v178 = DENDRO_836;
  v179 = dmul(v178, v177);
  v180 = DENDRO_262;
  v181 = DENDRO_93;
  v182 = dmul(v181, v180);
  v183 = dadd(v182, v179);
  v184 = dadd(v183, v176);
  v185 = dmul(v184, v173);
  v186 = dmul(v185, negone);
  v187 = DENDRO_164;
  v188 = DENDRO_270;
  v189 = dmul(v188, v187);
  v190 = dmul(v189, negone);
  v191 = DENDRO_229;
  v192 = DENDRO_744;
  v193 = DENDRO_146;
  v194 = DENDRO_525;
  v195 = dmul(v194, v193);
  v196 = DENDRO_221;
  v197 = DENDRO_522;
  v198 = dmul(v197, v196);
  v199 = DENDRO_296;
  v200 = DENDRO_527;
  v201 = dmul(v200, v199);
  v202 = dadd(v201, v198);
  v203 = dadd(v202, v195);
  v204 = dadd(v203, v192);
  v205 = dmul(v204, v191);
  v206 = dmul(v205, negone);
  v207 = DENDRO_308;
  v208 = DENDRO_67;
  v209 = dmul(v208, v207);
  v210 = dmul(v209, negone);
  v211 = DENDRO_357;
  v212 = DENDRO_766;
  v213 = DENDRO_767;
  v214 = dmul(v213, negone);
  v215 = dadd(v214, v212);
  v216 = dmul(v215, v211);
  v217 = dmul(v216, negone);
  v218 = DENDRO_357;
  v219 = DENDRO_782;
  v220 = DENDRO_828;
  v221 = DENDRO_150;
  v222 = DENDRO_179;
  v223 = dmul(v222, v221);
  v224 = dmul(v223, negone);
  v225 = dadd(v224, v220);
  v226 = dadd(v225, v219);
  v227 = dmul(v226, v218);
  v228 = dmul(v227, negone);
  v229 = DENDRO_357;
  v230 = DENDRO_786;
  v231 = DENDRO_835;
  v232 = DENDRO_197;
  v233 = DENDRO_358;
  v234 = dmul(v233, v232);
  v235 = dmul(v234, negone);
  v236 = dadd(v235, v231);
  v237 = dadd(v236, v230);
  v238 = dmul(v237, v229);
  v239 = dmul(v238, negone);
  v240 = DENDRO_733;
  v241 = DENDRO_837;
  v242 = dmul(v241, v240);
  v243 = dmul(v242, negone);
  v244 = DENDRO_818;
  v245 = DENDRO_823;
  v246 = DENDRO_140;
  v247 = DENDRO_173;
  v248 = dmul(v247, v246);
  v249 = dadd(v248, v245);
  v250 = dmul(v249, v244);
  v251 = dmul(v250, negone);
  v252 = DENDRO_818;
  v253 = DENDRO_820;
  v254 = DENDRO_197;
  v255 = DENDRO_240;
  v256 = dmul(v255, v254);
  v257 = DENDRO_205;
  v258 = DENDRO_238;
  v259 = dmul(v258, v257);
  v260 = dadd(v259, v256);
  v261 = dadd(v260, v253);
  v262 = dmul(v261, v252);
  v263 = dmul(v262, negone);
  v264 = DENDRO_124;
  v265 = DENDRO_54;
  v266 = DENDRO_824;
  v267 = DENDRO_825;
  v268 = DENDRO_826;
  v269 = dadd(v268, v267);
  v270 = dadd(v269, v266);
  v271 = dmul(v270, v265);
  v272 = dmul(v271, v264);
  v273 = dadd(v272, v263);
  v274 = dadd(v273, v251);
  v275 = dadd(v274, v243);
  v276 = dadd(v275, v239);
  v277 = dadd(v276, v228);
  v278 = dadd(v277, v217);
  v279 = dadd(v278, v210);
  v280 = dadd(v279, v206);
  v281 = dadd(v280, v190);
  v282 = dadd(v281, v186);
  v283 = dadd(v282, v172);
  v284 = dadd(v283, v159);
  v285 = dadd(v284, v150);
  v286 = dadd(v285, v141);
  v287 = dadd(v286, v132);
  v288 = dadd(v287, v124);
  v289 = dadd(v288, v120);
  v290 = dadd(v289, v108);
  v291 = dadd(v290, v101);
  v292 = dadd(v291, v94);
  v293 = dadd(v292, v89);
  v294 = dadd(v293, v82);
  v295 = dadd(v294, v70);
  v296 = dadd(v295, v57);
  v297 = dadd(v296, v46);
  v298 = dadd(v297, v38);
  v299 = dadd(v298, v37);
  v300 = dadd(v299, v36);
  v301 = dadd(v300, v35);
  v302 = dadd(v301, v34);
  v303 = dadd(v302, v33);
  v304 = dadd(v303, v32);
  v305 = dadd(v304, v31);
  v306 = dadd(v305, v30);
  v307 = dadd(v306, v29);
  v308 = dadd(v307, v28);
  v309 = dadd(v308, v27);
  v310 = dadd(v309, v26);
  v311 = dadd(v310, v25);
  v312 = dadd(v311, v24);
  v313 = dadd(v312, v23);
  v314 = dadd(v313, v22);
  v315 = dadd(v314, v21);
  v316 = dadd(v315, v20);
  v317 = dmul(v316, v19);
  v318 = DENDRO_733;
  v319 = DENDRO_804;
  v320 = dmul(v319, v318);
  v321 = DENDRO_816;
  v322 = DENDRO_738;
  v323 = DENDRO_742;
  v324 = DENDRO_739;
  v325 = dmul(v324, negone);
  v326 = DENDRO_740;
  v327 = dmul(v326, negone);
  v328 = dadd(v327, v325);
  v329 = dadd(v328, v323);
  v330 = dadd(v329, v322);
  v331 = dmul(v330, v321);
  v332 = DENDRO_817;
  v333 = DENDRO_737;
  v334 = DENDRO_60;
  v335 = DENDRO_509;
  v336 = DENDRO_65;
  v337 = DENDRO_733;
  v338 = dmul(v337, v336);
  v339 = dadd(v338, v335);
  v340 = dmul(v339, v334);
  v341 = dadd(v340, v333);
  v342 = dmul(v341, v332);
  v343 = 12.0;
  v344 = DENDRO_731;
  v345 = dmul(v344, v343);
  v346 = dmul(v345, negone);
  v347 = 6.00000000000000;
  v348 = DENDRO_732;
  v349 = DENDRO_140;
  v350 = DENDRO_37;
  v351 = dmul(v350, v349);
  v352 = DENDRO_143;
  v353 = DENDRO_54;
  v354 = dmul(v353, v352);
  v355 = dmul(v354, negone);
  v356 = DENDRO_46;
  v357 = DENDRO_69;
  v358 = dmul(v357, v356);
  v359 = dmul(v358, negone);
  v360 = DENDRO_60;
  v361 = DENDRO_87;
  v362 = gt1[pp];
  v363 = dmul(v362, v361);
  v364 = dmul(v363, v360);
  v365 = dadd(v364, v359);
  v366 = dadd(v365, v355);
  v367 = dadd(v366, v351);
  v368 = dmul(v367, v348);
  v369 = dmul(v368, v347);
  v370 = dadd(v369, v346);
  v371 = dadd(v370, v342);
  v372 = dadd(v371, v331);
  v373 = dadd(v372, v320);
  v374 = dadd(v373, v317);
  v375 = dmul(v374, v18);
  v376 = beta0[pp];
  v377 = *(agrad_0_At1+pp );
  v378 = dmul(v377, v376);
  v379 = beta1[pp];
  v380 = *(agrad_1_At1+pp );
  v381 = dmul(v380, v379);
  v382 = beta2[pp];
  v383 = *(agrad_2_At1+pp );
  v384 = dmul(v383, v382);
  v385 = At1[pp];
  v386 = DENDRO_7;
  v387 = dmul(v386, v385);
  v388 = dmul(v387, negone);
  v389 = alpha[pp];
  v390 = DENDRO_34;
  v391 = DENDRO_808;
  v392 = dmul(v391, v390);
  v393 = DENDRO_42;
  v394 = DENDRO_811;
  v395 = dmul(v394, v393);
  v396 = DENDRO_51;
  v397 = DENDRO_815;
  v398 = dmul(v397, v396);
  v399 = At1[pp];
  v400 = K[pp];
  v401 = dmul(v400, v399);
  v402 = dmul(v401, negone);
  v403 = dadd(v402, v398);
  v404 = dadd(v403, v395);
  v405 = dadd(v404, v392);
  v406 = dmul(v405, v389);
  v407 = dmul(v406, negone);
  v408 = dadd(v407, v388);
  v409 = dadd(v408, v384);
  v410 = dadd(v409, v381);
  v411 = dadd(v410, v378);
  v412 = dadd(v411, v375);
  v413 = dadd(v412, v17);
  v414 = dadd(v413, v14);
  v415 = dadd(v414, v11);
  v416 = dadd(v415, v8);
  v417 = dadd(v416, v5);
  v418 = dadd(v417, v2);
  At_rhs01[pp] = v418;
//--
  v0 = At3[pp];
  v1 = DENDRO_19;
  v2 = dmul(v1, v0);
  v3 = DENDRO_12;
  v4 = DENDRO_25;
  v5 = dmul(v4, v3);
  v6 = DENDRO_13;
  v7 = DENDRO_861;
  v8 = dmul(v7, v6);
  v9 = DENDRO_56;
  v10 = DENDRO_108;
  v11 = DENDRO_437;
  v12 = DENDRO_139;
  v13 = DENDRO_447;
  v14 = DENDRO_866;
  v15 = dmul(v14, negone);
  v16 = dadd(v15, v13);
  v17 = dmul(v16, v12);
  v18 = DENDRO_149;
  v19 = DENDRO_451;
  v20 = 1.00000000000000;
  v21 = DENDRO_865;
  v22 = dmul(v21, v20);
  v23 = dmul(v22, negone);
  v24 = dadd(v23, v19);
  v25 = dmul(v24, v18);
  v26 = DENDRO_155;
  v27 = DENDRO_462;
  v28 = DENDRO_464;
  v29 = dmul(v28, negone);
  v30 = dadd(v29, v27);
  v31 = dmul(v30, v26);
  v32 = DENDRO_155;
  v33 = DENDRO_867;
  v34 = DENDRO_456;
  v35 = DENDRO_819;
  v36 = dmul(v35, v34);
  v37 = dadd(v36, v33);
  v38 = dmul(v37, v32);
  v39 = DENDRO_155;
  v40 = DENDRO_869;
  v41 = DENDRO_71;
  v42 = DENDRO_868;
  v43 = dmul(v42, v41);
  v44 = dadd(v43, v40);
  v45 = dmul(v44, v39);
  v46 = DENDRO_155;
  v47 = DENDRO_156;
  v48 = DENDRO_821;
  v49 = dmul(v48, v47);
  v50 = DENDRO_179;
  v51 = DENDRO_456;
  v52 = dmul(v51, v50);
  v53 = dadd(v52, v49);
  v54 = dmul(v53, v46);
  v55 = DENDRO_161;
  v56 = DENDRO_452;
  v57 = 1.00000000000000;
  v58 = DENDRO_866;
  v59 = dmul(v58, v57);
  v60 = dmul(v59, negone);
  v61 = dadd(v60, v56);
  v62 = dmul(v61, v55);
  v63 = DENDRO_161;
  v64 = DENDRO_473;
  v65 = DENDRO_828;
  v66 = dadd(v65, v64);
  v67 = dmul(v66, v63);
  v68 = DENDRO_161;
  v69 = DENDRO_475;
  v70 = DENDRO_835;
  v71 = dadd(v70, v69);
  v72 = dmul(v71, v68);
  v73 = DENDRO_179;
  v74 = DENDRO_440;
  v75 = dmul(v74, v73);
  v76 = DENDRO_197;
  v77 = DENDRO_443;
  v78 = dmul(v77, v76);
  v79 = DENDRO_362;
  v80 = DENDRO_445;
  v81 = DENDRO_865;
  v82 = dmul(v81, negone);
  v83 = dadd(v82, v80);
  v84 = dmul(v83, v79);
  v85 = DENDRO_439;
  v86 = DENDRO_253;
  v87 = DENDRO_369;
  v88 = dadd(v87, v86);
  v89 = dmul(v88, v85);
  v90 = DENDRO_442;
  v91 = DENDRO_821;
  v92 = dmul(v91, v90);
  v93 = DENDRO_444;
  v94 = DENDRO_166;
  v95 = DENDRO_477;
  v96 = dadd(v95, v94);
  v97 = dmul(v96, v93);
  v98 = DENDRO_139;
  v99 = DENDRO_450;
  v100 = dmul(v99, negone);
  v101 = DENDRO_140;
  v102 = DENDRO_179;
  v103 = dmul(v102, v101);
  v104 = dadd(v103, v100);
  v105 = dmul(v104, v98);
  v106 = dmul(v105, negone);
  v107 = DENDRO_139;
  v108 = DENDRO_459;
  v109 = dmul(v108, negone);
  v110 = DENDRO_197;
  v111 = DENDRO_238;
  v112 = dmul(v111, v110);
  v113 = dadd(v112, v109);
  v114 = dmul(v113, v107);
  v115 = dmul(v114, negone);
  v116 = DENDRO_140;
  v117 = DENDRO_308;
  v118 = dmul(v117, v116);
  v119 = dmul(v118, negone);
  v120 = DENDRO_149;
  v121 = DENDRO_466;
  v122 = 1.00000000000000;
  v123 = DENDRO_468;
  v124 = dmul(v123, v122);
  v125 = dmul(v124, negone);
  v126 = dadd(v125, v121);
  v127 = dmul(v126, v120);
  v128 = dmul(v127, negone);
  v129 = DENDRO_149;
  v130 = DENDRO_471;
  v131 = DENDRO_469;
  v132 = DENDRO_819;
  v133 = dmul(v132, v131);
  v134 = dadd(v133, v130);
  v135 = dmul(v134, v129);
  v136 = dmul(v135, negone);
  v137 = DENDRO_229;
  v138 = DENDRO_476;
  v139 = DENDRO_232;
  v140 = DENDRO_284;
  v141 = dmul(v140, v139);
  v142 = DENDRO_234;
  v143 = DENDRO_266;
  v144 = dmul(v143, v142);
  v145 = DENDRO_235;
  v146 = DENDRO_304;
  v147 = dmul(v146, v145);
  v148 = dadd(v147, v144);
  v149 = dadd(v148, v141);
  v150 = dadd(v149, v138);
  v151 = dmul(v150, v137);
  v152 = dmul(v151, negone);
  v153 = DENDRO_238;
  v154 = DENDRO_288;
  v155 = dmul(v154, v153);
  v156 = dmul(v155, negone);
  v157 = DENDRO_257;
  v158 = DENDRO_270;
  v159 = dmul(v158, v157);
  v160 = dmul(v159, negone);
  v161 = DENDRO_336;
  v162 = gt3[pp];
  v163 = dmul(v162, v161);
  v164 = dmul(v163, negone);
  v165 = DENDRO_362;
  v166 = DENDRO_449;
  v167 = dmul(v166, negone);
  v168 = DENDRO_238;
  v169 = DENDRO_819;
  v170 = dmul(v169, v168);
  v171 = dadd(v170, v167);
  v172 = dmul(v171, v165);
  v173 = dmul(v172, negone);
  v174 = DENDRO_362;
  v175 = DENDRO_458;
  v176 = dmul(v175, negone);
  v177 = DENDRO_140;
  v178 = DENDRO_821;
  v179 = dmul(v178, v177);
  v180 = dadd(v179, v176);
  v181 = dmul(v180, v174);
  v182 = dmul(v181, negone);
  v183 = DENDRO_40;
  v184 = DENDRO_438;
  v185 = DENDRO_827;
  v186 = dmul(v185, v184);
  v187 = dmul(v186, v183);
  v188 = dadd(v187, v182);
  v189 = dadd(v188, v173);
  v190 = dadd(v189, v164);
  v191 = dadd(v190, v160);
  v192 = dadd(v191, v156);
  v193 = dadd(v192, v152);
  v194 = dadd(v193, v136);
  v195 = dadd(v194, v128);
  v196 = dadd(v195, v119);
  v197 = dadd(v196, v115);
  v198 = dadd(v197, v106);
  v199 = dadd(v198, v97);
  v200 = dadd(v199, v92);
  v201 = dadd(v200, v89);
  v202 = dadd(v201, v84);
  v203 = dadd(v202, v78);
  v204 = dadd(v203, v75);
  v205 = dadd(v204, v72);
  v206 = dadd(v205, v67);
  v207 = dadd(v206, v62);
  v208 = dadd(v207, v54);
  v209 = dadd(v208, v45);
  v210 = dadd(v209, v38);
  v211 = dadd(v210, v31);
  v212 = dadd(v211, v25);
  v213 = dadd(v212, v17);
  v214 = dadd(v213, v11);
  v215 = dmul(v214, v10);
  v216 = DENDRO_805;
  v217 = gt3[pp];
  v218 = dmul(v217, v216);
  v219 = DENDRO_83;
  v220 = DENDRO_284;
  v221 = DENDRO_422;
  v222 = DENDRO_85;
  v223 = DENDRO_86;
  v224 = DENDRO_84;
  v225 = dmul(v224, negone);
  v226 = dadd(v225, v223);
  v227 = dadd(v226, v222);
  v228 = dmul(v227, v221);
  v229 = dmul(v228, negone);
  v230 = dadd(v229, v220);
  v231 = dmul(v230, v219);
  v232 = DENDRO_864;
  v233 = DENDRO_304;
  v234 = DENDRO_422;
  v235 = DENDRO_341;
  v236 = DENDRO_342;
  v237 = DENDRO_103;
  v238 = dmul(v237, negone);
  v239 = dadd(v238, v236);
  v240 = dadd(v239, v235);
  v241 = dmul(v240, v234);
  v242 = dmul(v241, negone);
  v243 = dadd(v242, v233);
  v244 = dmul(v243, v232);
  v245 = 12.0;
  v246 = DENDRO_419;
  v247 = dmul(v246, v245);
  v248 = dmul(v247, negone);
  v249 = 12.0;
  v250 = DENDRO_58;
  v251 = DENDRO_429;
  v252 = DENDRO_60;
  v253 = DENDRO_424;
  v254 = DENDRO_425;
  v255 = DENDRO_65;
  v256 = dmul(v255, v254);
  v257 = dadd(v256, v253);
  v258 = dmul(v257, v252);
  v259 = dadd(v258, v251);
  v260 = dmul(v259, v250);
  v261 = dmul(v260, v249);
  v262 = dadd(v261, v248);
  v263 = dadd(v262, v244);
  v264 = dadd(v263, v231);
  v265 = dadd(v264, v218);
  v266 = dadd(v265, v215);
  v267 = dmul(v266, v9);
  v268 = beta0[pp];
  v269 = *(agrad_0_At3+pp );
  v270 = dmul(v269, v268);
  v271 = beta1[pp];
  v272 = *(agrad_1_At3+pp );
  v273 = dmul(v272, v271);
  v274 = beta2[pp];
  v275 = *(agrad_2_At3+pp );
  v276 = dmul(v275, v274);
  v277 = At3[pp];
  v278 = DENDRO_18;
  v279 = dmul(v278, v277);
  v280 = dmul(v279, negone);
  v281 = At3[pp];
  v282 = DENDRO_7;
  v283 = dmul(v282, v281);
  v284 = dmul(v283, negone);
  v285 = alpha[pp];
  v286 = DENDRO_34;
  v287 = DENDRO_811;
  v288 = dmul(v287, v286);
  v289 = DENDRO_808;
  v290 = DENDRO_862;
  v291 = dmul(v290, v289);
  v292 = DENDRO_815;
  v293 = DENDRO_863;
  v294 = dmul(v293, v292);
  v295 = At3[pp];
  v296 = K[pp];
  v297 = dmul(v296, v295);
  v298 = dmul(v297, negone);
  v299 = dadd(v298, v294);
  v300 = dadd(v299, v291);
  v301 = dadd(v300, v288);
  v302 = dmul(v301, v285);
  v303 = dmul(v302, negone);
  v304 = dadd(v303, v284);
  v305 = dadd(v304, v280);
  v306 = dadd(v305, v276);
  v307 = dadd(v306, v273);
  v308 = dadd(v307, v270);
  v309 = dadd(v308, v267);
  v310 = dadd(v309, v8);
  v311 = dadd(v310, v5);
  v312 = dadd(v311, v2);
  At_rhs11[pp] = v312;
//--
  v0 = beta0[pp];
  v1 = *(agrad_0_K+pp );
  v2 = dmul(v1, v0);
  v3 = beta1[pp];
  v4 = *(agrad_1_K+pp );
  v5 = dmul(v4, v3);
  v6 = beta2[pp];
  v7 = *(agrad_2_K+pp );
  v8 = dmul(v7, v6);
  v9 = 0.3333333333333333;
  v10 = alpha[pp];
  v11 = K[pp] * K[pp];
  v12 = At0[pp];
  v13 = DENDRO_893;
  v14 = DENDRO_898;
  v15 = dmul(v14, v13);
  v16 = dmul(v15, v12);
  v17 = At1[pp];
  v18 = DENDRO_903;
  v19 = DENDRO_914;
  v20 = dmul(v19, v18);
  v21 = dmul(v20, v17);
  v22 = At2[pp];
  v23 = DENDRO_903;
  v24 = DENDRO_904;
  v25 = dmul(v24, v23);
  v26 = dmul(v25, v22);
  v27 = At3[pp];
  v28 = DENDRO_893;
  v29 = DENDRO_901;
  v30 = dmul(v29, v28);
  v31 = dmul(v30, v27);
  v32 = At4[pp];
  v33 = DENDRO_903;
  v34 = DENDRO_924;
  v35 = dmul(v34, v33);
  v36 = dmul(v35, v32);
  v37 = At5[pp];
  v38 = DENDRO_893;
  v39 = DENDRO_902;
  v40 = dmul(v39, v38);
  v41 = dmul(v40, v37);
  v42 = dadd(v41, v36);
  v43 = dadd(v42, v31);
  v44 = dadd(v43, v26);
  v45 = dadd(v44, v21);
  v46 = dadd(v45, v16);
  v47 = dadd(v46, v11);
  v48 = dmul(v47, v10);
  v49 = dmul(v48, v9);
  v50 = DENDRO_886;
  v51 = chi[pp];
  v52 = DENDRO_637;
  v53 = dmul(v52, negone);
  v54 = DENDRO_888;
  v55 = DENDRO_645;
  v56 = DENDRO_33;
  v57 = DENDRO_653;
  v58 = dmul(v57, v56);
  v59 = DENDRO_33;
  v60 = DENDRO_654;
  v61 = dmul(v60, v59);
  v62 = DENDRO_60;
  v63 = DENDRO_61;
  v64 = DENDRO_649;
  v65 = dmul(v64, negone);
  v66 = dadd(v65, v63);
  v67 = dmul(v66, v62);
  v68 = dmul(v67, negone);
  v69 = dadd(v68, v61);
  v70 = dadd(v69, v58);
  v71 = dadd(v70, v55);
  v72 = dmul(v71, v54);
  v73 = DENDRO_889;
  v74 = DENDRO_642;
  v75 = DENDRO_643;
  v76 = DENDRO_33;
  v77 = DENDRO_655;
  v78 = dmul(v77, v76);
  v79 = DENDRO_60;
  v80 = DENDRO_62;
  v81 = DENDRO_641;
  v82 = dmul(v81, negone);
  v83 = dadd(v82, v80);
  v84 = dmul(v83, v79);
  v85 = dmul(v84, negone);
  v86 = dadd(v85, v78);
  v87 = dadd(v86, v75);
  v88 = dadd(v87, v74);
  v89 = dmul(v88, v73);
  v90 = DENDRO_887;
  v91 = DENDRO_94;
  v92 = DENDRO_656;
  v93 = DENDRO_882;
  v94 = gt4[pp];
  v95 = dmul(v94, v93);
  v96 = dadd(v95, v92);
  v97 = dmul(v96, v91);
  v98 = dmul(v97, v90);
  v99 = dadd(v98, v89);
  v100 = dadd(v99, v72);
  v101 = dadd(v100, v53);
  v102 = dmul(v101, v51);
  v103 = dmul(v102, v50);
  v104 = DENDRO_892;
  v105 = chi[pp];
  v106 = DENDRO_731;
  v107 = dmul(v106, negone);
  v108 = DENDRO_889;
  v109 = DENDRO_735;
  v110 = DENDRO_736;
  v111 = DENDRO_33;
  v112 = DENDRO_746;
  v113 = dmul(v112, v111);
  v114 = DENDRO_60;
  v115 = DENDRO_63;
  v116 = DENDRO_734;
  v117 = dmul(v116, negone);
  v118 = dadd(v117, v115);
  v119 = dmul(v118, v114);
  v120 = dmul(v119, negone);
  v121 = dadd(v120, v113);
  v122 = dadd(v121, v110);
  v123 = dadd(v122, v109);
  v124 = dmul(v123, v108);
  v125 = DENDRO_891;
  v126 = DENDRO_738;
  v127 = DENDRO_33;
  v128 = DENDRO_747;
  v129 = dmul(v128, v127);
  v130 = DENDRO_33;
  v131 = DENDRO_748;
  v132 = dmul(v131, v130);
  v133 = DENDRO_60;
  v134 = DENDRO_61;
  v135 = DENDRO_741;
  v136 = dmul(v135, negone);
  v137 = dadd(v136, v134);
  v138 = dmul(v137, v133);
  v139 = dmul(v138, negone);
  v140 = dadd(v139, v132);
  v141 = dadd(v140, v129);
  v142 = dadd(v141, v126);
  v143 = dmul(v142, v125);
  v144 = DENDRO_82;
  v145 = DENDRO_887;
  v146 = DENDRO_745;
  v147 = DENDRO_884;
  v148 = gt1[pp];
  v149 = dmul(v148, v147);
  v150 = dadd(v149, v146);
  v151 = dmul(v150, v145);
  v152 = dmul(v151, v144);
  v153 = dadd(v152, v143);
  v154 = dadd(v153, v124);
  v155 = dadd(v154, v107);
  v156 = dmul(v155, v105);
  v157 = dmul(v156, v104);
  v158 = DENDRO_879;
  v159 = chi[pp];
  v160 = DENDRO_337;
  v161 = dmul(v160, negone);
  v162 = DENDRO_638;
  v163 = DENDRO_413;
  v164 = DENDRO_338;
  v165 = DENDRO_882;
  v166 = dmul(v165, v164);
  v167 = dadd(v166, v163);
  v168 = dmul(v167, v162);
  v169 = DENDRO_82;
  v170 = DENDRO_33;
  v171 = DENDRO_404;
  v172 = dmul(v171, v170);
  v173 = DENDRO_33;
  v174 = DENDRO_407;
  v175 = dmul(v174, v173);
  v176 = DENDRO_33;
  v177 = DENDRO_410;
  v178 = dmul(v177, v176);
  v179 = DENDRO_60;
  v180 = DENDRO_348;
  v181 = DENDRO_351;
  v182 = dmul(v181, negone);
  v183 = dadd(v182, v180);
  v184 = dmul(v183, v179);
  v185 = dmul(v184, negone);
  v186 = dadd(v185, v178);
  v187 = dadd(v186, v175);
  v188 = dadd(v187, v172);
  v189 = dmul(v188, v169);
  v190 = DENDRO_880;
  v191 = DENDRO_412;
  v192 = DENDRO_338;
  v193 = DENDRO_881;
  v194 = dmul(v193, v192);
  v195 = dadd(v194, v191);
  v196 = dmul(v195, v190);
  v197 = dadd(v196, v189);
  v198 = dadd(v197, v168);
  v199 = dadd(v198, v161);
  v200 = dmul(v199, v159);
  v201 = dmul(v200, v158);
  v202 = dmul(v201, negone);
  v203 = DENDRO_883;
  v204 = chi[pp];
  v205 = DENDRO_419;
  v206 = dmul(v205, negone);
  v207 = DENDRO_58;
  v208 = DENDRO_428;
  v209 = DENDRO_33;
  v210 = DENDRO_480;
  v211 = dmul(v210, v209);
  v212 = DENDRO_33;
  v213 = DENDRO_481;
  v214 = dmul(v213, v212);
  v215 = DENDRO_60;
  v216 = DENDRO_423;
  v217 = DENDRO_427;
  v218 = dmul(v217, negone);
  v219 = dadd(v218, v216);
  v220 = dmul(v219, v215);
  v221 = dmul(v220, negone);
  v222 = dadd(v221, v214);
  v223 = dadd(v222, v211);
  v224 = dadd(v223, v208);
  v225 = dmul(v224, v207);
  v226 = DENDRO_638;
  v227 = DENDRO_482;
  v228 = DENDRO_421;
  v229 = DENDRO_882;
  v230 = dmul(v229, v228);
  v231 = dadd(v230, v227);
  v232 = dmul(v231, v226);
  v233 = DENDRO_732;
  v234 = DENDRO_479;
  v235 = DENDRO_421;
  v236 = DENDRO_884;
  v237 = dmul(v236, v235);
  v238 = dadd(v237, v234);
  v239 = dmul(v238, v233);
  v240 = dadd(v239, v232);
  v241 = dadd(v240, v225);
  v242 = dadd(v241, v206);
  v243 = dmul(v242, v204);
  v244 = dmul(v243, v203);
  v245 = dmul(v244, negone);
  v246 = DENDRO_885;
  v247 = chi[pp];
  v248 = DENDRO_57;
  v249 = dmul(v248, negone);
  v250 = DENDRO_732;
  v251 = DENDRO_499;
  v252 = DENDRO_59;
  v253 = DENDRO_884;
  v254 = dmul(v253, v252);
  v255 = dadd(v254, v251);
  v256 = dmul(v255, v250);
  v257 = DENDRO_880;
  v258 = DENDRO_500;
  v259 = DENDRO_59;
  v260 = DENDRO_881;
  v261 = dmul(v260, v259);
  v262 = dadd(v261, v258);
  v263 = dmul(v262, v257);
  v264 = DENDRO_94;
  v265 = DENDRO_100;
  v266 = DENDRO_33;
  v267 = DENDRO_501;
  v268 = dmul(v267, v266);
  v269 = DENDRO_33;
  v270 = DENDRO_502;
  v271 = dmul(v270, v269);
  v272 = DENDRO_60;
  v273 = DENDRO_101;
  v274 = DENDRO_106;
  v275 = dmul(v274, negone);
  v276 = dadd(v275, v273);
  v277 = dmul(v276, v272);
  v278 = dmul(v277, negone);
  v279 = dadd(v278, v271);
  v280 = dadd(v279, v268);
  v281 = dadd(v280, v265);
  v282 = dmul(v281, v264);
  v283 = dadd(v282, v263);
  v284 = dadd(v283, v256);
  v285 = dadd(v284, v249);
  v286 = dmul(v285, v247);
  v287 = dmul(v286, v246);
  v288 = dmul(v287, negone);
  v289 = DENDRO_890;
  v290 = chi[pp];
  v291 = DENDRO_503;
  v292 = dmul(v291, negone);
  v293 = DENDRO_888;
  v294 = DENDRO_508;
  v295 = DENDRO_33;
  v296 = DENDRO_523;
  v297 = dmul(v296, v295);
  v298 = DENDRO_33;
  v299 = DENDRO_524;
  v300 = dmul(v299, v298);
  v301 = DENDRO_60;
  v302 = DENDRO_511;
  v303 = DENDRO_63;
  v304 = dmul(v303, negone);
  v305 = dadd(v304, v302);
  v306 = dmul(v305, v301);
  v307 = dadd(v306, v300);
  v308 = dadd(v307, v297);
  v309 = dadd(v308, v294);
  v310 = dmul(v309, v293);
  v311 = DENDRO_891;
  v312 = DENDRO_516;
  v313 = DENDRO_33;
  v314 = DENDRO_528;
  v315 = dmul(v314, v313);
  v316 = DENDRO_33;
  v317 = DENDRO_529;
  v318 = dmul(v317, v316);
  v319 = DENDRO_60;
  v320 = DENDRO_518;
  v321 = DENDRO_62;
  v322 = dmul(v321, negone);
  v323 = dadd(v322, v320);
  v324 = dmul(v323, v319);
  v325 = dadd(v324, v318);
  v326 = dadd(v325, v315);
  v327 = dadd(v326, v312);
  v328 = dmul(v327, v311);
  v329 = DENDRO_58;
  v330 = DENDRO_887;
  v331 = DENDRO_526;
  v332 = DENDRO_881;
  v333 = gt2[pp];
  v334 = dmul(v333, v332);
  v335 = dadd(v334, v331);
  v336 = dmul(v335, v330);
  v337 = dmul(v336, v329);
  v338 = dadd(v337, v328);
  v339 = dadd(v338, v310);
  v340 = dadd(v339, v292);
  v341 = dmul(v340, v290);
  v342 = dmul(v341, v289);
  v343 = dmul(v342, negone);
  v344 = dadd(v343, v288);
  v345 = dadd(v344, v245);
  v346 = dadd(v345, v202);
  v347 = dadd(v346, v157);
  v348 = dadd(v347, v103);
  v349 = dadd(v348, v49);
  v350 = dadd(v349, v8);
  v351 = dadd(v350, v5);
  v352 = dadd(v351, v2);
  K_rhs[pp] = v352;
//--
  v0 = DENDRO_961;
  v1 = lambda[2];
  v2 = beta0[pp];
  v3 = *(agrad_0_B0+pp );
  v4 = dmul(v3, v2);
  v5 = beta1[pp];
  v6 = *(agrad_1_B0+pp );
  v7 = dmul(v6, v5);
  v8 = beta2[pp];
  v9 = *(agrad_2_B0+pp );
  v10 = dmul(v9, v8);
  v11 = dadd(v10, v7);
  v12 = dadd(v11, v4);
  v13 = dmul(v12, v1);
  v14 = B0[pp];
  v15 = eta;
  v16 = dmul(v15, v14);
  v17 = dmul(v16, negone);
  v18 = DENDRO_925;
  v19 = lambda[3];
  v20 = dmul(v19, v18);
  v21 = dmul(v20, negone);
  v22 = dadd(v21, v17);
  v23 = dadd(v22, v13);
  v24 = dadd(v23, v0);
  B_rhs0[pp] = v24;
//--
  v0 = DENDRO_10;
  v1 = DENDRO_4;
  v2 = DENDRO_977;
  v3 = dmul(v2, v1);
  v4 = DENDRO_940;
  v5 = DENDRO_972;
  v6 = dmul(v5, v4);
  v7 = DENDRO_941;
  v8 = DENDRO_975;
  v9 = dmul(v8, negone);
  v10 = DENDRO_952;
  v11 = DENDRO_971;
  v12 = dmul(v11, v10);
  v13 = dadd(v12, v9);
  v14 = dmul(v13, v7);
  v15 = DENDRO_941;
  v16 = DENDRO_976;
  v17 = dmul(v16, negone);
  v18 = DENDRO_943;
  v19 = DENDRO_972;
  v20 = dmul(v19, v18);
  v21 = dadd(v20, v17);
  v22 = dmul(v21, v15);
  v23 = DENDRO_948;
  v24 = DENDRO_971;
  v25 = dmul(v24, v23);
  v26 = DENDRO_827;
  v27 = DENDRO_946;
  v28 = dmul(v27, v26);
  v29 = dmul(v28, negone);
  v30 = DENDRO_960;
  v31 = DENDRO_977;
  v32 = dmul(v31, v30);
  v33 = dmul(v32, negone);
  v34 = DENDRO_971;
  v35 = DENDRO_974;
  v36 = dmul(v35, v34);
  v37 = dmul(v36, negone);
  v38 = DENDRO_972;
  v39 = DENDRO_973;
  v40 = dmul(v39, v38);
  v41 = dmul(v40, negone);
  v42 = DENDRO_124;
  v43 = DENDRO_16;
  v44 = DENDRO_287;
  v45 = dmul(v44, v43);
  v46 = dmul(v45, v42);
  v47 = DENDRO_124;
  v48 = DENDRO_307;
  v49 = DENDRO_9;
  v50 = dmul(v49, v48);
  v51 = dmul(v50, v47);
  v52 = dadd(v51, v46);
  v53 = dadd(v52, v41);
  v54 = dadd(v53, v37);
  v55 = dadd(v54, v33);
  v56 = dadd(v55, v29);
  v57 = dadd(v56, v25);
  v58 = dadd(v57, v22);
  v59 = dadd(v58, v14);
  v60 = dadd(v59, v6);
  v61 = dadd(v60, v3);
  v62 = dadd(v61, v0);
  Gt_rhs1[pp] = v62;
//--
  v0 = DENDRO_961;
  Gt_rhs0[pp] = v0;
//--
  v0 = DENDRO_980;
  Gt_rhs2[pp] = v0;
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 191
// Dendro: printing temp variables

// Dendro: printing variables
//--
chi_rhs[pp] = DENDRO_23*K[pp]*alpha[pp] - DENDRO_23*(DENDRO_2 + DENDRO_4 + DENDRO_6) + beta0[pp]*agrad(0, chi[pp]) + beta1[pp]*agrad(1, chi[pp]) + beta2[pp]*agrad(2, chi[pp]);
//--
gt_rhs02[pp] = -At2[pp]*DENDRO_0 + DENDRO_11*gt5[pp] + DENDRO_15*gt0[pp] + DENDRO_16*gt1[pp] + DENDRO_17*DENDRO_2 + DENDRO_17*DENDRO_6 - DENDRO_5*gt2[pp] + DENDRO_9*gt4[pp] + beta0[pp]*agrad(0, gt2[pp]) + beta1[pp]*agrad(1, gt2[pp]) + beta2[pp]*agrad(2, gt2[pp]);
//--
gt_rhs11[pp] = -At3[pp]*DENDRO_0 + DENDRO_12*DENDRO_8 + DENDRO_13*DENDRO_20 - DENDRO_18*gt3[pp] + DENDRO_19*gt3[pp] - DENDRO_7*gt3[pp] + beta0[pp]*agrad(0, gt3[pp]) + beta1[pp]*agrad(1, gt3[pp]) + beta2[pp]*agrad(2, gt3[pp]);
//--
gt_rhs12[pp] = -At4[pp]*DENDRO_0 + DENDRO_12*gt2[pp] + DENDRO_13*gt5[pp] + DENDRO_15*gt1[pp] + DENDRO_16*gt3[pp] - DENDRO_18*gt4[pp] + DENDRO_21*DENDRO_4 + DENDRO_21*DENDRO_6 + beta0[pp]*agrad(0, gt4[pp]) + beta1[pp]*agrad(1, gt4[pp]) + beta2[pp]*agrad(2, gt4[pp]);
//--
b_rhs2[pp] = B2[pp]*DENDRO_1 + lambda[1]*(beta0[pp]*agrad(0, beta2[pp]) + beta1[pp]*agrad(1, beta2[pp]) + beta2[pp]*agrad(2, beta2[pp]));
//--
b_rhs0[pp] = B0[pp]*DENDRO_1 + lambda[1]*(beta0[pp]*agrad(0, beta0[pp]) + beta1[pp]*agrad(1, beta0[pp]) + beta2[pp]*agrad(2, beta0[pp]));
//--
gt_rhs01[pp] = -At1[pp]*DENDRO_0 + DENDRO_11*gt4[pp] + DENDRO_12*gt0[pp] + DENDRO_13*gt2[pp] + DENDRO_14*DENDRO_2 + DENDRO_14*DENDRO_4 - DENDRO_7*gt1[pp] + DENDRO_9*gt3[pp] + beta0[pp]*agrad(0, gt1[pp]) + beta1[pp]*agrad(1, gt1[pp]) + beta2[pp]*agrad(2, gt1[pp]);
//--
b_rhs1[pp] = B1[pp]*DENDRO_1 + lambda[1]*(beta0[pp]*agrad(0, beta1[pp]) + beta1[pp]*agrad(1, beta1[pp]) + beta2[pp]*agrad(2, beta1[pp]));
//--
gt_rhs22[pp] = -At5[pp]*DENDRO_0 + DENDRO_10*DENDRO_15 + DENDRO_16*DENDRO_20 - DENDRO_18*gt5[pp] + DENDRO_22*gt5[pp] - DENDRO_5*gt5[pp] + beta0[pp]*agrad(0, gt5[pp]) + beta1[pp]*agrad(1, gt5[pp]) + beta2[pp]*agrad(2, gt5[pp]);
//--
gt_rhs00[pp] = -At0[pp]*DENDRO_0 + DENDRO_10*DENDRO_11 + DENDRO_3*gt0[pp] - DENDRO_5*gt0[pp] - DENDRO_7*gt0[pp] + DENDRO_8*DENDRO_9 + beta0[pp]*agrad(0, gt0[pp]) + beta1[pp]*agrad(1, gt0[pp]) + beta2[pp]*agrad(2, gt0[pp]);
//--
a_rhs[pp] = -DENDRO_0*K[pp] + lambda[0]*(beta0[pp]*agrad(0, alpha[pp]) + beta1[pp]*agrad(1, alpha[pp]) + beta2[pp]*agrad(2, alpha[pp]));
// Dendro: reduced ops: 191
// Dendro: }}} 
// Dendro vectorized code: {{{
//--
  double v0 = beta0[pp];
  double v1 = *(agrad_0_chi+pp );
  double v2 = dmul(v1, v0);
  double v3 = beta1[pp];
  double v4 = *(agrad_1_chi+pp );
  double v5 = dmul(v4, v3);
  double v6 = beta2[pp];
  double v7 = *(agrad_2_chi+pp );
  double v8 = dmul(v7, v6);
  double v9 = DENDRO_23;
  double v10 = DENDRO_2;
  double v11 = DENDRO_4;
  double v12 = DENDRO_6;
  double v13 = dadd(v12, v11);
  double v14 = dadd(v13, v10);
  double v15 = dmul(v14, v9);
  double v16 = dmul(v15, negone);
  double v17 = DENDRO_23;
  double v18 = K[pp];
  double v19 = alpha[pp];
  double v20 = dmul(v19, v18);
  double v21 = dmul(v20, v17);
  double v22 = dadd(v21, v16);
  double v23 = dadd(v22, v8);
  double v24 = dadd(v23, v5);
  double v25 = dadd(v24, v2);
  chi_rhs[pp] = v25;
//--
  v0 = DENDRO_11;
  v1 = gt5[pp];
  v2 = dmul(v1, v0);
  v3 = DENDRO_15;
  v4 = gt0[pp];
  v5 = dmul(v4, v3);
  v6 = DENDRO_16;
  v7 = gt1[pp];
  v8 = dmul(v7, v6);
  v9 = DENDRO_17;
  v10 = DENDRO_2;
  v11 = dmul(v10, v9);
  v12 = DENDRO_17;
  v13 = DENDRO_6;
  v14 = dmul(v13, v12);
  v15 = DENDRO_9;
  v16 = gt4[pp];
  v17 = dmul(v16, v15);
  v18 = beta0[pp];
  v19 = *(agrad_0_gt2+pp );
  v20 = dmul(v19, v18);
  v21 = beta1[pp];
  v22 = *(agrad_1_gt2+pp );
  v23 = dmul(v22, v21);
  v24 = beta2[pp];
  v25 = *(agrad_2_gt2+pp );
  double v26 = dmul(v25, v24);
  double v27 = At2[pp];
  double v28 = DENDRO_0;
  double v29 = dmul(v28, v27);
  double v30 = dmul(v29, negone);
  double v31 = DENDRO_5;
  double v32 = gt2[pp];
  double v33 = dmul(v32, v31);
  double v34 = dmul(v33, negone);
  double v35 = dadd(v34, v30);
  double v36 = dadd(v35, v26);
  double v37 = dadd(v36, v23);
  double v38 = dadd(v37, v20);
  double v39 = dadd(v38, v17);
  double v40 = dadd(v39, v14);
  double v41 = dadd(v40, v11);
  double v42 = dadd(v41, v8);
  double v43 = dadd(v42, v5);
  double v44 = dadd(v43, v2);
  gt_rhs02[pp] = v44;
//--
  v0 = DENDRO_12;
  v1 = DENDRO_8;
  v2 = dmul(v1, v0);
  v3 = DENDRO_13;
  v4 = DENDRO_20;
  v5 = dmul(v4, v3);
  v6 = DENDRO_19;
  v7 = gt3[pp];
  v8 = dmul(v7, v6);
  v9 = beta0[pp];
  v10 = *(agrad_0_gt3+pp );
  v11 = dmul(v10, v9);
  v12 = beta1[pp];
  v13 = *(agrad_1_gt3+pp );
  v14 = dmul(v13, v12);
  v15 = beta2[pp];
  v16 = *(agrad_2_gt3+pp );
  v17 = dmul(v16, v15);
  v18 = At3[pp];
  v19 = DENDRO_0;
  v20 = dmul(v19, v18);
  v21 = dmul(v20, negone);
  v22 = DENDRO_18;
  v23 = gt3[pp];
  v24 = dmul(v23, v22);
  v25 = dmul(v24, negone);
  v26 = DENDRO_7;
  v27 = gt3[pp];
  v28 = dmul(v27, v26);
  v29 = dmul(v28, negone);
  v30 = dadd(v29, v25);
  v31 = dadd(v30, v21);
  v32 = dadd(v31, v17);
  v33 = dadd(v32, v14);
  v34 = dadd(v33, v11);
  v35 = dadd(v34, v8);
  v36 = dadd(v35, v5);
  v37 = dadd(v36, v2);
  gt_rhs11[pp] = v37;
//--
  v0 = DENDRO_12;
  v1 = gt2[pp];
  v2 = dmul(v1, v0);
  v3 = DENDRO_13;
  v4 = gt5[pp];
  v5 = dmul(v4, v3);
  v6 = DENDRO_15;
  v7 = gt1[pp];
  v8 = dmul(v7, v6);
  v9 = DENDRO_16;
  v10 = gt3[pp];
  v11 = dmul(v10, v9);
  v12 = DENDRO_21;
  v13 = DENDRO_4;
  v14 = dmul(v13, v12);
  v15 = DENDRO_21;
  v16 = DENDRO_6;
  v17 = dmul(v16, v15);
  v18 = beta0[pp];
  v19 = *(agrad_0_gt4+pp );
  v20 = dmul(v19, v18);
  v21 = beta1[pp];
  v22 = *(agrad_1_gt4+pp );
  v23 = dmul(v22, v21);
  v24 = beta2[pp];
  v25 = *(agrad_2_gt4+pp );
  v26 = dmul(v25, v24);
  v27 = At4[pp];
  v28 = DENDRO_0;
  v29 = dmul(v28, v27);
  v30 = dmul(v29, negone);
  v31 = DENDRO_18;
  v32 = gt4[pp];
  v33 = dmul(v32, v31);
  v34 = dmul(v33, negone);
  v35 = dadd(v34, v30);
  v36 = dadd(v35, v26);
  v37 = dadd(v36, v23);
  v38 = dadd(v37, v20);
  v39 = dadd(v38, v17);
  v40 = dadd(v39, v14);
  v41 = dadd(v40, v11);
  v42 = dadd(v41, v8);
  v43 = dadd(v42, v5);
  v44 = dadd(v43, v2);
  gt_rhs12[pp] = v44;
//--
  v0 = B2[pp];
  v1 = DENDRO_1;
  v2 = dmul(v1, v0);
  v3 = lambda[1];
  v4 = beta0[pp];
  v5 = *(agrad_0_beta2+pp );
  v6 = dmul(v5, v4);
  v7 = beta1[pp];
  v8 = *(agrad_1_beta2+pp );
  v9 = dmul(v8, v7);
  v10 = beta2[pp];
  v11 = *(agrad_2_beta2+pp );
  v12 = dmul(v11, v10);
  v13 = dadd(v12, v9);
  v14 = dadd(v13, v6);
  v15 = dmul(v14, v3);
  v16 = dadd(v15, v2);
  b_rhs2[pp] = v16;
//--
  v0 = B0[pp];
  v1 = DENDRO_1;
  v2 = dmul(v1, v0);
  v3 = lambda[1];
  v4 = beta0[pp];
  v5 = *(agrad_0_beta0+pp );
  v6 = dmul(v5, v4);
  v7 = beta1[pp];
  v8 = *(agrad_1_beta0+pp );
  v9 = dmul(v8, v7);
  v10 = beta2[pp];
  v11 = *(agrad_2_beta0+pp );
  v12 = dmul(v11, v10);
  v13 = dadd(v12, v9);
  v14 = dadd(v13, v6);
  v15 = dmul(v14, v3);
  v16 = dadd(v15, v2);
  b_rhs0[pp] = v16;
//--
  v0 = DENDRO_11;
  v1 = gt4[pp];
  v2 = dmul(v1, v0);
  v3 = DENDRO_12;
  v4 = gt0[pp];
  v5 = dmul(v4, v3);
  v6 = DENDRO_13;
  v7 = gt2[pp];
  v8 = dmul(v7, v6);
  v9 = DENDRO_14;
  v10 = DENDRO_2;
  v11 = dmul(v10, v9);
  v12 = DENDRO_14;
  v13 = DENDRO_4;
  v14 = dmul(v13, v12);
  v15 = DENDRO_9;
  v16 = gt3[pp];
  v17 = dmul(v16, v15);
  v18 = beta0[pp];
  v19 = *(agrad_0_gt1+pp );
  v20 = dmul(v19, v18);
  v21 = beta1[pp];
  v22 = *(agrad_1_gt1+pp );
  v23 = dmul(v22, v21);
  v24 = beta2[pp];
  v25 = *(agrad_2_gt1+pp );
  v26 = dmul(v25, v24);
  v27 = At1[pp];
  v28 = DENDRO_0;
  v29 = dmul(v28, v27);
  v30 = dmul(v29, negone);
  v31 = DENDRO_7;
  v32 = gt1[pp];
  v33 = dmul(v32, v31);
  v34 = dmul(v33, negone);
  v35 = dadd(v34, v30);
  v36 = dadd(v35, v26);
  v37 = dadd(v36, v23);
  v38 = dadd(v37, v20);
  v39 = dadd(v38, v17);
  v40 = dadd(v39, v14);
  v41 = dadd(v40, v11);
  v42 = dadd(v41, v8);
  v43 = dadd(v42, v5);
  v44 = dadd(v43, v2);
  gt_rhs01[pp] = v44;
//--
  v0 = B1[pp];
  v1 = DENDRO_1;
  v2 = dmul(v1, v0);
  v3 = lambda[1];
  v4 = beta0[pp];
  v5 = *(agrad_0_beta1+pp );
  v6 = dmul(v5, v4);
  v7 = beta1[pp];
  v8 = *(agrad_1_beta1+pp );
  v9 = dmul(v8, v7);
  v10 = beta2[pp];
  v11 = *(agrad_2_beta1+pp );
  v12 = dmul(v11, v10);
  v13 = dadd(v12, v9);
  v14 = dadd(v13, v6);
  v15 = dmul(v14, v3);
  v16 = dadd(v15, v2);
  b_rhs1[pp] = v16;
//--
  v0 = DENDRO_10;
  v1 = DENDRO_15;
  v2 = dmul(v1, v0);
  v3 = DENDRO_16;
  v4 = DENDRO_20;
  v5 = dmul(v4, v3);
  v6 = DENDRO_22;
  v7 = gt5[pp];
  v8 = dmul(v7, v6);
  v9 = beta0[pp];
  v10 = *(agrad_0_gt5+pp );
  v11 = dmul(v10, v9);
  v12 = beta1[pp];
  v13 = *(agrad_1_gt5+pp );
  v14 = dmul(v13, v12);
  v15 = beta2[pp];
  v16 = *(agrad_2_gt5+pp );
  v17 = dmul(v16, v15);
  v18 = At5[pp];
  v19 = DENDRO_0;
  v20 = dmul(v19, v18);
  v21 = dmul(v20, negone);
  v22 = DENDRO_18;
  v23 = gt5[pp];
  v24 = dmul(v23, v22);
  v25 = dmul(v24, negone);
  v26 = DENDRO_5;
  v27 = gt5[pp];
  v28 = dmul(v27, v26);
  v29 = dmul(v28, negone);
  v30 = dadd(v29, v25);
  v31 = dadd(v30, v21);
  v32 = dadd(v31, v17);
  v33 = dadd(v32, v14);
  v34 = dadd(v33, v11);
  v35 = dadd(v34, v8);
  v36 = dadd(v35, v5);
  v37 = dadd(v36, v2);
  gt_rhs22[pp] = v37;
//--
  v0 = DENDRO_10;
  v1 = DENDRO_11;
  v2 = dmul(v1, v0);
  v3 = DENDRO_3;
  v4 = gt0[pp];
  v5 = dmul(v4, v3);
  v6 = DENDRO_8;
  v7 = DENDRO_9;
  v8 = dmul(v7, v6);
  v9 = beta0[pp];
  v10 = *(agrad_0_gt0+pp );
  v11 = dmul(v10, v9);
  v12 = beta1[pp];
  v13 = *(agrad_1_gt0+pp );
  v14 = dmul(v13, v12);
  v15 = beta2[pp];
  v16 = *(agrad_2_gt0+pp );
  v17 = dmul(v16, v15);
  v18 = At0[pp];
  v19 = DENDRO_0;
  v20 = dmul(v19, v18);
  v21 = dmul(v20, negone);
  v22 = DENDRO_5;
  v23 = gt0[pp];
  v24 = dmul(v23, v22);
  v25 = dmul(v24, negone);
  v26 = DENDRO_7;
  v27 = gt0[pp];
  v28 = dmul(v27, v26);
  v29 = dmul(v28, negone);
  v30 = dadd(v29, v25);
  v31 = dadd(v30, v21);
  v32 = dadd(v31, v17);
  v33 = dadd(v32, v14);
  v34 = dadd(v33, v11);
  v35 = dadd(v34, v8);
  v36 = dadd(v35, v5);
  v37 = dadd(v36, v2);
  gt_rhs00[pp] = v37;
//--
  v0 = lambda[0];
  v1 = beta0[pp];
  v2 = *(agrad_0_alpha+pp );
  v3 = dmul(v2, v1);
  v4 = beta1[pp];
  v5 = *(agrad_1_alpha+pp );
  v6 = dmul(v5, v4);
  v7 = beta2[pp];
  v8 = *(agrad_2_alpha+pp );
  v9 = dmul(v8, v7);
  v10 = dadd(v9, v6);
  v11 = dadd(v10, v3);
  v12 = dmul(v11, v0);
  v13 = DENDRO_0;
  v14 = K[pp];
  v15 = dmul(v14, v13);
  v16 = dmul(v15, negone);
  v17 = dadd(v16, v12);
  a_rhs[pp] = v17;
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
// Dendro: {{{ 
// Dendro: original ops: 0
// Dendro: printing temp variables

// Dendro: printing variables
// Dendro: reduced ops: 0
// Dendro: }}} 
// Dendro vectorized code: {{{
// Dendro vectorized code: }}} 
