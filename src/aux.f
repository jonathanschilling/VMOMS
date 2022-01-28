      SUBROUTINE AUX(YP, Y, X, PARAM)                                   ABSH0206
C***********************************************************************ABSH0207
C*****AUX EVALUATES THE DERIVATIVE OF Y(X) AND TRANSFORMS THE SYSTEM OF*ABSH0208
C*****SECOND ORDER ODES INTO A SYSTEM OF FIRST ORDER ODES - AL*YP=BL.  *ABSH0209
C*****REFERENCES:                                                      *ABSH0210
C*****L.L.LAO,S.P.HIRSHMAN,R.M.WIELAND,ORNL/TM-7616 (1981).            *ABSH0211
C*****                                ,PHYS FLUIDS 24 (1981) AUG       *ABSH0212
C*****L.L.LAO,R.M.WIELAND,W.A.HOULBERG,S.P.HIRSHMAN,ORNL/TM-7871 (1981)*ABSH0213
C*****LAST REVISION: 6/81 L.L.LAO, R.M.WIELAND, AND W.A.HOULBERG ORNL. *ABSH0214
C*****CALCULATED PARAMETERS:                                           *ABSH0215
C*****YP(I)-DERIVATIVE OF Y(I) WITH RESPECT TO X.                      *ABSH0216
C*****INPUT PARAMETERS:                                                *ABSH0217
C*****Y(I)-SOLUTION OF FIRST ORDER ODE I AT X.                         *ABSH0218
C*****X-ABSCISSA AT WHICH YP(I) ARE TO BE EVALUATED.                   *ABSH0219
C*****OTHER COMMENTS:                                                  *ABSH0220
C*****AUX IS REQUIRED BY D02AGF ODE SOLVER IN NAG LIBRARY.             *ABSH0221
C***********************************************************************ABSH0222
      PARAMETER MNMX=61                                                 ABSH0223
      PARAMETER MTHETA=11,MM1=10,MNMOMS=3                               ABSH0224
      PARAMETER MNEQ=6                                                  ABSH0225
      COMMON /CFEVAL/ NFL                                               ABSH0226
      COMMON /LOCMEQ/ XL0, DX, R0, A0, BT0, P0, AI0, E0, E1,            ABSH0227
     *     D0, ARC, SLPE, SLPR0, BETAP0, KMHD, IUNDER                   ABSH0228
      COMMON /NIO/ NIN, NOUT, NPLOT, NDEBUG, NTTY                       ABSH0229
      COMMON /SCFUNC/ S1TH(MTHETA), S2TH(MTHETA), C1TH(MTHETA),         ABSH0230
     *     C2TH(MTHETA)                                                 ABSH0231
      COMMON /SWORK2/ BSIFT(MM1), CSIFT(MM1), DSIFT(MM1),               ABSH0232
     *     BELLI(MM1), CELLI(MM1), DELLI(MM1), BTRIA(MM1),              ABSH0233
     *     CTRIA(MM1), DTRIA(MM1), R(MTHETA), RX(MTHETA),               ABSH0234
     *     RT(MTHETA), RXT(MTHETA), RTT(MTHETA), Z(MTHETA),             ABSH0235
     *     ZX(MTHETA), ZT(MTHETA), ZXT(MTHETA), ZTT(MTHETA),            ABSH0236
     *     GTT(MTHETA), TAU(MTHETA), GSQRT(MTHETA),                     ABSH0237
     *     GTTX(MTHETA), ZMOM1(MTHETA), RMOM0(MTHETA),                  ABSH0238
     *     TMOM(MTHETA), RMOM2(MTHETA), RTAU2(MTHETA),                  ABSH0239
     *     RTAU3(MTHETA), GS10(MTHETA), GS21(MTHETA),                   ABSH0240
     *     GS4(MTHETA), GS5(MTHETA), GS6(MTHETA), GS7(MTHETA),          ABSH0241
     *     GS8(MTHETA), GS2(MTHETA), GS12(MTHETA),                      ABSH0242
     *     GS22(MTHETA), CS13N(MTHETA), SS20(MTHETA),                   ABSH0243
     *     CS1N(MTHETA)                                                 ABSH0244
      COMMON /WEIGHT/ WFCN(MTHETA), CNORM                               ABSH0245
      DIMENSION YP(MNEQ), Y(MNEQ), PARAM(MNEQ)                          ABSH0246
      DATA NFL /0/                                                      ABSH0247
      NFL = NFL + 1                                                     ABSH0248
      NTHETA = MTHETA                                                   ABSH0249
      CALL PRSFUN(X, PX, PPX)                                           ABSH0250
      CALL CURFUN(X, AIX, AIPX)                                         ABSH0251
      AIPOI = AIPX/AIX                                                  ABSH0252
      AI = BETAP0/2.0*PPX/AIX**2                                        ABSH0253
C*****VECTORIZE THE COEFFICIENTS IN THE MOMENT EQUATIONS.               ABSH0254
      DO 10 I=1,NTHETA                                                  ABSH0255
           R(I) = Y(1) - X*C1TH(I) + Y(5)*C2TH(I)                       ABSH0256
           RX(I) = Y(2) - C1TH(I) + Y(6)*C2TH(I)                        ABSH0257
           RT(I) = X*S1TH(I) - 2.0*Y(5)*S2TH(I)                         ABSH0258
           RXT(I) = S1TH(I) - 2.0*Y(6)*S2TH(I)                          ABSH0259
           RTT(I) = X*C1TH(I) - 4.0*Y(5)*C2TH(I)                        ABSH0260
           Z(I) = Y(3)*(X*S1TH(I)+Y(5)*S2TH(I))                         ABSH0261
           ZX(I) = Y(4)*(X*S1TH(I)+Y(5)*S2TH(I)) +                      ABSH0262
     *          Y(3)*(S1TH(I)+Y(6)*S2TH(I))                             ABSH0263
           ZT(I) = Y(3)*(X*C1TH(I)+2.0*Y(5)*C2TH(I))                    ABSH0264
           ZXT(I) = Y(4)*(X*C1TH(I)+2.0*Y(5)*C2TH(I)) +                 ABSH0265
     *          Y(3)*(C1TH(I)+2.0*Y(6)*C2TH(I))                         ABSH0266
           ZTT(I) = -Y(3)*(X*S1TH(I)+4.0*Y(5)*S2TH(I))                  ABSH0267
           GTT(I) = RT(I)**2 + ZT(I)**2                                 ABSH0268
           TAU(I) = RT(I)*ZX(I) - RX(I)*ZT(I)                           ABSH0269
           GSQRT(I) = R(I)*TAU(I)                                       ABSH0270
           RTAU2(I) = GSQRT(I)*TAU(I)                                   ABSH0271
           RTAU3(I) = RTAU2(I)*TAU(I)                                   ABSH0272
           GS10(I) = ZT(I)*GTT(I)/RTAU3(I)                              ABSH0273
           GS12(I) = GS10(I)*C2TH(I)                                    ABSH0274
           GS2(I) = RT(I)*GTT(I)/RTAU3(I)                               ABSH0275
           GS21(I) = GS2(I)*S1TH(I)                                     ABSH0276
           GS22(I) = GS2(I)*S2TH(I)                                     ABSH0277
           GS4(I) = 1.0/R(I)                                            ABSH0278
           GS5(I) = ZT(I)/R(I)/GSQRT(I)                                 ABSH0279
           CS13N(I) = GTT(I)*(RX(I)*ZXT(I)-ZX(I)*RXT(I))                ABSH0280
           GTTX(I) = RT(I)*RXT(I) + ZT(I)*ZXT(I)                        ABSH0281
           GS6(I) = (CS13N(I)+(RX(I)*RT(I)+ZX(I)*ZT(I))*(RTT(I)*        ABSH0282
     *          ZX(I)+RT(I)*ZXT(I)-RXT(I)*ZT(I)-RX(I)*ZTT(I)))/         ABSH0283
     *          RTAU3(I)                                                ABSH0284
           GS7(I) = (GTTX(I)-RX(I)*RTT(I)-ZX(I)*ZTT(I))/RTAU2(I)        ABSH0285
           GS8(I) = GTT(I)/RTAU2(I)                                     ABSH0286
           SS20(I) = GS5(I) + GS6(I) + GS7(I)                           ABSH0287
           CS1N(I) = (2.0*GTTX(I)-GTT(I)*RX(I)/R(I)+CS13N(I)/           ABSH0288
     *          TAU(I))/GSQRT(I)                                        ABSH0289
           RMOM0(I) = ZT(I)*WFCN(I)                                     ABSH0290
           TMOM(I) = TAU(I)*WFCN(I)                                     ABSH0291
           RMOM2(I) = (Y(3)*RT(I)*S2TH(I)-ZT(I)*C2TH(I))*WFCN(I)        ABSH0292
C          ZMOM1(I) = (-X*S1TH(I)+Y(5)*S2TH(I))*RT(I)*WFCN(I)           ABSH0293
           ZMOM1(I) = (+X*S1TH(I)+Y(5)*S2TH(I))*RT(I)*WFCN(I)           ABSH002A
   10 CONTINUE                                                          ABSH0294
      YP(1) = Y(2)                                                      ABSH0295
      YP(3) = Y(4)                                                      ABSH0296
      YP(5) = Y(6)                                                      ABSH0297
C*****R0 MOMENT TERMS - SHIFT.                                          ABSH0298
      CR010 = SDOT(NTHETA,RMOM0,1,GS10,1)                               ABSH0299
      CR012 = SDOT(NTHETA,RMOM0,1,GS12,1)                               ABSH0300
      CR021 = SDOT(NTHETA,RMOM0,1,GS21,1)                               ABSH0301
      CR022 = SDOT(NTHETA,RMOM0,1,GS22,1)                               ABSH0302
      CR03 = SDOT(NTHETA,RMOM0,1,R,1)                                   ABSH0303
      CR04 = SDOT(NTHETA,RMOM0,1,GS4,1)                                 ABSH0304
      CR08 = SDOT(NTHETA,RMOM0,1,GS8,1)                                 ABSH0305
      CR0SS2 = SDOT(NTHETA,RMOM0,1,SS20,1)                              ABSH0306
      CPSI10 = SDOT(NTHETA,TMOM,1,GS10,1)                               ABSH0307
      CPSI12 = SDOT(NTHETA,TMOM,1,GS12,1)                               ABSH0308
      CPSI21 = SDOT(NTHETA,TMOM,1,GS21,1)                               ABSH0309
      CPSI22 = SDOT(NTHETA,TMOM,1,GS22,1)                               ABSH0310
      CPSI3 = SDOT(NTHETA,TMOM,1,R,1)                                   ABSH0311
      CPSI4 = SDOT(NTHETA,TMOM,1,GS4,1)                                 ABSH0312
      CPSI8 = SDOT(NTHETA,TMOM,1,GS8,1)                                 ABSH0313
      CSISS2 = SDOT(NTHETA,TMOM,1,SS20,1)                               ABSH0314
      CCS1 = SDOT(NTHETA,WFCN,1,CS1N,1)/CPSI8                           ABSH0315
      CR0GS4 = CR04/CPSI4                                               ABSH0316
      CR0GS8 = CR08/CPSI8                                               ABSH0317
      DR010 = CR010 - CR0GS8*CPSI10                                     ABSH0318
      DR012 = CR012 - CR0GS8*CPSI12                                     ABSH0319
      DR021 = CR021 - CR0GS8*CPSI21                                     ABSH0320
      DR022 = CR022 - CR0GS8*CPSI22                                     ABSH0321
      DR03 = CPSI8**2/CNORM**2*(CR03-CR0GS4*CPSI3)                      ABSH0322
      DR040 = CR0SS2 - CR0GS4*CSISS2                                    ABSH0323
      DR041 = CR08 - CR0GS4*CPSI8                                       ABSH0324
      DR05 = CCS1*(CR08-CR0GS4*CPSI8)                                   ABSH0325
      AL11 = DR010                                                      ABSH0326
      AL12 = -(X*DR021+Y(5)*DR022)                                      ABSH0327
      AL13 = DR012 - Y(3)*DR022                                         ABSH0328
      BL1 = -AI*DR03 - DR040 - DR041*AIPOI + DR05 +                     ABSH0329
     *     2.0*Y(4)*(DR021+Y(6)*DR022)                                  ABSH0330
      IF (KMHD.GT.1) GO TO 20                                           ABSH0331
C*****SOLUTION OF MATRIX FOR KMHD=1.                                    ABSH0332
      B3 = 2.0*D0                                                       ABSH0333
      B2 = 2.0*(E0-E1)                                                  ABSH0334
      B1 = (BL1-AL13*BL3-AL12*BL2)/AL11                                 ABSH0335
      GO TO 40                                                          ABSH0336
C*****Z0 MOMENT TERMS - ELONGATION.                                     ABSH0337
   20 CZ110 = SDOT(NTHETA,ZMOM1,1,GS10,1)                               ABSH0338
      CZ112 = SDOT(NTHETA,ZMOM1,1,GS12,1)                               ABSH0339
      CZ121 = SDOT(NTHETA,ZMOM1,1,GS21,1)                               ABSH0340
      CZ122 = SDOT(NTHETA,ZMOM1,1,GS22,1)                               ABSH0341
      CZ13 = SDOT(NTHETA,ZMOM1,1,R,1)                                   ABSH0342
      CZ14 = SDOT(NTHETA,ZMOM1,1,GS4,1)                                 ABSH0343
      CZ18 = SDOT(NTHETA,ZMOM1,1,GS8,1)                                 ABSH0344
      CZ1SS2 = SDOT(NTHETA,ZMOM1,1,SS20,1)                              ABSH0345
      CZ1GS4 = CZ14/CPSI4                                               ABSH0346
      CZ1GS8 = CZ18/CPSI8                                               ABSH0347
      DZ110 = CZ110 - CZ1GS8*CPSI10                                     ABSH0348
      DZ112 = CZ112 - CZ1GS8*CPSI12                                     ABSH0349
      DZ121 = CZ121 - CZ1GS8*CPSI21                                     ABSH0350
      DZ122 = CZ122 - CZ1GS8*CPSI22                                     ABSH0351
      DZ13 = CPSI8**2/CNORM**2*(CZ13-CZ1GS4*CPSI3)                      ABSH0352
      DZ140 = CZ1SS2 - CZ1GS4*CSISS2                                    ABSH0353
      DZ141 = CZ18 - CZ1GS4*CPSI8                                       ABSH0354
      DZ15 = CCS1*(CZ18-CZ1GS4*CPSI8)                                   ABSH0355
      AL21 = DZ110                                                      ABSH0356
      AL22 = -(X*DZ121+Y(5)*DZ122)                                      ABSH0357
      AL23 = DZ112 - Y(3)*DZ122                                         ABSH0358
      BL2 = -AI*DZ13 - DZ140 - DZ141*AIPOI + DZ15 +                     ABSH0359
     *     2.0*Y(4)*(DZ121+Y(6)*DZ122)                                  ABSH0360
      IF (KMHD.EQ.3) GO TO 30                                           ABSH0361
C*****SOLUTION OF MATRIX FOR KMHD=2.                                    ABSH0362
      B3 = 2.0*D0                                                       ABSH0363
      BL1 = BL1 - AL13*BL3                                              ABSH0364
      BL2 = BL2 - AL23*BL3                                              ABSH0365
      DET = AL11*AL22 - AL12*AL21                                       ABSH0366
      IF (DET.EQ.0.0) GO TO 50                                          ABSH0367
      B1 = (BL1*AL22-AL12*BL2)/DET                                      ABSH0368
      B2 = (AL11*BL2-BL1*AL21)/DET                                      ABSH0369
      GO TO 40                                                          ABSH0370
C*****R2 MOMENT TERMS - TRIANGULARITY.                                  ABSH0371
   30 CR210 = SDOT(NTHETA,RMOM2,1,GS10,1)                               ABSH0372
      CR212 = SDOT(NTHETA,RMOM2,1,GS12,1)                               ABSH0373
      CR221 = SDOT(NTHETA,RMOM2,1,GS21,1)                               ABSH0374
      CR222 = SDOT(NTHETA,RMOM2,1,GS22,1)                               ABSH0375
      CR23 = SDOT(NTHETA,RMOM2,1,R,1)                                   ABSH0376
      CR24 = SDOT(NTHETA,RMOM2,1,GS4,1)                                 ABSH0377
      CR28 = SDOT(NTHETA,RMOM2,1,GS8,1)                                 ABSH0378
      CR2SS2 = SDOT(NTHETA,RMOM2,1,SS20,1)                              ABSH0379
      CR2GS4 = CR24/CPSI4                                               ABSH0380
      CR2GS8 = CR28/CPSI8                                               ABSH0381
      DR210 = CR210 - CR2GS8*CPSI10                                     ABSH0382
      DR212 = CR212 - CR2GS8*CPSI12                                     ABSH0383
      DR221 = CR221 - CR2GS8*CPSI21                                     ABSH0384
      DR222 = CR222 - CR2GS8*CPSI22                                     ABSH0385
      DR23 = CPSI8**2/CNORM**2*(CR23-CR2GS4*CPSI3)                      ABSH0386
      DR240 = CR2SS2 - CR2GS4*CSISS2                                    ABSH0387
      DR241 = CR28 - CR2GS4*CPSI8                                       ABSH0388
      DR25 = CCS1*(CR28-CR2GS4*CPSI8)                                   ABSH0389
      AL31 = DR210                                                      ABSH0390
      AL32 = -(X*DR221+Y(5)*DR222)                                      ABSH0391
      AL33 = DR212 - Y(3)*DR222                                         ABSH0392
      BL3 = -AI*DR23 - DR240 - DR241*AIPOI + DR25 +                     ABSH0393
     *     2.0*Y(4)*(DR221+Y(6)*DR222)                                  ABSH0394
C*****SOLUTION OF MATRIX FOR KMHD=3.                                    ABSH0395
      DET1 = AL22*AL33 - AL23*AL32                                      ABSH0396
      DET2 = AL21*AL33 - AL23*AL31                                      ABSH0397
      DET3 = AL21*AL32 - AL22*AL31                                      ABSH0398
      DET = AL11*DET1 - AL12*DET2 + AL13*DET3                           ABSH0399
      IF (DET.EQ.0.0) GO TO 50                                          ABSH0400
      B1 = BL1*DET1 - AL12*(BL2*AL33-AL23*BL3) +                        ABSH0401
     *     AL13*(BL2*AL32-AL22*BL3)                                     ABSH0402
      B1 = B1/DET                                                       ABSH0403
      B2 = AL11*(BL2*AL33-AL23*BL3) - BL1*DET2 +                        ABSH0404
     *     AL13*(AL21*BL3-BL2*AL31)                                     ABSH0405
      B2 = B2/DET                                                       ABSH0406
      B3 = AL11*(AL22*BL3-BL2*AL32) - AL12*(AL21*BL3-BL2*AL31)          ABSH0407
     *     + BL1*DET3                                                   ABSH0408
      B3 = B3/DET                                                       ABSH0409
   40 YP(2) = B1                                                        ABSH0410
      YP(4) = B2                                                        ABSH0411
      YP(6) = B3                                                        ABSH0412
      IUNDER = 0                                                        ABSH0413
      RETURN                                                            ABSH0414
C*****SINGULAR MATRIX--USUALLY UNDERFLOW IN DET.                        ABSH0415
C*****SOMETIMES OCCURS WHEN THE SOLUTION OF THE MOMENT EQUATIONS IS     ABSH0416
C*****ABOUT TO FAIL.                                                    ABSH0417
   50 CONTINUE                                                          ABSH0418
C--      WRITE(NTTY,1010) X,Y                                           ABSH0419
C-- 1010 FORMAT(38H SINGULAR SYSTEM IN "AUX" AT X= ,F5.3,/,4H Y= ,6E10.2ABSH0420
      IUNDER = 1                                                        ABSH0421
      RETURN                                                            ABSH0422
      END                                                               ABSH0423
