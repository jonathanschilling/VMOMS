      SUBROUTINE GPASMA                                                 ABSH0538
C***********************************************************************ABSH0539
C*****GPASMA EVALUATES THE METRIC ELEMENTS FOR THE MHD EQUILIBRIUM.    *ABSH0540
C*****REFERENCES:                                                      *ABSH0541
C*****L.L.LAO,S.P.HIRSHMAN,R.M.WIELAND,ORNL/TM-7616 (1981).            *ABSH0542
C*****                                ,PHYS FLUIDS 24 (1981) AUG       *ABSH0543
C*****L.L.LAO,R.M.WIELAND,W.A.HOULBERG,S.P.HIRSHMAN,ORNL/TM-7871 (1981)*ABSH0544
C*****LAST REVISION: 6/81 L.L.LAO, R.M.WIELAND, AND W.A.HOULBERG ORNL. *ABSH0545
C*****OTHER COMMENTS:                                                  *ABSH0546
C*****MANY OF THE INPUT PARAMETERS ARE OBTAINED THROUGH COMMON BLOCKS  *ABSH0547
C*****FROM MHDEQ.                                                      *ABSH0548
C*****FOR EXPLANATION OF OTHER CALCULATED PARAMETERS SEE THE ABOVE     *ABSH0549
C*****REFERENCES AND THE OUTPUT ROUTINES.                              *ABSH0550
C***********************************************************************ABSH0551
      PARAMETER MTHETA=11,MM1=10,MNMOMS=3                               ABSH0552
      PARAMETER MNEQ=6                                                  ABSH0553
      PARAMETER MM1T=20                                                 ABSH0554
      COMMON /CD02/ XP(MM1), C(MM1,MNEQ), PARAM(MNEQ)                   ABSH0555
      COMMON /LOCGPA/ BETAT0, F0, FWALL, PSIC0, Q0, ARGEOM,             ABSH0556
     *     BETAPB, BETATB                                               ABSH0557
      COMMON /LOCMEQ/ XL0, DX, R0, A0, BT0, P0, AI0, E0, E1,            ABSH0558
     *     D0, ARC, SLPE, SLPR0, BETAP0, KMHD, IUNDER                   ABSH0559
      COMMON /SCFUNC/ S1TH(MTHETA), S2TH(MTHETA), C1TH(MTHETA),         ABSH0560
     *     C2TH(MTHETA)                                                 ABSH0561
      COMMON /SWORK2/ BSIFT(MM1), CSIFT(MM1), DSIFT(MM1),               ABSH0562
     *     BELLI(MM1), CELLI(MM1), DELLI(MM1), BTRIA(MM1),              ABSH0563
     *     CTRIA(MM1), DTRIA(MM1), R(MTHETA), RX(MTHETA),               ABSH0564
     *     RT(MTHETA), RXT(MTHETA), RTT(MTHETA), Z(MTHETA),             ABSH0565
     *     ZX(MTHETA), ZT(MTHETA), ZXT(MTHETA), ZTT(MTHETA),            ABSH0566
     *     GTT(MTHETA), TAU(MTHETA), GSQRT(MTHETA),                     ABSH0567
     *     GTTX(MTHETA), ZMOM1(MTHETA), RMOM0(MTHETA),                  ABSH0568
     *     TMOM(MTHETA), RMOM2(MTHETA), RTAU2(MTHETA),                  ABSH0569
     *     RTAU3(MTHETA), GS10(MTHETA), GS21(MTHETA),                   ABSH0570
     *     GS4(MTHETA), GS5(MTHETA), GS6(MTHETA), GS7(MTHETA),          ABSH0571
     *     GS8(MTHETA), GS2(MTHETA), GS12(MTHETA),                      ABSH0572
     *     GS22(MTHETA), CS13N(MTHETA), SS20(MTHETA),                   ABSH0573
     *     CS1N(MTHETA)                                                 ABSH0574
      COMMON /SWORK3/ PH(MM1), PPH(MM1), AIH(MM1), AIPH(MM1),           ABSH0575
     *     PSI(MM1), FPOL(MM1), SMQ(MM1), AGSQRT(MM1),                  ABSH0576
     *     AGTTGS(MM1), ATAUOR(MM1), AGRHO(MM1), ARIF(MM1),             ABSH0577
     *     AR2IF(MM1), FFP(MM1), AJORF(MM1), AJTM(MM1T),                ABSH0578
     *     RMID(MM1T), GTTGSQ(MTHETA), TAUOR(MTHETA),                   ABSH0579
     *     GRXX(MTHETA), R2I(MTHETA)                                    ABSH0580
      COMMON /WEIGHT/ WFCN(MTHETA), CNORM                               ABSH0581
      DATA CVEL /3.0E+10/, STATA /3.0E+9/, EVERG /1.6022E-12/           ABSH0582
      PI = ACOS(-1.0)                                                   ABSH0583
      M1 = MM1                                                          ABSH0584
      M1M1 = M1 - 1                                                     ABSH0585
      NTHETA = MTHETA                                                   ABSH0586
      ARGEOM = R0/A0                                                    ABSH0587
      FWALL = R0*BT0                                                    ABSH0588
      F0 = CVEL/2.0*FWALL/STATA                                         ABSH0589
      AJT0 = (F0*STATA)**2/(2.0*PI*A0**2*ARGEOM*(AI0*STATA))            ABSH0590
      BETAT0 = (P0*EVERG)*4.0*PI/BT0**2                                 ABSH0591
      Q0 = F0/AI0/ARGEOM**2                                             ABSH0592
C*****GET THE METRIC ELEMENTS AND VARIOUS AVERAGES.                     ABSH0593
      DO 10 J=1,M1                                                      ABSH0594
           CALL PRSFUN(XP(J), PH(J), PPH(J))                            ABSH0595
           CALL CURFUN(XP(J), AIH(J), AIPH(J))                          ABSH0596
   10 CONTINUE                                                          ABSH0597
      DO 30 J=1,M1                                                      ABSH0598
           R00 = C(J,1)                                                 ABSH0599
           R1 = -XP(J)                                                  ABSH0600
           R2 = C(J,5)                                                  ABSH0601
           R0X = C(J,2)                                                 ABSH0602
           R1X = -1.0                                                   ABSH0603
           R2X = C(J,6)                                                 ABSH0604
           Z1 = XP(J)*C(J,3)                                            ABSH0605
           Z2 = C(J,3)*C(J,5)                                           ABSH0606
           Z1X = C(J,3) + XP(J)*C(J,4)                                  ABSH0607
           Z2X = C(J,3)*R2X + C(J,4)*R2                                 ABSH0608
           DO 20 I=1,NTHETA                                             ABSH0609
                R(I) = R00 + R1*C1TH(I) + R2*C2TH(I)                    ABSH0610
                RX(I) = R0X + R1X*C1TH(I) + R2X*C2TH(I)                 ABSH0611
                RT(I) = -(R1*S1TH(I)+2.0*R2*S2TH(I))                    ABSH0612
                RXT(I) = -(R1X*S1TH(I)+2.0*R2X*S2TH(I))                 ABSH0613
                RTT(I) = -(R1*C1TH(I)+4.0*R2*C2TH(I))                   ABSH0614
                Z(I) = Z1*S1TH(I) + Z2*S2TH(I)                          ABSH0615
                ZX(I) = Z1X*S1TH(I) + Z2X*S2TH(I)                       ABSH0616
                ZT(I) = Z1*C1TH(I) + 2.0*Z2*C2TH(I)                     ABSH0617
                ZXT(I) = Z1X*C1TH(I) + 2.0*Z2X*C2TH(I)                  ABSH0618
                ZTT(I) = -(Z1*S1TH(I)+4.0*Z2*S2TH(I))                   ABSH0619
                TAU(I) = RT(I)*ZX(I) - RX(I)*ZT(I)                      ABSH0620
                GSQRT(I) = R(I)*TAU(I)                                  ABSH0621
                GTT(I) = RT(I)**2 + ZT(I)**2                            ABSH0622
                GTTGSQ(I) = GTT(I)/GSQRT(I)                             ABSH0623
                TAUOR(I) = TAU(I)/R(I)                                  ABSH0624
                GXT = RX(I)*RT(I) + ZX(I)*ZT(I)                         ABSH0625
                GXX = RX(I)**2 + ZX(I)**2                               ABSH0626
                GRXX(I) = GTT(I)/(GXX*GTT(I)-GXT**2)*WFCN(I)            ABSH0627
                R2I(I) = WFCN(I)/R(I)**2                                ABSH0628
   20      CONTINUE                                                     ABSH0629
           AGSQRT(J) = SDOT(NTHETA,WFCN,1,GSQRT,1)/CNORM                ABSH0630
           AGTTGS(J) = SDOT(NTHETA,WFCN,1,GTTGSQ,1)/CNORM               ABSH0631
           ATAUOR(J) = SDOT(NTHETA,WFCN,1,TAUOR,1)/CNORM                ABSH0632
           AGRHO(J) = SDOT(NTHETA,GSQRT,1,GRXX,1)/CNORM/                ABSH0633
     *          AGSQRT(J)                                               ABSH0634
           ARIF(J) = SDOT(NTHETA,TAU,1,WFCN,1)/CNORM/AGSQRT(J)          ABSH0635
           AR2IF(J) = SDOT(NTHETA,GSQRT,1,R2I,1)/CNORM/AGSQRT(J)        ABSH0636
           FFP(J) = -(BETAP0*AGSQRT(J)*PPH(J)+2.0*AIH(J)*AIPH(J)        ABSH0637
     *          /AGTTGS(J))/2.0/Q0**2/ARGEOM**4/ATAUOR(J)               ABSH0638
           AJTM(J) = -AGTTGS(J)/AIH(J)*(BETAT0/ARGEOM*PPH(J)*           ABSH0639
     *          R(1)+FFP(J)*ARGEOM/R(1))*AJT0                           ABSH0640
           AJTM(M1+J) = -AGTTGS(J)/AIH(J)*(BETAT0/ARGEOM*PPH(J)*        ABSH0641
     *          R(NTHETA)+FFP(J)*ARGEOM/R(NTHETA))*AJT0                 ABSH0642
           AJORF(J) = -AGTTGS(J)/AIH(J)*(BETAT0/ARGEOM*PPH(J)           ABSH0643
     *          +FFP(J)*ARGEOM*AR2IF(J))*AJT0/ARIF(J)                   ABSH0644
           RMID(J) = R(1)                                               ABSH0645
           RMID(M1+J) = R(NTHETA)                                       ABSH0646
   30 CONTINUE                                                          ABSH0647
C*****PERFORM INTEGRATIONS.                                             ABSH0648
      SUMF2 = 0.0                                                       ABSH0649
C     SUMG = AGSQRT(1)/2.0                                              ABSH0650
C     SUMPSI = AIH(1)/AGTTGS(1)/2.0                                     ABSH0651
C     SUMP = PH(1)*AGSQRT(1)/2.0                                        ABSH0652
      SUMG = AGSQRT(1)/2.0*XL0/DX                                       ABSH004A
      SUMPSI = AIH(1)/AGTTGS(1)/2.0*XL0/DX                              ABSH005A
      SUMP = PH(1)*AGSQRT(1)/2.0*XL0/DX                                 ABSH006A
      PSI(1) = SUMPSI*DX                                                ABSH0653
      DO 40 I=1,M1M1                                                    ABSH0654
           SUMG = SUMG + (AGSQRT(I)+AGSQRT(I+1))/2.0                    ABSH0655
           SUMF2 = SUMF2 + (FFP(M1-I)+FFP(M1-I+1))/2.0                  ABSH0656
           SUMP = SUMP + (PH(I)*AGSQRT(I)+PH(I+1)*AGSQRT(I+1))/         ABSH0657
     *          2.0                                                     ABSH0658
           SUMPSI = SUMPSI + (AIH(I)/AGTTGS(I)+AIH(I+1)/                ABSH0659
     *          AGTTGS(I+1))/2.0                                        ABSH0660
           FPOL(M1-I) = SQRT(1.0-SUMF2*2.0*DX)                          ABSH0661
           SMQ(M1-I) = FPOL(M1-I)/AIH(M1-I)*Q0*ARGEOM**2*               ABSH0662
     *          AGTTGS(M1-I)*ATAUOR(M1-I)                               ABSH0663
           PSI(I+1) = SUMPSI*DX                                         ABSH0664
   40 CONTINUE                                                          ABSH0665
      PSIC0 = 4.0*PI*A0*(AI0*STATA)/CVEL*SUMPSI*DX                      ABSH0666
      FPOL(M1) = 1.0                                                    ABSH0667
      SMQ(M1) = Q0*ARGEOM**2*AGTTGS(M1)*ATAUOR(M1)                      ABSH0668
      BETAPB = SUMP/SUMG*BETAP0*AGSQRT(M1)*AGTTGS(M1)                   ABSH0669
      BETATB = SUMP/SUMG*PI*8.0/BT0**2*(P0*EVERG)                       ABSH0670
      DO 50 I=1,M1                                                      ABSH0671
           PSI(I) = PSI(I)/PSI(M1)                                      ABSH0672
   50 CONTINUE                                                          ABSH0673
      RETURN                                                            ABSH0674
      END                                                               ABSH0675
