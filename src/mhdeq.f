      SUBROUTINE MHDEQ(KMHDEQ, IPARAM, NRHO, XR0, XBT0, XP0,            ABSH0676
     *     XAI0, XE0, XE1, XD0, RHO, SHIFT, ELONG, TRIANG,              ABSH0677
     *     IFAIL)                                                       ABSH0678
C***********************************************************************ABSH0679
C*****MHDEQ SOLVES FOR THE MOMENTS OF THE GRAD-SHAFRANOV EQUATION.     *ABSH0680
C*****REFERENCES:                                                      *ABSH0681
C*****L.L.LAO,S.P.HIRSHMAN,R.M.WIELAND,ORNL/TM-7616 (1981).            *ABSH0682
C*****                                ,PHYS FLUIDS 24 (1981) AUG       *ABSH0683
C*****L.L.LAO,R.M.WIELAND,W.A.HOULBERG,S.P.HIRSHMAN,ORNL/TM-7871 (1981)*ABSH0684
C*****LAST REVISION: 6/81 L.L.LAO, R.M.WIELAND, AND W.A.HOULBERG ORNL. *ABSH0685
C*****CALCULATED PARAMETERS:                                           *ABSH0686
C*****SHIFT(I)-SHIFT OF SURFACE I RELATIVE TO MAG AXIS-(CM).           *ABSH0687
C*****ELONG(I)-ELONGATION OF SURFACE I-(DIMENSIONLESS).                *ABSH0688
C*****TRIANG(I)-TRIANGULARITY OF SURFACE I-(DIMENSIONLESS).            *ABSH0689
C*****INPUT PARAMETERS:                                                *ABSH0690
C*****KMHDEQ-NUMBER OF EQUATIONS TO BE SOLVED.                         *ABSH0691
C*****      =1 SOLVE SHIFT, APPROXIMATE ELONGATION AND TRIANGULARITY.  *ABSH0692
C*****      =2 SOLVE SHIFT AND ELONGATION, APPROXIMATE TRIANGULARITY.  *ABSH0693
C*****      =3 SOLVE SHIFT, ELONGATION, AND TRIANGULARITY.             *ABSH0694
C*****IPARAM-SWITCH FOR INITIAL GUESS AT SOLUTION.                     *ABSH0695
C*****      =0 CONSTRUCT GUESS FROM XE0, XE1, AND XD0.                 *ABSH0696
C*****      >0 USE SOLUTION FROM PREVIOUS CALL IF AVAILABLE.           *ABSH0697
C*****NRHO-NUMBER OF RADIAL GRID POINTS.                               *ABSH0698
C*****R0-MAJOR RADIUS TO GEOMETRIC CENTER OF OUTERMOST SURFACE-(CM).   *ABSH0699
C*****BT0-VACUUM TOROIDAL FIELD AT R0-(GAUSS).                         *ABSH0700
C*****P0-PLASMA PRESSURE AT MAGNETIC AXIS-(EV/CM3).                    *ABSH0701
C*****AI0-TOTAL TOROIDAL CURRENT-(AMPS).                               *ABSH0702
C*****XE0-ELONGATION OF OUTERMOST FLUX SURFACE-(DIMENSIONLESS).        *ABSH0703
C*****XE1-ELONGATION ON AXIS FOR APPROXIMATE SOLUTION-(DIMENSIONLESS). *ABSH0704
C*****XD0-TRIANGULARITY OF OUTERMOST FLUX SURFACE-(DIMENSIONLESS).     *ABSH0705
C*****RHO(I)-HALF-DIAMETER OF SURFACE I IN MIDPLANE-(CM).              *ABSH0706
C*****OTHER COMMENTS:                                                  *ABSH0707
C*****LET RG(I) BE THE MAJOR RADIUS TO THE GEOMETRIC CENTER OF FLUX    *ABSH0708
C*****SURFACE I                                                        *ABSH0709
C*****THEN SHIFT(I)=RMAG-R0.                                           *ABSH0710
C*****LET ZMAX(I) BE THE MAXIMUM HEIGHT OF FLUX SURFACE I, R(ZMAX(I))  *ABSH0711
C*****THE DISTANCE FROM THIS POINT TO THE MAJOR AXIS THEN              *ABSH0712
C*****ELONG(I)=ZMAX(I)/RHO(I).                                         *ABSH0713
C*****TRIANG(I)=(RG(I)-R(ZMAX(I)))/RHO(I).                             *ABSH0714
C***********************************************************************ABSH0715
      PARAMETER MNMX=61                                                 ABSH0716
      PARAMETER MTHETA=11,MM1=10,MNMOMS=3                               ABSH0717
      PARAMETER MNEQ=6                                                  ABSH0718
      COMMON /CD02/ XP(MM1), C(MM1,MNEQ), PARAM(MNEQ)                   ABSH0719
      COMMON /CINIT/ INITM, INITC, INITP, ITPR                          ABSH0720
      COMMON /LOCD02/ ERROR(MNEQ), PARERR(MNEQ),                        ABSH0721
     *     MAT(MNEQ,MNEQ), WSPACE(MNEQ,9), WSPAC1(MNEQ),                ABSH0722
     *     WSPAC2(MNEQ), COPY(MNEQ,MNEQ), HH, HH0                       ABSH0723
      COMMON /LOCMEQ/ XL0, DX, R0, A0, BT0, P0, AI0, E0, E1,            ABSH0724
     *     D0, ARC, SLPE, SLPR0, BETAP0, KMHD, IUNDER                   ABSH0725
      COMMON /SCFUNC/ S1TH(MTHETA), S2TH(MTHETA), C1TH(MTHETA),         ABSH0726
     *     C2TH(MTHETA)                                                 ABSH0727
      COMMON /SWORK2/ BSIFT(MM1), CSIFT(MM1), DSIFT(MM1),               ABSH0728
     *     BELLI(MM1), CELLI(MM1), DELLI(MM1), BTRIA(MM1),              ABSH0729
     *     CTRIA(MM1), DTRIA(MM1), R(MTHETA), RX(MTHETA),               ABSH0730
     *     RT(MTHETA), RXT(MTHETA), RTT(MTHETA), Z(MTHETA),             ABSH0731
     *     ZX(MTHETA), ZT(MTHETA), ZXT(MTHETA), ZTT(MTHETA),            ABSH0732
     *     GTT(MTHETA), TAU(MTHETA), GSQRT(MTHETA),                     ABSH0733
     *     GTTX(MTHETA), ZMOM1(MTHETA), RMOM0(MTHETA),                  ABSH0734
     *     TMOM(MTHETA), RMOM2(MTHETA), RTAU2(MTHETA),                  ABSH0735
     *     RTAU3(MTHETA), GS10(MTHETA), GS21(MTHETA),                   ABSH0736
     *     GS4(MTHETA), GS5(MTHETA), GS6(MTHETA), GS7(MTHETA),          ABSH0737
     *     GS8(MTHETA), GS2(MTHETA), GS12(MTHETA),                      ABSH0738
     *     GS22(MTHETA), CS13N(MTHETA), SS20(MTHETA),                   ABSH0739
     *     CS1N(MTHETA)                                                 ABSH0740
      COMMON /WEIGHT/ WFCN(MTHETA), CNORM                               ABSH0741
      DIMENSION RHO(NRHO), SHIFT(NRHO), ELONG(NRHO),                    ABSH0742
     *     TRIANG(NRHO)                                                 ABSH0743
      EXTERNAL AUX, BCAUX, RAAUX, PRSOL                                 ABSH0744
      DATA ERRY /1.E-02/, ERRP /1.E-03/, SLPR0 /1.0E-4/, SLPE           ABSH0745
     *     /1.0E-4/                                                     ABSH0746
      DATA XL00 /0.03/, HH0 /0.01/                                      ABSH0747
      DATA CVEL /3.0E+10/, STATA /3.0E+9/, EVERG /1.6022E-12/           ABSH0748
C*****INITIALIZATION.                                                   ABSH0749
      IF (INITM.EQ.1) GO TO 30                                          ABSH0750
      M1 = MM1                                                          ABSH0751
      MTM2 = MTHETA - 2                                                 ABSH0752
      NEQ = MNEQ                                                        ABSH0753
C*****SET UP SINE AND COSINE ARRAYS.                                    ABSH0754
      CALL SINCOS                                                       ABSH0755
C*****SET UP THE INTEGRATION WEIGHTING ARRAY.                           ABSH0756
      WFCN(1) = 1.0                                                     ABSH0757
      DO 10 I=2,MTM2,2                                                  ABSH0758
           WFCN(I) = 4.0                                                ABSH0759
           WFCN(I+1) = 2.0                                              ABSH0760
   10 CONTINUE                                                          ABSH0761
      WFCN(MTHETA-1) = 4.0                                              ABSH0762
      WFCN(MTHETA) = 1.0                                                ABSH0763
      CNORM = 3.0*(MTHETA-1)                                            ABSH0764
      DO 20 I=1,6                                                       ABSH0765
           ERROR(I) = ERRY                                              ABSH0766
           PARERR(I) = ERRP                                             ABSH0767
   20 CONTINUE                                                          ABSH0768
      INITM = 1                                                         ABSH0769
   30 KMHD = KMHDEQ                                                     ABSH0770
      R0 = XR0                                                          ABSH0771
      BT0 = XBT0                                                        ABSH0772
      P0 = XP0                                                          ABSH0773
      AI0 = XAI0                                                        ABSH0774
      HH = HH0                                                          ABSH0775
      INITC = 0                                                         ABSH0776
      INITP = 0                                                         ABSH0777
      ITPR = 0                                                          ABSH0778
C*****SET UP DIMENSIONLESS RADIAL GRID.                                 ABSH0779
      A0 = RHO(NRHO)                                                    ABSH0780
      XL0 = XL00                                                        ABSH0781
      XL0 = AMIN1(XL0,RHO(2)/A0)                                        ABSH0782
      DX = (1.0-XL0)/(M1-1)                                             ABSH0783
      DO 40 I=1,M1                                                      ABSH0784
           XP(I) = XL0 + (I-1)*DX                                       ABSH0785
   40 CONTINUE                                                          ABSH0786
      PI = ACOS(-1.0)                                                   ABSH0787
      BETAP0 = 2.0*PI*CVEL**2*(P0*EVERG)*A0**2/(AI0*STATA)**2           ABSH0788
C*****CONVERT TO MOMENTS REPRESENTATION FOR ELONGATION AND TRIANGULARITYABSH0789
      CALL GEOTRN(1, 1.0, XE0, XD0, E0, D0)                             ABSH0790
      E1 = XE1                                                          ABSH0791
      ARC = (R0-D0*A0)/A0                                               ABSH0792
      IF (IPARAM.NE.0 .AND. PARAM(3).NE.0.0) GO TO 50                   ABSH0793
C*****SET UP THE INITIAL GUESS FOR THE PARAMETERS.                      ABSH0794
      PARAM(1) = ARC + 0.1                                              ABSH0795
      PARAM(2) = -0.4                                                   ABSH0796
      PARAM(3) = E1                                                     ABSH0797
      IF (KMHD.GT.1) PARAM(3) = 1.0 + (E0-1.0)/2.0                      ABSH0798
      PARAM(4) = 0.0                                                    ABSH0799
      PARAM(5) = 1.0                                                    ABSH0800
      PARAM(6) = 2.0*D0                                                 ABSH0801
C*****CALL THE ODE SOLVER.                                              ABSH0802
   50 IFAIL = 1                                                         ABSH0803
      CALL D02AGF(HH, ERROR, PARERR, PARAM, C, NEQ, 2*KMHD, M1,         ABSH0804
     *     AUX, BCAUX, RAAUX, PRSOL, MAT, COPY, WSPACE, WSPAC1,         ABSH0805
     *     WSPAC2, IFAIL)                                               ABSH0806
      IF (IUNDER.EQ.1) IFAIL = 8                                        ABSH0807
      IF (IFAIL.GT.0) RETURN                                            ABSH0808
C*****INTERPOLATE ONTO THE DRIVER GRID.                                 ABSH0809
      CALL SPLAAN(M1, XP, C(1,1), BSIFT, CSIFT, DSIFT)                  ABSH0810
      CALL SPLAAN(M1, XP, C(1,3), BELLI, CELLI, DELLI)                  ABSH0811
      CALL SPLAAN(M1, XP, C(1,5), BTRIA, CTRIA, DTRIA)                  ABSH0812
      DO 60 I=2,NRHO                                                    ABSH0813
           X = RHO(I)/A0                                                ABSH0814
           SHIFT(I) = SEVAL(M1,X,XP,C(1,1),BSIFT,CSIFT,DSIFT)           ABSH0815
           EMOM = SEVAL(M1,X,XP,C(1,3),BELLI,CELLI,DELLI)               ABSH0816
           TRMOM = SEVAL(M1,X,XP,C(1,5),BTRIA,CTRIA,DTRIA)              ABSH0817
           SHIFT(I) = (SHIFT(I)+TRMOM)*A0 - R0                          ABSH0818
           CALL GEOTRN(2, X, EMOM, TRMOM, ELONG(I), TRIANG(I))          ABSH0819
   60 CONTINUE                                                          ABSH0820
      SHIFT(1) = PARAM(1)*A0 - R0                                       ABSH0821
      ELONG(1) = PARAM(3)                                               ABSH0822
      TRIANG(1) = 0.0                                                   ABSH0823
      RETURN                                                            ABSH0824
      END                                                               ABSH0825
