C***********************************************************************ABSH0001
C*****MDRIV IS THE MAIN DRIVING FOR SOLVING THE MOMENT EQUATIONS OF THE*ABSH0002
C*****GRAD-SHAFRANOV EQUATION USING A VARIATIONAL METHOD.              *ABSH0003
C***********************************************************************ABSH0004
      PROGRAM VMOMS                                                     ABSH0005
      PARAMETER MTHETA=11,MM1=10,MNMOMS=3                               ABSH0006
      PARAMETER MNEQ=6,MNMX=61                                          ABSH0007
      COMMON /NIO/ NIN, NOUT, NPLOT, NDEBUG, NTTY                       ABSH0008
      COMMON /SWORK1/ PRES(MNMX), CURT(MNMX), RHO(MNMX),                ABSH0009
     *     BPRES(MNMX), CPRES(MNMX), DPRES(MNMX), BCURT(MNMX),          ABSH0010
     *     CCURT(MNMX), DCURT(MNMX), NRHO                               ABSH0011
      COMMON /VINPUT/ PA, PB, AIA, AIB, AIC, XE0, XE1, XD0              ABSH0012
      COMMON /VOUPUT/ SHIFT(MNMX), ELONG(MNMX), TRIANG(MNMX)            ABSH0013
      NAMELIST/MOMIN/KMHD,IPARAM,R0,A0,BT0,P0,PA,PB,AI0,AIA,XE0,        ABSH0014
     *XE1,XD0                                                           ABSH0015
      DATA IPARAM /0/, PA /2.0/, PB /1.0/, AIA /2.0/                    ABSH0016
      DATA R0 /94./, A0 /26./, BT0 /9100./, P0 /1.32E+17/, AI0          ABSH0017
     *     /2.E+5/                                                      ABSH0018
      DATA XE0 /1.55/, XE1 /1.0/, XD0 /0.18/                            ABSH0019
      NIN = 10                                                          ABSH0020
      NOUT = 11                                                         ABSH0021
      NPLOT = 12                                                        ABSH0022
      NTTY = 5                                                          ABSH0023
      NDEBUG = 13                                                       ABSH0024
      READ (NIN,MOMIN)                                                  ABSH0025
      NRHO = 25                                                         ABSH0026
      DXX = A0/FLOAT(NRHO-1)                                            ABSH0027
      DO 10 I=1,NRHO                                                    ABSH0028
           RHO(I) = (I-1)*DXX                                           ABSH0029
           X = RHO(I)/A0                                                ABSH0030
           PRES(I) = P0*(1.0-X**PA)**PB                                 ABSH0031
           XSQ = X**2                                                   ABSH0032
           AIB = 3.0 - 2.0*AIA                                          ABSH0033
           AIC = AIA - 2.0                                              ABSH0034
           CURT(I) = AI0*XSQ*(AIA+XSQ*(AIB+XSQ*AIC))                    ABSH0035
   10 CONTINUE                                                          ABSH0036
      IFAIL = 1                                                         ABSH0037
      CALL MHDEQ(KMHD, IPARAM, NRHO, R0, BT0, P0, AI0, XE0,             ABSH0038
     *     XE1, XD0, RHO, SHIFT, ELONG, TRIANG, IFAIL)                  ABSH0039
      CALL GPASMA                                                       ABSH0040
      CALL PRTMEQ                                                       ABSH0041
      CALL PLTMEQ                                                       ABSH0042
      STOP                                                              ABSH0043
      END                                                               ABSH0044
      SUBROUTINE PLTMEQ                                                 ABSH0045
C***********************************************************************ABSH0046
C*****PLTMEQ GENERATES A DATA FILE FOR THE EQUILIBRIUM QUANTITIES.     *ABSH0047
C***********************************************************************ABSH0048
      PARAMETER MNMX=61                                                 ABSH0049
      PARAMETER MTHETA=11,MM1=10,MNMOMS=3                               ABSH0050
      PARAMETER MNEQ=6                                                  ABSH0051
      PARAMETER MM1T=20                                                 ABSH0052
      COMMON /CD02/ XP(MM1), C(MM1,MNEQ), PARAM(MNEQ)                   ABSH0053
      COMMON /LOCGPA/ BETAT0, F0, FWALL, PSIC0, Q0, ARGEOM,             ABSH0054
     *     BETAPB, BETATB                                               ABSH0055
      COMMON /LOCMEQ/ XL0, DX, R0, A0, BT0, P0, AI0, E0, E1,            ABSH0056
     *     D0, ARC, SLPE, SLPR0, BETAP0, KMHD, IUNDER                   ABSH0057
      COMMON /NIO/ NIN, NOUT, NPLOT, NDEBUG, NTTY                       ABSH0058
      COMMON /SWORK1/ PRES(MNMX), CURT(MNMX), RHO(MNMX),                ABSH0059
     *     BPRES(MNMX), CPRES(MNMX), DPRES(MNMX), BCURT(MNMX),          ABSH0060
     *     CCURT(MNMX), DCURT(MNMX), NRHO                               ABSH0061
      COMMON /SWORK3/ PH(MM1), PPH(MM1), AIH(MM1), AIPH(MM1),           ABSH0062
     *     PSI(MM1), FPOL(MM1), SMQ(MM1), AGSQRT(MM1),                  ABSH0063
     *     AGTTGS(MM1), ATAUOR(MM1), AGRHO(MM1), ARIF(MM1),             ABSH0064
     *     AR2IF(MM1), FFP(MM1), AJORF(MM1), AJTM(MM1T),                ABSH0065
     *     RMID(MM1T), GTTGSQ(MTHETA), TAUOR(MTHETA),                   ABSH0066
     *     GRXX(MTHETA), R2I(MTHETA)                                    ABSH0067
      COMMON /VOUPUT/ SHIFT(MNMX), ELONG(MNMX), TRIANG(MNMX)            ABSH0068
      M1 = MM1                                                          ABSH0069
      WRITE (NPLOT,99999) M1, NRHO                                      ABSH0070
      WRITE (NPLOT,99998) ARGEOM, A0, D0, PSIC0, PARAM                  ABSH0071
      DO 10 I=1,MM1                                                     ABSH0072
           WRITE (NPLOT,99998) XP(I), C(I,1), C(I,3), C(I,5)            ABSH0073
           WRITE (NPLOT,99998) PSI(I), AGSQRT(I), SMQ(I),               ABSH0074
     *          AJORF(I), FFP(I), AJTM(I), AJTM(MM1+I),                 ABSH0075
     *          FPOL(I), RMID(I), RMID(MM1+I), AIH(I)                   ABSH0076
           WRITE (NPLOT,99998) PH(I), PPH(I)                            ABSH0077
   10 CONTINUE                                                          ABSH0078
      DO 20 I=1,NRHO                                                    ABSH0079
           WRITE (NPLOT,99998) RHO(I), SHIFT(I), ELONG(I),              ABSH0080
     *          TRIANG(I)                                               ABSH0081
   20 CONTINUE                                                          ABSH0082
      RETURN                                                            ABSH0083
99999 FORMAT (1X, 2I5)                                                  ABSH0084
99998 FORMAT (1X, 10E12.4)                                              ABSH0085
      END                                                               ABSH0086
      SUBROUTINE PRTMEQ                                                 ABSH0087
C***********************************************************************ABSH0088
C*****PRTMEQ WRITES AN OUTPUT FILE SUMMARIZING THE RESULTS.            *ABSH0089
C***********************************************************************ABSH0090
      PARAMETER MNMX=61                                                 ABSH0091
      PARAMETER MTHETA=11,MM1=10,MNMOMS=3                               ABSH0092
      PARAMETER MNEQ=6                                                  ABSH0093
      PARAMETER MM1T=20                                                 ABSH0094
      COMMON /CD02/ XP(MM1), C(MM1,MNEQ), PARAM(MNEQ)                   ABSH0095
      COMMON /LOCD02/ ERROR(MNEQ), PARERR(MNEQ),                        ABSH0096
     *     MAT(MNEQ,MNEQ), WSPACE(MNEQ,9), WSPAC1(MNEQ),                ABSH0097
     *     WSPAC2(MNEQ), COPY(MNEQ,MNEQ), HH, HH0                       ABSH0098
      COMMON /LOCGPA/ BETAT0, F0, FWALL, PSIC0, Q0, ARGEOM,             ABSH0099
     *     BETAPB, BETATB                                               ABSH0100
      COMMON /LOCMEQ/ XL0, DX, R0, A0, BT0, P0, AI0, E0, E1,            ABSH0101
     *     D0, ARC, SLPE, SLPR0, BETAP0, KMHD, IUNDER                   ABSH0102
      COMMON /NIO/ NIN, NOUT, NPLOT, NDEBUG, NTTY                       ABSH0103
      COMMON /SWORK1/ PRES(MNMX), CURT(MNMX), RHO(MNMX),                ABSH0104
     *     BPRES(MNMX), CPRES(MNMX), DPRES(MNMX), BCURT(MNMX),          ABSH0105
     *     CCURT(MNMX), DCURT(MNMX), NRHO                               ABSH0106
      COMMON /SWORK3/ PH(MM1), PPH(MM1), AIH(MM1), AIPH(MM1),           ABSH0107
     *     PSI(MM1), FPOL(MM1), SMQ(MM1), AGSQRT(MM1),                  ABSH0108
     *     AGTTGS(MM1), ATAUOR(MM1), AGRHO(MM1), ARIF(MM1),             ABSH0109
     *     AR2IF(MM1), FFP(MM1), AJORF(MM1), AJTM(MM1T),                ABSH0110
     *     RMID(MM1T), GTTGSQ(MTHETA), TAUOR(MTHETA),                   ABSH0111
     *     GRXX(MTHETA), R2I(MTHETA)                                    ABSH0112
      COMMON /VINPUT/ PA, PB, AIA, AIB, AIC, XE0, XE1, XD0              ABSH0113
      COMMON /VOUPUT/ SHIFT(MNMX), ELONG(MNMX), TRIANG(MNMX)            ABSH0114
      DIMENSION Y(MNEQ), YP(MNEQ), CPP(MM1,MNMOMS)                      ABSH0115
      PI = ACOS(-1.0)                                                   ABSH0116
      WRITE (NOUT,99999)                                                ABSH0117
      WRITE (NOUT,99998) R0, A0, XE0, XE1, XD0, BT0, P0, PA, PB         ABSH0118
      WRITE (NOUT,99997) AI0, AIA                                       ABSH0119
      WRITE (NOUT,99996) KMHD, HH0, XL0, PARERR, ERROR                  ABSH0120
      WRITE (NOUT,99995) BETAP0, F0, Q0, PSIC0, E0, E1, D0              ABSH0121
      DO 20 I=1,MM1                                                     ABSH0122
           DO 10 J=1,MNEQ                                               ABSH0123
                Y(J) = C(I,J)                                           ABSH0124
   10      CONTINUE                                                     ABSH0125
           X = XP(I)                                                    ABSH0126
           CALL AUX(YP, Y, X, PARAM)                                    ABSH0127
           CPP(I,1) = YP(2)                                             ABSH0128
           CPP(I,2) = YP(4)                                             ABSH0129
           CPP(I,3) = YP(6)                                             ABSH0130
   20 CONTINUE                                                          ABSH0131
      WRITE (NOUT,99994)                                                ABSH0132
      DO 30 I=1,MM1                                                     ABSH0133
           WRITE (NOUT,99993) I, XP(I), PH(I), PPH(I), AIH(I),          ABSH0134
     *          AIPH(I)                                                 ABSH0135
   30 CONTINUE                                                          ABSH0136
      WRITE (NOUT,99992)                                                ABSH0137
      DO 40 I=1,MM1                                                     ABSH0138
           WRITE (NOUT,99991) I, XP(I), C(I,1), C(I,2),                 ABSH0139
     *          CPP(I,1), C(I,3), C(I,4), CPP(I,2), C(I,5),             ABSH0140
     *          C(I,6), CPP(I,3)                                        ABSH0141
   40 CONTINUE                                                          ABSH0142
      WRITE (NOUT,99990)                                                ABSH0143
      DO 50 I=1,MM1                                                     ABSH0144
           WRITE (NOUT,99989) I, XP(I), SMQ(I), AJORF(I),               ABSH0145
     *          PSI(I), FPOL(I), FFP(I), AJTM(I), AJTM(MM1+I)           ABSH0146
   50 CONTINUE                                                          ABSH0147
      WRITE (NOUT,99988)                                                ABSH0148
      DO 60 I=1,MM1                                                     ABSH0149
           WRITE (NOUT,99989) I, XP(I), AGSQRT(I), AGTTGS(I),           ABSH0150
     *          ATAUOR(I), AGRHO(I), ARIF(I), AR2IF(I)                  ABSH0151
   60 CONTINUE                                                          ABSH0152
      WRITE (NOUT,99987)                                                ABSH0153
      WRITE (NOUT,99986)                                                ABSH0154
      DO 70 I=1,NRHO                                                    ABSH0155
           WRITE (NOUT,99985) I, RHO(I), PRES(I), CURT(I),              ABSH0156
     *          SHIFT(I), ELONG(I), TRIANG(I)                           ABSH0157
   70 CONTINUE                                                          ABSH0158
      RETURN                                                            ABSH0159
99999 FORMAT (1H1, //, 1X, 40(1H*), 24H  VMOMS  SAMPLE  OUTPUT ,        ABSH0160
     *     1H , 40(1H*), ///)                                           ABSH0161
99998 FORMAT (25X, 20(1H*), 20H PLASMA  PARAMETERS , 25(1H*),           ABSH0162
     *     //, 40X, 5HR0 = , F5.1, T60, 2HCM, /, 40X, 5HA0 = ,          ABSH0163
     *     F5.1, T60, 2HCM, /, 40X, 6HXE0 = , F5.2, T60,                ABSH0164
     *     9HGEOMETRIC, /, 40X, 6HXE1 = , F5.2, T60, 8HGEOMETRI,        ABSH0165
     *     1HC, /, 40X, 5HXD0 =, F5.2, T60, 9HGEOMETRIC, /,             ABSH0166
     *     40X, 6HBT0 = , F7.1, T60, 5HGAUSS, /, 40X, 5HP0 = ,          ABSH0167
     *     1PE9.2, T60, 8HEV/CM**3, /, 40X, 5HPA = , 0PF4.2, /,         ABSH0168
     *     40X, 5HPB = , F4.2)                                          ABSH0169
99997 FORMAT (40X, 8HAI0AM = , 1PE8.2, T60, 4HAMPS, /, 40X,             ABSH0170
     *     6HAIA = , 0PF4.2, ///)                                       ABSH0171
99996 FORMAT (25X, 20(1H*), 20H D02AGF  PARAMETERS , 20(1H*),           ABSH0172
     *     //, 40X, 9HKMHD   = , I4, /, 40X, 6HHH0 = , F7.5, /,         ABSH0173
     *     40X, 6HXL0 = , F7.5, /, 40X, 9HPARERR = , 1P6E10.2,          ABSH0174
     *     /, 40X, 9HERROR  = , 1P6E10.2, ///)                          ABSH0175
99995 FORMAT (25X, 20(1H*), 22H INTERNAL  PARAMETERS , 20(1H*),         ABSH0176
     *     //, 40X, 9HBETAP0 = , F6.2, /, 40X, 5HF0 = ,                 ABSH0177
     *     1PE10.2, 5H AMPS, /, 40X, 5HQ0 = , 1PE10.2, /, 40X,          ABSH0178
     *     8HPSIC0 = , 1PE10.2, 12H GAUSS-CM**2, /, 40X,                ABSH0179
     *     5HE0 = , 0PF4.2, T60, 7HMOMENTS, /, 40X, 5HE1 = ,            ABSH0180
     *     F4.2, T60, 7HMOMENTS, /, 40X, 5HD0 = , F5.3, T60,            ABSH0181
     *     7HMOMENTS)                                                   ABSH0182
99994 FORMAT (1H1, 25X, 20(1H*), 19H MOMENTS(INTERNAL) ,                ABSH0183
     *     20(1H*), //, 3H  I, 3X, 1HX, T15, 5HP-HAT, T27,              ABSH0184
     *     6HPP-HAT, T39, 5HI-HAT, T50, 6HIP-HAT)                       ABSH0185
99993 FORMAT (I3, F6.2, 1P4E12.3)                                       ABSH0186
99992 FORMAT (//, 3H  I, 3X, 1HX, T17, 2HR0, T29, 3HR0P, T39,           ABSH0187
     *     4HR0PP, T53, 2HE0, T65, 3HE0P, T76, 4HE0PP, T89,             ABSH0188
     *     2HD0, T101, 3HD0P, T112, 4HD0PP)                             ABSH0189
99991 FORMAT (I3, F6.2, 1P9E12.3)                                       ABSH0190
99990 FORMAT (//, 3H  I, 3X, 1HX, T16, 1HQ, T27, 5HJFLUX, T39,          ABSH0191
     *     3HPSI, T51, 1HF, T62, 3HFFP, T72, 9HJMDPL(IN), T84,          ABSH0192
     *     10HJMDPL(OUT), /, T25, 10H(STATAMPS), T72, 7H(STATAM,        ABSH0193
     *     3HPS), T84, 10H(STATAMPS))                                   ABSH0194
99989 FORMAT (I3, F6.2, 1P8E12.3)                                       ABSH0195
99988 FORMAT (//, 3H  I, 3X, 1HX, T14, 6H GSQRT, T26, 6H GTTGS,         ABSH0196
     *     T38, 6H GSGPP, T51, 5H GRHO, T62, 4HARIF, T74,               ABSH0197
     *     5HAR2IF)                                                     ABSH0198
99987 FORMAT (1H1, 25X, 20(1H*), 19H MOMENTS(EXTERNAL) ,                ABSH0199
     *     20(1H*), ////)                                               ABSH0200
99986 FORMAT (//, 3H  I, 3X, 1HR, T17, 1HP, T29, 2HAI, T39,             ABSH0201
     *     5HSHIFT, T51, 5HELONG, T63, 6HTRIANG, /, T6, 4H(CM),         ABSH0202
     *     T13, 10H(EV/CM**3), T27, 6H(AMPS), T39, 4H(CM), /)           ABSH0203
99985 FORMAT (I3, F6.2, 1P5E12.3)                                       ABSH0204
      END                                                               ABSH0205
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
      SUBROUTINE BCAUX(GL, GR, PARAM)                                   ABSH0424
C***********************************************************************ABSH0425
C*****BCAUX SETS UP THE BOUNDARY CONDITIONS FOR THE MHD MOMENTS.       *ABSH0426
C*****REFERENCES:                                                      *ABSH0427
C*****L.L.LAO,S.P.HIRSHMAN,R.M.WIELAND,ORNL/TM-7616 (1981).            *ABSH0428
C*****                                ,PHYS FLUIDS 24 (1981) AUG       *ABSH0429
C*****L.L.LAO,R.M.WIELAND,W.A.HOULBERG,S.P.HIRSHMAN,ORNL/TM-7871 (1981)*ABSH0430
C*****LAST REVISION: 6/81 L.L.LAO, R.M.WIELAND, AND W.A.HOULBERG ORNL. *ABSH0431
C*****CALCULATED PARAMETERS:                                           *ABSH0432
C*****GL(I)-VALUE OF Y(I) AT LEFT BOUNDARY.                            *ABSH0433
C*****GR(I)-VALUE OF Y(I) AT RIGHT BOUNDARY.                           *ABSH0434
C*****INPUT PARAMETERS:                                                *ABSH0435
C*****PARAM(I)-B.C. PARAMETER FOR FIRST ORDER ODE I.                   *ABSH0436
C*****XL0-LOWER BOUND OF ABSCISSA.                                     *ABSH0437
C*****SLPR0-ANALYTIC PARAMETER FOR LEFT B.C. AT X1.                    *ABSH0438
C*****SLPE-ANALYTIC PARAMETER FOR LEFT B.C. AT X1.                     *ABSH0439
C*****ARC-INVERSE ASPECT RATIO FOR OUTERMOST FLUX SURFACE.             *ABSH0440
C*****E0-MOMENTS REPRESENTATION OF ELONGATION AT X=1.                  *ABSH0441
C*****D0-MOMENTS REPRESENTATION OF TRIANGULARITY AT X=1.               *ABSH0442
C*****OTHER COMMENTS:                                                  *ABSH0443
C*****BCAUX IS REQUIRED BY D02AGF ODE SOLVER IN NAG LIBRARY.           *ABSH0444
C***********************************************************************ABSH0445
      PARAMETER MNMOMS=3,MNEQ=6                                         ABSH0446
      COMMON /LOCMEQ/ XL0, DX, R0, A0, BT0, P0, AI0, E0, E1,            ABSH0447
     *     D0, ARC, SLPE, SLPR0, BETAP0, KMHD, IUNDER                   ABSH0448
      DIMENSION GL(MNEQ), GR(MNEQ), PARAM(MNEQ)                         ABSH0449
C*****LEFT B.C.'S - CENTER.                                             ABSH0450
      GL(1) = PARAM(1) - SLPR0*XL0**2                                   ABSH0451
      GL(2) = -2.0*SLPR0*XL0                                            ABSH0452
      GL(3) = PARAM(3) + SLPE*XL0**2                                    ABSH0453
      GL(4) = 2.0*XL0*SLPE                                              ABSH0454
      GL(5) = 10.0*PARAM(5)*XL0**4                                      ABSH0455
      GL(6) = 40.0*PARAM(5)*XL0**3                                      ABSH0456
C*****RIGHT B.C.'S - EDGE.                                              ABSH0457
      GR(1) = ARC                                                       ABSH0458
      GR(2) = PARAM(2)                                                  ABSH0459
      GR(3) = E0                                                        ABSH0460
      GR(4) = PARAM(4)                                                  ABSH0461
      GR(5) = D0                                                        ABSH0462
      GR(6) = PARAM(6)                                                  ABSH0463
      RETURN                                                            ABSH0464
      END                                                               ABSH0465
      SUBROUTINE CURFUN(X, AI, AIP)                                     ABSH0466
C***********************************************************************ABSH0467
C*****CURFUN PROVIDES THE INTEGRATED TOROIDAL CURRENT PROFILE.         *ABSH0468
C*****REFERENCES:                                                      *ABSH0469
C*****L.L.LAO,R.M.WIELAND,W.A.HOULBERG,S.P.HIRSHMAN,ORNL/TM-7871 (1981)*ABSH0470
C*****LAST REVISION: 6/81 L.L.LAO, R.M.WIELAND, AND W.A.HOULBERG ORNL. *ABSH0471
C*****CALCULATED PARAMETERS:                                           *ABSH0472
C*****AI-NORMALIZED INTEGRATED TOROIDAL CURRENT AT X-(DIMENSIONLESS).  *ABSH0473
C*****AIP-DERIVATIVE OF AI WITH RESPECT TO X-(DIMENSIONLESS).          *ABSH0474
C*****INPUT PARAMETERS:                                                *ABSH0475
C*****NRHO-NUMBER OF RADIAL GRID POINTS.                               *ABSH0476
C*****X-RADIAL POSITION-(DIMENSIONLESS).                               *ABSH0477
C*****CURT(I)-TOROIDAL CURRENT INTEGRATED TO NODE I.                   *ABSH0478
C*****RHO(I)-RADIAL COORDINATE FOR NODE I.                             *ABSH0479
C***********************************************************************ABSH0480
      PARAMETER MNMX=61                                                 ABSH0481
      COMMON /CINIT/ INITM, INITC, INITP, ITPR                          ABSH0482
      COMMON /SWORK1/ PRES(MNMX), CURT(MNMX), RHO(MNMX),                ABSH0483
     *     BPRES(MNMX), CPRES(MNMX), DPRES(MNMX), BCURT(MNMX),          ABSH0484
     *     CCURT(MNMX), DCURT(MNMX), NRHO                               ABSH0485
      IF (INITC.NE.0) GO TO 10                                          ABSH0486
      CALL SPLAAN(NRHO, RHO, CURT, BCURT, CCURT, DCURT)                 ABSH0487
      INITC = 1                                                         ABSH0488
   10 R = X*RHO(NRHO)                                                   ABSH0489
      AI = SEVAL(NRHO,R,RHO,CURT,BCURT,CCURT,DCURT)/CURT(NRHO)          ABSH0490
      AIP = SPEVAL(NRHO,R,RHO,CURT,BCURT,CCURT,DCURT)*RHO(NRHO)/        ABSH0491
     *     CURT(NRHO)                                                   ABSH0492
      RETURN                                                            ABSH0493
      END                                                               ABSH0494
      SUBROUTINE GEOTRN(MODE, XIN, EIN, DIN, EOUT, DOUT)                ABSH0495
C***********************************************************************ABSH0496
C*****GEOTRN TRANSFORMS ELONGATION AND TRIANGULARITY FROM A MOMENTS TO *ABSH0497
C*****A GEOMETRICAL REPRESENTATION AND VICE VERSA.                     *ABSH0498
C*****REFERENCES:                                                      *ABSH0499
C*****L.L.LAO,S.P.HIRSHMAN,R.M.WIELAND,ORNL/TM-7616 (1981).            *ABSH0500
C*****                                ,PHYS FLUIDS 24 (1981) AUG       *ABSH0501
C*****L.L.LAO,R.M.WIELAND,W.A.HOULBERG,S.P.HIRSHMAN,ORNL/TM-7871 (1981)*ABSH0502
C*****LAST REVISION: 6/81 L.L.LAO, R.M.WIELAND, AND W.A.HOULBERG ORNL. *ABSH0503
C*****CALCULATED PARAMETERS:                                           *ABSH0504
C*****EOUT-OUTPUT ELONGATION.                                          *ABSH0505
C*****DOUT-OUTPUT TRIANGULARITY.                                       *ABSH0506
C*****INPUT PARAMETERS:                                                *ABSH0507
C*****MODE-DESIGNATES DIRECTION OF TRANSFORMATION.                     *ABSH0508
C*****    =1 GEOMETRICAL => MOMENTS.                                   *ABSH0509
C*****    =2 MOMENTS => GEOMETRICAL.                                   *ABSH0510
C*****XIN-INPUT REDUCED MINOR RADIUS                                   *ABSH0511
C*****EIN-INPUT ELONGATION.                                            *ABSH0512
C*****DIN-INPUT TRIANGULARITY.                                         *ABSH0513
C***********************************************************************ABSH0514
      IF (XIN.LE.0.0) GO TO 30                                          ABSH0515
      IF (MODE.EQ.2) GO TO 20                                           ABSH0516
C*****GEOMETRICAL ==> MOMENTS REPRESENTATION.                           ABSH0517
      DOUT = DIN/4.0                                                    ABSH0518
      CTC = 0.0                                                         ABSH0519
      DO 10 I=1,10                                                      ABSH0520
           CTC = 4.0*DOUT/(SQRT(XIN**2+32.0*DOUT**2)+XIN)               ABSH0521
           DOUT = XIN*DIN/(4.0-6.0*CTC**2)                              ABSH0522
   10 CONTINUE                                                          ABSH0523
      EOUT = XIN*EIN/(SQRT(1.0-CTC*CTC)*(XIN+2.0*DOUT*CTC))             ABSH0524
      RETURN                                                            ABSH0525
C*****MOMENTS ==> GEOMETRICAL REPRESENTATION.                           ABSH0526
   20 CTC = 4.0*DIN/(SQRT(XIN**2+32.0*DIN**2)+XIN)                      ABSH0527
      STC = SQRT(1.0-CTC*CTC)                                           ABSH0528
      S2TC = 2.0*STC*CTC                                                ABSH0529
      C2TC = 2.0*CTC*CTC - 1.0                                          ABSH0530
      EOUT = EIN*(XIN*STC+DIN*S2TC)/XIN                                 ABSH0531
      DOUT = DIN*(1.0-3.0*C2TC)/XIN                                     ABSH0532
      RETURN                                                            ABSH0533
   30 EOUT = EIN                                                        ABSH0534
      DOUT = DIN                                                        ABSH0535
      RETURN                                                            ABSH0536
      END                                                               ABSH0537
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
      SUBROUTINE PRSFUN(X, P, PP)                                       ABSH0826
C***********************************************************************ABSH0827
C*****PRSFUN PROVIDES THE PLASMA PRESSURE PROFILE.                     *ABSH0828
C*****REFERENCES:                                                      *ABSH0829
C*****L.L.LAO,R.M.WIELAND,W.A.HOULBERG,S.P.HIRSHMAN,ORNL/TM-7871 (1981)*ABSH0830
C*****LAST REVISION: 6/81 L.L.LAO, R.M.WIELAND, AND W.A.HOULBERG ORNL. *ABSH0831
C*****CALCULATED PARAMETERS:                                           *ABSH0832
C*****P-NORMALIZED PRESSURE PROFILE AT X-(DIMENSIONLESS).              *ABSH0833
C*****PP-DERIVATIVE OF P WITH RESPECT TO X-(DIMENSIONLESS).            *ABSH0834
C*****INPUT PARAMETERS:                                                *ABSH0835
C*****NRHO-NUMBER OF RADIAL GRID POINTS.                               *ABSH0836
C*****X-RADIAL POSITION-(DIMENSIONLESS).                               *ABSH0837
C*****PRES(I)-PLASMA PRESSURE AT NODE I.                               *ABSH0838
C*****RHO(I)-RADIAL COORDINATE FOR NODE I.                             *ABSH0839
C***********************************************************************ABSH0840
      PARAMETER MNMX=61                                                 ABSH0841
      COMMON /CINIT/ INITM, INITC, INITP, ITPR                          ABSH0842
      COMMON /SWORK1/ PRES(MNMX), CURT(MNMX), RHO(MNMX),                ABSH0843
     *     BPRES(MNMX), CPRES(MNMX), DPRES(MNMX), BCURT(MNMX),          ABSH0844
     *     CCURT(MNMX), DCURT(MNMX), NRHO                               ABSH0845
      IF (INITP.NE.0) GO TO 10                                          ABSH0846
      CALL SPLAAN(NRHO, RHO, PRES, BPRES, CPRES, DPRES)                 ABSH0847
      INITP = 1                                                         ABSH0848
   10 R = X*RHO(NRHO)                                                   ABSH0849
      P = SEVAL(NRHO,R,RHO,PRES,BPRES,CPRES,DPRES)/PRES(1)              ABSH0850
      PP = SPEVAL(NRHO,R,RHO,PRES,BPRES,CPRES,DPRES)*RHO(NRHO)/         ABSH0851
     *     PRES(1)                                                      ABSH0852
      RETURN                                                            ABSH0853
      END                                                               ABSH0854
      SUBROUTINE PRSOL(PARAM, RES, N1, ERR)                             ABSH0855
C***********************************************************************ABSH0856
C*****PRSOL PRINTS OUT THE ITERATION PARAMETERS FOR THE MHD MOMENTS.   *ABSH0857
C*****REFERENCES:                                                      *ABSH0858
C*****L.L.LAO,S.P.HIRSHMAN,R.M.WIELAND,ORNL/TM-7616 (1981).            *ABSH0859
C*****                                ,PHYS FLUIDS 24 (1981) AUG       *ABSH0860
C*****L.L.LAO,R.M.WIELAND,W.A.HOULBERG,S.P.HIRSHMAN,ORNL/TM-7871 (1981)*ABSH0861
C*****LAST REVISION: 6/81 L.L.LAO, R.M.WIELAND, AND W.A.HOULBERG ORNL. *ABSH0862
C*****INPUT PARAMETERS:                                                *ABSH0863
C*****PARAM(I)-B.C. PARAMETER FOR FIRST ORDER ODE I.                   *ABSH0864
C*****RES-SUM OF SQUARES OF ERR(I).                                    *ABSH0865
C*****N1-NUMBER OF FIRST ORDER ODE'S BEING SOLVED.                     *ABSH0866
C*****ERR(I)-DIFFERENCE BETWEEN THE TWO SOLUTIONS AT MATCHING POINT.   *ABSH0867
C*****OTHER COMMENTS:                                                  *ABSH0868
C*****PRSOL IS REQUIRED BY D02AGF ODE SOLVER IN NAG LIBRARY.           *ABSH0869
C*****SEE D02AGF WRITE-UP FOR OPTIONAL DIAGNOSTICS.                    *ABSH0870
C***********************************************************************ABSH0871
      PARAMETER MNMOMS=3,MNEQ=6                                         ABSH0872
      COMMON /CINIT/ INITM, INITC, INITP, ITPR                          ABSH0873
      COMMON /NIO/ NIN, NOUT, NPLOT, NDEBUG, NTTY                       ABSH0874
      DIMENSION PARAM(MNEQ), ERR(MNEQ)                                  ABSH0875
      IF (ITPR.GT.0) GO TO 10                                           ABSH0876
      WRITE (NOUT,99999)                                                ABSH0877
      WRITE (NOUT,99998)                                                ABSH0878
   10 ITPR = ITPR + 1                                                   ABSH0879
      WRITE (NOUT,99997) ITPR, (PARAM(I),I=1,N1)                        ABSH0880
      WRITE (NOUT,99996) (ERR(I),I=1,N1), RES                           ABSH0881
      RETURN                                                            ABSH0882
99999 FORMAT (/, 10X, 40(1H*), 23H   VMOMS SAMPLE OUTPUT ,              ABSH0883
     *     40(1H*), //)                                                 ABSH0884
99998 FORMAT (//, 10X, 40(1H*), 26H NEWTON ITERATION SUMMARY ,          ABSH0885
     *     40(1H*))                                                     ABSH0886
99997 FORMAT (//, 10X, 13H ITERATION = , I3, /, 10X, 8H PARAM =,        ABSH0887
     *     1H , 6(1PE10.3, 2X))                                         ABSH0888
99996 FORMAT (10X, 9H ERR   = , 6(1PE10.3, 2X), 12H  RESIDUE = ,        ABSH0889
     *     1PE10.3)                                                     ABSH0890
      END                                                               ABSH0891
      SUBROUTINE RAAUX(XL, XR, RMATCH, PARAM)                           ABSH0892
C***********************************************************************ABSH0893
C*****RAAUX SETS THE MATCHING POINT FOR THE MHD MOMENTS.               *ABSH0894
C*****REFERENCES:                                                      *ABSH0895
C*****L.L.LAO,S.P.HIRSHMAN,R.M.WIELAND,ORNL/TM-7616 (1981).            *ABSH0896
C*****                                ,PHYS FLUIDS 24 (1981) AUG       *ABSH0897
C*****L.L.LAO,R.M.WIELAND,W.A.HOULBERG,S.P.HIRSHMAN,ORNL/TM-7871 (1981)*ABSH0898
C*****LAST REVISION: 6/81 L.L.LAO, R.M.WIELAND, AND W.A.HOULBERG ORNL. *ABSH0899
C*****CALCULATED PARAMETERS:                                           *ABSH0900
C*****XL-LEFT (LOWER) LIMIT ON ABSCISSA.                               *ABSH0901
C*****XR-RIGHT (UPPER) LIMIT ON ABSCISSA.                              *ABSH0902
C*****RMATCH-MATCHING POINT FOR SHOOTING METHOD (RMATCH=1.0).          *ABSH0903
C*****INPUT PARAMETERS:                                                *ABSH0904
C*****XL0-LOWER BOUND OF ABSCISSA.                                     *ABSH0905
C*****OTHER COMMENTS:                                                  *ABSH0906
C*****RAAUX IS REQUIRED BY D02AGF ODE SOLVER IN NAG LIBRARY.           *ABSH0907
C***********************************************************************ABSH0908
      PARAMETER MNMOMS=3,MNEQ=6                                         ABSH0909
      COMMON /LOCMEQ/ XL0, DX, R0, A0, BT0, P0, AI0, E0, E1,            ABSH0910
     *     D0, ARC, SLPE, SLPR0, BETAP0, KMHD, IUNDER                   ABSH0911
      DIMENSION PARAM(MNEQ)                                             ABSH0912
      XL = XL0                                                          ABSH0913
      XR = 1.0                                                          ABSH0914
      RMATCH = 1.0                                                      ABSH0915
      RETURN                                                            ABSH0916
      END                                                               ABSH0917
      FUNCTION SDOT(N, SX, INCX, SY, INCY)                              ABSH0918
C***********************************************************************ABSH0919
C*****SDOT CALCULATES THE DOT PRODUCT OF THE VECTORS SX AND SY.        *ABSH0920
C*****LAST REVISION: 6/81 R.M.WIELAND AND W.A.HOULBERG ORNL.           *ABSH0921
C*****INPUT PARAMETERS:                                                *ABSH0922
C*****N-NUMBER OF ELEMENTS TO BE SUMMED.                               *ABSH0923
C*****SX,SY-VECTORS TO BE DOTTED.                                      *ABSH0924
C*****INCX,INCY-INCREMENTS IN SX,SY ARGUMENTS.                         *ABSH0925
C*****OTHER COMMENTS:                                                  *ABSH0926
C*****USE OMNILIB VERSION FOR CRAY OPTIMIZATION.                       *ABSH0927
C***********************************************************************ABSH0928
      DIMENSION SX(N), SY(N)                                            ABSH0929
      SUM = 0.0                                                         ABSH0930
      IX = 1                                                            ABSH0931
      IY = 1                                                            ABSH0932
      DO 10 I=1,N                                                       ABSH0933
           SUM = SUM + SX(IX)*SY(IY)                                    ABSH0934
           IX = IX + INCX                                               ABSH0935
           IY = IY + INCY                                               ABSH0936
   10 CONTINUE                                                          ABSH0937
      SDOT = SUM                                                        ABSH0938
      RETURN                                                            ABSH0939
      END                                                               ABSH0940
      FUNCTION SEVAL(N, U, X, Y, B, C, D)                               ABSH0941
C***********************************************************************ABSH0942
C*****SEVAL EVALUATES THE CUBIC SPLINE FUNCTI0N.                       *ABSH0943
C*****REFERENCES:                                                      *ABSH0944
C*****G.E.FORSYTHE,M.A.MALCOLM,C.B.MOLER, COMPUTER METHODS FOR MATHE-  *ABSH0945
C***** MATICAL COMPUTATIONS, PRENTICE-HALL, 1977, P.76.                *ABSH0946
C*****LAST REVISION: 6/81 R.M.WIELAND AND W.A.HOULBERG ORNL.           *ABSH0947
C*****INPUT PARAMETERS:                                                *ABSH0948
C*****N-NUMBER OF DATA POINTS.                                         *ABSH0949
C*****U-ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED.               *ABSH0950
C*****X-ARRAY CONTAINING THE DATA ABCISSAS.                            *ABSH0951
C*****Y-ARRAY CONTAINING THE DATA ORDINATES.                           *ABSH0952
C*****B,C,D-ARRAYS OF SPLINE COEFFICIENTS COMPUTED BY SPLINE.          *ABSH0953
C*****OTHER COMMENTS:                                                  *ABSH0954
C*****SEVAL=Y(I)+B(I)*(U-X(I))+C(I)*(U-X(I))**2+D(I)*(U-X(I))**3       *ABSH0955
C*****WHERE X(I).LT.U.LT.X(I+1), USING HORNER'S RULE.                  *ABSH0956
C*****IF U.LT.X(1) THEN I=1 IS USED.                                   *ABSH0957
C*****IF U.GE.X(N) THEN I=N IS USED.                                   *ABSH0958
C*****IF U IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL, THEN A    *ABSH0959
C*****BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER INTERVAL.     *ABSH0960
C***********************************************************************ABSH0961
      DIMENSION X(N), Y(N), B(N), C(N), D(N)                            ABSH0962
      DATA I /1/                                                        ABSH0963
      IF (I.GE.N) I = 1                                                 ABSH0964
      IF (U.LT.X(I)) GO TO 10                                           ABSH0965
      IF (U.LE.X(I+1)) GO TO 30                                         ABSH0966
C*****BINARY SEARCH.                                                    ABSH0967
   10 I = 1                                                             ABSH0968
      J = N + 1                                                         ABSH0969
   20 K = (I+J)/2                                                       ABSH0970
      IF (U.LT.X(K)) J = K                                              ABSH0971
      IF (U.GE.X(K)) I = K                                              ABSH0972
      IF (J.GT.I+1) GO TO 20                                            ABSH0973
C*****EVALUATE SPLINE.                                                  ABSH0974
   30 DX = U - X(I)                                                     ABSH0975
      SEVAL = Y(I) + DX*(B(I)+DX*(C(I)+DX*D(I)))                        ABSH0976
      RETURN                                                            ABSH0977
      END                                                               ABSH0978
      SUBROUTINE SINCOS                                                 ABSH0979
C***********************************************************************ABSH0980
C*****SINCOS SETS UP THE SINE AND COSINE ARRAYS FOR THE MHD MOMENTS.   *ABSH0981
C*****REFERENCES:                                                      *ABSH0982
C*****L.L.LAO,S.P.HIRSHMAN,R.M.WIELAND,ORNL/TM-7616 (1981).            *ABSH0983
C*****                                ,PHYS FLUIDS 24 (1981) AUG       *ABSH0984
C*****L.L.LAO,R.M.WIELAND,W.A.HOULBERG,S.P.HIRSHMAN,ORNL/TM-7871 (1981)*ABSH0985
C*****LAST REVISION: 6/81 L.L.LAO, R.M.WIELAND, AND W.A.HOULBERG ORNL. *ABSH0986
C***********************************************************************ABSH0987
      PARAMETER MTHETA=11,MM1=10,MNMOMS=3                               ABSH0988
      PARAMETER MNEQ=6                                                  ABSH0989
      COMMON /SCFUNC/ S1TH(MTHETA), S2TH(MTHETA), C1TH(MTHETA),         ABSH0990
     *     C2TH(MTHETA)                                                 ABSH0991
      PI = ACOS(-1.0)                                                   ABSH0992
      DTH = PI/(MTHETA-1)                                               ABSH0993
      DO 10 I=1,MTHETA                                                  ABSH0994
           X1 = (I-1)*DTH                                               ABSH0995
           X2 = 2.0*X1                                                  ABSH0996
           S1TH(I) = SIN(X1)                                            ABSH0997
           S2TH(I) = SIN(X2)                                            ABSH0998
           C1TH(I) = COS(X1)                                            ABSH0999
           C2TH(I) = COS(X2)                                            ABSH1000
   10 CONTINUE                                                          ABSH1001
      RETURN                                                            ABSH1002
      END                                                               ABSH1003
      FUNCTION SPEVAL(N, U, X, Y, B, C, D)                              ABSH1004
C***********************************************************************ABSH1005
C*****SPEVAL EVALUATES THE DERIVATIVE OF THE CUBIC SPLINE FUNCTI0N.    *ABSH1006
C*****REFERENCES:                                                      *ABSH1007
C*****G.E.FORSYTHE,M.A.MALCOLM,C.B.MOLER, COMPUTER METHODS FOR MATHE-  *ABSH1008
C***** MATICAL COMPUTATIONS, PRENTICE-HALL, 1977, P.76.                *ABSH1009
C*****LAST REVISION: 6/81 R.M.WIELAND AND W.A.HOULBERG ORNL.           *ABSH1010
C*****INPUT PARAMETERS:                                                *ABSH1011
C*****N-NUMBER OF DATA POINTS.                                         *ABSH1012
C*****U-ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED.               *ABSH1013
C*****X-ARRAY CONTAINING THE DATA ABCISSAS.                            *ABSH1014
C*****Y-ARRAY CONTAINING THE DATA ORDINATES.                           *ABSH1015
C*****B,C,D-ARRAYS OF SPLINE COEFFICIENTS COMPUTED BY SPLINE.          *ABSH1016
C*****OTHER COMMENTS:                                                  *ABSH1017
C*****SPEVAL=B(I)+2*C(I)*(U-X(I))+3*D(I)*(U-X(I))**2                   *ABSH1018
C*****WHERE X(I).LT.U.LT.X(I+1), USING HORNER'S RULE.                  *ABSH1019
C*****IF U.LT.X(1) THEN I=1 IS USED.                                   *ABSH1020
C*****IF U.GE.X(N) THEN I=N IS USED.                                   *ABSH1021
C*****IF U IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL, THEN A    *ABSH1022
C*****BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER INTERVAL.     *ABSH1023
C***********************************************************************ABSH1024
      DIMENSION X(N), Y(N), B(N), C(N), D(N)                            ABSH1025
      DATA I /1/                                                        ABSH1026
      IF (I.GE.N) I = 1                                                 ABSH1027
      IF (U.LT.X(I)) GO TO 10                                           ABSH1028
      IF (U.LE.X(I+1)) GO TO 30                                         ABSH1029
C*****BINARY SEARCH.                                                    ABSH1030
   10 I = 1                                                             ABSH1031
      J = N + 1                                                         ABSH1032
   20 K = (I+J)/2                                                       ABSH1033
      IF (U.LT.X(K)) J = K                                              ABSH1034
      IF (U.GE.X(K)) I = K                                              ABSH1035
      IF (J.GT.I+1) GO TO 20                                            ABSH1036
C*****EVALUATE SPLINE.                                                  ABSH1037
   30 DX = U - X(I)                                                     ABSH1038
      SPEVAL = B(I) + DX*(2.0*C(I)+3.0*DX*D(I))                         ABSH1039
      RETURN                                                            ABSH1040
      END                                                               ABSH1041
      SUBROUTINE SPLAAN(N, X, Y, B, C, D)                               ABSH1042
C***********************************************************************ABSH1043
C*****SPLAAN EVALUATES THE COEFFICIENTS FOR A CUBIC INTERPOLATING      *ABSH1044
C*****SPLINE FOR WHICH S'(X1)=0.                                       *ABSH1045
C*****REFERENCES:                                                      *ABSH1046
C*****MODIFICATION OF SPLINE ROUTINE IN BELOW REFERENCE BY R.M.WIELAND *ABSH1047
C*****FOR S'(X1)=0.                                                    *ABSH1048
C*****G.E.FORSYTHE,M.A.MALCOLM,C.B.MOLER, COMPUTER METHODS FOR MATHE-  *ABSH1049
C***** MATICAL COMPUTATIONS, PRENTICE-HALL, 1977, P.76.                *ABSH1050
C*****LAST REVISION: 6/81 R.M.WIELAND AND W.A.HOULBERG ORNL.           *ABSH1051
C*****CALCULATED PARAMETERS:                                           *ABSH1052
C*****B(I),C(I),D(I)-ARRAYS OF N SPLINE COEFFICIENTS.                  *ABSH1053
C*****Y(I)=S(X(I)).                                                    *ABSH1054
C*****B(I)=S'(X(I)).                                                   *ABSH1055
C*****C(I)=S''(X(I)).                                                  *ABSH1056
C*****D(I)=S'''(X(I)).                                                 *ABSH1057
C*****INPUT PARAMETERS:                                                *ABSH1058
C*****N-NUMBER OF DATA POINTS OR KNOTS-(N.GE.2).                       *ABSH1059
C*****X(I)-ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING ORDER.        *ABSH1060
C*****Y(I)-ORDINATES OF THE KNOTS.                                     *ABSH1061
C*****OTHER COMMENTS:                                                  *ABSH1062
C*****S(X)=Y(I)+B(I)*(X-X(I))+C(I)*(X-X(I))**2+D(I)*(X-X(I))**3        *ABSH1063
C*****FOR X(I).LE.X.LE.X(I+1)                                          *ABSH1064
C*****SEVAL CAN BE USED TO EVALUATE THE SPLINE.                        *ABSH1065
C*****SPEVAL CAN BE USED TO EVALUATE THE DERIVATIVE OF THE SPLINE.     *ABSH1066
C***********************************************************************ABSH1067
      DIMENSION X(N), Y(N), B(N), C(N), D(N)                            ABSH1068
      NM1 = N - 1                                                       ABSH1069
      IF (N.LT.2) RETURN                                                ABSH1070
      IF (N.LT.3) GO TO 60                                              ABSH1071
C*****SET UP TRIDIAGONAL SYSTEM:                                        ABSH1072
C*****B=DIAGONAL, D=OFFDIAGONAL, C=RIGHT HAND SIDE.                     ABSH1073
      D(1) = X(2) - X(1)                                                ABSH1074
      C(2) = (Y(2)-Y(1))/D(1)                                           ABSH1075
      DO 10 I=2,NM1                                                     ABSH1076
           D(I) = X(I+1) - X(I)                                         ABSH1077
           B(I) = 2.0*(D(I-1)+D(I))                                     ABSH1078
           C(I+1) = (Y(I+1)-Y(I))/D(I)                                  ABSH1079
           C(I) = C(I+1) - C(I)                                         ABSH1080
   10 CONTINUE                                                          ABSH1081
C*****END CONDITIONS:                                                   ABSH1082
C*****THIRD DERIVATIVES AT 1 AND N OBTAINED FROM DIVIDED DIFFERENCES.   ABSH1083
      B(1) = 2.0*D(1)                                                   ABSH1084
      B(N) = -D(N-1)                                                    ABSH1085
      C(1) = 0.0                                                        ABSH1086
      C(N) = 0.0                                                        ABSH1087
      IF (N.EQ.3) GO TO 20                                              ABSH1088
      C(1) = (Y(2)-Y(1))/D(1)                                           ABSH1089
      C(N) = C(N-1)/(X(N)-X(N-2)) - C(N-2)/(X(N-1)-X(N-3))              ABSH1090
      C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))                              ABSH1091
C*****FORWARD ELIMINATION.                                              ABSH1092
   20 DO 30 I=2,N                                                       ABSH1093
           T = D(I-1)/B(I-1)                                            ABSH1094
           B(I) = B(I) - T*D(I-1)                                       ABSH1095
           C(I) = C(I) - T*C(I-1)                                       ABSH1096
   30 CONTINUE                                                          ABSH1097
C*****BACK SUBSTITUTION.                                                ABSH1098
      C(N) = C(N)/B(N)                                                  ABSH1099
      DO 40 IB=1,NM1                                                    ABSH1100
           I = N - IB                                                   ABSH1101
           C(I) = (C(I)-D(I)*C(I+1))/B(I)                               ABSH1102
   40 CONTINUE                                                          ABSH1103
C*****C(I) IS NOW THE SIGMA(I) OF THE TEXT.                             ABSH1104
C*****COMPUTE POLYNOMIAL COEFFICIENTS.                                  ABSH1105
      B(N) = (Y(N)-Y(NM1))/D(NM1) + D(NM1)*(C(NM1)+2.0*C(N))            ABSH1106
      DO 50 I=1,NM1                                                     ABSH1107
           B(I) = (Y(I+1)-Y(I))/D(I) - D(I)*(C(I+1)+2.0*C(I))           ABSH1108
           D(I) = (C(I+1)-C(I))/D(I)                                    ABSH1109
           C(I) = 3.0*C(I)                                              ABSH1110
   50 CONTINUE                                                          ABSH1111
      C(N) = 3.0*C(N)                                                   ABSH1112
      D(N) = D(N-1)                                                     ABSH1113
      RETURN                                                            ABSH1114
C*****COEFFICIENTS FOR N=2.                                             ABSH1115
   60 B(1) = (Y(2)-Y(1))/(X(2)-X(1))                                    ABSH1116
      C(1) = 0.0                                                        ABSH1117
      D(1) = 0.0                                                        ABSH1118
      B(2) = B(1)                                                       ABSH1119
      C(2) = 0.0                                                        ABSH1120
      D(2) = 0.0                                                        ABSH1121
      RETURN                                                            ABSH1122
      END                                                               ABSH1123
      SUBROUTINE D02AGF(H, ERROR, PARERR, PARAM, C, N, N1, M1,          ABSH1124
     *     AUX, BCAUX, RAAUX, PRSOL, MAT, COPY, WSPACE, WSPAC1,         ABSH1125
     *     WSPAC2, IFAIL)                                               ABSH1126
C     NAG COPYRIGHT 1975                                                ABSH1127
C     EDITED BY JOYCE CLARKE OXFORD OEG NUCLEAR PHYSICS 11TH SEP 1976   ABSH1128
C                   FORTRAN MACRO VERSION FDIA24.TEC                    ABSH1129
C     MARK 4.5 REVISED                                                  ABSH1130
C     ADDITIONAL COMMENTS BY R. WIELAND // ORNL // JAN.,1981.           ABSH1131
C                                                                       ABSH1132
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC ABSH1133
C                                                                       ABSH1134
C   THE PARAMETERS ARE DEFINED AS FOLLOWS:                              ABSH1135
C                                                                       ABSH1136
C     H           INPUT ESTIMATE OF REQUIRED INTEGRATION STEP LENGTH.   ABSH1137
C                 CHANGED ON RETURN TO PROVIDE LAST STEP LENGTH USED.   ABSH1138
C                                                                       ABSH1139
C     ERROR(N)    INPUT A REAL ARRAY USED AS:                           ABSH1140
C                 A) A LOCAL ERROR BOUND FOR INTEGRATION                ABSH1141
C                 B) A CONVERGENCE TOLERANCE ON Y AT THE MATCHING POINT.ABSH1142
C                    ON RETURN, THE FINAL ERROR.                        ABSH1143
C                                                                       ABSH1144
C     PARERR(N1)  INPUT A REAL ARRAY USED AS:                           ABSH1145
C                 A) A CONVERGENCE TOLERANCE FOR THE PARAMETERS (P)     ABSH1146
C                 B) USED TO APPROXIMATE DELTA-P IN ESTIMATING THE      ABSH1147
C                    JACOBIAN. ON RETURN, THE FINAL ERROR.              ABSH1148
C                                                                       ABSH1149
C     PARAM(N1)   INPUT STARTING VALUES FOR P;                          ABSH1150
C                 ON RETURN, THE CORRECTED VALUES                       ABSH1151
C                                                                       ABSH1152
C     C(M1,N)     THE SOLUTION Y(I;J) OF THE J-TH COMPONENT OF Y        ABSH1153
C                 EVALUATED AT X(I) IS RETURNED IN C (I,J); THE X(I)    ABSH1154
C                 SPACING IS DETERMINED BY M1, BELOW.                   ABSH1155
C                                                                       ABSH1156
C     N           NUMBER OF ODE'S                                       ABSH1157
C                                                                       ABSH1158
C     N1          NUMBER OF PARAMETERS                                  ABSH1159
C                                                                       ABSH1160
C     M1          THE FINAL SOLUTION IS CALCULATED AT M1 EQUIDISTANT    ABSH1161
C                 POINTS.                                               ABSH1162
C                                                                       ABSH1163
C     AUX         A USER SUPPLIED SUBROUTINE:                           ABSH1164
C                 SUBROUTINE AUX (F,Y,X,PARAM)                          ABSH1165
C                 WHERE X,Y(N),PARAM(N1) ARE USED TO EVALUTE THE        ABSH1166
C                 N DERIVATIVES F(N) AT X.                              ABSH1167
C                                                                       ABSH1168
C     BCAUX       USER SUPPLIED SUBROUTINE:                             ABSH1169
C                 SUBROUTINE BCAUX(G,G1,PARAM)                          ABSH1170
C                 WHERE PARAM IS USED (IF NECESSARY) TO EVALUATE Y(N)   ABSH1171
C                 AT THE ENDPOINTS X AND X1 AND RETURN THEM IN G(N)     ABSH1172
C                 AND G1(N).                                            ABSH1173
C                                                                       ABSH1174
C     RAAUX       USER SUPPLIED SUBROUTINE:                             ABSH1175
C                 SUBROUTINE RAAUX(X,X1,R,PARAM)                        ABSH1176
C                 WHERE PARAM IS USED (IF NECESSARY) TO EVALUATE THE    ABSH1177
C                 ENDPOINTS X AND X1, AND THE MATCHING POINT R.         ABSH1178
C                                                                       ABSH1179
C     PRSOL       A USER SUPPLIED SUBROUTINE:                           ABSH1180
C                 SUBROUTINE PRSOL(PARAM,RES,N1,ERR)                    ABSH1181
C                 CALLED AT EACH NEWTON ITERATION; CAN BE USED TO OUTPUTABSH1182
C                 ANY OF THE PARAMETERS OF INTEREST; ERR(N) ARE THE     ABSH1183
C                 ERRORS AT R IN EACH COMPONENT Y, AND RES IS THE       ABSH1184
C                 SUM OF THE SQUARES OF THESE ERRORS.                   ABSH1185
C                                                                       ABSH1186
C     MAT(N1,N1)   A REAL WORK ARRAY                                    ABSH1187
C     COPY(N1,N1)  A REAL WORK ARRAY                                    ABSH1188
C     WSPACE(N,9) A REAL WORK ARRAY                                     ABSH1189
C     WSPAC1(N)   A REAL WORK ARRAY                                     ABSH1190
C     WSPAC2(N)   A REAL WORK ARRAY                                     ABSH1191
C                                                                       ABSH1192
C     IFAIL       AN ERROR FLAG; ON INPUT, ENTER 0, ON OUTPUT:          ABSH1193
C       IFAIL = 0   NORMAL RETURN                                       ABSH1194
C             = 1   N1 > N                                              ABSH1195
C             = 2   INTEGRATION FAILED TO CONVERGE WHILE CALC. JACOBIAN ABSH1196
C             = 3   THE CONDITION X <= R <= X1 DOES NOT HOLD.           ABSH1197
C             = 4   INTEGRATION FAILED TO CONVERGE                      ABSH1198
C             = 5   JACOBIAN IS SINGULAR                                ABSH1199
C             = 6   AFTER 3 ATTEMPTS AT HALVING DELTA-P, THE NEWTON-    ABSH1200
C                   RAPHSON LOOP STILL DOES NOT YIELD A DIMINISHING     ABSH1201
C                   RESIDUAL.                                           ABSH1202
C             = 7   AFTER 12 N-R ITERATIONS, THERE IS STILL NOT         ABSH1203
C                   SUFFICIENT CONVERGENCE.                             ABSH1204
C                                                                       ABSH1205
C    INTERNAL FLAG :                                                    ABSH1206
C    VRBOSE (LOGICAL)   .TRUE. TO GENERATE DEBUG OUTPUT ON NDEBUG       ABSH1207
C                                                                       ABSH1208
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC ABSH1209
      INTEGER COUNT, COUNT1, EM, ONE, CT, P01AAF, N, N1, M1,            ABSH1210
     *     IFAIL, M, I, K, ITEST, J                                     ABSH1211
      DOUBLE PRECISION SRNAME                                           ABSH1212
      REAL MAT(N1,N1), H, ERROR(N), PARERR(N1), PARAM(N1),              ABSH1213
     *     C(M1,N), COPY(N1,N1), WSPACE(N,9), WSPAC1(N),                ABSH1214
     *     WSPAC2(N), EPS, H0, X, X1, R, DUM, RESID, D, PERT,           ABSH1215
     *     OLDRES, DIST, C1, X02AAF                                     ABSH1216
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1217
      COMMON /CFEVAL/ NFL                                               ABSH1218
      COMMON /NIO/ NIN, NOUT, NPLOT, NDEBUG, NTTY                       ABSH1219
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1220
      EXTERNAL AUX                                                      ABSH1221
      LOGICAL VRBOSE                                                    ABSH1222
      DATA VRBOSE /.FALSE./                                             ABSH1223
      DATA SRNAME /8H D02AGF /                                          ABSH1224
C                                                                       ABSH1225
C     SOLVES A GENERAL BOUNDARY VALUE                                   ABSH1226
C     PROBLEM FOR N DIFFERENTIAL EQUATIONS                              ABSH1227
C     IN N1 PARAMETERS USING A SHOOTING                                 ABSH1228
C     AND MATCHING TECHNIQUE.  EPS IS THE                               ABSH1229
C     LARGEST REAL VARIABLE SUCH THAT 1+EPS=1                           ABSH1230
C     ALL IMPLICITLY DECLARED REALS MAY BE USED DOUBLE-LENGTH           ABSH1231
C     THE ARRAY COPY IS REDUNDANT                                       ABSH1232
      EPS = SPMPAR(1)                                                   ABSH1233
      M = M1 - 1                                                        ABSH1234
      IF (N1.LE.N) GO TO 10                                             ABSH1235
C                                                                       ABSH1236
C     *** IFAIL = 1 ***                                                 ABSH1237
      EM = 1                                                            ABSH1238
      GO TO 480                                                         ABSH1239
C                                                                       ABSH1240
C     *** SET NEWTON ITERATION COUNTER ***                              ABSH1241
   10 COUNT = 0                                                         ABSH1242
      H0 = H                                                            ABSH1243
      ONE = 1                                                           ABSH1244
      EM = -1                                                           ABSH1245
C                                                                       ABSH1246
C     FORMS THE RESIDUALS AT THE                                        ABSH1247
C     MATCHING POINT                                                    ABSH1248
C                                                                       ABSH1249
C     *** GET LEFT-HAND PT., RIGHT-HAND PT., MATCHING PT. ***           ABSH1250
   20 CALL RAAUX(X, X1, R, PARAM)                                       ABSH1251
      IF ((X-R)*(X1-R).LE.0.0) GO TO 30                                 ABSH1252
C                                                                       ABSH1253
C     *** IFAIL = 3 ***                                                 ABSH1254
      EM = 3                                                            ABSH1255
      GO TO 480                                                         ABSH1256
C                                                                       ABSH1257
   30 IF (H0*(X1-X).LT.0.0) H0 = -H0                                    ABSH1258
C     ### G(I)  ==> W1            ###                                   ABSH1259
C     ### G1(I) ==> W2 ==> W(I,8) ###                                   ABSH1260
      CALL BCAUX(WSPAC1, WSPAC2, PARAM)                                 ABSH1261
      H = H0                                                            ABSH1262
      DO 40 I=1,N                                                       ABSH1263
           WSPACE(I,8) = WSPAC2(I)                                      ABSH1264
   40 CONTINUE                                                          ABSH1265
      I = 1                                                             ABSH1266
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1267
      IF (VRBOSE) WRITE (NDEBUG,99999)                                  ABSH1268
      NFL = 0                                                           ABSH1269
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1270
C                                                                       ABSH1271
C     *** INTEGRATE FROM X ==> R // ON RETURN: Y(R-) ==> W1(I) ***      ABSH1272
      CALL D02AGZ(X, WSPAC1, ERROR, ONE, N, N1, I, R-X, H, AUX,         ABSH1273
     *     WSPACE, WSPAC2, PARAM)                                       ABSH1274
      IF (I.EQ.0) GO TO 50                                              ABSH1275
C                                                                       ABSH1276
C     *** IFAIL = 4 ***                                                 ABSH1277
      EM = 4                                                            ABSH1278
      GO TO 480                                                         ABSH1279
C                                                                       ABSH1280
C     ###  -Y(R-) ==> W(I,8) // G1(I) ==> W1(I) ###                     ABSH1281
   50 DO 60 I=1,N                                                       ABSH1282
           DUM = WSPACE(I,8)                                            ABSH1283
           WSPACE(I,8) = -WSPAC1(I)                                     ABSH1284
           WSPAC1(I) = DUM                                              ABSH1285
   60 CONTINUE                                                          ABSH1286
C                                                                       ABSH1287
C     *** INTEGRATE FROM X1 ==> R // ON RETURN: YR(+) ==> W1(I) ***     ABSH1288
      H = -H0                                                           ABSH1289
      I = 1                                                             ABSH1290
      CALL D02AGZ(X1, WSPAC1, ERROR, ONE, N, N1, I, R-X1, H,            ABSH1291
     *     AUX, WSPACE, WSPAC2, PARAM)                                  ABSH1292
      IF (I.EQ.0) GO TO 70                                              ABSH1293
C                                                                       ABSH1294
C     *** IFAIL = 4 ***                                                 ABSH1295
      EM = 4                                                            ABSH1296
      GO TO 480                                                         ABSH1297
C                                                                       ABSH1298
   70 RESID = 0.0                                                       ABSH1299
      CT = 0                                                            ABSH1300
C                                                                       ABSH1301
C     *** FORM ERROR RESIDUALS AT R:                  ***               ABSH1302
C     ###   R(I) = Y(R+) - Y(R-) ==> W(I,8) == S(I;P) ###               ABSH1303
C     ###   S(I;P) ==> W1(I)                          ###               ABSH1304
      DO 80 I=1,N1                                                      ABSH1305
           D = WSPAC1(I)                                                ABSH1306
           DUM = WSPACE(I,8)                                            ABSH1307
           WSPACE(I,8) = D + DUM                                        ABSH1308
           D = 1.0 + ABS(DUM) + ABS(D)                                  ABSH1309
           DUM = WSPACE(I,8)                                            ABSH1310
           IF (ABS(DUM).LT.ERROR(I)*D) CT = CT + 1                      ABSH1311
           WSPAC1(I) = DUM                                              ABSH1312
           RESID = RESID + DUM*DUM                                      ABSH1313
   80 CONTINUE                                                          ABSH1314
      CALL PRSOL(PARAM, RESID, N1, WSPAC1)                              ABSH1315
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1316
      IF (VRBOSE) WRITE (NDEBUG,99998) CT, NFL                          ABSH1317
      IF (VRBOSE) WRITE (NDEBUG,99997) RESID                            ABSH1318
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1319
C     *** EM = -1 FIRST TIME THRU ONLY ***                              ABSH1320
      IF (EM.NE.-1) GO TO 320                                           ABSH1321
   90 COUNT = COUNT + 1                                                 ABSH1322
      IF (COUNT.NE.12) GO TO 100                                        ABSH1323
C                                                                       ABSH1324
C     *** IFAIL = 7 ***                                                 ABSH1325
      EM = 7                                                            ABSH1326
      GO TO 480                                                         ABSH1327
C                                                                       ABSH1328
C     FORMS THE JACOBIAN BY NUMERICAL                                   ABSH1329
C     DIFFERENTIATION                                                   ABSH1330
C     *** DEL$P(K)$ ==> PERT ***                                        ABSH1331
C     *** P(K) + DEL$P(K)$ ==> P(K+) ***                                ABSH1332
  100 CONTINUE                                                          ABSH1333
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1334
      NFL = 0                                                           ABSH1335
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1336
      DO 180 K=1,N1                                                     ABSH1337
           PERT = 10.0*PARERR(K)*(1.0+ABS(PARAM(K)))                    ABSH1338
           PARAM(K) = PERT + PARAM(K)                                   ABSH1339
           CALL RAAUX(X, X1, R, PARAM)                                  ABSH1340
           IF ((X-R)*(X1-R).LE.0.0) GO TO 110                           ABSH1341
C                                                                       ABSH1342
C     *** IFAIL = 3 ***                                                 ABSH1343
           EM = 3                                                       ABSH1344
           GO TO 480                                                    ABSH1345
C                                                                       ABSH1346
  110      IF (H0*(X1-X).LT.0.0) H0 = -H0                               ABSH1347
           H = H0                                                       ABSH1348
           CALL BCAUX(WSPAC1, WSPAC2, PARAM)                            ABSH1349
C     ### G1(I)$P(K+)$ ==> W(I,7) ###                                   ABSH1350
           DO 120 I=1,N                                                 ABSH1351
                WSPACE(I,7) = WSPAC2(I)                                 ABSH1352
  120      CONTINUE                                                     ABSH1353
           I = 1                                                        ABSH1354
C                                                                       ABSH1355
C     *** INTEGRATE FROM X ==> R USING P(K+) ***                        ABSH1356
           CALL D02AGZ(X, WSPAC1, ERROR, ONE, N, N1, I, R-X, H,         ABSH1357
     *          AUX, WSPACE, WSPAC2, PARAM)                             ABSH1358
           IF (I.EQ.0) GO TO 130                                        ABSH1359
C                                                                       ABSH1360
C     *** IFAIL = 2 ***                                                 ABSH1361
           EM = 2                                                       ABSH1362
           GO TO 480                                                    ABSH1363
C                                                                       ABSH1364
C     ### M(I,K) = Y(R-;P(K+)) ###                                      ABSH1365
  130      DO 140 I=1,N1                                                ABSH1366
                MAT(I,K) = WSPAC1(I)                                    ABSH1367
  140      CONTINUE                                                     ABSH1368
           H = -H0                                                      ABSH1369
           I = 1                                                        ABSH1370
           DO 150 I=1,N                                                 ABSH1371
                WSPAC1(I) = WSPACE(I,7)                                 ABSH1372
  150      CONTINUE                                                     ABSH1373
C                                                                       ABSH1374
C     *** INTEGRATE FROM X1 ==> R USING P(K+) ***                       ABSH1375
           CALL D02AGZ(X1, WSPAC1, ERROR, ONE, N, N1, I, R-X1,          ABSH1376
     *          H, AUX, WSPACE, WSPAC2, PARAM)                          ABSH1377
           IF (I.EQ.0) GO TO 160                                        ABSH1378
C                                                                       ABSH1379
C     *** IFAIL = 2 ***                                                 ABSH1380
           EM = 2                                                       ABSH1381
           GO TO 480                                                    ABSH1382
C                                                                       ABSH1383
C     ### M(I,K) = (S(P(K+)-S(P(K))/DEL$P(K)$ ###                       ABSH1384
  160      DO 170 I=1,N1                                                ABSH1385
                MAT(I,K) = (MAT(I,K)-WSPAC1(I)+WSPACE(I,8))/PERT        ABSH1386
                IF (ABS(MAT(I,K)).LT.5.0*EPS*ABS(WSPACE(I,8))/          ABSH1387
     *               PERT) MAT(I,K) = 0.0                               ABSH1388
  170      CONTINUE                                                     ABSH1389
           PARAM(K) = PARAM(K) - PERT                                   ABSH1390
  180 CONTINUE                                                          ABSH1391
C                                                                       ABSH1392
C     *** NEW JACOBIAN ***                                              ABSH1393
      ITEST = 1                                                         ABSH1394
      EM = -3                                                           ABSH1395
C                                                                       ABSH1396
C     PERFORMS COLUMN SCALING ON THE JACOBIAN                           ABSH1397
C     AND FORMS A TRIANGULAR DECOMPOSITION                              ABSH1398
      DO 220 I=1,N1                                                     ABSH1399
           D = 0.0                                                      ABSH1400
           DO 190 J=1,N1                                                ABSH1401
                IF (ABS(MAT(J,I)).GT.D) D = ABS(MAT(J,I))               ABSH1402
  190      CONTINUE                                                     ABSH1403
           IF (D.NE.0.0) GO TO 200                                      ABSH1404
C                                                                       ABSH1405
C     *** IFAIL = 5 ***                                                 ABSH1406
           EM = 5                                                       ABSH1407
           GO TO 480                                                    ABSH1408
C                                                                       ABSH1409
C     *** NORMALIZE M(I,K) ***                                          ABSH1410
  200      DO 210 J=1,N1                                                ABSH1411
                MAT(J,I) = MAT(J,I)/D                                   ABSH1412
  210      CONTINUE                                                     ABSH1413
           WSPACE(I,7) = D                                              ABSH1414
  220 CONTINUE                                                          ABSH1415
      I = 1                                                             ABSH1416
C                                                                       ABSH1417
C     *** LU DECOMPOSITION OF M(I,K) // PIVOT ARRAY: W1(I) ***          ABSH1418
      CALL F03AFF(N1, EPS, MAT, N1, D, J, WSPAC1, I)                    ABSH1419
      IF (I.EQ.0) GO TO 230                                             ABSH1420
C                                                                       ABSH1421
C     *** IFAIL = 5 ***                                                 ABSH1422
      EM = 5                                                            ABSH1423
      GO TO 480                                                         ABSH1424
C                                                                       ABSH1425
C     ### PIVOT ARRAY ==> W(I,6) ###                                    ABSH1426
  230 CONTINUE                                                          ABSH1427
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1428
      IF (VRBOSE) WRITE (NDEBUG,99996) NFL                              ABSH1429
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1430
      DO 240 I=1,N1                                                     ABSH1431
           WSPACE(I,6) = WSPAC1(I)                                      ABSH1432
  240 CONTINUE                                                          ABSH1433
C                                                                       ABSH1434
C     USES A GENERALISED NEWTON RAPHSON                                 ABSH1435
C     TECHNIQUE TO SOLVE THE NONLINEAR                                  ABSH1436
C     EQUATIONS AT THE MATCHING POINT                                   ABSH1437
C     *** SOLVE THE SYSTEM :  M(I,K) * DEL$P(K)$ = S(I;P) ***           ABSH1438
  250 OLDRES = RESID                                                    ABSH1439
C                                                                       ABSH1440
C     *** SET NEWTON - RAPHSON COUNTER ***                              ABSH1441
      COUNT1 = 0                                                        ABSH1442
C                                                                       ABSH1443
C     ### S(P) ==> W1(I) ###                                            ABSH1444
      DO 260 I=1,N1                                                     ABSH1445
           WSPAC1(I) = WSPACE(I,6)                                      ABSH1446
           WSPACE(I,1) = WSPACE(I,8)                                    ABSH1447
  260 CONTINUE                                                          ABSH1448
C                                                                       ABSH1449
C     *** THE NEW DEL$P$ = W1 ***                                       ABSH1450
      CALL F04AJF(N1, ONE, MAT, N1, WSPAC1, WSPACE, N)                  ABSH1451
C                                                                       ABSH1452
C     *** NORMALIZE DEL$P$ / FORM NEW P ***                             ABSH1453
      DO 270 I=1,N1                                                     ABSH1454
           WSPACE(I,1) = WSPACE(I,1)/WSPACE(I,7)                        ABSH1455
           PARAM(I) = PARAM(I) + WSPACE(I,1)                            ABSH1456
  270 CONTINUE                                                          ABSH1457
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1458
      IF (VRBOSE) WRITE (NDEBUG,99995) PARAM                            ABSH1459
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1460
      IF (CT.LT.N1) GO TO 300                                           ABSH1461
      DO 280 I=1,N1                                                     ABSH1462
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1463
           IPP = I                                                      ABSH1464
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1465
           IF (ABS(WSPACE(I,1)).GT.PARERR(I)*(1.0+ABS(PARAM(I)))        ABSH1466
     *          ) GO TO 290                                             ABSH1467
  280 CONTINUE                                                          ABSH1468
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1469
      IF (VRBOSE) WRITE (NDEBUG,99994)                                  ABSH1470
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1471
C                                                                       ABSH1472
C     *** SOLUTION FOUND] ***                                           ABSH1473
      EM = -5                                                           ABSH1474
      GO TO 380                                                         ABSH1475
  290 CONTINUE                                                          ABSH1476
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1477
      PIP = PARERR(IPP)*(1.0+ABS(PARAM(IPP)))                           ABSH1478
      IF (VRBOSE) WRITE (NDEBUG,99993) IPP, WSPACE(IPP,1), PIP          ABSH1479
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1480
  300 DO 310 I=1,N1                                                     ABSH1481
           WSPACE(I,1) = -WSPACE(I,1)                                   ABSH1482
  310 CONTINUE                                                          ABSH1483
      GO TO 20                                                          ABSH1484
C                                                                       ABSH1485
C     ******************************************************            ABSH1486
C                                                                       ABSH1487
  320 CONTINUE                                                          ABSH1488
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1489
      IF (VRBOSE) WRITE (NDEBUG,99992) OLDRES                           ABSH1490
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1491
      IF (COUNT1.NE.0) GO TO 330                                        ABSH1492
      IF (RESID.GE.OLDRES/10.0) GO TO 330                               ABSH1493
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1494
      IF (VRBOSE) WRITE (NDEBUG,99991)                                  ABSH1495
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1496
C                                                                       ABSH1497
C     *** GO BACK AND FORM NEW DEL$P$ USING CURRENT JACOBIAN, NEW S *** ABSH1498
      EM = -2                                                           ABSH1499
      ITEST = 0                                                         ABSH1500
      GO TO 250                                                         ABSH1501
C                                                                       ABSH1502
C     *** FORM NEW JACOBIAN ***                                         ABSH1503
  330 CONTINUE                                                          ABSH1504
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1505
      IF (VRBOSE) WRITE (NDEBUG,99990)                                  ABSH1506
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1507
      IF (RESID.LT.OLDRES) GO TO 90                                     ABSH1508
      IF (COUNT1.NE.3) GO TO 360                                        ABSH1509
      IF (ITEST.EQ.0) GO TO 340                                         ABSH1510
C                                                                       ABSH1511
C     *** IFAIL = 6 ***                                                 ABSH1512
      EM = 6                                                            ABSH1513
      GO TO 480                                                         ABSH1514
  340 CONTINUE                                                          ABSH1515
      DO 350 I=1,N1                                                     ABSH1516
           PARAM(I) = PARAM(I) + WSPACE(I,1)                            ABSH1517
  350 CONTINUE                                                          ABSH1518
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1519
      IF (VRBOSE) WRITE (NDEBUG,99989) PARAM                            ABSH1520
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1521
      EM = -1                                                           ABSH1522
      GO TO 20                                                          ABSH1523
  360 COUNT1 = COUNT1 + 1                                               ABSH1524
      EM = -4                                                           ABSH1525
C                                                                       ABSH1526
C     *** SCALE-DOWN DEL$P$ ***                                         ABSH1527
      DO 370 I=1,N1                                                     ABSH1528
           WSPACE(I,1) = WSPACE(I,1)/2.0                                ABSH1529
           PARAM(I) = PARAM(I) + WSPACE(I,1)                            ABSH1530
  370 CONTINUE                                                          ABSH1531
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1532
      IF (VRBOSE) WRITE (NDEBUG,99988) PARAM                            ABSH1533
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1534
      GO TO 20                                                          ABSH1535
C                                                                       ABSH1536
C     ******************************************************            ABSH1537
C                                                                       ABSH1538
C     CALCULATES THE FINAL SOLUTION                                     ABSH1539
  380 CONTINUE                                                          ABSH1540
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1541
      IF (VRBOSE) WRITE (NDEBUG,99987)                                  ABSH1542
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1543
      IF (M.LE.0) GO TO 480                                             ABSH1544
      CALL RAAUX(X, X1, R, PARAM)                                       ABSH1545
      IF ((X-R)*(X1-R).LE.0.0) GO TO 390                                ABSH1546
      EM = 3                                                            ABSH1547
      GO TO 480                                                         ABSH1548
  390 IF (H0*(X1-X).LT.0.0) H0 = -H0                                    ABSH1549
      H = H0                                                            ABSH1550
      CALL BCAUX(WSPAC1, WSPAC2, PARAM)                                 ABSH1551
      DO 400 I=1,N                                                      ABSH1552
           WSPACE(I,7) = WSPAC2(I)                                      ABSH1553
  400 CONTINUE                                                          ABSH1554
      DIST = (X1-X)/FLOAT(M)                                            ABSH1555
      J = 1                                                             ABSH1556
      C1 = X                                                            ABSH1557
      K = 1                                                             ABSH1558
  410 DO 420 I=1,N                                                      ABSH1559
           C(J,I) = WSPAC1(I)                                           ABSH1560
  420 CONTINUE                                                          ABSH1561
      IF ((R-C1-0.25*DIST)*DIST.LE.0.0) GO TO 440                       ABSH1562
      I = 1                                                             ABSH1563
      CALL D02AGZ(C1, WSPAC1, ERROR, ONE, N, N1, I, DIST, H,            ABSH1564
     *     AUX, WSPACE, WSPAC2, PARAM)                                  ABSH1565
      IF (I.EQ.0) GO TO 430                                             ABSH1566
      EM = 4                                                            ABSH1567
      GO TO 480                                                         ABSH1568
  430 J = J + K                                                         ABSH1569
      GO TO 410                                                         ABSH1570
  440 IF (K.EQ.-1) GO TO 460                                            ABSH1571
      DIST = -DIST                                                      ABSH1572
      C1 = X1                                                           ABSH1573
      H = -H0                                                           ABSH1574
      J = M1                                                            ABSH1575
      K = -1                                                            ABSH1576
      DO 450 I=1,N                                                      ABSH1577
           WSPAC1(I) = WSPACE(I,7)                                      ABSH1578
  450 CONTINUE                                                          ABSH1579
      GO TO 410                                                         ABSH1580
  460 CALL BCAUX(WSPAC1, WSPAC2, PARAM)                                 ABSH1581
      DO 470 I=1,N                                                      ABSH1582
           C(1,I) = WSPAC1(I)                                           ABSH1583
           C(M1,I) = WSPAC2(I)                                          ABSH1584
  470 CONTINUE                                                          ABSH1585
  480 IF (EM.LE.0) IFAIL = 0                                            ABSH1586
      IF (EM.GT.0) IFAIL = P01AAF(IFAIL,EM,SRNAME)                      ABSH1587
      RETURN                                                            ABSH1588
C     END OF D02AGF                                                     ABSH1589
99999 FORMAT (10H FWD MARCH)                                            ABSH1590
99998 FORMAT (2X, I2, 2H Y13HS MATCH WITH , I6, 10H AUX EVALS)          ABSH1591
99997 FORMAT (16H PRSOL: RESID = , 1PE11.2)                             ABSH1592
99996 FORMAT (19H NEW JACOBIAN WITH , I6, 10H AUX EVALS)                ABSH1593
99995 FORMAT (17H SOLVE FOR NEW P:, /, 5X, 1P6E10.2)                    ABSH1594
99994 FORMAT (22H DEL P IS SMALL ENOUGH)                                ABSH1595
99993 FORMAT (10H  DEL P # , I1, 17H STILL TOO LARGE:, 1PE11.3,         ABSH1596
     *     2H >, E11.3)                                                 ABSH1597
99992 FORMAT (28H COMPARE AGAINST OLD RESID- , 1PE11.3)                 ABSH1598
99991 FORMAT (44H KEEP THE OLD JAC; GET NEW DEL P USING NEW S)          ABSH1599
99990 FORMAT (39H NEW RESID IS NOT A FACTOR OF 10 BETTER)               ABSH1600
99989 FORMAT (39H RESID IS WORSE & LIM=3; RESTORE OLD P:, /,            ABSH1601
     *     5X, 1P6E10.2)                                                ABSH1602
99988 FORMAT (32H RESID IS WORSE; BACK DOWN DEL P, /, 5X,               ABSH1603
     *     6HNEW P:, /, 5X, 1P6E10.2)                                   ABSH1604
99987 FORMAT (9H FINISHED)                                              ABSH1605
      END                                                               ABSH1606
      SUBROUTINE D02AGZ(X, Y, G, T, N, N1, M, H0, H, AUX,               ABSH1607
     *     WSPACE, WSPAC1, PARAM)                                       ABSH1608
C                                                                       ABSH1609
C     USES THE LOGIC OF NAG LIBRARY ROUTINE D02ABF                      ABSH1610
C     EPS AND DUM MAY BE USED DOUBLE LENGTH.  EPS IS THE                ABSH1611
C     SMALLEST REAL SUCH THAT 1+EPS>EPS.  SMAX                          ABSH1612
C     IS THE LARGEST INTEGER. DUM,ERR AND HS MAYBE DECLARED             ABSH1613
C     DOUBLE PRECISION                                                  ABSH1614
C     NAG COPYRIGHT 1975                                                ABSH1615
C     EDITED BY JOYCE CLARKE OXFORD OEG NUCLEAR PHYSICS 11TH SEP 1976   ABSH1616
C                   FORTRAN MACRO VERSION FDIA24.TEC                    ABSH1617
C     MARK 4.5 REVISED                                                  ABSH1618
C     ADDITIONAL COMMENTS BY R. WIELAND // ORNL // JAN.,1981.           ABSH1619
C                                                                       ABSH1620
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC             ABSH1621
C     THE SUBROUTINE PARAMETERS ARE DEFINED AS FOLLOWS:                 ABSH1622
C                                                                       ABSH1623
C     X       INPUT INITAL X VALUE; ON RETURN X=X+H0 IF SUCCESSFUL      ABSH1624
C                                                                       ABSH1625
C     Y(N)    INPUT INITIAL VALUES OF Y(X); ON RETURN Y(X+H0) IF        ABSH1626
C             SUCCESSFUL                                                ABSH1627
C                                                                       ABSH1628
C     G(N)    INPUT ERROR BOUNDS ON Y(N)                                ABSH1629
C                                                                       ABSH1630
C     T       INPUT TYPE OF ERROR TEST                                  ABSH1631
C                                                                       ABSH1632
C     N       NUMBER OF ODE'S IN THE SYSTEM                             ABSH1633
C                                                                       ABSH1634
C     M       INPUT ERROR FLAG AS 0; ON RETURN:                         ABSH1635
C             M = 0   NORMAL RETURN                                     ABSH1636
C               = 1   STEP LENGTH REPEATEDLY HALVED TO < 10**-4 THE     ABSH1637
C                     INITIALLY SUGGESTED STEP LENGTH                   ABSH1638
C               = 2   T .NE. 1 OR 2  OR 3                               ABSH1639
C               = 3   THE NUMBER OF STEPS REQUIRED EXCEEDS THE LARGEST  ABSH1640
C                     INTEGER ON THE MACHINE                            ABSH1641
C                                                                       ABSH1642
C     H0      INTERVAL OF INTEGRATION                                   ABSH1643
C                                                                       ABSH1644
C     H       INPUT AN ESTIMATE OF THE STEP LENGTH REQUIRED;            ABSH1645
C             ON RETURN, THE FINAL STEP LENGTH USED.                    ABSH1646
C                                                                       ABSH1647
C      THE REMAINING PARAMETERS ARE AS DEFINED IN 'D02AGF'              ABSH1648
C                                                                       ABSH1649
C     ON RETURN:                                                        ABSH1650
C       WSPACE(I,2) CONTAINS THE LOCAL ERROR AT X0+H FOR EACH Y(I)      ABSH1651
C       WSPACE(I,5) CONTAINS Y'(X+H0)                                   ABSH1652
C       WSPACE(I,3) CONTAINS Y'(X)                                      ABSH1653
C                                                                       ABSH1654
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCABSH1655
C                                                                       ABSH1656
      INTEGER D, S, SMAX, T, N, N1, M, I, J                             ABSH1657
      REAL X, Y(N), G(N), H0, H, WSPACE(N,9), WSPAC1(N),                ABSH1658
     *     PARAM(N1), EPS, DUM, P, Q, HS, ERR, X02AAF                   ABSH1659
      EXTERNAL AUX                                                      ABSH1660
C     *** H0 = INTERVAL LENGTH // H = STEP LENGTH ***                   ABSH1661
      EPS = SPMPAR(1)                                                   ABSH1662
      SMAX = I1MACH(9)                                                  ABSH1663
      M = 0                                                             ABSH1664
      DUM = AMAX1(ABS(X),ABS(X+H0))                                     ABSH1665
      IF (ABS(H0).LE.EPS*DUM) RETURN                                    ABSH1666
      DUM = H/H0                                                        ABSH1667
      IF (ABS(DUM).GT.1.0 .OR. H.EQ.0.0) GO TO 10                       ABSH1668
      DUM = ABS(H0/H+0.9)                                               ABSH1669
      IF (DUM.GT.FLOAT(SMAX)) GO TO 80                                  ABSH1670
      H = H0/AINT(DUM)                                                  ABSH1671
      GO TO 20                                                          ABSH1672
   10 H = H0                                                            ABSH1673
C                                                                       ABSH1674
C     *** T = 1 :: MIXED    ERROR TEST    ***                           ABSH1675
C     *** T = 2 :: ABSOLUTE ERROR TEST    ***                           ABSH1676
C     *** T = 3 :: RELATIVE ERROR TEST    ***                           ABSH1677
   20 P = 1.0                                                           ABSH1678
      IF (T.EQ.3) P = 0.0                                               ABSH1679
      Q = 1.0                                                           ABSH1680
      IF (T.EQ.2) Q = 0.0                                               ABSH1681
      DUM = H0/H + 0.1                                                  ABSH1682
C     *** S = # OF ANTICIPATED STEPS       ***                          ABSH1683
C     *** HS = MINIMUM STEP LENGTH ALLOWED ***                          ABSH1684
      S = INT(DUM)                                                      ABSH1685
      HS = 1.0E-4*ABS(H)                                                ABSH1686
   30 DO 40 I=1,N                                                       ABSH1687
           WSPACE(I,9) = Y(I)                                           ABSH1688
   40 CONTINUE                                                          ABSH1689
C                                                                       ABSH1690
C     *** INTEGRATE FROM X ==> X+H;                              ***    ABSH1691
C     ***   INPUT:  Y(X) ==> Y // X ==> X;                       ***    ABSH1692
C     ***   OUTPUT: Y(X+H) ==> Y // X+H ==> X;                   ***    ABSH1693
C     ***           LOCAL ERROR ==> W(I,2) // Y'(X+H) ==> W(I,5) ***    ABSH1694
C     ***           Y'(X) ==> W(I,3)                             ***    ABSH1695
      CALL D02AGY(Y, X, H, N, N1, AUX, WSPACE, WSPAC1, PARAM)           ABSH1696
      D = 0                                                             ABSH1697
      DO 50 I=1,N                                                       ABSH1698
           ERR = G(I)*(P+Q*ABS(Y(I)))                                   ABSH1699
           IF (WSPACE(I,2).GT.ERR) GO TO 60                             ABSH1700
           IF (40.0*WSPACE(I,2).GT.ERR) D = 1                           ABSH1701
   50 CONTINUE                                                          ABSH1702
      S = S - 1                                                         ABSH1703
C                                                                       ABSH1704
C     *** THE ONLY SUCCESSFUL WAY OUT ***                               ABSH1705
      IF (S.EQ.0) RETURN                                                ABSH1706
C                                                                       ABSH1707
      DUM = FLOAT(S)/2.0 + 0.1                                          ABSH1708
      J = INT(DUM)*2                                                    ABSH1709
      IF (D.NE.0 .OR. J.NE.S) GO TO 30                                  ABSH1710
      H = 2.0*H                                                         ABSH1711
      S = S/2                                                           ABSH1712
      GO TO 30                                                          ABSH1713
C                                                                       ABSH1714
C     *** LOCAL ERROR TOO LARGE ***                                     ABSH1715
   60 X = X - H                                                         ABSH1716
      DO 70 I=1,N                                                       ABSH1717
           Y(I) = WSPACE(I,9)                                           ABSH1718
   70 CONTINUE                                                          ABSH1719
      IF (S.GT.SMAX/2) GO TO 80                                         ABSH1720
      S = 2*S                                                           ABSH1721
      H = 0.5*H                                                         ABSH1722
      IF (ABS(H).GT.HS) GO TO 30                                        ABSH1723
   80 M = 1                                                             ABSH1724
      RETURN                                                            ABSH1725
      END                                                               ABSH1726
      SUBROUTINE D02AGY(Y, X, H, N, N1, AUX, WSPACE, WSPAC1,            ABSH1727
     *     PARAM)                                                       ABSH1728
C                                                                       ABSH1729
C     USES THE TECHNIQUE OF NAG LIBRARY                                 ABSH1730
C     PROCEDURE D02AAF WITH THE                                         ABSH1731
C     SPECIFICATION OF AUX CHANGED                                      ABSH1732
C     ALL IMPLICITLY DECLARED REALS MAY                                 ABSH1733
C     BE DECLARED DOUBLE PRECISION                                      ABSH1734
C     NAG COPYRIGHT 1975                                                ABSH1735
C     EDITED BY JOYCE CLARKE OXFORD OEG NUCLEAR PHYSICS 11TH SEP 1976   ABSH1736
C                   FORTRAN MACRO VERSION FDIA24.TEC                    ABSH1737
C     MARK 4.5 REVISED                                                  ABSH1738
C                                                                       ABSH1739
C   INTEGRATES A SYSTEM OF ODE'S USING MERSON'S METH.                   ABSH1740
C   CF.  LAMBERT, J.D., COMPUTATIONAL METHODS IN ORDINARY DIFFERENTIAL  ABSH1741
C        EQUATIONS, PP. 130-135, WILEY, 1973.                           ABSH1742
C                                                                       ABSH1743
      INTEGER N, N1, I                                                  ABSH1744
      REAL Y(N), X, H, WSPACE(N,9), WSPAC1(N), PARAM(N1), C1,           ABSH1745
     *     C2, U, V, W, DUM                                             ABSH1746
      CALL AUX(WSPAC1, Y, X, PARAM)                                     ABSH1747
      C1 = 1.0/3.0                                                      ABSH1748
      C2 = 1.0/6.0                                                      ABSH1749
      DO 10 I=1,N                                                       ABSH1750
           WSPACE(I,3) = WSPAC1(I)                                      ABSH1751
   10 CONTINUE                                                          ABSH1752
      U = C1*H                                                          ABSH1753
      DO 20 I=1,N                                                       ABSH1754
           WSPACE(I,4) = Y(I)                                           ABSH1755
           Y(I) = Y(I) + U*WSPACE(I,3)                                  ABSH1756
   20 CONTINUE                                                          ABSH1757
      CALL AUX(WSPAC1, Y, U+X, PARAM)                                   ABSH1758
      DO 30 I=1,N                                                       ABSH1759
           WSPACE(I,5) = WSPAC1(I)                                      ABSH1760
   30 CONTINUE                                                          ABSH1761
      V = H*C2                                                          ABSH1762
      DO 40 I=1,N                                                       ABSH1763
           Y(I) = WSPACE(I,4) + V*(WSPACE(I,3)+WSPACE(I,5))             ABSH1764
   40 CONTINUE                                                          ABSH1765
      CALL AUX(WSPAC1, Y, U+X, PARAM)                                   ABSH1766
      DO 50 I=1,N                                                       ABSH1767
           WSPACE(I,5) = WSPAC1(I)                                      ABSH1768
   50 CONTINUE                                                          ABSH1769
      U = H*0.125                                                       ABSH1770
      V = H*0.375                                                       ABSH1771
      DO 60 I=1,N                                                       ABSH1772
           Y(I) = WSPACE(I,4) + WSPACE(I,3)*U + WSPACE(I,5)*V           ABSH1773
   60 CONTINUE                                                          ABSH1774
      U = 0.5*H                                                         ABSH1775
      V = 1.5*H                                                         ABSH1776
      W = 2.0*H                                                         ABSH1777
      CALL AUX(WSPAC1, Y, U+X, PARAM)                                   ABSH1778
      DO 70 I=1,N                                                       ABSH1779
           Y(I) = WSPACE(I,4) + WSPACE(I,3)*U - WSPACE(I,5)*V +         ABSH1780
     *          WSPAC1(I)*W                                             ABSH1781
           WSPACE(I,5) = WSPAC1(I)                                      ABSH1782
   70 CONTINUE                                                          ABSH1783
      X = X + H                                                         ABSH1784
      CALL AUX(WSPAC1, Y, X, PARAM)                                     ABSH1785
      U = H*C2                                                          ABSH1786
      V = 2.0*H*C1                                                      ABSH1787
      DO 80 I=1,N                                                       ABSH1788
           W = WSPACE(I,4) + U*(WSPACE(I,3)+WSPAC1(I)) +                ABSH1789
     *          WSPACE(I,5)*V                                           ABSH1790
           DUM = W - Y(I)                                               ABSH1791
           WSPACE(I,2) = 0.2*ABS(DUM)                                   ABSH1792
           Y(I) = W                                                     ABSH1793
   80 CONTINUE                                                          ABSH1794
      RETURN                                                            ABSH1795
      END                                                               ABSH1796
      SUBROUTINE F03AFF(N, EPS, A, IA, D1, ID, P, IFAIL)                ABSH1797
C     MARK 2 RELEASE. NAG COPYRIGHT 1972                                ABSH1798
C     EDITED BY JOYCE CLARKE OXFORD OEG NUCLEAR PHYSICS 11TH SEP 1976   ABSH1799
C                   FORTRAN MACRO VERSION FDIA24.TEC                    ABSH1800
C     MARK 3 REVISED.                                                   ABSH1801
C     MARK 4.5 REVISED                                                  ABSH1802
C     SW3 SWITCH REVISION : R. WIELAND // ORNL //   JAN., 1981.         ABSH1803
C     .FALSE. GIVES NO EXTENDED PRECISION                               ABSH1804
C     .TRUE. GIVES EXTENDED PRECISION IF CONDITIONS COMMENTED           ABSH1805
C     ON IN X03AAZ ARE FULFILLED.                                       ABSH1806
C                                                                       ABSH1807
C     UNSYMDET                                                          ABSH1808
C     THE UNSYMMETRIC MATRIX, A, IS STORED IN THE N*N ARRAY A(I,J),     ABSH1809
C     I=1,N, J=1,N. THE DECOMPOSITION A=LU, WHERE L IS A                ABSH1810
C     LOWER TRIANGULAR MATRIX AND U IS A UNIT UPPER TRIANGULAR          ABSH1811
C     MATRIX, IS PERFORMED AND OVERWRITTEN ON A, OMITTING THE UNIT      ABSH1812
C     DIAGONAL OF U. A RECORD OF ANY INTERCHANGES MADE TO THE ROWS      ABSH1813
C     OF A IS KEPT IN P(I), I=1,N, SUCH THAT THE I-TH ROW AND           ABSH1814
C     THE P(I)-TH ROW WERE INTERCHANGED AT THE I-TH STEP. THE           ABSH1815
C     DETERMINANT, D1 * 2.0**ID, OF A IS ALSO COMPUTED. THE             ABSH1816
C     SUBROUTINE                                                        ABSH1817
C     WILL FAIL IF A, MODIFIED BY THE ROUNDING ERRORS, IS SINGULAR      ABSH1818
C     OR ALMOST SINGULAR. SETS IFAIL = 0 IF SUCCESSFUL ELSE IFAIL =     ABSH1819
C     1.                                                                ABSH1820
C     1ST DECEMBER 1971                                                 ABSH1821
C                                                                       ABSH1822
      INTEGER ISAVE, IFAIL, IFAIL1, I, N, IA, ID, K, L, K1, K2,         ABSH1823
     *     ISTART, J, P01AAF                                            ABSH1824
      DOUBLE PRECISION SRNAME                                           ABSH1825
      REAL Y, D2, D1, X, EPS, A(IA,N), P(N)                             ABSH1826
      LOGICAL SW3                                                       ABSH1827
      DATA SRNAME /8H F03AFF /                                          ABSH1828
      DATA SW3 /.FALSE./                                                ABSH1829
      ISAVE = IFAIL                                                     ABSH1830
      IFAIL1 = 0                                                        ABSH1831
      DO 10 I=1,N                                                       ABSH1832
           CALL X03AAF(A(I,1), IA*N+1-I, A(I,1), IA*N+1-I, N,           ABSH1833
     *          IA, IA, 0.0, 0.0, Y, D2, SW3, IFAIL1)                   ABSH1834
           IF (Y.LE.0.0) GO TO 100                                      ABSH1835
           P(I) = 1.0/SQRT(Y)                                           ABSH1836
   10 CONTINUE                                                          ABSH1837
      D1 = 1.0                                                          ABSH1838
      ID = 0                                                            ABSH1839
      DO 90 K=1,N                                                       ABSH1840
           L = K                                                        ABSH1841
           X = 0.0                                                      ABSH1842
           K1 = K - 1                                                   ABSH1843
           K2 = K + 1                                                   ABSH1844
           ISTART = K                                                   ABSH1845
           DO 20 I=ISTART,N                                             ABSH1846
                CALL X03AAF(A(I,1), N*IA+1-I, A(1,K),                   ABSH1847
     *               (N-K+1)*IA, K1, IA, 1, -A(I,K), 0.0, Y,            ABSH1848
     *               D2, SW3, IFAIL1)                                   ABSH1849
                A(I,K) = -Y                                             ABSH1850
                Y = ABS(Y*P(I))                                         ABSH1851
                IF (Y.LE.X) GO TO 20                                    ABSH1852
                X = Y                                                   ABSH1853
                L = I                                                   ABSH1854
   20      CONTINUE                                                     ABSH1855
           IF (L.EQ.K) GO TO 40                                         ABSH1856
           D1 = -D1                                                     ABSH1857
           DO 30 J=1,N                                                  ABSH1858
                Y = A(K,J)                                              ABSH1859
                A(K,J) = A(L,J)                                         ABSH1860
                A(L,J) = Y                                              ABSH1861
   30      CONTINUE                                                     ABSH1862
           P(L) = P(K)                                                  ABSH1863
   40      P(K) = L                                                     ABSH1864
           D1 = D1*A(K,K)                                               ABSH1865
           IF (X.LT.8.0*EPS) GO TO 100                                  ABSH1866
   50      IF (ABS(D1).LT.1.0) GO TO 60                                 ABSH1867
           D1 = D1*0.0625                                               ABSH1868
           ID = ID + 4                                                  ABSH1869
           GO TO 50                                                     ABSH1870
   60      IF (ABS(D1).GE.0.0625) GO TO 70                              ABSH1871
           D1 = D1*16.0                                                 ABSH1872
           ID = ID - 4                                                  ABSH1873
           GO TO 60                                                     ABSH1874
   70      X = -1.0/A(K,K)                                              ABSH1875
           IF (K2.GT.N) GO TO 90                                        ABSH1876
           DO 80 J=K2,N                                                 ABSH1877
                CALL X03AAF(A(K,1), N*IA+1-K, A(1,J),                   ABSH1878
     *               (N-J+1)*IA, K1, IA, 1, -A(K,J), 0.0, Y,            ABSH1879
     *               D2, SW3, IFAIL1)                                   ABSH1880
                A(K,J) = X*Y                                            ABSH1881
   80      CONTINUE                                                     ABSH1882
   90 CONTINUE                                                          ABSH1883
      IFAIL = 0                                                         ABSH1884
      RETURN                                                            ABSH1885
  100 IFAIL = P01AAF(ISAVE,1,SRNAME)                                    ABSH1886
      RETURN                                                            ABSH1887
      END                                                               ABSH1888
      SUBROUTINE F04AJF(N, IR, A, IA, P, B, IB)                         ABSH1889
C     MARK 2 RELEASE. NAG COPYRIGHT 1972                                ABSH1890
C     EDITED BY JOYCE CLARKE OXFORD OEG NUCLEAR PHYSICS 11TH SEP 1976   ABSH1891
C                   FORTRAN MACRO VERSION FDIA24.TEC                    ABSH1892
C     MARK 4 REVISED.                                                   ABSH1893
C     MARK 4.5 REVISED                                                  ABSH1894
C     SW3 SWITCH REVISION : R. WIELAND // ORNL //   JAN., 1981.         ABSH1895
C     .FALSE. GIVES NO EXTENDED PRECISION                               ABSH1896
C     .TRUE. GIVES EXTENDED PRECISION IF CONDITIONS COMMENTED           ABSH1897
C     ON IN X03AAZ ARE FULFILLED.                                       ABSH1898
C     UNSYMSOL                                                          ABSH1899
C     SOLVES AX=B, WHERE A IS AN UNSYMMETRIC MATRIX AND B IS AN         ABSH1900
C     N*IR                                                              ABSH1901
C     MATRIX OF IR RIGHT-HAND SIDES. THE SUBROUTINE F04AJF MUST BY      ABSH1902
C     PRECEDED BY F03AFF IN WHICH L AND U ARE PRODUCED IN A(I,J),       ABSH1903
C     FROM A, AND THE RECORD OF THE INTERCHANGES IS PRODUCED IN         ABSH1904
C     P(I). AX=B IS SOLVED IN THREE STEPS, INTERCHANGE THE              ABSH1905
C     ELEMENTS OF B, LY=B AND UX=Y. THE MATRICES Y AND THEN X ARE       ABSH1906
C     OVERWRITTEN ON B.                                                 ABSH1907
C     1ST AUGUST 1971                                                   ABSH1908
C                                                                       ABSH1909
      INTEGER IFAIL1, I, N, I1, K, IR, IA, IB, II                       ABSH1910
      REAL X, D1, D2, A(IA,N), P(N), B(IB,IR)                           ABSH1911
      LOGICAL SW3                                                       ABSH1912
      DATA SW3 /.FALSE./                                                ABSH1913
      IFAIL1 = 0                                                        ABSH1914
C     INTERCHANGING OF ELEMENTS OF B                                    ABSH1915
      DO 20 I=1,N                                                       ABSH1916
           I1 = P(I) + 0.5                                              ABSH1917
           IF (I1.EQ.I) GO TO 20                                        ABSH1918
           DO 10 K=1,IR                                                 ABSH1919
                X = B(I,K)                                              ABSH1920
                B(I,K) = B(I1,K)                                        ABSH1921
                B(I1,K) = X                                             ABSH1922
   10      CONTINUE                                                     ABSH1923
   20 CONTINUE                                                          ABSH1924
      DO 50 K=1,IR                                                      ABSH1925
C     SOLUTION OF LY= B                                                 ABSH1926
           DO 30 I=1,N                                                  ABSH1927
                I1 = I - 1                                              ABSH1928
                CALL X03AAF(A(I,1), N*IA-I+1, B(1,K),                   ABSH1929
     *               (IR-K+1)*IB, I1, IA, 1, B(I,K), 0.0, D1,           ABSH1930
     *               D2, SW3, IFAIL1)                                   ABSH1931
                B(I,K) = -D1/A(I,I)                                     ABSH1932
   30      CONTINUE                                                     ABSH1933
C     SOLUTION OF UX= Y                                                 ABSH1934
           B(N,K) = -B(N,K)                                             ABSH1935
           IF (N.EQ.1) GO TO 50                                         ABSH1936
           DO 40 II=2,N                                                 ABSH1937
                I = N - II + 1                                          ABSH1938
                I1 = I + 1                                              ABSH1939
                CALL X03AAF(A(I,I1), (N-I)*IA-I+1, B(I1,K),             ABSH1940
     *               (IR-K+1)*IB-I, N-I, IA, 1, B(I,K), 0.0,            ABSH1941
     *               D1, D2, SW3, IFAIL1)                               ABSH1942
                B(I,K) = -D1                                            ABSH1943
   40      CONTINUE                                                     ABSH1944
   50 CONTINUE                                                          ABSH1945
      RETURN                                                            ABSH1946
      END                                                               ABSH1947
      INTEGER FUNCTION P01AAF(IFAIL, ERROR, SRNAME)                     ABSH1948
C     MARK 1 RELEASE.  NAG COPYRIGHT 1971                               ABSH1949
C     MARK 3 REVISED                                                    ABSH1950
C     MARK 4A REVISED, IER-45                                           ABSH1951
C     MARK 4.5 REVISED                                                  ABSH1952
C     MARK 7 REVISED (DEC 1978) .... (APR 1979)                         ABSH1953
C     NOUT = NTTY :: REVISED BY R. WIELAND // ORNL // JAN., 1981.       ABSH1954
C     RETURNS THE VALUE OF ERROR OR TERMINATES THE PROGRAM.             ABSH1955
C     IF A HARD FAILURE OCCURS, THIS ROUTINE CALLS A FORTRAN AUXILIARY  ABSH1956
C     ROUTINE P01AAZ WHICH GIVES A TRACE, A FAILURE MESSAGE AND HALTS   ABSH1957
C     THE PROGRAM                                                       ABSH1958
      INTEGER ERROR, IFAIL, NOUT                                        ABSH1959
      DOUBLE PRECISION SRNAME                                           ABSH1960
      COMMON /NIO/ NIN, NOUT, NPLOT, NDEBUG, NTTY                       ABSH1961
C     TEST IF NO ERROR DETECTED                                         ABSH1962
      IF (ERROR.EQ.0) GO TO 20                                          ABSH1963
C     TEST FOR SOFT FAILURE                                             ABSH1964
      IF (MOD(IFAIL,10).EQ.1) GO TO 10                                  ABSH1965
C     HARD FAILURE                                                      ABSH1966
      WRITE (NTTY,99999) SRNAME, ERROR                                  ABSH1967
C     STOPPING MECHANISM MAY ALSO DIFFER                                ABSH1968
C      CALL P01AAZ (X)                                                  ABSH1969
C     STOP                                                              ABSH1970
C     SOFT FAIL                                                         ABSH1971
C     TEST IF ERROR MESSAGES SUPPRESSED                                 ABSH1972
   10 IF (MOD(IFAIL/10,10).EQ.0) GO TO 20                               ABSH1973
      WRITE (NTTY,99999) SRNAME, ERROR                                  ABSH1974
   20 P01AAF = ERROR                                                    ABSH1975
      RETURN                                                            ABSH1976
99999 FORMAT (1H0, 38HERROR DETECTED BY NAG LIBRARY ROUTINE ,           ABSH1977
     *     A8, 11H - IFAIL = , I5//)                                    ABSH1978
      END                                                               ABSH1979
C     AUTO EDIT 20 SEP 76                                               ABSH1980
C     AUTO EDIT 20 SEP 76                                               ABSH1981
      SUBROUTINE X03AAF(A, ISIZEA, B, ISIZEB, N, ISTEPA,                ABSH1982
     *     ISTEPB, C1, C2, D1, D2, SW, IFAIL)                           ABSH1983
C     NAG COPYRIGHT 1975                                                ABSH1984
C     EDITED BY JOYCE CLARKE OXFORD OEG NUCLEAR PHYSICS 11TH SEP 1976   ABSH1985
C                   FORTRAN MACRO VERSION FDIA24.TEC                    ABSH1986
C     MARK 4.5 RELEASE                                                  ABSH1987
C                                                                       ABSH1988
C     CALCULATES THE VALUE OF A SCALAR PRODUCT USING BASIC              ABSH1989
C     OR ADDITIONAL PRECISION AND ADDS IT TO A BASIC OR ADDITIONAL      ABSH1990
C     PRECISION INITIAL VALUE.                                          ABSH1991
C     X03AAF CALLS X03AAZ WHICH MAY BE IN ASSEMBLY LANGUAGE.            ABSH1992
C                                                                       ABSH1993
      INTEGER P01AAF, ISAVE, ISIZEA, ISIZEB, ISTEPA, ISTEPB,            ABSH1994
     *     IFAIL, IS, IT, N, I                                          ABSH1995
      DOUBLE PRECISION SUM, SRNAME                                      ABSH1996
      REAL A(ISIZEA), B(ISIZEB), C1, C2, D1, D2, X                      ABSH1997
      LOGICAL SW                                                        ABSH1998
      DATA SRNAME /8H X03AAF /                                          ABSH1999
      ISAVE = IFAIL                                                     ABSH2000
      IFAIL = 0                                                         ABSH2001
      IF (ISTEPA.GT.0 .AND. ISTEPB.GT.0) GO TO 10                       ABSH2002
      IFAIL = P01AAF(ISAVE,1,SRNAME)                                    ABSH2003
      RETURN                                                            ABSH2004
   10 IS = 1 - ISTEPA                                                   ABSH2005
      IT = 1 - ISTEPB                                                   ABSH2006
      IF (SW) GO TO 40                                                  ABSH2007
      X = 0.0                                                           ABSH2008
      IF (N.LT.1) GO TO 30                                              ABSH2009
      DO 20 I=1,N                                                       ABSH2010
           IS = IS + ISTEPA                                             ABSH2011
           IT = IT + ISTEPB                                             ABSH2012
           X = X + A(IS)*B(IT)                                          ABSH2013
   20 CONTINUE                                                          ABSH2014
   30 D1 = X + (C1+C2)                                                  ABSH2015
      D2 = 0.0                                                          ABSH2016
      RETURN                                                            ABSH2017
   40 SUM = 0.0D0                                                       ABSH2018
      IF (N.LT.1) GO TO 60                                              ABSH2019
      DO 50 I=1,N                                                       ABSH2020
           IS = IS + ISTEPA                                             ABSH2021
           IT = IT + ISTEPB                                             ABSH2022
           SUM = SUM + DBLE(A(IS))*B(IT)                                ABSH2023
   50 CONTINUE                                                          ABSH2024
   60 SUM = SUM + (DBLE(C1)+C2)                                         ABSH2025
      CALL X03AAZ(SUM, D1, D2)                                          ABSH2026
      RETURN                                                            ABSH2027
      END                                                               ABSH2028
      SUBROUTINE X03AAZ(DP, D1, D2)                                     ABSH2029
      DOUBLE PRECISION DP                                               ABSH2030
      REAL D1, D2                                                       ABSH2031
      D1 = DP                                                           ABSH2032
      D2 = DP - D1                                                      ABSH2033
      RETURN                                                            ABSH2034
C     RETURN                                                            ABSH2035
      END                                                               ABSH2036
      REAL FUNCTION SPMPAR(I)                                           ABSH2037
      INTEGER I                                                         ABSH2038
C     **********                                                        ABSH2039
C                                                                       ABSH2040
C     FUNCTION SPMPAR                                                   ABSH2041
C                                                                       ABSH2042
C     THIS FUNCTION PROVIDES SINGLE PRECISION MACHINE PARAMETERS        ABSH2043
C     WHEN THE APPROPRIATE SET OF DATA STATEMENTS IS ACTIVATED (BY      ABSH2044
C     REMOVING THE C FROM COLUMN 1) AND ALL OTHER DATA STATEMENTS ARE   ABSH2045
C     RENDERED INACTIVE. MOST OF THE PARAMETER VALUES WERE OBTAINED     ABSH2046
C     FROM THE CORRESPONDING BELL LABORATORIES PORT LIBRARY FUNCTION.   ABSH2047
C                                                                       ABSH2048
C     THE FUNCTION STATEMENT IS                                         ABSH2049
C                                                                       ABSH2050
C       REAL FUNCTION SPMPAR(I)                                         ABSH2051
C                                                                       ABSH2052
C     WHERE                                                             ABSH2053
C                                                                       ABSH2054
C       I IS AN INTEGER INPUT VARIABLE SET TO 1, 2, OR 3 WHICH          ABSH2055
C         SELECTS THE DESIRED MACHINE PARAMETER. IF THE MACHINE HAS     ABSH2056
C         T BASE B DIGITS AND ITS SMALLEST AND LARGEST EXPONENTS ARE    ABSH2057
C         EMIN AND EMAX, RESPECTIVELY, THEN THESE PARAMETERS ARE        ABSH2058
C                                                                       ABSH2059
C         SPMPAR(1) = B**(1 - T), THE MACHINE PRECISION,                ABSH2060
C                                                                       ABSH2061
C         SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,            ABSH2062
C                                                                       ABSH2063
C         SPMPAR(3) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.     ABSH2064
C                                                                       ABSH2065
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.         ABSH2066
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE             ABSH2067
C                                                                       ABSH2068
C     **********                                                        ABSH2069
      INTEGER MCHEPS(1)                                                 ABSH2070
      INTEGER MINMAG(1)                                                 ABSH2071
      INTEGER MAXMAG(1)                                                 ABSH2072
      REAL RMACH(3)                                                     ABSH2073
      EQUIVALENCE (RMACH(1),MCHEPS(1))                                  ABSH2074
      EQUIVALENCE (RMACH(2),MINMAG(1))                                  ABSH2075
      EQUIVALENCE (RMACH(3),MAXMAG(1))                                  ABSH2076
C                                                                       ABSH2077
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,                     ABSH2078
C     THE AMDAHL 470/V6, THE ICL 2900, THE ITEL AS/6,                   ABSH2079
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.                  ABSH2080
C                                                                       ABSH2081
C     DATA RMACH(1) / Z3C100000 /                                       ABSH2082
C     DATA RMACH(2) / Z00100000 /                                       ABSH2083
C     DATA RMACH(3) / Z7FFFFFFF /                                       ABSH2084
C                                                                       ABSH2085
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.              ABSH2086
C                                                                       ABSH2087
C     DATA RMACH(1) / O716400000000 /                                   ABSH2088
C     DATA RMACH(2) / O402400000000 /                                   ABSH2089
C     DATA RMACH(3) / O376777777777 /                                   ABSH2090
C                                                                       ABSH2091
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.                   ABSH2092
C                                                                       ABSH2093
C     DATA RMACH(1) / 16414000000000000000B /                           ABSH2094
C     DATA RMACH(2) / 00014000000000000000B /                           ABSH2095
C     DATA RMACH(3) / 37767777777777777777B /                           ABSH2096
C                                                                       ABSH2097
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).            ABSH2098
C                                                                       ABSH2099
C     DATA RMACH(1) / "147400000000 /                                   ABSH2100
C     DATA RMACH(2) / "000400000000 /                                   ABSH2101
C     DATA RMACH(3) / "377777777777 /                                   ABSH2102
C                                                                       ABSH2103
C     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING               ABSH2104
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).                 ABSH2105
C                                                                       ABSH2106
      DATA MCHEPS(1) /  889192448 /                                     ABSH2107
      DATA MINMAG(1) /    8388608 /                                     ABSH2108
      DATA MAXMAG(1) / 2147483647 /                                     ABSH2109
                                                                        ABSH2110
C      DATA RMACH(1) / O06500000000 /                                    ABSH2111
C      DATA RMACH(2) / O00040000000 /                                    ABSH2112
C      DATA RMACH(3) / O17777777777 /                                    ABSH2113
C                                                                       ABSH2114
C     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING               ABSH2115
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).                 ABSH2116
C                                                                       ABSH2117
C     DATA MCHEPS(1),MCHEPS(2) / 13568,     0 /                         ABSH2118
C     DATA MINMAG(1),MINMAG(2) /   128,     0 /                         ABSH2119
C     DATA MAXMAG(1),MAXMAG(2) / 32767,    -1 /                         ABSH2120
C                                                                       ABSH2121
C     DATA MCHEPS(1),MCHEPS(2) / O032400, O000000 /                     ABSH2122
C     DATA MINMAG(1),MINMAG(2) / O000200, O000000 /                     ABSH2123
C     DATA MAXMAG(1),MAXMAG(2) / O077777, O177777 /                     ABSH2124
C                                                                       ABSH2125
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS.       ABSH2126
C                                                                       ABSH2127
C     DATA RMACH(1) / O1301000000000000 /                               ABSH2128
C     DATA RMACH(2) / O1771000000000000 /                               ABSH2129
C     DATA RMACH(3) / O0777777777777777 /                               ABSH2130
C                                                                       ABSH2131
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.                  ABSH2132
C                                                                       ABSH2133
C     DATA RMACH(1) / Z4EA800000 /                                      ABSH2134
C     DATA RMACH(2) / Z400800000 /                                      ABSH2135
C     DATA RMACH(3) / Z5FFFFFFFF /                                      ABSH2136
C                                                                       ABSH2137
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.                     ABSH2138
C                                                                       ABSH2139
C     DATA RMACH(1) / O147400000000 /                                   ABSH2140
C     DATA RMACH(2) / O000400000000 /                                   ABSH2141
C     DATA RMACH(3) / O377777777777 /                                   ABSH2142
C                                                                       ABSH2143
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.             ABSH2144
C                                                                       ABSH2145
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -      ABSH2146
C     STATIC RMACH(3)                                                   ABSH2147
C                                                                       ABSH2148
C     DATA MINMAG/20K,0/,MAXMAG/77777K,177777K/                         ABSH2149
C     DATA MCHEPS/36020K,0/                                             ABSH2150
C                                                                       ABSH2151
C     MACHINE CONSTANTS FOR THE HARRIS 220.                             ABSH2152
C                                                                       ABSH2153
C     DATA MCHEPS(1) / '20000000, '00000353 /                           ABSH2154
C     DATA MINMAG(1) / '20000000, '00000201 /                           ABSH2155
C     DATA MAXMAG(1) / '37777777, '00000177 /                           ABSH2156
C                                                                       ABSH2157
C     MACHINE CONSTANTS FOR THE CRAY-1.                                 ABSH2158
C                                                                       ABSH2159
C      DATA RMACH(1) / 0377224000000000000000B /                        ABSH2160
C      DATA RMACH(2) / 0200034000000000000000B /                        ABSH2161
C      DATA RMACH(3) / 0577777777777777777776B /                        ABSH2162
C                                                                       ABSH2163
C     MACHINE CONSTANTS FOR THE PRIME 400.                              ABSH2164
C                                                                       ABSH2165
C     DATA MCHEPS(1) / :10000000153 /                                   ABSH2166
C     DATA MINMAG(1) / :10000000000 /                                   ABSH2167
C     DATA MAXMAG(1) / :17777777777 /                                   ABSH2168
C                                                                       ABSH2169
      SPMPAR = RMACH(I)                                                 ABSH2170
      RETURN                                                            ABSH2171
C                                                                       ABSH2172
C     LAST CARD OF FUNCTION SPMPAR.                                     ABSH2173
C                                                                       ABSH2174
      END                                                               ABSH2175
C                                                                       ABSH2176
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     ABSH2177
C                                                                       ABSH2178
      INTEGER FUNCTION I1MACH(I)                                        ABSH2179
C                                                                       ABSH2180
C  I/O UNIT NUMBERS.                                                    ABSH2181
C                                                                       ABSH2182
C    I1MACH( 1) = THE STANDARD INPUT UNIT.                              ABSH2183
C                                                                       ABSH2184
C    I1MACH( 2) = THE STANDARD OUTPUT UNIT.                             ABSH2185
C                                                                       ABSH2186
C    I1MACH( 3) = THE STANDARD PUNCH UNIT.                              ABSH2187
C                                                                       ABSH2188
C    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.                      ABSH2189
C                                                                       ABSH2190
C  WORDS.                                                               ABSH2191
C                                                                       ABSH2192
C    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.          ABSH2193
C                                                                       ABSH2194
C    I1MACH( 6) = THE NUMBER OF CHARACTERS PER INTEGER STORAGE UNIT.    ABSH2195
C                                                                       ABSH2196
C  INTEGERS.                                                            ABSH2197
C                                                                       ABSH2198
C    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM        ABSH2199
C                                                                       ABSH2200
C               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )          ABSH2201
C                                                                       ABSH2202
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.               ABSH2203
C                                                                       ABSH2204
C    I1MACH( 7) = A, THE BASE.                                          ABSH2205
C                                                                       ABSH2206
C    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.                       ABSH2207
C                                                                       ABSH2208
C    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.                      ABSH2209
C                                                                       ABSH2210
C  FLOATING-POINT NUMBERS.                                              ABSH2211
C                                                                       ABSH2212
C    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,      ABSH2213
C    BASE-B FORM                                                        ABSH2214
C                                                                       ABSH2215
C               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )            ABSH2216
C                                                                       ABSH2217
C               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,                 ABSH2218
C               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.                 ABSH2219
C                                                                       ABSH2220
C    I1MACH(10) = B, THE BASE.                                          ABSH2221
C                                                                       ABSH2222
C  SINGLE-PRECISION                                                     ABSH2223
C                                                                       ABSH2224
C    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.                       ABSH2225
C                                                                       ABSH2226
C    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.                        ABSH2227
C                                                                       ABSH2228
C    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.                         ABSH2229
C                                                                       ABSH2230
C  DOUBLE-PRECISION                                                     ABSH2231
C                                                                       ABSH2232
C    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.                       ABSH2233
C                                                                       ABSH2234
C    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.                        ABSH2235
C                                                                       ABSH2236
C    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.                         ABSH2237
C                                                                       ABSH2238
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,                 ABSH2239
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY            ABSH2240
C  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF                   ABSH2241
C  I1MACH(1) - I1MACH(4) SHOULD BE CHECKED FOR CONSISTENCY              ABSH2242
C  WITH THE LOCAL OPERATING SYSTEM.                                     ABSH2243
C                                                                       ABSH2244
      INTEGER IMACH(16), OUTPUT                                         ABSH2245
C                                                                       ABSH2246
      EQUIVALENCE (IMACH(4),OUTPUT)                                     ABSH2247
C                                                                       ABSH2248
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.                  ABSH2249
C                                                                       ABSH2250
C     DATA IMACH( 1) /    7 /                                           ABSH2251
C     DATA IMACH( 2) /    2 /                                           ABSH2252
C     DATA IMACH( 3) /    2 /                                           ABSH2253
C     DATA IMACH( 4) /    2 /                                           ABSH2254
C     DATA IMACH( 5) /   36 /                                           ABSH2255
C     DATA IMACH( 6) /    4 /                                           ABSH2256
C     DATA IMACH( 7) /    2 /                                           ABSH2257
C     DATA IMACH( 8) /   33 /                                           ABSH2258
C     DATA IMACH( 9) / Z1FFFFFFFF /                                     ABSH2259
C     DATA IMACH(10) /    2 /                                           ABSH2260
C     DATA IMACH(11) /   24 /                                           ABSH2261
C     DATA IMACH(12) / -256 /                                           ABSH2262
C     DATA IMACH(13) /  255 /                                           ABSH2263
C     DATA IMACH(14) /   60 /                                           ABSH2264
C     DATA IMACH(15) / -256 /                                           ABSH2265
C     DATA IMACH(16) /  255 /                                           ABSH2266
C                                                                       ABSH2267
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.                  ABSH2268
C                                                                       ABSH2269
C     DATA IMACH( 1) /   5 /                                            ABSH2270
C     DATA IMACH( 2) /   6 /                                            ABSH2271
C     DATA IMACH( 3) /   7 /                                            ABSH2272
C     DATA IMACH( 4) /   6 /                                            ABSH2273
C     DATA IMACH( 5) /  48 /                                            ABSH2274
C     DATA IMACH( 6) /   6 /                                            ABSH2275
C     DATA IMACH( 7) /   2 /                                            ABSH2276
C     DATA IMACH( 8) /  39 /                                            ABSH2277
C     DATA IMACH( 9) / O0007777777777777 /                              ABSH2278
C     DATA IMACH(10) /   8 /                                            ABSH2279
C     DATA IMACH(11) /  13 /                                            ABSH2280
C     DATA IMACH(12) / -50 /                                            ABSH2281
C     DATA IMACH(13) /  76 /                                            ABSH2282
C     DATA IMACH(14) /  26 /                                            ABSH2283
C     DATA IMACH(15) / -50 /                                            ABSH2284
C     DATA IMACH(16) /  76 /                                            ABSH2285
C                                                                       ABSH2286
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.            ABSH2287
C                                                                       ABSH2288
C     DATA IMACH( 1) /   5 /                                            ABSH2289
C     DATA IMACH( 2) /   6 /                                            ABSH2290
C     DATA IMACH( 3) /   7 /                                            ABSH2291
C     DATA IMACH( 4) /   6 /                                            ABSH2292
C     DATA IMACH( 5) /  48 /                                            ABSH2293
C     DATA IMACH( 6) /   6 /                                            ABSH2294
C     DATA IMACH( 7) /   2 /                                            ABSH2295
C     DATA IMACH( 8) /  39 /                                            ABSH2296
C     DATA IMACH( 9) / O0007777777777777 /                              ABSH2297
C     DATA IMACH(10) /   8 /                                            ABSH2298
C     DATA IMACH(11) /  13 /                                            ABSH2299
C     DATA IMACH(12) / -50 /                                            ABSH2300
C     DATA IMACH(13) /  76 /                                            ABSH2301
C     DATA IMACH(14) /  26 /                                            ABSH2302
C     DATA IMACH(15) / -32754 /                                         ABSH2303
C     DATA IMACH(16) /  32780 /                                         ABSH2304
C                                                                       ABSH2305
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.                   ABSH2306
C                                                                       ABSH2307
C     DATA IMACH( 1) /    5 /                                           ABSH2308
C     DATA IMACH( 2) /    6 /                                           ABSH2309
C     DATA IMACH( 3) /    7 /                                           ABSH2310
C     DATA IMACH( 4) /    6 /                                           ABSH2311
C     DATA IMACH( 5) /   60 /                                           ABSH2312
C     DATA IMACH( 6) /   10 /                                           ABSH2313
C     DATA IMACH( 7) /    2 /                                           ABSH2314
C     DATA IMACH( 8) /   48 /                                           ABSH2315
C     DATA IMACH( 9) / 00007777777777777777B /                          ABSH2316
C     DATA IMACH(10) /    2 /                                           ABSH2317
C     DATA IMACH(11) /   48 /                                           ABSH2318
C     DATA IMACH(12) / -974 /                                           ABSH2319
C     DATA IMACH(13) / 1070 /                                           ABSH2320
C     DATA IMACH(14) /   96 /                                           ABSH2321
C     DATA IMACH(15) / -927 /                                           ABSH2322
C     DATA IMACH(16) / 1070 /                                           ABSH2323
C                                                                       ABSH2324
C     MACHINE CONSTANTS FOR THE CRAY 1                                  ABSH2325
C                                                                       ABSH2326
C      DATA IMACH( 1) /   100 /                                         ABSH2327
C      DATA IMACH( 2) /   5 /                                           ABSH2328
C      DATA IMACH( 3) /   102 /                                         ABSH2329
C      DATA IMACH( 4) /   101 /                                         ABSH2330
C      DATA IMACH( 5) /    64 /                                         ABSH2331
C      DATA IMACH( 6) /     8 /                                         ABSH2332
C      DATA IMACH( 7) /     2 /                                         ABSH2333
C      DATA IMACH( 8) /    63 /                                         ABSH2334
C      DATA IMACH( 9) /  777777777777777777777B /                       ABSH2335
C      DATA IMACH(10) /     2 /                                         ABSH2336
C      DATA IMACH(11) /    48 /                                         ABSH2337
C      DATA IMACH(12) / -8192 /                                         ABSH2338
C      DATA IMACH(13) /  8191 /                                         ABSH2339
C      DATA IMACH(14) /    96 /                                         ABSH2340
C      DATA IMACH(15) / -8192 /                                         ABSH2341
C      DATA IMACH(16) /  8191 /                                         ABSH2342
C                                                                       ABSH2343
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200              ABSH2344
C                                                                       ABSH2345
C     DATA IMACH( 1) /   11 /                                           ABSH2346
C     DATA IMACH( 2) /   12 /                                           ABSH2347
C     DATA IMACH( 3) /    8 /                                           ABSH2348
C     DATA IMACH( 4) /   10 /                                           ABSH2349
C     DATA IMACH( 5) /   16 /                                           ABSH2350
C     DATA IMACH( 6) /    2 /                                           ABSH2351
C     DATA IMACH( 7) /    2 /                                           ABSH2352
C     DATA IMACH( 8) /   15 /                                           ABSH2353
C     DATA IMACH( 9) /32767 /                                           ABSH2354
C     DATA IMACH(10) /   16 /                                           ABSH2355
C     DATA IMACH(11) /    6 /                                           ABSH2356
C     DATA IMACH(12) /  -64 /                                           ABSH2357
C     DATA IMACH(13) /   63 /                                           ABSH2358
C     DATA IMACH(14) /   14 /                                           ABSH2359
C     DATA IMACH(15) /  -64 /                                           ABSH2360
C     DATA IMACH(16) /   63 /                                           ABSH2361
C                                                                       ABSH2362
C     MACHINE CONSTANTS FOR THE HARRIS 220                              ABSH2363
C                                                                       ABSH2364
C     DATA IMACH( 1) /       5 /                                        ABSH2365
C     DATA IMACH( 2) /       6 /                                        ABSH2366
C     DATA IMACH( 3) /       0 /                                        ABSH2367
C     DATA IMACH( 4) /       6 /                                        ABSH2368
C     DATA IMACH( 5) /      24 /                                        ABSH2369
C     DATA IMACH( 6) /       3 /                                        ABSH2370
C     DATA IMACH( 7) /       2 /                                        ABSH2371
C     DATA IMACH( 8) /      23 /                                        ABSH2372
C     DATA IMACH( 9) / 8388607 /                                        ABSH2373
C     DATA IMACH(10) /       2 /                                        ABSH2374
C     DATA IMACH(11) /      23 /                                        ABSH2375
C     DATA IMACH(12) /    -127 /                                        ABSH2376
C     DATA IMACH(13) /     127 /                                        ABSH2377
C     DATA IMACH(14) /      38 /                                        ABSH2378
C     DATA IMACH(15) /    -127 /                                        ABSH2379
C     DATA IMACH(16) /     127 /                                        ABSH2380
C                                                                       ABSH2381
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.              ABSH2382
C                                                                       ABSH2383
C     DATA IMACH( 1) /    5 /                                           ABSH2384
C     DATA IMACH( 2) /    6 /                                           ABSH2385
C     DATA IMACH( 3) /   43 /                                           ABSH2386
C     DATA IMACH( 4) /    6 /                                           ABSH2387
C     DATA IMACH( 5) /   36 /                                           ABSH2388
C     DATA IMACH( 6) /    6 /                                           ABSH2389
C     DATA IMACH( 7) /    2 /                                           ABSH2390
C     DATA IMACH( 8) /   35 /                                           ABSH2391
C     DATA IMACH( 9) / O377777777777 /                                  ABSH2392
C     DATA IMACH(10) /    2 /                                           ABSH2393
C     DATA IMACH(11) /   27 /                                           ABSH2394
C     DATA IMACH(12) / -127 /                                           ABSH2395
C     DATA IMACH(13) /  127 /                                           ABSH2396
C     DATA IMACH(14) /   63 /                                           ABSH2397
C     DATA IMACH(15) / -127 /                                           ABSH2398
C     DATA IMACH(16) /  127 /                                           ABSH2399
C                                                                       ABSH2400
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,                     ABSH2401
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.                  ABSH2402
C                                                                       ABSH2403
C     DATA IMACH( 1) /   5 /                                            ABSH2404
C     DATA IMACH( 2) /   6 /                                            ABSH2405
C     DATA IMACH( 3) /   7 /                                            ABSH2406
C     DATA IMACH( 4) /   6 /                                            ABSH2407
C     DATA IMACH( 5) /  32 /                                            ABSH2408
C     DATA IMACH( 6) /   4 /                                            ABSH2409
C     DATA IMACH( 7) /   2 /                                            ABSH2410
C     DATA IMACH( 8) /  31 /                                            ABSH2411
C     DATA IMACH( 9) / Z7FFFFFFF /                                      ABSH2412
C     DATA IMACH(10) /  16 /                                            ABSH2413
C     DATA IMACH(11) /   6 /                                            ABSH2414
C     DATA IMACH(12) / -64 /                                            ABSH2415
C     DATA IMACH(13) /  63 /                                            ABSH2416
C     DATA IMACH(14) /  14 /                                            ABSH2417
C     DATA IMACH(15) / -64 /                                            ABSH2418
C     DATA IMACH(16) /  63 /                                            ABSH2419
C                                                                       ABSH2420
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).                  ABSH2421
C                                                                       ABSH2422
C     DATA IMACH(1) /5/                                                 ABSH2423
C     DATA IMACH(2) /6/                                                 ABSH2424
C     DATA IMACH(3) /5/                                                 ABSH2425
C     DATA IMACH(4) /6/                                                 ABSH2426
C     DATA IMACH(5) /36/                                                ABSH2427
C     DATA IMACH(6) /5/                                                 ABSH2428
C     DATA IMACH(7) /2/                                                 ABSH2429
C     DATA IMACH(8) /35/                                                ABSH2430
C     DATA IMACH( 9) / "377777777777 /                                  ABSH2431
C     DATA IMACH(10) /2/                                                ABSH2432
C     DATA IMACH(11) /27/                                               ABSH2433
C     DATA IMACH(12) /-128/                                             ABSH2434
C     DATA IMACH(13) /127/                                              ABSH2435
C     DATA IMACH(14) /54/                                               ABSH2436
C     DATA IMACH(15) /-101/                                             ABSH2437
C     DATA IMACH(16) /127/                                              ABSH2438
C                                                                       ABSH2439
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).                  ABSH2440
C                                                                       ABSH2441
C     DATA IMACH( 1) /    5 /                                           ABSH2442
C     DATA IMACH( 2) /    6 /                                           ABSH2443
C     DATA IMACH( 3) /    5 /                                           ABSH2444
C     DATA IMACH( 4) /    6 /                                           ABSH2445
C     DATA IMACH( 5) /   36 /                                           ABSH2446
C     DATA IMACH( 6) /    5 /                                           ABSH2447
C     DATA IMACH( 7) /    2 /                                           ABSH2448
C     DATA IMACH( 8) /   35 /                                           ABSH2449
C     DATA IMACH( 9) / "377777777777 /                                  ABSH2450
C     DATA IMACH(10) /    2 /                                           ABSH2451
C     DATA IMACH(11) /   27 /                                           ABSH2452
C     DATA IMACH(12) / -128 /                                           ABSH2453
C     DATA IMACH(13) /  127 /                                           ABSH2454
C     DATA IMACH(14) /   62 /                                           ABSH2455
C     DATA IMACH(15) / -128 /                                           ABSH2456
C     DATA IMACH(16) /  127 /                                           ABSH2457
C                                                                       ABSH2458
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN'S SUPPORTING                 ABSH2459
C     32-BIT INTEGER ARITHMETIC.                                        ABSH2460
C                                                                       ABSH2461
      DATA IMACH( 1) /    5 /                                           ABSH2462
      DATA IMACH( 2) /    6 /                                           ABSH2463
      DATA IMACH( 3) /    5 /                                           ABSH2464
      DATA IMACH( 4) /    6 /                                           ABSH2465
      DATA IMACH( 5) /   32 /                                           ABSH2466
      DATA IMACH( 6) /    4 /                                           ABSH2467
      DATA IMACH( 7) /    2 /                                           ABSH2468
      DATA IMACH( 8) /   31 /                                           ABSH2469
      DATA IMACH( 9) / 2147483647 /                                     ABSH2470
      DATA IMACH(10) /    2 /                                           ABSH2471
      DATA IMACH(11) /   24 /                                           ABSH2472
      DATA IMACH(12) / -127 /                                           ABSH2473
      DATA IMACH(13) /  127 /                                           ABSH2474
      DATA IMACH(14) /   56 /                                           ABSH2475
      DATA IMACH(15) / -127 /                                           ABSH2476
      DATA IMACH(16) /  127 /                                           ABSH2477
C                                                                       ABSH2478
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN'S SUPPORTING                 ABSH2479
C     16-BIT INTEGER ARITHMETIC.                                        ABSH2480
C                                                                       ABSH2481
C     DATA IMACH( 1) /    5 /                                           ABSH2482
C     DATA IMACH( 2) /    6 /                                           ABSH2483
C     DATA IMACH( 3) /    5 /                                           ABSH2484
C     DATA IMACH( 4) /    6 /                                           ABSH2485
C     DATA IMACH( 5) /   16 /                                           ABSH2486
C     DATA IMACH( 6) /    2 /                                           ABSH2487
C     DATA IMACH( 7) /    2 /                                           ABSH2488
C     DATA IMACH( 8) /   15 /                                           ABSH2489
C     DATA IMACH( 9) / 32767 /                                          ABSH2490
C     DATA IMACH(10) /    2 /                                           ABSH2491
C     DATA IMACH(11) /   24 /                                           ABSH2492
C     DATA IMACH(12) / -127 /                                           ABSH2493
C     DATA IMACH(13) /  127 /                                           ABSH2494
C     DATA IMACH(14) /   56 /                                           ABSH2495
C     DATA IMACH(15) / -127 /                                           ABSH2496
C     DATA IMACH(16) /  127 /                                           ABSH2497
C                                                                       ABSH2498
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.                     ABSH2499
C                                                                       ABSH2500
C     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7            ABSH2501
C     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM.                   ABSH2502
C     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1.                   ABSH2503
C                                                                       ABSH2504
C     DATA IMACH( 1) /    5 /                                           ABSH2505
C     DATA IMACH( 2) /    6 /                                           ABSH2506
C     DATA IMACH( 3) /    7 /                                           ABSH2507
C     DATA IMACH( 4) /    6 /                                           ABSH2508
C     DATA IMACH( 5) /   36 /                                           ABSH2509
C     DATA IMACH( 6) /    6 /                                           ABSH2510
C     DATA IMACH( 7) /    2 /                                           ABSH2511
C     DATA IMACH( 8) /   35 /                                           ABSH2512
C     DATA IMACH( 9) / O377777777777 /                                  ABSH2513
C     DATA IMACH(10) /    2 /                                           ABSH2514
C     DATA IMACH(11) /   27 /                                           ABSH2515
C     DATA IMACH(12) / -128 /                                           ABSH2516
C     DATA IMACH(13) /  127 /                                           ABSH2517
C     DATA IMACH(14) /   60 /                                           ABSH2518
C     DATA IMACH(15) /-1024 /                                           ABSH2519
C     DATA IMACH(16) / 1023 /                                           ABSH2520
C                                                                       ABSH2521
      IF (I.LT.1 .OR. I.GT.16) GO TO 10                                 ABSH2522
C                                                                       ABSH2523
      I1MACH = IMACH(I)                                                 ABSH2524
      RETURN                                                            ABSH2525
C                                                                       ABSH2526
   10 WRITE (OUTPUT,99999)                                              ABSH2527
C                                                                       ABSH2528
C                                                                       ABSH2529
      STOP                                                              ABSH2530
C                                                                       ABSH2531
99999 FORMAT (39H1ERROR    1 IN I1MACH - I OUT OF BOUNDS)               ABSH2532
      END                                                               ABSH2533
