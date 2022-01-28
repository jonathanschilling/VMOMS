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
