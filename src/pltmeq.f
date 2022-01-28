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
