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
