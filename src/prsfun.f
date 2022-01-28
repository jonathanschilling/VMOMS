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
