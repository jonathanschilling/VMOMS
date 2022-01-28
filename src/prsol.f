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
