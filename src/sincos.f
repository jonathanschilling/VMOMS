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
