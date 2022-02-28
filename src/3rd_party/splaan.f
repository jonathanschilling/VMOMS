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
