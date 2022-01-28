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
