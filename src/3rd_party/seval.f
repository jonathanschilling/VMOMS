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
