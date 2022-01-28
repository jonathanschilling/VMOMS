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
