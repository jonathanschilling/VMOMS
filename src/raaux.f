      SUBROUTINE RAAUX(XL, XR, RMATCH, PARAM)                           ABSH0892
C***********************************************************************ABSH0893
C*****RAAUX SETS THE MATCHING POINT FOR THE MHD MOMENTS.               *ABSH0894
C*****REFERENCES:                                                      *ABSH0895
C*****L.L.LAO,S.P.HIRSHMAN,R.M.WIELAND,ORNL/TM-7616 (1981).            *ABSH0896
C*****                                ,PHYS FLUIDS 24 (1981) AUG       *ABSH0897
C*****L.L.LAO,R.M.WIELAND,W.A.HOULBERG,S.P.HIRSHMAN,ORNL/TM-7871 (1981)*ABSH0898
C*****LAST REVISION: 6/81 L.L.LAO, R.M.WIELAND, AND W.A.HOULBERG ORNL. *ABSH0899
C*****CALCULATED PARAMETERS:                                           *ABSH0900
C*****XL-LEFT (LOWER) LIMIT ON ABSCISSA.                               *ABSH0901
C*****XR-RIGHT (UPPER) LIMIT ON ABSCISSA.                              *ABSH0902
C*****RMATCH-MATCHING POINT FOR SHOOTING METHOD (RMATCH=1.0).          *ABSH0903
C*****INPUT PARAMETERS:                                                *ABSH0904
C*****XL0-LOWER BOUND OF ABSCISSA.                                     *ABSH0905
C*****OTHER COMMENTS:                                                  *ABSH0906
C*****RAAUX IS REQUIRED BY D02AGF ODE SOLVER IN NAG LIBRARY.           *ABSH0907
C***********************************************************************ABSH0908
      PARAMETER MNMOMS=3,MNEQ=6                                         ABSH0909
      COMMON /LOCMEQ/ XL0, DX, R0, A0, BT0, P0, AI0, E0, E1,            ABSH0910
     *     D0, ARC, SLPE, SLPR0, BETAP0, KMHD, IUNDER                   ABSH0911
      DIMENSION PARAM(MNEQ)                                             ABSH0912
      XL = XL0                                                          ABSH0913
      XR = 1.0                                                          ABSH0914
      RMATCH = 1.0                                                      ABSH0915
      RETURN                                                            ABSH0916
      END                                                               ABSH0917
