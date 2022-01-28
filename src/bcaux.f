      SUBROUTINE BCAUX(GL, GR, PARAM)                                   ABSH0424
C***********************************************************************ABSH0425
C*****BCAUX SETS UP THE BOUNDARY CONDITIONS FOR THE MHD MOMENTS.       *ABSH0426
C*****REFERENCES:                                                      *ABSH0427
C*****L.L.LAO,S.P.HIRSHMAN,R.M.WIELAND,ORNL/TM-7616 (1981).            *ABSH0428
C*****                                ,PHYS FLUIDS 24 (1981) AUG       *ABSH0429
C*****L.L.LAO,R.M.WIELAND,W.A.HOULBERG,S.P.HIRSHMAN,ORNL/TM-7871 (1981)*ABSH0430
C*****LAST REVISION: 6/81 L.L.LAO, R.M.WIELAND, AND W.A.HOULBERG ORNL. *ABSH0431
C*****CALCULATED PARAMETERS:                                           *ABSH0432
C*****GL(I)-VALUE OF Y(I) AT LEFT BOUNDARY.                            *ABSH0433
C*****GR(I)-VALUE OF Y(I) AT RIGHT BOUNDARY.                           *ABSH0434
C*****INPUT PARAMETERS:                                                *ABSH0435
C*****PARAM(I)-B.C. PARAMETER FOR FIRST ORDER ODE I.                   *ABSH0436
C*****XL0-LOWER BOUND OF ABSCISSA.                                     *ABSH0437
C*****SLPR0-ANALYTIC PARAMETER FOR LEFT B.C. AT X1.                    *ABSH0438
C*****SLPE-ANALYTIC PARAMETER FOR LEFT B.C. AT X1.                     *ABSH0439
C*****ARC-INVERSE ASPECT RATIO FOR OUTERMOST FLUX SURFACE.             *ABSH0440
C*****E0-MOMENTS REPRESENTATION OF ELONGATION AT X=1.                  *ABSH0441
C*****D0-MOMENTS REPRESENTATION OF TRIANGULARITY AT X=1.               *ABSH0442
C*****OTHER COMMENTS:                                                  *ABSH0443
C*****BCAUX IS REQUIRED BY D02AGF ODE SOLVER IN NAG LIBRARY.           *ABSH0444
C***********************************************************************ABSH0445
      PARAMETER MNMOMS=3,MNEQ=6                                         ABSH0446
      COMMON /LOCMEQ/ XL0, DX, R0, A0, BT0, P0, AI0, E0, E1,            ABSH0447
     *     D0, ARC, SLPE, SLPR0, BETAP0, KMHD, IUNDER                   ABSH0448
      DIMENSION GL(MNEQ), GR(MNEQ), PARAM(MNEQ)                         ABSH0449
C*****LEFT B.C.'S - CENTER.                                             ABSH0450
      GL(1) = PARAM(1) - SLPR0*XL0**2                                   ABSH0451
      GL(2) = -2.0*SLPR0*XL0                                            ABSH0452
      GL(3) = PARAM(3) + SLPE*XL0**2                                    ABSH0453
      GL(4) = 2.0*XL0*SLPE                                              ABSH0454
      GL(5) = 10.0*PARAM(5)*XL0**4                                      ABSH0455
      GL(6) = 40.0*PARAM(5)*XL0**3                                      ABSH0456
C*****RIGHT B.C.'S - EDGE.                                              ABSH0457
      GR(1) = ARC                                                       ABSH0458
      GR(2) = PARAM(2)                                                  ABSH0459
      GR(3) = E0                                                        ABSH0460
      GR(4) = PARAM(4)                                                  ABSH0461
      GR(5) = D0                                                        ABSH0462
      GR(6) = PARAM(6)                                                  ABSH0463
      RETURN                                                            ABSH0464
      END                                                               ABSH0465
