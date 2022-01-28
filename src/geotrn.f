      SUBROUTINE GEOTRN(MODE, XIN, EIN, DIN, EOUT, DOUT)                ABSH0495
C***********************************************************************ABSH0496
C*****GEOTRN TRANSFORMS ELONGATION AND TRIANGULARITY FROM A MOMENTS TO *ABSH0497
C*****A GEOMETRICAL REPRESENTATION AND VICE VERSA.                     *ABSH0498
C*****REFERENCES:                                                      *ABSH0499
C*****L.L.LAO,S.P.HIRSHMAN,R.M.WIELAND,ORNL/TM-7616 (1981).            *ABSH0500
C*****                                ,PHYS FLUIDS 24 (1981) AUG       *ABSH0501
C*****L.L.LAO,R.M.WIELAND,W.A.HOULBERG,S.P.HIRSHMAN,ORNL/TM-7871 (1981)*ABSH0502
C*****LAST REVISION: 6/81 L.L.LAO, R.M.WIELAND, AND W.A.HOULBERG ORNL. *ABSH0503
C*****CALCULATED PARAMETERS:                                           *ABSH0504
C*****EOUT-OUTPUT ELONGATION.                                          *ABSH0505
C*****DOUT-OUTPUT TRIANGULARITY.                                       *ABSH0506
C*****INPUT PARAMETERS:                                                *ABSH0507
C*****MODE-DESIGNATES DIRECTION OF TRANSFORMATION.                     *ABSH0508
C*****    =1 GEOMETRICAL => MOMENTS.                                   *ABSH0509
C*****    =2 MOMENTS => GEOMETRICAL.                                   *ABSH0510
C*****XIN-INPUT REDUCED MINOR RADIUS                                   *ABSH0511
C*****EIN-INPUT ELONGATION.                                            *ABSH0512
C*****DIN-INPUT TRIANGULARITY.                                         *ABSH0513
C***********************************************************************ABSH0514
      IF (XIN.LE.0.0) GO TO 30                                          ABSH0515
      IF (MODE.EQ.2) GO TO 20                                           ABSH0516
C*****GEOMETRICAL ==> MOMENTS REPRESENTATION.                           ABSH0517
      DOUT = DIN/4.0                                                    ABSH0518
      CTC = 0.0                                                         ABSH0519
      DO 10 I=1,10                                                      ABSH0520
           CTC = 4.0*DOUT/(SQRT(XIN**2+32.0*DOUT**2)+XIN)               ABSH0521
           DOUT = XIN*DIN/(4.0-6.0*CTC**2)                              ABSH0522
   10 CONTINUE                                                          ABSH0523
      EOUT = XIN*EIN/(SQRT(1.0-CTC*CTC)*(XIN+2.0*DOUT*CTC))             ABSH0524
      RETURN                                                            ABSH0525
C*****MOMENTS ==> GEOMETRICAL REPRESENTATION.                           ABSH0526
   20 CTC = 4.0*DIN/(SQRT(XIN**2+32.0*DIN**2)+XIN)                      ABSH0527
      STC = SQRT(1.0-CTC*CTC)                                           ABSH0528
      S2TC = 2.0*STC*CTC                                                ABSH0529
      C2TC = 2.0*CTC*CTC - 1.0                                          ABSH0530
      EOUT = EIN*(XIN*STC+DIN*S2TC)/XIN                                 ABSH0531
      DOUT = DIN*(1.0-3.0*C2TC)/XIN                                     ABSH0532
      RETURN                                                            ABSH0533
   30 EOUT = EIN                                                        ABSH0534
      DOUT = DIN                                                        ABSH0535
      RETURN                                                            ABSH0536
      END                                                               ABSH0537
