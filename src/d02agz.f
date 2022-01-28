      SUBROUTINE D02AGZ(X, Y, G, T, N, N1, M, H0, H, AUX,               ABSH1607
     *     WSPACE, WSPAC1, PARAM)                                       ABSH1608
C                                                                       ABSH1609
C     USES THE LOGIC OF NAG LIBRARY ROUTINE D02ABF                      ABSH1610
C     EPS AND DUM MAY BE USED DOUBLE LENGTH.  EPS IS THE                ABSH1611
C     SMALLEST REAL SUCH THAT 1+EPS>EPS.  SMAX                          ABSH1612
C     IS THE LARGEST INTEGER. DUM,ERR AND HS MAYBE DECLARED             ABSH1613
C     DOUBLE PRECISION                                                  ABSH1614
C     NAG COPYRIGHT 1975                                                ABSH1615
C     EDITED BY JOYCE CLARKE OXFORD OEG NUCLEAR PHYSICS 11TH SEP 1976   ABSH1616
C                   FORTRAN MACRO VERSION FDIA24.TEC                    ABSH1617
C     MARK 4.5 REVISED                                                  ABSH1618
C     ADDITIONAL COMMENTS BY R. WIELAND // ORNL // JAN.,1981.           ABSH1619
C                                                                       ABSH1620
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC             ABSH1621
C     THE SUBROUTINE PARAMETERS ARE DEFINED AS FOLLOWS:                 ABSH1622
C                                                                       ABSH1623
C     X       INPUT INITAL X VALUE; ON RETURN X=X+H0 IF SUCCESSFUL      ABSH1624
C                                                                       ABSH1625
C     Y(N)    INPUT INITIAL VALUES OF Y(X); ON RETURN Y(X+H0) IF        ABSH1626
C             SUCCESSFUL                                                ABSH1627
C                                                                       ABSH1628
C     G(N)    INPUT ERROR BOUNDS ON Y(N)                                ABSH1629
C                                                                       ABSH1630
C     T       INPUT TYPE OF ERROR TEST                                  ABSH1631
C                                                                       ABSH1632
C     N       NUMBER OF ODE'S IN THE SYSTEM                             ABSH1633
C                                                                       ABSH1634
C     M       INPUT ERROR FLAG AS 0; ON RETURN:                         ABSH1635
C             M = 0   NORMAL RETURN                                     ABSH1636
C               = 1   STEP LENGTH REPEATEDLY HALVED TO < 10**-4 THE     ABSH1637
C                     INITIALLY SUGGESTED STEP LENGTH                   ABSH1638
C               = 2   T .NE. 1 OR 2  OR 3                               ABSH1639
C               = 3   THE NUMBER OF STEPS REQUIRED EXCEEDS THE LARGEST  ABSH1640
C                     INTEGER ON THE MACHINE                            ABSH1641
C                                                                       ABSH1642
C     H0      INTERVAL OF INTEGRATION                                   ABSH1643
C                                                                       ABSH1644
C     H       INPUT AN ESTIMATE OF THE STEP LENGTH REQUIRED;            ABSH1645
C             ON RETURN, THE FINAL STEP LENGTH USED.                    ABSH1646
C                                                                       ABSH1647
C      THE REMAINING PARAMETERS ARE AS DEFINED IN 'D02AGF'              ABSH1648
C                                                                       ABSH1649
C     ON RETURN:                                                        ABSH1650
C       WSPACE(I,2) CONTAINS THE LOCAL ERROR AT X0+H FOR EACH Y(I)      ABSH1651
C       WSPACE(I,5) CONTAINS Y'(X+H0)                                   ABSH1652
C       WSPACE(I,3) CONTAINS Y'(X)                                      ABSH1653
C                                                                       ABSH1654
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCABSH1655
C                                                                       ABSH1656
      INTEGER D, S, SMAX, T, N, N1, M, I, J                             ABSH1657
      REAL X, Y(N), G(N), H0, H, WSPACE(N,9), WSPAC1(N),                ABSH1658
     *     PARAM(N1), EPS, DUM, P, Q, HS, ERR, X02AAF                   ABSH1659
      EXTERNAL AUX                                                      ABSH1660
C     *** H0 = INTERVAL LENGTH // H = STEP LENGTH ***                   ABSH1661
      EPS = SPMPAR(1)                                                   ABSH1662
      SMAX = I1MACH(9)                                                  ABSH1663
      M = 0                                                             ABSH1664
      DUM = AMAX1(ABS(X),ABS(X+H0))                                     ABSH1665
      IF (ABS(H0).LE.EPS*DUM) RETURN                                    ABSH1666
      DUM = H/H0                                                        ABSH1667
      IF (ABS(DUM).GT.1.0 .OR. H.EQ.0.0) GO TO 10                       ABSH1668
      DUM = ABS(H0/H+0.9)                                               ABSH1669
      IF (DUM.GT.FLOAT(SMAX)) GO TO 80                                  ABSH1670
      H = H0/AINT(DUM)                                                  ABSH1671
      GO TO 20                                                          ABSH1672
   10 H = H0                                                            ABSH1673
C                                                                       ABSH1674
C     *** T = 1 :: MIXED    ERROR TEST    ***                           ABSH1675
C     *** T = 2 :: ABSOLUTE ERROR TEST    ***                           ABSH1676
C     *** T = 3 :: RELATIVE ERROR TEST    ***                           ABSH1677
   20 P = 1.0                                                           ABSH1678
      IF (T.EQ.3) P = 0.0                                               ABSH1679
      Q = 1.0                                                           ABSH1680
      IF (T.EQ.2) Q = 0.0                                               ABSH1681
      DUM = H0/H + 0.1                                                  ABSH1682
C     *** S = # OF ANTICIPATED STEPS       ***                          ABSH1683
C     *** HS = MINIMUM STEP LENGTH ALLOWED ***                          ABSH1684
      S = INT(DUM)                                                      ABSH1685
      HS = 1.0E-4*ABS(H)                                                ABSH1686
   30 DO 40 I=1,N                                                       ABSH1687
           WSPACE(I,9) = Y(I)                                           ABSH1688
   40 CONTINUE                                                          ABSH1689
C                                                                       ABSH1690
C     *** INTEGRATE FROM X ==> X+H;                              ***    ABSH1691
C     ***   INPUT:  Y(X) ==> Y // X ==> X;                       ***    ABSH1692
C     ***   OUTPUT: Y(X+H) ==> Y // X+H ==> X;                   ***    ABSH1693
C     ***           LOCAL ERROR ==> W(I,2) // Y'(X+H) ==> W(I,5) ***    ABSH1694
C     ***           Y'(X) ==> W(I,3)                             ***    ABSH1695
      CALL D02AGY(Y, X, H, N, N1, AUX, WSPACE, WSPAC1, PARAM)           ABSH1696
      D = 0                                                             ABSH1697
      DO 50 I=1,N                                                       ABSH1698
           ERR = G(I)*(P+Q*ABS(Y(I)))                                   ABSH1699
           IF (WSPACE(I,2).GT.ERR) GO TO 60                             ABSH1700
           IF (40.0*WSPACE(I,2).GT.ERR) D = 1                           ABSH1701
   50 CONTINUE                                                          ABSH1702
      S = S - 1                                                         ABSH1703
C                                                                       ABSH1704
C     *** THE ONLY SUCCESSFUL WAY OUT ***                               ABSH1705
      IF (S.EQ.0) RETURN                                                ABSH1706
C                                                                       ABSH1707
      DUM = FLOAT(S)/2.0 + 0.1                                          ABSH1708
      J = INT(DUM)*2                                                    ABSH1709
      IF (D.NE.0 .OR. J.NE.S) GO TO 30                                  ABSH1710
      H = 2.0*H                                                         ABSH1711
      S = S/2                                                           ABSH1712
      GO TO 30                                                          ABSH1713
C                                                                       ABSH1714
C     *** LOCAL ERROR TOO LARGE ***                                     ABSH1715
   60 X = X - H                                                         ABSH1716
      DO 70 I=1,N                                                       ABSH1717
           Y(I) = WSPACE(I,9)                                           ABSH1718
   70 CONTINUE                                                          ABSH1719
      IF (S.GT.SMAX/2) GO TO 80                                         ABSH1720
      S = 2*S                                                           ABSH1721
      H = 0.5*H                                                         ABSH1722
      IF (ABS(H).GT.HS) GO TO 30                                        ABSH1723
   80 M = 1                                                             ABSH1724
      RETURN                                                            ABSH1725
      END                                                               ABSH1726
