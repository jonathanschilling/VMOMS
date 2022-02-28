      SUBROUTINE D02AGY(Y, X, H, N, N1, AUX, WSPACE, WSPAC1,            ABSH1727
     *     PARAM)                                                       ABSH1728
C                                                                       ABSH1729
C     USES THE TECHNIQUE OF NAG LIBRARY                                 ABSH1730
C     PROCEDURE D02AAF WITH THE                                         ABSH1731
C     SPECIFICATION OF AUX CHANGED                                      ABSH1732
C     ALL IMPLICITLY DECLARED REALS MAY                                 ABSH1733
C     BE DECLARED DOUBLE PRECISION                                      ABSH1734
C     NAG COPYRIGHT 1975                                                ABSH1735
C     EDITED BY JOYCE CLARKE OXFORD OEG NUCLEAR PHYSICS 11TH SEP 1976   ABSH1736
C                   FORTRAN MACRO VERSION FDIA24.TEC                    ABSH1737
C     MARK 4.5 REVISED                                                  ABSH1738
C                                                                       ABSH1739
C   INTEGRATES A SYSTEM OF ODE'S USING MERSON'S METH.                   ABSH1740
C   CF.  LAMBERT, J.D., COMPUTATIONAL METHODS IN ORDINARY DIFFERENTIAL  ABSH1741
C        EQUATIONS, PP. 130-135, WILEY, 1973.                           ABSH1742
C                                                                       ABSH1743
      INTEGER N, N1, I                                                  ABSH1744
      REAL Y(N), X, H, WSPACE(N,9), WSPAC1(N), PARAM(N1), C1,           ABSH1745
     *     C2, U, V, W, DUM                                             ABSH1746
      CALL AUX(WSPAC1, Y, X, PARAM)                                     ABSH1747
      C1 = 1.0/3.0                                                      ABSH1748
      C2 = 1.0/6.0                                                      ABSH1749
      DO 10 I=1,N                                                       ABSH1750
           WSPACE(I,3) = WSPAC1(I)                                      ABSH1751
   10 CONTINUE                                                          ABSH1752
      U = C1*H                                                          ABSH1753
      DO 20 I=1,N                                                       ABSH1754
           WSPACE(I,4) = Y(I)                                           ABSH1755
           Y(I) = Y(I) + U*WSPACE(I,3)                                  ABSH1756
   20 CONTINUE                                                          ABSH1757
      CALL AUX(WSPAC1, Y, U+X, PARAM)                                   ABSH1758
      DO 30 I=1,N                                                       ABSH1759
           WSPACE(I,5) = WSPAC1(I)                                      ABSH1760
   30 CONTINUE                                                          ABSH1761
      V = H*C2                                                          ABSH1762
      DO 40 I=1,N                                                       ABSH1763
           Y(I) = WSPACE(I,4) + V*(WSPACE(I,3)+WSPACE(I,5))             ABSH1764
   40 CONTINUE                                                          ABSH1765
      CALL AUX(WSPAC1, Y, U+X, PARAM)                                   ABSH1766
      DO 50 I=1,N                                                       ABSH1767
           WSPACE(I,5) = WSPAC1(I)                                      ABSH1768
   50 CONTINUE                                                          ABSH1769
      U = H*0.125                                                       ABSH1770
      V = H*0.375                                                       ABSH1771
      DO 60 I=1,N                                                       ABSH1772
           Y(I) = WSPACE(I,4) + WSPACE(I,3)*U + WSPACE(I,5)*V           ABSH1773
   60 CONTINUE                                                          ABSH1774
      U = 0.5*H                                                         ABSH1775
      V = 1.5*H                                                         ABSH1776
      W = 2.0*H                                                         ABSH1777
      CALL AUX(WSPAC1, Y, U+X, PARAM)                                   ABSH1778
      DO 70 I=1,N                                                       ABSH1779
           Y(I) = WSPACE(I,4) + WSPACE(I,3)*U - WSPACE(I,5)*V +         ABSH1780
     *          WSPAC1(I)*W                                             ABSH1781
           WSPACE(I,5) = WSPAC1(I)                                      ABSH1782
   70 CONTINUE                                                          ABSH1783
      X = X + H                                                         ABSH1784
      CALL AUX(WSPAC1, Y, X, PARAM)                                     ABSH1785
      U = H*C2                                                          ABSH1786
      V = 2.0*H*C1                                                      ABSH1787
      DO 80 I=1,N                                                       ABSH1788
           W = WSPACE(I,4) + U*(WSPACE(I,3)+WSPAC1(I)) +                ABSH1789
     *          WSPACE(I,5)*V                                           ABSH1790
           DUM = W - Y(I)                                               ABSH1791
           WSPACE(I,2) = 0.2*ABS(DUM)                                   ABSH1792
           Y(I) = W                                                     ABSH1793
   80 CONTINUE                                                          ABSH1794
      RETURN                                                            ABSH1795
      END                                                               ABSH1796
