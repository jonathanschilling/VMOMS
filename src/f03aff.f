      SUBROUTINE F03AFF(N, EPS, A, IA, D1, ID, P, IFAIL)                ABSH1797
C     MARK 2 RELEASE. NAG COPYRIGHT 1972                                ABSH1798
C     EDITED BY JOYCE CLARKE OXFORD OEG NUCLEAR PHYSICS 11TH SEP 1976   ABSH1799
C                   FORTRAN MACRO VERSION FDIA24.TEC                    ABSH1800
C     MARK 3 REVISED.                                                   ABSH1801
C     MARK 4.5 REVISED                                                  ABSH1802
C     SW3 SWITCH REVISION : R. WIELAND // ORNL //   JAN., 1981.         ABSH1803
C     .FALSE. GIVES NO EXTENDED PRECISION                               ABSH1804
C     .TRUE. GIVES EXTENDED PRECISION IF CONDITIONS COMMENTED           ABSH1805
C     ON IN X03AAZ ARE FULFILLED.                                       ABSH1806
C                                                                       ABSH1807
C     UNSYMDET                                                          ABSH1808
C     THE UNSYMMETRIC MATRIX, A, IS STORED IN THE N*N ARRAY A(I,J),     ABSH1809
C     I=1,N, J=1,N. THE DECOMPOSITION A=LU, WHERE L IS A                ABSH1810
C     LOWER TRIANGULAR MATRIX AND U IS A UNIT UPPER TRIANGULAR          ABSH1811
C     MATRIX, IS PERFORMED AND OVERWRITTEN ON A, OMITTING THE UNIT      ABSH1812
C     DIAGONAL OF U. A RECORD OF ANY INTERCHANGES MADE TO THE ROWS      ABSH1813
C     OF A IS KEPT IN P(I), I=1,N, SUCH THAT THE I-TH ROW AND           ABSH1814
C     THE P(I)-TH ROW WERE INTERCHANGED AT THE I-TH STEP. THE           ABSH1815
C     DETERMINANT, D1 * 2.0**ID, OF A IS ALSO COMPUTED. THE             ABSH1816
C     SUBROUTINE                                                        ABSH1817
C     WILL FAIL IF A, MODIFIED BY THE ROUNDING ERRORS, IS SINGULAR      ABSH1818
C     OR ALMOST SINGULAR. SETS IFAIL = 0 IF SUCCESSFUL ELSE IFAIL =     ABSH1819
C     1.                                                                ABSH1820
C     1ST DECEMBER 1971                                                 ABSH1821
C                                                                       ABSH1822
      INTEGER ISAVE, IFAIL, IFAIL1, I, N, IA, ID, K, L, K1, K2,         ABSH1823
     *     ISTART, J, P01AAF                                            ABSH1824
      DOUBLE PRECISION SRNAME                                           ABSH1825
      REAL Y, D2, D1, X, EPS, A(IA,N), P(N)                             ABSH1826
      LOGICAL SW3                                                       ABSH1827
      DATA SRNAME /8H F03AFF /                                          ABSH1828
      DATA SW3 /.FALSE./                                                ABSH1829
      ISAVE = IFAIL                                                     ABSH1830
      IFAIL1 = 0                                                        ABSH1831
      DO 10 I=1,N                                                       ABSH1832
           CALL X03AAF(A(I,1), IA*N+1-I, A(I,1), IA*N+1-I, N,           ABSH1833
     *          IA, IA, 0.0, 0.0, Y, D2, SW3, IFAIL1)                   ABSH1834
           IF (Y.LE.0.0) GO TO 100                                      ABSH1835
           P(I) = 1.0/SQRT(Y)                                           ABSH1836
   10 CONTINUE                                                          ABSH1837
      D1 = 1.0                                                          ABSH1838
      ID = 0                                                            ABSH1839
      DO 90 K=1,N                                                       ABSH1840
           L = K                                                        ABSH1841
           X = 0.0                                                      ABSH1842
           K1 = K - 1                                                   ABSH1843
           K2 = K + 1                                                   ABSH1844
           ISTART = K                                                   ABSH1845
           DO 20 I=ISTART,N                                             ABSH1846
                CALL X03AAF(A(I,1), N*IA+1-I, A(1,K),                   ABSH1847
     *               (N-K+1)*IA, K1, IA, 1, -A(I,K), 0.0, Y,            ABSH1848
     *               D2, SW3, IFAIL1)                                   ABSH1849
                A(I,K) = -Y                                             ABSH1850
                Y = ABS(Y*P(I))                                         ABSH1851
                IF (Y.LE.X) GO TO 20                                    ABSH1852
                X = Y                                                   ABSH1853
                L = I                                                   ABSH1854
   20      CONTINUE                                                     ABSH1855
           IF (L.EQ.K) GO TO 40                                         ABSH1856
           D1 = -D1                                                     ABSH1857
           DO 30 J=1,N                                                  ABSH1858
                Y = A(K,J)                                              ABSH1859
                A(K,J) = A(L,J)                                         ABSH1860
                A(L,J) = Y                                              ABSH1861
   30      CONTINUE                                                     ABSH1862
           P(L) = P(K)                                                  ABSH1863
   40      P(K) = L                                                     ABSH1864
           D1 = D1*A(K,K)                                               ABSH1865
           IF (X.LT.8.0*EPS) GO TO 100                                  ABSH1866
   50      IF (ABS(D1).LT.1.0) GO TO 60                                 ABSH1867
           D1 = D1*0.0625                                               ABSH1868
           ID = ID + 4                                                  ABSH1869
           GO TO 50                                                     ABSH1870
   60      IF (ABS(D1).GE.0.0625) GO TO 70                              ABSH1871
           D1 = D1*16.0                                                 ABSH1872
           ID = ID - 4                                                  ABSH1873
           GO TO 60                                                     ABSH1874
   70      X = -1.0/A(K,K)                                              ABSH1875
           IF (K2.GT.N) GO TO 90                                        ABSH1876
           DO 80 J=K2,N                                                 ABSH1877
                CALL X03AAF(A(K,1), N*IA+1-K, A(1,J),                   ABSH1878
     *               (N-J+1)*IA, K1, IA, 1, -A(K,J), 0.0, Y,            ABSH1879
     *               D2, SW3, IFAIL1)                                   ABSH1880
                A(K,J) = X*Y                                            ABSH1881
   80      CONTINUE                                                     ABSH1882
   90 CONTINUE                                                          ABSH1883
      IFAIL = 0                                                         ABSH1884
      RETURN                                                            ABSH1885
  100 IFAIL = P01AAF(ISAVE,1,SRNAME)                                    ABSH1886
      RETURN                                                            ABSH1887
      END                                                               ABSH1888
