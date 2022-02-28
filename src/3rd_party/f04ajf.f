      SUBROUTINE F04AJF(N, IR, A, IA, P, B, IB)                         ABSH1889
C     MARK 2 RELEASE. NAG COPYRIGHT 1972                                ABSH1890
C     EDITED BY JOYCE CLARKE OXFORD OEG NUCLEAR PHYSICS 11TH SEP 1976   ABSH1891
C                   FORTRAN MACRO VERSION FDIA24.TEC                    ABSH1892
C     MARK 4 REVISED.                                                   ABSH1893
C     MARK 4.5 REVISED                                                  ABSH1894
C     SW3 SWITCH REVISION : R. WIELAND // ORNL //   JAN., 1981.         ABSH1895
C     .FALSE. GIVES NO EXTENDED PRECISION                               ABSH1896
C     .TRUE. GIVES EXTENDED PRECISION IF CONDITIONS COMMENTED           ABSH1897
C     ON IN X03AAZ ARE FULFILLED.                                       ABSH1898
C     UNSYMSOL                                                          ABSH1899
C     SOLVES AX=B, WHERE A IS AN UNSYMMETRIC MATRIX AND B IS AN         ABSH1900
C     N*IR                                                              ABSH1901
C     MATRIX OF IR RIGHT-HAND SIDES. THE SUBROUTINE F04AJF MUST BY      ABSH1902
C     PRECEDED BY F03AFF IN WHICH L AND U ARE PRODUCED IN A(I,J),       ABSH1903
C     FROM A, AND THE RECORD OF THE INTERCHANGES IS PRODUCED IN         ABSH1904
C     P(I). AX=B IS SOLVED IN THREE STEPS, INTERCHANGE THE              ABSH1905
C     ELEMENTS OF B, LY=B AND UX=Y. THE MATRICES Y AND THEN X ARE       ABSH1906
C     OVERWRITTEN ON B.                                                 ABSH1907
C     1ST AUGUST 1971                                                   ABSH1908
C                                                                       ABSH1909
      INTEGER IFAIL1, I, N, I1, K, IR, IA, IB, II                       ABSH1910
      REAL X, D1, D2, A(IA,N), P(N), B(IB,IR)                           ABSH1911
      LOGICAL SW3                                                       ABSH1912
      DATA SW3 /.FALSE./                                                ABSH1913
      IFAIL1 = 0                                                        ABSH1914
C     INTERCHANGING OF ELEMENTS OF B                                    ABSH1915
      DO 20 I=1,N                                                       ABSH1916
           I1 = P(I) + 0.5                                              ABSH1917
           IF (I1.EQ.I) GO TO 20                                        ABSH1918
           DO 10 K=1,IR                                                 ABSH1919
                X = B(I,K)                                              ABSH1920
                B(I,K) = B(I1,K)                                        ABSH1921
                B(I1,K) = X                                             ABSH1922
   10      CONTINUE                                                     ABSH1923
   20 CONTINUE                                                          ABSH1924
      DO 50 K=1,IR                                                      ABSH1925
C     SOLUTION OF LY= B                                                 ABSH1926
           DO 30 I=1,N                                                  ABSH1927
                I1 = I - 1                                              ABSH1928
                CALL X03AAF(A(I,1), N*IA-I+1, B(1,K),                   ABSH1929
     *               (IR-K+1)*IB, I1, IA, 1, B(I,K), 0.0, D1,           ABSH1930
     *               D2, SW3, IFAIL1)                                   ABSH1931
                B(I,K) = -D1/A(I,I)                                     ABSH1932
   30      CONTINUE                                                     ABSH1933
C     SOLUTION OF UX= Y                                                 ABSH1934
           B(N,K) = -B(N,K)                                             ABSH1935
           IF (N.EQ.1) GO TO 50                                         ABSH1936
           DO 40 II=2,N                                                 ABSH1937
                I = N - II + 1                                          ABSH1938
                I1 = I + 1                                              ABSH1939
                CALL X03AAF(A(I,I1), (N-I)*IA-I+1, B(I1,K),             ABSH1940
     *               (IR-K+1)*IB-I, N-I, IA, 1, B(I,K), 0.0,            ABSH1941
     *               D1, D2, SW3, IFAIL1)                               ABSH1942
                B(I,K) = -D1                                            ABSH1943
   40      CONTINUE                                                     ABSH1944
   50 CONTINUE                                                          ABSH1945
      RETURN                                                            ABSH1946
      END                                                               ABSH1947
