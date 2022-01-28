      INTEGER FUNCTION P01AAF(IFAIL, ERROR, SRNAME)                     ABSH1948
C     MARK 1 RELEASE.  NAG COPYRIGHT 1971                               ABSH1949
C     MARK 3 REVISED                                                    ABSH1950
C     MARK 4A REVISED, IER-45                                           ABSH1951
C     MARK 4.5 REVISED                                                  ABSH1952
C     MARK 7 REVISED (DEC 1978) .... (APR 1979)                         ABSH1953
C     NOUT = NTTY :: REVISED BY R. WIELAND // ORNL // JAN., 1981.       ABSH1954
C     RETURNS THE VALUE OF ERROR OR TERMINATES THE PROGRAM.             ABSH1955
C     IF A HARD FAILURE OCCURS, THIS ROUTINE CALLS A FORTRAN AUXILIARY  ABSH1956
C     ROUTINE P01AAZ WHICH GIVES A TRACE, A FAILURE MESSAGE AND HALTS   ABSH1957
C     THE PROGRAM                                                       ABSH1958
      INTEGER ERROR, IFAIL, NOUT                                        ABSH1959
      DOUBLE PRECISION SRNAME                                           ABSH1960
      COMMON /NIO/ NIN, NOUT, NPLOT, NDEBUG, NTTY                       ABSH1961
C     TEST IF NO ERROR DETECTED                                         ABSH1962
      IF (ERROR.EQ.0) GO TO 20                                          ABSH1963
C     TEST FOR SOFT FAILURE                                             ABSH1964
      IF (MOD(IFAIL,10).EQ.1) GO TO 10                                  ABSH1965
C     HARD FAILURE                                                      ABSH1966
      WRITE (NTTY,99999) SRNAME, ERROR                                  ABSH1967
C     STOPPING MECHANISM MAY ALSO DIFFER                                ABSH1968
C      CALL P01AAZ (X)                                                  ABSH1969
C     STOP                                                              ABSH1970
C     SOFT FAIL                                                         ABSH1971
C     TEST IF ERROR MESSAGES SUPPRESSED                                 ABSH1972
   10 IF (MOD(IFAIL/10,10).EQ.0) GO TO 20                               ABSH1973
      WRITE (NTTY,99999) SRNAME, ERROR                                  ABSH1974
   20 P01AAF = ERROR                                                    ABSH1975
      RETURN                                                            ABSH1976
99999 FORMAT (1H0, 38HERROR DETECTED BY NAG LIBRARY ROUTINE ,           ABSH1977
     *     A8, 11H - IFAIL = , I5//)                                    ABSH1978
      END                                                               ABSH1979
C     AUTO EDIT 20 SEP 76                                               ABSH1980
C     AUTO EDIT 20 SEP 76                                               ABSH1981
