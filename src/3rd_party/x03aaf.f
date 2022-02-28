      SUBROUTINE X03AAF(A, ISIZEA, B, ISIZEB, N, ISTEPA,                ABSH1982
     *     ISTEPB, C1, C2, D1, D2, SW, IFAIL)                           ABSH1983
C     NAG COPYRIGHT 1975                                                ABSH1984
C     EDITED BY JOYCE CLARKE OXFORD OEG NUCLEAR PHYSICS 11TH SEP 1976   ABSH1985
C                   FORTRAN MACRO VERSION FDIA24.TEC                    ABSH1986
C     MARK 4.5 RELEASE                                                  ABSH1987
C                                                                       ABSH1988
C     CALCULATES THE VALUE OF A SCALAR PRODUCT USING BASIC              ABSH1989
C     OR ADDITIONAL PRECISION AND ADDS IT TO A BASIC OR ADDITIONAL      ABSH1990
C     PRECISION INITIAL VALUE.                                          ABSH1991
C     X03AAF CALLS X03AAZ WHICH MAY BE IN ASSEMBLY LANGUAGE.            ABSH1992
C                                                                       ABSH1993
      INTEGER P01AAF, ISAVE, ISIZEA, ISIZEB, ISTEPA, ISTEPB,            ABSH1994
     *     IFAIL, IS, IT, N, I                                          ABSH1995
      DOUBLE PRECISION SUM, SRNAME                                      ABSH1996
      REAL A(ISIZEA), B(ISIZEB), C1, C2, D1, D2, X                      ABSH1997
      LOGICAL SW                                                        ABSH1998
      DATA SRNAME /8H X03AAF /                                          ABSH1999
      ISAVE = IFAIL                                                     ABSH2000
      IFAIL = 0                                                         ABSH2001
      IF (ISTEPA.GT.0 .AND. ISTEPB.GT.0) GO TO 10                       ABSH2002
      IFAIL = P01AAF(ISAVE,1,SRNAME)                                    ABSH2003
      RETURN                                                            ABSH2004
   10 IS = 1 - ISTEPA                                                   ABSH2005
      IT = 1 - ISTEPB                                                   ABSH2006
      IF (SW) GO TO 40                                                  ABSH2007
      X = 0.0                                                           ABSH2008
      IF (N.LT.1) GO TO 30                                              ABSH2009
      DO 20 I=1,N                                                       ABSH2010
           IS = IS + ISTEPA                                             ABSH2011
           IT = IT + ISTEPB                                             ABSH2012
           X = X + A(IS)*B(IT)                                          ABSH2013
   20 CONTINUE                                                          ABSH2014
   30 D1 = X + (C1+C2)                                                  ABSH2015
      D2 = 0.0                                                          ABSH2016
      RETURN                                                            ABSH2017
   40 SUM = 0.0D0                                                       ABSH2018
      IF (N.LT.1) GO TO 60                                              ABSH2019
      DO 50 I=1,N                                                       ABSH2020
           IS = IS + ISTEPA                                             ABSH2021
           IT = IT + ISTEPB                                             ABSH2022
           SUM = SUM + DBLE(A(IS))*B(IT)                                ABSH2023
   50 CONTINUE                                                          ABSH2024
   60 SUM = SUM + (DBLE(C1)+C2)                                         ABSH2025
      CALL X03AAZ(SUM, D1, D2)                                          ABSH2026
      RETURN                                                            ABSH2027
      END                                                               ABSH2028
