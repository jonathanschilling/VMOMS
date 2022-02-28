      REAL FUNCTION SPMPAR(I)                                           ABSH2037
      INTEGER I                                                         ABSH2038
C     **********                                                        ABSH2039
C                                                                       ABSH2040
C     FUNCTION SPMPAR                                                   ABSH2041
C                                                                       ABSH2042
C     THIS FUNCTION PROVIDES SINGLE PRECISION MACHINE PARAMETERS        ABSH2043
C     WHEN THE APPROPRIATE SET OF DATA STATEMENTS IS ACTIVATED (BY      ABSH2044
C     REMOVING THE C FROM COLUMN 1) AND ALL OTHER DATA STATEMENTS ARE   ABSH2045
C     RENDERED INACTIVE. MOST OF THE PARAMETER VALUES WERE OBTAINED     ABSH2046
C     FROM THE CORRESPONDING BELL LABORATORIES PORT LIBRARY FUNCTION.   ABSH2047
C                                                                       ABSH2048
C     THE FUNCTION STATEMENT IS                                         ABSH2049
C                                                                       ABSH2050
C       REAL FUNCTION SPMPAR(I)                                         ABSH2051
C                                                                       ABSH2052
C     WHERE                                                             ABSH2053
C                                                                       ABSH2054
C       I IS AN INTEGER INPUT VARIABLE SET TO 1, 2, OR 3 WHICH          ABSH2055
C         SELECTS THE DESIRED MACHINE PARAMETER. IF THE MACHINE HAS     ABSH2056
C         T BASE B DIGITS AND ITS SMALLEST AND LARGEST EXPONENTS ARE    ABSH2057
C         EMIN AND EMAX, RESPECTIVELY, THEN THESE PARAMETERS ARE        ABSH2058
C                                                                       ABSH2059
C         SPMPAR(1) = B**(1 - T), THE MACHINE PRECISION,                ABSH2060
C                                                                       ABSH2061
C         SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,            ABSH2062
C                                                                       ABSH2063
C         SPMPAR(3) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.     ABSH2064
C                                                                       ABSH2065
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.         ABSH2066
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE             ABSH2067
C                                                                       ABSH2068
C     **********                                                        ABSH2069
      INTEGER MCHEPS(1)                                                 ABSH2070
      INTEGER MINMAG(1)                                                 ABSH2071
      INTEGER MAXMAG(1)                                                 ABSH2072
      REAL RMACH(3)                                                     ABSH2073
      EQUIVALENCE (RMACH(1),MCHEPS(1))                                  ABSH2074
      EQUIVALENCE (RMACH(2),MINMAG(1))                                  ABSH2075
      EQUIVALENCE (RMACH(3),MAXMAG(1))                                  ABSH2076
C                                                                       ABSH2077
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,                     ABSH2078
C     THE AMDAHL 470/V6, THE ICL 2900, THE ITEL AS/6,                   ABSH2079
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.                  ABSH2080
C                                                                       ABSH2081
C     DATA RMACH(1) / Z3C100000 /                                       ABSH2082
C     DATA RMACH(2) / Z00100000 /                                       ABSH2083
C     DATA RMACH(3) / Z7FFFFFFF /                                       ABSH2084
C                                                                       ABSH2085
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.              ABSH2086
C                                                                       ABSH2087
C     DATA RMACH(1) / O716400000000 /                                   ABSH2088
C     DATA RMACH(2) / O402400000000 /                                   ABSH2089
C     DATA RMACH(3) / O376777777777 /                                   ABSH2090
C                                                                       ABSH2091
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.                   ABSH2092
C                                                                       ABSH2093
C     DATA RMACH(1) / 16414000000000000000B /                           ABSH2094
C     DATA RMACH(2) / 00014000000000000000B /                           ABSH2095
C     DATA RMACH(3) / 37767777777777777777B /                           ABSH2096
C                                                                       ABSH2097
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).            ABSH2098
C                                                                       ABSH2099
C     DATA RMACH(1) / "147400000000 /                                   ABSH2100
C     DATA RMACH(2) / "000400000000 /                                   ABSH2101
C     DATA RMACH(3) / "377777777777 /                                   ABSH2102
C                                                                       ABSH2103
C     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING               ABSH2104
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).                 ABSH2105
C                                                                       ABSH2106
      DATA MCHEPS(1) /  889192448 /                                     ABSH2107
      DATA MINMAG(1) /    8388608 /                                     ABSH2108
      DATA MAXMAG(1) / 2147483647 /                                     ABSH2109
                                                                        ABSH2110
C      DATA RMACH(1) / O06500000000 /                                    ABSH2111
C      DATA RMACH(2) / O00040000000 /                                    ABSH2112
C      DATA RMACH(3) / O17777777777 /                                    ABSH2113
C                                                                       ABSH2114
C     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING               ABSH2115
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).                 ABSH2116
C                                                                       ABSH2117
C     DATA MCHEPS(1),MCHEPS(2) / 13568,     0 /                         ABSH2118
C     DATA MINMAG(1),MINMAG(2) /   128,     0 /                         ABSH2119
C     DATA MAXMAG(1),MAXMAG(2) / 32767,    -1 /                         ABSH2120
C                                                                       ABSH2121
C     DATA MCHEPS(1),MCHEPS(2) / O032400, O000000 /                     ABSH2122
C     DATA MINMAG(1),MINMAG(2) / O000200, O000000 /                     ABSH2123
C     DATA MAXMAG(1),MAXMAG(2) / O077777, O177777 /                     ABSH2124
C                                                                       ABSH2125
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS.       ABSH2126
C                                                                       ABSH2127
C     DATA RMACH(1) / O1301000000000000 /                               ABSH2128
C     DATA RMACH(2) / O1771000000000000 /                               ABSH2129
C     DATA RMACH(3) / O0777777777777777 /                               ABSH2130
C                                                                       ABSH2131
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.                  ABSH2132
C                                                                       ABSH2133
C     DATA RMACH(1) / Z4EA800000 /                                      ABSH2134
C     DATA RMACH(2) / Z400800000 /                                      ABSH2135
C     DATA RMACH(3) / Z5FFFFFFFF /                                      ABSH2136
C                                                                       ABSH2137
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.                     ABSH2138
C                                                                       ABSH2139
C     DATA RMACH(1) / O147400000000 /                                   ABSH2140
C     DATA RMACH(2) / O000400000000 /                                   ABSH2141
C     DATA RMACH(3) / O377777777777 /                                   ABSH2142
C                                                                       ABSH2143
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.             ABSH2144
C                                                                       ABSH2145
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -      ABSH2146
C     STATIC RMACH(3)                                                   ABSH2147
C                                                                       ABSH2148
C     DATA MINMAG/20K,0/,MAXMAG/77777K,177777K/                         ABSH2149
C     DATA MCHEPS/36020K,0/                                             ABSH2150
C                                                                       ABSH2151
C     MACHINE CONSTANTS FOR THE HARRIS 220.                             ABSH2152
C                                                                       ABSH2153
C     DATA MCHEPS(1) / '20000000, '00000353 /                           ABSH2154
C     DATA MINMAG(1) / '20000000, '00000201 /                           ABSH2155
C     DATA MAXMAG(1) / '37777777, '00000177 /                           ABSH2156
C                                                                       ABSH2157
C     MACHINE CONSTANTS FOR THE CRAY-1.                                 ABSH2158
C                                                                       ABSH2159
C      DATA RMACH(1) / 0377224000000000000000B /                        ABSH2160
C      DATA RMACH(2) / 0200034000000000000000B /                        ABSH2161
C      DATA RMACH(3) / 0577777777777777777776B /                        ABSH2162
C                                                                       ABSH2163
C     MACHINE CONSTANTS FOR THE PRIME 400.                              ABSH2164
C                                                                       ABSH2165
C     DATA MCHEPS(1) / :10000000153 /                                   ABSH2166
C     DATA MINMAG(1) / :10000000000 /                                   ABSH2167
C     DATA MAXMAG(1) / :17777777777 /                                   ABSH2168
C                                                                       ABSH2169
      SPMPAR = RMACH(I)                                                 ABSH2170
      RETURN                                                            ABSH2171
C                                                                       ABSH2172
C     LAST CARD OF FUNCTION SPMPAR.                                     ABSH2173
C                                                                       ABSH2174
      END                                                               ABSH2175
C                                                                       ABSH2176
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     ABSH2177
C                                                                       ABSH2178
