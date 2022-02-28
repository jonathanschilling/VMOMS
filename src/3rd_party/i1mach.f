      INTEGER FUNCTION I1MACH(I)                                        ABSH2179
C                                                                       ABSH2180
C  I/O UNIT NUMBERS.                                                    ABSH2181
C                                                                       ABSH2182
C    I1MACH( 1) = THE STANDARD INPUT UNIT.                              ABSH2183
C                                                                       ABSH2184
C    I1MACH( 2) = THE STANDARD OUTPUT UNIT.                             ABSH2185
C                                                                       ABSH2186
C    I1MACH( 3) = THE STANDARD PUNCH UNIT.                              ABSH2187
C                                                                       ABSH2188
C    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.                      ABSH2189
C                                                                       ABSH2190
C  WORDS.                                                               ABSH2191
C                                                                       ABSH2192
C    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.          ABSH2193
C                                                                       ABSH2194
C    I1MACH( 6) = THE NUMBER OF CHARACTERS PER INTEGER STORAGE UNIT.    ABSH2195
C                                                                       ABSH2196
C  INTEGERS.                                                            ABSH2197
C                                                                       ABSH2198
C    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM        ABSH2199
C                                                                       ABSH2200
C               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )          ABSH2201
C                                                                       ABSH2202
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.               ABSH2203
C                                                                       ABSH2204
C    I1MACH( 7) = A, THE BASE.                                          ABSH2205
C                                                                       ABSH2206
C    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.                       ABSH2207
C                                                                       ABSH2208
C    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.                      ABSH2209
C                                                                       ABSH2210
C  FLOATING-POINT NUMBERS.                                              ABSH2211
C                                                                       ABSH2212
C    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,      ABSH2213
C    BASE-B FORM                                                        ABSH2214
C                                                                       ABSH2215
C               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )            ABSH2216
C                                                                       ABSH2217
C               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,                 ABSH2218
C               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.                 ABSH2219
C                                                                       ABSH2220
C    I1MACH(10) = B, THE BASE.                                          ABSH2221
C                                                                       ABSH2222
C  SINGLE-PRECISION                                                     ABSH2223
C                                                                       ABSH2224
C    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.                       ABSH2225
C                                                                       ABSH2226
C    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.                        ABSH2227
C                                                                       ABSH2228
C    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.                         ABSH2229
C                                                                       ABSH2230
C  DOUBLE-PRECISION                                                     ABSH2231
C                                                                       ABSH2232
C    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.                       ABSH2233
C                                                                       ABSH2234
C    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.                        ABSH2235
C                                                                       ABSH2236
C    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.                         ABSH2237
C                                                                       ABSH2238
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,                 ABSH2239
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY            ABSH2240
C  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF                   ABSH2241
C  I1MACH(1) - I1MACH(4) SHOULD BE CHECKED FOR CONSISTENCY              ABSH2242
C  WITH THE LOCAL OPERATING SYSTEM.                                     ABSH2243
C                                                                       ABSH2244
      INTEGER IMACH(16), OUTPUT                                         ABSH2245
C                                                                       ABSH2246
      EQUIVALENCE (IMACH(4),OUTPUT)                                     ABSH2247
C                                                                       ABSH2248
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.                  ABSH2249
C                                                                       ABSH2250
C     DATA IMACH( 1) /    7 /                                           ABSH2251
C     DATA IMACH( 2) /    2 /                                           ABSH2252
C     DATA IMACH( 3) /    2 /                                           ABSH2253
C     DATA IMACH( 4) /    2 /                                           ABSH2254
C     DATA IMACH( 5) /   36 /                                           ABSH2255
C     DATA IMACH( 6) /    4 /                                           ABSH2256
C     DATA IMACH( 7) /    2 /                                           ABSH2257
C     DATA IMACH( 8) /   33 /                                           ABSH2258
C     DATA IMACH( 9) / Z1FFFFFFFF /                                     ABSH2259
C     DATA IMACH(10) /    2 /                                           ABSH2260
C     DATA IMACH(11) /   24 /                                           ABSH2261
C     DATA IMACH(12) / -256 /                                           ABSH2262
C     DATA IMACH(13) /  255 /                                           ABSH2263
C     DATA IMACH(14) /   60 /                                           ABSH2264
C     DATA IMACH(15) / -256 /                                           ABSH2265
C     DATA IMACH(16) /  255 /                                           ABSH2266
C                                                                       ABSH2267
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.                  ABSH2268
C                                                                       ABSH2269
C     DATA IMACH( 1) /   5 /                                            ABSH2270
C     DATA IMACH( 2) /   6 /                                            ABSH2271
C     DATA IMACH( 3) /   7 /                                            ABSH2272
C     DATA IMACH( 4) /   6 /                                            ABSH2273
C     DATA IMACH( 5) /  48 /                                            ABSH2274
C     DATA IMACH( 6) /   6 /                                            ABSH2275
C     DATA IMACH( 7) /   2 /                                            ABSH2276
C     DATA IMACH( 8) /  39 /                                            ABSH2277
C     DATA IMACH( 9) / O0007777777777777 /                              ABSH2278
C     DATA IMACH(10) /   8 /                                            ABSH2279
C     DATA IMACH(11) /  13 /                                            ABSH2280
C     DATA IMACH(12) / -50 /                                            ABSH2281
C     DATA IMACH(13) /  76 /                                            ABSH2282
C     DATA IMACH(14) /  26 /                                            ABSH2283
C     DATA IMACH(15) / -50 /                                            ABSH2284
C     DATA IMACH(16) /  76 /                                            ABSH2285
C                                                                       ABSH2286
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.            ABSH2287
C                                                                       ABSH2288
C     DATA IMACH( 1) /   5 /                                            ABSH2289
C     DATA IMACH( 2) /   6 /                                            ABSH2290
C     DATA IMACH( 3) /   7 /                                            ABSH2291
C     DATA IMACH( 4) /   6 /                                            ABSH2292
C     DATA IMACH( 5) /  48 /                                            ABSH2293
C     DATA IMACH( 6) /   6 /                                            ABSH2294
C     DATA IMACH( 7) /   2 /                                            ABSH2295
C     DATA IMACH( 8) /  39 /                                            ABSH2296
C     DATA IMACH( 9) / O0007777777777777 /                              ABSH2297
C     DATA IMACH(10) /   8 /                                            ABSH2298
C     DATA IMACH(11) /  13 /                                            ABSH2299
C     DATA IMACH(12) / -50 /                                            ABSH2300
C     DATA IMACH(13) /  76 /                                            ABSH2301
C     DATA IMACH(14) /  26 /                                            ABSH2302
C     DATA IMACH(15) / -32754 /                                         ABSH2303
C     DATA IMACH(16) /  32780 /                                         ABSH2304
C                                                                       ABSH2305
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.                   ABSH2306
C                                                                       ABSH2307
C     DATA IMACH( 1) /    5 /                                           ABSH2308
C     DATA IMACH( 2) /    6 /                                           ABSH2309
C     DATA IMACH( 3) /    7 /                                           ABSH2310
C     DATA IMACH( 4) /    6 /                                           ABSH2311
C     DATA IMACH( 5) /   60 /                                           ABSH2312
C     DATA IMACH( 6) /   10 /                                           ABSH2313
C     DATA IMACH( 7) /    2 /                                           ABSH2314
C     DATA IMACH( 8) /   48 /                                           ABSH2315
C     DATA IMACH( 9) / 00007777777777777777B /                          ABSH2316
C     DATA IMACH(10) /    2 /                                           ABSH2317
C     DATA IMACH(11) /   48 /                                           ABSH2318
C     DATA IMACH(12) / -974 /                                           ABSH2319
C     DATA IMACH(13) / 1070 /                                           ABSH2320
C     DATA IMACH(14) /   96 /                                           ABSH2321
C     DATA IMACH(15) / -927 /                                           ABSH2322
C     DATA IMACH(16) / 1070 /                                           ABSH2323
C                                                                       ABSH2324
C     MACHINE CONSTANTS FOR THE CRAY 1                                  ABSH2325
C                                                                       ABSH2326
C      DATA IMACH( 1) /   100 /                                         ABSH2327
C      DATA IMACH( 2) /   5 /                                           ABSH2328
C      DATA IMACH( 3) /   102 /                                         ABSH2329
C      DATA IMACH( 4) /   101 /                                         ABSH2330
C      DATA IMACH( 5) /    64 /                                         ABSH2331
C      DATA IMACH( 6) /     8 /                                         ABSH2332
C      DATA IMACH( 7) /     2 /                                         ABSH2333
C      DATA IMACH( 8) /    63 /                                         ABSH2334
C      DATA IMACH( 9) /  777777777777777777777B /                       ABSH2335
C      DATA IMACH(10) /     2 /                                         ABSH2336
C      DATA IMACH(11) /    48 /                                         ABSH2337
C      DATA IMACH(12) / -8192 /                                         ABSH2338
C      DATA IMACH(13) /  8191 /                                         ABSH2339
C      DATA IMACH(14) /    96 /                                         ABSH2340
C      DATA IMACH(15) / -8192 /                                         ABSH2341
C      DATA IMACH(16) /  8191 /                                         ABSH2342
C                                                                       ABSH2343
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200              ABSH2344
C                                                                       ABSH2345
C     DATA IMACH( 1) /   11 /                                           ABSH2346
C     DATA IMACH( 2) /   12 /                                           ABSH2347
C     DATA IMACH( 3) /    8 /                                           ABSH2348
C     DATA IMACH( 4) /   10 /                                           ABSH2349
C     DATA IMACH( 5) /   16 /                                           ABSH2350
C     DATA IMACH( 6) /    2 /                                           ABSH2351
C     DATA IMACH( 7) /    2 /                                           ABSH2352
C     DATA IMACH( 8) /   15 /                                           ABSH2353
C     DATA IMACH( 9) /32767 /                                           ABSH2354
C     DATA IMACH(10) /   16 /                                           ABSH2355
C     DATA IMACH(11) /    6 /                                           ABSH2356
C     DATA IMACH(12) /  -64 /                                           ABSH2357
C     DATA IMACH(13) /   63 /                                           ABSH2358
C     DATA IMACH(14) /   14 /                                           ABSH2359
C     DATA IMACH(15) /  -64 /                                           ABSH2360
C     DATA IMACH(16) /   63 /                                           ABSH2361
C                                                                       ABSH2362
C     MACHINE CONSTANTS FOR THE HARRIS 220                              ABSH2363
C                                                                       ABSH2364
C     DATA IMACH( 1) /       5 /                                        ABSH2365
C     DATA IMACH( 2) /       6 /                                        ABSH2366
C     DATA IMACH( 3) /       0 /                                        ABSH2367
C     DATA IMACH( 4) /       6 /                                        ABSH2368
C     DATA IMACH( 5) /      24 /                                        ABSH2369
C     DATA IMACH( 6) /       3 /                                        ABSH2370
C     DATA IMACH( 7) /       2 /                                        ABSH2371
C     DATA IMACH( 8) /      23 /                                        ABSH2372
C     DATA IMACH( 9) / 8388607 /                                        ABSH2373
C     DATA IMACH(10) /       2 /                                        ABSH2374
C     DATA IMACH(11) /      23 /                                        ABSH2375
C     DATA IMACH(12) /    -127 /                                        ABSH2376
C     DATA IMACH(13) /     127 /                                        ABSH2377
C     DATA IMACH(14) /      38 /                                        ABSH2378
C     DATA IMACH(15) /    -127 /                                        ABSH2379
C     DATA IMACH(16) /     127 /                                        ABSH2380
C                                                                       ABSH2381
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.              ABSH2382
C                                                                       ABSH2383
C     DATA IMACH( 1) /    5 /                                           ABSH2384
C     DATA IMACH( 2) /    6 /                                           ABSH2385
C     DATA IMACH( 3) /   43 /                                           ABSH2386
C     DATA IMACH( 4) /    6 /                                           ABSH2387
C     DATA IMACH( 5) /   36 /                                           ABSH2388
C     DATA IMACH( 6) /    6 /                                           ABSH2389
C     DATA IMACH( 7) /    2 /                                           ABSH2390
C     DATA IMACH( 8) /   35 /                                           ABSH2391
C     DATA IMACH( 9) / O377777777777 /                                  ABSH2392
C     DATA IMACH(10) /    2 /                                           ABSH2393
C     DATA IMACH(11) /   27 /                                           ABSH2394
C     DATA IMACH(12) / -127 /                                           ABSH2395
C     DATA IMACH(13) /  127 /                                           ABSH2396
C     DATA IMACH(14) /   63 /                                           ABSH2397
C     DATA IMACH(15) / -127 /                                           ABSH2398
C     DATA IMACH(16) /  127 /                                           ABSH2399
C                                                                       ABSH2400
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,                     ABSH2401
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.                  ABSH2402
C                                                                       ABSH2403
C     DATA IMACH( 1) /   5 /                                            ABSH2404
C     DATA IMACH( 2) /   6 /                                            ABSH2405
C     DATA IMACH( 3) /   7 /                                            ABSH2406
C     DATA IMACH( 4) /   6 /                                            ABSH2407
C     DATA IMACH( 5) /  32 /                                            ABSH2408
C     DATA IMACH( 6) /   4 /                                            ABSH2409
C     DATA IMACH( 7) /   2 /                                            ABSH2410
C     DATA IMACH( 8) /  31 /                                            ABSH2411
C     DATA IMACH( 9) / Z7FFFFFFF /                                      ABSH2412
C     DATA IMACH(10) /  16 /                                            ABSH2413
C     DATA IMACH(11) /   6 /                                            ABSH2414
C     DATA IMACH(12) / -64 /                                            ABSH2415
C     DATA IMACH(13) /  63 /                                            ABSH2416
C     DATA IMACH(14) /  14 /                                            ABSH2417
C     DATA IMACH(15) / -64 /                                            ABSH2418
C     DATA IMACH(16) /  63 /                                            ABSH2419
C                                                                       ABSH2420
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).                  ABSH2421
C                                                                       ABSH2422
C     DATA IMACH(1) /5/                                                 ABSH2423
C     DATA IMACH(2) /6/                                                 ABSH2424
C     DATA IMACH(3) /5/                                                 ABSH2425
C     DATA IMACH(4) /6/                                                 ABSH2426
C     DATA IMACH(5) /36/                                                ABSH2427
C     DATA IMACH(6) /5/                                                 ABSH2428
C     DATA IMACH(7) /2/                                                 ABSH2429
C     DATA IMACH(8) /35/                                                ABSH2430
C     DATA IMACH( 9) / "377777777777 /                                  ABSH2431
C     DATA IMACH(10) /2/                                                ABSH2432
C     DATA IMACH(11) /27/                                               ABSH2433
C     DATA IMACH(12) /-128/                                             ABSH2434
C     DATA IMACH(13) /127/                                              ABSH2435
C     DATA IMACH(14) /54/                                               ABSH2436
C     DATA IMACH(15) /-101/                                             ABSH2437
C     DATA IMACH(16) /127/                                              ABSH2438
C                                                                       ABSH2439
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).                  ABSH2440
C                                                                       ABSH2441
C     DATA IMACH( 1) /    5 /                                           ABSH2442
C     DATA IMACH( 2) /    6 /                                           ABSH2443
C     DATA IMACH( 3) /    5 /                                           ABSH2444
C     DATA IMACH( 4) /    6 /                                           ABSH2445
C     DATA IMACH( 5) /   36 /                                           ABSH2446
C     DATA IMACH( 6) /    5 /                                           ABSH2447
C     DATA IMACH( 7) /    2 /                                           ABSH2448
C     DATA IMACH( 8) /   35 /                                           ABSH2449
C     DATA IMACH( 9) / "377777777777 /                                  ABSH2450
C     DATA IMACH(10) /    2 /                                           ABSH2451
C     DATA IMACH(11) /   27 /                                           ABSH2452
C     DATA IMACH(12) / -128 /                                           ABSH2453
C     DATA IMACH(13) /  127 /                                           ABSH2454
C     DATA IMACH(14) /   62 /                                           ABSH2455
C     DATA IMACH(15) / -128 /                                           ABSH2456
C     DATA IMACH(16) /  127 /                                           ABSH2457
C                                                                       ABSH2458
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN'S SUPPORTING                 ABSH2459
C     32-BIT INTEGER ARITHMETIC.                                        ABSH2460
C                                                                       ABSH2461
      DATA IMACH( 1) /    5 /                                           ABSH2462
      DATA IMACH( 2) /    6 /                                           ABSH2463
      DATA IMACH( 3) /    5 /                                           ABSH2464
      DATA IMACH( 4) /    6 /                                           ABSH2465
      DATA IMACH( 5) /   32 /                                           ABSH2466
      DATA IMACH( 6) /    4 /                                           ABSH2467
      DATA IMACH( 7) /    2 /                                           ABSH2468
      DATA IMACH( 8) /   31 /                                           ABSH2469
      DATA IMACH( 9) / 2147483647 /                                     ABSH2470
      DATA IMACH(10) /    2 /                                           ABSH2471
      DATA IMACH(11) /   24 /                                           ABSH2472
      DATA IMACH(12) / -127 /                                           ABSH2473
      DATA IMACH(13) /  127 /                                           ABSH2474
      DATA IMACH(14) /   56 /                                           ABSH2475
      DATA IMACH(15) / -127 /                                           ABSH2476
      DATA IMACH(16) /  127 /                                           ABSH2477
C                                                                       ABSH2478
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN'S SUPPORTING                 ABSH2479
C     16-BIT INTEGER ARITHMETIC.                                        ABSH2480
C                                                                       ABSH2481
C     DATA IMACH( 1) /    5 /                                           ABSH2482
C     DATA IMACH( 2) /    6 /                                           ABSH2483
C     DATA IMACH( 3) /    5 /                                           ABSH2484
C     DATA IMACH( 4) /    6 /                                           ABSH2485
C     DATA IMACH( 5) /   16 /                                           ABSH2486
C     DATA IMACH( 6) /    2 /                                           ABSH2487
C     DATA IMACH( 7) /    2 /                                           ABSH2488
C     DATA IMACH( 8) /   15 /                                           ABSH2489
C     DATA IMACH( 9) / 32767 /                                          ABSH2490
C     DATA IMACH(10) /    2 /                                           ABSH2491
C     DATA IMACH(11) /   24 /                                           ABSH2492
C     DATA IMACH(12) / -127 /                                           ABSH2493
C     DATA IMACH(13) /  127 /                                           ABSH2494
C     DATA IMACH(14) /   56 /                                           ABSH2495
C     DATA IMACH(15) / -127 /                                           ABSH2496
C     DATA IMACH(16) /  127 /                                           ABSH2497
C                                                                       ABSH2498
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.                     ABSH2499
C                                                                       ABSH2500
C     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7            ABSH2501
C     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM.                   ABSH2502
C     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1.                   ABSH2503
C                                                                       ABSH2504
C     DATA IMACH( 1) /    5 /                                           ABSH2505
C     DATA IMACH( 2) /    6 /                                           ABSH2506
C     DATA IMACH( 3) /    7 /                                           ABSH2507
C     DATA IMACH( 4) /    6 /                                           ABSH2508
C     DATA IMACH( 5) /   36 /                                           ABSH2509
C     DATA IMACH( 6) /    6 /                                           ABSH2510
C     DATA IMACH( 7) /    2 /                                           ABSH2511
C     DATA IMACH( 8) /   35 /                                           ABSH2512
C     DATA IMACH( 9) / O377777777777 /                                  ABSH2513
C     DATA IMACH(10) /    2 /                                           ABSH2514
C     DATA IMACH(11) /   27 /                                           ABSH2515
C     DATA IMACH(12) / -128 /                                           ABSH2516
C     DATA IMACH(13) /  127 /                                           ABSH2517
C     DATA IMACH(14) /   60 /                                           ABSH2518
C     DATA IMACH(15) /-1024 /                                           ABSH2519
C     DATA IMACH(16) / 1023 /                                           ABSH2520
C                                                                       ABSH2521
      IF (I.LT.1 .OR. I.GT.16) GO TO 10                                 ABSH2522
C                                                                       ABSH2523
      I1MACH = IMACH(I)                                                 ABSH2524
      RETURN                                                            ABSH2525
C                                                                       ABSH2526
   10 WRITE (OUTPUT,99999)                                              ABSH2527
C                                                                       ABSH2528
C                                                                       ABSH2529
      STOP                                                              ABSH2530
C                                                                       ABSH2531
99999 FORMAT (39H1ERROR    1 IN I1MACH - I OUT OF BOUNDS)               ABSH2532
      END                                                               ABSH2533
