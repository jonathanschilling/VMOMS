      SUBROUTINE D02AGF(H, ERROR, PARERR, PARAM, C, N, N1, M1,          ABSH1124
     *     AUX, BCAUX, RAAUX, PRSOL, MAT, COPY, WSPACE, WSPAC1,         ABSH1125
     *     WSPAC2, IFAIL)                                               ABSH1126
C     NAG COPYRIGHT 1975                                                ABSH1127
C     EDITED BY JOYCE CLARKE OXFORD OEG NUCLEAR PHYSICS 11TH SEP 1976   ABSH1128
C                   FORTRAN MACRO VERSION FDIA24.TEC                    ABSH1129
C     MARK 4.5 REVISED                                                  ABSH1130
C     ADDITIONAL COMMENTS BY R. WIELAND // ORNL // JAN.,1981.           ABSH1131
C                                                                       ABSH1132
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC ABSH1133
C                                                                       ABSH1134
C   THE PARAMETERS ARE DEFINED AS FOLLOWS:                              ABSH1135
C                                                                       ABSH1136
C     H           INPUT ESTIMATE OF REQUIRED INTEGRATION STEP LENGTH.   ABSH1137
C                 CHANGED ON RETURN TO PROVIDE LAST STEP LENGTH USED.   ABSH1138
C                                                                       ABSH1139
C     ERROR(N)    INPUT A REAL ARRAY USED AS:                           ABSH1140
C                 A) A LOCAL ERROR BOUND FOR INTEGRATION                ABSH1141
C                 B) A CONVERGENCE TOLERANCE ON Y AT THE MATCHING POINT.ABSH1142
C                    ON RETURN, THE FINAL ERROR.                        ABSH1143
C                                                                       ABSH1144
C     PARERR(N1)  INPUT A REAL ARRAY USED AS:                           ABSH1145
C                 A) A CONVERGENCE TOLERANCE FOR THE PARAMETERS (P)     ABSH1146
C                 B) USED TO APPROXIMATE DELTA-P IN ESTIMATING THE      ABSH1147
C                    JACOBIAN. ON RETURN, THE FINAL ERROR.              ABSH1148
C                                                                       ABSH1149
C     PARAM(N1)   INPUT STARTING VALUES FOR P;                          ABSH1150
C                 ON RETURN, THE CORRECTED VALUES                       ABSH1151
C                                                                       ABSH1152
C     C(M1,N)     THE SOLUTION Y(I;J) OF THE J-TH COMPONENT OF Y        ABSH1153
C                 EVALUATED AT X(I) IS RETURNED IN C (I,J); THE X(I)    ABSH1154
C                 SPACING IS DETERMINED BY M1, BELOW.                   ABSH1155
C                                                                       ABSH1156
C     N           NUMBER OF ODE'S                                       ABSH1157
C                                                                       ABSH1158
C     N1          NUMBER OF PARAMETERS                                  ABSH1159
C                                                                       ABSH1160
C     M1          THE FINAL SOLUTION IS CALCULATED AT M1 EQUIDISTANT    ABSH1161
C                 POINTS.                                               ABSH1162
C                                                                       ABSH1163
C     AUX         A USER SUPPLIED SUBROUTINE:                           ABSH1164
C                 SUBROUTINE AUX (F,Y,X,PARAM)                          ABSH1165
C                 WHERE X,Y(N),PARAM(N1) ARE USED TO EVALUTE THE        ABSH1166
C                 N DERIVATIVES F(N) AT X.                              ABSH1167
C                                                                       ABSH1168
C     BCAUX       USER SUPPLIED SUBROUTINE:                             ABSH1169
C                 SUBROUTINE BCAUX(G,G1,PARAM)                          ABSH1170
C                 WHERE PARAM IS USED (IF NECESSARY) TO EVALUATE Y(N)   ABSH1171
C                 AT THE ENDPOINTS X AND X1 AND RETURN THEM IN G(N)     ABSH1172
C                 AND G1(N).                                            ABSH1173
C                                                                       ABSH1174
C     RAAUX       USER SUPPLIED SUBROUTINE:                             ABSH1175
C                 SUBROUTINE RAAUX(X,X1,R,PARAM)                        ABSH1176
C                 WHERE PARAM IS USED (IF NECESSARY) TO EVALUATE THE    ABSH1177
C                 ENDPOINTS X AND X1, AND THE MATCHING POINT R.         ABSH1178
C                                                                       ABSH1179
C     PRSOL       A USER SUPPLIED SUBROUTINE:                           ABSH1180
C                 SUBROUTINE PRSOL(PARAM,RES,N1,ERR)                    ABSH1181
C                 CALLED AT EACH NEWTON ITERATION; CAN BE USED TO OUTPUTABSH1182
C                 ANY OF THE PARAMETERS OF INTEREST; ERR(N) ARE THE     ABSH1183
C                 ERRORS AT R IN EACH COMPONENT Y, AND RES IS THE       ABSH1184
C                 SUM OF THE SQUARES OF THESE ERRORS.                   ABSH1185
C                                                                       ABSH1186
C     MAT(N1,N1)   A REAL WORK ARRAY                                    ABSH1187
C     COPY(N1,N1)  A REAL WORK ARRAY                                    ABSH1188
C     WSPACE(N,9) A REAL WORK ARRAY                                     ABSH1189
C     WSPAC1(N)   A REAL WORK ARRAY                                     ABSH1190
C     WSPAC2(N)   A REAL WORK ARRAY                                     ABSH1191
C                                                                       ABSH1192
C     IFAIL       AN ERROR FLAG; ON INPUT, ENTER 0, ON OUTPUT:          ABSH1193
C       IFAIL = 0   NORMAL RETURN                                       ABSH1194
C             = 1   N1 > N                                              ABSH1195
C             = 2   INTEGRATION FAILED TO CONVERGE WHILE CALC. JACOBIAN ABSH1196
C             = 3   THE CONDITION X <= R <= X1 DOES NOT HOLD.           ABSH1197
C             = 4   INTEGRATION FAILED TO CONVERGE                      ABSH1198
C             = 5   JACOBIAN IS SINGULAR                                ABSH1199
C             = 6   AFTER 3 ATTEMPTS AT HALVING DELTA-P, THE NEWTON-    ABSH1200
C                   RAPHSON LOOP STILL DOES NOT YIELD A DIMINISHING     ABSH1201
C                   RESIDUAL.                                           ABSH1202
C             = 7   AFTER 12 N-R ITERATIONS, THERE IS STILL NOT         ABSH1203
C                   SUFFICIENT CONVERGENCE.                             ABSH1204
C                                                                       ABSH1205
C    INTERNAL FLAG :                                                    ABSH1206
C    VRBOSE (LOGICAL)   .TRUE. TO GENERATE DEBUG OUTPUT ON NDEBUG       ABSH1207
C                                                                       ABSH1208
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC ABSH1209
      INTEGER COUNT, COUNT1, EM, ONE, CT, P01AAF, N, N1, M1,            ABSH1210
     *     IFAIL, M, I, K, ITEST, J                                     ABSH1211
      DOUBLE PRECISION SRNAME                                           ABSH1212
      REAL MAT(N1,N1), H, ERROR(N), PARERR(N1), PARAM(N1),              ABSH1213
     *     C(M1,N), COPY(N1,N1), WSPACE(N,9), WSPAC1(N),                ABSH1214
     *     WSPAC2(N), EPS, H0, X, X1, R, DUM, RESID, D, PERT,           ABSH1215
     *     OLDRES, DIST, C1, X02AAF                                     ABSH1216
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1217
      COMMON /CFEVAL/ NFL                                               ABSH1218
      COMMON /NIO/ NIN, NOUT, NPLOT, NDEBUG, NTTY                       ABSH1219
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1220
      EXTERNAL AUX                                                      ABSH1221
      LOGICAL VRBOSE                                                    ABSH1222
      DATA VRBOSE /.FALSE./                                             ABSH1223
      DATA SRNAME /8H D02AGF /                                          ABSH1224
C                                                                       ABSH1225
C     SOLVES A GENERAL BOUNDARY VALUE                                   ABSH1226
C     PROBLEM FOR N DIFFERENTIAL EQUATIONS                              ABSH1227
C     IN N1 PARAMETERS USING A SHOOTING                                 ABSH1228
C     AND MATCHING TECHNIQUE.  EPS IS THE                               ABSH1229
C     LARGEST REAL VARIABLE SUCH THAT 1+EPS=1                           ABSH1230
C     ALL IMPLICITLY DECLARED REALS MAY BE USED DOUBLE-LENGTH           ABSH1231
C     THE ARRAY COPY IS REDUNDANT                                       ABSH1232
      EPS = SPMPAR(1)                                                   ABSH1233
      M = M1 - 1                                                        ABSH1234
      IF (N1.LE.N) GO TO 10                                             ABSH1235
C                                                                       ABSH1236
C     *** IFAIL = 1 ***                                                 ABSH1237
      EM = 1                                                            ABSH1238
      GO TO 480                                                         ABSH1239
C                                                                       ABSH1240
C     *** SET NEWTON ITERATION COUNTER ***                              ABSH1241
   10 COUNT = 0                                                         ABSH1242
      H0 = H                                                            ABSH1243
      ONE = 1                                                           ABSH1244
      EM = -1                                                           ABSH1245
C                                                                       ABSH1246
C     FORMS THE RESIDUALS AT THE                                        ABSH1247
C     MATCHING POINT                                                    ABSH1248
C                                                                       ABSH1249
C     *** GET LEFT-HAND PT., RIGHT-HAND PT., MATCHING PT. ***           ABSH1250
   20 CALL RAAUX(X, X1, R, PARAM)                                       ABSH1251
      IF ((X-R)*(X1-R).LE.0.0) GO TO 30                                 ABSH1252
C                                                                       ABSH1253
C     *** IFAIL = 3 ***                                                 ABSH1254
      EM = 3                                                            ABSH1255
      GO TO 480                                                         ABSH1256
C                                                                       ABSH1257
   30 IF (H0*(X1-X).LT.0.0) H0 = -H0                                    ABSH1258
C     ### G(I)  ==> W1            ###                                   ABSH1259
C     ### G1(I) ==> W2 ==> W(I,8) ###                                   ABSH1260
      CALL BCAUX(WSPAC1, WSPAC2, PARAM)                                 ABSH1261
      H = H0                                                            ABSH1262
      DO 40 I=1,N                                                       ABSH1263
           WSPACE(I,8) = WSPAC2(I)                                      ABSH1264
   40 CONTINUE                                                          ABSH1265
      I = 1                                                             ABSH1266
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1267
      IF (VRBOSE) WRITE (NDEBUG,99999)                                  ABSH1268
      NFL = 0                                                           ABSH1269
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1270
C                                                                       ABSH1271
C     *** INTEGRATE FROM X ==> R // ON RETURN: Y(R-) ==> W1(I) ***      ABSH1272
      CALL D02AGZ(X, WSPAC1, ERROR, ONE, N, N1, I, R-X, H, AUX,         ABSH1273
     *     WSPACE, WSPAC2, PARAM)                                       ABSH1274
      IF (I.EQ.0) GO TO 50                                              ABSH1275
C                                                                       ABSH1276
C     *** IFAIL = 4 ***                                                 ABSH1277
      EM = 4                                                            ABSH1278
      GO TO 480                                                         ABSH1279
C                                                                       ABSH1280
C     ###  -Y(R-) ==> W(I,8) // G1(I) ==> W1(I) ###                     ABSH1281
   50 DO 60 I=1,N                                                       ABSH1282
           DUM = WSPACE(I,8)                                            ABSH1283
           WSPACE(I,8) = -WSPAC1(I)                                     ABSH1284
           WSPAC1(I) = DUM                                              ABSH1285
   60 CONTINUE                                                          ABSH1286
C                                                                       ABSH1287
C     *** INTEGRATE FROM X1 ==> R // ON RETURN: YR(+) ==> W1(I) ***     ABSH1288
      H = -H0                                                           ABSH1289
      I = 1                                                             ABSH1290
      CALL D02AGZ(X1, WSPAC1, ERROR, ONE, N, N1, I, R-X1, H,            ABSH1291
     *     AUX, WSPACE, WSPAC2, PARAM)                                  ABSH1292
      IF (I.EQ.0) GO TO 70                                              ABSH1293
C                                                                       ABSH1294
C     *** IFAIL = 4 ***                                                 ABSH1295
      EM = 4                                                            ABSH1296
      GO TO 480                                                         ABSH1297
C                                                                       ABSH1298
   70 RESID = 0.0                                                       ABSH1299
      CT = 0                                                            ABSH1300
C                                                                       ABSH1301
C     *** FORM ERROR RESIDUALS AT R:                  ***               ABSH1302
C     ###   R(I) = Y(R+) - Y(R-) ==> W(I,8) == S(I;P) ###               ABSH1303
C     ###   S(I;P) ==> W1(I)                          ###               ABSH1304
      DO 80 I=1,N1                                                      ABSH1305
           D = WSPAC1(I)                                                ABSH1306
           DUM = WSPACE(I,8)                                            ABSH1307
           WSPACE(I,8) = D + DUM                                        ABSH1308
           D = 1.0 + ABS(DUM) + ABS(D)                                  ABSH1309
           DUM = WSPACE(I,8)                                            ABSH1310
           IF (ABS(DUM).LT.ERROR(I)*D) CT = CT + 1                      ABSH1311
           WSPAC1(I) = DUM                                              ABSH1312
           RESID = RESID + DUM*DUM                                      ABSH1313
   80 CONTINUE                                                          ABSH1314
      CALL PRSOL(PARAM, RESID, N1, WSPAC1)                              ABSH1315
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1316
      IF (VRBOSE) WRITE (NDEBUG,99998) CT, NFL                          ABSH1317
      IF (VRBOSE) WRITE (NDEBUG,99997) RESID                            ABSH1318
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1319
C     *** EM = -1 FIRST TIME THRU ONLY ***                              ABSH1320
      IF (EM.NE.-1) GO TO 320                                           ABSH1321
   90 COUNT = COUNT + 1                                                 ABSH1322
      IF (COUNT.NE.12) GO TO 100                                        ABSH1323
C                                                                       ABSH1324
C     *** IFAIL = 7 ***                                                 ABSH1325
      EM = 7                                                            ABSH1326
      GO TO 480                                                         ABSH1327
C                                                                       ABSH1328
C     FORMS THE JACOBIAN BY NUMERICAL                                   ABSH1329
C     DIFFERENTIATION                                                   ABSH1330
C     *** DEL$P(K)$ ==> PERT ***                                        ABSH1331
C     *** P(K) + DEL$P(K)$ ==> P(K+) ***                                ABSH1332
  100 CONTINUE                                                          ABSH1333
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1334
      NFL = 0                                                           ABSH1335
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1336
      DO 180 K=1,N1                                                     ABSH1337
           PERT = 10.0*PARERR(K)*(1.0+ABS(PARAM(K)))                    ABSH1338
           PARAM(K) = PERT + PARAM(K)                                   ABSH1339
           CALL RAAUX(X, X1, R, PARAM)                                  ABSH1340
           IF ((X-R)*(X1-R).LE.0.0) GO TO 110                           ABSH1341
C                                                                       ABSH1342
C     *** IFAIL = 3 ***                                                 ABSH1343
           EM = 3                                                       ABSH1344
           GO TO 480                                                    ABSH1345
C                                                                       ABSH1346
  110      IF (H0*(X1-X).LT.0.0) H0 = -H0                               ABSH1347
           H = H0                                                       ABSH1348
           CALL BCAUX(WSPAC1, WSPAC2, PARAM)                            ABSH1349
C     ### G1(I)$P(K+)$ ==> W(I,7) ###                                   ABSH1350
           DO 120 I=1,N                                                 ABSH1351
                WSPACE(I,7) = WSPAC2(I)                                 ABSH1352
  120      CONTINUE                                                     ABSH1353
           I = 1                                                        ABSH1354
C                                                                       ABSH1355
C     *** INTEGRATE FROM X ==> R USING P(K+) ***                        ABSH1356
           CALL D02AGZ(X, WSPAC1, ERROR, ONE, N, N1, I, R-X, H,         ABSH1357
     *          AUX, WSPACE, WSPAC2, PARAM)                             ABSH1358
           IF (I.EQ.0) GO TO 130                                        ABSH1359
C                                                                       ABSH1360
C     *** IFAIL = 2 ***                                                 ABSH1361
           EM = 2                                                       ABSH1362
           GO TO 480                                                    ABSH1363
C                                                                       ABSH1364
C     ### M(I,K) = Y(R-;P(K+)) ###                                      ABSH1365
  130      DO 140 I=1,N1                                                ABSH1366
                MAT(I,K) = WSPAC1(I)                                    ABSH1367
  140      CONTINUE                                                     ABSH1368
           H = -H0                                                      ABSH1369
           I = 1                                                        ABSH1370
           DO 150 I=1,N                                                 ABSH1371
                WSPAC1(I) = WSPACE(I,7)                                 ABSH1372
  150      CONTINUE                                                     ABSH1373
C                                                                       ABSH1374
C     *** INTEGRATE FROM X1 ==> R USING P(K+) ***                       ABSH1375
           CALL D02AGZ(X1, WSPAC1, ERROR, ONE, N, N1, I, R-X1,          ABSH1376
     *          H, AUX, WSPACE, WSPAC2, PARAM)                          ABSH1377
           IF (I.EQ.0) GO TO 160                                        ABSH1378
C                                                                       ABSH1379
C     *** IFAIL = 2 ***                                                 ABSH1380
           EM = 2                                                       ABSH1381
           GO TO 480                                                    ABSH1382
C                                                                       ABSH1383
C     ### M(I,K) = (S(P(K+)-S(P(K))/DEL$P(K)$ ###                       ABSH1384
  160      DO 170 I=1,N1                                                ABSH1385
                MAT(I,K) = (MAT(I,K)-WSPAC1(I)+WSPACE(I,8))/PERT        ABSH1386
                IF (ABS(MAT(I,K)).LT.5.0*EPS*ABS(WSPACE(I,8))/          ABSH1387
     *               PERT) MAT(I,K) = 0.0                               ABSH1388
  170      CONTINUE                                                     ABSH1389
           PARAM(K) = PARAM(K) - PERT                                   ABSH1390
  180 CONTINUE                                                          ABSH1391
C                                                                       ABSH1392
C     *** NEW JACOBIAN ***                                              ABSH1393
      ITEST = 1                                                         ABSH1394
      EM = -3                                                           ABSH1395
C                                                                       ABSH1396
C     PERFORMS COLUMN SCALING ON THE JACOBIAN                           ABSH1397
C     AND FORMS A TRIANGULAR DECOMPOSITION                              ABSH1398
      DO 220 I=1,N1                                                     ABSH1399
           D = 0.0                                                      ABSH1400
           DO 190 J=1,N1                                                ABSH1401
                IF (ABS(MAT(J,I)).GT.D) D = ABS(MAT(J,I))               ABSH1402
  190      CONTINUE                                                     ABSH1403
           IF (D.NE.0.0) GO TO 200                                      ABSH1404
C                                                                       ABSH1405
C     *** IFAIL = 5 ***                                                 ABSH1406
           EM = 5                                                       ABSH1407
           GO TO 480                                                    ABSH1408
C                                                                       ABSH1409
C     *** NORMALIZE M(I,K) ***                                          ABSH1410
  200      DO 210 J=1,N1                                                ABSH1411
                MAT(J,I) = MAT(J,I)/D                                   ABSH1412
  210      CONTINUE                                                     ABSH1413
           WSPACE(I,7) = D                                              ABSH1414
  220 CONTINUE                                                          ABSH1415
      I = 1                                                             ABSH1416
C                                                                       ABSH1417
C     *** LU DECOMPOSITION OF M(I,K) // PIVOT ARRAY: W1(I) ***          ABSH1418
      CALL F03AFF(N1, EPS, MAT, N1, D, J, WSPAC1, I)                    ABSH1419
      IF (I.EQ.0) GO TO 230                                             ABSH1420
C                                                                       ABSH1421
C     *** IFAIL = 5 ***                                                 ABSH1422
      EM = 5                                                            ABSH1423
      GO TO 480                                                         ABSH1424
C                                                                       ABSH1425
C     ### PIVOT ARRAY ==> W(I,6) ###                                    ABSH1426
  230 CONTINUE                                                          ABSH1427
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1428
      IF (VRBOSE) WRITE (NDEBUG,99996) NFL                              ABSH1429
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1430
      DO 240 I=1,N1                                                     ABSH1431
           WSPACE(I,6) = WSPAC1(I)                                      ABSH1432
  240 CONTINUE                                                          ABSH1433
C                                                                       ABSH1434
C     USES A GENERALISED NEWTON RAPHSON                                 ABSH1435
C     TECHNIQUE TO SOLVE THE NONLINEAR                                  ABSH1436
C     EQUATIONS AT THE MATCHING POINT                                   ABSH1437
C     *** SOLVE THE SYSTEM :  M(I,K) * DEL$P(K)$ = S(I;P) ***           ABSH1438
  250 OLDRES = RESID                                                    ABSH1439
C                                                                       ABSH1440
C     *** SET NEWTON - RAPHSON COUNTER ***                              ABSH1441
      COUNT1 = 0                                                        ABSH1442
C                                                                       ABSH1443
C     ### S(P) ==> W1(I) ###                                            ABSH1444
      DO 260 I=1,N1                                                     ABSH1445
           WSPAC1(I) = WSPACE(I,6)                                      ABSH1446
           WSPACE(I,1) = WSPACE(I,8)                                    ABSH1447
  260 CONTINUE                                                          ABSH1448
C                                                                       ABSH1449
C     *** THE NEW DEL$P$ = W1 ***                                       ABSH1450
      CALL F04AJF(N1, ONE, MAT, N1, WSPAC1, WSPACE, N)                  ABSH1451
C                                                                       ABSH1452
C     *** NORMALIZE DEL$P$ / FORM NEW P ***                             ABSH1453
      DO 270 I=1,N1                                                     ABSH1454
           WSPACE(I,1) = WSPACE(I,1)/WSPACE(I,7)                        ABSH1455
           PARAM(I) = PARAM(I) + WSPACE(I,1)                            ABSH1456
  270 CONTINUE                                                          ABSH1457
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1458
      IF (VRBOSE) WRITE (NDEBUG,99995) PARAM                            ABSH1459
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1460
      IF (CT.LT.N1) GO TO 300                                           ABSH1461
      DO 280 I=1,N1                                                     ABSH1462
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1463
           IPP = I                                                      ABSH1464
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1465
           IF (ABS(WSPACE(I,1)).GT.PARERR(I)*(1.0+ABS(PARAM(I)))        ABSH1466
     *          ) GO TO 290                                             ABSH1467
  280 CONTINUE                                                          ABSH1468
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1469
      IF (VRBOSE) WRITE (NDEBUG,99994)                                  ABSH1470
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1471
C                                                                       ABSH1472
C     *** SOLUTION FOUND] ***                                           ABSH1473
      EM = -5                                                           ABSH1474
      GO TO 380                                                         ABSH1475
  290 CONTINUE                                                          ABSH1476
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1477
      PIP = PARERR(IPP)*(1.0+ABS(PARAM(IPP)))                           ABSH1478
      IF (VRBOSE) WRITE (NDEBUG,99993) IPP, WSPACE(IPP,1), PIP          ABSH1479
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1480
  300 DO 310 I=1,N1                                                     ABSH1481
           WSPACE(I,1) = -WSPACE(I,1)                                   ABSH1482
  310 CONTINUE                                                          ABSH1483
      GO TO 20                                                          ABSH1484
C                                                                       ABSH1485
C     ******************************************************            ABSH1486
C                                                                       ABSH1487
  320 CONTINUE                                                          ABSH1488
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1489
      IF (VRBOSE) WRITE (NDEBUG,99992) OLDRES                           ABSH1490
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1491
      IF (COUNT1.NE.0) GO TO 330                                        ABSH1492
      IF (RESID.GE.OLDRES/10.0) GO TO 330                               ABSH1493
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1494
      IF (VRBOSE) WRITE (NDEBUG,99991)                                  ABSH1495
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1496
C                                                                       ABSH1497
C     *** GO BACK AND FORM NEW DEL$P$ USING CURRENT JACOBIAN, NEW S *** ABSH1498
      EM = -2                                                           ABSH1499
      ITEST = 0                                                         ABSH1500
      GO TO 250                                                         ABSH1501
C                                                                       ABSH1502
C     *** FORM NEW JACOBIAN ***                                         ABSH1503
  330 CONTINUE                                                          ABSH1504
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1505
      IF (VRBOSE) WRITE (NDEBUG,99990)                                  ABSH1506
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1507
      IF (RESID.LT.OLDRES) GO TO 90                                     ABSH1508
      IF (COUNT1.NE.3) GO TO 360                                        ABSH1509
      IF (ITEST.EQ.0) GO TO 340                                         ABSH1510
C                                                                       ABSH1511
C     *** IFAIL = 6 ***                                                 ABSH1512
      EM = 6                                                            ABSH1513
      GO TO 480                                                         ABSH1514
  340 CONTINUE                                                          ABSH1515
      DO 350 I=1,N1                                                     ABSH1516
           PARAM(I) = PARAM(I) + WSPACE(I,1)                            ABSH1517
  350 CONTINUE                                                          ABSH1518
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1519
      IF (VRBOSE) WRITE (NDEBUG,99989) PARAM                            ABSH1520
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1521
      EM = -1                                                           ABSH1522
      GO TO 20                                                          ABSH1523
  360 COUNT1 = COUNT1 + 1                                               ABSH1524
      EM = -4                                                           ABSH1525
C                                                                       ABSH1526
C     *** SCALE-DOWN DEL$P$ ***                                         ABSH1527
      DO 370 I=1,N1                                                     ABSH1528
           WSPACE(I,1) = WSPACE(I,1)/2.0                                ABSH1529
           PARAM(I) = PARAM(I) + WSPACE(I,1)                            ABSH1530
  370 CONTINUE                                                          ABSH1531
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1532
      IF (VRBOSE) WRITE (NDEBUG,99988) PARAM                            ABSH1533
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1534
      GO TO 20                                                          ABSH1535
C                                                                       ABSH1536
C     ******************************************************            ABSH1537
C                                                                       ABSH1538
C     CALCULATES THE FINAL SOLUTION                                     ABSH1539
  380 CONTINUE                                                          ABSH1540
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1541
      IF (VRBOSE) WRITE (NDEBUG,99987)                                  ABSH1542
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                                  ABSH1543
      IF (M.LE.0) GO TO 480                                             ABSH1544
      CALL RAAUX(X, X1, R, PARAM)                                       ABSH1545
      IF ((X-R)*(X1-R).LE.0.0) GO TO 390                                ABSH1546
      EM = 3                                                            ABSH1547
      GO TO 480                                                         ABSH1548
  390 IF (H0*(X1-X).LT.0.0) H0 = -H0                                    ABSH1549
      H = H0                                                            ABSH1550
      CALL BCAUX(WSPAC1, WSPAC2, PARAM)                                 ABSH1551
      DO 400 I=1,N                                                      ABSH1552
           WSPACE(I,7) = WSPAC2(I)                                      ABSH1553
  400 CONTINUE                                                          ABSH1554
      DIST = (X1-X)/FLOAT(M)                                            ABSH1555
      J = 1                                                             ABSH1556
      C1 = X                                                            ABSH1557
      K = 1                                                             ABSH1558
  410 DO 420 I=1,N                                                      ABSH1559
           C(J,I) = WSPAC1(I)                                           ABSH1560
  420 CONTINUE                                                          ABSH1561
      IF ((R-C1-0.25*DIST)*DIST.LE.0.0) GO TO 440                       ABSH1562
      I = 1                                                             ABSH1563
      CALL D02AGZ(C1, WSPAC1, ERROR, ONE, N, N1, I, DIST, H,            ABSH1564
     *     AUX, WSPACE, WSPAC2, PARAM)                                  ABSH1565
      IF (I.EQ.0) GO TO 430                                             ABSH1566
      EM = 4                                                            ABSH1567
      GO TO 480                                                         ABSH1568
  430 J = J + K                                                         ABSH1569
      GO TO 410                                                         ABSH1570
  440 IF (K.EQ.-1) GO TO 460                                            ABSH1571
      DIST = -DIST                                                      ABSH1572
      C1 = X1                                                           ABSH1573
      H = -H0                                                           ABSH1574
      J = M1                                                            ABSH1575
      K = -1                                                            ABSH1576
      DO 450 I=1,N                                                      ABSH1577
           WSPAC1(I) = WSPACE(I,7)                                      ABSH1578
  450 CONTINUE                                                          ABSH1579
      GO TO 410                                                         ABSH1580
  460 CALL BCAUX(WSPAC1, WSPAC2, PARAM)                                 ABSH1581
      DO 470 I=1,N                                                      ABSH1582
           C(1,I) = WSPAC1(I)                                           ABSH1583
           C(M1,I) = WSPAC2(I)                                          ABSH1584
  470 CONTINUE                                                          ABSH1585
  480 IF (EM.LE.0) IFAIL = 0                                            ABSH1586
      IF (EM.GT.0) IFAIL = P01AAF(IFAIL,EM,SRNAME)                      ABSH1587
      RETURN                                                            ABSH1588
C     END OF D02AGF                                                     ABSH1589
99999 FORMAT (10H FWD MARCH)                                            ABSH1590
99998 FORMAT (2X, I2, 2H Y13HS MATCH WITH , I6, 10H AUX EVALS)          ABSH1591
99997 FORMAT (16H PRSOL: RESID = , 1PE11.2)                             ABSH1592
99996 FORMAT (19H NEW JACOBIAN WITH , I6, 10H AUX EVALS)                ABSH1593
99995 FORMAT (17H SOLVE FOR NEW P:, /, 5X, 1P6E10.2)                    ABSH1594
99994 FORMAT (22H DEL P IS SMALL ENOUGH)                                ABSH1595
99993 FORMAT (10H  DEL P # , I1, 17H STILL TOO LARGE:, 1PE11.3,         ABSH1596
     *     2H >, E11.3)                                                 ABSH1597
99992 FORMAT (28H COMPARE AGAINST OLD RESID- , 1PE11.3)                 ABSH1598
99991 FORMAT (44H KEEP THE OLD JAC; GET NEW DEL P USING NEW S)          ABSH1599
99990 FORMAT (39H NEW RESID IS NOT A FACTOR OF 10 BETTER)               ABSH1600
99989 FORMAT (39H RESID IS WORSE & LIM=3; RESTORE OLD P:, /,            ABSH1601
     *     5X, 1P6E10.2)                                                ABSH1602
99988 FORMAT (32H RESID IS WORSE; BACK DOWN DEL P, /, 5X,               ABSH1603
     *     6HNEW P:, /, 5X, 1P6E10.2)                                   ABSH1604
99987 FORMAT (9H FINISHED)                                              ABSH1605
      END                                                               ABSH1606
