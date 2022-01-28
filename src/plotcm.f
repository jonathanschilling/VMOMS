CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                ABSH2535
C                                                                       ABSH2536
C  THIS PROGRAM GENERATES PLOTS OF THE 2D AND 1D                        ABSH2537
C  SOLUTIONS OF GRAD-SHAFRANOF EQUATION                                 ABSH2538
C  ORNL, 09/80                                                          ABSH2539
C                                                                       ABSH2540
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                ABSH2541
      PROGRAM PLOTCM                                                    ABSH2542
      PARAMETER M2P=10,M2PT=20,MNMX=61                                  ABSH2543
      PARAMETER NTHETA=101                                              ABSH2544
      DIMENSION RMOM(MNMX), PARAM(6)                                    ABSH2545
      DIMENSION XMOM(M2P), SMOM(M2P), EMOM(M2P), DMOM(M2P)              ABSH2546
      DIMENSION RHO(MNMX), SHIFT(MNMX), ELONG(MNMX),                    ABSH2547
     *     TRIANG(MNMX)                                                 ABSH2548
      DIMENSION RTHETA(NTHETA), ZTHETA(NTHETA), PSIMOM(M2P)             ABSH2549
      DIMENSION DVDX(M2P), SMQ(M2P), AJORF(M2P), FFP(M2P),              ABSH2550
     *     AJTM(M2PT), RMID(M2PT), FPOL(M2P), PXP(M2P),                 ABSH2551
     *     PPXP(M2P), AISPL(M2P)                                        ABSH2552
      DIMENSION YRFL(M2PT)                                              ABSH2553
      DATA TWOPI /6.28319/                                              ABSH2554
      DATA PI /3.141593/, RADIAN /.017453/, STATA /3.0E+09/             ABSH2555
C                                                                       ABSH2556
      NPLOT = 12                                                        ABSH2557
C                                                                       ABSH2558
      READ (NPLOT,99999) M2, NRHO                                       ABSH2559
      READ (NPLOT,99998) ARC, A0, DGEOM, PSIC0, PARAM                   ABSH2560
C                                                                       ABSH2561
      DO 10 I=1,M2                                                      ABSH2562
           READ (NPLOT,99998) XMOM(I), SMOM(I), EMOM(I), DMOM(I)        ABSH2563
           READ (NPLOT,99998) PSIMOM(I), DVDX(I), SMQ(I),               ABSH2564
     *          AJORF(I), FFP(I), AJTM(M2-I+1), AJTM(M2+I),             ABSH2565
     *          FPOL(I), RMID(M2-I+1), RMID(M2+I), AISPL(I)             ABSH2566
           READ (NPLOT,99998) PXP(I), PPXP(I)                           ABSH2567
   10 CONTINUE                                                          ABSH2568
C                                                                       ABSH2569
      DO 20 I=1,M2                                                      ABSH2570
           PSIMOM(I) = PSIC0*PSIMOM(I)/(2.*PI)                          ABSH2571
           AJORF(I) = AJORF(I)/STATA                                    ABSH2572
   20 CONTINUE                                                          ABSH2573
C                                                                       ABSH2574
      DO 30 I=1,NRHO                                                    ABSH2575
           READ (NPLOT,99998) RHO(I), SHIFT(I), ELONG(I),               ABSH2576
     *          TRIANG(I)                                               ABSH2577
           RMOM(I) = RHO(I)/A0                                          ABSH2578
   30 CONTINUE                                                          ABSH2579
C                                                                       ABSH2580
      CALL VRSTEC                                                       ABSH2581
      CALL BGNPL(-1)                                                    ABSH2582
      CALL SETDEV(6, 6)                                                 ABSH2583
      CALL NOBRDR                                                       ABSH2584
      CALL CROSS                                                        ABSH2585
      CALL XTICKS(2)                                                    ABSH2586
      CALL YTICKS(2)                                                    ABSH2587
C                                                                       ABSH2588
C  PLOT 1-D CONTOURS                                                    ABSH2589
C                                                                       ABSH2590
C                                                                       ABSH2591
      EOUT = EMOM(M2)                                                   ABSH2592
      X1L = 5.0                                                         ABSH2593
      Y1L = X1L*EOUT                                                    ABSH2594
      CALL TITLE(9HPSI(R,Z)$, +100, 2HR$, 100, 2HZ$, 100, X1L,          ABSH2595
     *     Y1L)                                                         ABSH2596
      RWIDE = 2.0*A0*1.1                                                ABSH2597
      ZWIDE = EOUT*RWIDE/2.0                                            ABSH2598
      R00 = (ARC+DGEOM-1.1)*A0                                          ABSH2599
      CALL GRAF(R00, 5HSCALE, R00+RWIDE, -ZWIDE, 5HSCALE, ZWIDE)        ABSH2600
      IPDEL = M2/8                                                      ABSH2601
      IB = IPDEL + 1                                                    ABSH2602
      DO 50 I=IB,M2,IPDEL                                               ABSH2603
           R0 = SMOM(I)                                                 ABSH2604
           R1 = -XMOM(I)                                                ABSH2605
           R2 = DMOM(I)                                                 ABSH2606
           Z1 = EMOM(I)*XMOM(I)                                         ABSH2607
           Z2 = EMOM(I)*DMOM(I)                                         ABSH2608
           DO 40 J=1,NTHETA                                             ABSH2609
                THETA = (J-1)*360./FLOAT(NTHETA-1)*RADIAN               ABSH2610
                RTHETA(J) = A0*(R0+R1*COS(THETA)+R2*COS(2.*             ABSH2611
     *               THETA))                                            ABSH2612
                ZTHETA(J) = A0*(Z1*SIN(THETA)+Z2*SIN(2.*THETA))         ABSH2613
   40      CONTINUE                                                     ABSH2614
           CALL CURVE(RTHETA, ZTHETA, NTHETA, 0)                        ABSH2615
   50 CONTINUE                                                          ABSH2616
C                                                                       ABSH2617
C                                                                       ABSH2618
C  PLOT ZORNOC CONTOURS                                                 ABSH2619
C  PLOT CONSTANT THETA CONTOURS                                         ABSH2620
C                                                                       ABSH2621
      NTHP = 16                                                         ABSH2622
      NTHPM2 = NTHP - 2                                                 ABSH2623
      DTHETA = TWOPI/(NTHP-1)                                           ABSH2624
      M2P1 = M2 + 1                                                     ABSH2625
      RTHETA(1) = A0*PARAM(1)                                           ABSH2626
      ZTHETA(1) = 0.0                                                   ABSH2627
      DO 70 J=1,NTHPM2                                                  ABSH2628
           TH = J*DTHETA                                                ABSH2629
           DO 60 I=1,M2                                                 ABSH2630
                R0 = SMOM(I)                                            ABSH2631
                R1 = -XMOM(I)                                           ABSH2632
                R2 = DMOM(I)                                            ABSH2633
                Z1 = EMOM(I)*XMOM(I)                                    ABSH2634
                Z2 = EMOM(I)*DMOM(I)                                    ABSH2635
                RTHETA(I+1) = A0*(R0+R1*COS(TH)+R2*COS(2.0*TH))         ABSH2636
                ZTHETA(I+1) = A0*(Z1*SIN(TH)+Z2*SIN(2.0*TH))            ABSH2637
   60      CONTINUE                                                     ABSH2638
           CALL CURVE(RTHETA, ZTHETA, M2P1, 0)                          ABSH2639
   70 CONTINUE                                                          ABSH2640
      CALL ENDGR(0)                                                     ABSH2641
      CALL ENDPL(0)                                                     ABSH2642
C                                                                       ABSH2643
C                                                                       ABSH2644
C                                                                       ABSH2645
      XL = 2.6                                                          ABSH2646
      YL = 2.6                                                          ABSH2647
      X61 = 0.6                                                         ABSH2648
      X62 = 4.3                                                         ABSH2649
      Y61 = 0.5                                                         ABSH2650
      Y62 = 4.0                                                         ABSH2651
      Y63 = 7.5                                                         ABSH2652
      X6L = 3.0                                                         ABSH2653
      Y6L = 2.7                                                         ABSH2654
C                                                                       ABSH2655
C  PLOT THE SHIFT                                                       ABSH2656
C                                                                       ABSH2657
      YMAX = SHIFT(1)                                                   ABSH2658
      YMIN = 0.                                                         ABSH2659
      CALL PHYSOR(X61, Y63)                                             ABSH2660
      CALL TITLE(6HSHIFT$, 100, 2H $, 100, 6HRGEOM$, 100, X6L,          ABSH2661
     *     Y6L)                                                         ABSH2662
      CALL GRAF(0., .2, 1., YMIN, 5HSCALE, YMAX)                        ABSH2663
      CALL CURVE(RMOM, SHIFT, NRHO, 0)                                  ABSH2664
      CALL ENDGR(0)                                                     ABSH2665
C                                                                       ABSH2666
C  PLOT ELONGATION                                                      ABSH2667
C                                                                       ABSH2668
      CALL PHYSOR(X62, Y63)                                             ABSH2669
      CALL TITLE(11HELONGATION$, 100, 2H $, 100, 2HE$, 100,             ABSH2670
     *     X6L, Y6L)                                                    ABSH2671
      YMIN = AMIN1(ELONG(1),ELONG(NRHO))                                ABSH2672
      YMAX = AMAX1(ELONG(1),ELONG(NRHO))                                ABSH2673
      YMINE = (IFIX(YMIN*10.)-1.0)/10.0                                 ABSH2674
      YMAXE = (IFIX(YMAX*10.)+1.0)/10.0                                 ABSH2675
      CALL GRAF(0., .2, 1., YMINE, 0.1, YMAXE)                          ABSH2676
      CALL CURVE(RMOM, ELONG, NRHO, 0)                                  ABSH2677
      CALL ENDGR(0)                                                     ABSH2678
C                                                                       ABSH2679
C  PLOT TRIANGULARITY                                                   ABSH2680
C                                                                       ABSH2681
      CALL PHYSOR(X61, Y62)                                             ABSH2682
      CALL TITLE(15H TRIANGULARITY$, 100, 2H $, 100, 2HD$, 100,         ABSH2683
     *     X6L, Y6L)                                                    ABSH2684
      YMAX = 1.1*TRIANG(NRHO)                                           ABSH2685
      YMIN = 0.                                                         ABSH2686
      CALL GRAF(0., .2, 1., YMIN, 5HSCALE, YMAX)                        ABSH2687
      CALL CURVE(RMOM, TRIANG, NRHO, 0)                                 ABSH2688
      CALL ENDGR(0)                                                     ABSH2689
C                                                                       ABSH2690
C  PLOT VP(X)                                                           ABSH2691
C                                                                       ABSH2692
      CALL PHYSOR(X62, Y62)                                             ABSH2693
      CALL TITLE(7H VP(X)$, 100, 2H $, 100, 3HVP$, 100, X6L,            ABSH2694
     *     Y6L)                                                         ABSH2695
      YMIN = 0.0                                                        ABSH2696
      YMAX = DVDX(1)                                                    ABSH2697
      DO 80 I=2,M2                                                      ABSH2698
           YMAX = AMAX1(YMAX,DVDX(I))                                   ABSH2699
   80 CONTINUE                                                          ABSH2700
      YMAX = 1.1*YMAX                                                   ABSH2701
      CALL GRAF(0., .2, 1., YMIN, 5HSCALE, YMAX)                        ABSH2702
      CALL CURVE(XMOM, DVDX, M2, 0)                                     ABSH2703
      CALL ENDGR(0)                                                     ABSH2704
C                                                                       ABSH2705
C  PLOT Q(X)                                                            ABSH2706
C                                                                       ABSH2707
      CALL PHYSOR(X61, Y61)                                             ABSH2708
      CALL TITLE(6H Q(X)$, 100, 2HX$, 100, 2HQ$, 100, X6L, Y6L)         ABSH2709
      YMIN = 0.                                                         ABSH2710
      YMAX = IFIX(SMQ(M2)) + 1                                          ABSH2711
      CALL YINTAX                                                       ABSH2712
      CALL GRAF(0., .2, 1., YMIN, 1., YMAX)                             ABSH2713
      CALL CURVE(XMOM, SMQ, M2, 0)                                      ABSH2714
      CALL RESET(6HYINTAX)                                              ABSH2715
      CALL ENDGR(0)                                                     ABSH2716
C                                                                       ABSH2717
C  PLOT THE POLOIDAL FLUX                                               ABSH2718
C                                                                       ABSH2719
      CALL PHYSOR(X62, Y61)                                             ABSH2720
      CALL TITLE(8H PSI(X)$, 100, 2HX$, 100, 4HPSI$, 100, X6L,          ABSH2721
     *     Y6L)                                                         ABSH2722
      YMIN = PSIMOM(1)                                                  ABSH2723
      YMAX = 1.1*PSIMOM(M2)                                             ABSH2724
      CALL GRAF(0., .2, 1., YMIN, 5HSCALE, YMAX)                        ABSH2725
      CALL CURVE(XMOM, PSIMOM, M2, 0)                                   ABSH2726
      CALL ENDPL(0)                                                     ABSH2727
C                                                                       ABSH2728
C  PLOT JFLUX                                                           ABSH2729
C                                                                       ABSH2730
      CALL PHYSOR(X61, Y63)                                             ABSH2731
      CALL TITLE(7H JFLUX$, 100, 2HX$, 100, 2HJ$, 100, X6L, Y6L)        ABSH2732
      YMIN = 0.0                                                        ABSH2733
      YMAX = AJORF(1)*1.1                                               ABSH2734
      CALL GRAF(0., .2, 1., YMIN, 5HSCALE, YMAX)                        ABSH2735
      CALL CURVE(XMOM, AJORF, M2, 0)                                    ABSH2736
      CALL ENDGR(0)                                                     ABSH2737
C                                                                       ABSH2738
C  PLOT MID-PLANE J(I,J)                                                ABSH2739
C                                                                       ABSH2740
      CALL PHYSOR(X62, Y63)                                             ABSH2741
      CALL TITLE(6H JMID$, 100, 5HRMID$, 100, 2HJ$, 100, X6L,           ABSH2742
     *     Y6L)                                                         ABSH2743
      M2T = 2*M2                                                        ABSH2744
      XMIN = RMID(1)                                                    ABSH2745
      XMAX = RMID(1)                                                    ABSH2746
      YMIN = AJTM(1)                                                    ABSH2747
      YMAX = AJTM(1)                                                    ABSH2748
      DO 90 I=2,M2T                                                     ABSH2749
           XMIN = AMIN1(XMIN,RMID(I))                                   ABSH2750
           XMAX = AMAX1(XMAX,RMID(I))                                   ABSH2751
           YMAX = AMAX1(YMAX,AJTM(I))                                   ABSH2752
   90 CONTINUE                                                          ABSH2753
      XMIN = XMIN*A0                                                    ABSH2754
      XMAX = XMAX*A0                                                    ABSH2755
      DO 100 I=1,M2T                                                    ABSH2756
           AJTM(I) = AJTM(I)/YMAX                                       ABSH2757
           RMID(I) = RMID(I)*A0                                         ABSH2758
  100 CONTINUE                                                          ABSH2759
      YMIN = 0.                                                         ABSH2760
      YMAX = 1.                                                         ABSH2761
      CALL GRAF(XMIN, 5HSCALE, XMAX, YMIN, 5HSCALE, YMAX)               ABSH2762
      CALL CURVE(RMID, AJTM, M2T, 0)                                    ABSH2763
      GO TO 130                                                         ABSH2764
      DO 110 I=1,M2                                                     ABSH2765
           YRFL(M2+I) = PXP(I)                                          ABSH2766
           YRFL(I) = PXP(M2+1-I)                                        ABSH2767
  110 CONTINUE                                                          ABSH2768
      CALL DOT                                                          ABSH2769
      CALL YGRAXS(YMIN, 5HSCALE, YMAX, Y6L, 2H $, 0, X6L, 0.)           ABSH2770
      CALL CURVE(RMID, YRFL, M2T, 0)                                    ABSH2771
      CALL RESET(3HDOT)                                                 ABSH2772
C                                                                       ABSH2773
      DO 120 I=1,M2                                                     ABSH2774
           YRFL(M2+I) = AJORF(I)/AJORF(1)                               ABSH2775
           YRFL(I) = AJORF(M2+1-I)/AJORF(1)                             ABSH2776
  120 CONTINUE                                                          ABSH2777
      CALL DASH                                                         ABSH2778
      CALL YGRAXS(YMIN, 5HSCALE, YMAX, Y6L, 2H $, 0, X6L, 0.)           ABSH2779
      CALL CURVE(RMID, YRFL, M2T, 0)                                    ABSH2780
      CALL RESET(4HDASH)                                                ABSH2781
  130 CALL ENDGR(0)                                                     ABSH2782
C                                                                       ABSH2783
C  PLOT P(X)                                                            ABSH2784
C                                                                       ABSH2785
      CALL PHYSOR(X61, Y62-.2)                                          ABSH2786
      CALL TITLE(6H P(X)$, 100, 2HX$, 100, 2HP$, 100, X6L, Y6L)         ABSH2787
      YMIN = 0.0                                                        ABSH2788
      YMAX = 1.2                                                        ABSH2789
      CALL GRAF(0., .2, 1., YMIN, 5HSCALE, YMAX)                        ABSH2790
      CALL CURVE(XMOM, PXP, M2, 0)                                      ABSH2791
      GO TO 150                                                         ABSH2792
      YMIN = PPXP(1)                                                    ABSH2793
      YMAX = YMIN                                                       ABSH2794
      DO 140 I=1,M2                                                     ABSH2795
           YMIN = AMIN1(YMIN,PPXP(I))                                   ABSH2796
           YMAX = AMAX1(YMAX,PPXP(I))                                   ABSH2797
  140 CONTINUE                                                          ABSH2798
      CALL YGRAXS(YMIN, 5HSCALE, YMAX, Y6L, 2H $, -100, X6L, 0.)        ABSH2799
      CALL CURVE(XMOM, PPXP, M2, 0)                                     ABSH2800
  150 CALL ENDGR(0)                                                     ABSH2801
C                                                                       ABSH2802
C  PLOT I(X)                                                            ABSH2803
C                                                                       ABSH2804
      CALL PHYSOR(X62, Y62-.2)                                          ABSH2805
      CALL TITLE(6H I(X)$, 100, 2HX$, 100, 2HI$, 100, X6L, Y6L)         ABSH2806
      YMIN = 0.0                                                        ABSH2807
      YMAX = 1.0                                                        ABSH2808
      CALL GRAF(0., .2, 1., YMIN, 5HSCALE, YMAX)                        ABSH2809
      CALL CURVE(XMOM, AISPL, M2, 0)                                    ABSH2810
      CALL ENDGR(0)                                                     ABSH2811
C                                                                       ABSH2812
      CALL ENDPL(0)                                                     ABSH2813
C                                                                       ABSH2814
      CALL DONEPL                                                       ABSH2815
      STOP                                                              ABSH2816
99999 FORMAT (1X, 2I5)                                                  ABSH2817
99998 FORMAT (1X, 10E12.4)                                              ABSH2818
      END                                                               ABSH2819
