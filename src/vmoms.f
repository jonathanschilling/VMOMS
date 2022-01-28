C***********************************************************************ABSH0001
C*****MDRIV IS THE MAIN DRIVING FOR SOLVING THE MOMENT EQUATIONS OF THE*ABSH0002
C*****GRAD-SHAFRANOV EQUATION USING A VARIATIONAL METHOD.              *ABSH0003
C***********************************************************************ABSH0004
      PROGRAM VMOMS                                                     ABSH0005
      PARAMETER MTHETA=11,MM1=10,MNMOMS=3                               ABSH0006
      PARAMETER MNEQ=6,MNMX=61                                          ABSH0007
      COMMON /NIO/ NIN, NOUT, NPLOT, NDEBUG, NTTY                       ABSH0008
      COMMON /SWORK1/ PRES(MNMX), CURT(MNMX), RHO(MNMX),                ABSH0009
     *     BPRES(MNMX), CPRES(MNMX), DPRES(MNMX), BCURT(MNMX),          ABSH0010
     *     CCURT(MNMX), DCURT(MNMX), NRHO                               ABSH0011
      COMMON /VINPUT/ PA, PB, AIA, AIB, AIC, XE0, XE1, XD0              ABSH0012
      COMMON /VOUPUT/ SHIFT(MNMX), ELONG(MNMX), TRIANG(MNMX)            ABSH0013
      NAMELIST/MOMIN/KMHD,IPARAM,R0,A0,BT0,P0,PA,PB,AI0,AIA,XE0,        ABSH0014
     *XE1,XD0                                                           ABSH0015
      DATA IPARAM /0/, PA /2.0/, PB /1.0/, AIA /2.0/                    ABSH0016
      DATA R0 /94./, A0 /26./, BT0 /9100./, P0 /1.32E+17/, AI0          ABSH0017
     *     /2.E+5/                                                      ABSH0018
      DATA XE0 /1.55/, XE1 /1.0/, XD0 /0.18/                            ABSH0019
      NIN = 10                                                          ABSH0020
      NOUT = 11                                                         ABSH0021
      NPLOT = 12                                                        ABSH0022
      NTTY = 5                                                          ABSH0023
      NDEBUG = 13                                                       ABSH0024
      READ (NIN,MOMIN)                                                  ABSH0025
      NRHO = 25                                                         ABSH0026
      DXX = A0/FLOAT(NRHO-1)                                            ABSH0027
      DO 10 I=1,NRHO                                                    ABSH0028
           RHO(I) = (I-1)*DXX                                           ABSH0029
           X = RHO(I)/A0                                                ABSH0030
           PRES(I) = P0*(1.0-X**PA)**PB                                 ABSH0031
           XSQ = X**2                                                   ABSH0032
           AIB = 3.0 - 2.0*AIA                                          ABSH0033
           AIC = AIA - 2.0                                              ABSH0034
           CURT(I) = AI0*XSQ*(AIA+XSQ*(AIB+XSQ*AIC))                    ABSH0035
   10 CONTINUE                                                          ABSH0036
      IFAIL = 1                                                         ABSH0037
      CALL MHDEQ(KMHD, IPARAM, NRHO, R0, BT0, P0, AI0, XE0,             ABSH0038
     *     XE1, XD0, RHO, SHIFT, ELONG, TRIANG, IFAIL)                  ABSH0039
      CALL GPASMA                                                       ABSH0040
      CALL PRTMEQ                                                       ABSH0041
      CALL PLTMEQ                                                       ABSH0042
      STOP                                                              ABSH0043
      END                                                               ABSH0044
