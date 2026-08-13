C     ==================================================================
C     ALGORITHM AS 197, TRANSCRIBED FROM THE PUBLISHED ARTICLE
C     ==================================================================
C
C     G. Melard (1984) "Algorithm AS 197: A Fast Algorithm for the Exact
C     Likelihood of Autoregressive-Moving Average Models", Applied
C     Statistics 33(1), 104-114.  The listing is printed in full on
C     pages 110-114 and is transcribed here CHARACTER FOR CHARACTER from
C     `literature/as197.pdf`.  Do not tidy it: the point of this file is
C     that it is the paper and not our code.
C
C       FLIKAM   pages 110-112
C       TWACF    pages 112-114
C
C     Two things are NOT transcription and are recorded here:
C
C       * The listing declares REAL.  The paper's own "Precision" note
C         (page 108) says: "On machines with small word length, all the
C         real variables should be replaced by double precision
C         variables."  Rather than edit the source, the promotion is done
C         by the compiler with -fdefault-real-8, so this file stays
C         verbatim.  See tests/test_as197_published_fortran.py.
C
C       * The conditional-sum-of-squares variant that the paper prints as
C         a comment block inside FLIKAM (page 112) is left as a comment,
C         exactly as printed.  fue computes the exact likelihood.
C
C     Why this file exists: `docs/PROVENANCE.md` could say that
C     `usmelard.c` implements AS 197 and that the C and the Python port
C     agree, but both descend from Mauricio's code.  This is the
C     algorithm as its author published it, which is a different witness.
C     ==================================================================

      SUBROUTINE FLIKAM(P, MP, Q, MQ, W, E, N, SUMSQ, FACT, VW, VL,
     * MRP1, VK, MR, TOLER, IFAULT)
C
C        ALGORITHM AS 197  APPL. STATIST. (1984) VOL.33, NO.1
C
C        COMPUTES THE LIKELIHOOD FUNCTION OF AN AUTOREGRESSIVE-
C        MOVING AVERAGE PROCESS, EXPRESSED AS FACT*SUMSQ
C
      REAL P(MP), Q(MQ), W(N), E(N), VW(MRP1), VL(MRP1), VK(MR)
C
      REAL FACT, SUMSQ, TOLER, EPSIL1, ZERO, P0625, ONE, TWO, FOUR,
     * SIXTEN, A, ALF, AOR, DETCAR, DETMAN, FLJ, FN, R, VL1, VW1
C
      REAL ABS, SQRT
C
      DATA EPSIL1 /1.0E-10/
      DATA ZERO, P0625, ONE, TWO, FOUR, SIXTEN /0.0, 0.0625, 1.0, 2.0,
     * 4.0, 16.0/
C
      FACT = ZERO
      DETMAN = ONE
      DETCAR = ZERO
      SUMSQ = ZERO
      MXPQ = MAXO(MP, MQ)
      MXPQP1 = MXPQ + 1
      MQP1 = MQ + 1
      MPP1 = MP + 1
C
C        CALCULATION OF THE AUTOCOVARIANCE FUNCTION OF A PROCESS WITH
C        UNIT INNOVATION VARIANCE (VW) AND THE COVARIANCES BETWEEN THE
C        VARIABLE AND THE LAGGED INNOVATIONS (VL).
C
      CALL TWACF(P, MP, Q, MQ, VW, MXPQP1, VL, MXPQP1, VK, MXPQ, IFAULT)
      IF (MR .NE. MAXO(MP, MQP1)) IFAULT = 6
      IF (MRP1 .NE. MR + 1) IFAULT = 7
      IF (IFAULT .GT. 0) RETURN
C
C        COMPUTATION OF THE FIRST COLUMN OF MATRIX P (VK)
C
      VK(1) = VW(1)
      IF (MR .EQ. 1) GOTO 150
      DO 140 K = 2, MR
      VK(K) = ZERO
      IF (K .GT. MP) GOTO 120
      DO 110 J = K, MP
      JP2MK = J + 2 - K
      VK(K) = VK(K) + P(J) * VW(JP2MK)
  110 CONTINUE
  120 IF (K .GT. MQP1) GOTO 140
      DO 130 J = K, MQP1
      JP1MK = J + 1 - K
      VK(K) = VK(K) - Q(J - 1) * VL(JP1MK)
  130 CONTINUE
  140 CONTINUE
C
C        COMPUTATION OF THE INITIAL VECTORS L AND K (VL,VK).
C
  150 R = VK(1)
      VL(MR) = ZERO
      DO 160 J = 1, MR
      VW(J) = ZERO
      IF (J .NE. MR) VL(J) = VK(J + 1)
      IF (J .LE. MP) VL(J) = VL(J) + P(J) * R
      VK(J) = VL(J)
  160 CONTINUE
C
C        INITIALIZATION
C
      LAST = MPP1 - MQ
      LOOP = MP
      JFROM = MPP1
      VW(MPP1) = ZERO
      VL(MXPQP1) = ZERO
C
C        EXIT IF NO OBSERVATION, OTHERWISE LOOP ON TIME.
C
      IF (N .LE. 0) GOTO 500
      DO 290 I = 1, N
C
C        TEST FOR SKIPPED UPDATING
C
      IF (I .NE. LAST) GOTO 170
      LOOP = MINO(MP, MQ)
      JFROM = LOOP + 1
C
C        TEST FOR SWITCHING
C
      IF (MQ .LE. 0) GOTO 300
  170 IF (R .LE. EPSIL1) GOTO 400
      IF (ABS(R - ONE) .LT. TOLER .AND. I .GT. MXPQ) GOTO 300
C
C        UPDATING SCALARS
C
      DETMAN = DETMAN * R
  190 IF (ABS(DETMAN) .LT. ONE) GOTO 200
      DETMAN = DETMAN * P0625
      DETCAR = DETCAR + FOUR
      GOTO 190
  200 IF (ABS(DETMAN) .GE. P0625) GOTO 210
      DETMAN = DETMAN * SIXTEN
      DETCAR = DETCAR - FOUR
      GOTO 200
  210 VW1 = VW(1)
      A = W(I) - VW1
      E(I) = A / SQRT(R)
      AOR = A / R
      SUMSQ = SUMSQ + A * AOR
      VL1 = VL(1)
      ALF = VL1 / R
      R = R - ALF * VL1
      IF (LOOP .EQ. 0) GOTO 230
C
C        UPDATING VECTORS
C
      DO 220 J = 1, LOOP
      FLJ = VL(J + 1) + P(J) * VL1
      VW(J) = VW(J + 1) + P(J) * VW1 + AOR * VK(J)
      VL(J) = FLJ - ALF * VK(J)
      VK(J) = VK(J) - ALF * FLJ
  220 CONTINUE
  230 IF (JFROM .GT. MQ) GOTO 250
      DO 240 J = JFROM, MQ
      VW(J) = VW(J + 1) + AOR * VK(J)
      VL(J) = VL(J + 1) - ALF * VK(J)
      VK(J) = VK(J) - ALF * VL(J + 1)
  240 CONTINUE
  250 IF (JFROM .GT. MP) GOTO 270
      DO 260 J = JFROM, MP
  260 VW(J) = VW(J + 1) + P(J) * W(I)
  270 CONTINUE
  290 CONTINUE
      GOTO 390
C
C        QUICK RECURSIONS
C
  300 NEXTI = I
      IFAULT = -NEXTI
      DO 310 I = NEXTI, N
  310 E(I) = W(I)
      IF (MP .EQ. 0) GOTO 340
      DO 330 I = NEXTI, N
      DO 320 J = 1, MP
      IMJ = I - J
      E(I) = E(I) - P(J) * W(IMJ)
  320 CONTINUE
  330 CONTINUE
  340 IF (MQ .EQ. 0) GOTO 370
      DO 360 I = NEXTI, N
      DO 350 J = 1, MQ
      IMJ = I - J
      E(I) = E(I) + Q(J) * E(IMJ)
  350 CONTINUE
  360 CONTINUE
C
C        RETURN SUM OF SQUARES AND DETERMINANT
C
  370 DO 380 I = NEXTI, N
  380 SUMSQ = SUMSQ + E(I) * E(I)
C
C        CODE FOR CONDITIONAL SUM OF SQUARES
C        REPLACES ALL EXECUTABLE STATEMENTS UPTO AND
C        INCLUDING THAT LABELLED 380
C
C                FACT = ZERO
C                DETMAN = ONE
C                DETCAR = ZERO
C                SUMSQ = ZERO
C                MXPQ = MAXO(MP, MQ)
C                DO 380 I=MXPQ,N
C                E(I)=W(I)
C                IF (MP.LE.0) GOTO 340
C                DO 320 J=1,MP
C                IMJ=I-J
C                E(I)=E(I)-P(J)*W(IMJ)
C          320   CONTINUE
C          340   IF (MQ.LE.0) GOTO 380
C                DO 350 J=1,MQ
C                IMJ=I-J
C                E(I)=E(I)+Q(J)*E(IMJ)
C          350   CONTINUE
C          380   SUMSQ=SUMSQ+E(I)*E(I)
C
  390 FN = N
      FACT = DETMAN ** (ONE / FN) * TWO ** (DETCAR / FN)
      RETURN
C
C        EXECUTION ERRORS
C
  400 IFAULT = 8
      RETURN
  500 IFAULT = 9
      RETURN
      END

      SUBROUTINE TWACF(P, MP, Q, MQ, ACF, MA, CVLI, MXPQP1, ALPHA, MXPQ,
     * IFAULT)
C
C        ALGORITHM AS 197.1  APPL. STATIST. (1984) VOL.33, NO.1
C
C        IMPLEMENTATION OF THE ALGORITHM OF G. TUNNICLIFFE WILSON
C        (J. STATIST. COMPUT. SIMUL. 8, 1979, 301-309) FOR THE
C        COMPUTATION OF THE AUTOCOVARIANCE FUNCTION (ACF) OF AN ARMA
C        PROCESS OF ORDER (MP,MQ) AND UNIT INNOVATION VARIANCE.
C        THE AUTOREGRESSIVE AND MOVING AVERAGE COEFFICIENTS ARE STORED
C        IN VECTORS P AND Q, USING BOX AND JENKINS NOTATION. ON OUTPUT
C        VECTOR CVLI CONTAINS THE COVARIANCES BETWEEN THE VARIABLE AND
C        THE (K-1)-LAGGED INNOVATION, FOR K=1,...,MQ+1.
C
      REAL P(MP), Q(MQ), ACF(MA), CVLI(MXPQP1), ALPHA(MXPQ)
C
      REAL EPSIL2, ZERO, HALF, ONE, TWO, DIV
C
      DATA EPSIL2 /1.0E-10/
      DATA ZERO, HALF, ONE, TWO /0.0, 0.5, 1.0, 2.0/
C
      IFAULT = 0
      IF (MP .LT. 0 .OR. MQ .LT. 0) IFAULT = 1
      IF (MXPQ .NE. MAXO(MP, MQ)) IFAULT = 2
      IF (MXPQP1 .NE. MXPQ + 1) IFAULT = 3
      IF (MA .LT. MXPQP1) IFAULT = 4
      IF (IFAULT .GT. 0) RETURN
C
C        INITIALIZATION AND RETURN IF MP=MQ=0
C
      ACF(1) = ONE
      CVLI(1) = ONE
      IF (MA .EQ. 1) RETURN
      DO 10 I = 2, MA
   10 ACF(I) = ZERO
      IF (MXPQP1 .EQ. 1) RETURN
      DO 20 I = 2, MXPQP1
   20 CVLI(I) = ZERO
      DO 90 K = 1, MXPQ
   90 ALPHA(K) = ZERO
C
C        COMPUTATION OF THE A.C.F. OF THE MOVING AVERAGE PART,
C        STORED IN ACF.
C
      IF (MQ .EQ. 0) GOTO 180
      DO 130 K = 1, MQ
      CVLI(K + 1) = -Q(K)
      ACF(K + 1) = -Q(K)
      KC = MQ - K
      IF (KC .EQ. 0) GOTO 120
      DO 110 J = 1, KC
      JPK = J + K
      ACF(K + 1) = ACF(K + 1) + Q(J) * Q(JPK)
  110 CONTINUE
  120 ACF(1) = ACF(1) + Q(K) * Q(K)
  130 CONTINUE
C
C        INITIALIZATION OF CVLI = T.W.-S PHI -- RETURN IF MP=0.
C
  180 IF (MP .EQ. 0) RETURN
      DO 190 K = 1, MP
      ALPHA(K) = P(K)
      CVLI(K) = P(K)
  190 CONTINUE
C
C        COMPUTATION OF T.W.-S ALPHA AND DELTA
C        (DELTA STORED IN ACF WHICH IS GRADUALLY OVERWRITTEN)
C
      DO 290 K = 1, MXPQ
      KC = MXPQ - K
      IF (KC .GE. MP) GOTO 240
      DIV = ONE - ALPHA(KC + 1) * ALPHA(KC + 1)
      IF (DIV .LE. EPSIL2) GOTO 700
      IF (KC .EQ. 0) GOTO 290
      DO 230 J = 1, KC
      KCP1MJ = KC + 1 - J
      ALPHA(J) = (CVLI(J) + ALPHA(KC + 1) * CVLI(KCP1MJ)) / DIV
  230 CONTINUE
  240 IF (KC .GE. MQ) GOTO 260
      J1 = MAXO(KC + 1 - MP, 1)
      DO 250 J = J1, KC
      KCP1MJ = KC + 1 - J
      ACF(J + 1) = ACF(J + 1) + ACF(KC + 2) * ALPHA(KCP1MJ)
  250 CONTINUE
  260 IF (KC .GE. MP) GOTO 290
      DO 270 J = 1, KC
  270 CVLI(J) = ALPHA(J)
  290 CONTINUE
C
C        COMPUTATION OF T.W.-S NU
C        (NU IS STORED IN CVLI, COPIED INTO ACF)
C
      ACF(1) = HALF * ACF(1)
      DO 330 K = 1, MXPQ
      IF (K .GT. MP) GOTO 330
      KP1 = K + 1
      DIV = ONE - ALPHA(K) * ALPHA(K)
      DO 310 J = 1, KP1
      KP2MJ = K + 2 - J
      CVLI(J) = (ACF(J) + ALPHA(K) * ACF(KP2MJ)) / DIV
  310 CONTINUE
      DO 320 J = 1, KP1
  320 ACF(J) = CVLI(J)
  330 CONTINUE
C
C        COMPUTATION OF ACF (ACF IS GRADUALLY OVERWRITTEN)
C
      DO 430 I = 1, MA
      MIIM1P = MINO(I - 1, MP)
      IF (MIIM1P .EQ. 0) GOTO 430
      DO 420 J = 1, MIIM1P
      IMJ = I - J
      ACF(I) = ACF(I) + P(J) * ACF(IMJ)
  420 CONTINUE
  430 CONTINUE
      ACF(1) = ACF(1) * TWO
C
C        COMPUTATION OF CVLI -- RETURN WHEN MQ=0
C
      CVLI(1) = ONE
      IF (MQ .LE. 0) GOTO 600
      DO 530 K = 1, MQ
      CVLI(K + 1) = -Q(K)
      IF (MP .EQ. 0) GOTO 530
      MIKP = MINO(K, MP)
      DO 520 J = 1, MIKP
      KP1MJ = K + 1 - J
      CVLI(K + 1) = CVLI(K + 1) + P(J) * CVLI(KP1MJ)
  520 CONTINUE
  530 CONTINUE
C
  600 RETURN
C
C        EXECUTION ERROR DUE TO (NEAR) NON-STATIONARITY
C
  700 IFAULT = 5
      RETURN
      END
