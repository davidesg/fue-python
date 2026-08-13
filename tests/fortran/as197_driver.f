C     Driver for the transcribed AS 197 listing.  Not part of the paper.
C     Reads a model and a series from standard input and prints what
C     FLIKAM returns, at full precision:
C
C         mp mq n toler
C         p(1..mp)
C         q(1..mq)
C         w(1..n)
C
C     Output: IFAULT, SUMSQ, FACT.  The caller builds the log-likelihood.
      PROGRAM AS197D
      IMPLICIT NONE
      INTEGER MAXN, MAXR
      PARAMETER (MAXN = 5000, MAXR = 64)
      REAL P(MAXR), Q(MAXR), W(MAXN), E(MAXN)
      REAL VW(MAXR), VL(MAXR), VK(MAXR)
      REAL SUMSQ, FACT, TOLER
      INTEGER MP, MQ, N, MR, MRP1, IFAULT, I
C
      READ (*,*) MP, MQ, N, TOLER
      IF (MP .GT. 0) READ (*,*) (P(I), I = 1, MP)
      IF (MQ .GT. 0) READ (*,*) (Q(I), I = 1, MQ)
      READ (*,*) (W(I), I = 1, N)
C
      MR = MAX0(MP, MQ + 1)
      MRP1 = MR + 1
      IFAULT = 0
C
      CALL FLIKAM(P, MP, Q, MQ, W, E, N, SUMSQ, FACT, VW, VL,
     * MRP1, VK, MR, TOLER, IFAULT)
C
      WRITE (*,'(A,I6)')    'IFAULT ', IFAULT
      WRITE (*,'(A,E24.16)') 'SUMSQ  ', SUMSQ
      WRITE (*,'(A,E24.16)') 'FACT   ', FACT
      END

C     MAXO/MINO are the FORTRAN 66 spellings of MAX0/MIN0.  gfortran does
C     not know them, and the transcription is not to be edited, so they
C     are supplied here.  They are exactly integer max and min.
      INTEGER FUNCTION MAXO(I, J)
      INTEGER I, J
      MAXO = MAX0(I, J)
      RETURN
      END

      INTEGER FUNCTION MINO(I, J)
      INTEGER I, J
      MINO = MIN0(I, J)
      RETURN
      END
