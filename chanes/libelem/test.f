c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE TEST(M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/SUBC/   YY(15,15),JJ(15)
      COMMON/SUBS/   S(15,15)
      DOUBLE COMPLEX         S,A(24,24),B(576),Y(1),W,Z,D,YY,JJ
      DOUBLE COMPLEX         ZERO/(0.D0,0.D0)/,Y0/(0.02D0,0.D0)/
      INTEGER         N1(12),N2(12)
      EQUIVALENCE     (A(1,1),B(1)),(B(433),Y(1))

C     PRINT 678
C678  FORMAT(2X,' ####### TEST START ######## ')
      MM1=   M*M
      MM4=MM1/4
      W=(1.D0,0.D0)
      Z=(0.D0,0.D0)
      DO 10 I=1,M
      DO 10 J=1,M
   10 A(I,J)=(0.D0,0.D0)
      DO 15 I=1,576
   15 B(I)=(0.D0,0.D0)
C  FILLING OF THE MATRIX WITH COMPLEX NUMBERS
      DO 17 I=1,M
      DO 17 J=1,M
      A(I,J)=S(I,J)
   17 CONTINUE
C     PRINT 224
C 224 FORMAT(' ======    A OUTPUT.  =======')
C     PRINT 225,((A(I,J),J=1,M),I=1,M)
C 225 FORMAT(1X,4(1X,E13.6,1X,E13.6))
C
C     PRINT 226
C 226 FORMAT(/' ========  B OUTPUT. =======')
C     PRINT 23,(B(I23),I23=1,MM1)
C
      DO 1 J=1,M
      L=M*(J-1)
      DO 1 I=1,M
      K=432+L+I
      B(K)=A(I,J)
C
C     PRINT 240,J,L,I,K,B(K)
C 240 FORMAT(/,2X,'J=',I3,2X,'L=',I3,2X,'I=',I3,2X,'K=',I3,
C    *2X,'B(K)=4/4B',2E13.6)
      IF(I.EQ.J) B(K)=B(K)+W
C
C     PRINT 242,I,L,K,B(K)
C 242 FORMAT(/,2X,'I=',I3,2X,'L=',I3,2X,'K=',I3,
C    *2X,'B(K) ALONG THE DIAGONAL. +1=',2E13.6)
C
      K=288+L+I
      B(K)=-A(I,J)
      IF(I.EQ.J) B(K)=B(K)+W
C
C     PRINT 243,I,L,K,B(K)
C 243 FORMAT(/,2X,'I=',I3,2X,'L=',I3,2X,'K=',I3,
C    *2X,'3/4 B(K) MEASURES THE VALUE=',2E13.6)
C
    1 CONTINUE
C     PRINT 22
C  22 FORMAT('   ======  B 1  =====   ')
C
C     PRINT 23,(B(I23),I23=1,MM1)
C  23 FORMAT(2X,2(2X,E12.5,1X,E12.5))

C     PRINT 230
C 230 FORMAT('     =======   Y BEFORE THE INCREMENT OF YSMINV  =======')
C     PRINT 23,(Y(I23),I23=1,MM4)
C
C ***********************************
      CALL YSMINV(Y,M,D,N1,N2)
C ********************************
C
C     PRINT 221
C 221 FORMAT(4X,'   =======  Y AFTER THE INCREMENT OF YSMINV =====')
C     PRINT 23,(Y(I221),I221=1,MM4)
C
C
      IF(DSQRT(DREAL(D)**2+DIMAG(D)**2).LT.1.D-30) GO TO 444
      DO 2 I=1,M
      DO 2 J=1,M
      D=Z
      DO 3 K=1,M
      IK=288+(K-1)*M+I
      KJ=432+(J-1)*M+K
    3 D=D+B(IK)*B(KJ)
      IJ=144+(J-1)*M+I
    2 B(IJ)=D
C
C
C     PRINT 232
C 232 FORMAT(2X,' =====  B2  ======')
C     PRINT 23,(B(I23),I23=1,MM1)
C
      MM=M*M
      DO 4 I=1,MM
    4 B(288+I)=B(144+I)
C     PRINT 23,(B(I23),I23=1,MM)
      DO 5 J=1,M
      DO 5 I=1,M
    5 A(I,J)=B(288+(J-1)*M+I)
C
C     PRINT 270
C 270 FORMAT(2X,'==========   A RESULT.   ========')
C     PRINT 225,((A(I,J),J=1,M),I=1,M)
      GO TO 500

  444 PRINT 445
  445 FORMAT(2X,'    ABORTION OF THE OPERATION IN THE TEST PROGRAM   ')
C
  500 CONTINUE
C     OVERWRITE MATRIX A WITH MATRIX YY
      DO 460 I=1,M
      DO 460 J=1,M
      YY(I,J)=A(I,J)*Y0
  460 CONTINUE

C       WRITE(6, 9875)
C 9875 FORMAT(2X, ' #######  YY B TEST''E  ######  ')
C       DO 465 II=1,M
C       DO 465 JJJ=1,M
C       WRITE(6, 9876) II,JJJ,YY(II,JJJ)
C 9876 FORMAT(2X,' YY(',I3,',',I3,')=',2(2X,E13.6 ))
C  465 CONTINUE
C

      RETURN
C     DEBUG SUBTRACE
      END

      SUBROUTINE YSMINV(A,N,D,L,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  A(1),L(1),M(1)
      DOUBLE COMPLEX D,A,BIGA,HOLD
      COABS(D)=DSQRT(DREAL(D)**2+DIMAG(D)**2)
      D=(1.D0,0.D0)
      NK=-N
      DO 80 K=1,N
      NK=NK+N
      L(K)=K
      M(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
   10 IF(COABS(BIGA)-COABS(A(IJ)))15,20,20
   15 BIGA=A(IJ)
      L(K)=I
      M(K)=J
   20 CONTINUE
C  INTERCHANGE ROWS
      J=L(K)
      IF(J-K)35,35,25
   25 KI=K-N
      DO 30 I=1,N
      KI=KI+N
      HOLD=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
   30 A(JI)=HOLD
C  INTERCHANGE COLUMNS
   35 I=M(K)
      IF(I-K)45,45,38
   38 JP=N*(I-1)
      DO 40 J=1,N
      JK=NK+J
      JI=JP+J
      HOLD=-A(JK)
      A(JK)=A(JI)
   40 A(JI)=HOLD
   45 IF(COABS(BIGA)-1.D-20) 46,46,48
   46 D=(0.D0,0.D0)
      RETURN
   48 DO 55 I=1,N
      IF(I-K)50,55,50
   50 IK=NK+I
      A(IK)=A(IK)/(-BIGA)
   55 CONTINUE
C  REDUCE MATRIX
      DO 65 I=1,N
      IK=NK+I
      HOLD=A(IK)
      IJ=I-N
      DO 65 J=1,N
      IJ=IJ+N
      IF(I-K)60,65,60
   60 IF(J-K)62,65,62
   62 KJ=IJ-I+K
      A(IJ)=HOLD*A(KJ)+A(IJ)
   65 CONTINUE
      KJ=K-N
      DO 75 J=1,N
      KJ=KJ+N
      IF(J-K)70,75,70
   70 A(KJ)=A(KJ)/BIGA
   75 CONTINUE
      D=D*BIGA
      A(KK)=1.D0/BIGA
   80 CONTINUE
C  FINAL ROW AND COLUMN INTERCHANGE
      K=N
  100 K=K-1
      IF(K)150,150,105
  105 I=L(K)
      IF(I-K)120,120,108
  108 JQ=N*(K-1)
      JR=N*(I-1)
      DO 110 J=1,N
      JK=JQ+J
      HOLD=A(JK)
      JI=JR+J

      A(JK)=-A(JI)
  110 A(JI)=HOLD
  120 J=M(K)
      IF(J-K)100,100,125
  125 KI=K-N
      DO 130 I=1,N
      KI=KI+N

      HOLD=A(KI)
      JI=KI-K+J

      A(KI)=-A(JI)
  130 A(JI)=HOLD
      GO TO 100
  150 RETURN
C     DEBUG SUBTRACE,INIT
      END
