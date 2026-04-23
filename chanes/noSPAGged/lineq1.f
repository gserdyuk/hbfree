

      SUBROUTINE LINEQ1(A,NDIM,N,NDIM2,NPRZ,IP,NER,DET,ARGD)
C$LARGE: A,IP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NDIM,NDIM2),IP(NDIM)
      DOUBLE COMPLEX DET
C
C      PRINT 100, N, NDIM
C  100 FORMAT(2X,'LINEQ1 : N=',I5,' NDIM=',I5)

      IP(N)=1
      DO 6 K=1,N
      IF(K.EQ.N) GOTO5
      KP1=K+1
      M=K
      DO 1 I=KP1,N
      IF  (DABS(A(I,K)) .GT. DABS(A(M,K))) M=I
    1 CONTINUE
      IP(K)=M
      IF(M.NE.K) IP(N)=-IP(N)
      T=A(M,K)
      A(M,K)=A(K,K)
      A(K,K)=T
      IF (T.EQ.0.0D0) GO TO 5
      DO 2 I=KP1,N
    2 A(I,K)=-A(I,K)/T
      DO 4 J=KP1,N
      T=A(M,J)
      A(M,J)=A(K,J)
      A(K,J)=T
      IF (T.EQ.0.0D0) GO TO 4
      DO 3 I=KP1,N
    3 A(I,J)=A(I,J)+A(I,K)*T
    4 CONTINUE
    5 IF(DABS(A(K,K)).NE.0.D0) GO TO 6
      IP(N)=0
      GOTO20
    6 CONTINUE
   20 CONTINUE
      NER=0
      DET=(1.0D0,0.0D0)
   15 CONTINUE
      ARGD=1.D0
      IF(NPRZ.LT.1) RETURN
      IF(IP(N).EQ.0)GOTO16
      NP1=N+1
      IF(N.EQ.1)GOTO9
      NM1=N-1
      DO 7 K=1,NM1
      KP1=K+1
      M=IP(K)
      T=A(M,NP1)
      A(M,NP1)=A(K,NP1)
      A(K,NP1)=T
      DO 7 I=KP1,N
    7 A(I,NP1)=A(I,NP1)+A(I,K)*T
      DO 8 KB=1,NM1
      KM1=N-KB
      K=KM1+1
      A(K,NP1)=A(K,NP1)/A(K,K)
      T=-A(K,NP1)
      DO 8 I=1,KM1
    8 A(I,NP1)=A(I,NP1)+A(I,K)*T
    9 A(1,NP1)=A(1,NP1)/A(1,1)
      RETURN
   16 WRITE(6, 17) K
   17 FORMAT(///5X,'LINEQ1....'/15X,I6,'-TH COLUMN CONTAINS ONLY ZE',
     +     'ROES'/)
      NER=1
      RETURN
      END
