
C
C     ..................................................................
C
C        SUBROUTINE HARM
C
C        PURPOSE
C           PERFORMS DISCRETE COMPLEX FOURIER TRANSFORMS ON A COMPLEX
C           THREE DIMENSIONAL ARRAY
C
C        USAGE
C           CALL HARM (A,M,INV,S,IFSET,IFERR)
C
C        DESCRIPTION OF PARAMETERS
C           A     - AS INPUT, A CONTAINS THE COMPLEX, 3-DIMENSIONAL
C                   ARRAY TO BE TRANSFORMED.  THE REAL PART OF
C                   A(I1,I2,I3) IS STORED IN VECTOR FASHION IN A CELL
C                   WITH INDEX 2*(I3*N1*N2 + I2*N1 + I1) + 1 WHERE
C                   NI = 2**M(I), I=1,2,3 AND I1 = 0,1,...,N1-1 ETC.
C                   THE IMAGINARY PART IS IN THE CELL IMMEDIATELY
C                   FOLLOWING.  NOTE THAT THE SUBSCRIPT I1 INCREASES
C                   MOST RAPIDLY AND I3 INCREASES LEAST RAPIDLY.
C                   AS OUTPUT, A CONTAINS THE COMPLEX FOURIER
C                   TRANSFORM.  THE NUMBER OF CORE LOCATIONS OF
C                   ARRAY A IS 2*(N1*N2*N3)
C           M     - A THREE CELL VECTOR WHICH DETERMINES THE SIZES
C                   OF THE 3 DIMENSIONS OF THE ARRAY A.   THE SIZE,
C                   NI, OF THE I DIMENSION OF A IS 2**M(I), I = 1,2,3
C           INV   - A VECTOR WORK AREA FOR BIT AND INDEX MANIPULATION
C                   OF DIMENSION ONE EIGHTH THE NUMBER OF CORE
C                   LOCATIONS OF A, VIZ., (1/8)*2*N1*N2*N3
C           S     - A VECTOR WORK AREA FOR SINE TABLES WITH DIMENSION
C                   THE SAME AS INV
C           IFSET - AN OPTION PARAMETER WITH THE FOLLOWING SETTINGS
C                      0    SET UP SINE AND INV TABLES ONLY
C                      1    SET UP SINE AND INV TABLES ONLY AND
C                           CALCULATE FOURIER TRANSFORM
C                     -1    SET UP SINE AND INV TABLES ONLY AND
C                           CALCULATE INVERSE FOURIER TRANSFORM (FOR
C                           THE MEANING OF INVERSE SEE THE EQUATIONS
C                           UNDER METHOD BELOW)
C                      2    CALCULATE FOURIER TRANSFORM ONLY (ASSUME
C                           SINE AND INV TABLES EXIST)
C                     -2    CALCULATE INVERSE FOURIER TRANSFORM ONLY
C                           (ASSUME SINE AND INV TABLES EXIST)
C           IFERR - ERROR INDICATOR.   WHEN IFSET IS 0,+1,-1,
C                   IFERR = 1 MEANS THE MAXIMUM M(I) IS GREATER THAN
C                  20 , I=1,2,3   WHEN IFSET IS 2,-2 , IFERR = 1
C                   MEANS THAT THE SINE AND INV TABLES ARE NOT LARGE
C                   ENOUGH OR HAVE NOT BEEN COMPUTED .
C                   IF ON RETURN IFERR = 0 THEN NONE OF THE ABOVE
C                   CONDITIONS ARE PRESENT
C
C        REMARKS
C           THIS SUBROUTINE IS TO BE USED FOR COMPLEX, 3-DIMENSIONAL
C           ARRAYS IN WHICH EACH DIMENSION IS A POWER OF 2.  THE
C           MAXIMUM M(I) MUST NOT BE LESS THAN 3 OR GREATER THAN 20,
C           I = 1,2,3
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C         METHOD
C           FOR IFSET = +1, OR +2, THE FOURIER TRANSFORM OF COMPLEX
C           ARRAY A IS OBTAINED.
C
C                  N1-1   N2-1   N3-1                L1   L2   L3
C     X(J1,J2,J3)=SUM    SUM    SUM    A(K1,K2,K3)*W1  *W2  *W3
C                  K1=0   K2=0   K3=0
C
C                  WHERE WI IS THE N(I) ROOT OF UNITY AND L1=K1*J1,
C                        L2=K2*J2, L3=K3*J3
C
C
C           FOR IFSET = -1, OR -2, THE INVERSE FOURIER TRANSFORM A OF
C           COMPLEX ARRAY X IS OBTAINED.
C
C     A(K1,K2,K3)=
C               1      N1-1   N2-1   N3-1                -L1  -L2  -L3
C           -------- *SUM    SUM    SUM    X(J1,J2,J3)*W1  *W2  *W3
C           N1*N2*N3   J1=0   J2=0   J3=0
C
C
C           SEE J.W. COOLEY AND J.W. TUKEY, 'AN ALGORITHM FOR THE
C           MACHINE CALCULATION OF COMPLEX FOURIER SERIES',
C           MATHEMATICS OF COMPUTATIONS, VOL. 19 (APR. 1965), P. 297.
C
C     ..................................................................
C
      SUBROUTINE HARM(A,M,INV,S,IFSET, IFERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1),INV(1),S(1),N(3),M(3),NP(3),W(2),W2(2),W3(2)
      EQUIVALENCE (N1,N(1)),(N2,N(2)),(N3,N(3))
      SAVE MT,NT,MTV2
   10 IF( IABS(IFSET) - 1) 900,900,12
   12 MTT=MAX0(M(1),M(2),M(3)) -2
      ROOT2 = DSQRT(2.D0)
      IF (MTT-MT ) 14,14,13
   13 IFERR=1
      RETURN
   14 IFERR=0
      M1=M(1)
      M2=M(2)
      M3=M(3)
      N1=2**M1
      N2=2**M2
      N3=2**M3
   16 IF(IFSET) 18,18,20
   18 NX= N1*N2*N3
      FN = NX
      DO 19 I = 1,NX
      A(2*I-1) = A(2*I-1)/FN
   19 A(2*I) = -A(2*I)/FN
   20 NP(1)=N1*2
      NP(2)= NP(1)*N2
      NP(3)=NP(2)*N3
      DO 250 ID=1,3
      IL = NP(3)-NP(ID)
      IL1 = IL+1
      MI = M(ID)
      IF (MI)250,250,30
   30 IDIF=NP(ID)
      KBIT=NP(ID)
      MEV = 2*(MI/2)
      IF (MI - MEV )60,60,40
C
C     M IS ODD. DO L=1 CASE
   40 KBIT=KBIT/2
      KL=KBIT-2
      DO 50 I=1,IL1,IDIF
      KLAST=KL+I
      DO 50 K=I,KLAST,2
      KD=K+KBIT
C
C     DO ONE STEP WITH L=1,J=0
C     A(K)=A(K)+A(KD)
C     A(KD)=A(K)-A(KD)
C
      T=A(KD)
      A(KD)=A(K)-T
      A(K)=A(K)+T
      T=A(KD+1)
      A(KD+1)=A(K+1)-T
   50 A(K+1)=A(K+1)+T
      IF (MI - 1)250,250,52
   52 LFIRST =3
C
C     DEF - JLAST = 2**(L-2) -1
      JLAST=1
      GO TO 70
C
C     M IS EVEN
   60 LFIRST = 2
      JLAST=0
   70 DO 240 L=LFIRST,MI,2
      JJDIF=KBIT
      KBIT=KBIT/4
      KL=KBIT-2
C
C     DO FOR J=0
      DO 80 I=1,IL1,IDIF
      KLAST=I+KL
      DO 80 K=I,KLAST,2
      K1=K+KBIT
      K2=K1+KBIT
      K3=K2+KBIT
C
C     DO TWO STEPS WITH J=0
C     A(K)=A(K)+A(K2)
C     A(K2)=A(K)-A(K2)
C     A(K1)=A(K1)+A(K3)
C     A(K3)=A(K1)-A(K3)
C
C     A(K)=A(K)+A(K1)
C     A(K1)=A(K)-A(K1)
C     A(K2)=A(K2)+A(K3)*I
C     A(K3)=A(K2)-A(K3)*I
C
      T=A(K2)
      A(K2)=A(K)-T
      A(K)=A(K)+T
      T=A(K2+1)
      A(K2+1)=A(K+1)-T
      A(K+1)=A(K+1)+T
C
      T=A(K3)
      A(K3)=A(K1)-T
      A(K1)=A(K1)+T
      T=A(K3+1)
      A(K3+1)=A(K1+1)-T
      A(K1+1)=A(K1+1)+T
C
      T=A(K1)
      A(K1)=A(K)-T
      A(K)=A(K)+T
      T=A(K1+1)
      A(K1+1)=A(K+1)-T
      A(K+1)=A(K+1)+T
C
      R=-A(K3+1)
      T = A(K3)
      A(K3)=A(K2)-R
      A(K2)=A(K2)+R
      A(K3+1)=A(K2+1)-T
   80 A(K2+1)=A(K2+1)+T
      IF (JLAST) 235,235,82
   82 JJ=JJDIF   +1
C
C     DO FOR J=1
      ILAST= IL +JJ
      DO 85 I = JJ,ILAST,IDIF
      KLAST = KL+I
      DO 85 K=I,KLAST,2
      K1 = K+KBIT
      K2 = K1+KBIT
      K3 = K2+KBIT
C
C     LETTING W=(1+I)/ROOT2,W3=(-1+I)/ROOT2,W2=I,
C     A(K)=A(K)+A(K2)*I
C     A(K2)=A(K)-A(K2)*I
C     A(K1)=A(K1)*W+A(K3)*W3
C     A(K3)=A(K1)*W-A(K3)*W3
C
C     A(K)=A(K)+A(K1)
C     A(K1)=A(K)-A(K1)
C     A(K2)=A(K2)+A(K3)*I
C     A(K3)=A(K2)-A(K3)*I
C
      R =-A(K2+1)
      T = A(K2)
      A(K2) = A(K)-R
      A(K) = A(K)+R
      A(K2+1)=A(K+1)-T
      A(K+1)=A(K+1)+T
C
      AWR=A(K1)-A(K1+1)
      AWI = A(K1+1)+A(K1)
      R=-A(K3)-A(K3+1)
      T=A(K3)-A(K3+1)
      A(K3)=(AWR-R)/ROOT2
      A(K3+1)=(AWI-T)/ROOT2
      A(K1)=(AWR+R)/ROOT2
      A(K1+1)=(AWI+T)/ROOT2
      T= A(K1)
      A(K1)=A(K)-T
      A(K)=A(K)+T
      T=A(K1+1)
      A(K1+1)=A(K+1)-T
      A(K+1)=A(K+1)+T
      R=-A(K3+1)
      T=A(K3)
      A(K3)=A(K2)-R
      A(K2)=A(K2)+R
      A(K3+1)=A(K2+1)-T
   85 A(K2+1)=A(K2+1)+T
      IF(JLAST-1) 235,235,90
   90 JJ= JJ + JJDIF
C
C     NOW DO THE REMAINING J'S
      DO 230 J=2,JLAST
C
C     FETCH W'S
C     DEF- W=W**INV(J), W2=W**2, W3=W**3
   96 I=INV(J+1)
   98 IC=NT-I
      W(1)=S(IC)
      W(2)=S(I)
      I2=2*I
      I2C=NT-I2
      IF(I2C)120,110,100
C
C     2*I IS IN FIRST QUADRANT
  100 W2(1)=S(I2C)
      W2(2)=S(I2)
      GO TO 130
  110 W2(1)=0.D0
      W2(2)=1.D0
      GO TO 130
C
C     2*I IS IN SECOND QUADRANT
  120 I2CC = I2C+NT
      I2C=-I2C
      W2(1)=-S(I2C)
      W2(2)=S(I2CC)
  130 I3=I+I2
      I3C=NT-I3
      IF(I3C)160,150,140
C
C     I3 IN FIRST QUADRANT
  140 W3(1)=S(I3C)
      W3(2)=S(I3)
      GO TO 200
  150 W3(1)=0.D0
      W3(2)=1.D0
      GO TO 200
C
  160 I3CC=I3C+NT
      IF(I3CC)190,180,170
C
C     I3 IN SECOND QUADRANT
  170 I3C=-I3C
      W3(1)=-S(I3C)
      W3(2)=S(I3CC)
      GO TO 200
  180 W3(1)=-1.D0
      W3(2)=0.D0
      GO TO 200
C
C     3*I IN THIRD QUADRANT
  190 I3CCC=NT+I3CC
      I3CC = -I3CC
      W3(1)=-S(I3CCC)
      W3(2)=-S(I3CC)
  200 ILAST=IL+JJ
      DO 220 I=JJ,ILAST,IDIF
      KLAST=KL+I
      DO 220 K=I,KLAST,2
      K1=K+KBIT
      K2=K1+KBIT
      K3=K2+KBIT
C
C     DO TWO STEPS WITH J NOT 0
C     A(K)=A(K)+A(K2)*W2
C     A(K2)=A(K)-A(K2)*W2
C     A(K1)=A(K1)*W+A(K3)*W3
C     A(K3)=A(K1)*W-A(K3)*W3
C
C     A(K)=A(K)+A(K1)
C     A(K1)=A(K)-A(K1)
C     A(K2)=A(K2)+A(K3)*I
C     A(K3)=A(K2)-A(K3)*I
C
      R=A(K2)*W2(1)-A(K2+1)*W2(2)
      T=A(K2)*W2(2)+A(K2+1)*W2(1)
      A(K2)=A(K)-R
      A(K)=A(K)+R
      A(K2+1)=A(K+1)-T
      A(K+1)=A(K+1)+T
C
      R=A(K3)*W3(1)-A(K3+1)*W3(2)
      T=A(K3)*W3(2)+A(K3+1)*W3(1)
      AWR=A(K1)*W(1)-A(K1+1)*W(2)
      AWI=A(K1)*W(2)+A(K1+1)*W(1)
      A(K3)=AWR-R
      A(K3+1)=AWI-T
      A(K1)=AWR+R
      A(K1+1)=AWI+T
      T=A(K1)
      A(K1)=A(K)-T
      A(K)=A(K)+T
      T=A(K1+1)
      A(K1+1)=A(K+1)-T
      A(K+1)=A(K+1)+T
      R=-A(K3+1)
      T=A(K3)
      A(K3)=A(K2)-R
      A(K2)=A(K2)+R
      A(K3+1)=A(K2+1)-T
  220 A(K2+1)=A(K2+1)+T
C     END OF I AND K LOOPS
C
  230 JJ=JJDIF+JJ
C     END OF J-LOOP
C
  235 JLAST=4*JLAST+3
  240 CONTINUE
C     END OF  L  LOOP
C
  250 CONTINUE
C     END OF  ID  LOOP
C
C     WE NOW HAVE THE COMPLEX FOURIER SUMS BUT THEIR ADDRESSES ARE
C     BIT-REVERSED.  THE FOLLOWING ROUTINE PUTS THEM IN ORDER
      NTSQ=NT*NT
      M3MT=M3-MT
  350 IF(M3MT) 370,360,360
C
C     M3 GR. OR EQ. MT
  360 IGO3=1
      N3VNT=N3/NT
      MINN3=NT
      GO TO 380
C
C     M3 LESS THAN MT
  370 IGO3=2
      N3VNT=1
      NTVN3=NT/N3
      MINN3=N3
  380 JJD3 = NTSQ/N3
      M2MT=M2-MT
  450 IF (M2MT)470,460,460
C
C     M2 GR. OR EQ. MT
  460 IGO2=1
      N2VNT=N2/NT
      MINN2=NT
      GO TO 480
C
C     M2 LESS THAN MT
  470 IGO2 = 2
      N2VNT=1
      NTVN2=NT/N2
      MINN2=N2
  480 JJD2=NTSQ/N2
      M1MT=M1-MT
  550 IF(M1MT)570,560,560
C
C     M1 GR. OR EQ. MT
  560 IGO1=1
      N1VNT=N1/NT
      MINN1=NT
      GO TO 580
C
C     M1 LESS THAN MT
  570 IGO1=2
      N1VNT=1
      NTVN1=NT/N1
      MINN1=N1
  580 JJD1=NTSQ/N1
  600 JJ3=1
      J=1
      DO 880 JPP3=1,N3VNT
      IPP3=INV(JJ3)
      DO 870 JP3=1,MINN3
      GO TO (610,620),IGO3
  610 IP3=INV(JP3)*N3VNT
      GO TO 630
  620 IP3=INV(JP3)/NTVN3
  630 I3=(IPP3+IP3)*N2
  700 JJ2=1
      DO 870 JPP2=1,N2VNT
      IPP2=INV(JJ2)+I3
      DO 860 JP2=1,MINN2
      GO TO (710,720),IGO2
  710 IP2=INV(JP2)*N2VNT
      GO TO 730
  720 IP2=INV(JP2)/NTVN2
  730 I2=(IPP2+IP2)*N1
  800 JJ1=1
      DO 860 JPP1=1,N1VNT
      IPP1=INV(JJ1)+I2
      DO 850 JP1=1,MINN1
      GO TO (810,820),IGO1
  810 IP1=INV(JP1)*N1VNT
      GO TO 830
  820 IP1=INV(JP1)/NTVN1
  830 I=2*(IPP1+IP1)+1
      IF (J-I) 840,845,845
  840 T=A(I)
      A(I)=A(J)
      A(J)=T
      T=A(I+1)
      A(I+1)=A(J+1)
      A(J+1)=T
  845 CONTINUE
  850 J=J+2
  860 JJ1=JJ1+JJD1
C     END OF JPP1 AND JP2
C
  870 JJ2=JJ2+JJD2
C     END OF JPP2 AND JP3 LOOPS
C
  880 JJ3 = JJ3+JJD3
C     END OF JPP3 LOOP
C
  890 IF(IFSET)891,895,895
  891 DO 892 I = 1,NX
  892 A(2*I) = -A(2*I)
  895 RETURN
C
C     THE FOLLOWING PROGRAM COMPUTES THE SIN AND INV TABLES.
C
  900 MT=MAX0(M(1),M(2),M(3)) -2
      MT = MAX0(2,MT)
  904 IF (MT-20)906,906,905
  905 IFERR = 1
      GO TO 895
  906 IFERR=0
      NT=2**MT
      NTV2=NT/2
C
C     SET UP SIN TABLE
C     THETA=PIE/2**(L+1) FOR L=1
  910 THETA=.7853981634D0
C
C     JSTEP=2**(MT-L+1) FOR L=1
      JSTEP=NT
C
C     JDIF=2**(MT-L) FOR L=1
C     JDIF=2**(MT-L) FOR L=1
      JDIF=NTV2
      S(JDIF)=DSIN(THETA)
      DO 950 L=2,MT
      THETA=THETA/2.D0
      JSTEP2=JSTEP
      JSTEP=JDIF
      JDIF=JSTEP/2
      S(JDIF)=DSIN(THETA)
      JC1=NT-JDIF
      S(JC1)=DCOS(THETA)
      JLAST=NT-JSTEP2
      IF(JLAST - JSTEP) 950,920,920
  920 DO 940 J=JSTEP,JLAST,JSTEP
      JC=NT-J
      JD=J+JDIF
  940 S(JD)=S(J)*S(JC1)+S(JDIF)*S(JC)
  950 CONTINUE
C
C     SET UP INV(J) TABLE
C
  960 MTLEXP=NTV2
C
C     MTLEXP=2**(MT-L). FOR L=1
      LM1EXP=1
C
C     LM1EXP=2**(L-1). FOR L=1
      INV(1)=0
      DO 980 L=1,MT
      INV(LM1EXP+1) = MTLEXP
      DO 970 J=2,LM1EXP
      JJ=J+LM1EXP
  970 INV(JJ)=INV(J)+MTLEXP
      MTLEXP=MTLEXP/2
  980 LM1EXP=LM1EXP*2
  982 IF(IFSET)12,895,12
      END
