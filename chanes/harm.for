!*==HARM.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code

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
      SUBROUTINE harm(A,M,Inv,S,Ifset,Iferr)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION A , awi , awr , fn , r , root2 , S , t , theta ,
     &                 w , w2 , w3
      INTEGER i , i2 , i2c , i2cc , i3 , i3c , i3cc , i3ccc , ic , id ,
     &        idif , Iferr , Ifset , igo1 , igo2 , igo3 , il , il1 ,
     &        ilast , Inv
      INTEGER ip1 , ip2 , ip3 , ipp1 , ipp2 , ipp3 , j , jc , jc1 , jd ,
     &        jdif , jj , jj1 , jj2 , jj3 , jjd1 , jjd2 , jjd3 , jjdif ,
     &        jlast
      INTEGER jp1 , jp2 , jp3 , jpp1 , jpp2 , jpp3 , jstep , jstep2 ,
     &        k , k1 , k2 , k3 , kbit , kd , kl , klast , l , lfirst ,
     &        lm1exp , M
      INTEGER m1 , m1mt , m2 , m2mt , m3 , m3mt , mev , mi , minn1 ,
     &        minn2 , minn3 , mt , mtlexp , mtt , n , n1 , n1vnt , n2 ,
     &        n2vnt , n3
      INTEGER n3vnt , np , nt , ntsq , ntv2 , ntvn1 , ntvn2 , ntvn3 , nx
      DIMENSION A(1) , Inv(1) , S(1) , n(3) , M(3) , np(3) , w(2) ,
     &          w2(2) , w3(2)
      EQUIVALENCE (n1,n(1)) , (n2,n(2)) , (n3,n(3))
      SAVE mt , nt
      IF ( iabs(Ifset).LE.1 ) THEN
C
C     THE FOLLOWING PROGRAM COMPUTES THE SIN AND INV TABLES.
C
         mt = max0(M(1),M(2),M(3)) - 2
         mt = max0(2,mt)
         IF ( mt.LE.20 ) THEN
            Iferr = 0
            nt = 2**mt
            ntv2 = nt/2
C
C     SET UP SIN TABLE
C     THETA=PIE/2**(L+1) FOR L=1
            theta = .7853981634D0
C
C     JSTEP=2**(MT-L+1) FOR L=1
            jstep = nt
C
C     JDIF=2**(MT-L) FOR L=1
C     JDIF=2**(MT-L) FOR L=1
            jdif = ntv2
            S(jdif) = dsin(theta)
            DO l = 2 , mt
               theta = theta/2.D0
               jstep2 = jstep
               jstep = jdif
               jdif = jstep/2
               S(jdif) = dsin(theta)
               jc1 = nt - jdif
               S(jc1) = dcos(theta)
               jlast = nt - jstep2
               IF ( jlast.GE.jstep ) THEN
                  DO j = jstep , jlast , jstep
                     jc = nt - j
                     jd = j + jdif
                     S(jd) = S(j)*S(jc1) + S(jdif)*S(jc)
                  ENDDO
               ENDIF
            ENDDO
C
C     SET UP INV(J) TABLE
C
            mtlexp = ntv2
C
C     MTLEXP=2**(MT-L). FOR L=1
            lm1exp = 1
C
C     LM1EXP=2**(L-1). FOR L=1
            Inv(1) = 0
            DO l = 1 , mt
               Inv(lm1exp+1) = mtlexp
               DO j = 2 , lm1exp
                  jj = j + lm1exp
                  Inv(jj) = Inv(j) + mtlexp
               ENDDO
               mtlexp = mtlexp/2
               lm1exp = lm1exp*2
            ENDDO
            IF ( Ifset.EQ.0 ) GOTO 895
         ELSE
            Iferr = 1
            GOTO 895
         ENDIF
      ENDIF
      mtt = max0(M(1),M(2),M(3)) - 2
      root2 = dsqrt(2.D0)
      IF ( mtt.LE.mt ) THEN
         Iferr = 0
         m1 = M(1)
         m2 = M(2)
         m3 = M(3)
         n1 = 2**m1
         n2 = 2**m2
         n3 = 2**m3
         IF ( Ifset.LE.0 ) THEN
            nx = n1*n2*n3
            fn = nx
            DO i = 1 , nx
               A(2*i-1) = A(2*i-1)/fn
               A(2*i) = -A(2*i)/fn
            ENDDO
         ENDIF
         np(1) = n1*2
         np(2) = np(1)*n2
         np(3) = np(2)*n3
         DO id = 1 , 3
            il = np(3) - np(id)
            il1 = il + 1
            mi = M(id)
            IF ( mi.GT.0 ) THEN
               idif = np(id)
               kbit = np(id)
               mev = 2*(mi/2)
               IF ( mi.LE.mev ) THEN
C
C     M IS EVEN
                  lfirst = 2
                  jlast = 0
               ELSE
C
C     M IS ODD. DO L=1 CASE
                  kbit = kbit/2
                  kl = kbit - 2
                  DO i = 1 , il1 , idif
                     klast = kl + i
                     DO k = i , klast , 2
                        kd = k + kbit
C
C     DO ONE STEP WITH L=1,J=0
C     A(K)=A(K)+A(KD)
C     A(KD)=A(K)-A(KD)
C
                        t = A(kd)
                        A(kd) = A(k) - t
                        A(k) = A(k) + t
                        t = A(kd+1)
                        A(kd+1) = A(k+1) - t
                        A(k+1) = A(k+1) + t
                     ENDDO
                  ENDDO
                  IF ( mi.LE.1 ) GOTO 250
                  lfirst = 3
C
C     DEF - JLAST = 2**(L-2) -1
                  jlast = 1
               ENDIF
               DO l = lfirst , mi , 2
                  jjdif = kbit
                  kbit = kbit/4
                  kl = kbit - 2
C
C     DO FOR J=0
                  DO i = 1 , il1 , idif
                     klast = i + kl
                     DO k = i , klast , 2
                        k1 = k + kbit
                        k2 = k1 + kbit
                        k3 = k2 + kbit
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
                        t = A(k2)
                        A(k2) = A(k) - t
                        A(k) = A(k) + t
                        t = A(k2+1)
                        A(k2+1) = A(k+1) - t
                        A(k+1) = A(k+1) + t
C
                        t = A(k3)
                        A(k3) = A(k1) - t
                        A(k1) = A(k1) + t
                        t = A(k3+1)
                        A(k3+1) = A(k1+1) - t
                        A(k1+1) = A(k1+1) + t
C
                        t = A(k1)
                        A(k1) = A(k) - t
                        A(k) = A(k) + t
                        t = A(k1+1)
                        A(k1+1) = A(k+1) - t
                        A(k+1) = A(k+1) + t
C
                        r = -A(k3+1)
                        t = A(k3)
                        A(k3) = A(k2) - r
                        A(k2) = A(k2) + r
                        A(k3+1) = A(k2+1) - t
                        A(k2+1) = A(k2+1) + t
                     ENDDO
                  ENDDO
                  IF ( jlast.GT.0 ) THEN
                     jj = jjdif + 1
C
C     DO FOR J=1
                     ilast = il + jj
                     DO i = jj , ilast , idif
                        klast = kl + i
                        DO k = i , klast , 2
                           k1 = k + kbit
                           k2 = k1 + kbit
                           k3 = k2 + kbit
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
                           r = -A(k2+1)
                           t = A(k2)
                           A(k2) = A(k) - r
                           A(k) = A(k) + r
                           A(k2+1) = A(k+1) - t
                           A(k+1) = A(k+1) + t
C
                           awr = A(k1) - A(k1+1)
                           awi = A(k1+1) + A(k1)
                           r = -A(k3) - A(k3+1)
                           t = A(k3) - A(k3+1)
                           A(k3) = (awr-r)/root2
                           A(k3+1) = (awi-t)/root2
                           A(k1) = (awr+r)/root2
                           A(k1+1) = (awi+t)/root2
                           t = A(k1)
                           A(k1) = A(k) - t
                           A(k) = A(k) + t
                           t = A(k1+1)
                           A(k1+1) = A(k+1) - t
                           A(k+1) = A(k+1) + t
                           r = -A(k3+1)
                           t = A(k3)
                           A(k3) = A(k2) - r
                           A(k2) = A(k2) + r
                           A(k3+1) = A(k2+1) - t
                           A(k2+1) = A(k2+1) + t
                        ENDDO
                     ENDDO
                     IF ( jlast.GT.1 ) THEN
                        jj = jj + jjdif
C
C     NOW DO THE REMAINING J'S
                        DO j = 2 , jlast
C
C     FETCH W'S
C     DEF- W=W**INV(J), W2=W**2, W3=W**3
                           i = Inv(j+1)
                           ic = nt - i
                           w(1) = S(ic)
                           w(2) = S(i)
                           i2 = 2*i
                           i2c = nt - i2
                           IF ( i2c.LT.0 ) THEN
C
C     2*I IS IN SECOND QUADRANT
                              i2cc = i2c + nt
                              i2c = -i2c
                              w2(1) = -S(i2c)
                              w2(2) = S(i2cc)
                           ELSEIF ( i2c.EQ.0 ) THEN
                              w2(1) = 0.D0
                              w2(2) = 1.D0
                           ELSE
C
C     2*I IS IN FIRST QUADRANT
                              w2(1) = S(i2c)
                              w2(2) = S(i2)
                           ENDIF
                           i3 = i + i2
                           i3c = nt - i3
                           IF ( i3c.LT.0 ) THEN
C
                              i3cc = i3c + nt
                              IF ( i3cc.LT.0 ) THEN
C
C     3*I IN THIRD QUADRANT
                                 i3ccc = nt + i3cc
                                 i3cc = -i3cc
                                 w3(1) = -S(i3ccc)
                                 w3(2) = -S(i3cc)
                              ELSEIF ( i3cc.EQ.0 ) THEN
                                 w3(1) = -1.D0
                                 w3(2) = 0.D0
                              ELSE
C
C     I3 IN SECOND QUADRANT
                                 i3c = -i3c
                                 w3(1) = -S(i3c)
                                 w3(2) = S(i3cc)
                              ENDIF
                           ELSEIF ( i3c.EQ.0 ) THEN
                              w3(1) = 0.D0
                              w3(2) = 1.D0
                           ELSE
C
C     I3 IN FIRST QUADRANT
                              w3(1) = S(i3c)
                              w3(2) = S(i3)
                           ENDIF
                           ilast = il + jj
                           DO i = jj , ilast , idif
                              klast = kl + i
                              DO k = i , klast , 2
                                 k1 = k + kbit
                                 k2 = k1 + kbit
                                 k3 = k2 + kbit
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
                                 r = A(k2)*w2(1) - A(k2+1)*w2(2)
                                 t = A(k2)*w2(2) + A(k2+1)*w2(1)
                                 A(k2) = A(k) - r
                                 A(k) = A(k) + r
                                 A(k2+1) = A(k+1) - t
                                 A(k+1) = A(k+1) + t
C
                                 r = A(k3)*w3(1) - A(k3+1)*w3(2)
                                 t = A(k3)*w3(2) + A(k3+1)*w3(1)
                                 awr = A(k1)*w(1) - A(k1+1)*w(2)
                                 awi = A(k1)*w(2) + A(k1+1)*w(1)
                                 A(k3) = awr - r
                                 A(k3+1) = awi - t
                                 A(k1) = awr + r
                                 A(k1+1) = awi + t
                                 t = A(k1)
                                 A(k1) = A(k) - t
                                 A(k) = A(k) + t
                                 t = A(k1+1)
                                 A(k1+1) = A(k+1) - t
                                 A(k+1) = A(k+1) + t
                                 r = -A(k3+1)
                                 t = A(k3)
                                 A(k3) = A(k2) - r
                                 A(k2) = A(k2) + r
                                 A(k3+1) = A(k2+1) - t
                                 A(k2+1) = A(k2+1) + t
                              ENDDO
                           ENDDO
C     END OF I AND K LOOPS
C
                           jj = jjdif + jj
                        ENDDO
                     ENDIF
                  ENDIF
C     END OF J-LOOP
C
                  jlast = 4*jlast + 3
               ENDDO
            ENDIF
C     END OF  L  LOOP
C
 250     ENDDO
C     END OF  ID  LOOP
C
C     WE NOW HAVE THE COMPLEX FOURIER SUMS BUT THEIR ADDRESSES ARE
C     BIT-REVERSED.  THE FOLLOWING ROUTINE PUTS THEM IN ORDER
         ntsq = nt*nt
         m3mt = m3 - mt
         IF ( m3mt.LT.0 ) THEN
C
C     M3 LESS THAN MT
            igo3 = 2
            n3vnt = 1
            ntvn3 = nt/n3
            minn3 = n3
         ELSE
C
C     M3 GR. OR EQ. MT
            igo3 = 1
            n3vnt = n3/nt
            minn3 = nt
         ENDIF
         jjd3 = ntsq/n3
         m2mt = m2 - mt
         IF ( m2mt.LT.0 ) THEN
C
C     M2 LESS THAN MT
            igo2 = 2
            n2vnt = 1
            ntvn2 = nt/n2
            minn2 = n2
         ELSE
C
C     M2 GR. OR EQ. MT
            igo2 = 1
            n2vnt = n2/nt
            minn2 = nt
         ENDIF
         jjd2 = ntsq/n2
         m1mt = m1 - mt
         IF ( m1mt.LT.0 ) THEN
C
C     M1 LESS THAN MT
            igo1 = 2
            n1vnt = 1
            ntvn1 = nt/n1
            minn1 = n1
         ELSE
C
C     M1 GR. OR EQ. MT
            igo1 = 1
            n1vnt = n1/nt
            minn1 = nt
         ENDIF
         jjd1 = ntsq/n1
         jj3 = 1
         j = 1
         DO jpp3 = 1 , n3vnt
            ipp3 = Inv(jj3)
            DO jp3 = 1 , minn3
               IF ( igo3.EQ.2 ) THEN
                  ip3 = Inv(jp3)/ntvn3
               ELSE
                  ip3 = Inv(jp3)*n3vnt
               ENDIF
               i3 = (ipp3+ip3)*n2
               jj2 = 1
               DO jpp2 = 1 , n2vnt
                  ipp2 = Inv(jj2) + i3
                  DO jp2 = 1 , minn2
                     IF ( igo2.EQ.2 ) THEN
                        ip2 = Inv(jp2)/ntvn2
                     ELSE
                        ip2 = Inv(jp2)*n2vnt
                     ENDIF
                     i2 = (ipp2+ip2)*n1
                     jj1 = 1
                     DO jpp1 = 1 , n1vnt
                        ipp1 = Inv(jj1) + i2
                        DO jp1 = 1 , minn1
                           IF ( igo1.EQ.2 ) THEN
                              ip1 = Inv(jp1)/ntvn1
                           ELSE
                              ip1 = Inv(jp1)*n1vnt
                           ENDIF
                           i = 2*(ipp1+ip1) + 1
                           IF ( j.LT.i ) THEN
                              t = A(i)
                              A(i) = A(j)
                              A(j) = t
                              t = A(i+1)
                              A(i+1) = A(j+1)
                              A(j+1) = t
                           ENDIF
                           j = j + 2
                        ENDDO
                        jj1 = jj1 + jjd1
                     ENDDO
                  ENDDO
C     END OF JPP1 AND JP2
C
                  jj2 = jj2 + jjd2
               ENDDO
            ENDDO
C     END OF JPP2 AND JP3 LOOPS
C
            jj3 = jj3 + jjd3
         ENDDO
C     END OF JPP3 LOOP
C
         IF ( Ifset.LT.0 ) THEN
            DO i = 1 , nx
               A(2*i) = -A(2*i)
            ENDDO
         ENDIF
      ELSE
         Iferr = 1
         RETURN
      ENDIF
 895  RETURN
      END
