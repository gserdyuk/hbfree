!*==TEST.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE test(M)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , ij , ik , j , k , kj , l , M , mm , mm1 , mm4
      COMMON /subc  / Yy(15,15) , Jj(15)
      COMMON /subs  / S(15,15)
      DOUBLE COMPLEX S , a(24,24) , b(576) , y(1) , w , z , d , Yy , Jj
      DOUBLE COMPLEX y0/(0.02D0,0.D0)/
      INTEGER n1(12) , n2(12)
      EQUIVALENCE (a(1,1),b(1)) , (b(433),y(1))

C     PRINT 678
C678  FORMAT(2X,' ####### TEST START ######## ')
      mm1 = M*M
      mm4 = mm1/4
      w = (1.D0,0.D0)
      z = (0.D0,0.D0)
      DO i = 1 , M
         DO j = 1 , M
            a(i,j) = (0.D0,0.D0)
         ENDDO
      ENDDO
      DO i = 1 , 576
         b(i) = (0.D0,0.D0)
      ENDDO
C  FILLING OF THE MATRIX WITH COMPLEX NUMBERS
      DO i = 1 , M
         DO j = 1 , M
            a(i,j) = S(i,j)
         ENDDO
      ENDDO
C     PRINT 224
C 224 FORMAT(' ======    A OUTPUT.  =======')
C     PRINT 225,((A(I,J),J=1,M),I=1,M)
C 225 FORMAT(1X,4(1X,E13.6,1X,E13.6))
C
C     PRINT 226
C 226 FORMAT(/' ========  B OUTPUT. =======')
C     PRINT 23,(B(I23),I23=1,MM1)
C
      DO j = 1 , M
         l = M*(j-1)
         DO i = 1 , M
            k = 432 + l + i
            b(k) = a(i,j)
C
C     PRINT 240,J,L,I,K,B(K)
C 240 FORMAT(/,2X,'J=',I3,2X,'L=',I3,2X,'I=',I3,2X,'K=',I3,
C    *2X,'B(K)=4/4B',2E13.6)
            IF ( i.EQ.j ) b(k) = b(k) + w
C
C     PRINT 242,I,L,K,B(K)
C 242 FORMAT(/,2X,'I=',I3,2X,'L=',I3,2X,'K=',I3,
C    *2X,'B(K) ALONG THE DIAGONAL. +1=',2E13.6)
C
            k = 288 + l + i
            b(k) = -a(i,j)
            IF ( i.EQ.j ) b(k) = b(k) + w
C
C     PRINT 243,I,L,K,B(K)
C 243 FORMAT(/,2X,'I=',I3,2X,'L=',I3,2X,'K=',I3,
C    *2X,'3/4 B(K) MEASURES THE VALUE=',2E13.6)
C
         ENDDO
      ENDDO
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
      CALL ysminv(y,M,d,n1,n2)
C ********************************
C
C     PRINT 221
C 221 FORMAT(4X,'   =======  Y AFTER THE INCREMENT OF YSMINV =====')
C     PRINT 23,(Y(I221),I221=1,MM4)
C
C
      IF ( dsqrt(dreal(d)**2+dimag(d)**2).LT.1.D-30 ) THEN

         PRINT 445
 445     FORMAT (2X,
     &           '    ABORTION OF THE OPERATION IN THE TEST PROGRAM   ')
      ELSE
         DO i = 1 , M
            DO j = 1 , M
               d = z
               DO k = 1 , M
                  ik = 288 + (k-1)*M + i
                  kj = 432 + (j-1)*M + k
                  d = d + b(ik)*b(kj)
               ENDDO
               ij = 144 + (j-1)*M + i
               b(ij) = d
            ENDDO
         ENDDO
C
C
C     PRINT 232
C 232 FORMAT(2X,' =====  B2  ======')
C     PRINT 23,(B(I23),I23=1,MM1)
C
         mm = M*M
         DO i = 1 , mm
            b(288+i) = b(144+i)
         ENDDO
C     PRINT 23,(B(I23),I23=1,MM)
         DO j = 1 , M
            DO i = 1 , M
               a(i,j) = b(288+(j-1)*M+i)
            ENDDO
C
C     PRINT 270
C 270 FORMAT(2X,'==========   A RESULT.   ========')
C     PRINT 225,((A(I,J),J=1,M),I=1,M)
         ENDDO
      ENDIF
C
C     OVERWRITE MATRIX A WITH MATRIX YY
      DO i = 1 , M
         DO j = 1 , M
            Yy(i,j) = a(i,j)*y0
         ENDDO
      ENDDO

C       WRITE(6, 9875)
C 9875 FORMAT(2X, ' #######  YY B TEST''E  ######  ')
C       DO 465 II=1,M
C       DO 465 JJJ=1,M
C       WRITE(6, 9876) II,JJJ,YY(II,JJJ)
C 9876 FORMAT(2X,' YY(',I3,',',I3,')=',2(2X,E13.6 ))
C  465 CONTINUE
C

C     DEBUG SUBTRACE
      END

      SUBROUTINE ysminv(A,N,D,L,M)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION coabs
      INTEGER i , ij , ik , iz , j , ji , jk , jp , jq , jr , k , ki ,
     &        kj , kk , L , M , N , nk
      DIMENSION A(1) , L(1) , M(1)
      DOUBLE COMPLEX D , A , biga , hold
      coabs(D) = dsqrt(dreal(D)**2+dimag(D)**2)
      D = (1.D0,0.D0)
      nk = -N
      DO k = 1 , N
         nk = nk + N
         L(k) = k
         M(k) = k
         kk = nk + k
         biga = A(kk)
         DO j = k , N
            iz = N*(j-1)
            DO i = k , N
               ij = iz + i
               IF ( coabs(biga).LT.coabs(A(ij)) ) THEN
                  biga = A(ij)
                  L(k) = i
                  M(k) = j
               ENDIF
            ENDDO
         ENDDO
C  INTERCHANGE ROWS
         j = L(k)
         IF ( j.GT.k ) THEN
            ki = k - N
            DO i = 1 , N
               ki = ki + N
               hold = -A(ki)
               ji = ki - k + j
               A(ki) = A(ji)
               A(ji) = hold
            ENDDO
         ENDIF
C  INTERCHANGE COLUMNS
         i = M(k)
         IF ( i.GT.k ) THEN
            jp = N*(i-1)
            DO j = 1 , N
               jk = nk + j
               ji = jp + j
               hold = -A(jk)
               A(jk) = A(ji)
               A(ji) = hold
            ENDDO
         ENDIF
         IF ( coabs(biga).LE.1.D-20 ) THEN
            D = (0.D0,0.D0)
            RETURN
         ELSE
            DO i = 1 , N
               IF ( i.NE.k ) THEN
                  ik = nk + i
                  A(ik) = A(ik)/(-biga)
               ENDIF
            ENDDO
C  REDUCE MATRIX
            DO i = 1 , N
               ik = nk + i
               hold = A(ik)
               ij = i - N
               DO j = 1 , N
                  ij = ij + N
                  IF ( i.NE.k ) THEN
                     IF ( j.NE.k ) THEN
                        kj = ij - i + k
                        A(ij) = hold*A(kj) + A(ij)
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            kj = k - N
            DO j = 1 , N
               kj = kj + N
               IF ( j.NE.k ) A(kj) = A(kj)/biga
            ENDDO
            D = D*biga
            A(kk) = 1.D0/biga
         ENDIF
      ENDDO
C  FINAL ROW AND COLUMN INTERCHANGE
      k = N
      DO
         k = k - 1
         IF ( k.LE.0 ) RETURN
         i = L(k)
         IF ( i.GT.k ) THEN
            jq = N*(k-1)
            jr = N*(i-1)
            DO j = 1 , N
               jk = jq + j
               hold = A(jk)
               ji = jr + j

               A(jk) = -A(ji)
               A(ji) = hold
            ENDDO
         ENDIF
         j = M(k)
         IF ( j.GT.k ) THEN
            ki = k - N
            DO i = 1 , N
               ki = ki + N

               hold = A(ki)
               ji = ki - k + j

               A(ki) = -A(ji)
               A(ji) = hold
            ENDDO
         ENDIF
      ENDDO
C     DEBUG SUBTRACE,INIT
      END
