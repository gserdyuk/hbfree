!*==LINEQ1.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code


      SUBROUTINE lineq1(A,Ndim,N,Ndim2,Nprz,Ip,Ner,Det,Argd)
C$LARGE: A,IP
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION A , Argd , t
      INTEGER i , Ip , j , k , kb , km1 , kp1 , m , N , Ndim , Ndim2 ,
     &        Ner , nm1 , np1 , Nprz
      DIMENSION A(Ndim,Ndim2) , Ip(Ndim)
      DOUBLE COMPLEX Det
C
C      PRINT 100, N, NDIM
C  100 FORMAT(2X,'LINEQ1 : N=',I5,' NDIM=',I5)

      Ip(N) = 1
      DO k = 1 , N
         IF ( k.NE.N ) THEN
            kp1 = k + 1
            m = k
            DO i = kp1 , N
               IF ( dabs(A(i,k)).GT.dabs(A(m,k)) ) m = i
            ENDDO
            Ip(k) = m
            IF ( m.NE.k ) Ip(N) = -Ip(N)
            t = A(m,k)
            A(m,k) = A(k,k)
            A(k,k) = t
            IF ( t.NE.0.0D0 ) THEN
               DO i = kp1 , N
                  A(i,k) = -A(i,k)/t
               ENDDO
               DO j = kp1 , N
                  t = A(m,j)
                  A(m,j) = A(k,j)
                  A(k,j) = t
                  IF ( t.NE.0.0D0 ) THEN
                     DO i = kp1 , N
                        A(i,j) = A(i,j) + A(i,k)*t
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
         IF ( dabs(A(k,k)).EQ.0.D0 ) THEN
            Ip(N) = 0
            GOTO 20
         ENDIF
      ENDDO
 20   Ner = 0
      Det = (1.0D0,0.0D0)
      Argd = 1.D0
      IF ( Nprz.LT.1 ) RETURN
      IF ( Ip(N).EQ.0 ) THEN
         WRITE (6,17) k
 17      FORMAT (///5X,'LINEQ1....'/15X,I6,
     &           '-TH COLUMN CONTAINS ONLY ZE','ROES'/)
         Ner = 1
         RETURN
      ELSE
         np1 = N + 1
         IF ( N.NE.1 ) THEN
            nm1 = N - 1
            DO k = 1 , nm1
               kp1 = k + 1
               m = Ip(k)
               t = A(m,np1)
               A(m,np1) = A(k,np1)
               A(k,np1) = t
               DO i = kp1 , N
                  A(i,np1) = A(i,np1) + A(i,k)*t
               ENDDO
            ENDDO
            DO kb = 1 , nm1
               km1 = N - kb
               k = km1 + 1
               A(k,np1) = A(k,np1)/A(k,k)
               t = -A(k,np1)
               DO i = 1 , km1
                  A(i,np1) = A(i,np1) + A(i,k)*t
               ENDDO
            ENDDO
         ENDIF
      ENDIF
      A(1,np1) = A(1,np1)/A(1,1)
      RETURN
      END
