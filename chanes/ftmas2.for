!*==FTMAS2.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE ftmas2(Znn,Kr,Kc,Nnr,Knr,Knc,Kn,Nr,Nb,M1,B1,B2,Is,
     &                  Flgfft,*)
C
C-----------------------------------------------------------------------
C
C     THE SUBROUTINE IMPLEMENTS A TWO-DIMENSIONAL INVERSE FFT,
C     OUTPUT TO AN INTERMEDIATE NONLINEAR TRANSFORMATION,
C     AND A TWO-DIMENSIONAL DIRECT FFT OF TWO-DIMENSIONAL
C     SPARSE VECTORS.
C
C     IDENTIFIERS OF VARIABLES AND ARRAYS:
C
C          ZNN - PACKED ARRAY OF AMPLITUDES;
C          KNC - LENGTH OF THE ONE-DIMENSIONAL FFT;
C          NNR - ARRAY OF INDICES OF NONZERO ROWS;
C          KNR - NUMBER OF NONZERO ROWS;
C          KR  - ARRAY OF INDICES OF COLUMNS OF NONZERO
C                ELEMENTS IN EACH NONZERO ROW;
C          KC  - ARRAY CONTAINING THE NUMBER OF
C                NONZERO ELEMENTS IN EACH ROW;
C          B1  - WORKING ARRAY FOR STORING THE RESULTS
C                OF THE ONE-DIMENSIONAL FFT;
C          B2  - BUFFER FOR STORING TRANSFORMED ROWS;
C          NR  - NUMBER OF SIMULTANEOUSLY TRANSFORMED
C                AMPLITUDES IN ZNN DURING THE INVERSE FFT;
C          NB  - NUMBER OF SIMULTANEOUSLY TRANSFORMED
C                VECTORS OF TEMPORARY SAMPLES
C                DURING THE DIRECT FFT.
C
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , iferr , inv , Is , j , k , k1 , Kn , Knc , knc1 ,
     &        knc2 , Knr , l , lp , M1 , Mephf , n , n1 , n2 , Nb
      INTEGER Nr
      DOUBLE PRECISION s
      COMMON /mep   / Mephf , Flgmnw
      DOUBLE COMPLEX Znn(Kn,Nr) , B1(Knc,Nr) , B2(Knc,Knr,Nr) ,
     &               zero/(0.0D0,0.0D0)/
      INTEGER Flgmnw , Flgfft
      INTEGER Kr(1) , Kc(1) , Nnr(1)
C     INTEGER  INV(32),S(32),M(3)/0,0,0/
      DIMENSION inv(32) , s(32)
      INTEGER m(3)/0 , 0 , 0/
      SAVE inv , s
C  ARRAYS    INV(32), S(32), M(3) - WORKING ARRAYS FOR SUBROUTINE HARM
C     PRINT 976
C 976 FORMAT('    ZNN:')
C     PRINT 975,((ZNN(I,J),         J=1,NR), I=1,KN)
C 975 FORMAT(4(2X,E12.5))
C 977 FORMAT('    B1:')

      m(1) = M1
      IF ( Flgfft.EQ.0 ) THEN
         CALL harm(B1,m,inv,s,0,iferr)
         Flgfft = 1
      ENDIF
      knc1 = Knc/2
      knc2 = Knc + 2
C
C  FORMATION OF THE ROW ARRAY
C
      DO l = 1 , Nr
         n1 = 1
         n2 = 0
         DO k = 1 , Knr
            DO i = 1 , Knc
               B1(i,l) = zero
            ENDDO
            n2 = n2 + Kc(k)
            DO n = n1 , n2
               i = Kr(n)
               B1(i,l) = Znn(n,l)
            ENDDO
            n1 = n2 + 1
            IF ( Mephf.NE.0 ) THEN
               IF ( Nnr(k).EQ.1 ) THEN
                  DO i = 2 , knc1
                     k1 = knc2 - i
                     B1(k1,l) = dconjg(B1(i,l))
                  ENDDO
               ENDIF
C  FOURIER TRANSFORMATION (DIRECT) OF A ROW
               CALL harm(B1(1,l),m,inv,s,2,iferr)
               IF ( Mephf.NE.1 ) THEN
                  DO i = 1 , Knc
                     B2(i,k,l) = B1(i,l)
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDDO
C
C  PROCESSING OF COLUMNS FROM THE ROW ARRAY
C
      Is = Knc
      IF ( Mephf.NE.2 ) GOTO 1056
      Is = 1
 34   DO l = 1 , Nr
         DO n = 1 , Knc
            B1(n,l) = zero
         ENDDO
         DO k = 1 , Knr
            n = Nnr(k)
            B1(n,l) = B2(Is,k,l)
         ENDDO
         DO j = 2 , knc1
            k1 = knc2 - j
            B1(k1,l) = dconjg(B1(j,l))
         ENDDO
C  FOURIER TRANSFORMATION (DIRECT) OF A COLUMN
         CALL harm(B1(1,l),m,inv,s,2,iferr)
      ENDDO
 1056 RETURN 1
C-----------------------------------------------------------------------
      ENTRY ft2(Znn,Kr,Kc,Nnr,Knr,Knc,Kn,Nr,Nb,B1,B2,Is,Flgfft,*)
C
      IF ( Nb.GT.0 ) THEN
         IF ( Mephf.EQ.2 ) THEN
C  FOURIER TRANSFORMATION (INVERSE) OF A COLUMN
            DO lp = 1 , Nb
               CALL harm(B1(1,lp),m,inv,s,-2,iferr)
               DO k = 1 , Knr
                  n = Nnr(k)
                  B2(Is,k,lp) = B1(n,lp)
               ENDDO
            ENDDO
            Is = Is + 1
            IF ( Is.LE.Knc ) GOTO 34
         ENDIF
C
C  PROCESSING OF ROWS FROM THE COLUMN ARRAY
         DO lp = 1 , Nb
            n1 = 1
            n2 = 0
            DO k = 1 , Knr
               IF ( Mephf.NE.0 ) THEN
                  IF ( Mephf.NE.1 ) THEN
                     DO i = 1 , Knc
                        B1(i,lp) = B2(i,k,lp)
                     ENDDO
                  ENDIF
C  FOURIER TRANSFORMATION (INVERSE) OF ROWS
C     PRINT 977
C     PRINT 975,(B1(II,LP),II=1,KNC)
                  CALL harm(B1(1,lp),m,inv,s,-2,iferr)
               ENDIF
C     PRINT 977
C     PRINT 975,(B1(II,LP),II=1,KNC)
               n2 = n2 + Kc(k)
               DO n = n1 , n2
                  i = Kr(n)
                  Znn(n,lp) = B1(i,lp)
               ENDDO
               n1 = n2 + 1
            ENDDO
         ENDDO
      ELSE
         Is = Is + 1
         IF ( Is.LE.Knc ) GOTO 34
      ENDIF
C     PRINT 976
C     PRINT 975,((ZNN(I,J),J=1,NB),I=1,KN)
C     DEBUG SUBTRACE,INIT(N1,N2,IFERR,M,FLGFFT)
      END
