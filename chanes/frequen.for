!*==FREQUEN.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE frequen(Nf3,Kpr)
C ***
C *** FREQUENCY PROCESSING PROGRAM
C ***
C       LAST EDITING DATE  -  28.04.92
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION F , f1 , f1old , f2 , f2old , Pniff
      INTEGER ic , Iff , iffr , iffro , iffrp , iin , ik , inn1 , inn2 ,
     &        ip , irr , Kiff , Knniff , knpri , Kpr , nf1 , nf2 , Nf3 ,
     &        ng , Nniff
      INTEGER npri
      INCLUDE 'charint.i'
      INCLUDE 'circuit.i'

      COMMON /fre   / F(2)
      CHARACTER*4 is , Fne1 , Fne2
C
      COMMON /bliff / Iff(7,4) , Kiff , Nniff(4,8) , Knniff , Pniff(8) ,
     &                Fne1 , Fne2
      INTEGER iffron
C
      DATA ip/6/ , knpri/4/
      DATA iin/10/ , is/'    '/

C
C *** initialization of variables
      iffron = 0


C *** INITIAL STATE ASSIGNMENT
      f1 = 0.D0
      f2 = 0.D0
      nf1 = 0
      nf2 = 0
      Nf3 = 0
      Fne1 = is
      Fne2 = is
      inn1 = 0
      inn2 = 0

      DO ik = 1 , Knniff
         Nniff(4,ik) = 0
      ENDDO

C  STORE IN THE ARRAY PNIFF THE FREQUENCY VALUES RECORDED IN THE
C  ARRAY PARAM

      DO ic = 1 , Knniff
         Pniff(ic) = Param(Nniff(3,ic))
         IF ( Kpr.GE.2 ) WRITE (6,3004) ic , Pniff(ic) , Nniff(3,ic)
 3004    FORMAT (2X,'FREQUEN : PNIFF(',I2,')=',E13.6,' (PARAM(',I4,'))')
      ENDDO

C  NPRI - PRIORITY NUMBER OF THE FREQUENCY CORRESPONDING TO IFF(3,N) AND
C         NNIFF(2,K)
      npri = 1

C  IFFRP - CURRENT FREQUENCY NUMBER (POSITIVE)
C  IFFRO - CURRENT FREQUENCY NUMBER (ZERO)
      iffrp = 0
      iffro = 0
      DO
         DO ic = 1 , Knniff
C  DETERMINE THE ELEMENT WITH THE HIGHEST PRIORITY AND,
C  DEPENDING ON THE FREQUENCY VALUE, STORE THE CURRENT
C  FREQUENCY VALUE FOR THAT ELEMENT
            IF ( Nniff(2,ic).EQ.npri ) THEN
               IF ( Pniff(ic).EQ.0.D0 .AND. Nniff(4,ic).EQ.0 )
     &              iffro = iffro - 1
               IF ( Pniff(ic).EQ.0.D0 .AND. Nniff(4,ic).EQ.0 )
     &              Nniff(4,ic) = iffro
               IF ( Pniff(ic).NE.0.D0 .AND. Nniff(4,ic).EQ.0 )
     &              iffrp = iffrp + 1
               IF ( Pniff(ic).NE.0.D0 .AND. Nniff(4,ic).EQ.0 )
     &              Nniff(4,ic) = iffrp
               IF ( Kpr.GE.3 ) WRITE (6,3017) npri , ic , Nniff(4,ic)
 3017          FORMAT (2X,'NPRI=',I2,', NNIFF(4,',I2,')=',I4)

C  SCAN THE REMAINING ELEMENTS. IF AN ELEMENT WITH THE SAME
C  FREQUENCY VALUE IS FOUND, ASSIGN THE SAME CURRENT FREQUENCY
C  VALUE TO THAT ELEMENT.
               DO irr = 1 , Knniff
                  IF ( ic.NE.irr ) THEN
                     IF ( Nniff(4,irr).EQ.0 ) THEN
                        IF ( Pniff(ic).EQ.Pniff(irr) ) Nniff(4,irr)
     &                       = Nniff(4,ic)
                        IF ( Kpr.GE.3 ) WRITE (6,3019) irr , ic ,
     &                       Nniff(4,irr)
 3019                   FORMAT (2X,'IRR=',I2,', IC=',I2,
     &                          ', NNIFF(4,IRR)=',I3)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         npri = npri + 1
         IF ( npri.GT.knpri ) THEN

            DO ic = 1 , Knniff
               IF ( Kpr.GE.3 ) WRITE (6,3040) ic , (Nniff(ng,ic),ng=1,4)
 3040          FORMAT (2X,'NNIFF( ,',I3,')=',A4,',',I3,',',I3,',',I3)
            ENDDO

C  SELECTION OF THE 1ST FREQUENCY
            npri = 1
            DO
               DO ic = 1 , Knniff
                  IF ( Nniff(2,ic).EQ.npri ) THEN
                     IF ( Nniff(4,ic).EQ.1 ) GOTO 3070
                  ENDIF
               ENDDO
               npri = npri + 1
               IF ( npri.GT.knpri ) THEN

C  SELECTION OF THE 2ND FREQUENCY
                  npri = 1
                  GOTO 3077
               ENDIF
            ENDDO
 3070       f1 = Pniff(ic)
            Fne1 = i2c(Nniff(1,ic))
            inn1 = Nniff(3,ic)
            nf1 = 1

            IF ( Kpr.GE.3 ) WRITE (ip,3073) f1 , nf1 , inn1 , Fne1
 3073       FORMAT (2X,'F1=',E12.5,' NF1=',I3,' INN1=',I4,' FNE1=',A4)
            npri = 1
            GOTO 3077
         ENDIF
      ENDDO
 3077 DO
         DO ic = 1 , Knniff
            IF ( Nniff(2,ic).EQ.npri ) THEN
               IF ( Nniff(4,ic).EQ.2 ) GOTO 3090
            ENDIF
         ENDDO
         npri = npri + 1
         IF ( npri.GT.knpri ) THEN

C  PROCESSING OF ZERO FREQUENCY
            npri = 1
            GOTO 3107
         ENDIF
      ENDDO
 3090 f2 = Pniff(ic)
      Fne2 = i2c(Nniff(1,ic))
      inn2 = Nniff(3,ic)
      nf2 = 1

      IF ( Kpr.GE.3 ) WRITE (ip,3100) f2 , nf2 , inn2 , Fne2

      IF ( nf1.EQ.1 .AND. nf2.EQ.1 .AND. f1.NE.0 .AND. f2.NE.0.D0 )
     &     GOTO 3123
      npri = 1
 3107 DO
         iffro = -iffro
         DO iffr = 1 , iffro
            DO ic = 1 , Knniff
               IF ( Nniff(2,ic).EQ.npri ) THEN
                  IF ( Nniff(4,ic).EQ.-iffr ) GOTO 3120
               ENDIF
            ENDDO
         ENDDO
         npri = npri + 1
         IF ( npri.GT.knpri ) GOTO 3140
      ENDDO
 3120 f2 = Pniff(ic)
      Fne2 = i2c(Nniff(1,ic))
      inn2 = Nniff(3,ic)
      nf2 = 1
      Nniff(4,ic) = 2
      iffron = -iffr
      IF ( Kpr.GE.3 ) WRITE (ip,3100) f2 , nf2 , inn2 , Fne2

C  CHECK IF THERE ARE ANY MORE ZERO FREQUENCIES?
 3123 DO ic = 1 , Knniff
         IF ( Nniff(4,ic).LT.0 .AND. iffrp.LT.2 ) Nniff(4,ic) = 2
         IF ( Nniff(4,ic).LT.0 .AND. iffrp.GE.2 ) Nniff(4,ic) = iffron
         IF ( Kpr.GE.3 ) WRITE (6,3125) ic , Nniff(4,ic)
 3125    FORMAT (2X,'NNIFF(4,',I2,')=',I3)
      ENDDO

C  CHECK IF THERE ARE ANY FREQUENCIES OTHER THAN F1 AND F2?
 3140 DO ic = 1 , Knniff
         IF ( Nniff(4,ic).GT.2 ) Nf3 = 1
      ENDDO
C
      DO ic = 1 , Knniff
         IF ( Kpr.GE.3 ) WRITE (6,3160) ic , (Nniff(ng,ic),ng=1,4)
 3160    FORMAT (2X,'FREQUEN (OUT) : NNIFF( ,',I3,')=',A4,',',I3,',',I3,
     &           ',',I3)
      ENDDO
C
      DO ic = 1 , Knniff
         IF ( Kpr.GE.3 ) WRITE (6,3180) ic , Pniff(ic)
 3180    FORMAT (2X,'FREQUEN (OUT) : PNIFF(',I2,')=',E13.6)
      ENDDO
C
      f1old = f1
      f2old = f2

      F(1) = f1
      F(2) = f2
C
C

      IF ( Kpr.GE.2 ) WRITE (ip,4000) f1 , Fne1
 4000 FORMAT (2X,'FREQUEN: FREQUENCY VARIATION'/2X,'FREQUENCY 1 =',
     &        E12.5,' ( ELEMENT ',A4,' )')


      IF ( Kpr.GE.2 ) WRITE (ip,4010) f2 , Fne2
 4010 FORMAT (2X,'FREQUENCY 2 =',E12.5,' ( ELEMENT ',A4,' )')

      IF ( Kpr.GE.2 .AND. Nf3.EQ.1 ) WRITE (ip,4020)
 4020 FORMAT (2X,'ATTENTION: MORE THAN 2 FREQ. DEFINED!')
 3100 FORMAT (2X,'F2=',E12.5,' NF2=',I3,' INN2=',I4,' FNE2=',A4)

C
      END
