!*==INKOOR.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE inkoor(Maxsyn,Mnmax,F1,F2,*,*)
C
C  *****  ORDERING OF ARRAYS MN AND MN1,
C  *****  FILLING THE COORDINATE ARRAYS,
C  *****  OBTAINING THE EXPANDED FREQUENCY GRID,
C  *****  FILLING THE ARRAYS WR AND WS.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION alog2 , F1 , F2 , omega , omega1 , omega2 , pi
      INTEGER i , i1 , i2 , ipq1 , ipq2 , ir , ivald , ivald1 , j , j1 ,
     &        j2 , k , kk , Kn , Kn1 , Knc , kncmax , kncmin , Knr ,
     &        Knr1
      INTEGER krdchk , l , Maxsyn , Mephf , mnm , Mnmax , mp , ndim ,
     &        ndim1
      INCLUDE 'funcsize.i'
      COMMON /blw1  / W , W1
      DOUBLE PRECISION W(20) , W1(200)
      COMMON /blw2  / Wr , Ws
      INTEGER Wr(20,20) , Ws(20,20)
      COMMON /blk1  / Kr , Kr1 , Kc , Kc1 , Nnr , Nnr1 , Mn , Mn1
      INTEGER Kr(20) , Kc(10) , Nnr(10)
      INTEGER Kr1(200) , Kc1(20) , Nnr1(20)
      INTEGER Mn(2,20) , Mn1(2,200)
      COMMON /blk2  / Knc , Knr , Kn , Knr1 , Kn1
      COMMON /mep   / Mephf , Flgmnw
      INTEGER Flgmnw
      INTEGER b1dim , b2dim
      DATA ir/10/
C
C take sizes from 'funcsize.i'
      b1dim = B1_SIZE
      b2dim = B2_SIZE
      ndim = MAXKN
      ndim1 = MAXKN1

      pi = 4*datan(1.0D0)
      omega1 = 2.0D0*pi*F1
      omega2 = 2.0D0*pi*F2
      DO i = 1 , Kn
         IF ( Mn(1,i).LE.0 ) THEN
            IF ( Mn(1,i).NE.0 .OR. Mn(2,i).LT.0 ) THEN
               Mn(1,i) = -Mn(1,i)
               Mn(2,i) = -Mn(2,i)
            ENDIF
         ENDIF
      ENDDO
      Flgmnw = 0
      Mephf = 2
      IF ( Kn.NE.1 ) THEN
         CALL sort(Mn,Kn)




C  *****  CHANGE OF FREQUENCIES F1 AND F2, IF F2 = 0

         DO i = 1 , Kn
            IF ( Mn(2,i).NE.0 ) GOTO 1030
         ENDDO
         DO i = 1 , Kn
            mnm = Mn(1,i)
            Mn(1,i) = Mn(2,i)
            Mn(2,i) = mnm
         ENDDO
         omega = omega1
         omega1 = omega2
         omega2 = omega
         Mephf = 1
         Flgmnw = 1
C  *****  EXCLUSION OF IDENTICAL ROWS OF THE MATRIX
C  *****  AND COMPRESSION OF THE LIST MN AFTER EXCLUSION
C
 1030    i = 2
         DO WHILE ( i.LE.Kn )
            IF ( Mn(1,i).EQ.Mn(1,i-1) .AND. Mn(2,i).EQ.Mn(2,i-1) ) THEN
               Kn = Kn - 1
               IF ( i.GT.Kn ) GOTO 1065
               DO j = i , Kn
                  Mn(1,j) = Mn(1,j+1)
                  Mn(2,j) = Mn(2,j+1)
               ENDDO
            ELSE
               i = i + 1
            ENDIF
         ENDDO
      ENDIF
C  *****  ADDITION OF THE CONSTANT COMPONENT,
C  *****  IF IT IS ABSENT IN THE LIST MN
 1065 IF ( Mn(1,1).NE.0 .OR. Mn(2,1).NE.0 ) THEN
         DO i = 1 , Kn
            Mn(1,Kn-i+2) = Mn(1,Kn-i+1)
            Mn(2,Kn-i+2) = Mn(2,Kn-i+1)
         ENDDO
         Mn(1,1) = 0
         Mn(2,1) = 0
         Kn = Kn + 1
      ENDIF
      IF ( Kn.EQ.1 ) Mephf = 0
      DO i = 1 , Kn
         W(i) = Mn(1,i)*omega1 + Mn(2,i)*omega2
      ENDDO
      Kn1 = 0
      DO i = 1 , Kn
         i1 = Mn(1,i)
         i2 = Mn(2,i)
         DO j = 1 , Kn
            j1 = Mn(1,j)
            j2 = Mn(2,j)
            mp = 1
            DO l = 1 , 2
               ipq1 = i1 + mp*j1
               ipq2 = i2 + mp*j2
               IF ( ipq1.LT.0 ) THEN
                  ipq1 = -ipq1
                  ipq2 = -ipq2
               ELSEIF ( ipq1.EQ.0 ) THEN
                  IF ( ipq2.LT.0 ) ipq2 = -ipq2
               ENDIF
               IF ( Kn1.NE.0 ) THEN
                  DO k = 1 , Kn1
                     IF ( ipq1.EQ.Mn1(1,k) .AND. ipq2.EQ.Mn1(2,k) )
     &                    GOTO 60
                  ENDDO
               ENDIF
               Kn1 = Kn1 + 1
               IF ( Kn1.GT.ndim1 ) GOTO 200
               Mn1(1,Kn1) = ipq1
               Mn1(2,Kn1) = ipq2
 60            mp = -1
            ENDDO
         ENDDO
      ENDDO
      CALL sort(Mn1,Kn1)
      DO i = 1 , Kn1
         W1(i) = Mn1(1,i)*omega1 + Mn1(2,i)*omega2
      ENDDO
C  *****  CALCULATION OF THE NUMBER OF NONZERO ROWS
C  *****  ON THE EXPANDED FREQUENCY GRID - KNR1
      Knr1 = 1
      IF ( Kn.GE.2 ) THEN
         DO i = 2 , Kn1
            IF ( Mn1(1,i).NE.Mn1(1,i-1) ) Knr1 = Knr1 + 1
         ENDDO
      ENDIF
C  *****  SELECTION OF LIMITS MIN AND MAX FOR KNC
      Mnmax = 0
      alog2 = dlog10(2.D0)
      DO i = 1 , Kn
         Mnmax = max0(Mnmax,Mn(1,i),iabs(Mn(2,i)))
      ENDDO
      kncmin = max0(Mnmax*4,16)
      kncmin = 2**(idint(dlog10(kncmin-0.5D0)/alog2)+1)
      kncmax = min0(b1dim/Maxsyn,b2dim/(Maxsyn*Knr1),128)
      kncmax = 2**(idint(dlog10(kncmax+0.5D0)/alog2))
      IF ( kncmin.GT.kncmax ) THEN
C----------------------------------------------------------------------
C  *****  THE SPREAD OF THE TRANSFORMATION OF THE FOURIER
C  *****  EXCEEDS THE UPPER LIMIT VALUE
C
         Mnmax = kncmax/4

         RETURN 2
      ELSE
         IF ( Knc.GT.kncmin .AND. Knc.LT.kncmax )
     &        Knc = 2**(idint(dlog10(Knc-0.5D0)/alog2)+1)
         IF ( Knc.LE.kncmin ) Knc = kncmin
         IF ( Knc.GE.kncmax ) Knc = kncmax
C  *****  FILLING THE COORDINATE ARRAYS
         CALL koord(Mn,Kr,Kc,Nnr,Knr,Kn,Knc,ir,ndim)

C++++++++++++++++++++++++  CHANGES FROM 12.05.91  SERDYUK G.V.
C                          (TOTAL 3 CHANGES)
         ivald = krdchk(Mn,Kr,Kc,Nnr,Knr,Kn,Knc,ir,ndim)
         IF ( ivald.NE.0 ) THEN
            WRITE (6,5010) ivald
            PRINT 5010 , ivald
         ENDIF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C  *****  FILLING THE COORDINATE ARRAYS
C  *****  FOR THE EXPANDED FREQUENCY GRID
         CALL koord(Mn1,Kr1,Kc1,Nnr1,Knr1,Kn1,Knc,ndim,ndim1)
C++++++++++++++++++++++++  CHANGES FROM 12.05.91  SERDYUK G.V.
         ivald1 = krdchk(Mn1,Kr1,Kc1,Nnr1,Knr1,Kn1,Knc,ndim,ndim1)
         IF ( ivald1.NE.0 ) THEN
            WRITE (6,5020) ivald1
            PRINT 5020 , ivald1
         ENDIF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         IF ( Nnr(Knr).EQ.1 .AND. Mephf.EQ.2 ) Mephf = 1
C  *****  FILLING THE ARRAYS WR AND WS
         DO i = 1 , Kn
            i1 = Mn(1,i)
            i2 = Mn(2,i)
            DO j = 1 , Kn
               j1 = Mn(1,j)
               j2 = Mn(2,j)
               mp = 1
               DO l = 1 , 2
                  ipq1 = i1 + mp*j1
                  ipq2 = i2 + mp*j2
                  kk = 1
                  IF ( ipq1.LT.0 ) THEN
                     ipq1 = -ipq1
                  ELSEIF ( ipq1.EQ.0 ) THEN
                     IF ( ipq2.GE.0 ) GOTO 90
                  ELSE
                     GOTO 90
                  ENDIF
                  ipq2 = -ipq2
                  kk = -1
 90               DO k = 1 , Kn1
                     IF ( ipq1.EQ.Mn1(1,k) .AND. ipq2.EQ.Mn1(2,k) )
     &                    GOTO 100
                  ENDDO
 100              k = k*kk
                  IF ( l.EQ.1 ) Ws(i,j) = k
                  IF ( l.EQ.2 ) Wr(i,j) = k
                  mp = -1
               ENDDO
            ENDDO
         ENDDO
         RETURN
      ENDIF
C     DEBUG INIT(W,OMEGA1,OMEGA2,MN)
C----------------------------------------------------------------------
C  *****  THE SPREAD OF MN1 EXCEEDS THE MAXIMUM
 200  RETURN 1
C++++++++++++++++++++++++  CHANGES FROM 12.05.91  SERDYUK G.V.
 5010 FORMAT (2x,' INKOOR: ERROR OF REPRESENTATION OF MAIN GRID.',
     &        ' CODE =',i3)
 5020 FORMAT (2x,' INKOOR: ERROR OF REPRESENTATION OF AUX  GRID.',
     &        ' CODE =',i3)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      END
