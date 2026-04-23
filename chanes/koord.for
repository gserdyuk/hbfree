!*==KOORD.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE koord(Mn,Kr,Kc,Nnr,Knr,Kn,Knc,Ir1,Ir2)
C
C  *****  FILLING THE COORDINATE ARRAYS
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , Ir1 , Ir2 , j , Kn , Knc , Knr , m , m1 , mnr , ms , n
      INTEGER Kr(Ir2) , Kc(Ir1) , Nnr(Ir1)
      INTEGER Mn(2,Ir2)
C
      DO i = 1 , Ir1
         Kc(i) = 0
         Nnr(i) = 0
      ENDDO
      DO i = 1 , Ir2
         Kr(i) = 0
      ENDDO
      Knr = 0
      mnr = -1
      DO i = 1 , Kn
         n = Mn(2,i)
         IF ( n.LT.0 ) n = n + Knc
         Kr(i) = n + 1
         m = Mn(1,i)
         IF ( m.NE.mnr ) THEN
            Knr = Knr + 1
            Nnr(Knr) = m + 1
            mnr = m
         ENDIF
         IF ( i.LE.1 ) THEN
            ms = 1
            m1 = m
            j = 1
C ++++++++++++++++++++++++++++ CHANGE FROM 12.05.91
C                             TOTAL 1 PIECE.
            IF ( i.NE.Kn ) GOTO 40
         ELSEIF ( m.NE.m1 .AND. i.LT.Kn ) THEN
         ELSEIF ( .NOT.m.NE.m1 .AND. i.LT.Kn ) THEN
            ms = ms + 1
            GOTO 40
         ELSEIF ( m.NE.m1 .AND. .NOT.i.LT.Kn ) THEN
            Kc(j) = ms
            Kc(j+1) = 1
            GOTO 40
         ELSE
            ms = ms + 1
         ENDIF
C +++++++++++++++++++++++++++++++++++++++++++++++
         Kc(j) = ms
         ms = 1
         m1 = m
         j = j + 1
 40   ENDDO
      END
