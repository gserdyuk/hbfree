!*==TOPO.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE topo(S,Nd)
C
C     SUBROUTINE FOR SHIFTING N NODES (MOVE FROM THE POSITION S TO
C     THE VALUE ND OR ND+1)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , ie , ii , iii , iti , ityp , nat , Nd , nfir , nk ,
     &        nlen , nn
      INCLUDE 'circuit.i'
      INTEGER S , end
      IF ( Nd.EQ.0 ) RETURN
      S = S
      Nd = Nd
      end = Nmpnt*20
C     LOOP OVER ELEMENT TYPES
      DO i = 1 , end , 20
         nfir = Mpoint(i+9)
         nlen = Mpoint(i+6)
         ityp = Mpoint(i+5)
         ie = Mpoint(i+4)
C     LOOP OVER ELEMENT TYPES OF INTERNAL TYPE
         DO ii = 1 , ie
            nat = nfir + (ii-1)*nlen
            IF ( ityp.EQ.2 ) THEN
               nn = 4
               nk = Nodeel(nat+1) + 3
            ELSEIF ( ityp.EQ.3 ) THEN
               nn = 1
               nk = Mpoint(i+7)
            ELSE
C     INITIALIZATION OF GRID CHANGES INDICES OF ELEMENTS VECTOR N,
C     AND STORAGE OF NODE VALUES (FROM VARIABLE TYPE)
               nn = 1
               nk = 2
            ENDIF
            DO iii = nn , nk
C     CALCULATION OF NODE VALUE
               iti = Nodeel(nat+iii)
               IF ( iti.NE.0 ) THEN
C     SHIFTING THE NODE VALUE
                  IF ( S.LT.0 .AND. iti.GT.Nd ) THEN
                     iti = iti - Nd
                  ELSEIF ( S.LT.0 .AND. iti.LE.Nd ) THEN
                     iti = iti - Nd - 1
                  ELSEIF ( S.GT.0 .AND. iti.GT.0 ) THEN
                     iti = iti + Nd
                  ELSE
                     IF ( S.LE.0 .OR. iti.GE.0 ) GOTO 60
                     iti = iti + Nd + 1
                  ENDIF
C
                  Nodeel(nat+iii) = iti
               ENDIF
 60         ENDDO
         ENDDO
      ENDDO
C     DEBUG SUBTRACE,INIT(ITI,S,ND)
      END
