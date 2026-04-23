!*==RESTOR.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE restor(U,Y,Vectj,Isize_maxnode,S)

C     SUBROUTINE FOR CALCULATING THE CURRENTS OF NODES AND STRESSES IN THE EXTERNAL NODES,
C     AND ALSO THE CALCULATION OF THE STRESS CURRENTS. IN ACCORDANCE WITH THE USER NUMBER.

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER j
      COMMON /kolnal/ Kol(4) , Nal(4)
      COMMON /blk2  / Knc , Knr , Kn , Knr1 , Kn1
      LOGICAL Nal
      INTEGER Kol , k12 , Knc , Knr , Kn , Knr1 , Kn1 , irc
      INTEGER Isize_maxnode
      DOUBLE COMPLEX S(Isize_maxnode,20)
      DOUBLE PRECISION U(1)

      DOUBLE COMPLEX Y(Isize_maxnode,Isize_maxnode)
      DOUBLE COMPLEX Vectj(Isize_maxnode)
C      PRINT 3
C    3 FORMAT(2X,'ENTRY INTO THE SUBROUTINE RESTOR')

      k12 = Kol(1) + Kol(2)
      CALL topo(1,k12)
      DO j = 1 , Kn
         irc = j
C      print *, 'RESTOR: IRC: ', IRC
         CALL stback(U,irc,Y,Vectj,Isize_maxnode,S)
         CALL decode(S(1,irc),irc,Isize_maxnode)
      ENDDO

      END
