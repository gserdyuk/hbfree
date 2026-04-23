!*==MACHEP.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE machep(Epsim)
C
C   SUBROUTINE FOR COMPUTING MACHINE EPSILON.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Epsim
C      print 5
C    5 format(2x,'ROUTINE MACHEP')

      Epsim = 1.D0
      DO
         Epsim = Epsim/2.D0

C      print 15, epsim

         IF ( (Epsim+1.0).EQ.1. ) RETURN
      ENDDO

C      EPSIM=2.*EPSIM
C      EPSIM=.1084202172D-18

C      print 15, epsim
C   15 format(2x,'EPSIM=',E16.10)

C     DEBUG INIT(EPSIM),SUBTRACE
      END
