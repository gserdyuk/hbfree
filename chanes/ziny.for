!*==ZINY.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE ziny(Y,Vectj,Isize_maxnode)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , j
      INTEGER Isize_maxnode
      DOUBLE COMPLEX Vectj(Isize_maxnode)
      DOUBLE COMPLEX Y(Isize_maxnode,Isize_maxnode)
      DOUBLE COMPLEX zero/(0.0D0,0.0D0)/

C     ZERO OUT MATRIX Y AND VECTOR VEKTJ
C      print *, 'ZINY: ',ISIZE_MAXNODE
      DO i = 1 , Isize_maxnode
         Vectj(i) = zero
      ENDDO

      DO i = 1 , Isize_maxnode
         DO j = 1 , Isize_maxnode
            Y(j,i) = zero
         ENDDO
      ENDDO

      END
