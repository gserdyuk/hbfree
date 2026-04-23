!*==GRADIE.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE gradie(Ntot,N,Dj,F,Sf,Gr)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , j , N , Ntot
      DOUBLE PRECISION F(1) , Sf(1)
      DOUBLE PRECISION Gr(1)
c$LARGE: DJ
      DOUBLE PRECISION Dj(Ntot,1)


      DO i = 1 , N
         Gr(i) = 0.D0
         DO j = 1 , N
            Gr(i) = Gr(i) + Dj(j,i)*F(j)*Sf((j+1)/2)**2
         ENDDO
      ENDDO
C  THE SAME SCALE IS APPLIED TO TWO ADJACENT ELEMENTS,
C  SINCE RE AND IM ARE PARTS OF THE SAME NUMBER.

C     DEBUG SUBTRACE,INIT(GR)
      END
