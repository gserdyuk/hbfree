!*==REGUL.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE regul(Nk,Dfdx,F,N)
C   Subroutine for regulation of the Jacobian matrix.
C$LARGE:DFDX,F
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , j , N , Nk
      DOUBLE PRECISION pi
      DOUBLE PRECISION Dfdx , F
      DIMENSION Dfdx(Nk,Nk) , F(Nk)
      DO i = 1 , N
         pi = 0.D0
         DO j = 1 , N
            pi = pi + dabs(Dfdx(i,j))
         ENDDO
         IF ( pi.NE.0.D0 ) THEN
            DO j = 1 , N
               Dfdx(i,j) = Dfdx(i,j)/pi
            ENDDO
            F(i) = F(i)/pi
         ENDIF
         IF ( F(i).EQ.0.0D0 ) THEN
            Dfdx(i,i) = Dfdx(i,i)*1.01D0
            IF ( pi.EQ.0.D0 ) Dfdx(i,i) = 1.0D0
         ENDIF
      ENDDO

C      n1=n+1
C      do 20 i=1,n
C      print 10, F(I),I,(dfdx(i,j),j=1,n1)
C   10 format(2x,'REGUL : F=',E13.6,' DFDX(',I3,', J)=',(/3X,6(E13.6)))
C   20 continue

      END
