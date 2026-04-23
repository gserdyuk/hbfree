!*==NEINCK.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE neinck(N,Epsim,Idim,Scalef,Scaleu,U,Icode)

C     Performs validation and initial setup
C     of the parameters for the SOLVE program

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Epsim
      INTEGER i , Icode , Idim , jj , Limit , N , ndiv2
      COMMON /typval/ Typu , Typi
      DOUBLE PRECISION Typu , Typi
      DOUBLE PRECISION Scalef(1) , Scaleu(1)
      COMMON /newton/ Epssol , Epsdu , Epsmin , Maxdu , Limit
      DOUBLE PRECISION Epssol , Epsdu , Epsmin , Maxdu
      INTEGER limmax/1000/
      DOUBLE PRECISION U(1)
      DOUBLE PRECISION norstp , nor1

C  Validation of the value of N - the dimensionality of the system
      IF ( N.LT.1 ) THEN
         Icode = -1
         RETURN
      ELSEIF ( N.LE.Idim ) THEN
C Checking the scale settings
         IF ( Typu.LE.0.D0 ) Typu = 1.D0
         IF ( Typi.LE.0.D0 ) Typi = 1.D0
         ndiv2 = N/2
         DO i = 1 , ndiv2
            Scalef(i) = 1.D0/Typi
            Scaleu(i) = 1.D0/Typu
         ENDDO
C Setting the parameters of the method
C  EPSSOL
         IF ( Epssol.LE.0 ) Epssol = Epsim**(1.D0/3.D0)
C  EPSDU
         IF ( Epsdu.LE.0 ) Epsdu = Epsim**(2.D0/3.D0)
C  EPSMIN
         IF ( Epsmin.LE.0 ) Epsmin = Epsim**(2.D0/3.D0)
C  MAXDU
         IF ( Maxdu.LE.0.D0 ) THEN
            norstp = 0.D0
            nor1 = 0.D0
            DO jj = 1 , N
               norstp = norstp + (Scaleu((jj+1)/2)*dabs(U(jj)))**2
               nor1 = nor1 + (Scaleu((jj+1)/2))**2
            ENDDO
            norstp = dsqrt(norstp)
            nor1 = dsqrt(nor1)
            Maxdu = 1000.D0*dmax1(norstp,nor1)
         ENDIF
C  Limit on the number of iterations
         IF ( Limit.LE.0 .OR. Limit.GT.limmax ) Limit = limmax
      ELSE
         Icode = -2
         RETURN
      ENDIF
C     DEBUG SUBCHK,INIT
      END
