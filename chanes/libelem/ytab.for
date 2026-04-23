!*==YTAB.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE ytab(Om,P1,L1,P2,L2,P3,L3,N)
C *** SUBROUTINE FOR THE FORMATION OF THE Y-MATRIX FROM GIVEN
C     TASK PARAMETERS
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION df , df1 , f0 , fi , Om , pi
      INTEGER i , i10 , i40 , i50 , ib , ipi , ipr , j10 , kfr , L1 ,
     &        L2 , L3 , m , N
      COMMON /subc  / Y(15,15) , J(15)
      DOUBLE COMPLEX Y , J
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
C  DETERMINATION OF THE NUMBERS OF THE GIVEN DISCRETE POINTS
      m = L3/(2*N**2+1)
      DO i10 = 1 , N
         DO j10 = 1 , N
            Y(i10,j10) = dcmplx(0.D0,0.D0)
         ENDDO
      ENDDO
C  WE ASSUME THAT FI IS CLOSER TO F0
C  THAN ANY FJ FOR J = 1, M
      pi = 3.14159D0
      f0 = Om/(2*pi)
      df = f0
      kfr = 0
      DO i = 1 , m
         fi = P3(1+(2*N**2+1)*(i-1))
         df1 = dabs(f0-fi)
         IF ( df1.LE.df ) THEN
            df = df1
            kfr = i
         ENDIF
      ENDDO
C
C  CALCULATION OF IB (I - BASE)
C
      ib = (kfr-1)*(2*N**2+1) + 1
C
      DO i40 = 1 , N
         DO i50 = 1 , N

C  SEARCH FOR POSITIONS
            ipr = ib + (2*i50-1) + (i40-1)*2*N
            ipi = ipr + 1
C  FILLING OF THE Y-MATRIX
            Y(i50,i40) = dcmplx(P3(ipr),P3(ipi))
C      PRINT 25,I50,I40,Y(I50,I40)
C  25 FORMAT(2X,'Y(',I4,',',I4,')=',E13.6,2X,E13.6)
         ENDDO
      ENDDO
C     DEBUG SUBTRACE,INIT
      END
