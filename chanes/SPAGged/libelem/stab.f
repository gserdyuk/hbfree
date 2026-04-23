c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE STAB(OM,P1,L1,P2,L2,P3,L3,N)
C *** SUBROUTINE FOR FORMATION OF THE S-MATRIX FROM GIVEN
C     TASK PARAMETERS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/SUBS/     S(15,15)
      DOUBLE COMPLEX           S,J,Y
      DOUBLE PRECISION              P1(L1),P2(L2),P3(L3)

C  DETERMINATION OF THE NUMBERS OF THE GIVEN DISCRETE POINTS
      M=L3/(2*N**2+1)
      DO 10 I10=1,N
      DO 10 J10=1,N
   10 S(I10,J10)=DCMPLX(0.D0,0.D0)
C  WE ASSUME THAT FI IS CLOSER TO F0
C  THAN ANY FJ FOR J = 1, M
      PI=3.14159D0
      F0=OM/(2*PI)
      DF=F0
      KFR=0
      DO 20 I=1,M
      FI=P3(1+(2*N**2+1)*(I-1))
      DF1=DABS(F0-FI)
       IF(DF1.GT.DF) GO TO 20
       DF=DF1
       KFR=I
   20  CONTINUE
C
C  CALCULATION OF IB (I - BASE)
C
       IB=(KFR-1)*(2*N**2+1)+1
C
       DO 40 I40=1,N
       DO 50 I50=1,N

C  SEARCH FOR POSITIONS
       IPR=IB+(2*I50-1)+(I40-1)*2*N
       IPI=IPR+1
C  FILLING OF THE S-MATRIX
       S(I50,I40)=DCMPLX(P3(IPR),P3(IPI))
C      PRINT 25,I50,I40,S(I50,I40)
C  25 FORMAT(2X,'S(',I4,',',I4,')=',E13.6,2X,E13.6)
   50  CONTINUE
   40  CONTINUE
C
C
C *****************************
       CALL TEST(N)
C *****************************







       RETURN
C     DEBUG SUBTRACE,INIT
       END
