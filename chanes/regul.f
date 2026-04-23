c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE REGUL(NK,DFDX,F,N)
C   Subroutine for regulation of the Jacobian matrix.
C$LARGE:DFDX,F
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION  DFDX         , F
      DIMENSION DFDX(NK,NK),F(NK)
      DO 5 I=1,N
      PI=0.D0
      DO 1 J=1,N
    1 PI=PI+DABS(DFDX(I,J))
      IF(PI.EQ.0.D0) GO TO 3
      DO 2 J=1,N
    2 DFDX(I,J)=DFDX(I,J)/PI
      F(I)=F(I)/PI
    3 IF (F(I).NE.0.0D0) GO TO 5
      DFDX(I,I)=DFDX(I,I)*1.01D0
      IF(PI.EQ.0.D0) DFDX(I,I)=1.0D0
    5 CONTINUE

C      n1=n+1
C      do 20 i=1,n
C      print 10, F(I),I,(dfdx(i,j),j=1,n1)
C   10 format(2x,'REGUL : F=',E13.6,' DFDX(',I3,', J)=',(/3X,6(E13.6)))
C   20 continue

      RETURN
      END
