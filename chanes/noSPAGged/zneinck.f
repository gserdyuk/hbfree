c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE NEINCK(N,EPSIM,IDIM,SCALEF,SCALEU,U,ICODE)

C     Performs validation and initial setup
C     of the parameters for the SOLVE program

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/TYPVAL/   TYPU,TYPI
      DOUBLE PRECISION             TYPU,TYPI
      DOUBLE PRECISION             SCALEF(1),SCALEU(1)
      COMMON/NEWTON/   EPSSOL,EPSDU,EPSMIN,MAXDU,LIMIT
      DOUBLE PRECISION             EPSSOL,EPSDU,EPSMIN,MAXDU
      COMMON/PRINT /   KPRLEN,KPRSRT,KPRNKR,KPRLIN,KPRSOL,KPRVAR,
     +          KPRGRF,KPRQUP
      INTEGER          LIMMAX/1000/
      DOUBLE PRECISION             U(1)
      DOUBLE PRECISION             NORSTP,NOR1

C  Validation of the value of N - the dimensionality of the system
      IF (N.GE.1) GO TO 5
      ICODE=-1
      RETURN
    5 IF(N.LE.IDIM) GO TO 10
      ICODE=-2
      RETURN
   10 CONTINUE
C Checking the scale settings
      IF(TYPU.LE.0.D0) TYPU=1.D0
      IF(TYPI.LE.0.D0) TYPI=1.D0
      NDIV2=N/2
      DO 20 I=1,NDIV2
      SCALEF(I)=1.D0/TYPI
   20 SCALEU(I)=1.D0/TYPU
C Setting the parameters of the method
C  EPSSOL
      IF(EPSSOL.LE.0) EPSSOL=EPSIM**(1.D0/3.D0)
C  EPSDU
      IF(EPSDU.LE.0) EPSDU=EPSIM**(2.D0/3.D0)
C  EPSMIN
      IF(EPSMIN.LE.0) EPSMIN=EPSIM**(2.D0/3.D0)
C  MAXDU
      IF(MAXDU.GT.0.D0) GO TO 30
      NORSTP=0.D0
      NOR1=0.D0
      DO 40 JJ=1,N
      NORSTP=NORSTP+(SCALEU((JJ+1)/2)*DABS(U(JJ)))**2
      NOR1=NOR1+(SCALEU((JJ+1)/2))**2
   40 CONTINUE
      NORSTP=DSQRT(NORSTP)
      NOR1=DSQRT(NOR1)
      MAXDU=1000.D0*DMAX1(NORSTP,NOR1)
   30 CONTINUE
C  Limit on the number of iterations
      IF(LIMIT.LE.0.OR.LIMIT.GT.LIMMAX) LIMIT=LIMMAX
      RETURN
C     DEBUG SUBCHK,INIT
      END
