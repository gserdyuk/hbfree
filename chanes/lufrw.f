c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE LUFRW(ALU,B,NTOT,N,NF,NEND,FLAG)
C**************************************************************
C*       SUBROUTINE FORWARD ELIMINATION                      *
C*                                                            *
C*       B - VECTOR OF FREE TERMS (RIGHT-HAND SIDE)          *
C*                                                            *
C**************************************************************
C$LARGE: ALU,B
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX ALU(NTOT,NTOT),B(NTOT)
      INTEGER FLAG
C  CHECKING INPUT DATA
      FLAG=0
      IF((NF.GT.N).OR.(NEND.GT.N)) FLAG=5
      IF(NF.GT.NEND) FLAG=4
      IF((NF.GT.NTOT).OR.(NEND.GT.NTOT).OR.(N.GT.NTOT)) FLAG=3
      IF((NF.LT.1).OR.(NEND.LT.0).OR.(N.LT.1)) FLAG=2
      IF(NTOT.LT.1) FLAG=1
C  EXCEPTION: IF NEND = NF - 1, THEN THERE ARE NO PIVOT ELEMENTS
C  (NORMAL TERMINATION. RETURN)  (!) FLAG = 4 (!)
      IF(NEND-NF+1.EQ.0) RETURN
      IF(FLAG.NE.0) WRITE (6, 55) FLAG
      IF(FLAG.NE.0) RETURN
C  IF N = 1, THEN IT IS NECESSARY TO PERFORM ONLY  B(N) = B(N) / ALU(N,N)
      IF(N.EQ.1) GO TO 45
C  DIRECT ELIMINATION
C  IF NEND = N, THEN NEND = NEND - 1
      NE=NEND
      IF(NEND.EQ.N) NE=NEND-1
      DO 40 JCOL=NF,NE
      B(JCOL)=B(JCOL)/(ALU(JCOL,JCOL)+0.1D-30)
      JCOLP1=JCOL+1
      DO 50 IROW=JCOLP1,N
      B(IROW)=B(IROW)-ALU(IROW,JCOL)*B(JCOL)
   50 CONTINUE
   40 CONTINUE
C  THEREFORE, IN THE PREVIOUS COMMENT, WE USED NE INSTEAD OF NEND.
C  ALSO BECAUSE AN ERROR OCCURS: NEND IS DECREMENTED BY
C  ONE AND THEN PASSED FURTHER ON.
   45 IF(NEND.EQ.N) B(N)=B(N)/ALU(N,N)
      RETURN
   55 FORMAT('  LUFRW: '/'  INPUT DATA ERROR. FLAG=',I4)
      END
