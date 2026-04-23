c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE LUBACK(ALU,B,NTOT,N,NEND,FLAG)
C**************************************************************
C*             Subroutine for Path Traversal                  *
C*                                                            *
C*          Vector B - both input and output                  *
C**************************************************************
C$LARGE: ALU, B
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX ALU(NTOT,NTOT),B(NTOT)
      INTEGER FLAG
      FLAG=0
      IF(N.GT.NTOT) FLAG=3
      IF(N.LT.1) FLAG=2
      IF(NTOT.LT.1)  FLAG=1
      IF(FLAG.NE.0) WRITE(6, 65)   FLAG
      IF(FLAG.NE.0)RETURN
C If N=1, then path traversal exclusion is not necessary
      IF(N.EQ.1) RETURN
C Since U has a unique diagonal, U(N,N)*X(N) = B(N)
C  X(N)=B(N)
      NM1=N-1
      DO 60 JC=1,NM1
      JCOL=N-JC+1
C This changes the organization of the loop JCOL = N, 2
      IROWLM=JCOL-1
      IF(IROWLM.GT.NEND)IROWLM=NEND
      DO 70 IROW=1,IROWLM
C     WRITE(6, 64) IROW, B(IROW)
C  64 FORMAT(2X,'LUBAC : S(',I3,')=',E12.5,2X,E12.5)
      B(IROW)=B(IROW)-ALU(IROW,JCOL)*B(JCOL)
C     WRITE(6, 66) IROW,JCOL,ALU(IROW,JCOL),
C    *                  JCOL,B(JCOL),
C    *             IROW,     B(IROW)
C  66 FORMAT(2X,'        Y(',I3,I3,')=',E12.5,2X,E12.5,
C    *' S(',I3,')=',E12.5,2X,E12.5/
C    *2X,'         S(',I3,')=',E12.5,2X,E12.5)
  70  CONTINUE
  60  CONTINUE
      RETURN
  65  FORMAT('  LUBACK:'/'INPUT DATA ERROR ,FLAG=',I4)

      END
