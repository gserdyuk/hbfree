!*==LUBACK.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE luback(Alu,B,Ntot,N,Nend,Flag)
C**************************************************************
C*             Subroutine for Path Traversal                  *
C*                                                            *
C*          Vector B - both input and output                  *
C**************************************************************
C$LARGE: ALU, B
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER irow , irowlm , jc , jcol , N , Nend , nm1 , Ntot
      DOUBLE COMPLEX Alu(Ntot,Ntot) , B(Ntot)
      INTEGER Flag
      Flag = 0
      IF ( N.GT.Ntot ) Flag = 3
      IF ( N.LT.1 ) Flag = 2
      IF ( Ntot.LT.1 ) Flag = 1
      IF ( Flag.NE.0 ) WRITE (6,65) Flag
 65   FORMAT ('  LUBACK:'/'INPUT DATA ERROR ,FLAG=',I4)
      IF ( Flag.NE.0 ) RETURN
C If N=1, then path traversal exclusion is not necessary
      IF ( N.EQ.1 ) RETURN
C Since U has a unique diagonal, U(N,N)*X(N) = B(N)
C  X(N)=B(N)
      nm1 = N - 1
      DO jc = 1 , nm1
         jcol = N - jc + 1
C This changes the organization of the loop JCOL = N, 2
         irowlm = jcol - 1
         IF ( irowlm.GT.Nend ) irowlm = Nend
         DO irow = 1 , irowlm
C     WRITE(6, 64) IROW, B(IROW)
C  64 FORMAT(2X,'LUBAC : S(',I3,')=',E12.5,2X,E12.5)
            B(irow) = B(irow) - Alu(irow,jcol)*B(jcol)
C     WRITE(6, 66) IROW,JCOL,ALU(IROW,JCOL),
C    *                  JCOL,B(JCOL),
C    *             IROW,     B(IROW)
C  66 FORMAT(2X,'        Y(',I3,I3,')=',E12.5,2X,E12.5,
C    *' S(',I3,')=',E12.5,2X,E12.5/
C    *2X,'         S(',I3,')=',E12.5,2X,E12.5)
         ENDDO
      ENDDO

      END
