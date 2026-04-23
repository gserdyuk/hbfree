!*==LUCAN.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE lucan(Alu,Ntot,N,Nf,Flag)
C**************************************************************
C        SUBROUTINE FOR BRINGING MATRIX A TO A FORM           *
C     CANONICAL FOR THE Y-MATRIX.                             *
C**************************************************************
C
C$LARGE: ALU
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER irow , jcol , N , Nf , Ntot
      DOUBLE COMPLEX Alu(Ntot,Ntot)
      DOUBLE COMPLEX sum
      INTEGER Flag
C  CHECKING INPUT DATA
      Flag = 0
      IF ( Nf.GT.N ) Flag = 5
C  CORRECTIONS INTRODUCED COMPARED TO LUSLV
      IF ( (Nf.GT.Ntot) .OR. (N.GT.Ntot) ) Flag = 3
      IF ( (Nf.LT.1) .OR. (N.LT.1) ) Flag = 2
      IF ( Ntot.LT.1 ) Flag = 1
C      IF(FLAG.NE.0) WRITE(6, 95)  FLAG
      IF ( Flag.NE.0 ) RETURN
C  IF N=1, CANONICALIZATION IS CONSIDERED COMPLETE -
C  THE MATRIX CONTAINS ONLY 1 ELEMENT - ALU(1,1).
      IF ( N.EQ.1 ) RETURN
C  TRANSFER THE ACTIVE PART OF MATRIX A TO THE CANONICAL FORM.
C  DIAGONAL ELEMENTS WILL BE TRANSFERRED TO A.
      DO irow = Nf , N
         sum = 2.D0*Alu(irow,irow)
         DO jcol = Nf , N
            sum = sum - Alu(irow,jcol)
C     WRITE(6, 83) SUM, IROW, JCOL, ALU(IROW,JCOL)
C  83 FORMAT(2X,'LUCAN:SUM=',2(1X,D15.8),'ALU(',I3,I3,')=',2(1X,E12.5))
         ENDDO
         Alu(irow,irow) = sum
C     WRITE(6, 84) IROW, IROW, ALU(IROW,IROW)
C  84 FORMAT(2X,25X,'ALU(',I3,I3,')=',E12.5,2X,E12.5)
      ENDDO
C   95 FORMAT('   LUCAN:'/'   INPUT DATA ERROR. FLAG=',I4)
      END
