!*==LUFRW.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE lufrw(Alu,B,Ntot,N,Nf,Nend,Flag)
C**************************************************************
C*       SUBROUTINE FORWARD ELIMINATION                      *
C*                                                            *
C*       B - VECTOR OF FREE TERMS (RIGHT-HAND SIDE)          *
C*                                                            *
C**************************************************************
C$LARGE: ALU,B
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER irow , jcol , jcolp1 , N , ne , Nend , Nf , Ntot
      DOUBLE COMPLEX Alu(Ntot,Ntot) , B(Ntot)
      INTEGER Flag
C  CHECKING INPUT DATA
      Flag = 0
      IF ( (Nf.GT.N) .OR. (Nend.GT.N) ) Flag = 5
      IF ( Nf.GT.Nend ) Flag = 4
      IF ( (Nf.GT.Ntot) .OR. (Nend.GT.Ntot) .OR. (N.GT.Ntot) ) Flag = 3
      IF ( (Nf.LT.1) .OR. (Nend.LT.0) .OR. (N.LT.1) ) Flag = 2
      IF ( Ntot.LT.1 ) Flag = 1
C  EXCEPTION: IF NEND = NF - 1, THEN THERE ARE NO PIVOT ELEMENTS
C  (NORMAL TERMINATION. RETURN)  (!) FLAG = 4 (!)
      IF ( Nend-Nf+1.EQ.0 ) RETURN
      IF ( Flag.NE.0 ) WRITE (6,55) Flag
 55   FORMAT ('  LUFRW: '/'  INPUT DATA ERROR. FLAG=',I4)
      IF ( Flag.NE.0 ) RETURN
C  IF N = 1, THEN IT IS NECESSARY TO PERFORM ONLY  B(N) = B(N) / ALU(N,N)
      IF ( N.NE.1 ) THEN
C  DIRECT ELIMINATION
C  IF NEND = N, THEN NEND = NEND - 1
         ne = Nend
         IF ( Nend.EQ.N ) ne = Nend - 1
         DO jcol = Nf , ne
            B(jcol) = B(jcol)/(Alu(jcol,jcol)+0.1D-30)
            jcolp1 = jcol + 1
            DO irow = jcolp1 , N
               B(irow) = B(irow) - Alu(irow,jcol)*B(jcol)
            ENDDO
         ENDDO
      ENDIF
C  THEREFORE, IN THE PREVIOUS COMMENT, WE USED NE INSTEAD OF NEND.
C  ALSO BECAUSE AN ERROR OCCURS: NEND IS DECREMENTED BY
C  ONE AND THEN PASSED FURTHER ON.
      IF ( Nend.EQ.N ) B(N) = B(N)/Alu(N,N)
      END
