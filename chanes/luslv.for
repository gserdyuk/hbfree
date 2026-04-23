!*==LUSLV.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE luslv(Alu,Ntot,N,Nf,Nend,Flag)
C**************************************************************
C*          Subroutine for Complex LU Factorization           *
C*   Designed for LU Factorization of part of the matrix ALU. *
C*                                                            *
C*          FORMAL PARAMETERS:                                *
C*    ALU - OUTPUT MATRIX, HAS NORMAL (NON-SINGULAR) FORM FOR *
C*   Y-MATRIX, TYPE: THE OUTPUT IS: L-MATRIX, REPRESENTING A  *
C*   PREVIOUSLY SELECTED MATRIX, WITH L = L'*D, WHERE L' IS A *
C*   UNIT DIAGONAL, Y-MATRIX ALSO HAS A UNIT DIAGONAL, AND THE*
C*   ELEMENT A(NEND+1, NEND+1) BEGINS A DECOMPOSITION MATRIX, *
C*   WHICH CORRESPONDS TO THE GRAPH C WITH EXCLUDED NODES FROM*
C*   1..N.                                                    *
C*   DECOMPOSITION MATRIX MAINTAINS ITS NON-SINGULAR FORM.    *
C*   IN THE DECOMPOSED PART, THE DIAGONAL IS NON-VARIANT.     *
C*   NON-SINGULAR IS CONSIDERED A Y-MATRIX, WHERE DIAGONAL    *
C*   ELEMENTS ARE NON-FULL NODES OF CONDUCTIVITY, ONLY THE    *
C*   CONDUCTIVITY BETWEEN I-M AND ZERO NODES.                 *
C*    COMPLEX ALU(NTOT, NTOT)                                 *
C*      NTOT - FULL DIMENSION OF THE MATRIX ALU (BY THE       *
C*             DIMENSION DECLARED IN THE CALLING PROGRAM)     *
C*      N - DIMENSION OF THE PROCESSED (I.E. FILLED) PART OF ALU *
C*      NF - STARTING ROW, FROM WHICH THE FACTORIZATION SHOULD *
C*           BEGIN                                            *
C*      NEND - ENDING ROW, WHERE FACTORIZATION SHOULD BE ENDED*
C*      FLAG - ERROR FLAG:                                    *
C*        0 - NO ERRORS,                                      *
C*        1 - N, NEND OR NF LESS THAN 1                       *
C*        2 - N, NEND OR NF GREATER THAN NTOT                 *
C*        3 - NF GREATER THAN NEND                            *
C*        4 - NF OR NEND GREATER THAN N                       *
C*                                                            *
C**************************************************************

C$LARGE: ALU
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER irow , jcol , kpiv , kpivp1 , N , ne , Nend , nep1 , Nf ,
     &        Ntot
      DOUBLE COMPLEX Alu(Ntot,Ntot)
      DOUBLE COMPLEX diag(100) , aki , di , sum , diff
      INTEGER Flag
C  CHECKING THE INPUT DATA
C     DO 4 IV=1,N
C     DO 4 JV=1,N
C     WRITE(6, 3) IV, JV, ALU(IV,JV)
C   3 FORMAT(2X,'LUSLV: ALU(',I3,',',I3,')=',E12.5,2X,E12.5)
C   4 CONTINUE
      Flag = 0
      IF ( (Nf.GT.N) .OR. (Nend.GT.N) ) Flag = 5
      IF ( Nf.GT.Nend ) Flag = 4
      IF ( (Nf.GT.Ntot) .OR. (Nend.GT.Ntot) .OR. (N.GT.Ntot) ) Flag = 3
      IF ( (Nf.LT.1) .OR. (Nend.LT.0) .OR. (N.LT.1) ) Flag = 2
      IF ( Ntot.LT.1 ) Flag = 1
C  EXCLUSION: IF NEND=NF-1, THEN FUTURE ELEMENTS - 0
C  (NORMAL TERMINATION. RECOVERY)  (!) FLAG=4 (!)
      IF ( Nend-Nf+1.EQ.0 ) RETURN
      IF ( Flag.NE.0 ) WRITE (6,35) Flag
 35   FORMAT ('    LUSLV:'/'   ERR IN INPUT DATA.FLAG=',I4)
      IF ( Flag.NE.0 ) RETURN
C  IF N=1, THEN THE FACTORIZATION CAN BE CONSIDERED COMPLETED - THE MATRIX
C  CONTAINS ONE ELEMENT ALU(1,1)
      IF ( N.EQ.1 ) RETURN

C  TRANSFORMATION OF THE ACTIVE PART OF THE Y-MATRIX ALU TO A CANONICAL
C  FORM. DIAGONAL ELEMENTS WILL BE STORED IN A SEPARATE VECTOR OF DOUBLE
C  PRECISION. AFTER THE USE OF THE ELEMENT, THE DIAGONAL VECTOR WILL BE
C  REPLACED IN ALU. UNUSED ONES WILL BE CONVERTED TO CANONICAL FORM.
      DO irow = Nf , N
         sum = 2.D0*Alu(irow,irow)
         DO jcol = Nf , N
            sum = sum - Alu(irow,jcol)
         ENDDO
         diag(irow) = sum
      ENDDO
C  FACTORIZATION OF NEND-NF+1 FUTURE ELEMENTS
C  EXCLUSION: IF NEND = N, THEN NEND = NEND - 1
      ne = Nend
      IF ( Nend.EQ.N ) ne = Nend - 1
      DO kpiv = Nf , ne
C  ASSIGNMENT OF THE USED ELEMENT DIAG TO A(KPIV, KPIV)
C     WRITE(6,1111) NF,NE,KPIV,DIAG(KPIV)
C1111 FORMAT(2X,'LUSLV: NF=',I5,' NE=',I5,' KPIV=',I5,' DIAG=',2E12.5)
         Alu(kpiv,kpiv) = diag(kpiv)
C  NORMALIZATION FOR THE ELEMENTS OF THE NEXT ROW
         di = 1.D0/(diag(kpiv)+0.1D-30)
C  IN ALL COLUMNS OF THE ACTIVE SUBMATRIX:
         kpivp1 = kpiv + 1
         DO jcol = kpivp1 , N
C  NORMALIZATION THE ELEMENTS OF THE NEXT ROW
            aki = Alu(kpiv,jcol)*di
            Alu(kpiv,jcol) = aki
C  ASSIGN IT TO THE ELEMENT OF THE NEXT COLUMN
C  AND CALCULATE FROM THE CURRENT ELEMENT.
C  AND SO ON FOR ALL ELEMENTS OF THE COLUMN.
            DO irow = kpivp1 , N
               IF ( irow.NE.jcol ) Alu(irow,jcol) = Alu(irow,jcol)
     &              - aki*Alu(irow,kpiv)
               IF ( irow.EQ.jcol ) diag(irow) = diag(irow)
     &              - aki*Alu(irow,kpiv)
            ENDDO
         ENDDO
      ENDDO
C  TRANSFORMATION OF THE ACTIVE SUBMATRIX PART
C  TO A NORMALIZED FORM.
C
C  FOR ALL ROWS OF THE ACTIVE SUBMATRIX:
      nep1 = ne + 1
      DO irow = nep1 , N
         diff = diag(irow)
C  FOR ALL ELEMENTS OF THE ROW (INCLUDING THE DIAGONAL ELEMENT)
         DO jcol = nep1 , N
            IF ( irow.NE.jcol ) diff = diff + Alu(irow,jcol)
         ENDDO
         Alu(irow,irow) = diff
      ENDDO
      END
