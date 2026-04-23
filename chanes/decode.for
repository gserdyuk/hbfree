!*==DECODE.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE decode(U,Nohar,Isize_maxnode)
C
C  SUBROUTINE FOR PERMUTING ELEMENTS OF THE SOLUTION VECTOR
C  WITHIN THE VECTOR ITSELF TO MATCH THE NODE NUMBERING
C  WITH THE ORIGINAL NUMBERING.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , Nohar
      COMMON /kolnal/ Kol , Nal
      COMMON /nouzl / Nouz , Inouz1

      LOGICAL Nal(4)
      INTEGER Isize_maxnode
      INTEGER Kol(4) , ntot , Nouz(100) , Inouz1(100) , term
      DOUBLE COMPLEX U(1) , ut(Isize_maxnode)

C  TOTAL NUMBER OF NODES:
      ntot = Kol(1) + Kol(2) + Kol(3)
C     WRITE(6, 810) (IVA,NOUZ(IVA),IVA=1,NTOT)
C 810 FORMAT(2X,'DECODE: NOUZ(',I3,')=',I4)
C     WRITE(6, 815) NTOT
C 815 FORMAT(2X,'NTOT=',I5)

C  IN ALL NODES
      DO i = 1 , ntot
C  PERFORM PERMUTATION FROM U TO UT IN ACCORDANCE WITH
C  THE MAPPING VECTOR NOUZ    /SEE SUBROUTINE SORTUZ/
         term = Nouz(i)
         ut(term) = U(i)
      ENDDO
C  NOW COPY FROM THE TEMPORARY VECTOR UT
C  INTO THE SOLUTION VECTOR U:
      DO i = 1 , ntot
C     WRITE(6, 14) I,UT(I)
C  14 FORMAT(2X,'DECODE: U(',I3,')=',E12.5,2X,E12.5)
         U(i) = ut(i)
      ENDDO
C     DEBUG SUBTRACE,SUBCHK,INIT(I,TERM)
      END
