!*==STBACK.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE stback(U,Nrec,Y,Vj,Isize_maxnode,S)
C*********************************************************************
C     SUBROUTINE FOR CALCULATING THE SOLUTION OF THE SYSTEM OF EQUATIONS
C     FOR THE POTENTIAL FUNCTION AT THE BOUNDARY NODES
C          (ITERATIVE METHOD FOR LU-DECOMPOSITION TRANSFORMATION)
C*********************************************************************

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , in , k12 , k123 , k3
c$LARGE: VJ,Y

C      COMMON/FORMY/   VJ(70),Y(70,70)
      INTEGER Isize_maxnode , flag
      DOUBLE COMPLEX Y(Isize_maxnode,Isize_maxnode)
      DOUBLE COMPLEX Vj(Isize_maxnode)

      DOUBLE COMPLEX U(1)
C everywhere else U(1) is double precision , but here is double complex
C it is intentionally size is 1/2 of original (sure)

      DOUBLE COMPLEX S(Isize_maxnode,20)
      COMMON /kolnal/ Kol , Nal
      INTEGER Kol(4) , Nrec
      LOGICAL Nal(4)


      k12 = Kol(1) + Kol(2)
      k123 = k12 + Kol(3)

C      print *, 'STBACK: NREC: ', NREC

      CALL ziny(Y,Vj,Isize_maxnode)
C      print *, 'STBACK: ZINY PASSED :NREC ', NREC
C      stop
C      WRITE(6,111) K123,NREC,(II,VJ(II),II=1,K123)
c111   FORMAT(2X,'STBACK: K123=',I3,'  NREC=',I3/
c     *(2X,'VJ(',I3,')=',E12.5,2X,E12.5))
C     DECOMPOSITION AND CALCULATION OF THE REDUCED Y-MATRIX

C      print *, 'STBACK: NREC - just before DPACK: ', NREC
      CALL dpack2(Nrec,Y,Vj,Isize_maxnode)
c      WRITE (6, 111) K123,NREC, (II,VJ(II),II=1,K123)

C      write (6,*) 'Y matrix'
C      do ii=1,k123
C            write (6,120) (Y(ii,jj), jj=1,k123)
C      enddo
C120   format (2x,'(',1x,e12.5,1x,e12.5,')')

C  TRANSFORMATION OF THE POTENTIAL FUNCTION IN THE BOUNDARY SOLUTION S
      DO i = 1 , k123
         S(i,Nrec) = Vj(i)
      ENDDO
C      write (6,*) 'S: just filled:'
C      WRITE(6, 112) (II,NREC,S(II,NREC),II=1,K123)
C  112 FORMAT(2X,'S(',I3,I3,')=',E12.5,2X,E12.5)
C     TRANSFORMATION OF THE POTENTIAL FUNCTION AT THE BOUNDARY NODES IN 'XBOCT'
C               SOLUTION OF THE EQUATIONS

      IF ( Kol(3).NE.0 ) THEN
         k3 = Kol(3)
         DO in = 1 , k3
C     WRITE(6, 25) K12,IN,NREC,NREC,K3,IN, U((NREC-1)*K3+IN)
C  25 FORMAT(2X,'STBACK: K12,IN, NREC =',I4,I4,I4,
C    *          (2X,'U((',I3,'-1)*',I3,'+',I3,')=',
C    *          E12.5,2X,E12.5))
            S(k12+in,Nrec) = U((Nrec-1)*k3+in)
         ENDDO
      ENDIF
C      write (6,*) 'before LUBACK'
C      WRITE(6, 112) (II,NREC,S(II,NREC),II=1,K123)
C     ACCURATE CALCULATION OF THE SOLUTION
      CALL luback(Y,S(1,Nrec),Isize_maxnode,k123,k12,flag)
C      write (6,*) 'after LUBACK'
c      WRITE(6, 112) (II,NREC,S(II,NREC),II=1,K123)
C     DEBUG SUBTRACE
      END
