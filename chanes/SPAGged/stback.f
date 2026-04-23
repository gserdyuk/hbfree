c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE STBACK(U,NREC,Y,VJ,ISIZE_MAXNODE,S)
C*********************************************************************
C     SUBROUTINE FOR CALCULATING THE SOLUTION OF THE SYSTEM OF EQUATIONS
C     FOR THE POTENTIAL FUNCTION AT THE BOUNDARY NODES
C          (ITERATIVE METHOD FOR LU-DECOMPOSITION TRANSFORMATION)
C*********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$LARGE: VJ,Y

C      COMMON/FORMY/   VJ(70),Y(70,70)
      INTEGER ISIZE_MAXNODE,FLAG
      DOUBLE COMPLEX Y (ISIZE_MAXNODE,ISIZE_MAXNODE)
      DOUBLE COMPLEX VJ(ISIZE_MAXNODE)

      DOUBLE COMPLEX U(1)
C everywhere else U(1) is double precision , but here is double complex
C it is intentionally size is 1/2 of original (sure)

      DOUBLE COMPLEX S(ISIZE_MAXNODE,20)
      COMMON/KOLNAL/  KOL,NAL
      INTEGER         KOL(4),NREC
      LOGICAL         NAL(4)


      K12=KOL(1)+KOL(2)
      K123=K12+KOL(3)

C      print *, 'STBACK: NREC: ', NREC

      CALL ZINY(Y,VJ,ISIZE_MAXNODE)
C      print *, 'STBACK: ZINY PASSED :NREC ', NREC
C      stop
C      WRITE(6,111) K123,NREC,(II,VJ(II),II=1,K123)
c111   FORMAT(2X,'STBACK: K123=',I3,'  NREC=',I3/
c     *(2X,'VJ(',I3,')=',E12.5,2X,E12.5))
C     DECOMPOSITION AND CALCULATION OF THE REDUCED Y-MATRIX

C      print *, 'STBACK: NREC - just before DPACK: ', NREC
      CALL DPACK2(NREC,Y,VJ,ISIZE_MAXNODE)
c      WRITE (6, 111) K123,NREC, (II,VJ(II),II=1,K123)

C      write (6,*) 'Y matrix'
C      do ii=1,k123
C            write (6,120) (Y(ii,jj), jj=1,k123)
C      enddo
C120   format (2x,'(',1x,e12.5,1x,e12.5,')')

C  TRANSFORMATION OF THE POTENTIAL FUNCTION IN THE BOUNDARY SOLUTION S
      DO 40 I=1,K123
   40 S(I,NREC)=VJ(I)
C      write (6,*) 'S: just filled:'
C      WRITE(6, 112) (II,NREC,S(II,NREC),II=1,K123)
C  112 FORMAT(2X,'S(',I3,I3,')=',E12.5,2X,E12.5)
C     TRANSFORMATION OF THE POTENTIAL FUNCTION AT THE BOUNDARY NODES IN 'XBOCT'
C               SOLUTION OF THE EQUATIONS

      IF(KOL(3).EQ.0)GO TO 20
      K3=KOL(3)
      DO 30 IN=1,K3
C     WRITE(6, 25) K12,IN,NREC,NREC,K3,IN, U((NREC-1)*K3+IN)
C  25 FORMAT(2X,'STBACK: K12,IN, NREC =',I4,I4,I4,
C    *          (2X,'U((',I3,'-1)*',I3,'+',I3,')=',
C    *          E12.5,2X,E12.5))
   30 S(K12+IN,NREC)=U((NREC-1)*K3+IN)
C      write (6,*) 'before LUBACK'
C      WRITE(6, 112) (II,NREC,S(II,NREC),II=1,K123)
C     ACCURATE CALCULATION OF THE SOLUTION
   20 CALL LUBACK(Y,S(1,NREC),ISIZE_MAXNODE,K123,K12,FLAG)
C      write (6,*) 'after LUBACK'
c      WRITE(6, 112) (II,NREC,S(II,NREC),II=1,K123)
      RETURN
C     DEBUG SUBTRACE
      END
