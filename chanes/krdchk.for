!*==KRDCHK.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
C ++++++++++++++++++++++++++++++++++++++++++++++++
C  A POSSIBLE VARIANT OF THE SUBROUTINE/PROGRAM:
C                       (IF THERE IS TIME - TRY IT)
C
C      SUBROUTINE KOORD(MN,KR,KC,NNR,KNR,KN,KNC,IR1,IR2)
CC
CC  *****  FILLING THE COORDINATE ARRAYS
CC
C      INTEGER    KR(IR2),KC(IR1),NNR(IR1)
C      INTEGER    MN(2,IR2)
C
C      DO 10 I=1,IR1
C      KC(I)=0
C   10 NNR(I)=0
C      DO 15 I=1,IR2
C   15 KR(I)=0
C
C      KNR=1
C      MN11=MN11
C      IF(MN11.LT.0)MN11=MN11+KNC
C      NNR(1)=MN11+1
C
C      DO 40 I=1,KN
C
C      MN1I=MN(1,I)
C      MN2I=MN(2,I)
C      IF(MN1I.LT.0)MN1I=MN1I+KNC
C      IF(MN2I.LT.0)MN2I=MN2I+KNC
C
C      IF( NNR(KNR).NE.(MN1I+1)) THEN
C                                KNR=KNR+1
C                                NNR(KNR)=MN1I+1
C      ENDIF
C
C      KC(KNR)=KC(KNR)+1
C      KR(I)=MN2I+1
C 40   CONTINUE
C
C      RETURN
C      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                             Serdyuk G.V.   12.05.91.
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      INTEGER FUNCTION krdchk(Mn,Kr,Kc,Nnr,Knr,Kn,Knc,Ir1,Ir2)
C
C     WRITTEN ON 12.05.91    Serdyuk G.V.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , Ir1 , Ir2 , Kn , Knc , Knr
      INTEGER Kr(Ir2) , Kc(Ir1) , Nnr(Ir1)
      INTEGER Mn(2,Ir2)

      INTEGER err , sum1 , sum2 , sum3

C     CHECKING THE REPRESENTATION OF SPARSE MATRICES FOR
C     STORING RESULTS
C
C     LIST OF FORMAL PARAMETERS - SAME AS IN KOORD
C
C     IDENTIFIERS OF VARIABLES AND ARRAYS:
C
C          KNC - LENGTH OF ONE-DIMENSIONAL FFT;
C          NNR - ARRAY OF ROW NUMBERS OF NONZERO ELEMENTS;
C          KNR - NUMBER OF NONZERO ROWS;
C          KR  - ARRAY OF COLUMN NUMBERS OF NONZERO
C                ELEMENTS IN EACH NONZERO ROW;
C          KC  - ARRAY CONTAINING THE NUMBER OF
C                NONZERO ELEMENTS IN EACH ROW;

      err = 0


C    CHECKING THE NUMBER OF ROWS:
c           /                  \
c   i=IR1  |  1 if NNR(i) != 0  |
C    SUM  <                      >  == KNR  IF NOT, THEN err = err + 1
c    i=1   |  0 if NNR(i) == 0  |
c           \                  /
      sum1 = 0
      DO i = 1 , Ir1
         IF ( Nnr(i).NE.0 ) sum1 = sum1 + 1
      ENDDO

      IF ( sum1.NE.Knr ) err = err + 1

C     GENERAL CHECK OF ELEMENTS IN ROWS
c                        /                 \
c    i=IR1       i=IR2  |  1 if KR(i) != 0  |
C    SUM KC(i) == SUM   <                     >  IF NOT, THEN err = err + 2
c    i=1         i=1    |  0 if KR(i) == 0  |
c                        \                 /
      sum2 = 0
      DO i = 1 , Ir1
         sum2 = sum2 + Kc(i)
      ENDDO

      sum3 = 0
      DO i = 1 , Ir2
         IF ( Kr(i).NE.0 ) sum3 = sum3 + 1
      ENDDO

      IF ( sum2.NE.sum3 ) err = err + 2


C     CHECKING ORDERING IN ROWS (MAY NOT BE ENFORCED !!!)
c   for j=1 to KNR
c
C  1 + SUM {1}   =  KC(j)     IF NOT, THEN err = err + 4
c  KR(i+1) > KR(i)
c   KR(i) > 0
c
c  end_for
C                -------- !!!  NOT IMPLEMENTED AT THE MOMENT !!!

      krdchk = err
      END
