c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      INTEGER function krdchk(MN,KR,KC,NNR,KNR,KN,KNC,IR1,IR2)
C
C     WRITTEN ON 12.05.91    Serdyuk G.V.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER    KR(IR2),KC(IR1),NNR(IR1)
      INTEGER    MN(2,IR2)

      integer err,sum1,sum2,sum3

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

      err=0


C    CHECKING THE NUMBER OF ROWS:
c           /                  \
c   i=IR1  |  1 if NNR(i) != 0  |
C    SUM  <                      >  == KNR  IF NOT, THEN err = err + 1
c    i=1   |  0 if NNR(i) == 0  |
c           \                  /
      sum1=0
      do 10 i=1,IR1
 10   if (NNR(i).ne.0) sum1=sum1+1

      if (sum1.ne.KNR) err=err+1

C     GENERAL CHECK OF ELEMENTS IN ROWS
c                        /                 \
c    i=IR1       i=IR2  |  1 if KR(i) != 0  |
C    SUM KC(i) == SUM   <                     >  IF NOT, THEN err = err + 2
c    i=1         i=1    |  0 if KR(i) == 0  |
c                        \                 /
      sum2=0
      do 20 i=1,IR1
 20   sum2=sum2+KC(i)

      sum3=0
      do 30 i=1,IR2
 30   if( kr(i).ne.0) sum3=sum3+1

      if (sum2.ne.sum3) err=err+2


C     CHECKING ORDERING IN ROWS (MAY NOT BE ENFORCED !!!)
c   for j=1 to KNR
c
C  1 + SUM {1}   =  KC(j)     IF NOT, THEN err = err + 4
c  KR(i+1) > KR(i)
c   KR(i) > 0
c
c  end_for
C                -------- !!!  NOT IMPLEMENTED AT THE MOMENT !!!

      krdchk=err
      return
      end
