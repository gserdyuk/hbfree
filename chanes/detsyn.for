!*==DETSYN.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE detsyn(Maxsyn)
C
C     DETERMINATION OF THE MAXIMUM NUMBER OF INPUT/OUTPUT VALUES FOR FOURIER TRANSFORMATION
C                  (CALCULATIONS FOR INKOOR)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , iend , ity , Lenntp , Ln , Lp , Maxsyn , Mpoint , nb ,
     &        Nmpnt , Nn , Nnpr , Np , nr
      COMMON /point / Mpoint(1)
      COMMON /pointr/ Nmpnt , Nn , Np , Ln , Lp , Nnpr , Lenntp

C
      nr = 0
      nb = 0
      iend = Nmpnt*20
      DO i = 1 , iend , 20
         ity = Mpoint(i+5)
         IF ( ity.EQ.3 ) THEN
            nb = max0(Mpoint(i+17),nb)
            nr = max0(Mpoint(i+18),nr)
         ENDIF
      ENDDO
      Maxsyn = max0(nr,nb,1)
      END
