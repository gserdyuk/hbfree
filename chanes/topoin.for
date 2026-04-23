!*==TOPOIN.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE topoin(Kolpol,Nomer)
C
C     SUBROUTINE RECORDS IN THE ARRAY, STORING INFORMATION ABOUT THE SYSTEM,
C     INCLUDING DATA ON FORMATTING AND RECORDING TO THE MULTIPLEXER
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , j , k , Kolpol , l , Nomer
      INCLUDE 'circuit.i'
      INTEGER nam(6)/'NPOL' , 'REDU' , 'TERM' , 'NOMT' , 'SIMB' ,
     &        'ADDR'/

      i = (Nmpnt+Nomer-1)*20 + 1
C     STORAGE IN MPOINT
      Mpoint(i) = nam(1)
      Mpoint(i+1) = nam(2)
      Mpoint(i+2) = Nomer
      Mpoint(i+3) = Nomer
      Mpoint(i+4) = 1
      Mpoint(i+5) = 4
      Mpoint(i+6) = Kolpol + 8
      Mpoint(i+7) = 0
      Mpoint(i+8) = 0
      Mpoint(i+9) = Nnode + (Nomer-1)*Mpoint(Nmpnt*20+7)
      Mpoint(i+10) = 0
      Mpoint(i+11) = 0
      Mpoint(i+12) = 0
      Mpoint(i+13) = 0
      Mpoint(i+14) = Nomer - 1
      Mpoint(i+15) = 0
      Mpoint(i+16) = 0
      Mpoint(i+17) = 0
      Mpoint(i+18) = 0
      Mpoint(i+19) = 0

      k = Mpoint(i+9)
C     STORAGE IN NODEEL
      Nodeel(k) = nam(3)
      Nodeel(k+1) = Kolpol
      Nodeel(k+2) = Kolpol*(Kolpol+1)
      Nodeel(k+3) = nam(4)
      l = k + 3
      DO j = 1 , Kolpol
         Nodeel(l+j) = j
      ENDDO
      l = l + Kolpol
      Nodeel(l+1) = nam(5)
      Nodeel(l+2) = nam(6)
      Nodeel(l+3) = Nomer - 1
      Nodeel(l+4) = 0
C      DEBUG INIT,SUBTRACE
      END
