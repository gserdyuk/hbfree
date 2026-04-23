!*==STOP0.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE stop0(N,F,X,Sf,Sx,Icode,Kmaxdu)
C
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER j , Kprgrf , Kprlen , Kprlin , Kprnkr , Kprqup , Kprsol ,
     &        Kprsrt , Kprvar , Limit , N
      DOUBLE PRECISION F(1) , X(1)
      DOUBLE PRECISION Sf(1) , Sx(1)
      INTEGER Icode , Kmaxdu

      COMMON /newton/ Epssol , Epsdu , Epsmin , Maxdu , Limit
      DOUBLE PRECISION Epssol , Epsdu , Epsmin , Maxdu
      COMMON /print / Kprlen , Kprsrt , Kprnkr , Kprlin , Kprsol ,
     &                Kprvar , Kprgrf , Kprqup

      DOUBLE PRECISION conv

C     RESET THE COUNTER OF THE NUMBER
C  OF MAXIMUM STEPS TO 0
      Kmaxdu = 0

C     DETERMINATION OF THE BOUNDARY CONDITIONS AT THE INITIAL POINT
      Icode = 0
      conv = 0.0D0
      DO j = 1 , N , 1
         conv = dmax1(conv,Sf((j+1)/2)*dabs(F(j)))
      ENDDO
      IF ( conv.LE.1.D-02*Epssol ) Icode = 1

C  PRINT MESSAGES.
      IF ( Icode.EQ.1 .AND. Kprsol.GT.0 ) WRITE (6,100) conv
      IF ( Icode.EQ.1 ) PRINT 100 , conv

 100  FORMAT (2X,'     INITIAL POINT IS SOLUTION OF '/2X,
     &        ' SIMULTANEOUS EQUATIONS. TOLERANCE:',E13.6)
C     DEBUG SUBTRACE,INIT
      END
