!*==SVUTL.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE svutl(Om,P1,L1,P2,L2,P3,L3,N)
C
C
C     MODEL OF STRIPLINE - COUPLED LINES WITH LOSSES, MICRO-
C     STRIP LINES AND CONNECTION OF LINES
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION alfa , bett , cc , e0 , eps , er , fl , h , pi ,
     &                 rd , skl1 , skl2 , skl3 , Sv , t , tgd , w ,
     &                 zn0 , zo1 , zo2
      INTEGER i , L1 , L2 , L3 , N
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION Om , P1 , P2 , P3

      COMMON /subs  / Sv(15,15)
      DOUBLE PRECISION m0

C

C
C
      eps = P1(1)
      tgd = P1(2)
      rd = P1(4)
      h = P2(1)
      t = P2(2)
      w = P3(1)
      fl = P3(2)
      zo1 = P3(3)
      zo2 = P3(4)
      alfa = P3(5)
      zn0 = P3(6)
      skl1 = P3(7)
      skl2 = P3(8)
      skl3 = P3(9)
C
      er = eps
      e0 = 8.85418D-12
      pi = 3.1415926D0
      cc = 2.99D+8
      m0 = 1.256637D-06

C      SELECTION OF ELEMENTS
C   KZ-STL1=5;
C   XX-STL=4;
C   JUMP-STL=9;
C
C
C
      IF ( skl1.EQ.4.D0 ) THEN

C      OPEN CIRCUIT
         IF ( Om.EQ.0.D0 ) THEN
            Sv(1,1) = 0.0D0
            Sv(2,2) = 0.0D0
            Sv(1,2) = 0.99997D0
            Sv(2,1) = 0.99997D0
            N = 2
            CALL test(N)
         ELSE
            i = -1
            CALL bet11(Om,cc,bett,pi,w,h,er)
            CALL shleyf(i,zo1,bett,fl,w,h,er)
            N = 2
            CALL test(N)
         ENDIF
         RETURN
      ELSEIF ( skl1.EQ.5.D0 ) THEN

C      SHORT-CIRCUIT LINE
         IF ( Om.EQ.0.D0 ) THEN
            Sv(1,1) = .99997D0
            Sv(1,2) = 0.0D0
            Sv(2,1) = 0.0D0
            Sv(2,2) = .99997D0
            N = 2
            CALL test(N)
         ELSE
            i = 1
            CALL bet11(Om,cc,bett,pi,w,h,er)
            CALL shleyf(i,zo1,bett,fl,w,h,er)
            N = 2
            CALL test(N)
         ENDIF
         RETURN
      ELSEIF ( skl1.EQ.9.D0 ) THEN

C      JUMP IN WAVE IMPEDANCE
         CALL stepsc(zo1,zo2)
         N = 2
         CALL test(N)
         RETURN
      ENDIF
      PRINT 40 , skl1
 40   FORMAT ('  INCORRECT ELEMENT CODE :',E12.5,'STOP')
      STOP
      END

      SUBROUTINE bet11(Om,Cc,Bett,Pi,W,H,Er)
C      IF(OM.EQ.0.0) THEN
C     BETT=0.
C      ELSE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Bett , Cc , Er , ere , H , olam , Om , Pi , W
      CALL ereff(W,H,Er,ere)
      olam = Cc*2.D0*Pi/Om
      Bett = 2.D0*Pi*dsqrt(ere)/olam
C      ENDIF
      END

C      CALCULATION OF THE LINE
      SUBROUTINE shleyf(I,Z0s,Bett,Fl,W,H,Er)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Bett , Er , Fl , H , Sv , W , z , Z0s
      INTEGER I
      COMMON /subs  / Sv(15,15)
      DOUBLE COMPLEX st , ds
      st = (0.D0,1.D0)
      CALL z0mpl(z,W,H,Er)
      ds = st*2.D0*z/(dtan(Bett*Fl))**I/Z0s - I
      Sv(1,1) = I/ds
      Sv(1,2) = (ds+I)/ds
      Sv(2,1) = Sv(1,2)
      Sv(2,2) = -Sv(1,1)
      END

C       CALCULATION OF WAVE IMPEDANCE
      SUBROUTINE z0mpl(Z,W,H,Er)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Er , ere , H , W , wh , Z
      wh = W/H
      CALL ereff(W,H,Er,ere)
      IF ( wh.GT.1.D0 ) THEN
         Z = 60.D0/dsqrt(ere)/(wh+1.393D0+.667D0*dlog(wh+1.444D0))
         RETURN
      ENDIF
      Z = 376.7D0/dsqrt(ere)*dlog(8.D0/wh+.25D0*wh)
      RETURN
      END

C      CALCULATION OF EFFECTIVE DIELECTRIC CONSTANT
      SUBROUTINE ereff(W,H,Er,Ere)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Er , Ere , H , W , wh
      wh = W/H
      Ere = (Er+1.D0)/2.D0 + (Er-1.D0)/2.D0/dsqrt(1.D0+10.D0*wh)
      END

C      JUMP IN WAVE IMPEDANCE
      SUBROUTINE stepsc(Z1,Z2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION ds , Sv , Z1 , Z2
      COMMON /subs  / Sv(15,15)
      ds = Z1 + Z2
      Sv(1,1) = (Z2-Z1)/ds
      Sv(1,2) = 2*dsqrt(Z1*Z2)/ds
      Sv(2,1) = Sv(1,2)
      Sv(2,2) = -Sv(1,1)
      END
