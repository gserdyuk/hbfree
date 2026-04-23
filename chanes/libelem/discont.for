!*==DISCONT.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE discont(Om,P1,L1,P2,L2,P3,L3,N)
C
C  ANALOG OF SUBROUTINE SM38 (FREQ, NOL)
C
C  Stripline with Step Width Junction discontinuities (STEP)
C
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION d1 , d2 , ep , epe1 , epe2 , eps , eps1 , eps2 ,
     &                 f , fep1 , fl , h , r , sm63 , t , w1 , w2 ,
     &                 wef1 , wef2 , wn
      DOUBLE PRECISION wn1 , x , z1 , z2
      INTEGER L1 , L2 , L3 , N
      DIMENSION P1(L1) , P2(L2) , P3(L3)

      DOUBLE PRECISION Om , P1 , P2 , P3

      COMMON /subs  / Sv(15,15)

      DOUBLE COMPLEX Sv
      DOUBLE COMPLEX par(4) , param(10) , zp

      LOGICAL*4 syms , syml , nol


      syml = .TRUE.

C     MODEL PARAMETERS: W1, W2, T, H, EPS
C               W1   -   WIDTH OF THE STRIP BEFORE THE STEP (m)
C               W2   -   WIDTH OF THE STRIP AFTER THE STEP (m)
C               T    -   THICKNESS OF THE MPL STRIP (m)
C               H    -   THICKNESS OF THE MPL SUBSTRATE (m)
C               EPS  -   RELATIVE DIELECTRIC PERMITTIVITY
C
C     PARAMETERS COMMON FOR THE SCHEME
      eps = P1(1)
C
C     PARAMETERS COMMON FOR THIS TYPE OF ELEMENT
      h = P2(1)*1.D+03
      t = P2(2)*1.D+03
C
C     INDIVIDUAL PARAMETERS
C     IF FL=1, THE STEP IS ASYMMETRIC,
C     IF FL=-1, THE STEP IS SYMMETRIC
      fl = P3(1)
      w1 = P3(2)*1.D+03
      w2 = P3(3)*1.D+03
      syms = .FALSE.
      IF ( fl.GT.0 ) syms = .TRUE.
      Om = Om*1.D-09


      WRITE (6,5) w1 , w2 , t , h , eps
 5    FORMAT (2X,' STEP. INITIAL DATA:'/2X,'W1=',E12.5,'  W2=',E12.5,2X,
     &        'T=',E12.5,3X,'H=',E12.5/2X,'EPS=',E12.5)
C     IF(IFR.GT.1) GOTO 1
C
C
C     MICROSTRIP LINE
C
      ep = eps

      CALL msl(w1,t,h,ep,wef1,epe1,z1)
      CALL msl(w2,t,h,ep,wef2,epe2,z2)
      WRITE (6,6) z1 , z2 , epe1 , epe2 , wef1 , wef2
 6    FORMAT (2X,' STEP: Z1=',E12.5,'  Z2=',E12.5/2X,' EPE1=',E12.5,
     &        '    EPE2=',E12.5,'  WEF1=',E12.5,'  WEF2=',E12.5)
C   S-MATRIX FOR STEADY-STATE CURRENT MODE
      IF ( Om.NE.0.0D0 ) THEN
         eps1 = dsqrt(epe1)
         eps2 = dsqrt(epe2)
         r = z2/z1
         d1 = sm63(z1,eps1,h)
         d2 = sm63(z2,eps2,h)
         f = 149.89623D0/dmax1(eps1*d1,eps2*d2)

C    WARNING - CHECK THE ASSIGNMENT OF THE NEXT OPERATOR.
C    I MEAN THE CALCULATION OF THE SYMMETRIC STEP.

         IF ( .NOT.syms ) f = 2.D0*f

         wn = .02095845D0*Om
C       IF(SYML) GO TO 70
         CALL sm21(ep,epe1,h,z1,wn,fep1)
         wn1 = wn*dsqrt(fep1)
         CALL sm32(wn1,d1,d2,syms,x,nol)
         IF ( .NOT.nol ) THEN
            zp = dcmplx(0.D0,x)
            CALL sm01(zp,r,1,par)
C
            param(1) = par(1)
            param(2) = par(2)
            param(3) = par(3)

C     FORMATION OF THE S-MATRIX
            Sv(1,1) = param(1)
            Sv(1,2) = param(2)
            Sv(2,1) = param(2)
            Sv(2,2) = param(3)

            WRITE (6,21) Sv(1,1) , Sv(1,2) , Sv(2,1) , Sv(2,2)
 21         FORMAT (2X,'SV(1,1)=',E12.5,1X,E12.5/2X,'SV(1,2)=',E12.5,1X,
     &              E12.5/2X,'SV(2,1)=',E12.5,1X,E12.5/2X,'SV(2,2)=',
     &              E12.5,1X,E12.5)
C
C     FORMATION OF THE Y-MATRIX
            CALL test(N)
            RETURN
         ENDIF
      ELSE
         Sv(1,1) = dcmplx(0.0D0,0.0D0)
         Sv(1,2) = dcmplx(0.99997D0,0.0D0)
         Sv(2,1) = dcmplx(0.99997D0,0.0D0)
         Sv(2,2) = dcmplx(0.0D0,0.0D0)
         RETURN
      ENDIF
      WRITE (6,11) f
 11   FORMAT (10X,'STEP:','F  MAX =',E12.5)
      RETURN
C


      END



      SUBROUTINE sm32(Wng,B1,B2,Asym,B,Nol)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION a , a1 , alfa , alfa1 , B , B1 , B2 , bw , bw1 ,
     &                 c , d , d1 , r , wlg , Wng
      LOGICAL Nol , Asym
      Nol = .FALSE.
      wlg = 6.2831853D0/Wng
      IF ( Asym ) wlg = wlg*.5D0
      IF ( (B1.LE.wlg) .AND. (B2.LE.wlg) ) THEN
         r = B2/B1
         IF ( r.LT.1.D0 ) THEN
            bw = B1
            bw1 = B2
         ELSEIF ( r.EQ.1.D0 ) THEN
            B = 0.D0
            RETURN
         ELSE
            bw = B2
            bw1 = B1
         ENDIF
         alfa = bw1/bw
         bw = bw/wlg
         bw1 = bw1/wlg
         d = (1.D0+alfa)/(1.D0-alfa)
         d1 = d**(1.D0/alfa)
         c = 4.D0*alfa
         d = d**alfa
         alfa = alfa*alfa
         alfa1 = 1.D0/(1.D0-alfa)
         a = dsqrt(1.D0-bw*bw)
         a = (a+1.D0)/(1.D0-a)
         a = a*d*d - (1.D0+3.D0*alfa)*alfa1
         a1 = dsqrt(1.D0-bw1*bw1)
         a1 = (1.D0+a1)/(1.D0-a1)
         a1 = a1*d1*d1 + (3.D0+alfa)*alfa1
         c = c*alfa1
         d1 = d*d1
         B = dlog((1.D0/c)*dsqrt(d1))
         c = c*c
         B = B + 2.D0*(a+a1+2.D0*c)/(a*a1-c*c)
         a1 = ((5.D0*alfa-1.D0)*alfa1+1.3333333D0*alfa*c/a)**2
         alfa = bw1/bw
         d = ((1.D0-alfa)/(1.D0+alfa))**(4.D0*alfa)
         B = B + .0625D0*bw**2*d*a1
         B = B*bw*2.D0
         IF ( r.GT.1.D0 ) B = B*alfa
      ELSE
         Nol = .TRUE.
         RETURN
      ENDIF
      END



      SUBROUTINE sm01(Z,Rat,It,Pars)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER It
      DOUBLE PRECISION Rat , t
      DOUBLE COMPLEX Z , a , d , Pars(3)
      t = dfloat(It)
      a = Z + Rat
      d = 1.D0/(1.D0+a)
      Pars(1) = t*d*(a-1.D0)
      Pars(2) = 2.D0*d*dsqrt(Rat)
      Pars(3) = t*d*(Z-Rat+1.D0)
      END



      SUBROUTINE msl(W,T,H,Er,We,Ee,Z)
C  KOCTYKEBICH A.B., DEC. 1979
C  PROCEEDINGS OF THE IEEE, 1977, V.65, N11, MTT-25
C  BAHL I.J., GARG R. AUTHORS OF THE ARTICLE.
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION a1 , a2 , a3 , a4 , c , ea , eb , ec , Ee , ef ,
     &                 Er , H , T , W , We , wh , Z
      DATA a1 , a2 , a3 , a4/376.9911D0 , .3978873D0 , 12.56637D0 ,
     &     .1591549D0/
      We = W
      wh = W/H
      c = (Er-1.D0)/4.6D0*T/(H*dsqrt(wh))
      ea = (Er+1.D0)*.5D0
      eb = (Er-1.D0)*.5D0
      ec = 1.D0/dsqrt(1.D0+12.D0/wh)
      IF ( T.GE.1.D-5 ) THEN
         ef = .5D0*H/T
         IF ( wh.GT.a4 ) ef = a3*W/T
         We = We + a2/H*T*(1.D0+dlog(ef))
      ENDIF
      IF ( wh.GT.1.D0 ) THEN
C     WIDE STRIP
         Ee = ea - c + eb*ec
         Z = a1/(dsqrt(Ee)*(We/H+1.393D0+.667D0*dlog(We/H+1.444D0)))
         RETURN
      ENDIF
C     NARROW STRIP
      Ee = ea - c + eb*(ec+.04D0*(1.D0-wh)**2)
      Z = 60.D0/dsqrt(Ee)*dlog(8.D0/We*H+.25D0/H*We)
      RETURN
      END




      SUBROUTINE sm21(Ep,Epe,H,Z,Wn,Fep)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Ep , Epe , Fep , g , H , Wn , wnr , Z
      wnr = Z/H/119.916984D0
      g = .6D0 + .009D0*Z
      Fep = Ep - (Ep-Epe)/(1.D0+g*(Wn/wnr)**2)
      END



      SUBROUTINE sm22(W,T,B,Ep,Z)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION a , ak , ak1 , B , e , e1 , Ep , f , f1 , p , T ,
     &                 W , Z
      a = dsin(1.57079633D0/B*T)
      p = .5D0*dlog((1.D0+a)/(1.D0-a))
      a = dexp(1.57079633D0*W/(B-T)+p)
      ak = 2.D0/(a+1.D0/a)
      ak1 = dsqrt(1.D0-ak**2)
C
      CALL sm10(ak1,e,f)
      CALL sm10(ak,e1,f1)
      Z = 94.2477795D0*f/(dsqrt(Ep)*f1)
      END




      SUBROUTINE sm10(Ak1,E,F)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Ak1 , E , F , t
      t = Ak1**2
C
C     ELLIPTIC INTEGRAL OF THE SECOND KIND
C
      F = ((.032024666D0*t+.054544409D0)*t+.097932891D0)
     &    *t + 1.3862944D0 -
     &    (((.010944912D0*t+.060118519D0)*t+.12475074D0)*t+.5D0)*dlog(t)
C
C     ELLIPTIC INTEGRAL OF THE FIRST KIND
C
      E = ((.040905094D0*t+.085099193D0)*t+.44479204D0)*t + 1.D0 -
     &    (((.01382999D0*t+.08150224D0)*t+.24969795D0)*t)*dlog(t)
      END




      DOUBLE PRECISION FUNCTION sm63(Z,Sep,H)
C       CHANGES MADE ON 7.02.91 BY N.K.
C     LOGICAL*4 SYML
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION H , p , Sep , Z
      p = .25D0
C     IF(SYML)P=1.
      sm63 = 376.687D0*p*H/(Z*Sep)
      END
