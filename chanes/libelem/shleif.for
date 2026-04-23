!*==SHLEIF.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE shleif(Om,P1,L1,P2,L2,P3,L3,N,Name)
C     7.02.91
C
C  FORMER SUBROUTINE SM61(FREQ,OUT)
C
C  SHORT-CIRCUITED LINE, XX
C
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION a , as , b , d1 , d2 , ep , epe1 , epe2 , eps ,
     &                 eps1 , eps2 , f , fep1 , fep2 , h , sigm , sm63 ,
     &                 t , tep , tgd
      DOUBLE PRECISION w1 , w2 , wef1 , wef2 , wn , wn1 , wn2 , xa ,
     &                 xb , xc , xd , z1 , z2 , zs2
      INTEGER j , L1 , L2 , L3 , N , nol
      DIMENSION P1(L1) , P2(L2) , P3(L3)

      DOUBLE PRECISION Om , P1 , P2 , P3
      DOUBLE PRECISION l , mu

      COMMON /subs  / Sv(15,15)

      DOUBLE COMPLEX Sv
      DOUBLE COMPLEX x1 , x2 , x3 , x4 , u , p , r , y , z , gam
      DOUBLE COMPLEX par(4) , param(10)


      CHARACTER*4 Name , nkz , nxx

      DOUBLE COMPLEX ccexp
      ccexp(z) = dexp(dreal(z))*dcmplx(dcos(dimag(z)),dsin(dimag(z)))
C
C
      DATA nkz/'KZ  '/ , nxx/'XX  '/

C
C     MODEL PARAMETERS: W1, W2, T, H, EPS, TGD, SIGMA, L
C               W1   -   LINE WIDTH (M)
C               W2   -   LINE WIDTH (M)
C               T    -   STRIP THICKNESS OF MICROSTRIP LINE (M)
C               H    -   SUBSTRATE THICKNESS OF MICROSTRIP LINE (M)
C               EPS  -   RELATIVE DIELECTRIC PERMITTIVITY
C               TGD  -   TANGENT OF DIELECTRIC LOSS ANGLE
C               SIGM -   SPECIFIC CONDUCTIVITY OF THE STRIP (S/M)
C               L    -   LINE LENGTH (M)
C
C     PARAMETERS COMMON TO THE ENTIRE CIRCUIT
      eps = P1(1)
      tgd = P1(2)
      mu = P1(3)
      sigm = P1(4)/1.D+03
C
C     PARAMETERS COMMON TO THIS TYPE OF ELEMENT
      h = P2(1)*1.D+03
      t = P2(2)*1.D+03
C
C     INDIVIDUAL PARAMETERS
      w1 = P3(1)*1.D+03
      w2 = P3(2)*1.D+03
      l = P3(3)*1.D+03
C
      Om = Om*1.D-09
      u = (1.D0,0.D0)
      r = (-0.93584D0,0.D0)


      WRITE (6,5) w1 , w2 , t , h , eps , tgd , sigm , l , N
 5    FORMAT (2X,' LINE. INITIAL DATA:'/2X,'W OCH=',E12.5,'  W SCHL=',
     &        E12.5,2X,'T=',E12.5,3X,'H=',E12.5/2X,'EPS=',E12.5,3X,
     &        'TANG.DELTA=',E14.6,/,3X,'SIGM=',E14.6,3X,'L ым =',E12.5,
     &        1X,'N=',I5)
C     IF(IFR.GT.1) GOTO 1
C
C
C     MICROWAVE TRANSMISSION LINE
C
      ep = eps

      CALL msl(w1,t,h,ep,wef1,epe1,z1)
      CALL msl(w2,t,h,ep,wef2,epe2,z2)
      a = wef1/h
      b = wef2/h
      CALL sm22(w2,t,2.D0*h,ep,zs2)
      WRITE (6,6) z1 , z2
 6    FORMAT (2X,' Z OCH=',E12.5,'  Z ым=',E12.5)
      eps1 = dsqrt(epe1)
      eps2 = dsqrt(epe2)
      d1 = sm63(z1,eps1,h)
      d2 = sm63(z2,eps2,h)
      f = 149.89623D0/dmax1(eps1*d1,eps2*d2)
C   S - MATRIX FOR DC MODE
      IF ( Om.NE.0.0D0 ) THEN


         wn = .02095845D0*Om*1.D-09
         CALL sm21(ep,epe1,h,z1,wn,fep1)
         CALL sm21(ep,epe2,h,z2,wn,fep2)
         WRITE (6,9) fep1 , wn
 9       FORMAT (2X,'SHLEIF : FEP1, WN',E12.5,1X,E12.5)
         wn1 = wn*dsqrt(fep1)
         wn2 = wn*dsqrt(fep2)
         CALL sm36(wn1,d1,d2,xa,xb,xc,xd,nol)
C
C     PER UNIT LENGTH REFLECTION PARAMETERS
C
C     FOR MICROWAVE TRANSMISSION LINE
C
         CALL sm84(wn,w2,t,h,ep,tep,epe2,z2,zs2,b,sigm,as,gam)
         IF ( Name.EQ.nxx ) CALL sm35(w2,h,fep2,z2,wn2,r)
         p = r*ccexp(-2.D0*gam)
         z = (u+p)/(u-p)
         x2 = dcmplx(0.D0,-xb)
         x3 = dcmplx(0.D0,-xc)
         x4 = dcmplx(0.D0,xd)
         z = z*z2/z1
         z = (z+x4)*x3/(z+x3+x4) + x2
C   NORMALIZATION OF Z-WAVES
         z = z*50.D0

         WRITE (6,999) x2 , x3 , x4 , z
 999     FORMAT (2X,'X2=',2(1X,E12.5)/,'  X3=',2(1X,E12.5)/'  X4=',
     &           2(1X,E12.5)/'   Z=',2(1X,E12.5))

C     PARAMETERS OF ELEMENTS OF THE CIRCUIT
C
         x1 = dcmplx(0.D0,xa)
         x1 = x1*50.D0
         y = u/z
C
C     S-PARAMETER
C
C
C   ===============================

         CALL sm03(x1,y,x1,1.D0,1,par)
         DO j = 1 , 3
            param(j) = par(j)
         ENDDO
C
C    FORMATION OF THE S-MATRIX
C      L=1
C      DO 20 I=1,2
C      DO 20 J=1,2
C      CALL SM16(2,2,I,J,IJ,OUT)
C      K=2*IJ
C      SV(L,1)=PARAM(K-1)
C      SV(L,2)=PARAM(K)
C      WRITE(6,19) I,J,IJ,K,SV(L,1),SV(L,2)
C   19 FORMAT(2X,'SHLEIF : I,J,IJ,K =',4I5/2X,
C     *          'SV1=',E12.5,1X,E12.5,'  SV2=',E12.5,1X,E12.5)
C      L=L+1
C   20 CONTINUE
C
         Sv(1,1) = param(1)
         Sv(1,2) = param(2)
         Sv(2,1) = param(2)
         Sv(2,2) = param(3)
C
         WRITE (6,21) Sv(1,1) , Sv(1,2) , Sv(2,1) , Sv(2,2)
 21      FORMAT (2X,'SV(1,1)=',E12.5,1X,E12.5/2X,'SV(1,2)=',E12.5,1X,
     &           E12.5/2X,'SV(2,1)=',E12.5,1X,E12.5/2X,'SV(2,2)=',E12.5,
     &           1X,E12.5)

C
C    CONVERSION FROM S-MATRIX TO Y-MATRIX

         CALL test(N)
         RETURN
      ENDIF
      Sv(1,1) = dcmplx(0.0D0,0.0D0)
      Sv(1,2) = dcmplx(0.99997D0,0.0D0)
      Sv(2,1) = dcmplx(0.99997D0,0.0D0)
      Sv(2,2) = dcmplx(0.0D0,0.0D0)
      RETURN
C



      END


C      SUBROUTINE MSL(W,T,H,ER,WE,EE,Z)
C   KOSTYUKOVICH A.B., DEK. 1979
C  PROCEEDINGS OF THE IEEE, 1977, V.65, N11, MTT-25
C  BAHL I.J., GARG R. AUTHORS OF THE ARTICLE.
C     DATA A1,A2,A3,A4/376.9911,.3978873,12.56637,.1591549/
C      WE=W
C      WH=W/H
C      C=(ER-1.)/4.6*T/(H*SQRT(WH))
C         EA=(ER+1.)*.5
C      EB=(ER-1.)*.5
C      EC=1./SQRT(1.+12./WH)
C      IF(T.LT.1.E-5)GO TO 2
C      EF=.5*H/T
C      IF(WH.GT.A4)EF=A3*W/T
C      WE=WE+A2/H*T*(1.+ALOG(EF))
C    2 IF(WH.GT.1.) GO TO 1
CC  NARROW STRIP
C      EE=EA-C+EB*(EC+.04*(1.-WH)**2)
C      Z=60./SQRT(EE)*ALOG(8./WE*H+.25/H*WE)
C      RETURN
CC  WIDE STRIP
C    1 EE=EA-C+EB*EC
C      Z=A1/(SQRT(EE)*(WE/H+1.393+.667*ALOG(WE/H+1.444)))
C      RETURN
C      END

C      SUBROUTINE SM21(EP,EPE,H,Z,WN,FEP)
C      WNR=Z/H/119.916984
c      G=.6+.009*Z
C      FEP=EP-(EP-EPE)/(1.+G*(WN/WNR)**2)
C      RETURN
C      END
      SUBROUTINE sm25(Zs,Zm,H,Wn,Zf)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION g , H , Wn , wnr , Zf , Zm , Zs , zt
      zt = 2.D0*Zs
      wnr = Zm/H/119.916984D0
      g = .6D0 + .009D0*Zm
      Zf = zt + (Zm-zt)/(1.D0+g*(Wn/wnr)**2)
      END


      SUBROUTINE sm26(W,T,B,Ep,Tep,Rs,Z0,Wn,Alfc,Alfd)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION a , Alfc , Alfd , B , do , dsum , Ep , g , pi ,
     &                 Rs , s , T , Tep , W , Wn , x , Z0
      DATA pi/3.14159265D0/
      Alfc = 0.D0
      IF ( T.NE.0 ) THEN
         g = W/(B-T)
         IF ( g.LT.0.35D0 ) THEN
            x = T/W
            do = ((((.5343325D0*x-1.817245D0)*x+2.318656D0)*x-1.410188D0
     &           )*x+1.059607D0)*x + .4999405D0
            do = do*W
            dsum = ((((((((-17.18875D0*x+108.9861D0)*x-292.5251D0)*x+
     &             437.1874D0)*x-402.3509D0)*x+238.2397D0)*x-92.2279D0)
     &             *x+23.3323D0)*x-3.820357D0)*x + 1.547924D0
            dsum = dsum*B
            Alfc = Rs*(1.D0+dsum/do)/(2.D0*pi*Z0*B)
         ELSE
            a = 1.D0/(1.D0-T/B)
            s = a**2
            Alfc = Rs*Ep*Z0/(35456.86D0*B)
     &             *(a+2.D0*s*W/B+s*(1.D0+T/B)/pi*dlog((a+1.D0)/(a-1.D0)
     &             ))
         ENDIF
      ENDIF
      Alfd = .5D0*Wn*dsqrt(Ep)*Tep
      END


      SUBROUTINE sm27(W,Rw,T,H,Ep,Epe,Tep,Rs,Z,Wn,Alfc,Alfd)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION a , Alfc , Alfd , Ep , Epe , H , p , q , Rs ,
     &                 Rw , T , Tep , W , Wn , Z
      a = dsqrt(Ep)
      Alfd = .5D0*Wn*a*Tep
      IF ( Ep.NE.1.D0 ) THEN
         q = (Epe-1.D0)/(Ep-1.D0)
         Alfd = a*q*Alfd/dsqrt(Epe)
      ENDIF
C
      Alfc = 0.D0
      IF ( T.EQ.0.D0 ) RETURN
      p = W/H
      IF ( p.LE..159154943D0 ) THEN
         a = dlog(12.5663706D0*W/T-T/W)
      ELSE
         a = dlog(2.D0*H/T-T/H)
      ENDIF
      q = 1.D0 + (1.D0+a/3.14159265D0)/Rw
      Alfc = .159154943D0*q*(1.D0-.0625D0*Rw**2)
      IF ( p.GT.2.D0 ) THEN
         a = .5D0*Rw + .94D0
         q = Rw*(1.D0+1.D0/(a*3.14159265D0))*q
         Alfc = q/(Rw+.63661972D0*dlog(17.0794684D0*a))**2
      ENDIF
      Alfc = Alfc*Rs/(H*Z)
      END

C      SUBROUTINE SM22(W,T,B,EP,Z)
C      A=SIN(1.57079633/B*T)
C      P=.5*ALOG((1.+A)/(1.-A))
C      A=EXP(1.57079633*W/(B-T)+P)
C      AK=2./(A+1./A)
C      AK1=SQRT(1.-AK**2)
C
C      CALL SM10(AK1,E,F)
C      CALL SM10(AK,E1,F1)
C      Z=94.2477795*F/(SQRT(EP)*F1)
C      RETURN
C      END


      SUBROUTINE sm36(Wng,B,B1,Ba,Bb,Bc,Bd,Nol)
C     LOGICAL NOL
C     NOL=.FALSE.
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION a1 , a2 , ak , an , B , B1 , Ba , Bb , Bc , Bd ,
     &                 c , h , pi , w , wlg , Wng , x , z
      INTEGER Nol
      pi = 3.14159265D0
      wlg = 2.D0*pi/Wng
C     IF((B.GT.0.5*WLG).OR.(B1.GT.0.5*WLG)) NOL=.TRUE.
      ak = B1/(2.D0*B)
      c = datan(ak)
      x = ak*ak
      an = 1.D0*x
      h = datan(1.D0/ak)
      z = c/ak
      w = -2.D0*z
      a1 = -(2.D0*ak/pi)*dexp(w)
     &     *(1.D0+(5.D0+x)/(4.D0*(1.D0+x))*dexp(w)+(4.D0/an+
     &     ((5.D0+x)/an)**2)*dexp(2.D0*w)/9.D0)
      a2 = 2.D0*(ak*h+z+dlog(an/(4.D0*ak))-pi*an/(6.D0*ak)) - a1
      Ba = 2.D0*B1/wlg*(h+dlog(dsqrt(an))/ak)
      Bb = .5D0*(Ba-(2.D0*B/wlg*(pi*ak/3.D0+a1)))
      Bc = wlg/(2.D0*pi*B1)
      Bd = B/wlg*(pi/(3.D0*ak)+a2)
      END


      SUBROUTINE sm84(Wn,W,T,H,Ep,Tep,Es,Z,Zs,B,Sig,Dl,Gam)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION ac , ad , alf , B , bet , Dl , Ep , Es , fp , H ,
     &                 rs , Sig , sm28 , T , Tep , W , Wn , Z , zf , Zs
      DOUBLE COMPLEX Gam
C
      CALL sm21(Ep,Es,H,Z,Wn,fp)
C
      CALL sm25(Zs,Z,H,Wn,zf)
      rs = sm28(Wn,1.D0,Sig)
C
      CALL sm27(W,B,T,H,Ep,fp,Tep,rs,zf,Wn,ac,ad)
      alf = (ac+ad)
      bet = Wn*dsqrt(fp)
      bet = bet + .5D0*(ac-ad)**2/bet
      Gam = dcmplx(alf*Dl,bet*Dl)
      END



      SUBROUTINE sm03(Z1,Z2,Z3,Rat,It,Pars)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER It
      DOUBLE PRECISION Rat , t
      DOUBLE COMPLEX Z1 , Z2 , Z3 , a , b , c , d , Pars(3)
      t = dfloat(It)
      a = Z1 + 1.D0
      b = Z3 + Rat
      c = a*Z2 + 1.D0
      d = 1.D0/(a+b*c)
      Pars(2) = 2.D0*d*dsqrt(Rat)
      b = Z3 - Rat
      Pars(3) = t*d*(a+b*c)
      a = Z1 - 1.D0
      b = Z3 + Rat
      c = a*Z2 + 1.D0
      Pars(1) = t*d*(a+b*c)
      END


C      SUBROUTINE SM10(AK1,E,F)
C      T=AK1**2
CC
CC     ELLIPTIC INTEGRAL OF THE SECOND KIND
CC
C      F=((.032024666*T+.054544409)*T+.097932891)*T+1.3862944-
C     *(((.010944912*T+.060118519)*T+.12475074)*T+.5)*ALOG(T)
CC
CC     ELLIPTIC INTEGRAL OF THE FIRST KIND
CC
C      E=((.040905094*T+.085099193)*T+.44479204)*T+1.-
C     *(((.01382999*T+.08150224)*T+.24969795)*T)*ALOG(T)
C      RETURN
C      END




      DOUBLE PRECISION FUNCTION sm28(Wn,Amu,Sigm)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Amu , Sigm , Wn
      sm28 = dsqrt(Wn*Amu*188.49559D0/Sigm)
      END



C      FUNCTION SM63(Z,SEP,H)
CC       CHANGES MADE ON 7.02.91 BY N.K.
CC     LOGICAL*4 SYML
C      P=.25
CC     IF(SYML)P=1.
C      SM63=376.687*P*H/(Z*SEP)
C      RETURN
C      END


      SUBROUTINE sm16(Kc,Ms,I,J,Inc,Out)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER I , Inc , ir , ir1 , J , jr , jr1 , k , Kc , Ms
      LOGICAL*4 Out
      Out = .FALSE.
      k = 0
      ir1 = (I-1)/Kc
      jr1 = (J-1)/Kc
      IF ( ir1.NE.jr1 ) GOTO 20
      IF ( J.GT.Kc ) k = k + jr1
      IF ( k.LE.0 ) THEN
         ir = I
         jr = J
         IF ( Ms.EQ.1 ) GOTO 13
         IF ( Ms.EQ.2 ) GOTO 14
         IF ( Ms.EQ.3 ) GOTO 16
      ENDIF
      ir = mod((I-1),Kc) + 1
      jr = mod((J-1),Kc) + 1
      IF ( Ms.EQ.2 ) GOTO 14
      IF ( Ms.EQ.3 ) GOTO 16
 13   Inc = k*Kc**2 + Kc*(ir-1) + jr
      RETURN
 14   IF ( I.GT.J ) THEN
         jr = mod((I-1),Kc) + 1
         ir = mod((J-1),Kc) + 1
      ENDIF
      Inc = k*(Kc*(Kc+1))/2 + (ir-1)*Kc + jr - ir*(ir-1)/2
      RETURN
 16   IF ( I.EQ.J ) THEN
         Inc = k*Kc + ir
         RETURN
      ENDIF
 20   Out = .TRUE.
      Inc = 0
      END



      SUBROUTINE sm35(W,H,Epe,Z,Wn,R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION cl , Epe , gr , H , sm77 , sm78 , W , Wn , Z
      DOUBLE COMPLEX R , c , coexp
      coexp(c) = dexp(dreal(c))*dcmplx(dcos(dimag(c)),dsin(dimag(c)))
      cl = sm77(W,2.D0*H,Wn)
      R = dcmplx(0.D0,-2.D0*Wn*cl)
      gr = sm78(W,Epe,Wn)*Z
      gr = (1.D0-gr)/(1.D0+gr)
      R = dcmplx(gr,0.D0)*coexp(R)
      END



      DOUBLE PRECISION FUNCTION sm77(W,B,Wn1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION B , c , W , Wn1 , x
      c = .220635D0*B
      x = Wn1*c
      x = 2.D0*(2.D0*c+W)/(c+2.D0*W)*dcos(x)/dsin(x)
      sm77 = (1.570796D0-datan(x))/Wn1
      END



      DOUBLE PRECISION FUNCTION sm78(W,Epe,Bet)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Bet , Epe , P , sm08 , W
      EXTERNAL sm79
      COMMON /gint  / P(2)
      P(1) = Bet
      P(2) = W
      sm78 = dsqrt(Epe)
     &       *(sm08(0.D0,1.5707D0,sm79)+sm08(1.5708D0,3.14159265D0,sm79)
     &       )/2368.79422D0
      END


      DOUBLE PRECISION FUNCTION sm08(A,B,fun)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION A , B , c , d , e , fun , w , x
      INTEGER i
      DIMENSION x(6) , w(6)
      DATA x/0.1252334085D0 , 0.3678314989D0 , 0.5873179542D0 ,
     &     0.7699026741D0 , 0.9041172563D0 , 0.9815606342D0/ ,
     &     w/0.2491470458D0 , 0.2334925365D0 , 0.2031674267D0 ,
     &     0.1600783285D0 , 0.1069393259D0 , 0.0471753363D0/
C
      c = 0
      d = (A+B)/2.D0
      DO i = 1 , 6
         e = (B-A)*x(i)/2.D0
         c = c + w(i)*(fun(d+e)+fun(d-e))
      ENDDO
      sm08 = c*(B-A)/2.D0
      END



      DOUBLE PRECISION FUNCTION sm79(T)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION P , T
      COMMON /gint  / P(2)
      sm79 = (dsin(5.D0*P(1)*P(2)*dcos(T)))**2*(dsin(T))**3/(dcos(T))**2
      END
