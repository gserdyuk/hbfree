!*==MP.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE mp(Om,P1,L1,P2,L2,P3,L3,N)
C
C     MODEL OF MICROSTRIP LINES
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION a , alfa , alfd , alfp , b , bett , cc , d1 ,
     &                 d11 , d2 , d22 , dgr , e0 , eer , eps , er , f ,
     &                 fl , fr , g
      DOUBLE PRECISION h , olam , pi , rm , ro , rs , skl1 , skl2 , t ,
     &                 tgd , w , weh , weh1 , yo , yrr , zexp , zo , zop
      INTEGER i , j , L1 , L2 , L3 , N
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION Om , P1 , P2 , P3

      COMMON /subc  / Yy(15,15) , Vj(15)
      DOUBLE COMPLEX g1 , g2 , ch1 , sh1 , y11 , y12 , y13
      DOUBLE COMPLEX Yy , Vj , cth1 , csch , yr
      DOUBLE PRECISION m0
C
C
      eps = P1(1)
      tgd = P1(2)
      h = P2(1)
      t = P2(2)
      w = P3(1)
      fl = P3(2)
      eer = P3(3)
      zo = P3(4)
      alfa = P3(5)
      skl1 = P3(6)
      skl2 = P3(7)
C
      er = eps
      e0 = 8.85418D-12
      pi = 3.1415926D0
      cc = 2.99D+8
      m0 = 1.256637D-06
      rm = 5.8D-07
C
      IF ( Om.EQ.0.0D0 ) THEN
C
C+++++++++++++++++++++++++++++++++++++++++++
C  CALCULATION OF THE Y-MATRIX FOR DC
C+++++++++++++++++++++++++++++++++++++++++++
         DO i = 1 , 4
            DO j = 1 , 4
               Yy(i,j) = dcmplx(0.0D0,0.0D0)
            ENDDO
         ENDDO
C
         yrr = (1.D0/rm)*fl
         yr = dcmplx(yrr,0.0D0)
         Yy(1,1) = yr
         Yy(2,2) = Yy(1,1)
         Yy(1,2) = -Yy(1,1)
         Yy(2,1) = Yy(1,2)
      ELSE
         f = Om/(2.D0*pi)
C
C
         IF ( eer.LE.0.0D0 ) CALL serf(w,h,t,eer,weh1,er)
         IF ( zo.LE.0.0D0 ) CALL sz0(w,h,t,eer,weh1,zo)
C
         IF ( skl1.EQ.1.D0 ) THEN
C
C               CALCULATION OF FREQUENCY DISPERSION
C  DISPERSION OF WAVE IMPEDANCE AND DIELECTRIC CONSTANT
C
            f = Om/6.2831852D0
            g = dsqrt((zo-5.D0)/60.D0) + 0.004D0*zo
            fr = 0.3976D0*zo/(h*1.D+03)
            g1 = 1.D0 + g*(f*1.D-09/fr)**2.D0
            CALL z00p(w,h,t,er,zop)
            zo = zop - (zop-zo)/g1
            eer = er - (er-eer)/g1
         ENDIF
C
         IF ( skl2.EQ.1.D0 ) THEN
            WRITE (6,100) eer , zo , f
 100        FORMAT (8X,'EER=',F5.2,4X,'Zo=',F6.2,' OHM',4X,'F=',E11.4,
     &              ' HZ')
         ENDIF
C
C
         bett = Om*dsqrt(eer)/cc
         olam = 6.2831852D0/bett
C
         dgr = 360.D0*fl/olam
         IF ( skl2.EQ.1.D0 ) THEN
            WRITE (6,102) dgr , f
 102        FORMAT (8X,'PHASE SHIFT OF THE LINE=',F6.2,' DEG',4X,'F=',
     &              E11.4,' HZ')
         ENDIF
C

         IF ( alfa.LE.0.D0 ) THEN
C
C      CALCULATION OF LOSSES
C  LOSSES IN A MICROSTRIP LINE (Np/m)
C

            IF ( t.EQ.0.D0 ) t = 0.1D-08
            ro = 5.8D-07
            f = Om/6.2831852D0
            rs = dsqrt(3.1415926D0*f*m0*ro)
            olam = 3.D+08/f
            IF ( (w/h).GE.0.15915494D0 ) THEN
               b = h
               weh = w/h + 0.3978873D0*t/h*(1.D0+dlog(2.D0*h/t))
            ELSE
               b = 6.2831852D0*w
               weh = w/h + 0.3978873D0*t/h*(1.D0+dlog(12.56637D0*w/t))
            ENDIF
            a = 1.D0 + (1.D0+0.31830989D0*dlog(2.D0*b/t))/weh
            IF ( (w/h).LE.1.D0 ) THEN
               alfp = 1.38D0*a*rs*(32.D0-weh*weh)/(h*zo*(32.D0+weh*weh))
            ELSE
               alfp = 6.1D-05*a*rs*zo*eer/h*
     &                (weh+0.667D0*weh/(weh+1.444D0))
            ENDIF
            alfd = 27.3D0*er*(eer-1.D0)*tgd/(er-1.D0)/dsqrt(eer)/olam
            alfa = (alfp+alfd)/8.686D0
         ENDIF
C
         IF ( skl2.EQ.1.D0 ) THEN
            WRITE (6,101) alfa , f
 101        FORMAT (8X,'ALFA=',E11.4,' Np/m',4X,'F=',E11.4,' HZ')
         ENDIF
C
C
         yo = 1.D0/zo
         d2 = bett
         d1 = alfa
         g1 = dcmplx(d1,d2)
         d11 = -d1
         d22 = -d2
         g2 = dcmplx(d11,d22)
         ch1 = 0.5D0*(zexp(g1*fl)+zexp(g2*fl))
         sh1 = 0.5D0*(zexp(g1*fl)-zexp(g2*fl))
         cth1 = ch1/sh1
         csch = 1.D0/sh1
C
         y11 = yo*cth1
         y12 = yo*csch
         y13 = dcmplx(0.0D0,0.0D0)
C
         Yy(1,1) = y11
         Yy(2,2) = y11
         Yy(1,2) = -y12
         Yy(2,1) = -y12
         Yy(1,3) = y13
         Yy(1,4) = y13
         Yy(2,3) = y13
         Yy(2,4) = y13
         Yy(3,1) = y13
         Yy(3,2) = y13
         Yy(3,3) = y13
         Yy(3,4) = y13
         Yy(4,1) = y13
         Yy(4,2) = y13
         Yy(4,3) = y13
C
         Yy(4,4) = y13
      ENDIF
C
      IF ( skl2.EQ.1.D0 ) THEN
         WRITE (6,103) f , Yy(1,1) , Yy(1,2) , Yy(1,3) , Yy(1,4) ,
     &                 Yy(2,1) , Yy(2,2) , Yy(2,3) , Yy(2,4) , Yy(3,1) ,
     &                 Yy(3,2) , Yy(3,3) , Yy(3,4) , Yy(4,1) , Yy(4,2) ,
     &                 Yy(4,3) , Yy(4,4)
 103     FORMAT (24X,'F=',E11.4,' HZ',/(E11.4,E11.4),'i',1X,
     &           (E11.4,E11.4),'i',1X,(E11.4,E11.4),'i',1X,(E11.4,E11.4)
     &           ,'i',/(E11.4,E11.4),'i',1X,(E11.4,E11.4),'i',1X,
     &           (E11.4,E11.4),'i',1X,(E11.4,E11.4),'i',/(E11.4,E11.4),
     &           'i',1X,(E11.4,E11.4),'i',1X,(E11.4,E11.4),'i',1X,
     &           (E11.4,E11.4),'i',/(E11.4,E11.4),'i',1X,(E11.4,E11.4),
     &           'i',1X,(E11.4,E11.4),'i',1X,(E11.4,E11.4),'i')
      ENDIF
      END

      SUBROUTINE z00p(W,H,T,Er,Zop)
C
C  CALCULATION OF THE WAVE IMPEDANCE OF A STRIPLINE
C  FOR THE DISPERSION OF THE DIELECTRIC CONSTANT AND
C  WAVE IMPEDANCE OF THE STRIPLINE.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION a1 , b , Er , H , sk , sk1 , skk , T , W , w1bt ,
     &                 w33 , x , xm , Zop
      b = 2.D0*H
      IF ( (T/H).LE.0.00005D0 ) THEN
         sk = dtanh(1.5707963D0*W/b)
         sk1 = dsqrt(1.D0-sk*sk)
         IF ( sk.GE.0.0D0 .AND. sk.LE.0.7D0 ) THEN
            skk = 3.183098D-01*dlog(2.D0*(1.D0+dsqrt(sk1))
     &            /(1.D0-dsqrt(sk1)))
         ELSE
            skk = 1.D0/
     &            (3.1830989D-01*dlog(2.D0*(1.D0+dsqrt(sk))/(1.D0-dsqrt
     &            (sk))))
         ENDIF
         Zop = 94.24778D0*skk/dsqrt(Er)
      ELSE
         x = T/b
         xm = 2.D0/(1.D0+0.666666D0*x/(1.D0-x))
         a1 = x/(3.1415926D0*(1.D0-x))
     &        *(1.D0-0.5D0*dlog((x/(2.D0-x))**2.D0+
     &        (0.0796D0*x/(W/b+1.1D0*x))**xm))
         w33 = W/(b-T)
         w1bt = w33 + a1
         Zop = 30.D0/dsqrt(Er)
     &         *dlog(1.D0+1.2732395D0/w1bt*(2.546479D0/w1bt+
     &         dsqrt((2.546479D0/w1bt)**2.D0+6.27D0)))
      ENDIF
      END

      SUBROUTINE smpl(Om,P1,L1,P2,L2,P3,L3,N)
C
C     MODEL OF COUPLED MICROSTRIP LINES
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION alfe , alfo , bett , bette , betto , cc , ce ,
     &                 ceb , cend , co , cob , dfle , dflo , e0 , eer ,
     &                 eps , er , eree , ereo , f
      DOUBLE PRECISION fl , fl4 , h , pi , rd , rm , s , skl1 , skl2 ,
     &                 t , tgd , w , weh1 , wh , z0 , z01 , zoe , zoo
      INTEGER i , j , L1 , L2 , L3 , N
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION Om , P1 , P2 , P3

      COMMON /subc  / Yy(15,15) , Vj(15)
      DOUBLE COMPLEX JJ , qe , qo , yr , sie , sio , coe , coo , tae ,
     &               tao
      DOUBLE COMPLEX Yy , Vj
      DOUBLE COMPLEX yy11 , yy12 , yy13 , yy14
C
      DOUBLE PRECISION m0
      PARAMETER (JJ=(0.D0,1.D0))

C
C
      eps = P1(1)
      tgd = P1(2)
      rd = P1(4)
      h = P2(1)
      t = P2(2)
      w = P3(1)
      s = P3(2)
      fl = P3(3)
      eree = P3(4)
      ereo = P3(5)
      zoe = P3(6)
      zoo = P3(7)
      alfe = P3(8)
      alfo = P3(9)
      skl1 = P3(10)
      skl2 = P3(11)
C
      er = eps
      e0 = 8.85418D-12
      pi = 3.1415926D0
      cc = 2.99D+8
      m0 = 1.256637D-06
      rm = 5.8D-07
      f = Om/(2.D0*pi)
C
      IF ( Om.EQ.0.0D0 ) THEN
C
         DO i = 1 , 4
            DO j = 1 , 4
               Yy(i,j) = dcmplx(0.0D0,0.0D0)
            ENDDO
         ENDDO
         yr = JJ*(1.D0/rm)*fl
         Yy(1,1) = yr
         Yy(2,2) = Yy(1,1)
         Yy(3,3) = Yy(1,1)
         Yy(4,4) = Yy(1,1)
         Yy(1,3) = -Yy(1,1)
         Yy(2,4) = Yy(1,3)
         Yy(3,1) = Yy(1,3)
         Yy(4,2) = Yy(1,3)
      ELSE
C
         IF ( P3(4).LE.0.0D0 .AND. P3(5).LE.0.0D0 .AND.
     &        P3(6).LE.0.0D0 .AND. P3(7).LE.0.D0 ) THEN
            CALL serf(w,h,t,eer,weh1,er)
            CALL sz0(w,h,t,eer,weh1,z01)
            CALL coce(e0,er,w,h,cc,s,t,eer,z01,co,ce)
            CALL serf(w,h,t,eer,weh1,1.D0)
            CALL sz0(w,h,t,eer,weh1,z01)
            CALL coce(e0,1.D0,w,h,cc,s,t,eer,z01,cob,ceb)
C
            eree = ce/ceb
            ereo = co/cob
C
            zoe = 1.D0/(cc*dsqrt(ce*ceb))
            zoo = 1.D0/(cc*dsqrt(co*cob))
C
            IF ( skl1.EQ.1.D0 ) CALL dispz(Om,zoe,zoo,er,eree,ereo,h,w,
     &           s,pi)
         ENDIF
         z0 = dsqrt(zoe*zoo)
         IF ( P3(6).LE.0.0D0 .AND. P3(7).LE.0.0D0 ) THEN
            wh = w/h
            IF ( wh.LE.0.37D0 .AND. w.LE.0.3D0 .AND. s.LE.0.078D0 ) THEN
               zoo = zoo - 4.5D0
               zoe = zoe - 20.D0
            ELSE
               zoo = zoo - 4.5D0
               zoe = zoe - 4.5D0
            ENDIF
         ENDIF
C
         IF ( skl2.EQ.1.0D0 ) THEN
            WRITE (6,1020) eree , ereo , zoe , zoo , z0 , f
 1020       FORMAT (/,4X,'EREE=',F5.2,2X,'EREO=',F5.2,4X,'ZOE=',F6.2,
     &              'OHM',' ZOO=',F6.2,'OHM',4X,'Z0=',F6.2,'OHM   F=',
     &              E11.4,' HZ',/)
         ENDIF
C
         CALL bet1(Om,cc,eree,ereo,bett,bette,betto,pi)
C
         IF ( skl2.EQ.1.0D0 ) THEN
            fl4 = 2.D0*pi/bett/4.D0
            WRITE (6,1520) fl4 , f
 1520       FORMAT (/,4X,'WAVELENGTH LAMBDA/4.=',E11.4,
     &              '(m)   FREQUENCY F=',E11.4,' HZ',/)
         ENDIF
C
C
         IF ( alfe.LE.0.0D0 .OR. alfo.LE.0.0D0 )
     &        CALL alf1(w,rm,zoe,zoo,pi,Om,er,eree,ereo,tgd,cc,alfe,
     &        alfo)
C
C
         IF ( skl2.EQ.1.0D0 ) THEN
            WRITE (6,107) alfe , alfo , f
 107        FORMAT (/,8X,'ALFE=',E11.4,' Np/m',4X,'ALFO=',E11.4,' Np/m',
     &              4X,'F=',E11.4,' HZ',/)
         ENDIF
C
         IF ( skl1.GT.1.01D0 ) THEN
            cend = w*er/(120.D0*pi*cc)
     &             *(1.35D0/dlog10(4.D0*h/t)+w*1.333313D0/
     &             (3.0D0*h*dsqrt(er)))
            dfle = 0.92D0*datan(Om*cend*zoe)/bette
            dflo = 0.95D0*datan(Om*cend*zoo)/betto
            IF ( skl2.EQ.1.0D0 ) THEN
               WRITE (6,4151) cend , dfle , dflo
 4151          FORMAT (/,1X,
     &      '(FOR OPEN CIRCUIT AT END OF THE LINE) EDGE CAPACITANCE CK='
     &      ,E11.4,/10X,'INCREASE IN LENGTH DFLE=',E11.4,1X,
     &      'INCREASE IN LENGTH DFLO=',E11.4)
            ENDIF
         ELSE
            dfle = 0.0D0
            dflo = 0.0D0
         ENDIF
C
         qe = bette*(fl+dfle) - JJ*alfe*(fl+dfle)
         qo = betto*(fl+dflo) - JJ*alfo*(fl+dflo)
         sie = sin(qe)
         sio = sin(qo)
         coe = cos(qe)
         coo = cos(qo)
         tae = sie/coe
         tao = sio/coo
C
         yy11 = -JJ/2.D0*((1.D0/tae)/zoe+(1.D0/tao)/zoo)
         yy14 = JJ/2.D0*((1.D0/sin(qe))/zoe+(1.D0/sin(qo))/zoo)
         yy13 = JJ/2.D0*((1.D0/sin(qe))/zoe-(1.D0/sin(qo))/zoo)
         yy12 = -JJ/2.D0*((1.D0/tae)/zoe-(1.D0/tao)/zoo)
C
         Yy(1,1) = yy11
         Yy(1,2) = yy12
         Yy(1,3) = yy14
         Yy(1,4) = yy13
         Yy(2,1) = yy12
         Yy(2,2) = yy11
         Yy(2,3) = yy13
         Yy(2,4) = yy14
         Yy(3,1) = yy14
         Yy(3,2) = yy13
         Yy(3,3) = yy11
         Yy(3,4) = yy12
         Yy(4,1) = yy13
         Yy(4,2) = yy14
         Yy(4,3) = yy12
C
         Yy(4,4) = yy11
      ENDIF
C
      IF ( skl2.EQ.1.0D0 ) THEN
         WRITE (6,103) f , Yy(1,1) , Yy(1,2) , Yy(1,3) , Yy(1,4) ,
     &                 Yy(2,1) , Yy(2,2) , Yy(2,3) , Yy(2,4) , Yy(3,1) ,
     &                 Yy(3,2) , Yy(3,3) , Yy(3,4) , Yy(4,1) , Yy(4,2) ,
     &                 Yy(4,3) , Yy(4,4)
 103     FORMAT (24X,'F=',E11.4,' HZ',/(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i',2X,(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i',/(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i',2X,(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i',/(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i',2X,(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i',/(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i',2X,(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i')
      ENDIF
C
      END

      SUBROUTINE bet1(Om,Cc,Eree,Ereo,Bett,Bette,Betto,Pi)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Bett , Bette , Betto , Cc , Eree , Ereo , olam ,
     &                 Om , Pi
      olam = Cc*2.D0*Pi/Om
      Bette = 2.D0*Pi*dsqrt(Eree)/olam
      Betto = 2.D0*Pi*dsqrt(Ereo)/olam
      Bett = 0.5D0*(Bette+Betto)
C
      END

      SUBROUTINE alf1(W,Rm,Zoe,Zoo,Pi,Om,Er,Eree,Ereo,Tgd,Cc,Alfe,Alfo)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Alfe , alfec , Alfo , alfoc , Cc , Er , Eree ,
     &                 Ereo , frq , Om , Pi , Rm , Tgd , W , Zoe , Zoo
      frq = Om/(Pi*2.D0)
      Alfe = 27.3D0*Er*(Eree-1.D0)*Tgd/((Er-1.D0)*dsqrt(Eree)*Cc/frq)
      Alfo = 27.3D0*Er*(Ereo-1.D0)*Tgd/((Er-1.D0)*dsqrt(Ereo)*Cc/frq)
      Alfe = Alfe/8.686D0
      Alfo = Alfo/8.686D0
      alfec = Pi*dsqrt(frq*Rm)/(W*Zoo)
      alfoc = Pi*dsqrt(frq*Rm)/(W*Zoe)
      Alfe = (Alfe+alfec)*1.D-04 + 0.25D0
      Alfo = (Alfo+alfoc)*1.D-04 - 0.035D0
      END

      SUBROUTINE coce(E0,Er,W1,H1,Cc,S1,T1,Eer,Z01,Co,Ce)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION a , Cc , Ce , cf1 , cfe , cfo , cga , cgd , Co ,
     &                 cpe , cpo , dth , dwh , E0 , Eer , Er , H1 , S1 ,
     &                 sk , skk1
      DOUBLE PRECISION T1 , t2 , W1 , wh1 , wh2 , Z01
      t2 = T1/2.D0
      wh1 = W1/H1
      wh2 = wh1
      IF ( S1.GE.t2 ) THEN
         dth = T1/(Er*S1)
         IF ( wh1.LE.0.159155D0 ) THEN
            dwh = 0.397887D0*T1*(1.D0+dlog(12.5663D0*W1/T1))/H1
         ELSE
            dwh = 0.397887D0*T1*(1.D0+dlog(2.D0*H1/T1))/H1
         ENDIF
         wh1 = W1/H1 + dwh*(1.D0-0.5D0*dexp(-0.69D0*dwh/dth))
         wh2 = wh1 + dth
      ENDIF
      cpe = E0*Er*wh1
      cpo = E0*Er*wh2
      cfe = 0.5D0*dsqrt(Eer)/(Cc*Z01) - 0.5D0*cpe
      cfo = 0.5D0*dsqrt(Eer)/(Cc*Z01) - 0.5D0*cpo
      a = dexp(-0.1D0*dexp(2.33D0-2.53D0*wh1))
      cf1 = cfe*dsqrt(Er/Eer)/(1.D0+a*(H1/S1)*dtanh(10.D0*S1/H1))
      Ce = cpe + cfe + cf1
      sk = (S1/H1)/(S1/H1+2.D0*wh2)
      CALL sk1k(sk,skk1)
      cga = E0*skk1
      cgd = E0*Er/3.14115926D0*dlog(1.D0/dtanh(S1*0.7854D0/H1))
     &      + 0.65D0*cfo*(0.02D0*dsqrt(Er)*H1/S1+1.D0-1.D0/(Er*Er))
      Co = cpo + cfo + cga + cgd
      END

      SUBROUTINE sk1k(Sk,Skk1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Sk , sk1 , skk , Skk1
      sk1 = dsqrt(1.D0-Sk)
      IF ( Sk.LE.0.7D0 ) THEN
         skk = dlog(2.D0*(1.D0+dsqrt(sk1))/(1.D0-dsqrt(sk1)))
     &         /3.1415926D0
      ELSE
         skk = 3.1415926D0/dlog(2.D0*(1.D0+dsqrt(Sk))/(1.D0-dsqrt(Sk)))
      ENDIF
      Skk1 = skk
      END

      SUBROUTINE dispz(Om,Zoe,Zoo,Er,Eree,Ereo,H,W,S,Pi)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Er , Eree , Ereo , fr1 , fre , fro , gre , gro ,
     &                 H , Om , Pi , pk1e , pk1o , pke , pko , S , W ,
     &                 Zoe , zoes , Zoo
      DOUBLE PRECISION zoos
      DOUBLE PRECISION kpke , kpko
      fr1 = Om/(2.D0*Pi)*1.D-09
      fro = 0.7952D0*Zoo/(H*1.D+03)
      fre = 0.1988D0*Zoe/(H*1.D+03)
      gro = 0.6D0 + 0.018D0*Zoo
      gre = 0.6D0 + 0.0045D0*Zoe
      Eree = Er - (Er-Eree)/(1.D0+(fr1/fre)**2.D0*gre)
      Ereo = Er - (Er-Ereo)/(1.D0+(fr1/fro)**2.D0*gro)
      pke = dtanh(Pi*0.25D0*W/H)*dtanh(Pi*0.25D0*(W+S)/H)
      pk1e = dsqrt(1.D0-pke*pke)
      pko = (dtanh(Pi*0.25D0*W/H))/(dtanh(Pi*0.25D0*(W+S)/H))
      pk1o = dsqrt(1.D0-pko*pko)
      IF ( pke.LE.0.7D0 ) THEN
         kpke = 1.D0/Pi*dlog(2.D0*(1.D0+dsqrt(pk1e))/(1.D0-dsqrt(pk1e)))
      ELSE
         kpke = 1.D0/
     &          (1.D0/Pi*dlog(2.D0*(1.D0+dsqrt(pke))/(1.D0-dsqrt(pke))))
      ENDIF
      IF ( pko.LE.07 ) THEN
         kpko = 1.D0/Pi*dlog(2.D0*(1.D0+dsqrt(pk1o))/(1.D0-dsqrt(pk1o)))
      ELSE
         kpko = 1.D0/
     &          (1.D0/Pi*dlog(2.D0*(1.D0+dsqrt(pko))/(1.D0-dsqrt(pko))))
      ENDIF
      zoes = 2.D0*30.D0*Pi*kpke/dsqrt(Er)
      zoos = 2.D0*30.D0*Pi*kpko/dsqrt(Er)
      Zoe = zoes - (zoes-Zoe)/(1.D0+(fr1/fre)**1.6D0*gre)
      Zoo = zoos - (zoos-Zoo)/(1.D0+(fr1/fro)**1.6D0*gro)
      END

      SUBROUTINE serf(W1,H1,T1,Eer,Weh1,Er)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION dwh , Eer , Er , fwh , H1 , q , T1 , th , W1 ,
     &                 Weh1 , wh
      wh = W1/H1
      th = T1/H1
      IF ( th.EQ.0.0D0 ) THEN
         Eer = (Er+1.D0)*0.5D0 + 0.5D0*(Er-1.D0)/dsqrt(1.D0+10.D0/wh)
      ELSE
         IF ( wh.LE.0.159155D0 ) THEN
            dwh = 0.39788736D0*th*(1.D0+dlog(12.56637D0*W1/T1))
         ELSE
            dwh = 0.39788736D0*th*(1.D0+dlog(2.D0/th))
         ENDIF
         Weh1 = wh + dwh
         fwh = 1.D0/dsqrt(1.D0+10.D0/wh)
         q = (Er-1.D0)*T1/(H1*4.6D0*dsqrt(wh))
         Eer = (Er+1.D0)*0.5D0 + 0.5D0*(Er-1.D0)*fwh - q
      ENDIF
      END

      SUBROUTINE sz0(W1,H1,T1,Eer,Weh1,Z01)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Eer , H1 , T1 , th , W1 , Weh1 , wh , Z01
      wh = W1/H1
      th = T1/H1
      IF ( th.LE.1.0D0 ) THEN
         IF ( wh.LE.1.0D0 ) THEN
            Z01 = 60.D0*dlog(8.D0/wh+0.25D0*wh)/dsqrt(Eer)
         ELSE
            Z01 = 376.7D0/(dsqrt(Eer)
     &            *(wh+1.393D0+0.667D0*dlog(wh+1.444D0)))
         ENDIF
      ELSEIF ( wh.LE.1.0D0 ) THEN
         Z01 = 60.D0*dlog(8.D0/Weh1+0.25D0*Weh1)/dsqrt(Eer)
      ELSE
         Z01 = 376.7D0/(dsqrt(Eer)
     &         *(Weh1+1.393D0+0.667D0*dlog(Weh1+1.444D0)))
      ENDIF
      END

      SUBROUTINE lang(Om,P1,L1,P2,L2,P3,L3,N)
C
C     MODEL OF COUPLED MICROSTRIP LINES
C         ( CONNECTION OF LINES )
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION alfe , alfo , b01 , b02 , b1 , b11 , b2 , b22 ,
     &                 b3 , b4 , b44 , b5 , b55 , bett , bette , betto ,
     &                 cc , ce , ceb , co
      DOUBLE PRECISION cob , e0 , eer , enn , eps , er , eree , ereo ,
     &                 f , fl , fl4 , h , pi , rd , rm , s , skl1 ,
     &                 skl2 , t , tgd
      DOUBLE PRECISION w , weh1 , wh , z0 , z01 , zoe , zoo
      INTEGER i , j , L1 , L2 , L3 , N
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION Om , P1 , P2 , P3

      COMMON /subc  / Yy(15,15) , Vj(15)
      DOUBLE COMPLEX JJ , qe , qo , yr , sie , sio , coe , coo , tae ,
     &               tao
      DOUBLE COMPLEX Yy , Vj
      DOUBLE COMPLEX yy11 , yy12 , yy13 , yy14
C
      DOUBLE PRECISION m0
      PARAMETER (JJ=(0.D0,1.D0))

C
C
      eps = P1(1)
      tgd = P1(2)
      rd = P1(4)
      h = P2(1)
      t = P2(2)
      w = P3(1)
      s = P3(2)
      fl = P3(3)
      eree = P3(4)
      ereo = P3(5)
      zoe = P3(6)
      zoo = P3(7)
      alfe = P3(8)
      alfo = P3(9)
      enn = P3(10)
      skl1 = P3(11)
      skl2 = P3(12)
C
      er = eps
      e0 = 8.85418D-12
      pi = 3.1415926D0
      cc = 2.99D+8
      m0 = 1.256637D-06
      rm = 5.8D-07
      f = Om/(2.D0*pi)
      wh = w/h
C
      IF ( Om.EQ.0.0D0 ) THEN
C
         DO i = 1 , 4
            DO j = 1 , 4
               Yy(i,j) = dcmplx(0.0D0,0.0D0)
            ENDDO
         ENDDO
         yr = JJ*(1.D0/rm)*fl
         Yy(1,1) = yr
         Yy(2,2) = Yy(1,1)
         Yy(3,3) = Yy(1,1)
         Yy(4,4) = Yy(1,1)
         Yy(1,4) = -Yy(1,1)
         Yy(2,3) = Yy(1,4)
         Yy(4,1) = Yy(1,4)
         Yy(3,2) = Yy(1,4)
      ELSE
C
         IF ( P3(4).LE.0.0D0 .AND. P3(5).LE.0.0D0 .AND.
     &        P3(6).LE.0.0D0 .AND. P3(7).LE.0.D0 ) THEN
            CALL serf(w,h,t,eer,weh1,er)
            CALL sz0(w,h,t,eer,weh1,z01)
            CALL coce(e0,er,w,h,cc,s,t,eer,z01,co,ce)
            CALL serf(w,h,t,eer,weh1,1.D0)
            CALL sz0(w,h,t,eer,weh1,z01)
            CALL coce(e0,1.D0,w,h,cc,s,t,eer,z01,cob,ceb)
            eree = ce/ceb
            ereo = co/cob
C
            b1 = enn*wh
            b11 = enn*wh
            b2 = wh/(3.D0*dsqrt(er))
            b22 = wh/(3.D0*dsqrt(er))
            b3 = 1.35D0/(dlog10(4.D0*h/t))
            b4 = (2.D0*enn-1.D0)*w*s/(3.D0*h*(3.4D0*w+s)*dsqrt(er))
            b44 = 4.D0*(2.D0*enn-1.D0)
     &            *3.4D0*w/(3.D0*(s+3.4D0*w)*dsqrt(er))
            b5 = 1.35D0*(2.D0*enn-1.D0)*s/(dlog10(4.D0*h/t)*(3.4D0*w+s))
            b55 = 1.35D0*(2.D0*enn-1.D0)
     &            /dlog10(4.D0*s*dtanh(4.D0*h/s)/pi/t)
            b01 = 376.9D0/dsqrt(er)
            b02 = 378.D0/dsqrt(er)
C
            zoe = b01/(b1+b2+b3+b4+b5)
            zoe = zoe*1.0039843D0
C
            zoo = b02/(b11+b22+b3+b44+b55)
            zoo = zoo*1.0440309D0
            IF ( skl1.NE.1.D0 ) CALL dispz(Om,zoe,zoo,er,eree,ereo,h,w,
     &           s,pi)
         ENDIF
         z0 = dsqrt(zoe*zoo)
C
         IF ( skl2.EQ.1.D0 ) THEN
            WRITE (6,1020) eree , ereo , zoe , zoo , z0 , f
 1020       FORMAT (4X,'EREE=',F5.2,2X,'EREO=',F5.2,4X,'ZOE=',F6.2,
     &              'OHM',' ZOO=',F6.2,'OHM',4X,'Z0=',F6.2,'OHM   F=',
     &              E11.4,' HZ')
         ENDIF
C
C
         CALL bet1(Om,cc,eree,ereo,bett,bette,betto,pi)
C
         IF ( skl2.EQ.1.D0 ) THEN
            fl4 = 2.D0*pi/bett/4.D0
            WRITE (6,1520) fl4 , f
 1520       FORMAT (4X,'WAVELENGTH LAMBDA/4.=',E11.4,
     &              '(M)   FREQUENCY F=',E11.4,' HZ')
         ENDIF
C
C
         IF ( alfe.LE.0.D0 .OR. alfo.LE.0.D0 )
     &        CALL alf1(w,rm,zoe,zoo,pi,Om,er,eree,ereo,tgd,cc,alfe,
     &        alfo)
C
C
         IF ( skl2.EQ.1.D0 ) THEN
            WRITE (6,107) alfe , alfo , f
 107        FORMAT (8X,'ALFE=',E11.4,' Np/m',4X,'ALFO=',E11.4,' Np/m',
     &              4X,'F=',E11.4,' HZ')
         ENDIF
C
         qe = bette*fl - JJ*alfe*fl
         qo = betto*fl - JJ*alfo*fl
         sie = sin(qe)
         sio = sin(qo)
         coe = cos(qe)
         coo = cos(qo)
         tae = sie/coe
         tao = sio/coo
C
C
         yy11 = -JJ/2.D0*((1.D0/tae)/zoe+(1.D0/tao)/zoo)
         yy14 = JJ/2.D0*((1.D0/sin(qe))/zoe+(1.D0/sin(qo))/zoo)
         yy13 = JJ/2.D0*((1.D0/sin(qe))/zoe-(1.D0/sin(qo))/zoo)
         yy12 = -JJ/2.D0*((1.D0/tae)/zoe-(1.D0/tao)/zoo)
C

C
         Yy(1,1) = yy11
         Yy(1,2) = yy12
         Yy(1,4) = yy14
         Yy(1,3) = yy13
C
         Yy(2,1) = yy12
         Yy(2,2) = yy11
         Yy(2,4) = yy13
         Yy(2,3) = yy14
C
         Yy(4,1) = yy14
         Yy(4,2) = yy13
         Yy(4,4) = yy11
         Yy(4,3) = yy12
C
         Yy(3,1) = yy13
         Yy(3,2) = yy14
         Yy(3,4) = yy12
C
         Yy(3,3) = yy11
      ENDIF
C
      IF ( skl2.EQ.1.D0 ) THEN
         WRITE (6,103) f , Yy(1,1) , Yy(1,2) , Yy(1,3) , Yy(1,4) ,
     &                 Yy(2,1) , Yy(2,2) , Yy(2,3) , Yy(2,4) , Yy(3,1) ,
     &                 Yy(3,2) , Yy(3,3) , Yy(3,4) , Yy(4,1) , Yy(4,2) ,
     &                 Yy(4,3) , Yy(4,4)
 103     FORMAT (24X,'F=',E11.4,' HZ',/(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i',2X,(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i',/(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i',2X,(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i',/(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i',2X,(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i',/(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i',2X,(E11.4,1X,E11.4),'i',2X,
     &           (E11.4,1X,E11.4),'i')
      ENDIF
C
      END
