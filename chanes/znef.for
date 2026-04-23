!*==NEF.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE nef(Ntot,U,F,Sf,Fnor,Flag,Flgfft)
C******************************************************************
C                                                                 *
C     Step for calculating the i-th iteration of the vector of constraints (F) *
C                                                                 *
C                                                                 *
C******************************************************************
C                                                                 *
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , i1 , i1end , i2 , i2e , idu , idubeg , if , ifind1 ,
     &        ifind2 , ind1 , ind2 , ind3 , iterm , iu , iubeg , iw ,
     &        j , jw , jwbeg
      INTEGER k , k3 , kdun , ke , Kn , Kn1 , Knc , knc2 , kni1 ,
     &        knodes , Knr , Knr1 , Kol , kop , kou , kst , kunel ,
     &        kur , kw , kznn
      INTEGER l1 , l2 , l3 , le , m , m1 , mf , n0 , n1 , n2 , n3 , na ,
     &        nadru , nb , nb1 , nf , nr , nr1 , Ntot , nu
      INTEGER nui21 , nui22
      INCLUDE 'circuit.i'
      INCLUDE 'funcsize.i'
      EQUIVALENCE (k3,Kol(3))
      COMMON /kolnal/ Kol(4) , Nal(4)
C  Modification from 30.01.91    Modified by Serdyuk G.V.     See MAIN
      COMMON /maty  / Buffer(6000) , Buflen
      DOUBLE COMPLEX Buffer
      INTEGER*4 Buflen , Flgfft

      LOGICAL Nal
      COMMON /blk1  / Kr , Kr1 , Kc , Kc1 , Nnr , Nnr1 , Mn , Mn1
      COMMON /blk2  / Knc , Knr , Kn , Knr1 , Kn1
      COMMON /blw1  / W , W1
      DOUBLE PRECISION W(20) , W1(200)
      INTEGER Kr(20) , Kc(10) , Nnr(10) , Mn(2,20)
      INTEGER Kr1(200) , Kc1(20) , Nnr1(20) , Mn1(2,200)
      INTEGER name(4) , noi(5,2) , nou(5,2) , nu1(5,2) , nu2(5,2)
      INTEGER kouv(5) , kopv(5) , koi , nr1v(5) , nb1v(5) , ng , indep
      DOUBLE PRECISION Sf(1) , Fnor
      DOUBLE COMPLEX F(1) , U(1)
      DOUBLE COMPLEX unel(MAXU_MD*MAXKN) , dundt(MAXD_MD*MAXKN)
      DOUBLE COMPLEX b1(B1_SIZE) , b2(B2_SIZE)
      DOUBLE COMPLEX znn((MAXU_MD+MAXD_MD)*MAXKN)
C..............Dimensions - Specify!!!

C     COMPLEX VJ,YR   **** Old. 30.01.91  Modified by Serdyuk G.V.  See MAIN

      DOUBLE COMPLEX im/(0.0D0,1.0D0)/ , zero/(0.0D0,0.0D0)/
      DOUBLE COMPLEX u1 , u2
      LOGICAL Flag , exs(5,2) , exist(5,2)

C  Modification from 30.01.91    Modified by Serdyuk G.V.     See MAIN
C   Built-in functions for calculating addresses in BUFFER
      ifind1(i,m,nu) = i + (m-1)*nu
      ifind2(i,j,m,nu,mf) = nu*mf + i + (j-1)*nu + (m-1)*nu*nu

C..............Warning! MAX Dimension F(600=20*30).....
C  Initialization of the vector

      kur = Kn*k3
      DO j = 1 , kur
         F(j) = zero
      ENDDO
      n1 = Nnetpr
      l1 = Lenntp
C  Number of iterations = Number of frequencies + Number of nodes
C  If iterations exceed the size of the vector F, jump to label ..600.

      IF ( kur.GT.Ntot ) THEN

C  SMALL ERROR F
         PRINT 48 , k3 , Kn , Ntot
         WRITE (6,48) k3 , Kn , Ntot
         STOP
      ELSE

C  If the number of non-zero rows is greater than 10, jump to label ..500.

         IF ( Knr.LE.10 ) THEN
C  Define the path for 2-K iterations in KNC. KNC is the size of the
C  transformation matrix.
            k = 8
            DO i = 3 , 7
               IF ( Knc.LE.k ) GOTO 25
               k = k + k
            ENDDO
         ENDIF
         GOTO 500
C............This is the embedded /by-default/ size B1...
 25      Knc = k
         m1 = i
         knc2 = Knc + Knc
C  LOOP TO SEARCH FOR NEW ELEMENTS
C  CALCULATING STEP FOR SEARCHING NEW ELEMENTS.
C  NMPNT - NUMBER OF TYPES OF ELEMENTS.
         i1end = 20*Nmpnt
         DO i1 = 1 , i1end , 20
            IF ( Mpoint(i1+5).EQ.3 ) THEN
C  IF A NEW ELEMENT IS FOUND, JUMP TO LABEL ..100..
               i2e = Mpoint(i1+4)
C  I2E - NUMBER OF ELEMENTS OF THIS TYPE.
               IF ( i2e.NE.0 ) THEN
                  nf = Mpoint(i1+9)
                  le = Mpoint(i1+6)
                  n2 = Mpoint(i1+11)
C  L2 - NUMBER OF PARAMETERS OF THIS TYPE OF ELEMENT.
                  l2 = Mpoint(i1+10)
C  DEFINING THE TYPE IDENTIFIER OF THE ELEMENT.
                  name(1) = Mpoint(i1)
                  name(2) = Mpoint(i1+1)
                  name(3) = Mpoint(i1+2)
                  name(4) = Mpoint(i1+3)
                  CALL libmd5(name,noi,nou,exist,koi,kouv,kopv,nr1v,
     &                        nb1v)
C  LOOP OVER ELEMENT TYPES
                  DO i2 = 1 , i2e
                     l3 = Mpoint(i1+12)
C     WRITE(6,1002) I1,L3
C1002 FORMAT(2X,'ZNEF: MPOINT(',I3,'+12)=',I4,' (=L3)')
C  N3=MPOINT(I1+13)+(I2-1)*L3 -THIS CANNOT BE DONE!!! IT MUST BE DONE THIS WAY:
                     knodes = Mpoint(i1+7)
                     na = nf + (i2-1)*le
                     n3 = Nodeel(na+knodes+3)
C     WRITE(6,1003) NA,KNODES,N3
C1003 FORMAT(2X,'      NODEEL(',I3,'+',I3,'+3)=N3=',I4)
C  LOOP OVER INDEPENDENT (BETWEEN EACH OTHER) SOURCES OF CURRENT
                     nadru = 1
                     DO indep = 1 , koi
                        ng = indep

                        kou = kouv(ng)
                        kop = kopv(ng)
                        nr = kou + kop
                        nb = 1
                        nr1 = nr1v(ng)
                        nb1 = nb1v(ng)
C  DETERMINATION OF THE NODE OF THE APPLICATION CONTROLLING CURRENT
C  (NU1) AND CURRENT (NU2)
C  AND TRANSFER OF COMMANDS FROM EXIST TO EXS, ADJUSTING TO THE CURRENT
C  GROUP
                        DO i = 1 , kou
C     WRITE(6,1004) NADRU,NOU(NADRU,1),NA,NODEEL(NA+NOU(NADRU,1))
C1004 FORMAT(2X,'ZNEF: NOU(',I3,',1)=',I5,' NA=',I3,' NU1(I,1)=',I3)

                           nu1(i,1) = Nodeel(na+nou(nadru,1))
C     WRITE(6,1005) NADRU,NOU(NADRU,2),NA,NODEEL(NA+NOU(NADRU,2))
C1005 FORMAT(2X,'ZNEF: NOU(',I3,',2)=',I5,' NA=',I3,' NU1(I,2)=',I3)
                           nu1(i,2) = Nodeel(na+nou(nadru,2))
                           exs(i,1) = exist(nadru,1)
                           exs(i,2) = exist(nadru,2)
                           nadru = nadru + 1
                        ENDDO
                        nu2(1,1) = Nodeel(na+noi(indep,1))
                        nu2(1,2) = Nodeel(na+noi(indep,2))

C  FILLING ZNN, UNEL, DUNDT
C          ALL: UNEL(KN*KOU)
C               DUNDT(KN*KOP)
C               ZNN(KN*(KOU+KOP))
C     MPOINT(I1+17)=KOU+KOP
C
                        kunel = Kn*kou
                        kdun = Kn*kop
                        kznn = (kou+kop)*Kn
                        DO i = 1 , kunel
                           unel(i) = zero
                        ENDDO
                        IF ( kdun.NE.0 ) THEN
                           DO i = 1 , kdun
                              dundt(i) = zero
                           ENDDO
                        ENDIF
C  FILLING UNEL, DUNDT, ZNN
C  KOU 'BLOCKOB' IS FILLED ...
                        idu = 0
                        DO iu = 1 , kou
                           iubeg = (iu-1)*Kn
C  FOR KN ELEMENTS
                           DO jw = 1 , Kn
                              jwbeg = (jw-1)*k3
                              u1 = zero
                              u2 = zero
C     WRITE(6, 41) IU,NU1(IU,1),JWBEG
CC!  *             NU1(IU,1),U(JWBEG+NU1(IU,1))

C  41 FORMAT(2X,'ZNEF: NU1(',I3,',1)=',I5,'  JWBEG=',I5)
CC!  * 'U(',I3,'+',I3,')=',E12.5,2X,E12.5)

                              IF ( nu1(iu,1).NE.0 )
     &                             u1 = U(jwbeg+nu1(iu,1))
                              IF ( nu1(iu,2).NE.0 )
     &                             u2 = U(jwbeg+nu1(iu,2))
                              unel(iubeg+jw) = u1 - u2
                              znn(iubeg+jw) = unel(iubeg+jw)
                           ENDDO
                           IF ( exs(iu,2) ) THEN
                              idu = idu + 1
                              DO jw = 1 , Kn
                                 idubeg = (idu-1)*Kn
                                 dundt(idubeg+jw) = im*W(jw)
     &                              *unel(iubeg+jw)
                                 iterm = kunel + idubeg
                                 znn(iterm+jw) = dundt(idubeg+jw)
                              ENDDO
                           ENDIF
                        ENDDO
C  FOURIER TRANSFORMATION - COMPUTATION OF IHE
                        if = 0
C     WRITE(6,1000) N1,L1,N2,L2,N3,L3
                        CALL ftmas2(znn,Kr,Kc,Nnr,Knr,Knc,Kn,nr,nb,m1,
     &                              b1,b2,if,Flgfft,*50)

C1000 FORMAT(2X,'ZNEF: N1=',I3,' L1=',I3,' N2=',I3,' L2=',I3,
C    *         ' N3=',I5,' L3=',I3)

 50                     CALL libmd3(name,ng,n1,l1,n2,l2,n3,l3,b1,knc2,
     &                              nr,*300)

C     PRINT 1232,(B1(II),II=1,KNC)
                        CALL ft2(znn,Kr,Kc,Nnr,Knr,Knc,Kn,nr,nb,b1,b2,
     &                           if,Flgfft,*50)
C !       ðPé ðOCTAHOBKE VFFT-BCTABéTø M1         $
C     WRITE(6, 1232) (II,ZNN(II),II=1,KZNN)
C  FORMATION OF THE BOUNDARY VECTOR - ENTRY OF ZNN INTO F
                        i = 1
                        nui21 = nu2(i,1)
                        nui22 = nu2(i,2)
                        kni1 = Kn*(i-1)
                        DO j = 1 , Kn
                           n0 = (j-1)*k3
                           ind1 = n0 + nui21
                           ind2 = n0 + nui22
                           ind3 = j + kni1
                           IF ( nui21.NE.0 ) F(ind1) = F(ind1)
     &                          + znn(ind3)
                           IF ( nui22.NE.0 ) F(ind2) = F(ind2)
     &                          - znn(ind3)
                        ENDDO
C     PRINT 6543,IM
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDDO

C  ENTRY OF THE LINEAR PART OF THE VECTORS AND THE INDEPENDENT SOURCES
C  INTO THE BOUNDARY VECTOR.
         DO iw = 1 , Kn
            kw = (iw-1)*k3
            DO kst = 1 , k3
               F(kw+kst) = F(kw+kst) - Buffer(ifind1(kst,iw,k3))
               DO ke = 1 , k3
C     WRITE(6,3111) KW,KST,F(KW+KST),KE,KST,IW,
C    *		BUFFER( IFIND2(KE,KST,IW,K3,KN) ),KE,KW,
C    *              U(KE+KW)
C3111 FORMAT(2X,'F(',I3,'+',I3,')=',E12.5,2X,E12.5/
C    *       2X,'YR(',3(I5),')=',E12.5,2X,E12.5/
C    *       2X,'U(',I3,'+',I3,')=',E12.5,2X,E12.5)
                  F(kw+kst) = F(kw+kst)
     &                        + Buffer(ifind2(ke,kst,iw,k3,Kn))*U(ke+kw)
               ENDDO
            ENDDO
         ENDDO

C  CALCULATE FNOR=1/2*FVEC_TRANSP*FVEC, (FVEC=F)
         Fnor = 0.0D0
         DO i = 1 , kur
C     WRITE(6,155) FNOR,I,SF(I),F(I)
C 155 FORMAT(2X,'ZNEF : FNOR,I,SF(I),F(I) = ',E12.5,I5,2X,E12.5,
C    *       5X,E12.5,2X,E12.5)
            Fnor = Fnor + Sf(i)*Sf(i)*F(i)*dconjg(F(i))
         ENDDO
         Fnor = Fnor/2.D0

C  F        FOURIER TRANSFORMATION
         Flag = .TRUE.
         RETURN
      ENDIF
C     DEBUG SUBTRACE
C  EXIT INTO THE RESTRICTED AREA
 300  Flag = .FALSE.
      STOP

C  ERROR KN AND/OR KNC
 500  PRINT 49 , Knr , Knc
      WRITE (6,49) Knr , Knc
      STOP
C
C***********************************************************************
 49   FORMAT (/3X,'INPUT DATA ERROR:'/3X,'KNR=',I5,5X,'KNC=',I5)

 48   FORMAT (//9X,' CIRCUIT CAN NOT BE ANALYSED:'/3X,
     &        ' NUMBER OF EQUATIONS EXCEEDS MAXIMUM.'/3X,
     &        ' NUMBER OF NODES=',I5/3X,' NUMBER OF FREQUENCIES=',I5/3X,
     &        ' MAX SIZE=',I5)
C1232 FORMAT(2X,'ZNN(',I3,')=',E13.6,2X,E13.6)
      END
