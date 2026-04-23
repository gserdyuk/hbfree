!*==NEJAC.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE nejac(Ntot,U,Dfdx,Flag,Flgfft)
C******************************************************************
C                                                                 *
C     SUBROUTINE COMPUTES THE JACOBIAN DFDX                       *
C     AT THE I-TH ITERATION                                       *
C                                                                 *
C******************************************************************
C                                                                 *
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , i1 , i160 , i1end , i2 , i2e , ib1 , ib2 , ibl ,
     &        idiag , idr , ids , idu , idubeg , if , ifind2 , indzn1 ,
     &        ir , irb , is
      INTEGER is1 , is2 , isb , itr , its , iu , iubeg , j , j1kn ,
     &        j1kn1 , jj , jw , jwbeg , k , k3 , kdun , kki , kku , Kn ,
     &        Kn1
      INTEGER Knc , knc2 , knodes , Knr , Knr1 , Kol , kop , kou ,
     &        kunel , kunel1 , kur , kur1 , kznn1 , l1 , l2 , l3 , le ,
     &        ls , m , m1
      INTEGER mf , n0 , n001 , n002 , n01 , n02 , n1 , n2 , n3 , na ,
     &        nadru , nb , nb1 , nf , nr , nr1 , Ntot , nu
      EQUIVALENCE (k3,Kol(3))
      INCLUDE 'circuit.i'
      INCLUDE 'funcsize.i'
      COMMON /kolnal/ Kol(4) , Nal(4)

      COMMON /maty  / Buffer(6000) , Buflen
      DOUBLE COMPLEX Buffer
      INTEGER*4 Buflen
      INTEGER Flgfft

      LOGICAL Nal
      COMMON /blk1  / Kr , Kr1 , Kc , Kc1 , Nnr , Nnr1 , Mn , Mn1
      COMMON /blk2  / Knc , Knr , Kn , Knr1 , Kn1
      COMMON /blw1  / W , W1
      COMMON /blw2  / Wr , Ws
      DOUBLE PRECISION W(20) , W1(200)
      INTEGER Kr(20) , Kc(10) , Nnr(10) , Mn(2,20)
      INTEGER Kr1(200) , Kc1(20) , Nnr1(20) , Mn1(2,200)
      INTEGER Wr(20,20) , Ws(20,20)
C  DATA TYPE OF NAME AND PARTIALLY MPOINT MUST BE CHARACTER*4
      INTEGER name(4) , noi(5,2) , nou(5,2) , nu1(5,2) , nu2(5,2)
      INTEGER kouv(5) , kopv(5) , koi , nr1v(5) , nb1v(5) , ng , indep
      DOUBLE COMPLEX U(1) , Dfdx(Ntot,2,1)
      DOUBLE COMPLEX b1(B1_SIZE) , b2(B2_SIZE)
      DOUBLE COMPLEX unel(MAXU_MD*MAXKN1) , dundt(MAXD_MD*MAXKN1)
      DOUBLE COMPLEX znn1((MAXU_MD+MAXD_MD)*MAXKN1)

C ..............DIMENSIONS - TO BE SPECIFIED!!!.....
      DOUBLE COMPLEX ew , g , gi
      DOUBLE COMPLEX sum1 , sum2 , diff1 , diff2
      DOUBLE COMPLEX im/(0.0D0,1.0D0)/ , zero/(0.0D0,0.0D0)/
      DOUBLE COMPLEX u1 , u2
      LOGICAL Flag , exs(5,2) , exist(5,2)
C ..............ATTENTION! MAX DIMENSION F(600=20*30).....

C  CHANGE FROM 30.01.91    MODIFIED BY SERDYUK G.V.   SEE MAIN
C     BUILT-IN FUNCTIONS FOR INDEXING:
c      IFIND1(I,M,NU)=(I+(M-1)*NU)
      ifind2(i,j,m,nu,mf) = nu*mf + i + (j-1)*nu + (m-1)*nu*nu


      n1 = Nnetpr
      l1 = Lenntp
C  NUMBER OF EQUATIONS = NUMBER OF FREQUENCIES + NUMBER OF NODES
      kur = Kn*k3
C  ZEROING THE JACOBIAN (DFDX)
      kur1 = kur + 1
      DO i = 1 , kur
         DO j = 1 , kur
            Dfdx(j,1,i) = zero
            Dfdx(j,2,i) = zero
         ENDDO
      ENDDO
C  IF THE NUMBER OF EQUATIONS IS GREATER THAN THE SYSTEM DIMENSION
C  GO TO LABEL ..600..
      IF ( kur.GT.Ntot ) THEN

C  SMALL DIMENSION OF DFDX
         WRITE (6,48) k3 , Kn , Ntot
         PRINT 48 , k3 , Kn , Ntot
         STOP
      ELSE
C  IF THE NUMBER OF NON-ZERO ROWS IS GREATER THAN 10,
C           GO TO LABEL ..500..
         IF ( Knr.LE.10 ) THEN

C  DETERMINATION OF MULTIPLICITY 2-KE KNC.
C  KNC - DIMENSION OF FOURIER TRANSFORMATION.
            k = 8
            DO i = 3 , 7
               IF ( Knc.LE.k ) GOTO 25
               k = k + k
            ENDDO
         ENDIF
         GOTO 500
C ............THIS IS DUE /IN MY OPINION/ TO THE DIMENSION OF B1...
 25      Knc = k
         m1 = i
         knc2 = Knc + Knc
C  LOOP FOR SEARCHING NONLINEAR ELEMENTS
C  CALCULATION OF STEP SIZE FOR SEARCHING NONLINEAR ELEMENTS.
C  NMPNT - NUMBER OF TYPES OF ELEMENTS.
         i1end = 20*Nmpnt
         DO i1 = 1 , i1end , 20
            IF ( Mpoint(i1+5).EQ.3 ) THEN
C  IF A NONLINEAR ELEMENT IS FOUND, GO TO LABEL ..100..
               i2e = Mpoint(i1+4)
C  I2E - NUMBER OF ELEMENTS OF THE GIVEN TYPE.
               IF ( i2e.NE.0 ) THEN
                  nf = Mpoint(i1+9)
                  le = Mpoint(i1+6)
                  n2 = Mpoint(i1+11)
C  L2 - NUMBER OF PARAMETERS OF THE GIVEN TYPE OF ELEMENTS.
                  l2 = Mpoint(i1+10)
C  DETERMINATION OF THE DESIGNATION OF THE ELEMENT TYPE
                  name(1) = Mpoint(i1)
                  name(2) = Mpoint(i1+1)
                  name(3) = Mpoint(i1+2)
                  name(4) = Mpoint(i1+3)
                  CALL libmd5(name,noi,nou,exist,koi,kouv,kopv,nr1v,
     &                        nb1v)
C  CYCLE OVER ELEMENTS OF THE TYPE
                  DO i2 = 1 , i2e
                     l3 = Mpoint(i1+12)
C  N3=MPOINT(I1+13)+(I2-1)*L3 - THIS CANNOT BE DONE!!! IT MUST BE DONE LIKE THIS:
                     knodes = Mpoint(i1+7)
                     na = nf + (i2-1)*le
                     n3 = Nodeel(na+knodes+3)
C  CYCLE OVER INDEPENDENT (INTERRELATED) CURRENT SOURCES
                     nadru = 1
                     DO indep = 1 , koi
                        ng = indep

                        kou = kouv(ng)
                        kop = kopv(ng)
                        nr = kou + kop
                        nb = 1
                        nr1 = nr1v(ng)
                        nb1 = nb1v(ng)
C  DETERMINATION OF THE NODE OF THE APPLICATION CONTROLLING
C  THE STRESS (NU1) AND THE POINT (NU2)
C  AND TRANSFER OF FLAGS FROM EXIST TO EXS, CORRESPONDING TO THE CURRENT
C  GROUP
                        DO i = 1 , kou
                           nu1(i,1) = Nodeel(na+nou(nadru,1))
                           nu1(i,2) = Nodeel(na+nou(nadru,2))
                           exs(i,1) = exist(nadru,1)
                           exs(i,2) = exist(nadru,2)
                           nadru = nadru + 1
                        ENDDO
                        nu2(1,1) = Nodeel(na+noi(indep,1))
                        nu2(1,2) = Nodeel(na+noi(indep,2))

C     PRINT 6543,IM
C  FILLING OF ZNN, UNEL, DUNDT
C       ENTIRE: UNEL(KN*KOU)
C               DUNDT(KN*KOP)
C               ZNN(KN*(KOU+KOP))
C     MPOINT(I1+17)=KOU+KOP
C
                        kunel = Kn*kou
                        kdun = Kn*kop
C     KZNN=(KOU+KOP)*KN
                        DO i = 1 , kunel
                           unel(i) = zero
                        ENDDO
                        IF ( kdun.NE.0 ) THEN
                           DO i = 1 , kdun
                              dundt(i) = zero
                           ENDDO
                        ENDIF
C  FILLING OF UNEL, DUNDT, ZNN
C  KOU 'BLOCKS' ARE FILLED ...
                        idu = 0
                        DO iu = 1 , kou
                           iubeg = (iu-1)*Kn
C  рO KN ьм-TOB
                           DO jw = 1 , Kn
                              jwbeg = (jw-1)*k3
                              u1 = zero
                              u2 = zero
                              IF ( nu1(iu,1).NE.0 )
     &                             u1 = U(jwbeg+nu1(iu,1))
                              IF ( nu1(iu,2).NE.0 )
     &                             u2 = U(jwbeg+nu1(iu,2))
                              unel(iubeg+jw) = u1 - u2
                           ENDDO
                           IF ( exs(iu,2) ) THEN
                              idu = idu + 1
                              DO jw = 1 , Kn
                                 idubeg = (idu-1)*Kn
                                 dundt(idubeg+jw) = im*W(jw)
     &                              *unel(iubeg+jw)
                              ENDDO
                           ENDIF
                        ENDDO

C  CALCULATION OF UNEL AND DUNDT ON THE DESIRED GRID
C  !! THE VALUE NR1 ASSUMES THE ROLE OF KOU+KOP
                        kznn1 = nr1*Kn1
                        DO j = 1 , kznn1
                           znn1(j) = zero
                        ENDDO
C  LENGTH OF THE ZNN1 SECTION, WHERE THE UNEL STRESSES ARE CONTAINED
                        kunel1 = Kn1*kou
C  COMPLETE REFERENCE OF UNEL AND DUNDT IN ZNN1
                        jj = 0
                        DO j = 1 , kou
                           IF ( exs(j,2) ) jj = jj + 1
                           j1kn1 = (j-1)*Kn1
                           j1kn = (j-1)*Kn
                           DO i = 1 , Kn
                              ls = Ws(i,1)
                              znn1(j1kn1+ls) = unel(j1kn+i)
                              IF ( exs(j,2) ) THEN
                                 indzn1 = kunel1 + ls + (jj-1)*Kn1
                                 znn1(indzn1) = dundt(i+j1kn)
                              ENDIF
                           ENDDO
                        ENDDO

C     TRANSFORMATION OF THE FUNCTION - CALCULATION OF DI/DU

                        if = 0
                        CALL ftmas2(znn1,Kr1,Kc1,Nnr1,Knr1,Knc,Kn1,nr1,
     &                              nb1,m1,b1,b2,if,Flgfft,*70)
 70                     CALL libmd4(name,ng,n1,l1,n2,l2,n3,l3,b1,knc2,
     &                              nr1,*300)
                        CALL ft2(znn1,Kr1,Kc1,Nnr1,Knr1,Knc,Kn1,nr1,nb1,
     &                           b1,b2,if,Flgfft,*70)
C !    SEE WARNINGS IN THE FIRST LIST

C  FILLING OF DF/DX
                        n0 = Kn1*kou
C  $-  LENGTH OF THE ZNN1 SECTION, CONTAINING DI/DU
C  THE NEXT PART CONTAINS DI/D(DU/DT)
                        n01 = 0
                        n02 = 0
                        kki = 1
C     INDEP, NG     - SELECTED VARIABLES.
C
                        DO kku = 1 , kou
                           n01 = n01 + 1
                           IF ( exs(kku,2) ) n02 = n02 + 1
C THAT IS, IF A BLOCK OF KN1 ELEMENTS PASSES - PROCEED TO THE NEXT.
                           isb = 0
                           DO is = 1 , Kn
                              ew = im*W(is)
                              irb = 0
                              DO ir = 1 , Kn
                                 sum1 = zero
                                 sum2 = zero
                                 diff1 = zero
                                 diff2 = zero

                                 n001 = (n01-1)*Kn1 + 1
                                 CALL sumdif(znn1(n001),ir,is,sum1,
     &                              diff1)
                                 IF ( exs(kku,2) ) THEN
                                    n002 = n0 + (n02-1)*Kn1 + 1
                                    CALL sumdif(znn1(n002),ir,is,sum2,
     &                                 diff2)
                                 ENDIF
C      WRITE(6, 120) SUM1,EW,DIFF2,G,
C     *              IM,DIFF1,IS,W(IS),SUM2,GI
C  120 FORMAT(2X,'ZNEJAC: SUM1=',E13.6,',',E13.6,' EW=',E13.6,',',E13.6/
C     *       2X,'        DIFF2=',E13.6,',',E13.6,' G=',E13.6,',',E13.6/
C     *       2X,'       IM=',E13.6,',',E13.6,' DIFF1=',E13.6,',',E13.6/
C     *       2X,'       W(',I5,')=',E13.6,' SUM2=',E13.6,',',E13.6/
C     *       2X,'       GI=',E13.6,',',E13.6)
                                 g = sum1 + ew*diff2
                                 gi = im*diff1 - W(is)*sum2
                                 ib1 = irb + nu2(kki,1)
                                 ib2 = irb + nu2(kki,2)
                                 is1 = isb + nu1(kku,1)
                                 is2 = isb + nu1(kku,2)
C
C  FILLING OF Y-MATRIX. IF THE NODE IS ALREADY CONNECTED - NONZERO.
                                 IF ( nu1(kku,1).NE.0 .AND. nu2(kki,1)
     &                                .NE.0 ) Dfdx(ib1,1,is1)
     &                                = Dfdx(ib1,1,is1) + g
                                 IF ( nu1(kku,2).NE.0 .AND. nu2(kki,2)
     &                                .NE.0 ) Dfdx(ib2,1,is2)
     &                                = Dfdx(ib2,1,is2) + g
                                 IF ( nu1(kku,2).NE.0 .AND. nu2(kki,1)
     &                                .NE.0 ) Dfdx(ib1,1,is2)
     &                                = Dfdx(ib1,1,is2) - g
                                 IF ( nu1(kku,1).NE.0 .AND. nu2(kki,2)
     &                                .NE.0 ) Dfdx(ib2,1,is1)
     &                                = Dfdx(ib2,1,is1) - g

                                 IF ( nu1(kku,1).NE.0 .AND. nu2(kki,1)
     &                                .NE.0 ) Dfdx(ib1,2,is1)
     &                                = Dfdx(ib1,2,is1) + gi
                                 IF ( nu1(kku,2).NE.0 .AND. nu2(kki,2)
     &                                .NE.0 ) Dfdx(ib2,2,is2)
     &                                = Dfdx(ib2,2,is2) + gi
                                 IF ( nu1(kku,2).NE.0 .AND. nu2(kki,1)
     &                                .NE.0 ) Dfdx(ib1,2,is2)
     &                                = Dfdx(ib1,2,is2) - gi
                                 IF ( nu1(kku,1).NE.0 .AND. nu2(kki,2)
     &                                .NE.0 ) Dfdx(ib2,2,is1)
     &                                = Dfdx(ib2,2,is1) - gi

                                 irb = irb + k3
                              ENDDO
                              isb = isb + k3
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDDO

C  DESIGNATION OF Y-MATRIX
         idiag = 0
         DO ibl = 1 , Kn
            DO itr = 1 , k3
               idr = idiag + itr
               DO its = 1 , k3
                  ids = idiag + its

C  CHANGE FROM 30.01.91    MODIFIED BY SERDYUK G.V.
                  Dfdx(ids,1,idr) = Dfdx(ids,1,idr)
     &                              + Buffer(ifind2(its,itr,ibl,k3,Kn))
                  IF ( idiag.NE.0 ) Dfdx(ids,2,idr) = Dfdx(ids,2,idr)
     &                 + im*Buffer(ifind2(its,itr,ibl,k3,Kn))

               ENDDO
            ENDDO
            idiag = idiag + k3
         ENDDO

         DO i160 = 1 , k3
            Dfdx(i160,2,i160) = im
         ENDDO

C  F & DFDX FORMATION
         Flag = .TRUE.
         RETURN
      ENDIF
C     DEBUG SUBTRACE,INIT(I160,IB1,IB2,IS1,IS2,IDS,IDR)
C  EXIT TO A RESTRICTED AREA
 300  Flag = .FALSE.
      STOP

C  ERROR IN KN AND/OR KNC
 500  WRITE (6,49) Knr , Knc
      PRINT 49 , Knr , Knc
      STOP
C
C***********************************************************************
 49   FORMAT (/3X,'INPUT DATA ERROR:'/3X,'KNR=',I5,5X,'KNC=',I5)

 48   FORMAT (//9X,' CIRCUIT CAN NOT BE ANALYSED:'/3X,
     &        ' NUMBER OF EQUATIONS EXCEEDS MAXIMUM.'/3X,
     &        ' NUMBER OF NODES=',I5/3X,' NUMBER OF FREQUENCIES=',I5/3X,
     &        ' MAX SIZE=',I5)
C     DEBUG SUBTRACE,INIT(IB1,IB2,IS1,IS2,IRB,ISB,NU1,NU2)
C     DEBUG SUBTRACE,SUBCHK
      END
