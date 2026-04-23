!*==SUPRIM.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE suprim(T,U,Deltu)
C                                                                 *
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , i1 , i1end , i2 , i2e , idu , idubeg , iterm , iu ,
     &        iubeg , jw , jwbeg , k3 , kdun , Kn , Kn1 , Knc , knodes ,
     &        Knr , Knr1
      INTEGER Kol , kop , kou , kunel , kur , kval , l1 , l2 , l3 , le ,
     &        n1 , n2 , n3 , na , nadru , nb , nb1 , nf , nr , nr1
      EQUIVALENCE (k3,Kol(3))
      INCLUDE 'circuit.i'

      COMMON /kolnal/ Kol(4) , Nal(4)
      LOGICAL Nal
      COMMON /blk2  / Knc , Knr , Kn , Knr1 , Kn1
      COMMON /blw1  / W , W1
      DOUBLE PRECISION W(20) , W1(200)
      INTEGER name(4) , noi(5,2) , nou(5,2) , nu1(5,2) , nu2(5,2)
      INTEGER kouv(5) , kopv(5) , koi , nr1v(5) , nb1v(5) , ng , indep
      DOUBLE PRECISION T , t1
      DOUBLE COMPLEX U , Deltu
      DIMENSION U(1) , Deltu(1)
      DOUBLE COMPLEX val(150) , dval(150)
      DOUBLE COMPLEX im/(0.0D0,1.0D0)/ , zero/(0.0D0,0.0D0)/
      DOUBLE COMPLEX u1 , u2 , du1 , du2
      LOGICAL exs(5,2) , exist(5,2)
C     PRINT 26,(U(JJJ),DELTU(JJJ),JJJ=1,2 )
C  26 FORMAT('   U:',2E14.6,'  DELTU:',2E14.6 )
      n1 = Nnetpr
      l1 = Lenntp
      kur = Kn*k3
      T = 1.D0
C  LOOP FOR SEARCHING NONLINEAR ELEMENTS
C  COMPUTATION OF STEP FOR SEARCHING NONLINEAR ELEMENTS.
C  NMPNT - NUMBER OF TYPES OF ELEMENTS.
      i1end = 20*Nmpnt
      DO i1 = 1 , i1end , 20
         IF ( Mpoint(i1+5).EQ.3 ) THEN
C  IF A LINEAR ELEMENT IS FOUND, GO TO LABEL ..100..
            i2e = Mpoint(i1+4)
C  I2E - NUMBER OF ELEMENTS OF THIS TYPE.
            IF ( i2e.NE.0 ) THEN
               nf = Mpoint(i1+9)
               le = Mpoint(i1+6)
               n2 = Mpoint(i1+11)
C  L2-NUMBER OF PARAMETERS OF THIS TYPE OF ELEMENTS.
               l2 = Mpoint(i1+10)
C  DETERMINATION OF THE NAME OF THE ELEMENT TYPE
               name(1) = Mpoint(i1)
               name(2) = Mpoint(i1+1)
               name(3) = Mpoint(i1+2)
               name(4) = Mpoint(i1+3)
               CALL libmd5(name,noi,nou,exist,koi,kouv,kopv,nr1v,nb1v)
C  LOOP OVER ELEMENTS OF THE TYPE
               DO i2 = 1 , i2e
                  l3 = Mpoint(i1+12)
                  knodes = Mpoint(i1+7)
                  na = nf + (i2-1)*le
                  n3 = Nodeel(na+knodes+3)
C  LOOP OVER INDEPENDENT GROUPS OF SOURCES
                  nadru = 1
                  DO indep = 1 , koi
                     ng = indep
                     koi = 1
                     kou = kouv(ng)
                     kop = kopv(ng)
                     nr = kou + kop
                     nb = koi
                     nr1 = nr1v(ng)
                     nb1 = nb1v(ng)
C  DETERMINATION OF THE NODE NUMBER FOR THE APPLICATION OF CONTROL
C  VOLTAGE (NU1) AND CURRENT (NU2)
C  AND TRANSFER OF FLAGS FROM EXIST TO EXS RELATED TO THE CURRENT
C   GROUP
                     DO i = 1 , kou
                        nu1(i,1) = Nodeel(na+nou(nadru,1))
                        nu1(i,2) = Nodeel(na+nou(nadru,2))
                        exs(i,1) = exist(nadru,1)
                        exs(i,2) = exist(nadru,2)
                        nadru = nadru + 1
                     ENDDO
                     nu2(indep,1) = Nodeel(na+noi(indep,1))
                     nu2(indep,2) = Nodeel(na+noi(indep,2))

C  FILLING VAL, DVAL
C        TOTAL: VAL(KN*(KOU+KOP))
C              DVAL(KN*(KOU+KOP))
C
C     MPOINT(I1+17)=KOU+KOP
                     kunel = Kn*kou
                     kdun = Kn*kop
                     kval = nr*Kn
                     DO i = 1 , kval
                        val(i) = zero
                        dval(i) = zero
                     ENDDO
C  FILLING VAL, DVAL
C  FILLED WITH KOU 'BLOCKS' ...
                     idu = 0
                     DO iu = 1 , kou
                        iubeg = (iu-1)*Kn
C  BY KN ELEMENTS
                        DO jw = 1 , Kn
                           jwbeg = (jw-1)*k3
                           u1 = zero
                           u2 = zero
                           du1 = zero
                           du2 = zero
                           IF ( nu1(iu,1).NE.0 ) THEN
                              u1 = U(jwbeg+nu1(iu,1))
                              du1 = Deltu(jwbeg+nu1(iu,1))
                           ENDIF
                           IF ( nu1(iu,2).NE.0 ) THEN
                              u2 = U(jwbeg+nu1(iu,2))
                              du2 = Deltu(jwbeg+nu1(iu,2))
                           ENDIF
                           val(iubeg+jw) = u1 - u2
                           dval(iubeg+jw) = du1 - du2
                        ENDDO
                        IF ( exs(iu,2) ) THEN
                           idu = idu + 1
                           DO jw = 1 , Kn
                              idubeg = (idu-1)*Kn
                              iterm = kunel + idubeg
                              val(iterm+jw) = im*W(jw)*val(iubeg+jw)
                              dval(iterm+jw) = im*W(jw)*dval(iubeg+jw)
                           ENDDO
                        ENDIF
                     ENDDO
C  IN THE LIST OF PARAMETERS LIBMD5, A CHECK IS MADE: 'IS THERE A CONTROL STEP?'
                     t1 = T
                     CALL libmd6(name,ng,n1,l1,n2,l2,n3,l3,val,dval,Kn,
     &                           nr,t1)
                     T = dmin1(T,t1)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
      ENDDO
C  MULTIPLY Y BY T
      DO i = 1 , kur
         Deltu(i) = Deltu(i)*T
      ENDDO
C     DEBUG SUBTRACE,INIT(T,T1)
      END
