!*==SOLVE.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE solve(U)
C  DRIVER FOR THE SOLUTION OF NONLINEAR EQUATIONS SYSTEMS
C          BY THE NEWTON METHOD.
C
C  Ref.:     Dennis J., M.L., Schnabel P.
C       Numerical Methods for Unconstrained Minimization
C       and Solution of Nonlinear Equations: Per. and Engl.
C       M.:Mir, 1988.
C
C         PARAMETERS AND CHANGES:
C  PLEASE UPDATE!
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION cosga , dgamma , epsim , estim , ga , grnor ,
     &                 slope , ynor
      INTEGER i , Iapr , ii , irang1 , irang2 , iretcd , iter , Itermc ,
     &        k0 , k3 , Kn , Kn1 , Knc , Knr , Knr1 , ko , Kprgrf ,
     &        Kprlen , Kprlin , Kprnkr
      INTEGER Kprqup , Kprsol , Kprsrt , Kprvar , kur , Mglob , n , n1
      EQUIVALENCE (k3,Kol(3)) , (kur,k0,ko)
      DOUBLE PRECISION U(1)

      COMMON /ter   / Itermc

      DOUBLE PRECISION f(2*k3*Kn) , y(2*k3*Kn) , gr(2*k3*Kn)
      DOUBLE PRECISION un(2*k3*Kn) , fn(2*k3*Kn)
      DOUBLE PRECISION dfdx(2*k3*Kn,2*k3*Kn+1)
C     IRANG - RANGE OF THE PROBLEM SIZE.
      INTEGER irang
      INTEGER index(2*k3*Kn)
C SF and SX are hafl size - "scalers"
      DOUBLE PRECISION sf(k3*Kn) , sx(k3*Kn)


C  THESE 3 LINES - FOR LINEQ1
      INTEGER ner
      DOUBLE COMPLEX det
      DOUBLE PRECISION argd

      COMMON /kolnal/ Kol , Nal
      INTEGER Kol(4)
      LOGICAL Nal(4)
      COMMON /blk2  / Knc , Knr , Kn , Knr1 , Kn1
      COMMON /modglb/ Mglob , Iapr
      COMMON /print / Kprlen , Kprsrt , Kprnkr , Kprlin , Kprsol ,
     &                Kprvar , Kprgrf , Kprqup
      LOGICAL flag

C  P-P CONNECTION WITH IRANG AND PABEH POLYGON IRANG.
      INTEGER termcd , kmaxdu
      INTEGER flgfft
      DOUBLE PRECISION ts , tls , fnor , fnorn
      LOGICAL maxtkn

C INITIALIZATION

      irang = 2*k3*Kn

      irang2 = irang/2
      irang1 = irang + 1
      n = 2*Kn*k3
      n1 = n + 1
      iter = 0
      flgfft = 0
      ts = 1.D0
      tls = 1.D0
      termcd = 0

C this is temporary - should be value of UIN assigned
C      DO I=1, N
C        UIN(I)=0.
C      ENDDO

      IF ( Kprsol.GE.2 ) WRITE (6,33) n
 33   FORMAT (2X,'ZNEWSOL (÷èïä) :  N=',I4)
      IF ( Kprsol.GE.2 ) WRITE (6,2) (U(i),i=1,n)


C  CALCULATE MACHEPS
      CALL machep(epsim)
C  SETUP AND CHECK OF INPUT INFORMATION
      CALL neinck(n,epsim,irang,sf,sx,U,termcd)
      IF ( termcd.GE.0 ) THEN

C  CALCULATE VECTOR FUNCTION
         CALL nef(irang2,U,fn,sf,fnor,flag,flgfft)
C    FLAG - ERROR INDICATION IN THE MODEL FUNCTION

C  CHECK IF SOLUTION ALREADY EXISTS
         CALL stop0(n,fn,U,sf,sx,termcd,kmaxdu)
         IF ( termcd.EQ.1 ) RETURN
C  CALCULATE JACOBIAN
         CALL nejac(irang2,U,dfdx,flag,flgfft)

C  CALCULATE GRADIENT
C   FOR THE STANDARD NORM OF THE VECTOR ACCORDING TO THE PARAMETER
         CALL gradie(irang,n,dfdx,fn,sf,gr)
C F <- FN
         DO i = 1 , n
            f(i) = fn(i)
         ENDDO
C JACOBIAN NOW ANALYTICAL


C ITERATIONS
         DO WHILE ( termcd.EQ.0 )
            iter = iter + 1
C PRINTING
            IF ( Kprsol.GE.2 ) WRITE (6,5) iter
            PRINT 5 , iter
            IF ( Kprsol.GE.3 ) WRITE (6,2) (U(i),i=1,n)
            IF ( Kprsol.GE.3 ) WRITE (6,1) (f(i),i=1,n)
 1          FORMAT (2X,'  F= ..'/(3X,6(E12.5)))

C SOLUTION OF THE AFFINE MODEL
            CALL regul(irang,dfdx,f,n)
            IF ( Kprsol.GE.3 ) WRITE (6,4900) n , n1
            IF ( Kprsol.GE.3 ) PRINT 4900 , n , n1
            IF ( Kprsol.GE.3 ) WRITE (6,5000) (dfdx(i,n1),i=1,n)
            IF ( Kprsol.GE.3 ) PRINT 5000 , (dfdx(i,n1),i=1,n)
            IF ( Kprsol.GE.3 ) WRITE (6,5005) (f(i),i=1,n)
            IF ( Kprsol.GE.3 ) PRINT 5005 , (f(i),i=1,n)

            DO i = 1 , n
               dfdx(i,n1) = f(i)
            ENDDO
            IF ( Kprsol.GE.3 ) WRITE (6,5100) (dfdx(i,n1),i=1,n)
            IF ( Kprsol.GE.3 ) PRINT 5100 , (dfdx(i,n1),i=1,n)
C 5100 FORMAT(2X,'ZNEWSOL (40)    : DFDX=',3X,E13.6)
C   40 CONTINUE

            CALL lineq1(dfdx,irang,n,irang1,1,index,ner,det,argd)
C  INTEGER INDEX(IRANG),COMPLEX DET,REAL ARGD,INTEGER NER  õöE
C DECLARATIONS. NER=1 - EVERYTHING IS GOOD, NER=0 - ERROR IN COLUMN ALIGNMENT.
            IF ( Kprsol.GE.3 ) WRITE (6,5200) (dfdx(i,n1),i=1,n)
            IF ( Kprsol.GE.3 ) PRINT 5200 , (dfdx(i,n1),i=1,n)

            DO i = 1 , n
               y(i) = dfdx(i,n1)
            ENDDO
            IF ( Kprsol.GE.3 ) WRITE (6,4000) (y(i),i=1,n)
            IF ( Kprsol.GE.3 ) PRINT 4000 , (y(i),i=1,n)

C  EVALUATION OF THE MINIMUM PARTS OF THE TRANSFER AT THE ZERO FREQUENCY
            CALL otchm(y,k3)
            IF ( Kprsol.GE.3 ) WRITE (6,4100) (y(i),i=1,n)
            IF ( Kprsol.GE.3 ) PRINT 4100 , (y(i),i=1,n)

C   CALCULATION OF THE ANGLE BETWEEN GRADIENT AND PROJECTION. ESTIMATION OF THE MODEL JACOBIAN.
C   L2 NORM OF THE GRADIENT AND PROJECTION
            slope = 0.D0
            grnor = 0.D0
            ynor = 0.D0
            DO ii = 1 , n
               slope = slope - gr(ii)*y(ii)
               grnor = grnor + gr(ii)*gr(ii)
               ynor = ynor + y(ii)*y(ii)
            ENDDO
            ynor = dsqrt(ynor)
            grnor = dsqrt(grnor)
C  COS(GAMMA)
            cosga = slope/(grnor*ynor)
            IF ( Kprsol.GE.3 ) WRITE (6,3001) slope , grnor , ynor ,
     &                                cosga
 3001       FORMAT (2X,'ZNEWSOL : SLOPE=',E12.5,' GRNOR=',E12.5,
     &              ' YNOR=',E12.5/12X,'COSGA=',E21.14)
            IF ( cosga.GT.1.000000000D0 ) cosga = 1.000000000D0
            IF ( cosga.LT.-1.000000000D0 ) cosga = -1.000000000D0
            dgamma = dacos(cosga)
            ga = (dgamma/3.14D0)*180.D0
            IF ( Kprsol.GE.2 ) WRITE (6,2003) ga
            PRINT 2003 , ga
            estim = -1.D0/cosga
            IF ( Kprsol.GE.2 ) WRITE (6,2004) estim
            PRINT 2004 , estim

C  GLOBALIZATION:
C   .1. APPROXIMATION
            IF ( Iapr.EQ.1 ) THEN
               CALL suprim(ts,U,y)
               IF ( Kprsol.GE.3 ) WRITE (6,4200) (y(i),i=1,n)
               IF ( Kprsol.GE.3 ) PRINT 4200 , (y(i),i=1,n)
C
C     APPROXIMATION
               DO i = 1 , n
                  un(i) = U(i) - y(i)
               ENDDO
               CALL nef(irang2,un,fn,sf,fnorn,flag,flgfft)
            ENDIF

C .2. NONLINEAR SEARCH.
            IF ( Mglob.EQ.1 ) CALL lserch(irang2,n,U,fnor,gr,y,sx,sf,
     &           iretcd,maxtkn,un,fn,fnorn,tls,flag,flgfft)
            IF ( Kprsol.GE.3 ) WRITE (6,4300) (y(i),i=1,n)
            IF ( Kprsol.GE.3 ) PRINT 4300 , (y(i),i=1,n)

C  IF GLOBALIZATION DID NOT OCCUR, THEN IT IS NECESSARY TO CALCULATE UN, FN, FNORN
            IF ( Mglob.NE.1 .AND. Iapr.NE.1 ) THEN
               DO i = 1 , n
                  un(i) = U(i) - y(i)
               ENDDO
               CALL nef(irang2,U,fn,sf,fnorn,flag,flgfft)
            ENDIF

C  PRINT AGAIN
            IF ( Kprsol.GE.3 ) WRITE (6,3) (y(i),i=1,n)
            IF ( Kprsol.GE.3 ) PRINT 3 , (y(i),i=1,n)
C  STOP THE ITERATION...
C    JACOBIAN:
            CALL nejac(irang2,un,dfdx,flag,flgfft)
C    GRADIENT:
            CALL gradie(irang,n,dfdx,fn,sf,gr)
C  CHECKING THE CONDITION OF THE SYSTEM
            CALL stop(n,U,y,fn,fnorn,gr,sx,sf,iretcd,iter,maxtkn,kmaxdu,
     &                termcd)
            IF ( Kprsol.GE.2 ) WRITE (6,4400) (y(i),i=1,n)
            IF ( Kprsol.GE.3 ) PRINT 4400 , (y(i),i=1,n)


            IF ( Kprsol.GE.3 ) WRITE (6,2) (U(i),i=1,n)


            Itermc = termcd

C  AND PRINT AGAIN
            IF ( Kprsol.GE.3 ) WRITE (6,4) fnor , ts , tls
 4          FORMAT (2X,'FNOR=',E12.5,2X,'TS=',E12.5,2X,'TLS=',E12.5)


C  U <- UN , F <- FN , FNOR <- FNORN
            DO i = 1 , n
               U(i) = un(i)
               f(i) = fn(i)
            ENDDO
            fnor = fnorn

            IF ( Kprsol.GE.2 ) WRITE (6,2) (U(i),i=1,n)
         ENDDO
C___E_N_D__________________________________________________________
         IF ( Kprsol.GE.2 ) WRITE (6,1200)
 1200    FORMAT (2X,'ZNEWSOL (EXIT)  : ')
         IF ( Kprsol.GE.2 ) WRITE (6,2) (U(i),i=1,n)
      ELSE
         IF ( Kprsol.GT.0 ) WRITE (6,1010) termcd
         PRINT 1010 , termcd
         STOP
      ENDIF
 4900 FORMAT (2X,'ZNEWSOL        : N=',I5,'  N1=',I5)
 5000 FORMAT (2X,'ZNEWSOL (REGUL): DFDX='/(3X,6(E13.6)))
 5005 FORMAT (2X,'                 F   ='/(3X,6(E13.6)))
 5100 FORMAT (2X,'ZNEWSOL (40)    : DFDX='/(3X,6(E13.6)))
 5200 FORMAT (2X,'ZNEWSOL (LINEQ1): DFDX='/(3X,6(E13.6)))
 4000 FORMAT (2X,'ZNEWSOL (50)    : Y='/(3X,6(E13.6)))
 4100 FORMAT (2X,'ZNEWSOL (OTCHM) : Y='/(3X,6(E13.6)))
 2003 FORMAT (' *** IN THIS POINT ANGLE BTWN GRAD. AND CORRECTION=',
     &        F8.4)
 2004 FORMAT (' *** LOWER ESTIM. OF CONDITION NUMBER    =',E14.7)
 4200 FORMAT (2X,'ZNEWSOL (SUPRIM): Y='/(3X,6(E13.6)))
 4300 FORMAT (2X,'ZNEWSOL (LSERCH): Y='/(3X,6(E13.6)))
 4400 FORMAT (2X,'ZNEWSOL (STOP): Y='/(3X,6(E13.6)))
 5    FORMAT (2X,'ITERATION ',I4)
 2    FORMAT (2X,'  U= ..'/(3X,6(E12.5)))
 3    FORMAT (2X,' DU= ..'/(3X,6(E12.5)))


C
C
 1010 FORMAT (2X,'  ##### ZNEWSOL:  A T T E N T I O N ! !    #######'/2X
     &        ,'     FATAL ERROR IN INPUT DATA.'/2X,
     &        '            RETURN CODE =',I2)
      END
