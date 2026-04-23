!*==LSERCH.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE lserch(Irang2,N,X,Fnor,G,Y,Sx,Sf,Iretcd,Maxtkn,Xn,Fn,
     &                  Fnor1,Lambda,Flag,Flgfft)
C
C  LINE SEARCH.   MODIFIED GOLDSTEIN ALGORITHM - ARMIJO
C
C$LARGE:Y
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION a , alpha , b , div , Fnor , Fnor1 , fprev1 ,
     &                 rel , rellen , sinv , slope , v1 , v2
      INTEGER i , Irang2 , Kprgrf , Kprlen , Kprlin , Kprnkr , Kprqup ,
     &        Kprsol , Kprsrt , Kprvar , Limit , N
      DOUBLE PRECISION X(1) , G(1) , Y(1) , Sx(1) , Sf(1) , Xn(1) ,
     &                 Fn(1)
      INTEGER Iretcd
      LOGICAL Maxtkn
      COMMON /newton/ Epssol , Epsdu , Epsmin , Maxdu , Limit
      DOUBLE PRECISION Maxdu , Epssol , Epsdu , Epsmin
      COMMON /print / Kprlen , Kprsrt , Kprnkr , Kprlin , Kprsol ,
     &                Kprvar , Kprgrf , Kprqup
      DOUBLE PRECISION newtln , minlbd , Lambda , ltemp , lprev
C   THIS IS FOR NEF
      INTEGER Flgfft
      LOGICAL Flag

C INITIALIZATION OF THE ALGORITHM
      Maxtkn = .FALSE.
      Iretcd = 2
      alpha = 1.D-4

C  CALCULATION OF THE L-2 NORM OF VECTOR DU
      newtln = 0.D0
      DO i = 1 , N
         newtln = newtln + (Y(i)*Sx((i+1)/2))**2
      ENDDO
      newtln = dsqrt(newtln)
      IF ( newtln.GT.Maxdu ) THEN

C NEWTON STEP IS GREATER THAN MAXDU
         rel = Maxdu/newtln
         DO i = 1 , N
            Y(i) = rel*Y(i)
         ENDDO
         newtln = Maxdu
      ENDIF

C  CALCULATION OF DECREASE RATE
      slope = 0.D0
      DO i = 1 , N
         slope = slope - G(i)*Y(i)
      ENDDO

C  RELATIVE STEP LENGTH
      rellen = 0.D0
      DO i = 1 , N
         sinv = 1/Sx((i+1)/2)
         rellen = dmax1(rellen,dabs(Y(i))/dmax1(dabs(X(i)),sinv))
      ENDDO

C  MINIMALLY ALLOWABLE STEP LENGTH
      minlbd = Epsdu/rellen
C INITIALIZATION LAMBDA
      Lambda = 1.0D0
      DO

C_ COMPUTE LAMBDA_____________________________________________________
         DO i = 1 , N
            Xn(i) = X(i) - Lambda*Y(i)
         ENDDO

         CALL nef(Irang2,Xn,Fn,Sf,Fnor1,Flag,Flgfft)
C  COMPUTATION OF L-2 NORM OF F+ IS DONE IN NEF
C  CHECK: FT <= FC + ALPHA * LAMBDA * SLOPE
         IF ( (Fnor1-Fnor).LE.(alpha*Lambda*slope) ) THEN
C  THERE IS A GOOD POINT
            Iretcd = 0.D0
C MAXTKN=?
            IF ( Lambda.EQ.1.0D0 .AND. newtln.GT.0.99D0*Maxdu )
     &           Maxtkn = .TRUE.
C_RETURN______________________________________________________________
            RETURN

         ELSEIF ( Lambda.GE.minlbd ) THEN


C  MESSAGE ABOUT THE START OF A SINGLE SEARCH
C   (AT THE FIRST STEP, WHEN LAMBDA=1)
            IF ( Kprsol.GE.2 .AND. Lambda.EQ.1.D0 ) WRITE (6,501)
            IF ( Lambda.EQ.1.D0 ) PRINT 501

C  DECREASE LAMBDA
            IF ( Lambda.LT.1 ) THEN

C  CUBIC INTERPOLATION
               div = 1/(Lambda-lprev)
               v1 = Fnor1 - Fnor - Lambda*slope
               v2 = fprev1 - Fnor - lprev*slope
               a = div*(v1/(Lambda**2)-v2/(lprev**2))
               b = div*(-v1*lprev/(Lambda**2)+v2*Lambda/(lprev**2))
C     DISC=B*B-3.*A*SLOPE
C  IF A=0 - CUBIC INTERPOLATION DEGENERATES INTO QUADRATIC.
               IF ( a.EQ.0.D0 ) ltemp = -slope/(2.D0*b)
C  DEGENERATE INTERPOLATION
               IF ( a.NE.0.D0 ) ltemp = -b/(3.D0*a)
     &                                  + dsqrt((b/(3.D0*a))**2-
     &                                  slope/(3.D0*a))
C  CHECK: LTEMP > 0.5 * LAMBDA
               IF ( ltemp.GT.Lambda/2 ) ltemp = Lambda/2.D0
            ELSE

C  FIRST FRACTION. QUADRATIC INTERPOLATION
               ltemp = -slope/(2*(Fnor1-Fnor-slope))
            ENDIF

C  UPDATE: LTEMP CALCULATED.
            lprev = Lambda
            fprev1 = Fnor1

C  CHECK: LTEMP <= 0.1 * LAMBDA
            IF ( ltemp.LE.0.1D0*Lambda ) ltemp = 0.1D0*Lambda
            Lambda = ltemp

C  MESSAGE ABOUT THE START OF THE SEARCH
            IF ( Kprsol.GE.2 ) WRITE (6,502) Fnor1 , Lambda
            PRINT 502 , Fnor1 , Lambda
C  DECREASE LAMBDA OF THE CURRENT STEP
            IF ( Iretcd.LT.2 ) RETURN
         ELSE
C  GOOD POINT WAS NOT FOUND
            Iretcd = 1.D0
            RETURN
         ENDIF
      ENDDO
C     DEBUG SUBTRACE,INIT(NEWTLN,REL,MINLBD,RELLEN,FNOR1,FNOR,SLOPE,
C    *              DISC, LTEMP,LAMBDA,LPREV,DIV,A,B)
 501  FORMAT (15X,'  ONE-DIM SEARCH  : ')
 502  FORMAT (15X,'   FNOR=',E13.6,',  LAMBDA=',E13.6)

      END
