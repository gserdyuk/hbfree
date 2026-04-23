!*==STOP.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE stop(N,X,Dx,F,Fnor,G,Sx,Sf,Iretcd,Iter,Maxtkn,Kmaxdu,
     &                Termcd)
C
C   SUBROUTINE FOR DETERMINING THE REASON FOR STOPPING.
C     * TERMCD=0 - NO STOP;
C       ...

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION conv , dunor , Fnor , sinv , term
      INTEGER i , Iretcd , Iter , Kmaxdu , Kprgrf , Kprlen , Kprlin ,
     &        Kprnkr , Kprqup , Kprsol , Kprsrt , Kprvar , Limit , N
      COMMON /print / Kprlen , Kprsrt , Kprnkr , Kprlin , Kprsol ,
     &                Kprvar , Kprgrf , Kprqup
      COMMON /newton/ Epssol , Epsdu , Epsmin , Maxdu , Limit
      DOUBLE PRECISION Epssol , Epsdu , Epsmin , Maxdu
      INTEGER Termcd
      LOGICAL Maxtkn
      DOUBLE PRECISION X(1) , Dx(1) , F(1) , G(1)
      DOUBLE PRECISION Sx(1) , Sf(1)

      Termcd = 0

      IF ( Iretcd.NE.1 ) THEN

C  MAXIMUM SCALED RESIDUAL TERM:
         conv = 0.D0
         DO i = 1 , N
            conv = dmax1(conv,Sf((i+1)/2)*dabs(F(i)))
         ENDDO
         IF ( conv.GT.Epssol ) THEN

C  MAXIMUM SCALED CORRECTION TERM
            dunor = 0.D0
            DO i = 1 , N
               sinv = 1/Sx((i+1)/2)
               dunor = dmax1(dunor,dabs(Dx(i)/dmax1(dabs(X(i)),sinv)))
            ENDDO
            IF ( dunor.LE.Epsdu ) THEN
C  CORRECTION IS LESS THAN THE MINIMUM ALLOWED.
               Termcd = 2
               IF ( Kprsol.GT.0 ) WRITE (6,102) Iter , dunor , Epsdu ,
     &              conv , Fnor
               PRINT 102 , Iter , dunor , Epsdu , conv , Fnor
               RETURN

            ELSEIF ( Iter.LT.Limit ) THEN

               IF ( Maxtkn ) THEN
C  A STEP OF LENGTH MAXDU WAS TAKEN
                  Kmaxdu = Kmaxdu + 1
                  IF ( Kmaxdu.GE.5 ) THEN
C  FIVE STEPS OF LENGTH MAXDU WERE TAKEN
                     Termcd = 5
                     IF ( Kprsol.GT.0 ) WRITE (6,105) Iter , Maxdu ,
     &                    conv , Epssol
                     PRINT 105 , Iter , Maxdu , conv , Epssol
                     RETURN
                  ENDIF
               ENDIF

               Kmaxdu = 0
C   DEFINITION OF THE REQUIRED DECREASE IN FNOR
               term = 0.D0
               DO i = 1 , N
                  sinv = 1/Sx((i+1)/2)
                  term = dmax1(term,dabs(G(i))*dmax1(X(i),sinv)
     &                   /dmax1(Fnor,dfloat(N)/2.D0))
C     REMOVED PRINT 222 C WITH FORMAT 8.12.89. KOSCHMANOVA N.W. *****
               ENDDO
               IF ( term.GT.Epsmin ) THEN
C   NOTHING HAPPENED
                  IF ( Kprsol.GE.2 ) WRITE (6,107) conv , dunor , term
                  PRINT 107 , conv , dunor , term
                  RETURN
               ENDIF
            ELSE
C  ITERATION LIMIT EXHAUSTED
               Termcd = 4
               IF ( Kprsol.GT.0 ) WRITE (6,104) Iter , conv , Epssol ,
     &              dunor , Epsdu
               PRINT 104 , Iter , conv , Epssol , dunor , Epsdu
               RETURN
            ENDIF
         ELSE
C  WE HAVE AN APPROXIMATE SOLUTION
C   (IF EPSSOL IS NOT VERY SMALL)
            Termcd = 1
            IF ( Kprsol.GT.0 ) WRITE (6,101) Iter , conv , Epssol , Fnor
            PRINT 101 , Iter , conv , Epssol , Fnor
            RETURN
         ENDIF
      ELSE
C  IN LSERCH WE FAILED TO MAKE A SATISFACTORY STEP
C        SCHA
         Termcd = 3
         IF ( Kprsol.GT.0 ) WRITE (6,103) Iter , Epsdu , Fnor
         PRINT 103 , Iter , Epsdu , Fnor
         RETURN
      ENDIF
C   SO, WE ARE AT THE LOCAL MINIMUM OF FNOR.
      Termcd = 6
      IF ( Kprsol.GT.0 ) WRITE (6,106) Iter , term , Epsmin , Fnor ,
     &                                 conv , Epssol
      PRINT 106 , Iter , term , Epsmin , Fnor , conv , Epssol
      RETURN
 101  FORMAT ('     AT ',I5,
     &        ' -TH ITERATION CONVERGED TO '/'  SOLUTION WITH ERROR <= '
     &        ,E13.6,' ( < ',E13.6,
     &        ' )'/'  1/2 OF SQUARED L-2 NORM OF ERROR =',E13.6)
 102  FORMAT ('     AT ',I5,' -T ITERATION MAX STEP=',E13.6,' ( <',
     &        E13.6,')'/'  ERROR =',
     &        E13.6/'  1/2 OF SQUARED L-2 NORM OF ERROR =',E13.6)
 103  FORMAT ('     AT ',I5,
     &        ' -TH ITERATION CAN NOT MAKE GOOD STEP'/'  > ',
     &        E13.6/'  1/2 OF SQUARED L-2 NORM OF ERROR =',E13.6)
 104  FORMAT ('     ITERATION LIMIT:',I4,' IS REACHED.'/'  ERROR =',
     &        E13.6,'( > ',E13.6,' );'/'  STEP =',E13.6,'( > ',E13.6,
     &        ' ).')
 105  FORMAT ('     AT ',I4,
     &        ' ITERATION , IT IS MADE'/'  5 STEPS OF LENGTH ',
     &        E13.6/'  ERROR =',E13.6,' ( > ',E13.6,' ).')
 106  FORMAT ('  IT IS LOCAL MINIMUM :  '/
     &        '              ITERATION                        ',
     &        I4/'              GRADIENT OF ERROR NORM        ',E13.6,
     &        '(<',E13.6,
     &        ')'/'              1/2 SQUARED  L-2 ERROR NORM   ',
     &        E13.6/'               ERROR                        ',
     &        E13.6/'               REQUIRED                     ',
     &        E13.6)

 107  FORMAT (' ERROR  =',E13.6,', STEP=',E13.6,
     &        ', SPEED OF DECREASING=',E13.6)



C     DEBUG SUBTRACE,INIT
      END
