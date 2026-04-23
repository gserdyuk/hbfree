!*==INITI.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE initi(U,Yy,Vectj,Isize_maxnode,Ierr)

C*********************************************************************
C  SUBROUTINE 'START' - CALLING FOR THE PROCESSING, CALCULATION OF
C  THE FREQUENCY GRID AND FORMATION OF Y AND J ON THE FREQUENCY GRID.
C  IN CASE OF A LINEAR SCHEME
C*                    - CALCULATES THE SOLUTION OF THE EQUATIONS
C*********************************************************************

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER ii , io , ioo , irr , irro , Isize_maxnode , k12 , k23 ,
     &        k3 , Kmni , Kn , Kn1 , Knc , Knr , Knr1 , maxsyn , Mni ,
     &        mnmax
C   CHANGES MADE SIMILAR TO MAIN 30.01.91 (SERDYUK G.V.)
C $LARGE: J,Y
C      COMMON/MATY/    J(15,20),Y(15,15,20)
C     COMPLEX          Y,J

      DOUBLE PRECISION U(1)

      DOUBLE COMPLEX Yy(Isize_maxnode,Isize_maxnode)
      DOUBLE COMPLEX Vectj(Isize_maxnode)

      COMMON /blw1  / W , W1
      COMMON /fre   / F(2)
      DOUBLE PRECISION W(20) , W1(200) , F
      COMMON /blk1  / Kr , Kr1 , Kc , Kc1 , Nnr , Nnr1 , Mn , Mn1
      COMMON /blk2  / Knc , Knr , Kn , Knr1 , Kn1
      INTEGER Kr(20) , Kc(10) , Nnr(10) , Kr1(200) , Kc1(20) , Nnr1(20)
      INTEGER Mn(2,20) , Mn1(2,200)
      COMMON /kolnal/ Kol , Nal
      INTEGER Kol(4) , ip/6/
      LOGICAL Nal(4)
      COMMON /blmni / Mni(2,20) , Kmni
      INTEGER Ierr



C     CALL OF THE CORRECTION SUBROUTINE AND RECALCULATION OF NODE COORDINATES
      CALL sortuz

C     CHECKING FOR THE PRESENCE OF BOUNDARY ELEMENTS
C     IF(KOL(3)-KOL(4).EQ.0) WRITE(IP, 1020)

C     SUBROUTINE FOR DETERMINING THE MAXIMUM NUMBER OF INPUT/OUTPUT VARIABLES FOR
C           TRANSFORMATION FUNCTION
      CALL detsyn(maxsyn)

C     CALCULATION OF FREQUENCY GRIDS AND PR. (SEE COMMENTS IN SUBROUTINE)
      DO ii = 1 , Kn
         Mn(1,ii) = Mni(1,ii)
         Mn(2,ii) = Mni(2,ii)
      ENDDO

      CALL inkoor(maxsyn,mnmax,F(1),F(2),*100,*101)


C     FORMATION OF Y AND J AT EACH FREQUENCY GRID POINT (WITH W ACCURACY)

      DO io = 1 , Kn
         ioo = io
         CALL stepfr(W(ioo),ioo,Yy,Vectj,Isize_maxnode)
      ENDDO
C     TRANSFORMATION OF NODE CURRENTS
      k3 = Kol(3)
      k23 = Kol(2) + k3
      k12 = Kol(1) + Kol(2)
      CALL topo(-1,k12)
C     ASSIGNMENT OF TOPO-X-K FOR MULTIPLEXER IN ACCURACY MPOINT & NODEEL
      CALL topoin(k23,1)
      CALL topoin(k3,2)

C     IF THERE ARE NO BOUNDARY ELEMENTS AND GRAPHICAL NODES
C            - CALCULATION OF THE SOLUTION OF THE EQUATIONS
      IF ( Nal(3) .OR. Nal(4) ) THEN
         Ierr = 0
         RETURN
      ENDIF

C      WRITE(IP, 1020)
      DO irro = 1 , Kn
         irr = irro
         CALL stback(U,irr)
      ENDDO

      Ierr = 1
      RETURN

C***** LARGE GRID MN1 (MORE THAN 200) *****************************
 100  WRITE (ip,1000)

C**********************************************************************
C      ENTRY AGAIN
C this entry should not used any more
C
C     ENTRY FOR REPEATED CALL OF THE FORMATION OF CHANGING ELEMENTS
C
C
C     CALCULATION OF NODE CURRENTS
C      K12=KOL(1)+KOL(2)
C      CALL TOPO(1,K12)
C
C
C     REPEATED FORMATION OF THE BOUNDARY PART
C      DO 60 IO=1,KN
C      IOO=IO
C   60 CALL DOUBLE(W(IOO),IOO)
C
C     TRANSFORMATION OF NODE FUNCTIONS /AGAIN/
C      CALL TOPO(-1,K12)
C
C     ASSIGNMENT OF TOPO-X-K
C      K3=KOL(3)
C      CALL TOPOIN(K3,2)
C
C     IF THERE ARE NO BOUNDARY ELEMENTS AND GRAPHICAL NODES
C            - CALCULATION OF THE SOLUTION OF THE EQUATIONS
C      IF(NAL(3).OR.NAL(4)) RETURN
C      WRITE(IP, 1020)
C      DO 70 IRRO=1,KN
C      IRR= IRRO
C      CALL STBACK(U,IRR)
C   70 CONTINUE
C
C
C      RETURN 1

 1000 FORMAT (/10X,'SIZE OF AUX. FREQUENCY GRID IS MORE THAN MAX=200.')

      Ierr = 2
      RETURN

C***** KNC > KNCMAX - EXIT INTO THE FORBIDDEN AREA *******************
 101  WRITE (ip,1010) mnmax
 1010 FORMAT (/10X,'CANT COMPUTE HARMONICS HIGHER THAN ',I4)

      Ierr = 3
C
      END
