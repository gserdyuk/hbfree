!*==CURT1.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE curt1(Ivar)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER Ivar
      INTEGER*4 Mt
      INTEGER*4 klc , klv , knl
C
      COMMON /mdla  / Mt(15)
C
      Mt(1) = Ivar
      Mt(2) = Ivar
      Mt(3) = Ivar
      Mt(4) = 2
      Mt(5) = 2
      Mt(6) = 2
C
      klc = (1-Ivar)*4
      klv = 4 - klc
      knl = 4
C
      END

      SUBROUTINE curt2(Om,P1,L1,P2,L2,P3,L3)
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION c12 , c13 , g1 , g2 , g23 , g3 , g31 , Om , x
      INTEGER i , k , L1 , L2 , L3
      DOUBLE PRECISION P1 , P2 , P3
      DOUBLE COMPLEX Y , J , zero/(0.0D0,0.0D0)/ , cr , ci
C
      COMMON /subc  / Y(15,15) , J(15)
C
      DIMENSION P1(L1) , P2(L2) , P3(L3)
C
      cr(x) = dcmplx(x,0.0D0)
      ci(x) = dcmplx(0.0D0,Om*x)
C
      g1 = 1.D0/P3(1)
      g2 = 1.D0/P3(2)
      g3 = 1.D0/P3(3)
      IF ( Om.EQ.0.D0 ) THEN
         g31 = 0.000001D0
         g23 = 0.000001D0
      ELSE
         g31 = 1.D0/P3(14)
         g23 = 1.D0/P3(15)
      ENDIF
      c12 = P3(4)
      c13 = P3(5)
C
      DO i = 1 , 8
         DO k = 1 , 8
            Y(k,i) = zero
         ENDDO
      ENDDO
C
      Y(1,1) = cr(g2)
      Y(1,5) = -cr(g2)
      Y(2,2) = cr(g1)
      Y(2,4) = -cr(g1)
      Y(3,3) = cr(g3)
      Y(3,6) = -cr(g3)
      Y(4,2) = Y(2,4)
      Y(4,4) = cr(g1) + ci(c12) + ci(c13) + cr(g31)
      Y(4,5) = -ci(c12)
      Y(4,6) = -ci(c13) - cr(g31)
      Y(5,1) = Y(1,5)
      Y(5,4) = Y(4,5)
      Y(5,5) = cr(g2) + ci(c12) + cr(g23)
      Y(5,6) = -cr(g23)
      Y(6,5) = Y(5,6)
      Y(6,3) = Y(3,6)
      Y(6,4) = Y(4,6)
      Y(6,6) = cr(g3) + ci(c13)

C
      END

      SUBROUTINE curt3(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION al , as , bt , c230 , th , u , u1 , u2 , un1 ,
     &                 un2 , un3 , vbi , vbi0 , vt
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DIMENSION B1(Knc2,Nr)
      DIMENSION P1(L1) , P2(L2) , P3(L3)
C
      DOUBLE PRECISION P1 , P2 , P3
      DOUBLE PRECISION is , ld
      DOUBLE PRECISION B1
      DOUBLE PRECISION ueps/0.05D0/ , argmax/40.0D0/
C
      th(u) = dtanh(al*u)
      un1(u) = is*(dexp(as*u)-1.D0)
      un2(u) = c230/dsqrt(1.0D0-u/vbi)
      un3(u1,u2) = bt*(u1+vt)*(u1+vt)*(1.D0+ld*u2)*th(u2)
C
      al = P3(6)
      bt = P3(7)
      ld = P3(8)
      vt = P3(9)
      is = P3(10)
      as = P3(11)
      c230 = P3(12)
      vbi = P3(13)
C     UEPS=0.05
      vbi0 = vbi - ueps
C
      IF ( Ng.NE.1 ) THEN
C
         DO k = 1 , Knc2 , 2
            u1 = B1(k,1)
            u2 = B1(k,2)
C
            IF ( (u1+vt).LE.0.0D0 ) THEN
               B1(k,1) = 0.0D0
            ELSE
               B1(k,1) = un3(u1,u2)
            ENDIF
            B1(k+1,1) = 0.0D0
         ENDDO
         RETURN
      ENDIF
C
      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         IF ( u*as.GT.argmax ) u = argmax/as
         IF ( u.GT.vbi0 ) THEN
C
C     IF U > VBI0 = VBI - UEPS, THEN THE CAPACITANCE IS CALCULATED BY
C     THE FORMULA C(U) = C(VBI0) + DC/DU(U=VBI0) * (U - VBI0)
C
            B1(k,1) = un1(u) + un2(vbi0)*(1.D0+(u-vbi0)/(2.D0*ueps))
     &                *B1(k,2)
         ELSE
C     THE DISPLACEMENT CURRENT IS DEFINED AS C(U)*DU/DT AND NOT
C     D(C(U)*U)/DT=(C(U)+U*DC/DU)*DU/DT
C
            B1(k,1) = un1(u) + un2(u)*B1(k,2)
         ENDIF
C
         B1(k+1,1) = 0.0D0

      ENDDO
      RETURN
C
      END

      SUBROUTINE curt4(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION al , as , bt , c230 , dun1 , dun2 , dun3 , th ,
     &                 u , u1 , u2 , un2 , un3 , vbi , vbi0 , vt
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DIMENSION B1(Knc2,Nr)
C
      DOUBLE PRECISION P1 , P2 , P3
      DOUBLE PRECISION B1
      DOUBLE PRECISION is , ld
      DOUBLE PRECISION ueps/0.05D0/ , argmin/ - 120.0D0/ ,
     &                 argmax/40.0D0/
C
      th(u) = dtanh(al*u)
      dun1(u) = is*as*dexp(as*u)
      un2(u) = c230/dsqrt(1.0D0-u/vbi)
      un3(u) = un2(u)/2*(vbi-u)
      dun2(u1,u2) = 2.D0*bt*(u1+vt)*(1.D0+ld*u2)*th(u2)
      dun3(u1,u2) = bt*(u1+vt)*(u1+vt)
     &              *((1.D0+ld*u2)*al*(1.D0-th(u2)*th(u2))+ld*th(u2))
C
      al = P3(6)
      bt = P3(7)
      ld = P3(8)
      vt = P3(9)
      is = P3(10)
      as = P3(11)
      c230 = P3(12)
      vbi = P3(13)
C     UEPS=0.05
      vbi0 = vbi - ueps
C
      IF ( Ng.NE.1 ) THEN
C
         DO k = 1 , Knc2 , 2
            u1 = B1(k,1)
            u2 = B1(k,2)
            IF ( (u1+vt).LE.0.0D0 ) THEN
C
               B1(k,1) = 0.0D0
               B1(k,2) = 0.0D0
            ELSE
               B1(k,1) = dun2(u1,u2)
               B1(k,2) = dun3(u1,u2)
            ENDIF
C
            B1(k+1,1) = 0.0D0
            B1(k+1,2) = 0.0D0
         ENDDO
         RETURN
      ENDIF
      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         IF ( u*as.GT.argmax ) u = argmax/as
         IF ( u.GT.vbi0 ) THEN
C  U IS CLOSE TO VBI. UN2(U) -> K SINGULARITY
            B1(k,1) = dun1(u) + un2(vbi0)/(2.D0*ueps)*B1(k,2)
            B1(k,2) = un2(vbi0)*(1.D0+(u-vbi0)/(2.D0*ueps))
         ELSEIF ( u*as.LT.argmin ) THEN
C  U IS TOO SMALL ( EXP(AS*U) <= 1.E-60 )
            B1(k,1) = un3(u)*B1(k,2)
            B1(k,2) = un2(u)
         ELSE
C  NORMAL CASE:
            B1(k,1) = dun1(u) + un3(u)*B1(k,2)
            B1(k,2) = un2(u)
         ENDIF
C
         B1(k+1,1) = 0.0D0
         B1(k+1,2) = 0.0D0
C
      ENDDO
      RETURN
C
C
      END

      SUBROUTINE curt5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL Exist(5,2)
      INTEGER Noi(5,2) , Nou(5,2)
      INTEGER Koi , Kouv(5) , Kopv(5) , Nr1v(5) , Nb1v(5)
C
      Koi = 2
C
      Noi(1,1) = 5
      Noi(1,2) = 6
C
      Nou(1,1) = 5
      Nou(1,2) = 6
C
      Exist(1,1) = .TRUE.
      Exist(1,2) = .TRUE.
C
      Kouv(1) = 1
      Kopv(1) = 1
C
      Nr1v(1) = 2
      Nb1v(1) = 2
C
      Noi(2,1) = 4
      Noi(2,2) = 6
C
      Nou(2,1) = 5
      Nou(2,2) = 6
      Nou(3,1) = 4
      Nou(3,2) = 6
C
      Kouv(2) = 2
      Kopv(2) = 0
C
      Nr1v(2) = 2
      Nb1v(2) = 2
C
      Exist(2,1) = .TRUE.
      Exist(2,2) = .FALSE.
      Exist(3,1) = .TRUE.
      Exist(3,2) = .FALSE.
C
      END

      SUBROUTINE curt6(Ng,P1,L1,P2,L2,P3,L3,Val,Dval,Kn,Nr,T)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , Kn , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION T , unew1 , uold1 , zabs
      DOUBLE PRECISION P1 , P2 , P3
      DOUBLE COMPLEX Val , Dval
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DIMENSION Val(Kn,Nr) , Dval(Kn,Nr)
      DOUBLE PRECISION uold , unew , ustep
      DOUBLE PRECISION delta , ubound , stepm

      delta = 0.01D0
C      UBOUND=0.6
C      STEPM =0.3
      ubound = P3(16)
      stepm = P3(17)
C      UBOUND י STEPM גשלי 0.5 י 0.1 .
      IF ( Ng.NE.1 ) RETURN
      uold = 0.0D0
      unew = 0.0D0
      IF ( Kn.GE.2 ) THEN
         DO i = 2 , Kn
            uold = uold + zabs(Val(i,1))
            unew = unew + zabs(Val(i,1)+Dval(i,1))
         ENDDO
      ENDIF

      uold = uold + uold + dreal(Val(1,1))
      uold1 = uold
      unew = unew + unew + dreal(Val(1,1)+Dval(1,1))
      unew1 = unew
      ustep = unew - uold
      IF ( unew.LE.ubound ) RETURN
      IF ( uold.GT.ubound ) THEN
         IF ( ustep.GT.stepm ) T = stepm/ustep
         RETURN
      ENDIF
      T = (ubound+delta-uold)/ustep
      RETURN
C
      END
