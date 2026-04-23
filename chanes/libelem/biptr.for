!*==BIPTR1.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE biptr1(Ivar)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER Mt
      INTEGER*4 klc , klv , knl , Ivar
      COMMON /mdla  / Mt(15)
      Mt(1) = Ivar
      Mt(2) = Ivar
      Mt(3) = Ivar
      Mt(4) = Ivar
      Mt(5) = 2
      Mt(6) = 2
      Mt(7) = 2
      klc = 4*(1-Ivar)
      klv = 4*Ivar
      knl = 3
      END



      SUBROUTINE biptr2(Om,P1,L1,P2,L2,P3,L3)
C
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION alb , ale , alk , c14 , c24 , c34 , g1 , g2 ,
     &                 g3 , Om , x , zb , ze , zk
      INTEGER i , k , L1 , L2 , L3
      COMMON /subc  / Y(15,15) , J(15)
      DOUBLE PRECISION P1(L1) , P2(L2) , P3(L3)
      DOUBLE COMPLEX Y , J , ci , gb , gk , ge , zero


c      CR(X) = DCMPLX(X,   0.0D0)
      ci(x) = dcmplx(0.0D0,Om*x)

      zero = dcmplx(0.0D0,0.0D0)

      g1 = 1.D0/P3(1)
      g2 = 1.D0/P3(2)
      g3 = 1.D0/P3(3)

      c14 = P3(4)
      c24 = P3(5)
      c34 = P3(6)

      ale = P3(7)
      alk = P3(8)
      alb = P3(9)

      zb = P3(1)**2 + (Om*alb)**2
      zk = P3(2)**2 + (Om*alk)**2
      ze = P3(3)**2 + (Om*ale)**2

      gb = dcmplx(g1/zb,-Om*alb/zb)
      gk = dcmplx(g2/zk,-Om*alk/zk)
      ge = dcmplx(g3/ze,-Om*alk/ze)


      DO i = 1 , 7
         DO k = 1 , 7
            Y(i,k) = zero
         ENDDO
      ENDDO

      Y(1,1) = gb + ci(c14)
      Y(2,2) = gk + ci(c24)
      Y(3,3) = ge + ci(c34)

      Y(1,5) = -gb
      Y(5,1) = Y(1,5)

      Y(1,4) = -ci(c14)
      Y(4,1) = Y(1,4)

      Y(2,6) = -gk
      Y(6,2) = Y(2,6)

      Y(2,4) = -ci(c24)
      Y(4,2) = Y(2,4)

      Y(3,7) = -ge
      Y(7,3) = Y(3,7)

      Y(3,4) = -ci(c34)
      Y(4,3) = Y(3,4)

      END



      SUBROUTINE biptr3(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)
C
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION alfi , alfn , cbe0 , cbk0 , fi0e , fi0k , taue ,
     &                 tauk , tete , tetk , u , un1 , un10 , un2 , un3 ,
     &                 un4 , un5 , un6 , un7 , un8
      DOUBLE PRECISION un9
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1(L1) , P2(L2) , P3(L3) , B1(Knc2,Nr) , ie0 ,
     &                 ik0 , ne , nk , argmax
      DATA argmax/40.0D0/
C
C
      un1(u) = ie0*(dexp(tete*u)-1.D0)
      un2(u) = ik0*(dexp(tetk*u)-1.D0)

      un3(u) = cbe0/((1.D0+u/fi0e)**ne)
      un4(u) = cbe0*(1.D0+ne*u/fi0e)

      un5(u) = taue*tete*(ie0*(dexp(tete*u)-1.D0)+ie0)

      un6(u) = cbk0/((1.D0+u/fi0k)**nk)
      un7(u) = cbk0*(1.D0+nk*u/fi0k)

      un8(u) = tauk*tetk*(ik0*(dexp(tetk*u)-1.D0)+ik0)

      un9(u) = alfn*ie0*(dexp(tete*u)-1.D0)
      un10(u) = alfi*ik0*(dexp(tetk*u)-1.D0)

      ie0 = P3(10)
      tete = P3(11)
      ik0 = P3(12)
      tetk = P3(13)
      cbe0 = P3(14)
      fi0e = P3(15)
      ne = P3(16)
      taue = P3(17)
      cbk0 = P3(18)
      fi0k = P3(19)
      nk = P3(20)
      tauk = P3(21)
      alfn = P3(22)
      alfi = P3(23)

C -------- NG=1 (CONTROL VOLTAGE UBE)
C -------- NG=2 (CONTROL VOLTAGE UBC)

      IF ( Ng.NE.1 ) THEN
C
C
         DO k = 1 , Knc2
            u = B1(k,1)
            IF ( u*tetk.GT.argmax ) u = argmax/tetk
            IF ( u.GE.0.0D0 ) THEN
C
C -------- UBK>=0.0
               B1(k,1) = un2(u) + un7(u)*B1(k,2) + un8(u)*B1(k,2)
     &                   + un10(u)
            ELSE
C
C -------- UBK<0.0
               B1(k,1) = un2(u) + un6(u)*B1(k,2) + un8(u)*B1(k,2)
     &                   + un10(u)
            ENDIF

            B1(k+1,1) = 0.0D0

         ENDDO
         RETURN
      ENDIF

      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         IF ( u*tete.GT.argmax ) u = argmax/tete
         IF ( u.GE.0.0D0 ) THEN
C
C -------- UBE>=0.0
            B1(k,1) = un1(u) + un4(u)*B1(k,2) + un5(u)*B1(k,2) + un9(u)
         ELSE
C
C -------- UBE<0.0
            B1(k,1) = un1(u) + un3(u)*B1(k,2) + un5(u)*B1(k,2) + un9(u)
         ENDIF

         B1(k+1,1) = 0.0D0

      ENDDO
      RETURN



      END



      SUBROUTINE biptr4(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)
C
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION alfi , alfn , cbe0 , cbk0 , dun1 , dun10 , dun2 ,
     &                 dun3 , dun4 , dun5 , dun6 , dun7 , dun8 , dun9 ,
     &                 fi0e , fi0k , taue , tauk , tete , tetk
      DOUBLE PRECISION u , un1 , un2 , un3 , un4 , un5 , un6 , un7 , un8
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1(L1) , P2(L2) , P3(L3) , B1(Knc2,Nr) , ie0 ,
     &                 ik0 , ne , nk , argmax , argmin
      DATA argmax/40.0D0/ , argmin/ - 120.0D0/
C
C
C
C
      dun1(u) = ie0*tete*dexp(tete*u)
      dun2(u) = ik0*tetk*dexp(tetk*u)
      dun3(u) = (-ne*cbe0*(1.D0+u/fi0e)**(-(ne+1.D0)))/fi0e
      dun4(u) = cbe0*ne/fi0e
      dun5(u) = taue*tete*dun1(u)
      dun6(u) = (-nk*cbk0*(1.D0+u/fi0k)**(-(nk+1)))/fi0k
      dun7(u) = cbk0*nk*u/fi0k
      dun8(u) = tauk*tetk*dun2(u)
      dun9(u) = alfn*dun1(u)
      dun10(u) = alfi*dun2(u)
C
C
      un1(u) = ie0*(dexp(tete*u)-1.D0)
      un2(u) = ik0*(dexp(tetk*u)-1.D0)
      un3(u) = cbe0/((1.D0+u/fi0e)**ne)
      un4(u) = cbe0*(1.D0+ne*u/fi0e)
      un5(u) = taue*tete*(un1(u)+ie0)
      un6(u) = cbk0/((1.D0+u/fi0k)**nk)
      un7(u) = cbk0*(1.D0+nk*u/fi0k)
      un8(u) = tauk*tetk*(un2(u)+ik0)
C
C
      ie0 = P3(10)
      tete = P3(11)
      ik0 = P3(12)
      tetk = P3(13)
      cbe0 = P3(14)
      fi0e = P3(15)
      ne = P3(16)
      taue = P3(17)
      cbk0 = P3(18)
      fi0k = P3(19)
      nk = P3(20)
      tauk = P3(21)
      alfn = P3(22)
      alfi = P3(23)
C
C -------- NG=1 (CONTROL VOLTAGE UBE)
C -------- NG=2 (CONTROL VOLTAGE UBC)

      IF ( Ng.NE.1 ) THEN

         DO k = 1 , Knc2 , 2
            u = B1(k,1)
            IF ( u*tetk.GT.argmax ) u = argmax/tetk
            IF ( u.GE.0.0D0 ) THEN

               IF ( u*tetk.LT.argmin ) THEN
C
C -------- U TOO HIGH
                  B1(k,1) = dun7(u)*B1(k,2) + dun8(u)*B1(k,2) + dun10(u)
                  B1(k,2) = un7(u) + un8(u)
               ELSE
C
C -------- NORMAL CASE
                  B1(k,1) = dun2(u) + dun7(u)*B1(k,2) + dun8(u)*B1(k,2)
     &                      + dun10(u)
                  B1(k,2) = un7(u) + un8(u)
               ENDIF
            ELSEIF ( u*tetk.LT.argmin ) THEN
C
C -------- U TOO HIGH
               B1(k,1) = dun6(u)*B1(k,2) + dun8(u)*B1(k,2) + dun10(u)
               B1(k,2) = un6(u) + un8(u)
            ELSE
C
C -------- NORMAL CASE
               B1(k,1) = dun2(u) + dun6(u)*B1(k,2) + dun8(u)*B1(k,2)
     &                   + dun10(u)
               B1(k,2) = un6(u) + un8(u)
            ENDIF

            B1(k+1,1) = 0.0D0
            B1(k+1,2) = 0.0D0

         ENDDO
         RETURN
      ENDIF

      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         IF ( u*tete.GT.argmax ) u = argmax/tete
         IF ( u.GE.0.0D0 ) THEN

            IF ( u*tete.LT.argmin ) THEN
C
C -------- U TOO HIGH
               B1(k,1) = dun4(u)*B1(k,2) + dun5(u)*B1(k,2) + dun9(u)
               B1(k,2) = un4(u) + un5(u)
            ELSE
C
C -------- NORMAL CASE
               B1(k,1) = dun1(u) + dun4(u)*B1(k,2) + dun5(u)*B1(k,2)
     &                   + dun9(u)
               B1(k,2) = un4(u) + un5(u)
            ENDIF
         ELSEIF ( u*tete.LT.argmin ) THEN
C
C -------- U TOO HIGH
            B1(k,1) = dun3(u)*B1(k,2) + dun5(u)*B1(k,2) + dun9(u)
            B1(k,2) = un3(u) + un5(u)
         ELSE
C
C -------- NORMAL CASE
            B1(k,1) = dun1(u) + dun3(u)*B1(k,2) + dun5(u)*B1(k,2)
     &                + dun9(u)
            B1(k,2) = un3(u) + un5(u)
         ENDIF

         B1(k+1,1) = 0.0D0
         B1(k+1,2) = 0.0D0
      ENDDO
      RETURN

      END



      SUBROUTINE biptr5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nriv,Nbiv)
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER Noi(5,2) , Nou(5,2)
      LOGICAL Exist(5,2)
      INTEGER Koi , Kouv(5) , Kopv(5) , Nriv(5) , Nbiv(5)
C
C
      Koi = 4
C
C -------- CGK, CBK, I2
      Nou(1,1) = 2
      Nou(1,2) = 1
      Noi(1,1) = 6
      Noi(1,2) = 5

      Exist(1,1) = .TRUE.
      Exist(1,2) = .TRUE.

      Kouv(1) = 1
      Kopv(1) = 1
      Nriv(1) = 2
      Nbiv(1) = 2
C
C -------- CGE, CBE, I1
      Nou(2,1) = 3
      Nou(2,2) = 1
      Noi(2,1) = 7
      Noi(2,2) = 5

      Exist(2,1) = .TRUE.
      Exist(1,2) = .TRUE.

      Kouv(2) = 1
      Kopv(2) = 1
      Nriv(2) = 2
      Nbiv(2) = 2
C
C -------- ALFN*I1
      Nou(3,1) = 2
      Nou(3,2) = 1
      Noi(3,1) = 6
      Noi(3,2) = 5

      Exist(3,1) = .TRUE.
      Exist(3,2) = .TRUE.

      Kouv(3) = 1
      Kopv(3) = 0
      Nriv(3) = 1
      Nbiv(3) = 1
C
C -------- ALFI*I2
      Nou(4,1) = 3
      Nou(4,2) = 1
      Noi(4,1) = 7
      Noi(4,2) = 5

      Exist(4,1) = .TRUE.
      Exist(4,2) = .TRUE.

      Kouv(4) = 1
      Kopv(4) = 0
      Nriv(4) = 1
      Nbiv(4) = 1

      END
