c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE BIPTR1(IVAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER*4 KLC,KLV,KNL,IVAR
      COMMON/MDLA/MT(15)
      MT(1)=IVAR
      MT(2)=IVAR
      MT(3)=IVAR
      MT(4)=IVAR
      MT(5)=2
      MT(6)=2
      MT(7)=2
      KLC=4*(1-IVAR)
      KLV=4*IVAR
      KNL=3
      RETURN
      END



      SUBROUTINE BIPTR2(OM,P1,L1,P2,L2,P3,L3)
C
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/SUBC/      Y(15,15),J(15)
      DOUBLE PRECISION              P1(L1),P2(L2),P3(L3)
      DOUBLE COMPLEX           Y,J,                  CR,CI,GB,GK,GE,
     +              ZERO


c      CR(X) = DCMPLX(X,   0.0D0)
      CI(X) = DCMPLX(0.0D0,OM*X)

      ZERO = DCMPLX(0.0D0,0.0D0)

      G1 =1.D0/P3(1)
      G2 =1.D0/P3(2)
      G3 =1.D0/P3(3)

      C14=P3(4)
      C24=P3(5)
      C34=P3(6)

      ALE=P3(7)
      ALK=P3(8)
      ALB=P3(9)

      ZB =P3(1)**2+(OM*ALB)**2
      ZK =P3(2)**2+(OM*ALK)**2
      ZE =P3(3)**2+(OM*ALE)**2

      GB =DCMPLX(G1/ZB,-OM*ALB/ZB)
      GK =DCMPLX(G2/ZK,-OM*ALK/ZK)
      GE =DCMPLX(G3/ZE,-OM*ALK/ZE)


      DO 10  I=1,7
      DO 10  K=1,7
      Y(I,K)=ZERO
   10 CONTINUE

      Y(1,1) = GB+CI(C14)
      Y(2,2) = GK+CI(C24)
      Y(3,3) = GE+CI(C34)

      Y(1,5) = -GB
      Y(5,1) = Y(1,5)

      Y(1,4) = -CI(C14)
      Y(4,1) = Y(1,4)

      Y(2,6) = -GK
      Y(6,2) = Y(2,6)

      Y(2,4) = -CI(C24)
      Y(4,2) = Y(2,4)

      Y(3,7) = -GE
      Y(7,3) = Y(3,7)

      Y(3,4) = -CI(C34)
      Y(4,3) = Y(3,4)

      RETURN
      END



      SUBROUTINE BIPTR3(NG,P1,L1,P2,L2,P3,L3,B1,KNC2,NR,*)
C
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION              P1(L1),P2(L2),P3(L3),
     +   B1(KNC2,NR),                  IE0,IK0,NE,NK,ARGMAX
      DATA              ARGMAX/40.0D0/
C
C
      UN1(U) = IE0*(DEXP(TETE*U)-1.D0)
      UN2(U) = IK0*(DEXP(TETK*U)-1.D0)

      UN3(U) = CBE0/((1.D0+U/FI0E)**NE)
      UN4(U) = CBE0*(1.D0+NE*U/FI0E)

      UN5(U) = TAUE*TETE*(IE0*(DEXP(TETE*U)-1.D0)+IE0)

      UN6(U) = CBK0/((1.D0+U/FI0K)**NK)
      UN7(U) = CBK0*(1.D0+NK*U/FI0K)

      UN8(U) = TAUK*TETK*(IK0*(DEXP(TETK*U)-1.D0)+IK0)

      UN9(U) = ALFN*IE0*(DEXP(TETE*U)-1.D0)
      UN10(U)= ALFI*IK0*(DEXP(TETK*U)-1.D0)

      IE0  =P3(10)
      TETE =P3(11)
      IK0  =P3(12)
      TETK =P3(13)
      CBE0 =P3(14)
      FI0E =P3(15)
      NE   =P3(16)
      TAUE =P3(17)
      CBK0 =P3(18)
      FI0K =P3(19)
      NK   =P3(20)
      TAUK =P3(21)
      ALFN =P3(22)
      ALFI =P3(23)

C -------- NG=1 (CONTROL VOLTAGE UBE)
C -------- NG=2 (CONTROL VOLTAGE UBC)

      IF(NG.NE.1) GOTO 100

      DO 10  K=1,KNC2,2
      U =B1(K,1)
      IF(U*TETE.GT.ARGMAX) U =ARGMAX/TETE
      IF(U.GE.0.0D0) GOTO 5
C
C -------- UBE<0.0
      B1(K,1) =UN1(U)+UN3(U)*B1(K,2)+UN5(U)*B1(K,2)+UN9(U)
      GOTO 9
C
C -------- UBE>=0.0
    5 CONTINUE
      B1(K,1) =UN1(U)+UN4(U)*B1(K,2)+UN5(U)*B1(K,2)+UN9(U)

    9 CONTINUE
      B1(K+1,1) =0.0D0

   10 CONTINUE
      RETURN
C
C
  100 CONTINUE
      DO 30 K=1,KNC2
      U =B1(K,1)
      IF(U*TETK.GT.ARGMAX) U =ARGMAX/TETK
      IF(U.GE.0.0D0) GOTO 20
C
C -------- UBK<0.0
      B1(K,1) =UN2(U)+UN6(U)*B1(K,2)+UN8(U)*B1(K,2)+UN10(U)
      GOTO 25
C
C -------- UBK>=0.0
   20 CONTINUE
      B1(K,1) =UN2(U)+UN7(U)*B1(K,2)+UN8(U)*B1(K,2)+UN10(U)

   25 CONTINUE
      B1(K+1,1) =0.0D0

   30 CONTINUE


      RETURN

      END



      SUBROUTINE BIPTR4(NG,P1,L1,P2,L2,P3,L3,B1,KNC2,NR,*)
C
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION              P1(L1),P2(L2),P3(L3),
     +   B1(KNC2,NR),                  IE0,IK0,NE,NK,ARGMAX,ARGMIN
      DATA              ARGMAX/40.0D0/,ARGMIN/-120.0D0/
C
C
C
C
      DUN1(U) = IE0*TETE*DEXP(TETE*U)
      DUN2(U) = IK0*TETK*DEXP(TETK*U)
      DUN3(U) = (-NE*CBE0*(1.D0+U/FI0E)**(-(NE+1.D0)))/FI0E
      DUN4(U) = CBE0*NE/FI0E
      DUN5(U) = TAUE*TETE*DUN1(U)
      DUN6(U) = (-NK*CBK0*(1.D0+U/FI0K)**(-(NK+1)))/FI0K
      DUN7(U) = CBK0*NK*U/FI0K
      DUN8(U) = TAUK*TETK*DUN2(U)
      DUN9(U) = ALFN*DUN1(U)
      DUN10(U)= ALFI*DUN2(U)
C
C
      UN1(U) = IE0*(DEXP(TETE*U)-1.D0)
      UN2(U) = IK0*(DEXP(TETK*U)-1.D0)
      UN3(U) = CBE0/((1.D0+U/FI0E)**NE)
      UN4(U) = CBE0*(1.D0+NE*U/FI0E)
      UN5(U) = TAUE*TETE*(UN1(U)+IE0)
      UN6(U) = CBK0/((1.D0+U/FI0K)**NK)
      UN7(U) = CBK0*(1.D0+NK*U/FI0K)
      UN8(U) = TAUK*TETK*(UN2(U)+IK0)
C
C
      IE0  = P3(10)
      TETE = P3(11)
      IK0  = P3(12)
      TETK = P3(13)
      CBE0 = P3(14)
      FI0E = P3(15)
      NE   = P3(16)
      TAUE = P3(17)
      CBK0 = P3(18)
      FI0K = P3(19)
      NK   = P3(20)
      TAUK = P3(21)
      ALFN = P3(22)
      ALFI = P3(23)
C
C -------- NG=1 (CONTROL VOLTAGE UBE)
C -------- NG=2 (CONTROL VOLTAGE UBC)

      IF(NG.NE.1) GOTO 100

      DO 10  K=1,KNC2,2
      U =B1(K,1)
      IF(U*TETE.GT.ARGMAX) U =ARGMAX/TETE
      IF(U.GE.0.0D0) GOTO 5
      IF(U*TETE.LT.ARGMIN) GOTO 3
C
C -------- NORMAL CASE
      B1(K,1) =DUN1(U)+DUN3(U)*B1(K,2)+DUN5(U)*B1(K,2)+DUN9(U)
      B1(K,2) =UN3(U)+UN5(U)
      GOTO 8
C
C -------- U TOO HIGH
    3 CONTINUE
      B1(K,1) =DUN3(U)*B1(K,2)+DUN5(U)*B1(K,2)+DUN9(U)
      B1(K,2) =UN3(U)+UN5(U)
      GOTO 8

    5 CONTINUE
      IF(U*TETE.LT.ARGMIN) GOTO 7
C
C -------- NORMAL CASE
      B1(K,1) =DUN1(U)+DUN4(U)*B1(K,2)+DUN5(U)*B1(K,2)+DUN9(U)
      B1(K,2) =UN4(U)+UN5(U)
      GOTO 8
C
C -------- U TOO HIGH
    7 CONTINUE
      B1(K,1) =DUN4(U)*B1(K,2)+DUN5(U)*B1(K,2)+DUN9(U)
      B1(K,2) =UN4(U)+UN5(U)

    8 CONTINUE
      B1(K+1,1) =0.0D0
      B1(K+1,2) =0.0D0
  10  CONTINUE
      RETURN

  100 CONTINUE
      DO 30 K=1,KNC2,2
      U =B1(K,1)
      IF(U*TETK.GT.ARGMAX) U =ARGMAX/TETK
      IF(U.GE.0.0D0) GOTO 25
      IF(U*TETK.LT.ARGMIN) GOTO 23
C
C -------- NORMAL CASE
      B1(K,1) =DUN2(U)+DUN6(U)*B1(K,2)+DUN8(U)*B1(K,2)+DUN10(U)
      B1(K,2) =UN6(U)+UN8(U)
      GOTO 28
C
C -------- U TOO HIGH
   23 CONTINUE
      B1(K,1) =DUN6(U)*B1(K,2)+DUN8(U)*B1(K,2)+DUN10(U)
      B1(K,2) =UN6(U)+UN8(U)
      GOTO 28

   25 CONTINUE
      IF(U*TETK.LT.ARGMIN) GOTO 27
C
C -------- NORMAL CASE
      B1(K,1) =DUN2(U)+DUN7(U)*B1(K,2)+DUN8(U)*B1(K,2)+DUN10(U)
      B1(K,2) =UN7(U)+UN8(U)
      GOTO 28
C
C -------- U TOO HIGH
   27 CONTINUE
      B1(K,1) =DUN7(U)*B1(K,2)+DUN8(U)*B1(K,2)+DUN10(U)
      B1(K,2) =UN7(U)+UN8(U)

   28 CONTINUE
      B1(K+1,1)=0.0D0
      B1(K+1,2)=0.0D0

   30 CONTINUE

      RETURN
      END



      SUBROUTINE BIPTR5(NOI,NOU,EXIST,KOI,KOUV,KOPV,NRIV,NBIV)
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NOI(5,2),NOU(5,2)
      LOGICAL EXIST(5,2)
      INTEGER KOI,KOUV(5),KOPV(5),NRIV(5),NBIV(5)
C
C
      KOI=4
C
C -------- CGK, CBK, I2
      NOU(1,1) = 2
      NOU(1,2) = 1
      NOI(1,1) = 6
      NOI(1,2) = 5

      EXIST(1,1) = .TRUE.
      EXIST(1,2) = .TRUE.

      KOUV(1) = 1
      KOPV(1) = 1
      NRIV(1) = 2
      NBIV(1) = 2
C
C -------- CGE, CBE, I1
      NOU(2,1) = 3
      NOU(2,2) = 1
      NOI(2,1) = 7
      NOI(2,2) = 5

      EXIST(2,1) = .TRUE.
      EXIST(1,2) = .TRUE.

      KOUV(2) = 1
      KOPV(2) = 1
      NRIV(2) = 2
      NBIV(2) = 2
C
C -------- ALFN*I1
      NOU(3,1) = 2
      NOU(3,2) = 1
      NOI(3,1) = 6
      NOI(3,2) = 5

      EXIST(3,1) = .TRUE.
      EXIST(3,2) = .TRUE.

      KOUV(3) = 1
      KOPV(3) = 0
      NRIV(3) = 1
      NBIV(3) = 1
C
C -------- ALFI*I2
      NOU(4,1) = 3
      NOU(4,2) = 1
      NOI(4,1) = 7
      NOI(4,2) = 5

      EXIST(4,1) = .TRUE.
      EXIST(4,2) = .TRUE.

      KOUV(4) = 1
      KOPV(4) = 0
      NRIV(4) = 1
      NBIV(4) = 1

      RETURN
      END
