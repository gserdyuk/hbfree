c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



       SUBROUTINE LINE1(OM,P1,L1,P2,L2,P3,L3)

C
C MPOINT=
C 'LIB0','LL0 ','    ','    ',?,2,12,4,0,?,1,?,2,?,?,0,0,0,0,0
C NODEEL=
C 'NAME', 4, 2, 1, NODE1, NODE2, NODE3, NODE4, '    ', ?, ?, ?
C PARAM = R K.3., (- COMMON FOR TYPE)
C        Z0,LENTH

C      TRANSMISSION LINE MODEL WITHOUT LOSS WITH NON-BOUNDARY Y1 0-I-------I-0 Y2
C                 WITH COMMON CURRENT FLOW
C                                                             Y3 0-I-------I-0 Y4
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION OM,P1,P2,P3
      DOUBLE PRECISION RSH,Z0,LENTH,TETA,C
      DIMENSION P1(L1),P2(L2),P3(L3)
      COMMON/SUBC/YL(15,15),VJ(15)
      DOUBLE COMPLEX YL,VJ
      DOUBLE COMPLEX A,B
C     TRANSMISSION LINE WITHOUT LOSS. PARAMETERS: R K.3., Z0, LENGTH
C         C = CONSTANT CBETA.
      RSH=P2(1)
      Z0=P3(1)
      LENTH=P3(2)
      C=300 000 000.0D0

      TETA=OM*LENTH/C
      Z0SIN=Z0*DSIN(TETA)
      IF(DABS(Z0SIN).LT.RSH)B=-1.0D0/DCMPLX(RSH,Z0SIN)
      IF(DABS(Z0SIN).GE.RSH)B=-1.0D0/DCMPLX(0.0D0,Z0SIN)
      A= DCMPLX(DCOS(TETA),0.0D0)*(-B)

      YL(1,1)= A
      YL(1,2)= B
      YL(1,3)=-A
      YL(1,4)=-B
      YL(2,1)= B
      YL(2,2)= A
      YL(2,3)=-B
      YL(2,4)=-A
      YL(3,1)=-A
      YL(3,2)=-B
      YL(3,3)= A
      YL(3,4)= B
      YL(4,1)=-B
      YL(4,2)=-A
      YL(4,3)= B
      YL(4,4)= A
C
C   POSSIBILITY OF CHANGE A-B = -SIN(TETA/2.) / IM * Z0 * COS(TETA/2.)
C
      RETURN
C     DEBUG INIT
      END

      SUBROUTINE LINE2(OM,P1,L1,P2,L2,P3,L3)

C MPOINT=
C 'LIB1','LL0 ','    ','    ',?,2,10,2,0,?,1,?,2,?,?,0,0,0,0,0
C NODEEL=
C 'NAME', 4, 2, 1, NODE1, NODE2, '    ', ?, ?, ?
C PARAM = R K.3., (- COMMON FOR TYPE)
C        Z0,LENTH

C      TRANSMISSION LINE MODEL WITHOUT LOSS WITH BOUNDARY Y1 0-I-------I-0 Y2
C                 WITH COMMON CURRENT FLOW
C                                                            !-I-------I-!
C                                                           _!_         _!_
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION OM,P1,P2,P3
      DOUBLE PRECISION RSH,Z0,LENTH,TETA,C
      DIMENSION P1(L1),P2(L2),P3(L3)
      COMMON/SUBC/YL(15,15),VJ(15)
      DOUBLE COMPLEX YL,VJ
      DOUBLE COMPLEX A,B
C     TRANSMISSION LINE WITHOUT LOSS. PARAMETERS: R K.3., Z0, LENGTH
C         C = CONSTANT CBETA.

      RSH=P2(1)
      Z0=P3(1)
      LENTH=P3(2)
      C=300 000 000.0D0

      TETA= OM*LENTH/C
      Z0SIN=Z0*DSIN(TETA)
      IF(DABS(Z0SIN).LT.RSH)B=-1.0D0/DCMPLX(RSH,Z0SIN)
      IF(DABS(Z0SIN).GE.RSH)B=-1.0D0/DCMPLX(0.0D0,Z0SIN)
      A= DCMPLX(DCOS(TETA),0.0D0)*(-B)

      YL(1,1)= A
      YL(1,2)= B
      YL(2,1)= B
      YL(2,2)= A


      RETURN
C     DEBUG INIT
      END

