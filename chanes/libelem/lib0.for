!*==LINE1.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE line1(Om,P1,L1,P2,L2,P3,L3)

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
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER L1 , L2 , L3
      DOUBLE PRECISION z0sin
      DOUBLE PRECISION Om , P1 , P2 , P3
      DOUBLE PRECISION rsh , z0 , lenth , teta , c
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      COMMON /subc  / Yl(15,15) , Vj(15)
      DOUBLE COMPLEX Yl , Vj
      DOUBLE COMPLEX a , b
C     TRANSMISSION LINE WITHOUT LOSS. PARAMETERS: R K.3., Z0, LENGTH
C         C = CONSTANT CBETA.
      rsh = P2(1)
      z0 = P3(1)
      lenth = P3(2)
      c = 300000000.0D0

      teta = Om*lenth/c
      z0sin = z0*dsin(teta)
      IF ( dabs(z0sin).LT.rsh ) b = -1.0D0/dcmplx(rsh,z0sin)
      IF ( dabs(z0sin).GE.rsh ) b = -1.0D0/dcmplx(0.0D0,z0sin)
      a = dcmplx(dcos(teta),0.0D0)*(-b)

      Yl(1,1) = a
      Yl(1,2) = b
      Yl(1,3) = -a
      Yl(1,4) = -b
      Yl(2,1) = b
      Yl(2,2) = a
      Yl(2,3) = -b
      Yl(2,4) = -a
      Yl(3,1) = -a
      Yl(3,2) = -b
      Yl(3,3) = a
      Yl(3,4) = b
      Yl(4,1) = -b
      Yl(4,2) = -a
      Yl(4,3) = b
      Yl(4,4) = a
C
C   POSSIBILITY OF CHANGE A-B = -SIN(TETA/2.) / IM * Z0 * COS(TETA/2.)
C
C     DEBUG INIT
      END

      SUBROUTINE line2(Om,P1,L1,P2,L2,P3,L3)

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
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER L1 , L2 , L3
      DOUBLE PRECISION z0sin
      DOUBLE PRECISION Om , P1 , P2 , P3
      DOUBLE PRECISION rsh , z0 , lenth , teta , c
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      COMMON /subc  / Yl(15,15) , Vj(15)
      DOUBLE COMPLEX Yl , Vj
      DOUBLE COMPLEX a , b
C     TRANSMISSION LINE WITHOUT LOSS. PARAMETERS: R K.3., Z0, LENGTH
C         C = CONSTANT CBETA.

      rsh = P2(1)
      z0 = P3(1)
      lenth = P3(2)
      c = 300000000.0D0

      teta = Om*lenth/c
      z0sin = z0*dsin(teta)
      IF ( dabs(z0sin).LT.rsh ) b = -1.0D0/dcmplx(rsh,z0sin)
      IF ( dabs(z0sin).GE.rsh ) b = -1.0D0/dcmplx(0.0D0,z0sin)
      a = dcmplx(dcos(teta),0.0D0)*(-b)

      Yl(1,1) = a
      Yl(1,2) = b
      Yl(2,1) = b
      Yl(2,2) = a


C     DEBUG INIT
      END
