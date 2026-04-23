!*==INDSV.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE indsv(Om,P1,L1,P2,L2,P3,L3)
C
C     MODEL OF THREE MUTUALLY INDUCTIVELY COUPLED
C     TWO-PORT NETWORKS (LINEAR)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION g1 , g2 , g3
      INTEGER L1 , L2 , L3 , n1 , n2
      DOUBLE PRECISION Om , P1 , P2 , P3
      DOUBLE PRECISION r1 , r2 , r3
      DOUBLE PRECISION k1 , k2 , k3
      DOUBLE PRECISION ll1 , ll2 , ll3
      DOUBLE PRECISION m1 , m2 , m3
      DOUBLE PRECISION flag1 , flag2 , flag3
      DOUBLE PRECISION za1 , za2 , za3 , zb1 , zb2 , zb3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      COMMON /subc  / Yl(15,15) , Vj(15)
      DOUBLE COMPLEX Yl , Vj
      DOUBLE COMPLEX a , b , c , d , e , f
      DOUBLE COMPLEX z1 , z2 , z3 , zm1 , zm2 , zm3 , z

C     PARAMETERS COMMON TO THE SCHEME: G1, G2, G3
C     INDIVIDUAL PARAMETERS:
C                Z1 = R1 + j*OM*LL1 - COMPLEX TOTAL
C     IMPEDANCE OF THE BRANCH CONTAINING TWO-PORT Z1.
C                Z2 = R2 + j*OM*LL2 - COMPLEX TOTAL
C     IMPEDANCE OF THE BRANCH CONTAINING TWO-PORT Z2.
C                Z3 = R3 + j*OM*LL3 - COMPLEX TOTAL
C     IMPEDANCE OF THE BRANCH CONTAINING TWO-PORT Z3.
C                M1 - MUTUAL INDUCTANCE OF TWO-PORTS
C                    Z1 AND Z2.
C                M2 - MUTUAL INDUCTANCE OF TWO-PORTS
C                    Z2 AND Z3.
C                M3 - MUTUAL INDUCTANCE OF TWO-PORTS
C                    Z3 AND Z1.
C                K1 - COUPLING COEFFICIENT OF TWO-PORTS
C                    Z1 AND Z2.
C                K2 - COUPLING COEFFICIENT OF TWO-PORTS
C                    Z2 AND Z3.
C                K3 - COUPLING COEFFICIENT OF TWO-PORTS
C                    Z3 AND Z1.
C            FLAG1 = 0. - TWO-PORTS Z1 AND Z2 CONNECTED IN PHASE
C                    1. - TWO-PORTS Z1 AND Z2 CONNECTED OUT OF PHASE
C            FLAG2 = 0. - TWO-PORTS Z2 AND Z3 CONNECTED IN PHASE
C                    1. - TWO-PORTS Z2 AND Z3 CONNECTED OUT OF PHASE
C            FLAG3 = 0. - TWO-PORTS Z3 AND Z1 CONNECTED IN PHASE
C                    1. - TWO-PORTS Z3 AND Z1 CONNECTED OUT OF PHASE
C

      g1 = P2(1)
      g2 = P2(1)
      g3 = P2(1)
      IF ( P3(1).EQ.0.0D0 ) P3(1) = 1.D-11
      r1 = P3(1)
      ll1 = P3(2)
      IF ( P3(3).EQ.0.0D0 ) P3(3) = 1.D-11
      r2 = P3(3)
      ll2 = P3(4)
      IF ( P3(5).EQ.0.0D0 ) P3(5) = 1.D-11
      r3 = P3(5)
      ll3 = P3(6)
      k1 = P3(7)
      k2 = P3(8)
      k3 = P3(9)
      flag1 = P3(10)
      flag2 = P3(11)
      flag3 = P3(12)

      m1 = k1*dsqrt(ll1*ll2)
      m2 = k2*dsqrt(ll2*ll3)
      m3 = k3*dsqrt(ll3*ll1)

      za1 = Om*ll1
      z1 = dcmplx(r1,za1)
      za2 = Om*ll2
      z2 = dcmplx(r2,za2)
      za3 = Om*ll3
      z3 = dcmplx(r3,za3)

      zb1 = Om*m1
      zb2 = Om*m2
      zb3 = Om*m3

      IF ( flag1.EQ.0.0D0 ) zm1 = dcmplx(0.0D0,zb1)
      IF ( flag1.EQ.1.0D0 ) zm1 = -dcmplx(0.0D0,zb1)
      IF ( flag1.EQ.0.0D0 .OR. flag1.EQ.1.0D0 ) THEN

         IF ( flag2.EQ.0.0D0 ) zm2 = dcmplx(0.0D0,zb2)
         IF ( flag2.EQ.1.0D0 ) zm2 = -dcmplx(0.0D0,zb2)
         IF ( flag2.EQ.0.0D0 .OR. flag2.EQ.1.0D0 ) THEN

            IF ( flag3.EQ.0.0D0 ) zm3 = dcmplx(0.0D0,zb3)
            IF ( flag3.EQ.1.0D0 ) zm3 = -dcmplx(0.0D0,zb3)
            IF ( flag3.EQ.0.0D0 .OR. flag3.EQ.1.0D0 ) THEN

               WRITE (6,10) z1 , z2 , z3 , zm1 , zm2 , zm3
 10            FORMAT (2X,'Z1 = ',E12.5,' j',E12.5/2X,'Z2 = ',E12.5,
     &                 ' j',E12.5/2X,'Z3 = ',E12.5,' j',E12.5/2X,
     &                 'ZM1 = ',E12.5,' j',E12.5/2X,'ZM2 = ',E12.5,' j',
     &                 E12.5/2X,'ZM3 = ',E12.5,' j',E12.5)

               z = z1*z2*z3 + 2*zm1*zm2*zm3 - z2*zm3**2 - z3*zm1**2 -
     &             z1*zm2**2

               a = (z2*z3-zm2**2)/z
               b = (-zm1*z3+zm2*zm3)/z
               c = (zm1*zm2-zm3*z2)/z
               d = (z1*z3-zm3**2)/z
               e = (-zm2*z1+zm1*zm3)/z
               f = (z1*z2-zm1**2)/z

               WRITE (6,15) z , a , b , c , d , e , f
 15            FORMAT (2X,'Z = ',E12.5,' j',E12.5/2X,'A = ',E12.5,' j',
     &                 E12.5/2X,'B = ',E12.5,' j',E12.5/2X,'C = ',E12.5,
     &                 ' j',E12.5/2X,'D = ',E12.5,' j',E12.5/2X,'E = ',
     &                 E12.5,' j',E12.5/2X,'F = ',E12.5,' j',E12.5)

               Yl(1,1) = a + g1
               Yl(1,2) = -a
               Yl(1,3) = b - g1
               Yl(1,4) = -b
               Yl(1,5) = c
               Yl(1,6) = -c
               Yl(2,1) = -a
               Yl(2,2) = a + g3
               Yl(2,3) = -b
               Yl(2,4) = b
               Yl(2,5) = -c
               Yl(2,6) = c - g3
               Yl(3,1) = b - g1
               Yl(3,2) = -b
               Yl(3,3) = d + g1
               Yl(3,4) = -d
               Yl(3,5) = e
               Yl(3,6) = -e
               Yl(4,1) = -b
               Yl(4,2) = b
               Yl(4,3) = -d
               Yl(4,4) = d + g2
               Yl(4,5) = -e - g2
               Yl(4,6) = e
               Yl(5,1) = c
               Yl(5,2) = -c
               Yl(5,3) = e
               Yl(5,4) = -e - g2
               Yl(5,5) = f + g2
               Yl(5,6) = -f
               Yl(6,1) = -c
               Yl(6,2) = c - g3
               Yl(6,3) = -e
               Yl(6,4) = e
               Yl(6,5) = -f
               Yl(6,6) = f + g3

               WRITE (6,20) ((n1,n2,Yl(n1,n2),n2=1,6),n1=1,6)
 20            FORMAT (2X,'INDSV',2X,'YL(',I2,',',I2,')=',E12.5,2X,
     &                 E12.5)

               RETURN
            ENDIF
         ENDIF
      ENDIF

      WRITE (6,30)
 30   FORMAT (2X,52('$')/2X,
     &        '$$                E R R O R                        $$'/2X
     &        ,'$$    IN THE INPUT DATA. 10-12 INDIVIDUAL          $$'/2
     &        X,'$$    PARAMETERS FOR THE INDUCTIVELY-COUPLED       $$'/
     &        2X,'$$    TWO-PORT NETWORKS CAN ONLY HAVE              $$'
     &        /2X,
     &        '$$    THE FOLLOWING VALUES:                        $$'/2X
     &        ,'$$    0. - TWO-PORTS CONNECTED IN PHASE.           $$'/2
     &        X,'$$    1. - TWO-PORTS CONNECTED OUT OF PHASE.       $$'/
     &        2X,52('$'))


      END
