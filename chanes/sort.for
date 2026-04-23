!*==SORT.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE sort(Mn,Kn)
C
C  *****  SORTING OF THE MN ARRAY
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , i1 , j , j1 , k1 , Kn , krit , krit1 , m , m1 , mnm ,
     &        n , n1
      INTEGER Mn(2,Kn)
C
      IF ( Kn.LT.2 ) RETURN
      k1 = Kn - 1
      DO i = 1 , k1
         m = Mn(1,i)
         n = Mn(2,i)
         krit = 1024*m + n
         IF ( n.LT.0 ) krit = krit + 256
         j1 = i
         i1 = i + 1
         DO j = i1 , Kn
            m1 = Mn(1,j)
            n1 = Mn(2,j)
            krit1 = 1024*m1 + n1
            IF ( n1.LT.0 ) krit1 = krit1 + 256
            IF ( krit.GT.krit1 ) THEN
               krit = krit1
               j1 = j
            ENDIF
         ENDDO
         IF ( j1.NE.i ) THEN
            mnm = Mn(1,i)
            Mn(1,i) = Mn(1,j1)
            Mn(1,j1) = mnm
            mnm = Mn(2,i)
            Mn(2,i) = Mn(2,j1)
            Mn(2,j1) = mnm
         ENDIF
      ENDDO
      END
