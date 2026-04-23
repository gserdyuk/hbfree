!*==SUMDIF.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE sumdif(Q,Ir,Is,Sum,Diff)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER Ir , Is , izr , izs
      COMMON /blw1  / W , W1
      DOUBLE PRECISION W(20) , W1(200)
      COMMON /blw2  / Wr , Ws
      INTEGER Wr(20,20) , Ws(20,20)
      DOUBLE COMPLEX Q(1) , Sum , Diff , qizs
C
C
      izs = Ws(Ir,Is)
      qizs = Q(izs)
      izr = Wr(Ir,Is)
      IF ( izr.LT.0 ) THEN
         izr = -izr
         Sum = dconjg(Q(izr))
      ELSE
         Sum = Q(izr)
      ENDIF
      Diff = Sum - qizs
      IF ( W(Is).NE.0.D0 ) Sum = Sum + qizs
C     DEBUG INIT(Q,IR,IS,SUM,DIFF),SUBTRACE
      END
