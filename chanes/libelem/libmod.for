!*==LIBMD1.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE libmd1(Name,Ivar)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER Ivar , Kn , Knc2 , L1 , L2 , L3 , N1 , N2 , N3 , Ng , Nr
      DOUBLE PRECISION T
      COMMON /params/ P(1)
      LOGICAL Exist(5,2)
      INTEGER*4 Noi(5,2) , Nou(5,2)
      INTEGER*4 Koi , Kouv(5) , Kopv(5) , Nr1v(5) , Nb1v(5)
      DOUBLE PRECISION B1(1) , Om , P
      DOUBLE COMPLEX Val(1) , Dval(1)
      CHARACTER*4 Name(4) , names(17)
      DATA names/'GN  ' , 'CN  ' , 'VD  ' , 'ICU ' , 'FET ' , 'LIN ' ,
     &     'POLY' , 'JUNC' , 'BARR' , 'DIFF' , 'SCHT' , 'SP  ' ,
     &     'CURT' , 'TAJ ' , 'BIP ' , 'TR  ' , 'CUSD'/

C     WRITE(6,100) NAME,NAMES(11)
C 100 FORMAT(2X,'LIBMD1 NAME =',4A4,2X,4A4)

      IF ( Name(1).EQ.names(1) ) THEN
         IF ( Name(2).EQ.names(6) ) CALL lin1(Ivar)
         IF ( Name(2).EQ.names(7) ) CALL poly51(Ivar)
         IF ( Name(2).EQ.names(8) ) CALL junc1(Ivar)
      ENDIF

      IF ( Name(1).EQ.names(3) ) THEN
         IF ( Name(2).EQ.names(11) ) CALL mdsch1(Ivar)
         IF ( Name(2).EQ.names(8) ) CALL junc1(Ivar)
      ENDIF

      IF ( Name(1).EQ.names(5) ) THEN
         IF ( Name(2).EQ.names(13) ) CALL curt1(Ivar)
      ENDIF
C     IF(NAME(2).EQ.NAMES(14)) CALL TAJ1  (IVAR)

      IF ( Name(1).EQ.names(2) ) THEN
         IF ( Name(2).EQ.names(6) ) CALL clin1(Ivar)
         IF ( Name(2).EQ.names(7) ) CALL cpoly1(Ivar)
         IF ( Name(2).EQ.names(9) ) CALL cbarr1(Ivar)
         IF ( Name(2).EQ.names(10) ) CALL cdiff1(Ivar)
      ENDIF

      IF ( Name(1).EQ.names(4) ) THEN
C     IF(NAME(2).EQ.NAMES( 6)) CALL ICULN1(IVAR)
         IF ( Name(2).EQ.names(7) ) CALL icupl1(Ivar)
C     IF(NAME(2).EQ.NAMES(12)) CALL ICUSP1(IVAR)
         IF ( Name(2).EQ.names(08) ) CALL icuj1(Ivar)
      ENDIF
C      IF(NAME(2).EQ.NAMES(17)) CALL CUSD1 (IVAR)

      IF ( Name(1).EQ.names(15) ) THEN
         IF ( Name(2).EQ.names(16) ) CALL biptr1(Ivar)
      ENDIF

C  SCHT = SCHOTTKY DIODE
C  JUNC = DIODE WITHOUT NONLINEAR CAPACITANCE
      RETURN
C2222222222222222222222222222222222222222222222222222222222222222222222
      ENTRY libmd2(Name,Om,N1,L1,N2,L2,N3,L3)
C     WRITE(6,1000) N1,L1,N2,L2,N3,L3
C1000 FORMAT(2X,'LIBMD2: N1=',I3,' L1=',I3,
C    *         ' N2=',I3,' L2=',I3,
C    *         ' N3=',I5,' L3=',I3)

      IF ( Name(1).EQ.names(1) ) THEN
         IF ( Name(2).EQ.names(6) )
     &        CALL lin2(Om,P(N1),L1,P(N2),L2,P(N3),L3)
         IF ( Name(2).EQ.names(7) )
     &        CALL poly52(Om,P(N1),L1,P(N2),L2,P(N3),L3)
         IF ( Name(2).EQ.names(8) )
     &        CALL junc2(Om,P(N1),L1,P(N2),L2,P(N3),L3)
      ENDIF

      IF ( Name(1).EQ.names(3) ) THEN
         IF ( Name(2).EQ.names(11) )
     &        CALL mdsch2(Om,P(N1),L1,P(N2),L2,P(N3),L3)
         IF ( Name(2).EQ.names(8) )
     &        CALL junc2(Om,P(N1),L1,P(N2),L2,P(N3),L3)
      ENDIF

      IF ( Name(1).EQ.names(5) ) THEN
         IF ( Name(2).EQ.names(13) )
     &        CALL curt2(Om,P(N1),L1,P(N2),L2,P(N3),L3)
      ENDIF
C     IF(NAME(2).EQ.NAMES(14))CALL TAJ2  (OM,P(N1),L1,P(N2),
C    *                                           L2,P(N3),L3)

      IF ( Name(1).EQ.names(2) ) THEN
         IF ( Name(2).EQ.names(6) )
     &        CALL clin2(Om,P(N1),L1,P(N2),L2,P(N3),L3)
         IF ( Name(2).EQ.names(7) )
     &        CALL cpoly2(Om,P(N1),L1,P(N2),L2,P(N3),L3)
         IF ( Name(2).EQ.names(9) )
     &        CALL cbarr2(Om,P(N1),L1,P(N2),L2,P(N3),L3)
         IF ( Name(2).EQ.names(10) )
     &        CALL cdiff2(Om,P(N1),L1,P(N2),L2,P(N3),L3)
      ENDIF

      IF ( Name(1).EQ.names(4) ) THEN
C     IF(NAME(2).EQ.NAMES( 6))CALL ICULN2(OM,P(N1),L1,P(N2),
C    *                                           L2,P(N3),L3)
         IF ( Name(2).EQ.names(7) )
     &        CALL icupl2(Om,P(N1),L1,P(N2),L2,P(N3),L3)
C     IF(NAME(2).EQ.NAMES(12))CALL ICUSP2(OM,P(N1),L1,P(N2),
C    *                                           L2,P(N3),L3)
         IF ( Name(2).EQ.names(8) )
     &        CALL icuj2(Om,P(N1),L1,P(N2),L2,P(N3),L3)
      ENDIF
C      IF(NAME(2).EQ.NAMES(17))CALL CUSD2 (OM,P(N1),L1,P(N2),
C     +                                           L2,P(N3),L3)

      IF ( Name(1).EQ.names(15) ) THEN
         IF ( Name(2).EQ.names(16) )
     &        CALL biptr2(Om,P(N1),L1,P(N2),L2,P(N3),L3)
      ENDIF


      RETURN
C3333333333333333333333333333333333333333333333333333333333333333333333
      ENTRY libmd3(Name,Ng,N1,L1,N2,L2,N3,L3,B1,Knc2,Nr,*)
C     WRITE(6,1001) N1,L1,N2,L2,N3,L3
C1001 FORMAT(2X,'LIBMD3: N1=',I3,' L1=',I3,' N2=',I3,' L2=',I3,
C    *                 ' N3=',I5,' L3=',I3)
      IF ( Name(1).EQ.names(1) ) THEN
         IF ( Name(2).EQ.names(6) )
     &        CALL lin3(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
         IF ( Name(2).EQ.names(7) )
     &        CALL poly53(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
         IF ( Name(2).EQ.names(8) )
     &        CALL junc3(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
      ENDIF

      IF ( Name(1).EQ.names(3) ) THEN

         IF ( Name(2).EQ.names(11) )
     &        CALL mdsch3(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)

         IF ( Name(2).EQ.names(8) )
     &        CALL junc3(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
      ENDIF

      IF ( Name(1).EQ.names(5) ) THEN
         IF ( Name(2).EQ.names(13) )
     &        CALL curt3(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
      ENDIF
C     IF(NAME(2).EQ.NAMES(14))CALL TAJ3  (NG,P(N1),L1,P(N2),L2,P(N3),
C    *                                           L3,B1,KNC2,NR,*300)

      IF ( Name(1).EQ.names(2) ) THEN
         IF ( Name(2).EQ.names(6) )
     &        CALL clin3(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
         IF ( Name(2).EQ.names(7) )
     &        CALL cpoly3(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
         IF ( Name(2).EQ.names(9) )
     &        CALL cbarr3(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
         IF ( Name(2).EQ.names(10) )
     &        CALL cdiff3(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
      ENDIF

      IF ( Name(1).EQ.names(4) ) THEN
C     IF(NAME(2).EQ.NAMES( 6))CALL ICULN3(NG,P(N1),L1,P(N2),L2,P(N3),
C    *                                           L3,B1,KNC2,NR,*300)
         IF ( Name(2).EQ.names(7) )
     &        CALL icupl3(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
C     IF(NAME(2).EQ.NAMES(12))CALL ICUSP3(NG,P(N1),L1,P(N2),L2,P(N3),
C    *                                           L3,B1,KNC2,NR,*300)
         IF ( Name(2).EQ.names(8) )
     &        CALL icuj3(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
      ENDIF
C      IF(NAME(2).EQ.NAMES(17))CALL CUSD3 (NG,P(N1),L1,P(N2),L2,P(N3),
C     +                                           L3,B1,KNC2,NR,*300)

      IF ( Name(1).EQ.names(15) ) THEN
         IF ( Name(2).EQ.names(16) )
     &        CALL biptr3(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
      ENDIF

      RETURN
 300  RETURN 1
C4444444444444444444444444444444444444444444444444444444444444444444444
      ENTRY libmd4(Name,Ng,N1,L1,N2,L2,N3,L3,B1,Knc2,Nr,*)
      IF ( Name(1).EQ.names(1) ) THEN
         IF ( Name(2).EQ.names(6) )
     &        CALL lin4(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
         IF ( Name(2).EQ.names(7) )
     &        CALL poly54(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
         IF ( Name(2).EQ.names(8) )
     &        CALL junc4(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
      ENDIF

      IF ( Name(1).EQ.names(3) ) THEN
         IF ( Name(2).EQ.names(11) )
     &        CALL mdsch4(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
         IF ( Name(2).EQ.names(8) )
     &        CALL junc4(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
      ENDIF

      IF ( Name(1).EQ.names(5) ) THEN
         IF ( Name(2).EQ.names(13) )
     &        CALL curt4(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
      ENDIF
C     IF(NAME(2).EQ.NAMES(14))CALL TAJ4  (NG,P(N1),L1,P(N2),L2,P(N3),
C    *                                           L3,B1,KNC2,NR,*300)

      IF ( Name(1).EQ.names(2) ) THEN
         IF ( Name(2).EQ.names(6) )
     &        CALL clin4(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
         IF ( Name(2).EQ.names(7) )
     &        CALL cpoly4(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
         IF ( Name(2).EQ.names(9) )
     &        CALL cbarr4(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
         IF ( Name(2).EQ.names(10) )
     &        CALL cdiff4(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
      ENDIF

      IF ( Name(1).EQ.names(4) ) THEN
C     IF(NAME(2).EQ.NAMES( 6))CALL ICULN4(NG,P(N1),L1,P(N2),L2,P(N3),
C    *                                           L3,B1,KNC2,NR,*300)
         IF ( Name(2).EQ.names(7) )
     &        CALL icupl4(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
C     IF(NAME(2).EQ.NAMES(12))CALL ICUSP4(NG,P(N1),L1,P(N2),L2,P(N3),
C    *                                           L3,B1,KNC2,NR,*300)
         IF ( Name(2).EQ.names(8) )
     &        CALL icuj4(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
      ENDIF
C      IF(NAME(2).EQ.NAMES(17))CALL CUSD4 (NG,P(N1),L1,P(N2),L2,P(N3),
C     +                                           L3,B1,KNC2,NR,*300)

      IF ( Name(1).EQ.names(15) ) THEN
         IF ( Name(2).EQ.names(16) )
     &        CALL biptr4(Ng,P(N1),L1,P(N2),L2,P(N3),L3,B1,Knc2,Nr,*300)
      ENDIF
      RETURN
C5555555555555555555555555555555555555555555555555555555555555555555555
      ENTRY libmd5(Name,Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
      IF ( Name(1).EQ.names(1) ) THEN
         IF ( Name(2).EQ.names(6) )
     &        CALL lin5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
         IF ( Name(2).EQ.names(7) ) CALL poly55(Noi,Nou,Exist,Koi,Kouv,
     &        Kopv,Nr1v,Nb1v)
         IF ( Name(2).EQ.names(8) )
     &        CALL junc5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
      ENDIF

      IF ( Name(1).EQ.names(3) ) THEN
         IF ( Name(2).EQ.names(11) ) CALL mdsch5(Noi,Nou,Exist,Koi,Kouv,
     &        Kopv,Nr1v,Nb1v)
         IF ( Name(2).EQ.names(8) )
     &        CALL junc5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
      ENDIF

      IF ( Name(1).EQ.names(5) ) THEN
         IF ( Name(2).EQ.names(13) )
     &        CALL curt5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
      ENDIF
C     IF(NAME(2).EQ.NAMES(14))CALL TAJ5
C    *                 (NOI,NOU,EXIST,KOI,KOUV,KOPV,NR1V,NB1V)

      IF ( Name(1).EQ.names(2) ) THEN
         IF ( Name(2).EQ.names(6) )
     &        CALL clin5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
         IF ( Name(2).EQ.names(7) ) CALL cpoly5(Noi,Nou,Exist,Koi,Kouv,
     &        Kopv,Nr1v,Nb1v)
         IF ( Name(2).EQ.names(9) ) CALL cbarr5(Noi,Nou,Exist,Koi,Kouv,
     &        Kopv,Nr1v,Nb1v)
         IF ( Name(2).EQ.names(10) ) CALL cdiff5(Noi,Nou,Exist,Koi,Kouv,
     &        Kopv,Nr1v,Nb1v)
      ENDIF

      IF ( Name(1).EQ.names(4) ) THEN
C     IF(NAME(2).EQ.NAMES( 6))CALL ICULN5
C    *                 (NOI,NOU,EXIST,KOI,KOUV,KOPV,NR1V,NB1V)
         IF ( Name(2).EQ.names(7) ) CALL icupl5(Noi,Nou,Exist,Koi,Kouv,
     &        Kopv,Nr1v,Nb1v)
C     IF(NAME(2).EQ.NAMES(12))CALL ICUSP5
C    *                 (NOI,NOU,EXIST,KOI,KOUV,KOPV,NR1V,NB1V)
         IF ( Name(2).EQ.names(8) )
     &        CALL icuj5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
      ENDIF
C      IF(NAME(2).EQ.NAMES(17))CALL CUSD5
C     +                 (NOI,NOU,EXIST,KOI,KOUV,KOPV,NR1V,NB1V)

      IF ( Name(1).EQ.names(15) ) THEN
         IF ( Name(2).EQ.names(16) ) CALL biptr5(Noi,Nou,Exist,Koi,Kouv,
     &        Kopv,Nr1v,Nb1v)
      ENDIF
      RETURN
C6666666666666666666666666666666666666666666666666666666666666666666666
      ENTRY libmd6(Name,Ng,N1,L1,N2,L2,N3,L3,Val,Dval,Kn,Nr,T)

      IF ( Name(1).EQ.names(1) ) THEN
         IF ( Name(2).EQ.names(8) )
     &        CALL junc6(Ng,P(N1),L1,P(N2),L2,P(N3),L3,Val,Dval,Kn,Nr,T)
      ENDIF

      IF ( Name(1).EQ.names(3) ) THEN
         IF ( Name(2).EQ.names(11) )
     &        CALL mdsch6(Ng,P(N1),L1,P(N2),L2,P(N3),L3,Val,Dval,Kn,Nr,
     &        T)
         IF ( Name(2).EQ.names(8) )
     &        CALL junc6(Ng,P(N1),L1,P(N2),L2,P(N3),L3,Val,Dval,Kn,Nr,T)
      ENDIF

      IF ( Name(1).EQ.names(5) ) THEN
         IF ( Name(2).EQ.names(13) )
     &        CALL curt6(Ng,P(N1),L1,P(N2),L2,P(N3),L3,Val,Dval,Kn,Nr,T)
         IF ( Name(2).EQ.names(14) )
     &        CALL junc6(Ng,P(N1),L1,P(N2),L2,P(N3),L3,Val,Dval,Kn,Nr,T)
      ENDIF

      IF ( Name(1).EQ.names(2) ) THEN
         IF ( Name(2).EQ.names(10) )
     &        CALL cdiff6(Ng,P(N1),L1,P(N2),L2,P(N3),L3,Val,Dval,Kn,Nr,
     &        T)
      ENDIF

      IF ( Name(1).EQ.names(4) ) THEN
         IF ( Name(2).EQ.names(8) )
     &        CALL icuj6(Ng,P(N1),L1,P(N2),L2,P(N3),L3,Val,Dval,Kn,Nr,T)
      ENDIF
      END
