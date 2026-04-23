c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE LIBMD1(NAME,IVAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/PARAMS/    P(1)
      LOGICAL           EXIST(5,2)
      INTEGER*4         NOI(5,2),NOU(5,2)
      INTEGER*4         KOI,KOUV(5),KOPV(5),NR1V(5),NB1V(5)
      DOUBLE PRECISION              B1(1),OM,P
      DOUBLE COMPLEX           VAL(1),DVAL(1)
      CHARACTER*4       NAME(4),NAMES(17)
      DATA NAMES          /'GN  ','CN  ','VD  ','ICU ','FET ',
     + 'LIN ','POLY','JUNC','BARR','DIFF',           'SCHT','SP  ','CURT
     +','TAJ ','BIP ',           'TR  ','CUSD'/

C     WRITE(6,100) NAME,NAMES(11)
C 100 FORMAT(2X,'LIBMD1 NAME =',4A4,2X,4A4)

      IF(NAME(1).NE.NAMES(1))GOTO 110
      IF(NAME(2).EQ.NAMES(6)) CALL LIN1  (IVAR)
      IF(NAME(2).EQ.NAMES(7)) CALL POLY51(IVAR)
      IF(NAME(2).EQ.NAMES(8)) CALL JUNC1 (IVAR)
  110 CONTINUE

      IF(NAME(1).NE.NAMES(3))GOTO 120
      IF(NAME(2).EQ.NAMES(11)) CALL MDSCH1(IVAR)
      IF(NAME(2).EQ.NAMES( 8)) CALL JUNC1(IVAR)
  120 CONTINUE

      IF(NAME(1).NE.NAMES(5))GOTO 130
      IF(NAME(2).EQ.NAMES(13)) CALL CURT1 (IVAR)
C     IF(NAME(2).EQ.NAMES(14)) CALL TAJ1  (IVAR)
  130 CONTINUE

      IF(NAME(1).NE.NAMES(2))GOTO 140
      IF(NAME(2).EQ.NAMES( 6)) CALL CLIN1 (IVAR)
      IF(NAME(2).EQ.NAMES( 7)) CALL CPOLY1(IVAR)
      IF(NAME(2).EQ.NAMES( 9)) CALL CBARR1(IVAR)
      IF(NAME(2).EQ.NAMES(10)) CALL CDIFF1(IVAR)
  140 CONTINUE

      IF(NAME(1).NE.NAMES(4))GOTO 150
C     IF(NAME(2).EQ.NAMES( 6)) CALL ICULN1(IVAR)
      IF(NAME(2).EQ.NAMES( 7)) CALL ICUPL1(IVAR)
C     IF(NAME(2).EQ.NAMES(12)) CALL ICUSP1(IVAR)
      IF(NAME(2).EQ.NAMES(08)) CALL ICUJ1 (IVAR)
C      IF(NAME(2).EQ.NAMES(17)) CALL CUSD1 (IVAR)
  150 CONTINUE

      IF(NAME(1).NE.NAMES(15)) GOTO 160
      IF(NAME(2).EQ.NAMES(16)) CALL BIPTR1(IVAR)
  160 CONTINUE

C  SCHT = SCHOTTKY DIODE
C  JUNC = DIODE WITHOUT NONLINEAR CAPACITANCE
      RETURN
C2222222222222222222222222222222222222222222222222222222222222222222222
      ENTRY LIBMD2(NAME,OM,N1,L1,N2,L2,N3,L3)
C     WRITE(6,1000) N1,L1,N2,L2,N3,L3
C1000 FORMAT(2X,'LIBMD2: N1=',I3,' L1=',I3,
C    *         ' N2=',I3,' L2=',I3,
C    *         ' N3=',I5,' L3=',I3)

      IF(NAME(1).NE.NAMES(1))GOTO 210
      IF(NAME(2).EQ.NAMES(6))CALL LIN2  (OM,P(N1),L1,P(N2),
     +                              L2,P(N3),L3)
      IF(NAME(2).EQ.NAMES(7))CALL POLY52(OM,P(N1),L1,P(N2),
     +                              L2,P(N3),L3)
      IF(NAME(2).EQ.NAMES(8))CALL JUNC2 (OM,P(N1),L1,P(N2),
     +                              L2,P(N3),L3)
 210  CONTINUE

      IF(NAME(1).NE.NAMES(3))GOTO 220
      IF(NAME(2).EQ.NAMES(11))CALL MDSCH2(OM,P(N1),L1,P(N2),
     +                               L2,P(N3),L3)
      IF(NAME(2).EQ.NAMES( 8))CALL JUNC2(OM,P(N1),L1,P(N2),
     +                              L2,P(N3),L3)
 220  CONTINUE

      IF(NAME(1).NE.NAMES(5))GOTO 230
      IF(NAME(2).EQ.NAMES(13))CALL CURT2 (OM,P(N1),L1,P(N2),
     +                               L2,P(N3),L3)
C     IF(NAME(2).EQ.NAMES(14))CALL TAJ2  (OM,P(N1),L1,P(N2),
C    *                                           L2,P(N3),L3)
 230  CONTINUE

      IF(NAME(1).NE.NAMES(2))GOTO 240
      IF(NAME(2).EQ.NAMES( 6))CALL CLIN2 (OM,P(N1),L1,P(N2),
     +                               L2,P(N3),L3)
      IF(NAME(2).EQ.NAMES( 7))CALL CPOLY2(OM,P(N1),L1,P(N2),
     +                               L2,P(N3),L3)
      IF(NAME(2).EQ.NAMES( 9))CALL CBARR2(OM,P(N1),L1,P(N2),
     +                               L2,P(N3),L3)
      IF(NAME(2).EQ.NAMES(10))CALL CDIFF2(OM,P(N1),L1,P(N2),
     +                               L2,P(N3),L3)
 240  CONTINUE

      IF(NAME(1).NE.NAMES(4))GOTO 250
C     IF(NAME(2).EQ.NAMES( 6))CALL ICULN2(OM,P(N1),L1,P(N2),
C    *                                           L2,P(N3),L3)
      IF(NAME(2).EQ.NAMES( 7))CALL ICUPL2(OM,P(N1),L1,P(N2),
     +                               L2,P(N3),L3)
C     IF(NAME(2).EQ.NAMES(12))CALL ICUSP2(OM,P(N1),L1,P(N2),
C    *                                           L2,P(N3),L3)
      IF(NAME(2).EQ.NAMES( 8))CALL ICUJ2 (OM,P(N1),L1,P(N2),
     +                               L2,P(N3),L3)
C      IF(NAME(2).EQ.NAMES(17))CALL CUSD2 (OM,P(N1),L1,P(N2),
C     +                                           L2,P(N3),L3)
 250  CONTINUE

      IF(NAME(1).NE.NAMES(15)) GOTO 260
      IF(NAME(2).EQ.NAMES(16)) CALL BIPTR2(OM,P(N1),L1,P(N2),
     +                                 L2,P(N3),L3)
 260  CONTINUE


      RETURN
C3333333333333333333333333333333333333333333333333333333333333333333333
      ENTRY LIBMD3(NAME,NG,N1,L1,N2,L2,N3,L3,B1,KNC2,NR,*)
C     WRITE(6,1001) N1,L1,N2,L2,N3,L3
C1001 FORMAT(2X,'LIBMD3: N1=',I3,' L1=',I3,' N2=',I3,' L2=',I3,
C    *                 ' N3=',I5,' L3=',I3)
      IF(NAME(1).NE.NAMES(1))GOTO 310
      IF(NAME(2).EQ.NAMES(6))CALL LIN3  (NG,P(N1),L1,P(N2),L2,P(N3),L3,
     +                                             B1,KNC2,NR,*300)
      IF(NAME(2).EQ.NAMES(7))CALL POLY53(NG,P(N1),L1,P(N2),L2,P(N3),L3,
     +                                             B1,KNC2,NR,*300)
      IF(NAME(2).EQ.NAMES(8))CALL JUNC3 (NG,P(N1),L1,P(N2),L2,P(N3),L3,
     +                                             B1,KNC2,NR,*300)
  310 CONTINUE

      IF(NAME(1).NE.NAMES(3))GOTO 320

      IF(NAME(2).EQ.NAMES(11))CALL MDSCH3(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)

      IF(NAME(2).EQ.NAMES( 8))CALL JUNC3(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                       L3,B1,KNC2,NR,*300)
 320  CONTINUE

      IF(NAME(1).NE.NAMES(5))GOTO 330
      IF(NAME(2).EQ.NAMES(13))CALL CURT3 (NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
C     IF(NAME(2).EQ.NAMES(14))CALL TAJ3  (NG,P(N1),L1,P(N2),L2,P(N3),
C    *                                           L3,B1,KNC2,NR,*300)
 330  CONTINUE

      IF(NAME(1).NE.NAMES(2))GOTO 340
      IF(NAME(2).EQ.NAMES( 6))CALL CLIN3 (NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
      IF(NAME(2).EQ.NAMES( 7))CALL CPOLY3(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
      IF(NAME(2).EQ.NAMES( 9))CALL CBARR3(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
      IF(NAME(2).EQ.NAMES(10))CALL CDIFF3(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
  340 CONTINUE

      IF(NAME(1).NE.NAMES(4))GOTO 350
C     IF(NAME(2).EQ.NAMES( 6))CALL ICULN3(NG,P(N1),L1,P(N2),L2,P(N3),
C    *                                           L3,B1,KNC2,NR,*300)
      IF(NAME(2).EQ.NAMES( 7))CALL ICUPL3(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
C     IF(NAME(2).EQ.NAMES(12))CALL ICUSP3(NG,P(N1),L1,P(N2),L2,P(N3),
C    *                                           L3,B1,KNC2,NR,*300)
      IF(NAME(2).EQ.NAMES( 8))CALL ICUJ3 (NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
C      IF(NAME(2).EQ.NAMES(17))CALL CUSD3 (NG,P(N1),L1,P(N2),L2,P(N3),
C     +                                           L3,B1,KNC2,NR,*300)
  350 CONTINUE

      IF(NAME(1).NE.NAMES(15)) GOTO 360
      IF(NAME(2).EQ.NAMES(16)) CALL BIPTR3(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                          L3,B1,KNC2,NR,*300)
  360 CONTINUE

      RETURN
  300 RETURN 1
C4444444444444444444444444444444444444444444444444444444444444444444444
      ENTRY LIBMD4(NAME,NG,N1,L1,N2,L2,N3,L3,B1,KNC2,NR,*)
      IF(NAME(1).NE.NAMES(1))GOTO 410
      IF(NAME(2).EQ.NAMES(6))CALL LIN4  (NG,P(N1),L1,P(N2),L2,P(N3),L3,
     +                                             B1,KNC2,NR,*300)
      IF(NAME(2).EQ.NAMES(7))CALL POLY54(NG,P(N1),L1,P(N2),L2,P(N3),L3,
     +                                             B1,KNC2,NR,*300)
      IF(NAME(2).EQ.NAMES(8))CALL JUNC4 (NG,P(N1),L1,P(N2),L2,P(N3),L3,
     +                                             B1,KNC2,NR,*300)
  410 CONTINUE

      IF(NAME(1).NE.NAMES(3))GOTO 420
      IF(NAME(2).EQ.NAMES(11))CALL MDSCH4(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
      IF(NAME(2).EQ.NAMES( 8))CALL JUNC4(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                       L3,B1,KNC2,NR,*300)
 420  CONTINUE

      IF(NAME(1).NE.NAMES(5))GOTO 430
      IF(NAME(2).EQ.NAMES(13))CALL CURT4 (NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
C     IF(NAME(2).EQ.NAMES(14))CALL TAJ4  (NG,P(N1),L1,P(N2),L2,P(N3),
C    *                                           L3,B1,KNC2,NR,*300)
 430  CONTINUE

      IF(NAME(1).NE.NAMES(2))GOTO 440
      IF(NAME(2).EQ.NAMES( 6))CALL CLIN4 (NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
      IF(NAME(2).EQ.NAMES( 7))CALL CPOLY4(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
      IF(NAME(2).EQ.NAMES( 9))CALL CBARR4(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
      IF(NAME(2).EQ.NAMES(10))CALL CDIFF4(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
 440  CONTINUE

      IF(NAME(1).NE.NAMES(4))GOTO 450
C     IF(NAME(2).EQ.NAMES( 6))CALL ICULN4(NG,P(N1),L1,P(N2),L2,P(N3),
C    *                                           L3,B1,KNC2,NR,*300)
      IF(NAME(2).EQ.NAMES( 7))CALL ICUPL4(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
C     IF(NAME(2).EQ.NAMES(12))CALL ICUSP4(NG,P(N1),L1,P(N2),L2,P(N3),
C    *                                           L3,B1,KNC2,NR,*300)
      IF(NAME(2).EQ.NAMES( 8))CALL ICUJ4 (NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
C      IF(NAME(2).EQ.NAMES(17))CALL CUSD4 (NG,P(N1),L1,P(N2),L2,P(N3),
C     +                                           L3,B1,KNC2,NR,*300)
 450  CONTINUE

      IF(NAME(1).NE.NAMES(15))GOTO 460
      IF(NAME(2).EQ.NAMES(16))CALL BIPTR4(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                        L3,B1,KNC2,NR,*300)
 460  CONTINUE
      RETURN
C5555555555555555555555555555555555555555555555555555555555555555555555
      ENTRY LIBMD5(NAME,NOI,NOU,EXIST,KOI,KOUV,KOPV,NR1V,NB1V)
      IF(NAME(1).NE.NAMES(1))GOTO 510
      IF(NAME(2).EQ.NAMES(6))CALL LIN5                 (NOI,NOU,EXIST,KO
     +I,KOUV,KOPV,NR1V,NB1V)
      IF(NAME(2).EQ.NAMES(7))CALL POLY55                 (NOI,NOU,EXIST,
     +KOI,KOUV,KOPV,NR1V,NB1V)
      IF(NAME(2).EQ.NAMES(8))CALL JUNC5                 (NOI,NOU,EXIST,K
     +OI,KOUV,KOPV,NR1V,NB1V)
 510  CONTINUE

      IF(NAME(1).NE.NAMES(3))GOTO 520
      IF(NAME(2).EQ.NAMES(11))CALL MDSCH5                 (NOI,NOU,EXIST
     +,KOI,KOUV,KOPV,NR1V,NB1V)
      IF(NAME(2).EQ.NAMES( 8))CALL JUNC5                 (NOI,NOU,EXIST,
     +KOI,KOUV,KOPV,NR1V,NB1V)
 520  CONTINUE

      IF(NAME(1).NE.NAMES(5))GOTO 530
      IF(NAME(2).EQ.NAMES(13))CALL CURT5                 (NOI,NOU,EXIST,
     +KOI,KOUV,KOPV,NR1V,NB1V)
C     IF(NAME(2).EQ.NAMES(14))CALL TAJ5
C    *                 (NOI,NOU,EXIST,KOI,KOUV,KOPV,NR1V,NB1V)
 530  CONTINUE

      IF(NAME(1).NE.NAMES(2))GOTO 540
      IF(NAME(2).EQ.NAMES( 6))CALL CLIN5                 (NOI,NOU,EXIST,
     +KOI,KOUV,KOPV,NR1V,NB1V)
      IF(NAME(2).EQ.NAMES( 7))CALL CPOLY5                 (NOI,NOU,EXIST
     +,KOI,KOUV,KOPV,NR1V,NB1V)
      IF(NAME(2).EQ.NAMES( 9))CALL CBARR5                 (NOI,NOU,EXIST
     +,KOI,KOUV,KOPV,NR1V,NB1V)
      IF(NAME(2).EQ.NAMES(10))CALL CDIFF5                 (NOI,NOU,EXIST
     +,KOI,KOUV,KOPV,NR1V,NB1V)
  540 CONTINUE

      IF(NAME(1).NE.NAMES(4))GOTO 550
C     IF(NAME(2).EQ.NAMES( 6))CALL ICULN5
C    *                 (NOI,NOU,EXIST,KOI,KOUV,KOPV,NR1V,NB1V)
      IF(NAME(2).EQ.NAMES( 7))CALL ICUPL5                 (NOI,NOU,EXIST
     +,KOI,KOUV,KOPV,NR1V,NB1V)
C     IF(NAME(2).EQ.NAMES(12))CALL ICUSP5
C    *                 (NOI,NOU,EXIST,KOI,KOUV,KOPV,NR1V,NB1V)
      IF(NAME(2).EQ.NAMES( 8))CALL ICUJ5                 (NOI,NOU,EXIST,
     +KOI,KOUV,KOPV,NR1V,NB1V)
C      IF(NAME(2).EQ.NAMES(17))CALL CUSD5
C     +                 (NOI,NOU,EXIST,KOI,KOUV,KOPV,NR1V,NB1V)
  550 CONTINUE

      IF(NAME(1).NE.NAMES(15))GOTO 560
      IF(NAME(2).EQ.NAMES(16))CALL BIPTR5                 (NOI,NOU,EXIST
     +,KOI,KOUV,KOPV,NR1V,NB1V)
 560  CONTINUE
      RETURN
C6666666666666666666666666666666666666666666666666666666666666666666666
      ENTRY LIBMD6(NAME,NG,N1,L1,N2,L2,N3,L3,VAL,DVAL,KN,NR,T)

      IF(NAME(1).NE.NAMES(1))GOTO 610
      IF(NAME(2).EQ.NAMES(8))CALL JUNC6(NG,P(N1),L1,P(N2),L2,P(N3),
     +                             L3,VAL,DVAL,KN,NR,T)
 610  CONTINUE

      IF(NAME(1).NE.NAMES(3))GOTO 620
      IF(NAME(2).EQ.NAMES(11))CALL MDSCH6(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                 L3,VAL,DVAL,KN,NR,T)
      IF(NAME(2).EQ.NAMES( 8))CALL  JUNC6(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                 L3,VAL,DVAL,KN,NR,T)
 620  CONTINUE

      IF(NAME(1).NE.NAMES(5))GOTO 630
      IF(NAME(2).EQ.NAMES(13))CALL CURT6 (NG,P(N1),L1,P(N2),L2,P(N3),
     +                                 L3,VAL,DVAL,KN,NR,T)
      IF(NAME(2).EQ.NAMES(14))CALL JUNC6 (NG,P(N1),L1,P(N2),L2,P(N3),
     +                                 L3,VAL,DVAL,KN,NR,T)
  630 CONTINUE

      IF(NAME(1).NE.NAMES(2))GOTO 650
      IF(NAME(2).EQ.NAMES(10))CALL CDIFF6(NG,P(N1),L1,P(N2),L2,P(N3),
     +                                 L3,VAL,DVAL,KN,NR,T)
  650 CONTINUE

      IF(NAME(1).NE.NAMES(4))GOTO 640
      IF(NAME(2).EQ.NAMES(8))CALL ICUJ6 (NG,P(N1),L1,P(N2),L2,P(N3),L3,
     +                                  VAL,DVAL,KN,NR,T)
  640 RETURN
      END
