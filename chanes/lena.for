!*==LENA.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE lena
C ***
C *** THE PROGRAM CONTAINS THE ARRAY "A"
C ***
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER Iff , Ka , Kiff , Knniff , Nniff
      DOUBLE PRECISION Pniff
      COMMON /blka  / A , Ka
      COMMON /bliff / Iff(7,4) , Kiff , Nniff(4,8) , Knniff , Pniff(8) ,
     &                Fne1 , Fne2
      CHARACTER*4 Fne1 , Fne2
C
C
C      ARRAY "A" IS A PREPARATION FOR THE FORM "MPOINT"
      INTEGER*4 A(20,50)
C                         KA - MAX. NUMBER OF TYPES OF ELEMENTS
      DATA A/'R   ' , '    ' , '    ' , '    ' , 0 , 1 , 5 , 2 , 0 , 0 ,
     &     1 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 'L   ' , '    ' ,
     &     '    ' , '    ' , 0 , 1 , 5 , 2 , 0 , 0 , 1 , 0 , 1 , 0 , 0 ,
     &     0 , 0 , 0 , 0 , 0 , 'C   ' , '    ' , '    ' , '    ' , 0 ,
     &     1 , 5 , 2 , 0 , 0 , 1 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     &     'E   ' , '    ' , '    ' , '    ' , 0 , 2 , 10 , 2 , 0 , 0 ,
     &     1 , 0 , 4 , 0 , 0 , 0 , 0 , 0 , 0 , 2 , 'J   ' , '    ' ,
     &     '    ' , '    ' , 0 , 2 , 10 , 2 , 0 , 0 , 1 , 0 , 4 , 0 ,
     &     0 , 0 , 0 , 0 , 0 , 2 , 'P2  ' , '    ' , '    ' , '    ' ,
     &     0 , 2 , 10 , 2 , 0 , 0 , 1 , 0 , 7 , 0 , 0 , 0 , 0 , 0 , 0 ,
     &     2 , 'VD  ' , 'SCHT' , '    ' , '    ' , 0 , 3 , 10 , 3 , 1 ,
     &     0 , 2 , 0 , 8 , 0 , 0 , 0 , 0 , 1 , 1 , 0 , 'GN  ' , 'LIN ' ,
     &     '    ' , '    ' , 0 , 3 , 9 , 2 , 0 , 0 , 1 , 0 , 1 , 0 , 0 ,
     &     0 , 0 , 1 , 1 , 0 , 'GN  ' , 'POLY' , '    ' , '    ' , 0 ,
     &     3 , 9 , 2 , 0 , 0 , 1 , 0 , 7 , 0 , 0 , 0 , 0 , 1 , 1 , 0 ,
     &     'GN  ' , 'JUNC' , '    ' , '    ' , 0 , 3 , 9 , 2 , 0 , 0 ,
     &     1 , 0 , 4 , 0 , 0 , 0 , 0 , 1 , 1 , 0 , 'FET ' , 'CURT' ,
     &     '    ' , '    ' , 0 , 3 , 13 , 6 , 3 , 0 , 1 , 0 , 17 , 0 ,
     &     0 , 0 , 0 , 2 , 2 , 0 , 'ICU ' , 'JUNC' , '    ' , '    ' ,
     &     0 , 3 , 11 , 4 , 0 , 0 , 1 , 0 , 5 , 0 , 0 , 0 , 0 , 1 , 1 ,
     &     0 , 'ICU ' , 'POLY' , '    ' , '    ' , 0 , 3 , 11 , 4 , 0 ,
     &     0 , 1 , 0 , 7 , 0 , 0 , 0 , 0 , 1 , 1 , 0 , 'CN  ' , 'DIFF' ,
     &     '    ' , '    ' , 0 , 3 , 9 , 2 , 0 , 0 , 1 , 0 , 5 , 0 , 0 ,
     &     0 , 0 , 2 , 2 , 0 , 'CN  ' , 'BARR' , '    ' , '    ' , 0 ,
     &     3 , 9 , 2 , 0 , 0 , 1 , 0 , 4 , 0 , 0 , 0 , 0 , 2 , 2 , 0 ,
     &     'CN  ' , 'LIN ' , '    ' , '    ' , 0 , 3 , 9 , 2 , 0 , 0 ,
     &     1 , 0 , 1 , 0 , 0 , 0 , 0 , 2 , 2 , 0 , 'CN  ' , 'POLY' ,
     &     '    ' , '    ' , 0 , 3 , 9 , 2 , 0 , 0 , 1 , 0 , 7 , 0 , 0 ,
     &     0 , 0 , 2 , 2 , 0 , 'LIB0' , 'LL0 ' , '    ' , '    ' , 0 ,
     &     2 , 12 , 4 , 0 , 0 , 1 , 0 , 2 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     &     'YTAB' , '    ' , '    ' , '    ' , 0 , 5 , 0 , 0 , 0 , 0 ,
     &     1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 'STAB' , '    ' ,
     &     '    ' , '    ' , 0 , 5 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 ,
     &     0 , 0 , 0 , 0 , 1 , 'LIB ' , 'MPL ' , '    ' , '    ' , 0 ,
     &     2 , 12 , 4 , 0 , 0 , 2 , 0 , 7 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
     &     'LIB ' , 'LANG' , '    ' , '    ' , 0 , 2 , 12 , 4 , 0 , 0 ,
     &     2 , 0 , 12 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 'LIB ' , 'SMPL' ,
     &     '    ' , '    ' , 0 , 2 , 12 , 4 , 0 , 0 , 2 , 0 , 11 , 0 ,
     &     0 , 0 , 0 , 0 , 0 , 0 , 'BIP ' , 'TR  ' , '    ' , '    ' ,
     &     0 , 3 , 14 , 7 , 4 , 0 , 1 , 0 , 23 , 0 , 0 , 0 , 0 , 2 , 2 ,
     &     0 , 'SHL ' , 'KZ  ' , '    ' , '    ' , 0 , 2 , 10 , 2 , 0 ,
     &     0 , 2 , 0 , 3 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 'SHL ' , 'XX  ' ,
     &     '    ' , '    ' , 0 , 2 , 10 , 2 , 0 , 0 , 2 , 0 , 3 , 0 ,
     &     0 , 0 , 0 , 0 , 0 , 0 , 'DISC' , '    ' , '    ' , '    ' ,
     &     0 , 2 , 10 , 2 , 0 , 0 , 2 , 0 , 3 , 0 , 0 , 0 , 0 , 0 , 0 ,
     &     0 , 'INDS' , '    ' , '    ' , '    ' , 0 , 2 , 14 , 6 , 0 ,
     &     0 , 1 , 0 , 12 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 440*0/
C
C      WHEN ENTERING A NEW TYPE OF ELEMENTS, IT IS NECESSARY TO MAKE
C      ADDITIONS TO THE INSTRUCTION FILE INSTIL.TXT (SUBDIRECTORY YINSTR)
C      AND TO CHAPTER5.REP (SUBDIRECTORY REPORT)
C
C      ELEMENT TYPES '4HBIP TR', 'SHL KZ', 'SHL XX'
C      ARE NOT USED (OR NO LONGER USED).
C
      Ka = 28
C
C   FILLING THE ARRAY 'A'
C     A(1) ... A(4) - ELEMENT TYPE NAME;
C     A(6) - FLAG THAT DEFINES THE ELEMENT TYPE:
C              LINEAR BIPOLAR ELEMENT - 1,
C              LINEAR MULTIPOLAR ELEMENT - 2,
C              NONLINEAR ELEMENT - 3,
C              Y-, S-MATRICES - 5.
C     A(7) - DESCRIPTION LINE LENGTH OF ONE ELEMENT
C            OF THIS TYPE IN "NODEEL":
C              FOR LINEAR BIPOLAR ELEMENTS - 5,
C              FOR LINEAR MULTIPOLAR ELEMENTS - 8 + NUMBER OF EXTERNAL NODES,
C              FOR NONLINEAR ELEMENTS - 7 + NUMBER OF EXTERNAL NODES,
C              FOR Y-, S-MATRICES - 0.
C     A(8) - NUMBER OF NODES IN THE MATHEMATICAL MODEL
C            (EXTERNAL + INTERNAL):
C              FOR LINEAR BIPOLAR ELEMENTS - 2,
C              FOR LINEAR MULTIPOLAR ELEMENTS - 0, AS THIS VALUE
C              IS AUTOMATICALLY DETERMINED FOR EACH LINEAR MULTIPOLAR ELEMENT,
C              FOR NONLINEAR ELEMENTS - NUMBER OF NODES,
C              FOR Y-, S-MATRICES - 0.
C    A(9) - NUMBER OF INTERNAL NODES, AUTOMATICALLY NUMBERED
C           WHEN DESCRIBING NODES OF THE ELEMENT.
C    A(11) - NUMBER OF PARAMETERS COMMON TO THIS TYPE (IF PARAMETERS
C            ARE ABSENT, ENTER 1).
C    A(13) - NUMBER OF INDIVIDUAL PARAMETERS DESCRIBING
C            THE MATHEMATICAL MODEL OF THE ELEMENT.
C    A(18) - FOR NONLINEAR ELEMENTS, ENTER THE MAXIMUM
C            NUMBER OF INPUT VARIABLES IN THE CURRENT CALCULATION.
C    A(19) - FOR NONLINEAR ELEMENTS, ENTER THE MAXIMUM
C            NUMBER OF OUTPUT VARIABLES IN THE CURRENT CALCULATION.
C
C
C  ELEMENT TYPE "R"        SUBROUTINE LIN2P         (file FORMER)
C  ELEMENT TYPE "L"        SUBROUTINE LIN2P         (file FORMER)
C  ELEMENT TYPE "C"        SUBROUTINE LIN2P         (file FORMER)
C  ELEMENT TYPE "E"        SUBROUTINE EMF           (file LIBLIN)
C  ELEMENT TYPE "J"        SUBROUTINE JDRFE         (file LIBLIN)
C  ELEMENT TYPE "P2"       SUBROUTINE RFFIPN        (file LIBLIN)
C  ELEMENT TYPE "VD SCHT"  SUBROUTINE MDSCH1-MDSCH6 (file MDSCH)
C  ELEMENT TYPE "GN LLIN"  SUBROUTINE LIN1-LIN5     (file LIN)
C  ELEMENT TYPE "GN POLY"  SUBROUTINE POLY51-POLY55 (file POLY5)
C  ELEMENT TYPE "GN JUNC"  SUBROUTINE JUNC1-JUNC6   (file JUNC)
C  ELEMENT TYPE "FET CURT" SUBROUTINE CURT1-CURT6   (file CURT)
C  ELEMENT TYPE "ICU JUNC" SUBROUTINE ICUJ1-ICUJ6   (file ICUJUNC)
C  ELEMENT TYPE "ICU POLY" SUBROUTINE ICUPL1-ICUPL5 (file ICUPOLY)
C  ELEMENT TYPE "CN DIFF"  SUBROUTINE CDIFF1-CDIFF6 (file CDIFF)
C  ELEMENT TYPE "CN BARR"  SUBROUTINE CBARR1-CBARR5 (file CBARR)
C  ELEMENT TYPE "CN LIN"   SUBROUTINE CLIN1-CLIN5   (file CLIN)
C  ELEMENT TYPE "CN POLY"  SUBROUTINE CPOLY1-CPOLY5 (file CPOLY)
C  ELEMENT TYPE "LIB0 LL0" SUBROUTINE LINE1,LINE2   (file LIB0)
C  ELEMENT TYPE "YTAB"     SUBROUTINE YTAB          (file YTAB)
C  ELEMENT TYPE "STAB"     SUBROUTINE STAB          (file STAB)
C  ELEMENT TYPE "LIB MPL"  SUBROUTINE MP            (file MPL)
C  ELEMENT TYPE "LIB LANG" SUBROUTINE LANG          (file LANGE)
C  ELEMENT TYPE "LIB SMPL" SUBROUTINE SMPL          (file SMPL)
C  ELEMENT TYPE "BIP TR"   SUBROUTINE BIPTR1-BIPTR5 (file BIPTR)
C  ELEMENT TYPE "SHL KZ"   SUBROUTINE SHLEIF        (file SHLEIF)
C  ELEMENT TYPE "SHL XX"   SUBROUTINE SHLEIF        (file SHLEIF)
C  ELEMENT TYPE "DISC"     SUBROUTINE DISCONT       (file DISCONT)
C  ELEMENT TYPE "INDS"     SUBROUTINE INDSV         (file INDSV)
C
C
      DATA Iff/'P2  ' , '    ' , '    ' , '    ' , 2 , 0 , 0 , 'P2  ' ,
     &     '    ' , '    ' , '    ' , 5 , 0 , 0 , 'E   ' , '    ' ,
     &     '    ' , '    ' , 2 , 0 , 0 , 'J   ' , '    ' , '    ' ,
     &     '    ' , 2 , 0 , 0/

      DATA Kiff/4/ , Knniff/0/
C
      END
