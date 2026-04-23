c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE STEPFR (OM,NREC,Y,VJ,ISIZE_MAXNODE)
C*********************************************************************
C* Subroutine for forming and reduction of multi-channel linear      *
C* subsystem and organizing their I/O on the DA.                     *
C*                                                                   *
C*********************************************************************

C Change from 30.01.91. See MAIN
C$LARGE: BUFFER
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/MATY/  BUFFER (6000),BUFLEN
      DOUBLE COMPLEX       BUFFER
      INTEGER*4     BUFLEN
      COMMON/BLK2/  KNC,KNR,KN,KNR1,KN1

      INTEGER ISIZE_MAXNODE
      DOUBLE COMPLEX VJ(ISIZE_MAXNODE)
      DOUBLE COMPLEX Y(ISIZE_MAXNODE,ISIZE_MAXNODE)

      COMMON/KOLNAL/ KOL,NAL
      INTEGER        KOL(4),DY,FLAG
      LOGICAL        NAL(4)
      DOUBLE PRECISION           OM

C Change from 30.01.91. See MAIN
C Built-in functions for indexing:
      IFIND1(I,M,NU)=(I+(M-1)*NU)
      IFIND2(I,J,M,NU,MF)=NU*MF+I+(J-1)*NU+(M-1)*NU*NU


      DY=ISIZE_MAXNODE
C @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      KKK=KOL(1)+KOL(2)+KOL(3)

C Initialization of the matrices Y and the vector VJ
      CALL ZINY(Y,VJ,ISIZE_MAXNODE)
C      print *,'after ZINY'
C      print *,(VJ(ikkk),ikkk=1,20)

      IF(.NOT.NAL(1)) GO TO 20

C     LIN.CONST. - FORMATION
      CALL LNFORM (OM,Y,VJ,ISIZE_MAXNODE)

C   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C      WRITE(6,2)((III,JJJ,Y(III,JJJ),JJJ=1,KKK),III=1,KKK)
C   1  FORMAT(2X,' Y   STEPFR ',I2,' KOL(1),KOL(2),KOL(3),KKK=',4I3)
C   2  FORMAT(2X,'Y(',I3,',',I3,')=',1X,E13.6,1X,E13.6)
C      WRITE(6,3) (III, VJ(III),III=1,KKK)
C   3  FORMAT(2X,'STEPFR VJ(',I3,')=',E13.6,',',E13.6)

      IF(KOL(1).EQ.0) GO TO 10
      NF=1
      NEND=KOL(1)
      N=KOL(1)+KOL(2)+KOL(3)
C  Reduction Y
      CALL LUSLV (Y,DY,N,NF,NEND,FLAG)

C  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C  Reduction VJ
      CALL LUFRW (Y,VJ,DY,N,NF,NEND,FLAG)

C  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C      write (6,*) 'Y matrix'
C      do ii=1,n
C            write (6,120) (Y(ii,jj), jj=1,n)
C      enddo
C120   format (2x,'(',1x,e12.5,1x,e12.5,')')


C Packing and writing of matrix Y and vector VJ to DA
   10 CALL PACK1 (NREC,Y,VJ,ISIZE_MAXNODE)
      GO TO 20
C**********************************************************************
      ENTRY DOUBLE(OM,NREC)
C**********************************************************************
C*                                                                    *
C* Entry in subroutine for reformation (when modifying elements)      *
C*                                                                    *
C**********************************************************************

C Initialization of Y
      CALL ZINY(Y,VJ,ISIZE_MAXNODE)
C Reading from DA and unpacking the formatted Y and VJ from constant elements
      CALL DPACK1(NREC,Y,VJ,ISIZE_MAXNODE)



   20 IF(.NOT.NAL(2)) GO TO 30

C Linear Variable - Initialization
      CALL VLFORM(OM,Y,VJ,ISIZE_MAXNODE)

C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      IF(KOL(2).EQ.0) GO TO 30
      NF=KOL(1)+1
      NEND=KOL(1)+KOL(2)
      N=KOL(1)+KOL(2)+KOL(3)
C Reduction of the full matrix LIN.-subsystem
      CALL LUSLV (Y,DY,N,NF,NEND,FLAG)

C @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C Reduction of the full vector of the given currents in the LIN.-subsystem
      CALL LUFRW (Y,VJ,DY,N,NF,NEND,FLAG)


C @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C Transformation of matrix Y to canonical form
   30 CONTINUE
      N=KOL(1)+KOL(2)+KOL(3)
      NF=KOL(1)+KOL(2)+1
      CALL LUCAN(Y,DY,N,NF,FLAG)

C @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C Execution of the operation for setting the input values
C (IF(OMEGA.NE.0.)). Linked with the DPF.
C
C
      K123=KOL(1)+KOL(2)+KOL(3)
      IF(NREC.EQ.1)GOTO 80
      DO 90 JJ=1,K123
   90 VJ(JJ)=VJ(JJ)/2.D0
   80 CONTINUE

C Packing and writing Y to MD
      CALL PACK2 (NREC,Y,VJ,ISIZE_MAXNODE)

C Transfer of result to the function and reductions in the /MATY/ block,
C storing Y and J. Reduced linear subckt of the entire frequency range.
      K3=KOL(3)
      K12=KOL(1)+KOL(2)
      DO 70 KI=1,K3
      IND1=KI+K12

C Change from 30.01.91, modified by Serdyuk G.V.
C      JR(KI,NREC)=VJ(IND1)       **** OLD VERSION
      BUFFER( IFIND1(KI,NREC,K3) )=VJ(IND1)

      DO 70 KJ=1,K3
      IND2=K12+KJ

C Change from 30.01.91, modified by Serdyuk G.V.
C      YR(KJ,KI,NREC)=Y(IND2,IND1)   **** OLD VERSION
      ILOCAL=IFIND2(KJ,KI,NREC,K3,KN)
      IF (ILOCAL.GT.BUFLEN) GOTO 100
      BUFFER( ILOCAL )=Y(IND2,IND1)

C  NO    WRITE(6,3000) KJ,KI,NREC,YR(KJ,KI,NREC)
C  NO 3000 FORMAT(2X,'STEPFR,YR(',I3,',',I3,',',I3,')=',E12.5,2X,E12.5)
C      WRITE(6,3003) IND2,IND1,Y(IND2,IND1)
C 3003 FORMAT(2X,'        Y(',I3,',',I3,')=     ',E12.5,2X,E12.5)
   70 CONTINUE
      RETURN
100   WRITE(6,2000)
      STOP
C     DEBUG SUBTRACE,INIT(K12,K3,KI,KJ,IND1,IND2)

2000  FORMAT (2X,' BUFFER SIZE FOR REDUCED MATRICES EXCEDED'/        2X,
     +' ABNORMAL TERMINATION.')



      END
