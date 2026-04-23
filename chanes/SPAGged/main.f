c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c

      PROGRAM HBL

c*********************************************************************
c expected command line options
c   usage: hb[.exe] <infile> [<outfile> -f<freq-table> -p<p-freq-table> -u
c     <init> -a]
c          infile         - HBP input file
c          outfile        - HBP output file, optional, stdout if not set
c          freq-table     - output variables in frequency domain
c          p-freq-table   - "pulsed" output variables in frequency domain
c          uinit          - file of input variables (not supported now)
c          -a             - about
c
c
c
c
C*********************************************************************
C*                                                                   *
C*   MAIN PROGRAMM = SPEC'HOB'( Version 2.0 ),    11.11.87           *
C*                                                                   *
C*      DATE LAST EDITED  -  28.04.92                                *
C*      last edited August-2002                                      *
C*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      include 'circuit.i'
      COMMON/PRINT /   KPRLEN,KPRSRT,KPRNKR,KPRLIN,KPRSOL,KPRVAR,
     +          KPRGRF,KPRQUP
      COMMON/NEWTON/   EPSSOL,EPSDU,EPSMIN,MAXDU,LIMIT
      DOUBLE PRECISION             EPSSOL,EPSDU,EPSMIN,MAXDU
      COMMON/TYPVAL/   TYPU,TYPI
      DOUBLE PRECISION             TYPU,TYPI
      COMMON/FRE/      F(2)
      COMMON/KOLNAL/   KOL,NAL
      COMMON/MEP/      MEPHF,FLGMNW

      INTEGER ISIZE_UIN
      PARAMETER (ISIZE_MAXVAR = 1000)
      DOUBLE PRECISION UIN(ISIZE_MAXVAR)
      INTEGER ISIZE_MAXNODE
      PARAMETER (ISIZE_MAXNODE = 200)
      INTEGER NOUZ(ISIZE_MAXNODE), INOUZ1(ISIZE_MAXNODE)
C ISIZE_MAXVAR and ISIZE_MAXNODE are really different :
C ISIZE_MAXVAR - represents maximum number of variables (reduced number of nodes
C                X number of frequencies
C ISIZE_MAXNODE - number of nodes in circuit before reduction
C these numbers are not checked yet for circuit ... :-(
      DOUBLE COMPLEX VECTJ(ISIZE_MAXNODE)
      DOUBLE COMPLEX YY(ISIZE_MAXNODE,ISIZE_MAXNODE)

      DOUBLE COMPLEX S(ISIZE_MAXNODE,20)


      COMMON/BLIFF/IFF(7,4),KIFF,NNIFF(4,8),KNNIFF,PNIFF(8),FNE1,FNE2

C  CHANGE FROM 01/30/91 CHANGED BY G.V. SERDIUK
C$LARGE: BUFFER
      COMMON/MATY/     BUFFER (6000),BUFLEN
      DOUBLE COMPLEX          BUFFER
      INTEGER*4        BUFLEN


      COMMON/BLW1/     W,W1
      COMMON/BLW2/     WR,WS
      COMMON/BLK1/     KR,KR1,KC,KC1,NNR,NNR1,MN,MN1
      COMMON/BLK2/     KNC,KNR,KN,KNR1,KN1
      COMMON/BLMNI/    MNI(2,20),KMNI
      COMMON/MODGLB/   MGLOB,IAPR
      COMMON/SERV/     EPSIW,LIMERR,KITU
      COMMON/KOLUZ/    KMOD,KSCH
      COMMON/TER/      ITERMC

      DOUBLE PRECISION             PI/3.14159D0/,WIRCPI,WI12PI

      LOGICAL          NAL(4)
      INTEGER          KOL(4)
      INTEGER          WR(20,20),WS(20,20),FLGMNW,MN(2,20),MN1(2,200)
      INTEGER          KR(20),KC(10),NNR(10),KR1(200),KC1(20),NNR1(20)
      DOUBLE PRECISION W(20),W1(200)

      CHARACTER*4  NAM(13)/'  KR','  KC',' NNR',' KR1',' KC1','NNR1',
     +           '  MN',' MN1','  WR','  WS','   W','  W1','   F'/
      CHARACTER*4      FNE1,FNE2,   WW3
      CHARACTER*12     WWW2,FS1
      INTEGER IERR

      DATA             IIN/10/,IOUT/6/

      NAMELIST/UINIT/  UIN

      integer n_comline_arg, n_arg
      character*256 infile, outfile
      character*256 cl_param

C  SEE CORRECTION ABOVE FROM 30.01.91 BY SERDYUK G.V.
      IFIND1(I,M,NU)=(I+(M-1)*NU)
      IFIND2(I,J,M,NU,MF)=NU*MF+I+(J-1)*NU+(M-1)*NU*NU
      BUFLEN =6000
C
c =========================================================
c parsing command line
      n_comline_arg=iargc()
      if (n_comline_arg.eq.0) then
c  output usage message
      print *,'usage: hb[.exe] <infile> [<outfile> -f<freq-table>',
     +        '                -p<p-freq-table> -u<init> -a]'
      print *,' infile       - HBP input file; '
      print *,' outfile      - HBP output file, optional, stdout',
     +                                              ' if not set'
      print *,' freq-table   - output variables in frequency domain'
      print *,' p-freq-table - "pulsed" output variables in frequency',
     +                                              ' domain'
      print *,' uinit        - file of input variables (not',
     +                                               ' supported now)'
      print *,' -a           - about'
      stop
      endif

c  start parsing
c  assign defaults
      outfile=' '
      infile=' '
      do n_arg=1,n_comline_arg
        call getarg(n_arg, cl_param)
        if(cl_param(1:1).ne.'-') then
c this is not a parameter. store it: first infile, then - outfile
            if(infile(1:1).eq.' ') then
                infile=cl_param
                print *, 'in:',infile
            elseif(outfile(1:1).eq.' ') then
                outfile=cl_param
                print *, 'out:',outfile
            else
            print *, 'unknown string',cl_param
            endif
        else
c need not to open channel 6 for aut - it is stdout. parse switches
            if (cl_param(2:2).eq.'f') then
                print *,'f parameter'
            elseif(cl_param(2:2).eq.'p') then
                print *,'p parameter'
            elseif(cl_param(2:2).eq.'u') then
                print *,'u parameter'
            elseif(cl_param(2:2).eq.'a') then
                print *,'HArmonic BALAnce simulator, ',
     +                  '(c) Gennady Serdyuk, 1989-2002',
     +                  ' gserdyuk@mail.ru'
c                return
            else
                print *,'unknown parameter', cl_param
            endif
        endif
      enddo

      if(outfile(1:1).ne.' ') then
c not default value - open file
        OPEN(6,FILE=cl_param)
      endif
      if(infile(1:1).ne.' ') then
        OPEN(10,FILE=infile)
      endif


C =========================================================
C ##!##  HERE IN PACKET U=UMAX/2 KPOME U WHEN OMEGA=0. ANALOGICAL VECTJ.
C
C ***READING THE INTRODUCTION TASK
      WRITE (6,460)
      PRINT    460
C
C
C *********************************************************
      CALL LEN_S
C *********************************************************

      WRITE (6, 460)
C
C    npass - pass number through the program from 2000 CONTINUE
C                                              to  300 CONTINUE
C                                          ..... GO TO 2000
C                                        28.02.92
      npass=0
2000  CONTINUE
      npass=npass+1

C *********************************************************
      CALL INITI(UIN,YY,VECTJ,ISIZE_MAXNODE,IERR)
      IF (IERR.eq.1) then
        goto 1000
      elseif (IERR.eq.2) then
        goto 1010
      elseif (IERR.eq.3) then
        goto 1020
      endif
C else - just continue
C *********************************************************

      IF(KPRNKR.LE.2) GOTO 1800
      WRITE (6, 320) KNC,KN
      WRITE (6, 220) ((NAM(7),I,J,MN(I,J),I=1,2),J=1,KN)
      WRITE (6, 270)
      WRITE (6, 210) (NAM(1),I,KR(I),I=1,KN)
      WRITE (6, 270)
      WRITE (6, 210) (NAM(2),I,KC(I),I=1,KNR)
      WRITE (6, 270)
      WRITE (6, 210) (NAM(3),I,NNR(I),I=1,KNR)
      WRITE (6, 270)
      WRITE (6, 210) (NAM(4),I,KR1(I),I=1,KN1)
      WRITE (6, 270)
      WRITE (6, 210) (NAM(5),I,KC1(I),I=1,KNR1)
      WRITE (6, 270)
      WRITE (6, 210) (NAM(6),I,NNR1(I),I=1,KNR1)
      WRITE (6, 270)
      WRITE (6, 220) ((NAM(8),I,J,MN1(I,J),I=1,2),J=1,KN1)

      WRITE (6, 290) NAM(9),(I,I=1,KN)
      WRITE (6, 270)
      DO 501 I=1,KN
  501 WRITE (6, 350) I,(WR(I,J),J=1,KN)

      WRITE (6, 290) NAM(10),(I,I=1,KN)
      WRITE (6, 270)
      DO 502 I=1,KN
  502 WRITE (6, 350) I,(WS(I,J),J=1,KN)

      WRITE (6, 270)
      WRITE (6, 230) (NAM(11),I,W(I),I=1,KN)
      WRITE (6, 270)
      WRITE (6, 230) (NAM(12),I,W1(I),I=1,KN1)
      WRITE (6, 280)
 1800 CONTINUE
C
      K1=KOL(1)
      K2=KOL(2)
      K3=KOL(3)
      K12=K1+K2
      K23=K2+K3
      K123=K12+K3
      KO0=K3*KN
C

C
      IF(KPRLIN.GE.2) WRITE (6,400)
                      PRINT    400
      DO 50 I1=1,KN
      WI12PI=W(I1)/(2.D0*PI)
      IF(KPRLIN.GE.2) WRITE (6, 410) WI12PI
                      PRINT     410, WI12PI
      DO 50 I2=1,K3
      DO 49 I3=1,K3
      IF(KPRLIN.GE.2)         WRITE (6, 421) I3,I2,I1,BUFFER(IFIND2(I3,I
     +2,I1,K3,KN))
               PRINT     421, I3,I2,I1,BUFFER(IFIND2(I3,I2,I1,K3,KN))

 421  FORMAT(2X,'Y(',I3,',',I3,',',I3,')=',E12.5,2X,E12.5)
  49  CONTINUE



      IF(KPRLIN.GE.2) WRITE (6, 415) I2,I1, BUFFER(IFIND1(I2,I1,K3))
                      PRINT     415, I2,I1, BUFFER(IFIND1(I2,I1,K3))
  50  CONTINUE
C
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     READING AND WRITING OF INITIAL STRESS APPROXIMATIONS
C     FROM A FILE WITH THE EXTENSION ".UIN"
      IF(npass.ne.1.or.kitu.ne.3) GOTO 60
C
C     ASSIGNING THE SYMBOLIC VARIABLE WWW2 THE NAME OF
C     THE INPUT FILE.
C      INQUIRE (10, NAME=WWW2)
C      WW3='.UIN'
C      CALL NACH(WW3, FS1, WWW2)
C     THE SUBROUTINE NACH ASSIGNS THE NAME OF THE INPUT
C     FILE WITH THE EXTENSION .UIN TO THE SYMBOLIC VARIABLE FS1
C    so far - initial aproximation file wil have name UIN.UIN.
C we will return to this later
C NACH is removed
      FS1='UIN.UIN'

      OPEN (12, FILE=FS1)
      READ(12,UINIT,ERR=55)
      CLOSE (12, STATUS='DELETE')
      GOTO 80
   55 CONTINUE
      CLOSE (12, STATUS='DELETE')
   60 CONTINUE

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      if (npass.gt.1.and.kitu.ge.2)GOTO 80
      KO02=2*KO0
      DO 70 IDINT=1,KO02
   70 UIN(IDINT)=0.0D0
   80 CONTINUE
C
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IF(KITU.EQ.1) READ(IIN,UINIT)
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C
C
C
C *********************************************************
      CALL SOLVE(UIN)
C *********************************************************
C
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C       WRITING THE FILE OF INITIAL STRESS APPROXIMATIONS
C                  WITH THE EXTENSION "*UIN"
      IF(NPASS.EQ.1.AND.KITU.EQ.3) OPEN(12,file=fs1)
      IF(NPASS.EQ.1.AND.KITU.EQ.3) WRITE(12,UINIT)
      IF(NPASS.EQ.1.AND.KITU.EQ.3) KITU=2
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C *********************************************************
 1000 CALL RESTOR(UIN,YY,VECTJ,ISIZE_MAXNODE,S)
C *********************************************************
C
C
C
      K123=KOL(1)+KOL(2)+KOL(3)
      IF(KN.LT.2) GOTO 110
      DO 90 IRC=2,KN
      DO 90 I=1,K123
      S(I,IRC)=2.D0*S(I,IRC)
C    BRING TO THE ACTUAL VALUE
   90 S(I,IRC)=S(I,IRC)/DSQRT(2.D0)
C  90 CONTINUE
  110 CONTINUE
C
C     IF(KPRSOL.LE.1) GOTO 150
      DO 130 IRC=1,KN
      WIRCPI=W(IRC)/(2.D0*PI)
      IF(KPRSOL.GT.2) WRITE (6, 430)  IRC,MN(1,IRC),MN(2,IRC),WIRCPI
      PRINT     430,  IRC,MN(1,IRC),MN(2,IRC),WIRCPI
      DO 129 I=1,K123
      IF(KPRSOL.GT.2) WRITE (6, 440)  I,S(I,IRC)
      PRINT     440,  I,S(I,IRC)
  129 CONTINUE
  130 CONTINUE

      CALL wrtraw(ISIZE_MAXNODE,k123,Kn,W,s,pi,infile)

 1010 CONTINUE
 1020 CONTINUE
C
C
  150 CONTINUE
C
C    IF THE ITERATION LIMIT IS EXCEEDED OR THE METHOD
C    DOES NOT CONVERGE, WE DO NOT CALCULATE QUALITY INDICATORS
      IF(ITERMC.EQ.3.OR.ITERMC.EQ.4) GOTO 1050
C
C
C *********************************************************
C    instead of call QUPP - just rudimentary
C       NQUPV=0
C      CALL QUPP
C *********************************************************

 1050 CONTINUE
C ...
  300 CONTINUE
C      IF(KOLVAR.EQ.0)
      WRITE (6, 460)
      WRITE (6, 450)
      WRITE (6, 460)

      STOP

  210 FORMAT((1X,7(A5,'(',I2,')=',I3)))
  220 FORMAT((1X,6(A5,'(',I2,',',I2,')=',I3)))
  230 FORMAT((1X,5(A5,'(',I2,')=',E12.5)))
  270 FORMAT(1X)
  280 FORMAT(//1X,120('*'))
  290 FORMAT(//10X,'MATRIX ',A4/10X,12('_')/10X,20I5)
  320 FORMAT(/3X,'KNC=',I3,5X,'KN=',I3)
  350 FORMAT(I8,2X,20I5)
  400 FORMAT(10X,'RESULTS OF Y-MATRICES FORMING   ')
  410 FORMAT(10X,'FREQUENCY',E12.5)
  415 FORMAT(2X,'VECTJ(',I3,I4,')=',E12.5,2X,E12.5)
  420 FORMAT(10E12.5)
  430 FORMAT(5X,' FREQUENCY',I3,' COMBINATION(',I3,',',I3,')',
     + ' VALUE    ',E13.6)
  440 FORMAT(2X,'U(',I3,')=',E13.6,1X,E13.6)
  450 FORMAT(20X,'***************** END *************************')
  460 FORMAT(1X,79('*')/1X,7('*** hbfree '),'**'/1X,79('*'))
  480 FORMAT(10X,'SOLUTION OF NLN EQ. STARTED ')

      END

      SUBROUTINE wrtraw(Isize_maxnode,K123,Kn,W,S,Pi,Infile)
      INTEGER irc, Isize_maxnode , K123 , Kn
      DOUBLE PRECISION W(20) , wircpi , Pi
      DOUBLE COMPLEX S(Isize_maxnode,20)
      CHARACTER*64 Infile, line
      INTEGER i , jrc , inend , ifend, hblnode, lenname , nvar
      CHARACTER*9 rawfile / '         '/
      CHARACTER*15 nodefile / '               '/
      CHARACTER*6 spicenode(20)
      INTEGER dt(8)
      CHARACTER*10 b(3)

      inend = 0
      ifend = 0
      DO i = 1 , 12
         IF ( Infile(i:i).EQ.'.' ) THEN
            inend = i
            EXIT
         ENDIF
         IF ( Infile(i:i).NE.' ' ) ifend = i
      ENDDO
      IF ( inend.EQ.0 ) THEN
         IF (ifend.GT.5) STOP ' NAME TOO LONG'
         rawfile(1:ifend) = Infile(1:ifend)
         rawfile(ifend+1:ifend+4) = '.raw'
         nodefile(1:ifend) = Infile(1:ifend)
         nodefile(ifend+1:ifend+10) = '.ckt.nodes'
      ELSE
         IF (inend.GT.6) STOP ' NAME TOO LONG'
         rawfile(1:inend-1) = Infile(1:inend-1)
         rawfile(inend:inend+4) = '.raw'
         nodefile(1:inend-1) = Infile(1:inend-1)
         nodefile(inend:inend+9) = '.ckt.nodes'
      ENDIF

      OPEN (15, FILE=trim(nodefile),STATUS='OLD')
      nvar = 0
      DO i = 1 , K123+5
         READ (15,'(A20)',END=1500) line
         IF (i.GT.5) THEN
            nvar = nvar + 1
            READ (line,'(I1,TR1,A)') hblnode, spicenode(nvar)
            IF (hblnode .eq. 0) THEN
C              SKIP GND NODE
               nvar = nvar - 1
            ENDIF
         ENDIF
      ENDDO
 1500 CONTINUE
 
      call date_and_time(b(1), b(2), b(3), dt)

      OPEN (16, FILE=trim(rawfile))
 
      WRITE (16,'(A)') 'Title: spice test'
      WRITE (16, 170) dt(3), dt(2), dt(1), dt(5), dt(6), dt(7)
      WRITE (16,'(A)') 'Plotname:  Harmonic Balance Simulation'
      WRITE (16,'(A)') 'Flags: complex'
      WRITE (16,'(A,I4)') 'No. Variables: ' , nvar+1
      WRITE (16,'(A,I4)') 'No. Points: ' , 3*Kn
      WRITE (16,'(A)') 'Command:  version 3f5'
 
      WRITE (16,'(A)') 'Variables:'
      WRITE (16,*) '     0    frequency    frequency'
      DO i = 1 , nvar
         lenname = index(spicenode(i),achar(9))-1
         WRITE (16,140) i , spicenode(i)(1:lenname)
 140     FORMAT (3X,I4,'    V(', A, ')    voltage')
      ENDDO
 
      WRITE (16,'(A)') 'Values:'
      jrc = 0
      DO irc = 1 , Kn
         wircpi = W(irc)/(2.D0*Pi)
         WRITE (16,150) jrc , wircpi , 0.0
         DO i = 1 , nvar
            WRITE (16,160) 0.0 , 0.0
         ENDDO
         jrc = jrc + 1
         WRITE (16,150) jrc , wircpi , 0.0
         DO i = 1 , nvar
            WRITE (16,160) dreal(S(i,irc)) , dimag(S(i,irc))
         ENDDO
         jrc = jrc + 1
         WRITE (16,150) jrc , wircpi , 0.0
         DO i = 1 , nvar
            WRITE (16,160) 0.0 , 0.0
         ENDDO
         jrc = jrc + 1
      ENDDO
 
      CLOSE (15)
      CLOSE (16)
 150  FORMAT (1X,I4,4X,SP,E20.13,',',SP,E20.13)
 160  FORMAT (9X,SP,E20.13,',',SP,E20.13)
 170  FORMAT ('Date: ',I2.2,'.',I2.2,'.',I4, ' ',I2.2, ':', I2.2, ':',
     & I2.2)

      END
