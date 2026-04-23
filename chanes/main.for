!*==HBL.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c

      PROGRAM hbl

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
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Epsiw , F , Pniff
      INTEGER i , i1 , i2 , i3 , Iapr , iargc , idint , Iff , ifind1 ,
     &        ifind2 , iin , iout , irc , ISIZE_MAXVAR , Itermc , j ,
     &        k1 , k12 , k123 , k2
      INTEGER k23 , k3 , Kiff , Kitu , Kmni , Kmod , Kn , Kn1 , Knc ,
     &        Knniff , Knr , Knr1 , ko0 , ko02 , Kprgrf , Kprlen ,
     &        Kprlin , Kprnkr , Kprqup , Kprsol
      INTEGER Kprsrt , Kprvar , Ksch , Limerr , Limit , m , Mephf , mf ,
     &        Mglob , Mni , Nniff , npass , nu
c      include 'circuit.i'
      COMMON /print / Kprlen , Kprsrt , Kprnkr , Kprlin , Kprsol ,
     &                Kprvar , Kprgrf , Kprqup
      COMMON /newton/ Epssol , Epsdu , Epsmin , Maxdu , Limit
      DOUBLE PRECISION Epssol , Epsdu , Epsmin , Maxdu
      COMMON /typval/ Typu , Typi
      DOUBLE PRECISION Typu , Typi
      COMMON /fre   / F(2)
      COMMON /kolnal/ Kol , Nal
      COMMON /mep   / Mephf , Flgmnw

      PARAMETER (ISIZE_MAXVAR=1000)
      DOUBLE PRECISION uin(ISIZE_MAXVAR)
      INTEGER ISIZE_MAXNODE
      PARAMETER (ISIZE_MAXNODE=200)
C ISIZE_MAXVAR and ISIZE_MAXNODE are really different :
C ISIZE_MAXVAR - represents maximum number of variables (reduced number of nodes
C                X number of frequencies
C ISIZE_MAXNODE - number of nodes in circuit before reduction
C these numbers are not checked yet for circuit ... :-(
      DOUBLE COMPLEX vectj(ISIZE_MAXNODE)
      DOUBLE COMPLEX yy(ISIZE_MAXNODE,ISIZE_MAXNODE)

      DOUBLE COMPLEX s(ISIZE_MAXNODE,20)


      COMMON /bliff / Iff(7,4) , Kiff , Nniff(4,8) , Knniff , Pniff(8) ,
     &                Fne1 , Fne2

C  CHANGE FROM 01/30/91 CHANGED BY G.V. SERDIUK
C$LARGE: BUFFER
      COMMON /maty  / Buffer(6000) , Buflen
      DOUBLE COMPLEX Buffer
      INTEGER*4 Buflen


      COMMON /blw1  / W , W1
      COMMON /blw2  / Wr , Ws
      COMMON /blk1  / Kr , Kr1 , Kc , Kc1 , Nnr , Nnr1 , Mn , Mn1
      COMMON /blk2  / Knc , Knr , Kn , Knr1 , Kn1
      COMMON /blmni / Mni(2,20) , Kmni
      COMMON /modglb/ Mglob , Iapr
      COMMON /serv  / Epsiw , Limerr , Kitu
      COMMON /koluz / Kmod , Ksch
      COMMON /ter   / Itermc

      DOUBLE PRECISION pi/3.14159D0/ , wircpi , wi12pi

      LOGICAL Nal(4)
      INTEGER Kol(4)
      INTEGER Wr(20,20) , Ws(20,20) , Flgmnw , Mn(2,20) , Mn1(2,200)
      INTEGER Kr(20) , Kc(10) , Nnr(10) , Kr1(200) , Kc1(20) , Nnr1(20)
      DOUBLE PRECISION W(20) , W1(200)

      CHARACTER*4 nam(13)/'  KR' , '  KC' , ' NNR' , ' KR1' , ' KC1' ,
     &            'NNR1' , '  MN' , ' MN1' , '  WR' , '  WS' , '   W' ,
     &            '  W1' , '   F'/
      CHARACTER*4 Fne1 , Fne2
      CHARACTER*12 fs1
      INTEGER ierr

      DATA iin/10/ , iout/6/

      NAMELIST /uinit / uin

      INTEGER n_comline_arg , n_arg
      CHARACTER*256 infile , outfile
      CHARACTER*256 cl_param

C  SEE CORRECTION ABOVE FROM 30.01.91 BY SERDYUK G.V.
      ifind1(i,m,nu) = (i+(m-1)*nu)
      ifind2(i,j,m,nu,mf) = nu*mf + i + (j-1)*nu + (m-1)*nu*nu
      Buflen = 6000
C
c =========================================================
c parsing command line
      n_comline_arg = iargc()
      IF ( n_comline_arg.EQ.0 ) THEN
c  output usage message
         PRINT * , 'usage: hb[.exe] <infile> [<outfile> -f<freq-table>'
     &         , '                -p<p-freq-table> -u<init> -a]'
         PRINT * , ' infile       - HBP input file; '
         PRINT * , ' outfile      - HBP output file, optional, stdout' ,
     &         ' if not set'
         PRINT * ,
     &         ' freq-table   - output variables in frequency domain'
         PRINT * ,
     &         ' p-freq-table - "pulsed" output variables in frequency'
     &         , ' domain'
         PRINT * , ' uinit        - file of input variables (not' ,
     &         ' supported now)'
         PRINT * , ' -a           - about'
         STOP
      ENDIF

c  start parsing
c  assign defaults
      outfile = ' '
      infile = ' '
      DO n_arg = 1 , n_comline_arg
         CALL getarg(n_arg,cl_param)
         IF ( cl_param(1:1).NE.'-' ) THEN
c this is not a parameter. store it: first infile, then - outfile
            IF ( infile(1:1).EQ.' ' ) THEN
               infile = cl_param
               PRINT * , 'in:' , infile
            ELSEIF ( outfile(1:1).EQ.' ' ) THEN
               outfile = cl_param
               PRINT * , 'out:' , outfile
            ELSE
               PRINT * , 'unknown string' , cl_param
            ENDIF
c need not to open channel 6 for aut - it is stdout. parse switches
         ELSEIF ( cl_param(2:2).EQ.'f' ) THEN
            PRINT * , 'f parameter'
         ELSEIF ( cl_param(2:2).EQ.'p' ) THEN
            PRINT * , 'p parameter'
         ELSEIF ( cl_param(2:2).EQ.'u' ) THEN
            PRINT * , 'u parameter'
         ELSEIF ( cl_param(2:2).EQ.'a' ) THEN
            PRINT * , 'HArmonic BALAnce simulator, ' ,
     &            '(c) Gennady Serdyuk, 1989-2002' , ' gserdyuk@mail.ru'
c                return
         ELSE
            PRINT * , 'unknown parameter' , cl_param
         ENDIF
      ENDDO

c not default value - open file
      IF ( outfile(1:1).NE.' ' ) OPEN (6,FILE=cl_param)
      IF ( infile(1:1).NE.' ' ) OPEN (10,FILE=infile)


C =========================================================
C ##!##  HERE IN PACKET U=UMAX/2 KPOME U WHEN OMEGA=0. ANALOGICAL VECTJ.
C
C ***READING THE INTRODUCTION TASK
      WRITE (6,460)
      PRINT 460
C
C
C *********************************************************
      CALL len_s
C *********************************************************

      WRITE (6,460)
C
C    npass - pass number through the program from 2000 CONTINUE
C                                              to  300 CONTINUE
C                                          ..... GO TO 2000
C                                        28.02.92
      npass = 0
      npass = npass + 1

C *********************************************************
      CALL initi(uin,yy,vectj,ISIZE_MAXNODE,ierr)
      IF ( ierr.EQ.1 ) GOTO 1000
      IF ( ierr.EQ.2 ) GOTO 1010
      IF ( ierr.EQ.3 ) GOTO 1010
C else - just continue
C *********************************************************

      IF ( Kprnkr.GT.2 ) THEN
         WRITE (6,320) Knc , Kn
 320     FORMAT (/3X,'KNC=',I3,5X,'KN=',I3)
         WRITE (6,220) ((nam(7),i,j,Mn(i,j),i=1,2),j=1,Kn)
         WRITE (6,270)
         WRITE (6,210) (nam(1),i,Kr(i),i=1,Kn)
         WRITE (6,270)
         WRITE (6,210) (nam(2),i,Kc(i),i=1,Knr)
         WRITE (6,270)
         WRITE (6,210) (nam(3),i,Nnr(i),i=1,Knr)
         WRITE (6,270)
         WRITE (6,210) (nam(4),i,Kr1(i),i=1,Kn1)
         WRITE (6,270)
         WRITE (6,210) (nam(5),i,Kc1(i),i=1,Knr1)
         WRITE (6,270)
         WRITE (6,210) (nam(6),i,Nnr1(i),i=1,Knr1)
         WRITE (6,270)
         WRITE (6,220) ((nam(8),i,j,Mn1(i,j),i=1,2),j=1,Kn1)

         WRITE (6,290) nam(9) , (i,i=1,Kn)
         WRITE (6,270)
         DO i = 1 , Kn
            WRITE (6,350) i , (Wr(i,j),j=1,Kn)
         ENDDO

         WRITE (6,290) nam(10) , (i,i=1,Kn)
         WRITE (6,270)
         DO i = 1 , Kn
            WRITE (6,350) i , (Ws(i,j),j=1,Kn)
         ENDDO

         WRITE (6,270)
         WRITE (6,230) (nam(11),i,W(i),i=1,Kn)
         WRITE (6,270)
         WRITE (6,230) (nam(12),i,W1(i),i=1,Kn1)
         WRITE (6,280)
 280     FORMAT (//1X,120('*'))
      ENDIF
C
      k1 = Kol(1)
      k2 = Kol(2)
      k3 = Kol(3)
      k12 = k1 + k2
      k23 = k2 + k3
      k123 = k12 + k3
      ko0 = k3*Kn
C

C
      IF ( Kprlin.GE.2 ) WRITE (6,400)
      PRINT 400
      DO i1 = 1 , Kn
         wi12pi = W(i1)/(2.D0*pi)
         IF ( Kprlin.GE.2 ) WRITE (6,410) wi12pi
         PRINT 410 , wi12pi
         DO i2 = 1 , k3
            DO i3 = 1 , k3
               IF ( Kprlin.GE.2 ) WRITE (6,421) i3 , i2 , i1 ,
     &              Buffer(ifind2(i3,i2,i1,k3,Kn))
               PRINT 421 , i3 , i2 , i1 , Buffer(ifind2(i3,i2,i1,k3,Kn))
            ENDDO



            IF ( Kprlin.GE.2 ) WRITE (6,415) i2 , i1 ,
     &                                Buffer(ifind1(i2,i1,k3))
            PRINT 415 , i2 , i1 , Buffer(ifind1(i2,i1,k3))
         ENDDO
      ENDDO
C
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     READING AND WRITING OF INITIAL STRESS APPROXIMATIONS
C     FROM A FILE WITH THE EXTENSION ".UIN"
      IF ( npass.NE.1 .OR. Kitu.NE.3 ) GOTO 60
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
      fs1 = 'UIN.UIN'

      OPEN (12,FILE=fs1)
      READ (12,uinit,ERR=55)
      CLOSE (12,STATUS='DELETE')
      GOTO 80
 55   CLOSE (12,STATUS='DELETE')

C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
 60   IF ( npass.LE.1 .OR. Kitu.LT.2 ) THEN
         ko02 = 2*ko0
         DO idint = 1 , ko02
            uin(idint) = 0.0D0
         ENDDO
      ENDIF
C
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 80   IF ( Kitu.EQ.1 ) READ (iin,uinit)
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C
C
C
C *********************************************************
      CALL solve(uin)
C *********************************************************
C
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C       WRITING THE FILE OF INITIAL STRESS APPROXIMATIONS
C                  WITH THE EXTENSION "*UIN"
      IF ( npass.EQ.1 .AND. Kitu.EQ.3 ) OPEN (12,FILE=fs1)
      IF ( npass.EQ.1 .AND. Kitu.EQ.3 ) WRITE (12,uinit)
      IF ( npass.EQ.1 .AND. Kitu.EQ.3 ) Kitu = 2
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C *********************************************************
 1000 CALL restor(uin,yy,vectj,ISIZE_MAXNODE,s)
C *********************************************************
C
C
C
      k123 = Kol(1) + Kol(2) + Kol(3)
      IF ( Kn.GE.2 ) THEN
         DO irc = 2 , Kn
            DO i = 1 , k123
               s(i,irc) = 2.D0*s(i,irc)
C    BRING TO THE ACTUAL VALUE
               s(i,irc) = s(i,irc)/dsqrt(2.D0)
            ENDDO
         ENDDO
      ENDIF
C  90 CONTINUE
C
C     IF(KPRSOL.LE.1) GOTO 150
      DO irc = 1 , Kn
         wircpi = W(irc)/(2.D0*pi)
         IF ( Kprsol.GT.2 ) WRITE (6,430) irc , Mn(1,irc) , Mn(2,irc) ,
     &                             wircpi
         PRINT 430 , irc , Mn(1,irc) , Mn(2,irc) , wircpi
         DO i = 1 , k123
            IF ( Kprsol.GT.2 ) WRITE (6,440) i , s(i,irc)
            PRINT 440 , i , s(i,irc)
         ENDDO
      ENDDO

      CALL wrtraw(ISIZE_MAXNODE,k123,Kn,W,s,pi,infile)

C
C
C
C    IF THE ITERATION LIMIT IS EXCEEDED OR THE METHOD
C    DOES NOT CONVERGE, WE DO NOT CALCULATE QUALITY INDICATORS
 1010 IF ( Itermc.EQ.3 .OR. Itermc.EQ.4 ) THEN
      ENDIF
C
C
C *********************************************************
C    instead of call QUPP - just rudimentary
C       NQUPV=0
C      CALL QUPP
C *********************************************************

C ...
C      IF(KOLVAR.EQ.0)
      WRITE (6,460)
      WRITE (6,450)
 450  FORMAT (20X,'***************** END *************************')
      WRITE (6,460)

 421  FORMAT (2X,'Y(',I3,',',I3,',',I3,')=',E12.5,2X,E12.5)


 210  FORMAT ((1X,7(A5,'(',I2,')=',I3)))
 220  FORMAT ((1X,6(A5,'(',I2,',',I2,')=',I3)))
 230  FORMAT ((1X,5(A5,'(',I2,')=',E12.5)))
 270  FORMAT (1X)
 290  FORMAT (//10X,'MATRIX ',A4/10X,12('_')/10X,20I5)
 350  FORMAT (I8,2X,20I5)
 400  FORMAT (10X,'RESULTS OF Y-MATRICES FORMING   ')
 410  FORMAT (10X,'FREQUENCY',E12.5)
 415  FORMAT (2X,'VECTJ(',I3,I4,')=',E12.5,2X,E12.5)
 420  FORMAT (10E12.5)
 430  FORMAT (5X,' FREQUENCY',I3,' COMBINATION(',I3,',',I3,')',
     &        ' VALUE    ',E13.6)
 440  FORMAT (2X,'U(',I3,')=',E13.6,1X,E13.6)
 460  FORMAT (1X,79('*')/1X,7('*** hbfree '),'**'/1X,79('*'))
 480  FORMAT (10X,'SOLUTION OF NLN EQ. STARTED ')

      END

      SUBROUTINE wrtraw(Isize_maxnode,K123,Kn,W,S,Pi,Infile)
      INTEGER irc , Isize_maxnode , K123 , Kn
      DOUBLE PRECISION W(20) , wircpi , Pi
      DOUBLE COMPLEX S(Isize_maxnode,20)
      CHARACTER*64 Infile , line
      INTEGER i , jrc , inend , ifend , hblnode , lenname , nvar
      CHARACTER*9 rawfile/'         '/
      CHARACTER*13 nodefile/'             '/
      CHARACTER*6 spicenode(20)
      INTEGER dt(8)
      CHARACTER*10 b(3)

      inend = 0
      ifend = 0
      DO i = 1 , 12
         IF ( Infile(i:i).EQ.'.' ) THEN
            inend = i
            GOTO 100
         ENDIF
         IF ( Infile(i:i).NE.' ' ) ifend = i
      ENDDO
 100  IF ( inend.EQ.0 ) THEN
         IF ( ifend.GT.8 ) STOP ' NAME TOO LONG'
         rawfile(1:ifend) = Infile(1:ifend)
         rawfile(ifend+1:ifend+4) = '.raw'
         nodefile(1:ifend) = Infile(1:ifend)
         nodefile(ifend+1:ifend+10) = '.ckt.nodes'
      ELSE
         IF ( inend.GT.9 ) STOP ' NAME TOO LONG'
         rawfile(1:inend-1) = Infile(1:inend-1)
         rawfile(inend:inend+4) = '.raw'
         nodefile(1:inend-1) = Infile(1:inend-1)
         nodefile(inend:inend+9) = '.ckt.nodes'
      ENDIF

      OPEN (15,FILE=trim(nodefile),STATUS='OLD')
      nvar = 0
      DO i = 1 , K123 + 5
         READ (15,'(A20)',END=1500) line
         IF ( i.GT.5 ) THEN
            nvar = nvar + 1
            READ (line,'(I1,TR1,A)') hblnode , spicenode(nvar)
C              SKIP GND NODE
            IF ( hblnode.EQ.0 ) nvar = nvar - 1
         ENDIF
      ENDDO
C  ADD UNMAPPED S ENTRIES AS INTERNAL VARIABLES
      DO i = nvar + 1 , K123
         nvar = nvar + 1
         WRITE (spicenode(nvar),'(A,I2.2)') '#int', i
      ENDDO

 1500 CALL date_and_time(b(1),b(2),b(3),dt)

      OPEN (16,FILE=trim(rawfile))

      WRITE (16,'(A)') 'Title: spice test'
      WRITE (16,170) dt(3) , dt(2) , dt(1) , dt(5) , dt(6) , dt(7)
 170  FORMAT ('Date: ',I2.2,'.',I2.2,'.',I4,' ',I2.2,':',I2.2,':',I2.2)
      WRITE (16,'(A)') 'Plotname:  Harmonic Balance Simulation'
      WRITE (16,'(A)') 'Flags: complex'
      WRITE (16,'(A,I4)') 'No. Variables: ' , nvar + 1
      WRITE (16,'(A,I4)') 'No. Points: ' , 3*Kn
      WRITE (16,'(A)') 'Command:  version hbfree'

      WRITE (16,'(A)') 'Variables:'
      WRITE (16,*) '     0    frequency    frequency'
      DO i = 1 , nvar
         lenname = index(spicenode(i),achar(9)) - 1
         WRITE (16,140) i , spicenode(i)(1:lenname)
 140     FORMAT (3X,I4,'    V(',A,')    voltage')
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

      END
