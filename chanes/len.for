!*==LEN_S.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE len_s
C ***
C *** PROGRAM FOR READING INPUT DATA
C ***
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION b , Epsiw , F , f1 , f1old , f2 , f2old , fs ,
     &                 Pniff
      INTEGER i , ia , iadr , Iapr , idop , iex , Iff , ii , iin , ij ,
     &        ik , ikn , ikol , ikv , ikv1 , ikv2 , il , ip , ipp , ipr
      INTEGER itn , iz , j118 , jf , jj , jjf , jjf1 , kf , Kiff ,
     &        Kitu , Kmni , Kn , Kn1 , Knc , Knniff , Knr , Knr1 , kol ,
     &        kolc , kolpar
      INTEGER kolpol , koluz , kolvn , kpr , Kprgrf , Kprlen , Kprlin ,
     &        Kprnkr , Kprqup , Kprsol , Kprsrt , Kprvar , lenmpo ,
     &        Limerr , Limit , Mglob , Mni , n , n1 , n20
      INTEGER nadr , ndl , nf3 , ni , niadr , nkolc , nm , nm1 , nm20 ,
     &        nmp , nmpoin , nn , nnd , nnda , Nniff , nnn , nnnnn ,
     &        nnnp , nnod , nnode1
      INTEGER nnode3 , nnode4 , nnode5 , nnode6 , nnode7 , nnode8 ,
     &        nnodek , nnodel , nnp , npara1 , npkolc
      INCLUDE 'charint.i'
      INCLUDE 'circuit.i'

      COMMON /fre   / F(2)
      COMMON /blk1  / Kr , Kr1 , Kc , Kc1 , Nnr , Nnr1 , Mn , Mn1
      INTEGER Kr(20) , Kr1(200) , Kc(10) , Kc1(20)
      INTEGER Nnr(10) , Nnr1(20) , Mn(2,20) , Mn1(2,200)
      COMMON /blk2  / Knc , Knr , Kn , Knr1 , Kn1
      COMMON /serv  / Epsiw , Limerr , Kitu
      COMMON /newton/ Epssol , Epsdu , Epsmin , Maxdu , Limit
      DOUBLE PRECISION Epssol , Epsdu , Epsmin , Maxdu
      COMMON /modglb/ Mglob , Iapr/blmni / Mni(2,20) , Kmni
      COMMON /typval/ Typu , Typi
      DOUBLE PRECISION Typu , Typi
      COMMON /print / Kprlen , Kprsrt , Kprnkr , Kprlin , Kprsol ,
     &                Kprvar , Kprgrf , Kprqup
      DOUBLE PRECISION eps , mu , ro , tgd , p(50) , par(1500) ,
     &                 ys(1500)
      COMMON /blka  / A , Ka
      INTEGER*4 A(20,50) , Ka , knot(50)


      CHARACTER*4 it(4) , ne , koh , itob , ipar , is , istr , Fne1 ,
     &            Fne2 , graf , yes , no
      CHARACTER*60 name , nnnn
C
C
      COMMON /bliff / Iff(7,4) , Kiff , Nniff(4,8) , Knniff , Pniff(8) ,
     &                Fne1 , Fne2
C
      DATA it/'    ' , '    ' , '    ' , '    '/ , ip/6/
      DATA nnnn/
     &    '                                                            '
     &    /
      DATA iin/10/ , koh/'END '/ , itob/'    '/ , iex/'EXTR'/ ,
     &     ipar/'    '/
      DATA is/'    '/ , yes/'YES '/ , no/'NO  '/

      NAMELIST /circom/ eps , mu , ro , tgd
      NAMELIST /typ   / it , kol , p
      NAMELIST /elem  / ne , knot , ipr , istr , idop , par , ipar
      NAMELIST /frequ / f1 , f2 , Mn , Kn
      NAMELIST /serv  / name , graf , Epsiw , Limerr , Kitu , Epssol ,
     &   Epsdu , Epsmin , Maxdu , Limit , Kprlen , Kprsrt , Kprnkr ,
     &   Kprlin , Kprsol , Kprvar , Kprqup , Mglob , Iapr , Knc
C
C     ASSIGNMENT OF INITIAL VALUES TO VARIABLES
C
      lenmpo = 500
C             LENMPO - MAXIMUM LENGTH OF THE "MPOINT" ARRAY
      Lennod = 500
C             LENNOD - MAXIMUM LENGTH OF THE "NODEEL" ARRAY
      Lenpar = 2000
C             LENPAR - MAXIMUM LENGTH OF THE "PARAM" ARRAY
      Nnetpr = 1
C             NNETPR - ADDRESS OF THE START OF COMMON PARAMETERS
C                      FOR THE ENTIRE CIRCUIT IN "PARAM"
      Lenntp = 4
C             LENNTP - NUMBER OF PARAMETERS COMMON TO
C                      THE ENTIRE CIRCUIT
C
C     INITIAL "ZEROING" OF ARRAYS MPOINT, NODEEL, PARAM
      DO ii = 1 , 500
         Mpoint(ii) = 0
         Nodeel(ii) = 0
      ENDDO
      DO ii = 1 , 2000
         Param(ii) = 0.0D0
      ENDDO
C
C     ASSIGNMENT OF INITIAL VALUES
      Nnode = 1
C     NNODE - INDEX OF THE FIRST FREE ELEMENT IN "NODEEL"
      nmpoin = 1
C     NMPOIN - INDEX OF THE FIRST FREE ELEMENT IN "MPOINT"
      Nparam = 1
C     NPARAM - INDEX OF THE FIRST FREE ELEMENT IN "PARAM"
      Nmpnt = 0
C     NMPNT - NUMBER OF PROCESSED DIRECT-TYPE ELEMENTS,
C            I.E., NUMBER OF FILLED ROWS IN "MPOINT"
C
C *** INPUT NAMELIST /SERV/
C            NAME   - NAME OF THE CIRCUIT
C            GRAF   - GRAPHIC MODE SETTING
C            LIMERR - MAXIMUM NUMBER OF ERRORS
C            EPSIW  - FREQUENCY COMPARISON ACCURACY
C            LIMIT  - MAXIMUM NUMBER OF ITERATIONS
C            KPRLEN - PRINT MODE FOR INPUT READING SUBROUTINE
C            KPRSRT - PRINT MODE FOR CIRCUIT NODE SORTING SUBROUTINE
C            KPRNKR - PRINT MODE FOR FREQUENCY MESH GENERATION SUBROUTINE
C            KPRLIN - PRINT MODE FOR PROCESSING LINE SUBCIRCUIT SUBROUTINE
C            KPRSOL - PRINT MODE FOR SOLVING NONLINEAR SYSTEM SUBROUTINE
C            KPRVAR - PRINT MODE FOR PARAMETER VARIATION
C            KPRQUP - PRINT MODE FOR QUALITY ASSESSMENT OF INPUT DATA
C            EPSSOL - REQUIRED SOLUTION ACCURACY
C            EPSDU  - MINIMUM POSSIBLE STEP SIZE
C            EPSMIN - USED FOR LOCAL MINIMUM DETECTION,
C                       NOT CONSIDERED A SOLUTION
C            MAXDU  - DETERMINES MAXIMUM STEP SIZE
C            KITU   - SWITCH FOR INITIAL APPROXIMATION TYPE
C            KNC    - NUMBER OF POINTS FOR UNIFORM DFT EXPANSION
C            MGLOB  - DETERMINES GLOBALIZATION METHOD
C                     0 - NO GLOBALIZATION
C                     1 - LINE SEARCH IN NEWTON DIRECTION
C            IAPR   - ENABLE/DISABLE A PRIORI STEP SIZE LIMITATION METHOD
C                     0 - DISABLED / 1 - ENABLED
C     DEFAULT VALUE ASSIGNMENT
      Limerr = 3
      Epsiw = 1.0D-5
      Kitu = 0

      Kprsrt = 0
      Kprnkr = 0
      Kprlin = 0
      Kprsol = 0
      Kprlen = 1
      Kprvar = 1
      Kprqup = 1
      Knc = 32

      name = nnnn
      graf = no

C   DEFAULT SETTINGS SEE SUBROUTINE NEINCK
      Epssol = 1.D-6
      Epsdu = 1.0D-6
      Epsmin = 1.0D-7
      Maxdu = 0.0D0
      Limit = 500

      Iapr = 0
      Mglob = 0

      Typi = 0.D0
      Typu = 0.D0

C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      READ (iin,serv)
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      IF ( Kprlen.GT.0 ) WRITE (ip,6)
 6    FORMAT (//5X,'C I R C U I T  DESCRIPTION '/5X,
     &        '--------------------------'/)
      IF ( Kprlen.GT.0 .AND. name.NE.nnnn ) WRITE (ip,7) name
 7    FORMAT (5X,A60//)
C     *         5X,'--------------------------'//)

      Kprgrf = 0
      IF ( graf.EQ.yes ) Kprgrf = 1

C *** INPUT "NAMELIST/CIRKOM/
C             EPS - DIELECTRIC PERMITTIVITY OF THE SUBSTRATE
C             MU  - MAGNETIC PERMEABILITY OF THE CONDUCTOR
C                   MATERIAL
C             RO  - SPECIFIC RESISTANCE OF THE CONDUCTOR
C                   MATERIAL
C             TGD - TANGENT OF THE LOSS ANGLE IN THE DIELECTRIC
C                   SUBSTRATE
C
C
C     DEFAULT VALUE SETTINGS
C
C
      eps = 9.6D0
      mu = 1.0D0
      ro = 5.7D+07
      tgd = 1.D-04
C
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      READ (iin,circom)
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C
C AT THIS STAGE, 4 GENERAL PARAMETERS (EPS, MU, RO, TGD)
C ARE USED FOR THE ENTIRE CIRCUIT, I.E., LENNTP=4. IN THIS CASE, NNETPR=1.
C
      Param(Nnetpr) = eps
      Param(Nnetpr+1) = tgd
      Param(Nnetpr+2) = mu
      Param(Nnetpr+3) = ro
C
C
      IF ( Kprlen.GE.0 ) WRITE (ip,8) eps , mu , ro , tgd
 8    FORMAT (/2X,'   COMMON PARAMETERS OF CKT :'/5X,
     &        'DIELECTRIC PERMITTIVITY OF SUBSTRATE    ',E12.5/5X,
     &        'MAGNITE PERMITTIVITY OF CONDUCTORS            = ',
     &        E12.5/5X,'SPESIFIC RESISTIVITY OF CONDUCTORS           = '
     &        ,E12.5/5X,'tg(delta) - LOSSES IN DIELECTRIC          = ',
     &        E12.5)
C
C *********************************************************
C *** CALLING THE SUBROUTINE CONTAINING THE ARRAY "A"
      CALL lena
C *********************************************************
C
C
C
C     INCREASE THE COUNTER NPARAM BY THE NUMBER OF COMMON PARAMETERS
      Nparam = Nparam + Lenntp
C
C *** INPUT NAMELIST/TYP/ IT, KOL, P
C          IT - ELEMENT TYPE
C          KOL - NUMBER OF ELEMENTS OF THIS TYPE
C          P  - COMMON PARAMETERS OF THE TYPE (IF THEY EXIST)
C
 10   kol = 1
      DO ipp = 1 , 50
         p(ipp) = 0.D0
      ENDDO
      DO ipp = 1 , 4
         it(ipp) = itob
      ENDDO
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      READ (iin,typ)
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      IF ( it(1).EQ.koh ) THEN
C
C
C
C
C *********************************************************
C *** CALL TO THE FREQUENCY PROCESSING SUBROUTINE
         kpr = Kprlen
         CALL frequen(nf3,kpr)
C *********************************************************
C
         f1old = F(1)
         f2old = F(2)
         f1 = F(1)
         f2 = F(2)
C
C *** READING NAMELIST/FREQU/ F1,F2,MN,KN
C                     F1 - LOWER FREQUENCY
C                     F2 - UPPER FREQUENCY
C                     MN - COMBINATION COEFFICIENT MATRIX
C                     KN - NUMBER OF COMBINATIONS
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
         READ (iin,frequ)
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C
C     ASSIGNING THE VALUES OF F1 AND F2 TO THE ARRAY F
         F(1) = f1
         F(2) = f2
C
C     COPYING THE ARRAY MN TO THE ARRAY MNI FOR ITS RENAMING
         DO j118 = 1 , Kn
            Mni(1,j118) = Mn(1,j118)
            Mni(2,j118) = Mn(2,j118)
         ENDDO
         Kmni = Kn
C
C     PRINTING INFORMATION FROM NAMELIST/FRE/

         IF ( Kprlen.GE.1 .AND. f1old.NE.f1 ) WRITE (ip,120) f1
 120     FORMAT (2X,'FUNDAMENTAL FR.  F1 =',E12.5,',')

         IF ( Kprlen.GE.1 .AND. f1old.EQ.f1 ) WRITE (ip,121) f1 , Fne1
 121     FORMAT (2X,'FUNDAMENTAL FR.  F1 =',E12.5,' ( ELEMENT ',A4,' )')

         IF ( Kprlen.GE.1 .AND. f2old.NE.f2 ) WRITE (ip,122) f2

         IF ( Kprlen.GE.1 .AND. f2old.EQ.f2 .AND. Fne2.NE.is )
     &        WRITE (ip,123) f2 , Fne2
 123     FORMAT (2X,'FUNDAMENTAL FR.  F2 =',E12.5,' ( ELEMENT ',A4,' )')

         IF ( Kprlen.GE.1 .AND. f2old.EQ.f2 .AND. Fne2.EQ.is )
     &        WRITE (ip,122) f2
         IF ( Kprlen.GE.1 .AND. nf3.EQ.1 ) WRITE (ip,124)
 124     FORMAT (//2X,'A T T E N T  O N: MORE THAN 2 FREQUENCIES  !'//
     &           2X,' MAKE CHANGES IN INPUT TASK !')

         IF ( nf3.EQ.1 ) STOP

         IF ( Kprlen.GE.1 ) WRITE (ip,125)
 125     FORMAT (2X,'FREQUENCY GRID   ( F = mi*F1 + ni*F2 ) ')
         DO iz = 1 , Kn
            fs = Mn(1,iz)*F(1) + Mn(2,iz)*F(2)
            IF ( Kprlen.GE.1 ) WRITE (ip,126) iz , Mn(1,iz) , iz ,
     &                                Mn(2,iz) , fs
 126        FORMAT (12X,'m',I2,' =',I3,',   n',I2,' =',I3,'   F = ',
     &              E13.6)
         ENDDO

C     CONSIDERING THE CASE OF TWO ZERO FREQUENCIES
         IF ( F(1).EQ.0.D0 .AND. F(2).EQ.0.D0 ) THEN
            Mn(1,1) = 0
            Mn(2,1) = 0
            Mni(1,1) = 0
            Mni(2,1) = 0
            Kn = 1
            IF ( Kprlen.GE.1 ) WRITE (ip,128)
 128        FORMAT (/2X,'ALL FREQUENCIES ARE = TO ZERO '/2X,
     &              'ONLY DC ANALYSIS !'/)
         ENDIF
C
C     END OF NAMELIST/FRE/ INPUT
C
C
C *** PRINTING INFORMATION FROM NAMELIST/SERV/
         IF ( Kprlen.GE.1 ) WRITE (ip,130) Epsiw , Limit , Kprlen ,
     &                             Kprsrt , Kprnkr , Kprlin , Kprsol ,
     &                             Kprvar , Kprgrf , Kprqup , Epssol ,
     &                             Kitu , Knc
 130     FORMAT (2X,'SERVICE VALUES  :',20X,
     &           'FREQUENCIES CMP PRECISION = ',E12.5/20X,
     &           'MAX NUMBER OF ITERATIONS  = ',I3/20X,
     &           'PRINTINGS     = ',8(I2,',')/20X,
     &           'REQUIRED PREC. OF SOLUTION = ',E12.5/20X,
     &           'KIND OF INITIAL APPROXIMATION = ',I2/20X,
     &           'NUMBER OF POINTS IN FFT = ',I2/)


C ********** CONTROLLED PRINTING OF ARRAYS ******************

         IF ( Kprlen.GT.1 ) THEN

            WRITE (ip,*) 'PRINTS AT OUT OF SUBR LEN'

            WRITE (ip,610) Epsiw , Limit , Kprlen , Kprsrt , Kprnkr ,
     &                     Kprlin , Kprsol , Kprvar , Kprqup , Epssol ,
     &                     Epsdu , Epsmin , Maxdu , Kitu , Knc , Mglob ,
     &                     Iapr
 610        FORMAT (10X,' CONTROL VALUES :',/2X,
     &              'EPSIW (1.0E-5) - FREQUENCIES CMP PRECISION = ',
     &              E12.5/2X,'LIMIT (500) - MAX NUM OF ITERATIONS = ',
     &              I3/2X,
     &              'KPRLEN (1) - PRN LEVEL OF INPUT TASK READING = ',
     &              I2/2X,
     &         'KPRSRT (0) - PRN LEVEL PRG. SORTING OF CIRCUIT NODES = '
     &         ,I2/2X,
     &        'KPRNKR (0) - PRN LEVEL  CONSTRUCTION  FREQUENCY GRIDS = '
     &        ,I2/2X,
     &        'KPRLIN (0) - PRN LEVEL  PROCESSING  LINES SUBCIRCUITS = '
     &        ,I2/2X,
     &       'KPRSOL (0) - PRN LEVEL  SOL.  SYSTEMS LINEAR EQUATIONS = '
     &       ,I2/2X,
     &       'KPRVAR(1) - PRN LEVEL    FOR PARAMETER VARIATION = ',
     &       I2/2X,
     &       'KPRQUP (0) - PRN LEVEL READING INPUT ASSIG. FOR K.P. = ',
     &       I2/2X,'EPSSOL (1.E-6) - REQUIRED PREC. OF SOLUTION   = ',
     &       E12.5/2X,'EPSDU (1.E-6) - MINIMUM ALLOWED STEP = ',
     &       E12.5/2X,'EPSMIN (1.E-7) - USED TO DETECT LOCAL MINIMUM ',
     &       /2X,'                 WHICH NOT A SOLUTION = ',E12.5/2X,
     &       'MAXDU (0.0) - DEF. MAX ALLOWED STEP=',E12.5/2X,
     &       'KITU (0) - KIND OF INITIAL APPROXIMATION = ',I2/2X,
     &       'KNC (32) - NUMBER OF POINTS IN FFT = ',I2/2X,
     &       'MGLOB (0) - GLOBALIZATION APPROACH = ',I2/2X,
     &       '            0 - NO GLOBALIZATION',/2X,
     &       '            1 - LINEAR SEARCH IN NEWTON DIRECTION',/2X,
     &       'IAPR (0) - TURN ON A-PRIORY STEP LIMITATION  = ',I2/2X,
     &       '            0 - OFF 1 - ON   ',/)
C
C *** CONTROLLED PRINTING OF THE ARRAY "MPOINT"
            WRITE (ip,*) 'MPOINT = '
            DO jj = 1 , Nmpnt
               n1 = 20*(jj-1) + 1
               n20 = n1 + 20 - 1
               WRITE (ip,611) n1 , n20 , (Mpoint(n),n=n1,n20)
 611           FORMAT (2X,'(',I3,',',I3,')=',4A4,16(I3))
            ENDDO
C
C *** CONTROLLED PRINTING OF THE ARRAY "NODEEL"
            nnod = Nnode - 1
            DO jj = 1 , nnod
               WRITE (ip,620) jj , Nodeel(jj) , Nodeel(jj)
 620           FORMAT (2X,'NODEEL(',I3,')=',I4,2X,A4)
            ENDDO
            WRITE (ip,633) Nnode
 633        FORMAT (/2X,'FIRST FREE CELL IN  NODEEL - ',I4/)

C *** CONTROLLED PRINTING OF THE ARRAY "PARAM"
            nnp = Nparam - 1
            WRITE (ip,635) (nn,Param(nn),nn=1,nnp)
 635        FORMAT (2X,'PARAM(',I4,')=',E12.5)
            WRITE (ip,637) Nparam
 637        FORMAT (/2X,'FIRST FREE CELL IN PARAM - ',I4/)
C
C *** CONTROLLED PRINTING OF THE ARRAY "IFF"
            DO jf = 1 , Kiff
               WRITE (ip,640) jf , (Iff(jjf1,jf),jjf1=1,7)
 640           FORMAT (2X,'IFF(   ,',I3,')=',4A4,3(1X,I3))
            ENDDO
C
C *** CONTROLLED PRINTING OF THE ARRAY "NNIFF"
            DO ij = 1 , Knniff
               WRITE (ip,665) ij , (Nniff(i,ij),i=1,4)
 665           FORMAT (2X,'NNIFF( ,',I2,')=',A4,', ',I2,', ',I4,', ',I3)
            ENDDO

C *** CONTROLLED PRINTING OF VARIABLES
            WRITE (ip,676) Nmpnt , Nnode , Nparam , Lennod , Lenpar ,
     &                     Nnetpr , Lenntp
 676        FORMAT (2X,'NMPNT=',I2,' NNODE=',I3,1X,'NPARAM=',I4,
     &              ' LENNOD=',I3/2X,'LENPAR=',I4,' NNETPR=',I3,1X,
     &              'LENNTP=',I3)
         ENDIF
C


C *** END OF DIRECT INPUT *********************************************
C
C
         RETURN
      ELSE
         Nmpnt = Nmpnt + 1
C
C     OUTPUT OF INFORMATION FOR NAMELIST/TYP/
         IF ( Kprlen.GT.0 ) WRITE (ip,9) (it(itn),itn=1,4)
 9       FORMAT (/1X,'ELEMENT TYPE - ',4A4)
C
C     DETERMINATION OF THE ROW NUMBER (IA) IN ARRAY "A",
C     WHICH CORRESPONDS TO THE PROCESSED ELEMENT TYPE
         DO ia = 1 , Ka
            IF ( c2i(it(1)).EQ.A(1,ia) .AND. c2i(it(2)).EQ.A(2,ia) .AND.
     &           c2i(it(3)).EQ.A(3,ia) .AND. c2i(it(4)).EQ.A(4,ia) )
     &           GOTO 12
            IF ( ia.EQ.Ka ) GOTO 200
         ENDDO
C
C     IS THERE ENOUGH SPACE IN THE "MPOINT" ARRAY TO
C     ADD A NEW ELEMENT TYPE?
         nnnnn = nmpoin + 20
         IF ( nnnnn.GT.lenmpo ) THEN
C
            WRITE (ip,140)
            WRITE (ip,305)
 305        FORMAT (2X,'NO MORE ROOM FOR INPUT')
C IT IS NECESSARY TO INCREASE THE SIZE OF THE STRING "MPOINT" AND
C ASSIGN A NEW VALUE TO THE VARIABLE LENMPO, COMPARED TO LEN
            GOTO 500
         ENDIF
C
C *** FILLING THE "MPOINT" ARRAY  ************************************
C
 12      nm = nmpoin - 1
C       1-4 - NAME OF ELEMENT TYPE (NO MORE THAN 16 CHARACTERS;
C           THE FIRST FOUR CHARACTERS ARE KEY)
         Mpoint(nm+1) = A(1,ia)
         Mpoint(nm+2) = A(2,ia)
         Mpoint(nm+3) = A(3,ia)
         Mpoint(nm+4) = A(4,ia)
C       5 - NUMBER OF ELEMENTS OF THIS TYPE
         Mpoint(nm+5) = kol
C       6 - FLAG DETERMINING THE VARIETY OF THE ELEMENT
         Mpoint(nm+6) = A(6,ia)
C       7 - LENGTH OF THE DESCRIPTION STRING FOR A SINGLE ELEMENT
C           OF THIS TYPE IN "NODEEL"
         Mpoint(nm+7) = A(7,ia)
C       8 - NUMBER OF NODES IN THE MATHEMATICAL MODEL
C           (FOR LINEAR POLYNOMIALS, I.E., A(6,IA)=2,
C           MPOINT STORES 0 BECAUSE THIS VALUE IS DETERMINED
C           IN "NODEEL" SEPARATELY FOR EACH LINEAR POLYNOMIAL)
         IF ( A(6,ia).EQ.2 ) Mpoint(nm+8) = 0
         IF ( A(6,ia).NE.2 ) Mpoint(nm+8) = A(8,ia)
C       9 - NUMBER OF INTERNAL NODES AUTOMATICALLY NUMBERED
C           WHEN DESCRIBING THE CONNECTION NODES OF THE ELEMENT
         Mpoint(nm+9) = A(9,ia)
C      10 - ADDRESS OF THE START OF ELEMENT DESCRIPTIONS OF THIS
C           TYPE IN "NODEEL"
         Mpoint(nm+10) = Nnode
C      11 - NUMBER OF PARAMETERS COMMON TO THIS ELEMENT TYPE
         Mpoint(nm+11) = A(11,ia)
C      12 - ADDRESS OF THE START OF THE LIST OF PARAMETERS
C           COMMON TO THIS ELEMENT TYPE IN THE "PARAM" ARRAY
         Mpoint(nm+12) = Nparam
C      13 - NUMBER OF PARAMETERS DESCRIBING THE MATHEMATICAL
C           MODEL OF THE ELEMENT (ADDITIONAL ANALOGUE POSITION N8)
         IF ( A(6,ia).EQ.2 ) Mpoint(nm+13) = A(13,ia)
         IF ( A(6,ia).NE.2 ) Mpoint(nm+13) = A(13,ia)
C      14 - ADDRESS OF THE START OF THE PARAMETER LIST FOR
C           DIFFERENT VARIANTS OF MATHEMATICAL MODELS
C           OF THIS ELEMENT TYPE IN "PARAM"
         Mpoint(nm+14) = A(11,ia) + Nparam
C      15 - FLAG FOR PARAMETER VARIATION OF ELEMENTS
C           FILLED WHEN A VARIABLE ELEMENT ARRIVES.
C           INITIAL VALUE IS 0, STORED IN ARRAY "A"
         Mpoint(nm+15) = A(15,ia)
C      16 - ADDRESS OF THE LIST OF THE FULL NAME OF THE ELEMENT TYPE
C           AND ITS OUTPUTS (NOT USED YET)
         Mpoint(nm+16) = A(16,ia)
C      17 - ADDRESS OF THE START OF THE LIST OF FULL NAMES
C           OF ELEMENT PARAMETERS (NOT USED YET)
         Mpoint(nm+17) = A(17,ia)
C      18 - MAXIMUM NUMBER OF INPUT QUANTITIES WHEN CALCULATING
C           THE CURRENT OF A NONLINEAR ELEMENT
C           (FOR OTHER TYPES, THIS PARAMETER = 0)
         Mpoint(nm+18) = A(18,ia)
C      19 - MAXIMUM NUMBER OF OUTPUT CURRENTS DURING THE ANALYSIS
C           OF A NONLINEAR ELEMENT
C           (FOR OTHER TYPES, THIS PARAMETER = 0)
         Mpoint(nm+19) = A(19,ia)
C      20 - ADDITIONAL FLAG (RESERVED)
         Mpoint(nm+20) = A(20,ia)
C
C     CONTROL PRINTING OF THE "MPOINT" ARRAY ROW
         nm1 = nm + 1
         nm20 = nm + 20
         IF ( Kprlen.GT.2 ) WRITE (ip,14) (Mpoint(nn),nn=nm1,nm20)
 14      FORMAT ('  MPOINT= ',4A4,16(1X,I3)/)
C
         nmpoin = nm + 20 + 1
C
C *** FORMULATION OF THE "PARAM" ARRAY FOR GENERAL PARAMETERS TYPE *****
C
         kolc = A(11,ia)
         DO ikol = 1 , kolc
            nkolc = Nparam - 1 + ikol
            IF ( Kprlen.GE.3 ) WRITE (ip,1513) kolc , ikol , Nparam ,
     &                                nkolc , p(ikol)
 1513       FORMAT (2X,'KOLC=',I3,' IKOL=',I3,' NPARAM=',I4,' NKOLC=',
     &              I4,' P(IKOL)=',E12.5)
            Param(nkolc) = p(ikol)
            IF ( Kprlen.GT.3 ) THEN
               WRITE (ip,1515) nkolc , Param(nkolc)
 1515          FORMAT (2X,'COMM PARAM(',I4,')=',E12.5)
            ENDIF
         ENDDO
C     OUTPUT INFORMATION ABOUT GENERAL PARAMETERS TYPE
         npkolc = Nparam - 1 + kolc
         DO n = Nparam , npkolc
            b = 0.1D-15
            IF ( Param(n).GE.b .OR. Param(n).LE.-b ) GOTO 17
         ENDDO
         GOTO 19
      ENDIF
 17   IF ( Kprlen.GE.1 ) WRITE (ip,18) (Param(nkolc),nkolc=Nparam,npkolc
     &                                 )
 18   FORMAT (22X,'COMM PARAM=',3(E12.5,', ')/34X,3(E12.5,', ')/)
 19   Nparam = Nparam + kolc
C
C *** PROCESSING OF ALL ELEMENTS OF ONE TYPE ***************************
C
C     FIRST, WE CHECK KOL, GIVEN IN THE INPUT DATA
      IF ( kol.LT.1 ) THEN
C
         WRITE (ip,140)
         WRITE (ip,375) kol
 375     FORMAT (2X,'WRONG NUMBER OF ELEMENTS (KOL=',I4,')'/2X,
     &           'IN &TYP')
         GOTO 500
      ELSE
C
C
         DO ik = 1 , kol
C
C     INITIAL VALUE ASSIGNMENT FOR ELEMENTS
C     OF THE SAME TYPE
            DO i = 1 , 50
               knot(i) = -77777
            ENDDO
            DO i = 1 , 1500
               par(i) = 0.5D+35
               ys(i) = 0.5D+35
            ENDDO
            ipr = 0
            istr = is
            idop = 0
            ipar = itob
C
C *** INPUT NAMELIST /ELEM/ NE, KNOT, IPR, ISTR, IDOP, PAR, IPAR
C          NE   - SYMBOLIC DESIGNATION OF THE ELEMENT
C                 (4 CHARACTERS)
C          KNOT - NODE NUMBERS OF THE SCHEME
C          IPR  - PARAMETER VARIATION INDICATOR OF THE ELEMENT
C          ISTR - TYPE OF SEMICONDUCTOR STRUCTURE
C                 (ONLY FOR NONLINEAR ELEMENTS)
C          IDOP - GROUP MEMBERSHIP NUMBER
C          PAR  - PARAMETERS SPECIFIC TO A PARTICULAR ELEMENT
C          IPAR - DESIGNATION OF AN ELEMENT WHOSE
C                 PARAMETERS ARE IDENTICAL TO THOSE
C                 OF ELEMENT NE
C
C
C  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            READ (iin,elem)
C  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C
            IF ( Kprlen.GE.1 ) WRITE (ip,21) ne
 21         FORMAT (6X,'ELEMENT - ',A4,':')
C
C     CHECKING THE DESIGNATIONS OF THE PROCESSED ELEMENTS
            DO il = 1 , Nmpnt
               IF ( c2i(ne).EQ.Mpoint(1+(Nmpnt-1)*20) ) GOTO 240
            ENDDO
            DO i = 1 , Nnode
               IF ( c2i(ne).EQ.Nodeel(i) ) GOTO 240
            ENDDO
C
C     DETERMINATION OF THE NUMBER OF NODES (KOLPOL)
C
            DO i = 1 , 50
               IF ( knot(i).GT.-77770 ) kolpol = i
               IF ( kolpol.EQ.i .AND. Kprlen.GE.3 ) WRITE (ip,24) i ,
     &              knot(i) , kolpol
 24            FORMAT (2X,' KNOT(',I3,')=',I6,'       KOLPOL=',I3)
            ENDDO
C
C  &&&&  INPUT OF PARAMETERS SPECIFIED VIA CHANNEL 8  &&&&&&&&&&
            IF ( c2i(ipar).NE.iex ) GOTO 27
            OPEN (8,FILE='SNTPY',STATUS='OLD')
            READ (8,*,END=259,ERR=259) (ys(i),i=1,27)
 259        IF ( Kprlen.GE.3 ) WRITE (ip,25) (i,ys(i),i=1,27)
 25         FORMAT (2X,'YS(',I4,')=',E13.6)
            DO i = 1 , 1500
               par(i) = ys(i)
            ENDDO
C  &&&&  END OF INPUT  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     DETERMINATION OF THE NUMBER OF PARAMETERS (KOLPAR)
C     RECEIVED AT THE INPUT
 27         DO i = 1 , 1500
               IF ( par(i).LT.0.5D+34 ) kolpar = i
               IF ( kolpar.EQ.i .AND. Kprlen.GE.3 ) WRITE (ip,28) i ,
     &              par(i) , kolpar
 28            FORMAT (2X,'  PAR(',I4,')=',E12.5,' KOLPAR=',I4)
            ENDDO
C
C
C     DEPENDING ON THE TYPE OF ELEMENT, THE INITIAL
C     DATA IS PROCESSED DIFFERENTLY
            IF ( A(6,ia).NE.1 ) THEN
               IF ( A(6,ia).EQ.2 ) THEN
C
C === FORMATION OF THE "NODEEL" STRING FOR LINEAR
C === MULTIPLEXERS
C
C     IS THERE ENOUGH SPACE IN THE "NODEEL" ARRAY?
                  nnnnn = Nnode + 8 + kolpol
                  IF ( nnnnn.GT.Lennod ) GOTO 350
                  Nodeel(Nnode) = c2i(ne)
C       NE - SYMBOLIC DESIGNATION OF THE ELEMENT
C     THE NUMBER OF TERMINALS IS STORED IN "A(8,IA)"
                  IF ( kolpol.NE.A(8,ia) ) GOTO 320
                  Nodeel(Nnode+1) = kolpol
C     THE NUMBER OF PARAMETERS DESCRIBING THE ELEMENT
C     AT A SINGLE FREQUENCY POINT IS STORED IN "A(13,IA)"
                  IF ( kolpar.NE.A(13,ia) ) GOTO 330
                  Nodeel(Nnode+2) = kolpar
C     INDICATOR DESCRIBING THE CLASS OF THE ELEMENT
C     (ACTIVE OR PASSIVE)
C     STORAGE LOCATION IN "A(20,IA)"
                  Nodeel(Nnode+3) = A(20,ia)
C     CONNECTION NODES
                  DO ikn = 1 , kolpol
                     Nodeel(Nnode+3+ikn) = knot(ikn)
                  ENDDO
                  ikn = kolpol
C     SYMBOLIC DESIGNATION OF THE PARAMETER SET
C     (AT THIS STAGE, THE SYMBOLIC DESIGNATION
C     OF THE ELEMENT IS RECORDED, I.E., "NE")
                  Nodeel(Nnode+4+ikn) = c2i(ne)
C     THE ADDRESS OF THE PARAMETER SET IN "PARAM" IS CALCULATED
                  iadr = Mpoint(nm+14) + (ik-1)*Mpoint(nm+13)
                  iadr = Nparam
                  Nodeel(Nnode+5+ikn) = iadr
                  Nodeel(Nnode+6+ikn) = ipr
C       IPR - VARIATION INDICATOR. MAY BE ABSENT
C            IN THE "NAMELIST" LIST
                  IF ( ipr.EQ.1 ) Mpoint(nm+15) = 1
                  Nodeel(Nnode+7+ikn) = idop
C       IDOP - GROUP MEMBERSHIP NUMBER. MAY
C              BE ABSENT IN THE "NAMELIST" LIST.
C
C     PRINTING INFORMATION ABOUT "NAMELIST/ELEM/"
C
                  IF ( Kprlen.GE.1 ) WRITE (ip,44)
     &                 (knot(jj),jj=1,kolpol)
 44               FORMAT (22X,'NODES   .=',10(I2,','))
C
C     CONTROL PRINTING OF THE "NODEEL" STRING
                  IF ( Kprlen.GT.2 ) THEN
                     nnode3 = Nnode + 3 + kolpol
                     WRITE (ip,45) Nnode
 45                  FORMAT (2X,'NODEEL(',I3,')=')
                     WRITE (ip,46) (Nodeel(nn),nn=Nnode,nnode3)
 46                  FORMAT (15X,A4,10(1X,I4))
                     nnode4 = Nnode + 4 + kolpol
                     nnode7 = Nnode + 7 + kolpol
                     nnode8 = Nnode + 8 + kolpol
                     WRITE (ip,47) (Nodeel(nn),nn=nnode4,nnode7) ,
     &                             nnode8
 47                  FORMAT (17X,A4,3(1X,I4),'  NNODE CB.=',I3)
                  ENDIF
C
                  Nnode = Nnode + 8 + kolpol
C
                  GOTO 60
               ELSEIF ( A(6,ia).EQ.3 ) THEN
C
C
C
C
C === FORMATION OF THE "NODEEL" STRING FOR NONLINEAR
C === ELEMENTS
C
C
C     IS THERE ENOUGH SPACE IN THE ARRAY "NODEEL"
                  nnnnn = Nnode + 7 + kolpol + A(9,ia)
                  IF ( nnnnn.GT.Lennod ) GOTO 350
C
                  Nodeel(Nnode) = c2i(ne)
C       NE - SYMBOLIC DESIGNATION OF THE ELEMENT
C     THE TOTAL NUMBER OF NODES IS STORED IN "A(8,IA)"
C     THE NUMBER OF EXTERNAL NODES IS EQUAL TO THE DIFFERENCE
C     BETWEEN THE TOTAL NUMBER OF NODES AND THE NUMBER OF INTERNAL NODES
                  koluz = A(8,ia) - A(9,ia)
C     CHECKING THE NUMBER OF EXTERNAL NODES
                  IF ( kolpol.NE.koluz ) GOTO 320
C     NODE NUMBERS TO WHICH THE EXTERNAL
C     ELECTRODES OF THE ELEMENT ARE CONNECTED
                  DO ikn = 1 , koluz
                     Nodeel(Nnode+ikn) = knot(ikn)
                  ENDDO
C     THE NUMBER OF INTERNAL NODES IS STORED IN "A(9,IA)"
                  kolvn = A(9,ia)
C     IF KOLVN=0, THEN INTERNAL NODES ARE NOT CONSIDERED
                  IF ( kolvn.NE.0 ) THEN
C     NUMBERS OF INTERNAL NODES
                     ikv1 = 1 + koluz
                     ikv2 = koluz + kolvn
                     DO ikv = ikv1 , ikv2
                        Nodeel(Nnode+ikv) = 100 + ikv
                        knot(ikv) = 100 + ikv
                     ENDDO
                  ENDIF
                  IF ( kolvn.EQ.0 ) ikv = koluz
                  IF ( kolvn.NE.0 ) ikv = ikv2
C     SYMBOLIC DESIGNATION OF THE PARAMETER SET
C     (AT THIS STAGE, THE SYMBOLIC DESIGNATION
C     OF THE ELEMENT TYPE IS STORED, I.E., "IT(1)", "IT(2)")
                  Nodeel(Nnode+1+ikv) = c2i(it(1))
                  Nodeel(Nnode+2+ikv) = c2i(it(2))
C     THE ADDRESS OF THE PARAMETER SET IN "PARAM" IS CALCULATED
                  iadr = Nparam
                  Nodeel(Nnode+3+ikv) = iadr
                  Nodeel(Nnode+4+ikv) = ipr
C       IPR - VARIATION FLAG. MAY BE ABSENT IN
C             THE "NAMELIST" LIST
                  IF ( ipr.EQ.1 ) Mpoint(nm+15) = 1
C     IF(ISTR.EQ.IS) GOTO 260  *** No longer needed 3.11.89 ***
                  Nodeel(Nnode+5+ikv) = c2i(istr)
C       ISTR - TYPE OF SEMICONDUCTOR STRUCTURE
C              (ONLY FOR NONLINEAR ELEMENTS)
                  Nodeel(Nnode+6+ikv) = idop
C       IDOP - NODE GROUP MEMBERSHIP FLAG. MAY
C              BE ABSENT IN THE "NAMELIST" LIST
C      CHECKING THE NUMBER OF PARAMETERS
                  IF ( kolpar.NE.A(13,ia) ) GOTO 330
C
C
C     PRINTING INFORMATION FOR "NAMELIST/ELEM/"
C
                  IF ( Kprlen.GE.1 ) WRITE (ip,53) (knot(jj),jj=1,koluz)
 53               FORMAT (22X,'NODES    =',10(I2,','))
C
                  IF ( Kprlen.GT.1 .AND. kolvn.GT.0 ) WRITE (ip,52)
     &                 (knot(jj),jj=ikv1,ikv2)
 52               FORMAT (22X,'INTERN NODS=',10(I3,',')/34X,10(I3,','))
C
                  IF ( Kprlen.GE.1 .AND. istr.NE.is ) WRITE (ip,54) istr
 54               FORMAT (22X,'STRUCTURE - ',A4)
C
C     CONTROLLED PRINTING OF THE "NODEEL" LINE
                  IF ( Kprlen.GT.2 ) THEN
                     nnodek = Nnode + koluz
                     WRITE (ip,555) Nnode
 555                 FORMAT (2X,'NODEEL(',I3,')=')
                     WRITE (ip,57) (Nodeel(nn),nn=Nnode,nnodek)
 57                  FORMAT (15X,A4,10(1X,I4))
                     IF ( kolvn.NE.0 ) THEN
                        nnodek = Nnode + ikv1
                        nnodel = Nnode + ikv2
                        WRITE (ip,58) (Nodeel(nn),nn=nnodek,nnodel)
 58                     FORMAT (15X,10(1X,I4))
                     ENDIF
                     nnode1 = Nnode + 1 + ikv
                     nnode6 = Nnode + 6 + ikv
                     nnode7 = Nnode + 7 + ikv
                     WRITE (ip,59) (Nodeel(nn),nn=nnode1,nnode6) ,
     &                             nnode7
 59                  FORMAT (15X,A4,A4,2(1X,I4),1X,A4,1X,I4,
     &                       ' NNODE CB.=',I3)
                  ENDIF
C
                  Nnode = Nnode + 7 + ikv
                  GOTO 60
               ELSEIF ( A(6,ia).EQ.5 ) THEN
C
C
C === FORMATION OF THE "NODEEL" STRING FOR ELEMENTS
C === SPECIFIED USING Y- OR S-MATRICES
C
                  nnnnn = Nnode + 8 + kolpol
                  IF ( nnnnn.GT.Lennod ) GOTO 350
                  Nodeel(Nnode) = c2i(ne)
C       NE - SYMBOLIC DESIGNATION OF THE ELEMENT
C     THE NUMBER OF TERMINALS FOR EACH ELEMENT IS DETERMINED AUTOMATICALLY.
                  Nodeel(Nnode+1) = kolpol
C     THE NUMBER OF PARAMETERS IS DETERMINED AUTOMATICALLY FOR EACH ELEMENT
C     AUTOMATICALLY CHECKED
                  Nodeel(Nnode+2) = kolpar
C  ?? INDICATOR, DESCRIBING THE TYPE OF MATRIX (Y OR S) ??
                  Nodeel(Nnode+3) = A(20,ia)
C     CONNECTION NODES
                  DO ikn = 1 , kolpol
                     Nodeel(Nnode+3+ikn) = knot(ikn)
                  ENDDO
                  ikn = kolpol
C     SYMBOLIC DESIGNATION OF THE PARAMETER SET
C     (AT THIS STAGE, THE SYMBOLIC DESIGNATION
C      OF THE ELEMENT IS RECORDED, I.E., "NE")
                  Nodeel(Nnode+4+ikn) = c2i(ne)
C     THE ADDRESS OF THE PARAMETER SET IN "PARAM"
                  Nodeel(Nnode+5+ikn) = Nparam
C !!! VARIATION INDICATOR (OPTIONAL OPERATION)
                  Nodeel(Nnode+6+ikn) = ipr
C !!! IF A VARIATION INDICATOR IS SET, THEN IT IS NECESSARY
C !!! TO UPDATE ADDITIONAL INFORMATION IN "MPOINT"
                  IF ( ipr.EQ.1 ) Mpoint(nm+15) = 1
C !!! GROUP MEMBERSHIP NUMBER (OPTIONAL OPERATION)
                  Nodeel(Nnode+7+ikn) = idop
C
C     PRINTING INFORMATION ABOUT "NAMELIST/ELEM/"
C
                  IF ( Kprlen.GE.1 ) WRITE (ip,8044)
     &                 (knot(jj),jj=1,kolpol)
 8044             FORMAT (22X,'NODES    =',10(I2,','))
C
C     CONTROL PRINTING OF THE "NODEEL" STRING
                  IF ( Kprlen.GT.2 ) THEN
                     nnode3 = Nnode + 3 + kolpol
                     WRITE (ip,8045) Nnode
 8045                FORMAT (2X,'NODEEL(',I3,')=')
                     WRITE (ip,8046) (Nodeel(nn),nn=Nnode,nnode3)
 8046                FORMAT (15X,A4,10(1X,I3))
                     nnode4 = Nnode + 4 + kolpol
                     nnode7 = Nnode + 7 + kolpol
                     nnode8 = Nnode + 8 + kolpol
                     WRITE (ip,8047) (Nodeel(nn),nn=nnode4,nnode7) ,
     &                               nnode8
 8047                FORMAT (17X,A4,3(1X,I4),'  NNODE CB.=',I3)
                  ENDIF
C
                  Nnode = Nnode + 8 + kolpol
C
                  GOTO 60
               ENDIF
            ENDIF
C
C === FORMATION OF THE "NODEEL" STRING FOR LINEAR DIODES
C
C       WHETHER THE PLACE IN THE "NODEEL" STRING IS SUFFICIENT
C       FOR THE RECORD OF A NEW ELEMENT
            nnnnn = Nnode + 5
            IF ( nnnnn.GT.Lennod ) GOTO 350
C
C     CHECKING THE NUMBER OF CONNECTION NODES AND PARAMETERS
            IF ( kolpol.NE.2 ) GOTO 320
            IF ( kolpar.NE.A(13,ia) ) GOTO 330
C
            Nodeel(Nnode) = c2i(ne)
            IF ( nnnnn.GT.Lennod ) GOTO 350
            Nodeel(Nnode+1) = knot(1)
            Nodeel(Nnode+2) = knot(2)
C       KNOT(1),KNOT(2) - CONNECTED NODES
            Nodeel(Nnode+3) = ipr
C       IPR - VARIATION INDICATOR. MAY BE ABSENT
C             IN THE "NAMELIST" LIST
            IF ( ipr.EQ.1 ) Mpoint(nm+15) = 1
            Nodeel(Nnode+4) = idop
C       IDOP - GROUP MEMBERSHIP NUMBER. MAY
C              BE ABSENT IN THE "NAMELIST" LIST.
C
            IF ( Kprlen.GE.1 ) WRITE (ip,33) knot(1) , knot(2)
 33         FORMAT (22X,'NODES CONN.=',I2,',',I2)
C
C     CONTROL PRINTING OF THE "NODEEL" STRING
            IF ( Kprlen.GT.2 ) THEN
               nnode4 = Nnode + 4
               nnode5 = Nnode + 5
               WRITE (ip,34) Nnode
 34            FORMAT (2X,'NODEEL(',I3,')=')
               WRITE (ip,35) (Nodeel(n),n=Nnode,nnode4) , nnode5
 35            FORMAT (15X,A4,4(2X,I4),2X,'NNODE CB.=',I3)
            ENDIF
C
            Nnode = Nnode + 5
C
C
C     MESSAGE ABOUT GROUP MEMBERSHIP
C     (PRINTED IF THE ELEMENT IS EXTERNAL)
 60         IF ( Kprlen.GE.1 .AND. idop.EQ.1 ) WRITE (ip,660)
 660        FORMAT (22X,'ELEMENT -EXTERNAL')
C
C === IN CASE OF OUTPUTTING A LIST OF PARAMETERS AND
C === THE PRESENCE OF A LINK TO ONE OF THE PREVIOUS
C === ELEMENTS:
            IF ( ipar.NE.itob ) THEN
               IF ( Kprlen.GE.3 ) WRITE (ip,670) ipar , iex
 670           FORMAT (2X,'IPAR=',A4,' IEX=',A4)
               IF ( c2i(ipar).NE.iex ) THEN
C     (FOR NONLINEAR DOUBLE-CURRENT ELEMENTS, THIS IS NOT ALLOWED)
                  IF ( A(6,ia).EQ.1 ) GOTO 280
C     (SAME AS FOR ELEMENTS DEFINED BY THE MATRIX)
                  IF ( A(6,ia).EQ.4 ) GOTO 340
C
C     DETERMINING THE ADDRESS OF THE START OF THE DESCRIPTION OF
C     ELEMENTS OF THE GIVEN TYPE IN "NODEEL"
                  nmp = nmpoin - 21
                  nadr = Mpoint(nmp+10)
C     DETERMINING THE LENGTH OF THE DESCRIPTION LINE FOR ONE
C     ELEMENT OF THE GIVEN TYPE
                  ndl = A(7,ia)
C     CHECKING FOR THE PRESENCE OF THE ELEMENT, SPECIFIED
C     IN THE LINK
                  DO nn = 1 , kol
                     nnd = nadr + ndl*(nn-1)
                     IF ( c2i(ipar).EQ.Nodeel(nnd) ) GOTO 62
                     IF ( Kprlen.GE.3 ) WRITE (ip,680) nn , kol
 680                 FORMAT (2X,'NN=',I4,' KOL=',I4)
                     IF ( nn.EQ.kol ) GOTO 220
                  ENDDO
C     DETERMINING THE ADDRESS OF THE START OF THE PARAMETER LIST
C     IN THE ARRAY "PARAM", LINKED TO THE ARRAY
C     "NODEEL", FOR THE ELEMENT SPECIFIED IN THE LINK.
 62               IF ( A(6,ia).EQ.2 ) niadr = nnd + 5 + A(8,ia)
                  IF ( A(6,ia).EQ.3 ) niadr = nnd + 3 + A(8,ia)
                  iadr = Nodeel(niadr)
C     WRITING INTO THE ARRAY "NODEEL" THE NEW VALUE
C     OF THE ADDRESS FOR THE PROCESSED ELEMENT.
                  nnn = Nnode - ndl
                  IF ( A(6,ia).EQ.2 ) nnda = nnn + 5 + ikn
                  IF ( A(6,ia).EQ.3 ) nnda = nnn + 3 + ikv
                  Nodeel(nnda) = iadr
C     PRINT 6666,NND,NIADR,IADR,NNN,NNODE,IKV
C6666 FORMAT(2X,'NND,NIADR,IADR,NNN,NNODE,IKV=',6I4)
                  IF ( Kprlen.GE.2 ) WRITE (ip,6262) nnda , iadr
 6262             FORMAT (2X,'NEW VALUE  "NODEEL(',I4,')"=',I4)
                  IF ( Kprlen.GE.1 ) WRITE (ip,66) ipar
 66               FORMAT (22X,'PARAMETERS ARE IDENTICAL TO PAR',
     &                    'OF ELEM. ',A4)
                  GOTO 100
               ENDIF
            ENDIF
C
C *** FORMATTING THE ARRAY "PARAM"  ***************************************
C
C     PRELIMINARY CHECK. VERIFY IF THE LOCATION IN "PARAM" IS VALID

            nnnnn = Nparam + kolpar
            IF ( nnnnn.GT.Lenpar ) GOTO 360

            f1 = 0.D0
            f2 = 0.D0

            DO jf = 1 , Kiff
               IF ( Kprlen.GE.3 ) WRITE (ip,69) jf
 69            FORMAT (2X,'   JF=',I3)
               jjf = -1
               IF ( Iff(1,jf).EQ.c2i(it(1)) .AND. Iff(2,jf)
     &              .EQ.c2i(it(2)) .AND. Iff(3,jf).EQ.c2i(it(3)) .AND.
     &              Iff(4,jf).EQ.c2i(it(4)) ) THEN
                  Iff(6,jf) = jf
                  jjf = jf
                  IF ( Kprlen.GE.3 ) WRITE (ip,71) jf ,
     &                 (Iff(jjf1,jf),jjf1=1,7)
 71               FORMAT (2X,'IFF(   ,',I3,')=',4A4,3(1X,I3))
C      IF(IFF(1,JF).EQ.IT(1).AND.IFF(2,JF).EQ.IT(2).AND.
C     *   IFF(3,JF).EQ.IT(3).AND.IFF(4,JF).EQ.IT(4)) JJF=JF
               ENDIF

               DO ikol = 1 , kolpar
                  ni = Nparam - 1 + ikol
                  Param(ni) = par(ikol)
                  IF ( jjf.NE.-1 ) THEN
                     DO kf = 1 , Kiff
C     IF(KPRLEN.GT.3) WRITE(IP,76) IFF(6,KF),KF,IKOL,IFF(5,KF)
C  76 FORMAT(2X,'IFF(6,KF)=',I3,' KF=',I3,' IKOL=',I3,' IFF(5,KF)=',I3)
                        IF ( Iff(6,kf).EQ.jjf .AND. ikol.EQ.Iff(5,kf) )
     &                       THEN
                           Knniff = Knniff + 1
                           Nniff(1,Knniff) = c2i(ne)
                           Nniff(2,Knniff) = jjf
                           Nniff(3,Knniff) = ni
                           IF ( Kprlen.GE.3 ) WRITE (ip,80) Knniff ,
     &                          ne , Knniff , jjf , Knniff , ni , kf ,
     &                          ni , Param(ni)
 80                        FORMAT (2X,'NNIFF=(1',I3,')= NE  = ',A4/2X,
     &                             'NNIFF=(2',I3,')= JJF =',I3/2X,
     &                             'NNIFF=(3',I3,')= NI  =',I3/2X,
     &                             ' KF=',I3,' PARAM(',I4,')=',E12.5)

                           GOTO 85
                        ENDIF
                     ENDDO
                  ENDIF
 85            ENDDO
            ENDDO
C     CONTROLLED PRINTING OF "PARAM"
            npara1 = Nparam - 1 + kolpar
            IF ( Kprlen.GE.1 ) WRITE (ip,90)
     &                                (Param(nn),nn=Nparam,npara1)
 90         FORMAT (22X,'PARAMETRS=',4(E12.5,',')/(32X,4(E12.5,',')))
C    *      (32X,4(E12.5,','))/(32X,4(E12.5,','))/(32X,4(E12.5,',))/)
            IF ( Kprlen.GE.3 ) WRITE (ip,95) Nparam , npara1
 95         FORMAT (22X,'IN ARRAY  "PARAM" FILLED ','POS. SINCE',I4,
     &              ' TO ',I4)
C
            Nparam = Nparam + kolpar
            nnnp = Nparam
 100     ENDDO
C
C *** END OF NAMELIST/TYP/ AND NAMELIST/ELEM/ INPUT
C
         GOTO 10
      ENDIF
C
 200  WRITE (ip,140)
      WRITE (ip,210) it
 210  FORMAT (2X,'TYPE"',4A4,'" IS ABSENT IN  ','LIBRARY OF ELEMENTS')
      GOTO 500
C
 220  WRITE (ip,140)
      WRITE (ip,230) ipar
 230  FORMAT (2X,'ELEMENT "',A4,'" IS ABSENT IN CKT')
      GOTO 500
C
 240  WRITE (ip,140)
      WRITE (ip,250) ne
 250  FORMAT (2X,'NAME OF ELEMEN (NE="',A4,'") '/2X,
     &        'SHOULD NOT BE EQUAL TO NAME OF '/2X,
     &        'TYPE OR NAME OF OTHER ELEMENT ')
      GOTO 500
C
      WRITE (ip,140)
      WRITE (ip,270) ne
 270  FORMAT (2X,'FOR ELEMENT ',A4,' IS NOT DEF.'/2X,
     &        'SEMICONDUCTOR SRUCTURE ')
      GOTO 500
C
 280  WRITE (ip,140)
      WRITE (ip,290)
 290  FORMAT (2X,'LINEAR 2-POLES DO NOT ALLOW IPAR VALUE IN  &ELEM')
      GOTO 500
C
 320  WRITE (ip,140)
      WRITE (ip,325) ne , kolpol
 325  FORMAT (2X,'WRONG NODES NUMBER FOR ELEMENT ',A4,'(',I3,')')
      GOTO 500
C
 330  WRITE (ip,140)
      WRITE (ip,335) ne , kolpar
 335  FORMAT (2X,'WRONG PARAMETERS NUMBER FOR ELEMENT ',A4,'(',I4,')')
      GOTO 500
C
 340  WRITE (ip,140)
      WRITE (ip,345)
 345  FORMAT (2X,'CANT USE IPAR IN &ELEM FOR TABLE DEFINED COMPONENT')
      GOTO 500
C
 350  WRITE (ip,140)
      WRITE (ip,355)
 355  FORMAT (2X,'NO MORE ROOM IN NODEEL')
C IT IS NECESSARY TO INCREASE THE SIZE OF THE STRING NODEEL AND ASSIGN A NEW VALUE
C TO THE VARIABLE LENNOD, COMPARED TO LEN
      GOTO 500
C
 360  WRITE (ip,140)
      WRITE (ip,365)
C IT IS NECESSARY TO INCREASE THE SIZE OF THE STRING "PARAM" AND ASSIGN A NEW VALUE
C TO THE VARIABLE LENPAR
 365  FORMAT (2X,'NO MORE ROOM IN PARAM')
C
C
C
 500  WRITE (ip,140)
C
      STOP
 122  FORMAT (2X,'FUNDAMENTAL FR.  F2 =',E12.5,',')
C
C
C *** ERROR DIAGNOSTICS *********************************************
 140  FORMAT (/1X,78('-')//)
C
      END
