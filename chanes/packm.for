!*==PACK1.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c





      SUBROUTINE pack1(Nrec,Y,Vj,Isize_maxnode)
C
C     this routine just stores Y(part) and VJ arrays and allows to get
C     them when necesary. most probably it should be rewritten - but
C     now I left it as is (just simplify it a bit)
C     if stored via PACK1 - should be get via DPACK1 - same for PACK2/DPACK2
C
C***********************************************************************
C
C
C    AAAAA  !!     POSSIBLE INCREASE OF YEL AND NOY UP TO 900, BUT THEN
C   A    A  !!     DURING I/O THERE WILL BE EXCHANGE IN 2 RECORDS
C  A     A  !!     !!!-RECORD LENGTH <= TRACK LENGTH = 7040 BYTES
C  AAAAAAA  !!
C  A     A         ...'NREC...=>...'2*NREC-1...
C  A     A  !!
C
C
C        OPERATIONS VY = CABS... AND VJR = CABS... CAN BE EXCLUDED
C     BY DECLARING: EQUIVALENCE((Y(1,1),RY(1,1,1)),(VJ(1),RJR(1,1))
C                   REAL RY(2,100,100), RJR(2,100)
C
C    AND IN CODE: IF(RJR(1,I).EQ.0.0 .AND. RJR(2,I).EQ.0.0) GO TO 10
C                 IF(RY(1,II,I).EQ.0.0 .AND. RY(2,II,I).EQ.0.0) GO TO 20
C
c$LARGE: VJ,Y
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Epsiw
      INTEGER i , ie , if , ii , jj3 , k , Kitu , kl , Kpky , kym ,
     &        Limerr , Nrec

      INTEGER Isize_maxnode
      DOUBLE COMPLEX Y(Isize_maxnode,Isize_maxnode)
      DOUBLE COMPLEX Vj(Isize_maxnode)

      COMMON /elpack/ Kpky(2,20)
      COMMON /serv  / Epsiw , Limerr , Kitu
      COMMON /kolnal/ Kol , Nal
      DOUBLE COMPLEX yel(762)
      INTEGER Kol(4) , entry , noy(2,762) , kerr/0/
c$LARGE: RY
      LOGICAL Nal(4)
      DATA kym/762/

C      OPEN ( 1, FILE=  'TMP1_XX',
C    *          STATUS='NEW',
C     *          ACCESS='DIRECT',
C     *          RECL=   12992,
C    *          CHARE= 'COMPAT',
C     *          MODE=  'READWRITE'
C     *    )
C
C     DEFINE FILE 1(20,3248,U,ASS1)
C
C     OPEN ( 2, FILE=  'TMP2_XX',
C    *          STATUS='NEW',
C     *          ACCESS='DIRECT',
C     *          RECL=   12992,
C    *          CHARE= 'COMPAT',
C     *          MODE=  'READWRITE'
C     *   )
C
C     DEFINE FILE 2(20,3248,U,ASS2)

      entry = 1
      if = Kol(1) + 1
      ie = Kol(1) + Kol(2) + Kol(3)
      GOTO 50
C      write (6,*) 'Y matrix at ENTRY =1 '

C*********************************************************************
      ENTRY pack2(Nrec,Y,Vj,Isize_maxnode)
C
C     INPUT FOR PACKING AND WRITING MATRIX Y, VECTOR J WITHOUT PACKING
      entry = 2
      if = Kol(1) + Kol(2) + 1
      ie = Kol(1) + Kol(2) + Kol(3)

C      write (6,*) 'Y matrix (ENTRY 2)'
C      do ii=1,IE
C            write (6,1220) (Y(ii,jj), jj=1,IE)
C      enddo
C1220   format (2x,'(',1x,e12.5,1x,e12.5,')')


 50   k = 1
      DO i = 1 , ie
         jj3 = i
         IF ( jj3.GE.if ) jj3 = ie
         DO ii = 1 , jj3
C     VY=CABS(Y(II,I))
C     IF(VY     .EQ.0.0 ) GO TO 20
            IF ( dble(Y(ii,i)).NE.0.0D0 .OR. dimag(Y(ii,i)).NE.0.0D0 )
     &           THEN
               IF ( k.GT.kym ) GOTO 110
               yel(k) = Y(ii,i)
               noy(1,k) = ii
               noy(2,k) = i
               k = k + 1
            ENDIF
         ENDDO
      ENDDO

      Kpky(entry,Nrec) = k - 1
C      print *,'WR: ENTRY, NREC, KPKY =',ENTRY, NREC, KPKY

C     ACTUAL WRITE OPERATION

C      print *,'YEL, NOY, NOY: WRITING'
C      DO I=1,50
C        print *,I, YEL(I), NOY(1,I), NOY(2,I)
C      ENDDO


      IF ( entry.EQ.1 ) CALL savedata(1,Nrec,yel,noy,Vj,Isize_maxnode,
     &                                kym)

C      WRITE(1,REC=NREC) YEL, NOY, VJ
C     IF(ENTRY.EQ.1) WRITE(1'NREC)YEL,NOY,VJ  !! REC=20 !!

      IF ( entry.EQ.2 ) CALL savedata(2,Nrec,yel,noy,Vj,Isize_maxnode,
     &                                kym)
C      WRITE(2,REC=NREC) YEL, NOY, VJ
C     IF(ENTRY.EQ.2) WRITE(2'NREC)YEL,NOY,VJ  !! REC=20 !!

      RETURN

C********************************************************************
      ENTRY dpack1(Nrec,Y,Vj,Isize_maxnode)
C
C     ENTRY POINT FOR READING AND UNPACKING VJ & Y
C
      entry = 1
C     READ(1,REC=NREC,ERR=120) YEL,NOY,VJ
C 121 READ(1'NREC,ERR=120)YEL,NOY,VJ
 121  CALL getdata(1,Nrec,yel,noy,Vj,Isize_maxnode,kym)

      GOTO 70

C**********************************************************************
      ENTRY dpack2(Nrec,Y,Vj,Isize_maxnode)
C
C     READING Y AND UNPACKING
C
      entry = 2
C      print *, 'DPACK2: NREC: ', NREC
C     READ(2,REC=NREC,ERR=130) YEL,NOY,VJ
C 131 READ(2'NREC,ERR=130)YEL,NOY,VJ
 131  CALL getdata(2,Nrec,yel,noy,Vj,Isize_maxnode,kym)

C      print *,'YEL, NOY, NOY: READING'
C      DO I=1,50
C        print *,I, YEL(I), NOY(1,I), NOY(2,I)
C      ENDDO


 70   kl = Kpky(entry,Nrec)
C      print *,'READ: ENTRY, NREC, KPKY =',ENTRY, NREC, KPKY

      DO k = 1 , kl
         Y(noy(1,k),noy(2,k)) = yel(k)
      ENDDO

C      write (6,*) 'Y matrix (UNPACKED)'
C      do ii=1,4
C            write (6,1220) (Y(ii,jj), jj=1,4)
C      enddo

      RETURN

C
 110  WRITE (6,111) Nrec , kym
 111  FORMAT (10X,'MATRIX  Y AT FREQUENCY ',I4,'CONTAINS'/10X,
     &        'MORE THAN',I4,'ELEMENTS ')
      STOP

C
      kerr = kerr + 1
      IF ( kerr.LT.Limerr ) GOTO 121
      WRITE (6,122) kerr , entry
      STOP

C
      kerr = kerr + 1
      IF ( kerr.LT.Limerr ) GOTO 131
      WRITE (6,122) kerr , entry
      STOP
 122  FORMAT (10X,I2,'-TH FILE IO ERROR',I2)
C     DEBUG SUBTRACE,INIT(NREC,YEL,NOY,VJ,IF,IE,ENTRY,K)
      END



      SUBROUTINE savedata(Ir,Nrec,Yel,Noy,Vj,Size_vj,Size_yel)
C TO SAVE DATA IN FILE. NOW - EMULATE BY STATIC ARRAYS TO SIMPLIFY
c
C IR - FIRST OR SECOND SET OF MATRICES AFTER FISRT OF SECOND REDUCTION
C NREC- NUMBER OF FREQUENCY - SO FAR 1-20
C YEL - ELEMENTS OF Y -MATRIX
C VJ - VECTOR OF FREE SOURCES
C NOI - POSITIONS OF ELEMENTS IN Y-MATRIX
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i
      INTEGER Ir , Nrec , Size_vj , Size_yel
      DOUBLE COMPLEX Vj(Size_vj)
      DOUBLE COMPLEX Yel(Size_yel)
      INTEGER Noy(2,Size_yel)


      DOUBLE COMPLEX Yel_save(2,20,1000)
      DOUBLE COMPLEX Vj_save(2,20,200)
      INTEGER Noy_save(2,20,2000,2)
      COMMON /savedata_cb/ Yel_save , Vj_save , Noy_save

      DO i = 1 , Size_vj
         Vj_save(Ir,Nrec,i) = Vj(i)
      ENDDO

      DO i = 1 , Size_yel
         Yel_save(Ir,Nrec,i) = Yel(i)
         Noy_save(Ir,Nrec,i,1) = Noy(1,i)
         Noy_save(Ir,Nrec,i,2) = Noy(2,i)
      ENDDO

C      print *,'SAVE::: YEL, NOY, NOY'
C      DO I=1,50
C       print *,I, YEL_SAVE(IR,NREC,I), NOY_SAVE(IR,NREC,I,1),
C     +            NOY_SAVE(IR,NREC,I,2)
c      ENDDO



      END

      SUBROUTINE getdata(Ir,Nrec,Yel,Noy,Vj,Size_vj,Size_yel)
C TO GET DATA IN FROM FILE. NOW - EMULATE BY STATIC ARRAYS TO SIMPLIFY
c
C IR - FIRST OR SECOND SET OF MATRICES AFTER FISRT OF SECOND REDUCTION
C NREC- NUMBER OF FREQUENCY - SO FAR 1-20
C YEL - ELEMENTS OF Y -MATRIX
C VJ - VECTOR OF FREE SOURCES
C NOI - POSITIONS OF ELEMENTS IN Y-MATRIX
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i
      INTEGER Ir , Nrec , Size_vj , Size_yel
      DOUBLE COMPLEX Vj(Size_vj)
      DOUBLE COMPLEX Yel(Size_yel)
      INTEGER Noy(2,Size_yel)


      DOUBLE COMPLEX Yel_save(2,20,1000)
      DOUBLE COMPLEX Vj_save(2,20,200)
      INTEGER Noy_save(2,20,2000,2)
      COMMON /savedata_cb/ Yel_save , Vj_save , Noy_save

      DO i = 1 , Size_vj
         Vj(i) = Vj_save(Ir,Nrec,i)
      ENDDO

      DO i = 1 , Size_yel
         Yel(i) = Yel_save(Ir,Nrec,i)
         Noy(1,i) = Noy_save(Ir,Nrec,i,1)
         Noy(2,i) = Noy_save(Ir,Nrec,i,2)
      ENDDO

      END
