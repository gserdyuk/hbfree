!*==SORTUZ.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c





      SUBROUTINE sortuz

C  VARIANT WITHOUT FIXED. BOUNDARY NODES IN SYSTEM. UPPER.
C              STICKS WITH CURRENT TEXTS
C                 IN MODULES.
C
C

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i00 , i00end , i01 , i01end , i02 , i10 , i11 , i12 ,
     &        i20 , i20end , i40 , ide , idel , ie10 , ie11 , iin ,
     &        inuz , ip , ir , iterm2
      INTEGER itime , itime3 , ityp , itype , iva , iznn , j , kll ,
     &        Kprgrf , Kprlen , Kprlin , Kprnkr , Kprqup , Kprsol ,
     &        Kprsrt , Kprvar , luk0 , Mt , naaa , name
      INTEGER natadr , nb , nf , nfirst , nlen , np , nu
      INCLUDE 'circuit.i'
      LOGICAL Nal(4)
      INTEGER*4 nldim(100) , lcdim(100) , lvdim(100) , Kol(4)
      INTEGER*4 Nouz(100) , Inouz1(100)
      INTEGER*4 Kmod , Ksch , ivar , luk
      COMMON /kolnal/ Kol , Nal
      COMMON /koluz / Kmod , Ksch
      COMMON /mdla  / Mt(15)

      COMMON /nouzl / Nouz , Inouz1
      COMMON /print / Kprlen , Kprsrt , Kprnkr , Kprlin , Kprsol ,
     &                Kprvar , Kprgrf , Kprqup

      DIMENSION name(4)


C  BEGINNING. BOUNDARY CONDITIONS *****************************************

      Kol(1) = 0
      Kol(2) = 0
      Kol(3) = 0
      Kol(4) = 0
      Nal(1) = .FALSE.
      Nal(2) = .FALSE.
      Nal(3) = .FALSE.
      Nal(4) = .FALSE.
      luk0 = 100
      luk = luk0
      DO j = 1 , luk0
         nldim(j) = 0
         lvdim(j) = 0
         lcdim(j) = 0
         Nouz(j) = 0
         Inouz1(j) = 0
      ENDDO

C  DO FOR ALL TYPES OF ELEMENTS          *****     *****     *****

      ie10 = Nmpnt*20
      DO i10 = 1 , ie10 , 20
         IF ( Mpoint(i10+4).NE.0 ) THEN
            nfirst = Mpoint(i10+9)
            nlen = Mpoint(i10+6)
            np = nfirst

            name(1) = Mpoint(i10)
            name(2) = Mpoint(i10+1)
            name(3) = Mpoint(i10+2)
            name(4) = Mpoint(i10+3)
            ityp = Mpoint(i10+5)
C  DO FOR ELEMENTS WITHIN THE TYPE RANGE     *****     *****     *****
            ie11 = Mpoint(i10+4)
            DO i11 = 1 , ie11

C  ADDRESS OF ELEMENT IN NODEEL          *****     *****     *****
               natadr = nfirst + (i11-1)*nlen

C  IF THE ELEMENT IS LINEAR:             *****     *****     *****
               IF ( ityp.EQ.3 ) THEN

C  NON-LINEAR MULTIPOLE         *****     *****     *****
                  ide = Mpoint(i10+7)
                  ivar = Nodeel(natadr+ide+3+1)
                  CALL libmd1(name,ivar)
                  nb = 2 - 1
                  nf = Mpoint(i10+7) + 1 - 1

C  IF LINEAR MULTIPOLE                   *****     *****     *****
               ELSEIF ( ityp.EQ.2 ) THEN

C  LINEAR MULTIPOLE:          *****     *****     *****

                  ide = Nodeel(natadr+1)
                  ivar = Nodeel(natadr+ide+6)
                  itype = ivar
                  nb = 5 - 1
                  nf = Nodeel(natadr+1) + 4 - 1

C  IF Y-MATRIX                          *****     *****     *****
               ELSEIF ( ityp.EQ.5 ) THEN

C  Y-MATRIX                       *****     *****     ******

                  ide = Nodeel(np+1)
                  ivar = Nodeel(np+ide+6)
                  itype = ivar
                  nb = 5 - 1
                  nf = Nodeel(np+1) + 4 - 1
                  natadr = np
                  np = np + 8 + ide
               ELSE

C  IF LINEAR DOUBLE MULTIPOLE            *****     *****     *****
                  ivar = Nodeel(natadr+3)
                  itype = ivar
                  nb = 2 - 1
                  nf = 3 - 1
               ENDIF

C  ADDRESS OF FIRST AND LAST NODES BASED ON ITYP
C  PREVIOUS CHECK

C  LOOP THROUGH NODES IN ELEMENT:              *****     *****     *****
               DO i12 = nb , nf , 1
C  NODE HOMEP                            *****     *****     *****
                  inuz = Nodeel(natadr+i12)
                  IF ( inuz.NE.0 ) THEN
                     IF ( ityp.EQ.3 ) THEN

C  DETERMINE NODE TYPE                   *****     *****     *****

                        itype = Mt(i12)

C  CHECK IF NODES ARE WITHIN THE NUMBERED RANGE          *****     *****     *****

                        itime3 = nf - Mpoint(i10+8)
                        IF ( i12.GT.itime3 ) THEN

C  PROCESSING - NUMBERED RANGE           *****     *****     *****

                           luk = luk + 1
                           inuz = luk
                           Nodeel(natadr+i12) = luk
                        ENDIF
                     ENDIF

C  DETERMINE NODE TYPE AND CHECK IF THIS NODE IS ALREADY
C  FILLED IN ONE OF THE FILLED ARRAYS.
C  IF NOT - CONTINUE
                     IF ( itype.EQ.0 ) CALL find(Kol(1),lcdim,inuz)
                     IF ( itype.NE.0 ) THEN
                        IF ( itype.EQ.1 ) CALL find(Kol(2),lvdim,inuz)
                        IF ( itype.NE.1 ) THEN
                           IF ( itype.EQ.2 )
     &                          CALL find(Kol(3),nldim,inuz)
                           IF ( itype.EQ.2 ) THEN
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
C
               ENDDO
            ENDDO
         ENDIF
      ENDDO
C  DESTROY (SET TO 0) NODES IN LCDIM AND LVDIM IF THEY EXIST IN NLDIM
C  AND IN LCDIM IF THEY EXIST IN LVDIM

      CALL down(Kol(3),nldim,Kol(2),lvdim)
      CALL down(Kol(3),nldim,Kol(1),lcdim)
      CALL down(Kol(2),lvdim,Kol(1),lcdim)

C  REWRITE IN SEQUENCE FROM LCDIM, LVDIM, NLDIM TO INOUZ STARTING WITH
C  MINIMUM ELEMENTS OF THE ARRAYS

C  KSCN = NUMBER OF NODES IN CIRCUIT
C
      Ksch = 0
      CALL ranger(Ksch,Nouz,lcdim,Kol(1),Nal(1))
      CALL ranger(Ksch,Nouz,lvdim,Kol(2),Nal(2))
      CALL ranger(Ksch,Nouz,nldim,Kol(3),Nal(3))
      Kmod = luk - luk0
      Ksch = Ksch - Kmod
C  SHIFT NODES TO LUK0-KSCH POSITION
      idel = luk0 - Ksch
      DO i40 = 1 , luk0
         iznn = Nouz(i40)
         IF ( iznn.GT.luk0 ) Nouz(i40) = iznn - idel
      ENDDO
C CONVERSION OF NUMERATIONS
      DO nu = 1 , luk0
         i20end = Ksch + Kmod
         DO i20 = 1 , i20end
            iterm2 = Nouz(i20)
            IF ( iterm2.EQ.nu ) GOTO 80
         ENDDO
         GOTO 70
 80      Inouz1(nu) = i20
 70   ENDDO
C NOUZ(I)=NUMBER OF USE (BH.H.1)
C INOUZ1(I)=BH.H.1(USER TYPE)
C
C KMOD=NUMBER OF NODES IN MODULE.
C KSCH=NUMBER OF NODES IN CIRCUIT
C
C TRANSFER TOPOLOGICAL INFORMATION
C LOOP OVER ELEMENT TYPES
      i00end = Nmpnt*20
      DO i00 = 1 , i00end , 20
         nfirst = Mpoint(i00+9)
         nlen = Mpoint(i00+6)
         ityp = Mpoint(i00+5)
         i01end = Mpoint(i00+4)
         np = nfirst

         DO i01 = 1 , i01end
            natadr = nfirst + (i01-1)*nlen

            IF ( ityp.EQ.2 ) THEN
               nb = 4
               nf = Nodeel(natadr+1) + 3
            ELSEIF ( ityp.EQ.3 ) THEN

               nb = 1
               nf = Mpoint(i00+7)
            ELSEIF ( ityp.EQ.4 ) THEN
            ELSEIF ( ityp.EQ.5 ) THEN
               nb = 4
               nf = Nodeel(np+1) + 3
               natadr = np
               np = np + 8 + Nodeel(np+1)
            ELSE
               nb = 1
               nf = 2
            ENDIF
            DO i02 = nb , nf
               itime = Nodeel(natadr+i02)
               IF ( itime.NE.0 ) THEN
                  IF ( itime.GT.luk0 ) itime = itime - idel
                  Nodeel(natadr+i02) = Inouz1(itime)
                  naaa = natadr + i02
               ENDIF
C     WRITE(6,900) NATADR,I02,NODEEL(NAAA)
C 900 FORMAT(2X,'SORTUZ : NODEEL(',I3,'+',I3,')=',I5)
            ENDDO
         ENDDO
      ENDDO
      kll = Kol(1) + Kol(2) + Kol(3)
      IF ( Kprsrt.GE.2 ) THEN
         WRITE (6,800)
 800     FORMAT (10X,'RESULTS OF WORK: SORTUZ:')
         WRITE (6,810) (iva,Nouz(iva),iva=1,kll)
 810     FORMAT (2X,'NOUZ(',I3,')=',I4)
         WRITE (6,820)
 820     FORMAT (120(' '))
         WRITE (6,830) (iva,Inouz1(iva),iva=1,kll)
 830     FORMAT (2X,'INOUZ1(',I3,')=',I4)
         WRITE (6,833) (Kol(ip),ip=1,3)
 833     FORMAT (10X,'ARRAY KOL:'/6X,3(4X,I4))
         WRITE (6,836) (Nal(ir),ir=1,3)
 836     FORMAT (10X,'ARRAY NAL:'/6X,3(4X,L4))
         WRITE (6,840)
 840     FORMAT (10X,'ARRAY NODEEL AFTER RECODING')
         WRITE (6,850) (iin,Nodeel(iin),Nodeel(iin),iin=1,luk0)
 850     FORMAT (2X,I3,2X,I4,2X,A4)
      ENDIF

      END


      SUBROUTINE find(Nend,Nndim,Num)
C     IN THE ARRAY NNDIM, SEARCHES FOR THE VALUE NUM. IF IT IS NOT FOUND, WRITES IT.
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i
      INTEGER*4 Nend , Num , Nndim(100)
      IF ( Nend.NE.0 ) THEN
         DO i = 1 , Nend
            IF ( Nndim(i).EQ.Num ) RETURN
         ENDDO
      ENDIF
      Nend = Nend + 1
      Nndim(Nend) = Num
      END


      SUBROUTINE down(Ed,Domin,Ei,Insert)
C     DESTROYING ELEMENTS OF THE INSERT WHICH EXIST IN THE DOMAIN
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER id , ii
      INTEGER*4 Domin(100) , Insert(100) , Ed , Ei
      DO id = 1 , Ed
         IF ( Domin(id).NE.0 ) THEN
            DO ii = 1 , Ei
               IF ( Domin(id).EQ.Insert(ii) ) THEN
                  Insert(ii) = 0
                  GOTO 10
               ENDIF
            ENDDO
         ENDIF
 10   ENDDO
      END


      SUBROUTINE ranger(Ks,Inz,Ndim,End,Exist)
C     SELECTS ELEMENTS FROM NDIM TO INZ STARTING FROM MIN
C
C     END - IS REDEFINED. AT THE ENTRY, END - THE TOTAL NUMBER OF NODES
C  OF THE GIVEN TYPE BEFORE EXCLUSION BY DOWN, AT EXIT - THE ACTIVE
C  NUMBER OF NODES OF THE GIVEN TYPE /SEE. OUTPUT/
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , k , min , nmin
      LOGICAL Exist
      INTEGER*4 Inz(100) , Ndim(100) , Ks , End
      Exist = .FALSE.
      IF ( End.EQ.0 ) RETURN
      k = 0
      DO
         min = 9999
         DO i = 1 , End
            IF ( Ndim(i).NE.0 .AND. Ndim(i).LT.min ) THEN
               nmin = i
               min = Ndim(nmin)
            ENDIF
         ENDDO
         IF ( min.EQ.9999 ) THEN
            End = k
            Exist = .TRUE.
            RETURN
         ELSE
            k = k + 1
            Ks = Ks + 1
            Ndim(nmin) = 0
            Inz(Ks) = min
         ENDIF
      ENDDO
C     DEBUG INIT(KS,INZ,EXIST)
      END
