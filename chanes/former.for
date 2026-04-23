!*==LNFORM.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE lnform(Omega1,Y,Vj,Isize_maxnode)
C
C

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i10 , i10e , i11 , i11e , ibase , ide , iii , ityp ,
     &        lenbas , lenth , natadr , newnp , nfirst , nlen , np ,
     &        np10
      INCLUDE 'circuit.i'
      INCLUDE 'charint.i'

      INTEGER Isize_maxnode
      DOUBLE COMPLEX Y(Isize_maxnode,Isize_maxnode)
      DOUBLE COMPLEX Vj(Isize_maxnode)
C      COMMON/FORMY/VEKTJ(70),Y(70,70)

      LOGICAL dsign
      INTEGER varabl , adpar1 , adres1
      DOUBLE PRECISION omega , Omega1
      CHARACTER*4 name(4)


C      PRINT *, 'ENTER LNFORM'
C      print *,(vj(ikkk),ikkk=1,20)

      varabl = 0
      GOTO 5
C  ENTRY FOR THE FORMATION OF A NEW MODIFIED PART.
      ENTRY vlform(Omega1,Y,Vj,Isize_maxnode)
      varabl = 1
C  DEFINING THE FREQUENCY FLAG.
 5    dsign = .TRUE.
      IF ( Omega1.LT.0.0D0 ) dsign = .FALSE.
      omega = dabs(Omega1)

C     DO FOR TYPES OF ELEMENTS
      i10e = Nmpnt*20
      DO i10 = 1 , i10e , 20
C
C     IF NO ELEMENT OF THE GIVEN TYPE EXISTS - SKIP
         IF ( Mpoint(i10+4).NE.0 ) THEN
C
C     IT MAY BE THAT ONE OF THE ELEMENTS OF THIS TYPE DOES NOT CHANGE,
C     AND THE FORMATION - IS CHANGED. SKIP
            IF ( varabl.EQ.0 .OR. Mpoint(i10+14).NE.0 ) THEN

               nfirst = Mpoint(i10+9)
               nlen = Mpoint(i10+6)
               ityp = Mpoint(i10+5)
               Lenpar = Mpoint(i10+10)
               ibase = Mpoint(i10+13)
C     LENBAS=MPOINT(I10+12)-LENPAR-LENNTP
               lenbas = Mpoint(i10+12)
               adpar1 = Mpoint(i10+11)
               name(1) = i2c(Mpoint(i10))
C      WRITE(6,1110) NAME(1)
C 1110 FORMAT(2X,'LNFORM: NAME(1)=',A4)
C      WRITE(6,1111) NFIRST,NLEN,ITYP,LENPAR,IBASE,
C     *              LENNTP,LENBAS,ADPAR1
C 1111 FORMAT(2X,'LNFORM: NFIRST, NLEN, ITYP, LENPAR, IBASE, LENNTP=',
C     *  6I5,/10X,'LENBAS=',I3,' ADPAR1=',I3)
               name(2) = i2c(Mpoint(i10+1))
               name(3) = i2c(Mpoint(i10+2))
               name(4) = i2c(Mpoint(i10+3))
C
               np = nfirst
C
C     DO FOR ELEMENTS OF THIS TYPE
               i11e = Mpoint(i10+4)
               DO i11 = 1 , i11e
                  natadr = nfirst + (i11-1)*nlen


C ## LINEAR 2X-MULTIPLIER ############################################
                  IF ( ityp.NE.1 ) THEN

C ## LINEAR MULTIPLIER ###########################################
                     IF ( ityp.EQ.2 ) THEN
C  CHECK THE CORRESPONDENCE OF IVAR AND VARABL
                        ide = Nodeel(natadr+1)
                        IF ( Nodeel(natadr+ide+6).EQ.varabl ) THEN
C  START ADDRESS AND LENGTH OF THE INDIVIDUAL PARAMETER SUBSET.
                           adres1 = Nodeel(natadr+ide+5)
                           lenth = Nodeel(natadr+2)
C  CALLING THE ELEMENT MATRIX.
                           CALL liblin(name,omega,Nnetpr,Lenntp,adpar1,
     &                                 Lenpar,adres1,lenth,ide)
C  ENTRY INTO THE TOTAL MATRIX.
                           CALL linnp(omega,natadr,dsign,Y,Vj,
     &                                Isize_maxnode)
                        ENDIF

C ## Y-MATRIX ########################################################
                     ELSEIF ( ityp.EQ.5 ) THEN
                        np10 = np + 10
C     PRINT 355,NP,(NODEEL(JNP),JNP=NP,NP10)
C 355 FORMAT(2X,'NODEEL(',I3,')=',A4,9(1X,I5))

C  CHECK THE MATCHING OF IVAL AND VARABL
                        ide = Nodeel(np+1)
                        IF ( Nodeel(np+ide+6).EQ.varabl ) THEN
C  START ADDRESS AND LENGTH OF THE INDIVIDUAL PARAMETER SUBSET.
                           adres1 = Nodeel(np+ide+5)
                           lenth = Nodeel(np+2)
C  CALLING THE ELEMENT MATRIX.
                           CALL liblin(name,omega,Nnetpr,Lenntp,adpar1,
     &                                 Lenpar,adres1,lenth,ide)
C  ENTRY INTO THE TOTAL MATRIX.
                           newnp = np
                           CALL linnp(omega,newnp,dsign,Y,Vj,
     &                                Isize_maxnode)
                           np = np + 8 + ide
                        ENDIF
C
C
C ### NONLINEAR MULTIPOLE ########################################
                     ELSEIF ( ityp.EQ.3 ) THEN
C  CHECK IVAR - VARABL
                        ide = Mpoint(i10+7)
                        IF ( Nodeel(natadr+ide+4).EQ.varabl ) THEN
C  SOME VALUES (INCLUDING ADDRESS)
                           iii = i10
                           adres1 = Nodeel(natadr+ide+3)
C  Y-MATRIX OF THE LINEAR PART OF THE NONLINEAR ELEMENT...
                           lenbas = Mpoint(i10+12)
                           CALL libmd2(name,omega,Nnetpr,Lenntp,adpar1,
     &                                 Lenpar,adres1,lenbas)
C  ...IS ENTERED INTO THE TOTAL Y-MATRIX.
                           CALL nonlin(iii,natadr,dsign,Y,Vj,
     &                                 Isize_maxnode)
                        ENDIF
                     ENDIF
C  CHECK THE CORRESPONDENCE OF THE VALIDATION FLAG IVAR TO THE TYPE OF
C  FORMATION VARABL
                  ELSEIF ( Nodeel(natadr+3).EQ.varabl ) THEN
                     CALL lin2p(omega,name,Param(1),Param(adpar1),
     &                          Param(ibase+i11-1),dsign,natadr,Y,Vj,
     &                          Isize_maxnode)
                  ENDIF

C END OF DO - BY ELEMENTS
C      print *,'element  ',I11
C      print *,(vj(ikkk),ikkk=1,20)


               ENDDO
            ENDIF
         ENDIF
C END OF DO - BY ELEMENT TYPES.
      ENDDO
      END


      SUBROUTINE lin2p(Om,N,Pn,Pt,Pe,S,Nr,Y,J,Isize_maxnode)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER iy , jy , kk , Nr
      INCLUDE 'circuit.i'

      DOUBLE PRECISION Om
      LOGICAL S
      DIMENSION names(4) , N(4)
      DOUBLE PRECISION Pn , Pt , Pe , zn
      DOUBLE COMPLEX ys
C      COMMON/NODEL/NODEEL(1)

C      COMMON/FORMY/J(70),Y(70,70)
      INTEGER Isize_maxnode
      DOUBLE COMPLEX Y(Isize_maxnode,Isize_maxnode)
      DOUBLE COMPLEX J(Isize_maxnode)

      CHARACTER*4 names/'R   ' , 'L   ' , 'C   ' , 'G   '/ , N

      IF ( N(1).EQ.names(1) ) ys = dcmplx(1.D0/Pe,0.0D0)
      IF ( N(1).EQ.names(3) ) ys = dcmplx(Pt,Om*Pe)
      IF ( N(1).EQ.names(2) ) THEN
         zn = Pt**2 + (Pe*Om)**2
         ys = dcmplx(Pt/zn,-Pe*Om/zn)
      ENDIF
      IF ( N(1).EQ.names(4) ) ys = dcmplx(Pe,0.0D0)
      IF ( .NOT.S ) ys = dconjg(ys)
C  CALCULATE THE NODE NUMBERS OF THE ELEMENT
      iy = Nodeel(Nr+1)
      jy = Nodeel(Nr+2)
C     WRITE(6, 13) NR, IY, JY
C  13 FORMAT(2X,'LIN2P: NR, IY, JY =',3I5)
C  IF THE ELEMENT IS DISCONNECTED - RETURN
      IF ( iy.EQ.jy ) RETURN
C  IF CONNECTED TO GROUND - PROCESS BELOW
      IF ( iy.EQ.0 .OR. jy.EQ.0 ) THEN
C  GROUNDED ELEMENT
         kk = iy + jy
         IF ( kk.NE.0 ) Y(kk,kk) = Y(kk,kk) + ys
         RETURN
      ENDIF
C  AND IF INCLUDED IN THE CIRCUIT AND NOT CONNECTED TO GROUND:
      Y(iy,jy) = Y(iy,jy) - ys
      Y(jy,iy) = Y(jy,iy) - ys
C     WRITE(6, 15) IY,JY,Y(IY,JY),JY,IY,Y(JY,IY)
C  15 FORMAT(2X,'LIN2P:'/(2X,'Y(',I3,',',I3,')=',E12.5,2X,E12.5))
      RETURN
C     WRITE(6, 15) KK,KK,Y(KK,KK)
      END


      SUBROUTINE linnp(Omega,Nr,S,Y,Vectj,Isize_maxnode)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i20 , i30 , i40 , i50 , i60 , im , in , Is , iy , j , jn ,
     &        jy , Nr , nri20
      INCLUDE 'circuit.i'
      LOGICAL S , input
      INTEGER end

      DOUBLE PRECISION Omega
      DOUBLE COMPLEX zero/(0.0D0,0.0D0)/ , yd
      DOUBLE COMPLEX Suby , Subj

c      COMMON/POINT/MPOINT(1)/NODEL/NODEEL(1)
      COMMON /subc  / Suby(15,15) , Subj(15)

C      COMMON/FORMY/VEKTJ(70),Y(70,70)
      INTEGER Isize_maxnode
      DOUBLE COMPLEX Y(Isize_maxnode,Isize_maxnode)
      DOUBLE COMPLEX Vectj(Isize_maxnode)

C ENTRY FOR INSERTING INTO THE Y-MATRIX OF NONLINEAR MULTIPOLES...
      input = .FALSE.

      end = Nodeel(Nr+1)
      Nr = Nr + 3
C      WRITE(6,100) END, NR
C  100 FORMAT(2X,'LINNP: END=',I4,', NR=',I4)
C      WRITE(6, 2)((II,JJ,SUBY(II,JJ),JJ=1,END),II=1,END)
C    2 FORMAT (2X,'LINNP:'/('  SUBY(',I3,',',I3,')=',E13.6,2X,E13.6))
      GOTO 3
C ...OR Y-MATRIX OF THE LINEAR PART OF NONLINEAR ELEMENTS.
      ENTRY nonlin(Is,Nr,S,Y,Vectj,Isize_maxnode)
      end = Mpoint(Is+7)

      input = .TRUE.
C CALCULATE THE COMPLETE COMPRESSION MATRIX (IF THE FREQUENCY < 0.)
 3    IF ( .NOT.(S) ) THEN
         DO in = 1 , end
            DO im = 1 , end
               Suby(im,in) = dconjg(Suby(im,in))
            ENDDO
         ENDDO
      ENDIF
C      WRITE(6, 2) ((II,JJ,SUBY(II,JJ),JJ=1,END),II=1,END)

C TRANSFORM THE SUBMATRIX TO THE NECESSARY FORMAT:
C  Y(I,I)=Y(I,I)+Y(I,J),PRE J=1..N,(J.NE.I).                !
      DO in = 1 , end
         yd = zero
         DO jn = 1 , end
            yd = yd + Suby(in,jn)
         ENDDO
         Suby(in,in) = yd
      ENDDO
C      WRITE(6, 2) ((II,JJ,SUBY(II,JJ),JJ=1,END),II=1,END)

C ENTRY INTO THE DIAGONAL OF THE GROUND CONNECTION NODES:
      DO i20 = 1 , end
         nri20 = Nr + i20
C      WRITE(6,110) END, NR, NRI20, NODEEL(NRI20), I20
C  110 FORMAT(2X,'END=',I4,' NR=',I4,' NODEEL(',I4,')=',I4,' I20=',I4)
         IF ( Nodeel(Nr+i20).EQ.0 ) THEN
            DO i30 = 1 , end
C      WRITE(6,120) I30, I20, SUBY(I30,I30), SUBY(I20,I30)
C  120 FORMAT(2X,'I30',I4,' I20',I4,' SUBY=',2(1X,2(1X,E13.6)))
               Suby(i30,i30) = Suby(i30,i30) - Suby(i20,i30)
            ENDDO
         ENDIF
C      WRITE(6,130) SUBY(I30,I30)
C  130 FORMAT(2X,'SUBY-SUBY=',2(1X,E13.6))
      ENDDO
C      WRITE(6, 2) ((II,JJ,SUBY(II,JJ),JJ=1,END),II=1,END)

C FORMATION (APPROPRIATELY)
      DO i40 = 1 , end
         iy = Nodeel(Nr+i40)
         IF ( iy.NE.0 ) THEN

            DO i50 = 1 , end
               jy = Nodeel(Nr+i50)
               IF ( jy.NE.0 ) THEN

C ATTENTION! MINED 28.02.92

                  IF ( i40.NE.i50 .AND. iy.EQ.jy ) THEN
                     WRITE (6,43)
 43                  FORMAT (2x,
     &                      ' ####### LINNP : LOOP IS DETECTED  #######'
     &                      )
                     GOTO 50
                  ENDIF

                  Y(jy,iy) = Y(jy,iy) + Suby(i50,i40)
               ENDIF
C      WRITE(6,45) JY,IY,Y(JY,IY)
C   45 FORMAT(2X,'LINNP:'/ 2X,'Y(',I3,',',I3,')=',2(E12.5,2X))

 50         ENDDO
         ENDIF
      ENDDO


      IF ( input ) RETURN
      IF ( Nodeel(Nr).NE.2 ) RETURN
C ENTRY INTO VECTORS OF CURRENT SOURCES
C OF LINEAR ACTIVE N-RECEIVERS.
      DO i60 = 1 , end
         j = Nodeel(Nr+i60)
         IF ( j.NE.0 ) Vectj(j) = Vectj(j) + Subj(i60)
C      WRITE(6, 55) J,VEKTJ(J)
C   55 FORMAT(2X,'LINNP:'/ 2X,'VEKTJ(',I3,')=',2(E12.5,2X))
      ENDDO

      END
