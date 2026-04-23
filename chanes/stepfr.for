!*==STEPFR.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE stepfr(Om,Nrec,Y,Vj,Isize_maxnode)
C*********************************************************************
C* Subroutine for forming and reduction of multi-channel linear      *
C* subsystem and organizing their I/O on the DA.                     *
C*                                                                   *
C*********************************************************************

C Change from 30.01.91. See MAIN
C$LARGE: BUFFER
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , ifind1 , ifind2 , ilocal , ind1 , ind2 , j , jj ,
     &        k12 , k123 , k3 , ki , kj , kkk , Kn , Kn1 , Knc , Knr ,
     &        Knr1 , m
      INTEGER mf , n , nend , nf , Nrec , nu
      COMMON /maty  / Buffer(6000) , Buflen
      DOUBLE COMPLEX Buffer
      INTEGER*4 Buflen
      COMMON /blk2  / Knc , Knr , Kn , Knr1 , Kn1

      INTEGER Isize_maxnode
      DOUBLE COMPLEX Vj(Isize_maxnode)
      DOUBLE COMPLEX Y(Isize_maxnode,Isize_maxnode)

      COMMON /kolnal/ Kol , Nal
      INTEGER Kol(4) , dy , flag
      LOGICAL Nal(4)
      DOUBLE PRECISION Om

C Change from 30.01.91. See MAIN
C Built-in functions for indexing:
      ifind1(i,m,nu) = (i+(m-1)*nu)
      ifind2(i,j,m,nu,mf) = nu*mf + i + (j-1)*nu + (m-1)*nu*nu


      dy = Isize_maxnode
C @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      kkk = Kol(1) + Kol(2) + Kol(3)

C Initialization of the matrices Y and the vector VJ
      CALL ziny(Y,Vj,Isize_maxnode)
C      print *,'after ZINY'
C      print *,(VJ(ikkk),ikkk=1,20)

      IF ( Nal(1) ) THEN

C     LIN.CONST. - FORMATION
         CALL lnform(Om,Y,Vj,Isize_maxnode)

C   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C      WRITE(6,2)((III,JJJ,Y(III,JJJ),JJJ=1,KKK),III=1,KKK)
C   1  FORMAT(2X,' Y   STEPFR ',I2,' KOL(1),KOL(2),KOL(3),KKK=',4I3)
C   2  FORMAT(2X,'Y(',I3,',',I3,')=',1X,E13.6,1X,E13.6)
C      WRITE(6,3) (III, VJ(III),III=1,KKK)
C   3  FORMAT(2X,'STEPFR VJ(',I3,')=',E13.6,',',E13.6)

         IF ( Kol(1).NE.0 ) THEN
            nf = 1
            nend = Kol(1)
            n = Kol(1) + Kol(2) + Kol(3)
C  Reduction Y
            CALL luslv(Y,dy,n,nf,nend,flag)

C  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C  Reduction VJ
            CALL lufrw(Y,Vj,dy,n,nf,nend,flag)
         ENDIF

C  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C      write (6,*) 'Y matrix'
C      do ii=1,n
C            write (6,120) (Y(ii,jj), jj=1,n)
C      enddo
C120   format (2x,'(',1x,e12.5,1x,e12.5,')')


C Packing and writing of matrix Y and vector VJ to DA
         CALL pack1(Nrec,Y,Vj,Isize_maxnode)
      ENDIF
      GOTO 20
C**********************************************************************
      ENTRY double(Om,Nrec)
C**********************************************************************
C*                                                                    *
C* Entry in subroutine for reformation (when modifying elements)      *
C*                                                                    *
C**********************************************************************

C Initialization of Y
      CALL ziny(Y,Vj,Isize_maxnode)
C Reading from DA and unpacking the formatted Y and VJ from constant elements
      CALL dpack1(Nrec,Y,Vj,Isize_maxnode)



 20   IF ( Nal(2) ) THEN

C Linear Variable - Initialization
         CALL vlform(Om,Y,Vj,Isize_maxnode)

C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         IF ( Kol(2).NE.0 ) THEN
            nf = Kol(1) + 1
            nend = Kol(1) + Kol(2)
            n = Kol(1) + Kol(2) + Kol(3)
C Reduction of the full matrix LIN.-subsystem
            CALL luslv(Y,dy,n,nf,nend,flag)

C @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C Reduction of the full vector of the given currents in the LIN.-subsystem
            CALL lufrw(Y,Vj,dy,n,nf,nend,flag)
         ENDIF
      ENDIF


C @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C Transformation of matrix Y to canonical form
      n = Kol(1) + Kol(2) + Kol(3)
      nf = Kol(1) + Kol(2) + 1
      CALL lucan(Y,dy,n,nf,flag)

C @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

C Execution of the operation for setting the input values
C (IF(OMEGA.NE.0.)). Linked with the DPF.
C
C
      k123 = Kol(1) + Kol(2) + Kol(3)
      IF ( Nrec.NE.1 ) THEN
         DO jj = 1 , k123
            Vj(jj) = Vj(jj)/2.D0
         ENDDO
      ENDIF

C Packing and writing Y to MD
      CALL pack2(Nrec,Y,Vj,Isize_maxnode)

C Transfer of result to the function and reductions in the /MATY/ block,
C storing Y and J. Reduced linear subckt of the entire frequency range.
      k3 = Kol(3)
      k12 = Kol(1) + Kol(2)
      DO ki = 1 , k3
         ind1 = ki + k12

C Change from 30.01.91, modified by Serdyuk G.V.
C      JR(KI,NREC)=VJ(IND1)       **** OLD VERSION
         Buffer(ifind1(ki,Nrec,k3)) = Vj(ind1)

         DO kj = 1 , k3
            ind2 = k12 + kj

C Change from 30.01.91, modified by Serdyuk G.V.
C      YR(KJ,KI,NREC)=Y(IND2,IND1)   **** OLD VERSION
            ilocal = ifind2(kj,ki,Nrec,k3,Kn)
            IF ( ilocal.GT.Buflen ) GOTO 100
            Buffer(ilocal) = Y(ind2,ind1)

C  NO    WRITE(6,3000) KJ,KI,NREC,YR(KJ,KI,NREC)
C  NO 3000 FORMAT(2X,'STEPFR,YR(',I3,',',I3,',',I3,')=',E12.5,2X,E12.5)
C      WRITE(6,3003) IND2,IND1,Y(IND2,IND1)
C 3003 FORMAT(2X,'        Y(',I3,',',I3,')=     ',E12.5,2X,E12.5)
         ENDDO
      ENDDO
      RETURN
 100  WRITE (6,2000)
C     DEBUG SUBTRACE,INIT(K12,K3,KI,KJ,IND1,IND2)

 2000 FORMAT (2X,' BUFFER SIZE FOR REDUCED MATRICES EXCEDED'/2X,
     &        ' ABNORMAL TERMINATION.')
      STOP



      END
