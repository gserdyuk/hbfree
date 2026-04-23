!*==LIBLIN.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE liblin(Name,Om,Na1,Le1,Na2,Le2,Na3,Le3,Npol)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER Le1 , Le2 , Le3 , Na1 , Na2 , Na3 , Npol
      DOUBLE PRECISION Om
      CHARACTER*4 names(20) , Name(4)
      DOUBLE PRECISION Param
      COMMON /params/ Param(1)
      DATA names/'P2  ' , 'E   ' , 'LIB0' , 'LL0 ' , 'J   ' , 'YTAB' ,
     &     'STAB' , 'LIB ' , 'LANG' , 'SHL ' , 'KZ  ' , 'XX  ' ,
     &     'SMPL' , 'MPL ' , 'DISC' , 'INDS' , 4*'XXXX'/
C
C
C
C     WRITE(6, 5) NAME,OM,NA1,LE1,NA2,LE2,NA3,LE3,NPOL
C   5 FORMAT(2X,'LIBLIN : N=',4A4,' OM=',E12.5,' NA1=',I5,' LE1=',I5/
C    *    11X,'NA2=',I5,' LE2=',I5,' NA3=',I5,' LE3=',I5,' NPOL=',I5)
C
C
      IF ( Name(1).EQ.names(3) ) THEN
         IF ( Name(2).EQ.names(4) )
     &        CALL line1(Om,Param(Na1),Le1,Param(Na2),Le2,Param(Na3),
     &        Le3)
      ENDIF
      IF ( Name(1).EQ.names(1) )
     &     CALL rffipn(Om,Param(Na1),Le1,Param(Na2),Le2,Param(Na3),Le3)
      IF ( Name(1).EQ.names(2) ) CALL emf(Om,Param(Na1),Le1,Param(Na2),
     &     Le2,Param(Na3),Le3)
      IF ( Name(1).EQ.names(5) )
     &     CALL jdrive(Om,Param(Na1),Le1,Param(Na2),Le2,Param(Na3),Le3)
      IF ( Name(1).EQ.names(6) ) CALL ytab(Om,Param(Na1),Le1,Param(Na2),
     &     Le2,Param(Na3),Le3,Npol)
      IF ( Name(1).EQ.names(7) ) CALL stab(Om,Param(Na1),Le1,Param(Na2),
     &     Le2,Param(Na3),Le3,Npol)
      IF ( Name(1).EQ.names(8) .AND. Name(2).EQ.names(9) )
     &     CALL lang(Om,Param(Na1),Le1,Param(Na2),Le2,Param(Na3),Le3,
     &     Npol)
      IF ( Name(1).EQ.names(8) .AND. Name(2).EQ.names(13) )
     &     CALL smpl(Om,Param(Na1),Le1,Param(Na2),Le2,Param(Na3),Le3,
     &     Npol)
      IF ( Name(1).EQ.names(8) .AND. Name(2).EQ.names(14) )
     &     CALL mp(Om,Param(Na1),Le1,Param(Na2),Le2,Param(Na3),Le3,Npol)
C      IF(NAME(1).EQ.NAMES(10).AND.NAME(2).EQ.NAMES(11)) CALL SHLEIF(OM,
C     +                                    PARAM(NA1),LE1,PARAM(NA2),
C     +                                    LE2,PARAM(NA3),LE3,NPOL,
C     +                                    NAME(2))
C      IF(NAME(1).EQ.NAMES(10).AND.NAME(2).EQ.NAMES(12)) CALL SHLEIF(OM,
C     +                                    PARAM(NA1),LE1,PARAM(NA2),
C     +                                    LE2,PARAM(NA3),LE3,NPOL,
C     +                                    NAME(2))
C      IF(NAME(1).EQ.NAMES(15))CALL DISCONT(OM,PARAM(NA1),LE1,PARAM(NA2),
C     +                                    LE2,PARAM(NA3),LE3,NPOL)

      IF ( Name(1).EQ.names(16) )
     &     CALL indsv(Om,Param(Na1),Le1,Param(Na2),Le2,Param(Na3),Le3)
C
C     DEBUG SUBTRACE
      END




      SUBROUTINE rffipn(Om,P1,L1,P2,L2,P3,L3)
C     CONVERT THE VALUE TO EPSI / REDUCE TO 0.1E-05 /
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Epsi
      INTEGER i , ika , Kitu , L1 , L2 , L3 , Limerr , numer
      COMMON /serv  / Epsi , Limerr , Kitu
      DOUBLE PRECISION P1 , P2 , P3 , sq1 , sq2
      DOUBLE PRECISION Om , pi
      DOUBLE COMPLEX Suby , Subj

      COMMON /subc  / Suby(15,15) , Subj(15)
      DIMENSION P1(L1) , P2(L2) , P3(L3)

      pi = 4.D0*datan(1.D0)

      Suby(1,1) = dcmplx(1.D0/P3(1),0.0D0)
      Suby(2,2) = Suby(1,1)
      Suby(1,2) = -Suby(1,1)
      Suby(2,1) = Suby(1,2)
      Subj(1) = 0.0D0
      Subj(2) = 0.0D0

      numer = (L3-1)/3

      DO i = 1 , numer
         ika = 3*(i-1) + 2

         IF ( dabs(Om-2.D0*pi*P3(ika)).LE.Om*Epsi ) THEN

            IF ( Om.NE.0.D0 ) THEN
               sq1 = dcos(P3(ika+1))*dsqrt(8.0D0*P3(ika+2)/P3(1))
               sq2 = dsin(P3(ika+1))*dsqrt(8.0D0*P3(ika+2)/P3(1))
            ELSE
               sq1 = dsqrt(4.0D0*P3(ika+2)/P3(1))
               sq2 = 0.D0
            ENDIF

C  ACCORDING TO THE DEFINITION, A POWER SOURCE WITH P POWER DELIVERS P POWER TO THE LOAD
C  UNDER THE CONDITION THAT ITS IMPEDANCE EQUALS THE LOAD IMPEDANCE
C  FROM THIS, WE GET   J = 2 * SQRT(P / R_internal)
C  NEXT, FOR FREQUENCIES F != 0, WE MULTIPLY BY SQRT(2) TO GET
C  THE INSTANTANEOUS VALUE
C
C   INSTANTANEOUS CURRENT VALUES MUST BE PRESENT IN ALL MODELS !!!
C

            Subj(1) = dcmplx(sq1,sq2)
            Subj(2) = -Subj(1)
         ENDIF
      ENDDO
C     DEBUG                         SUBTRACE
      END

      SUBROUTINE emf(Om,P1,L1,P2,L2,P3,L3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Epsi , sq12
      INTEGER Kitu , L1 , L2 , L3 , Limerr
      COMMON /serv  / Epsi , Limerr , Kitu
      DOUBLE PRECISION P1 , P2 , P3 , pi , sq1 , sq2
      DOUBLE PRECISION Om
      DOUBLE COMPLEX Suby , Subj

      COMMON /subc  / Suby(15,15) , Subj(15)
      DIMENSION P1(L1) , P2(L2) , P3(L3)

      sq12 = dsqrt(2.D0)
      pi = 4.D0*datan(1.D0)

      Suby(1,1) = dcmplx(1.D0/P3(1),0.0D0)
      Suby(2,2) = Suby(1,1)
      Suby(1,2) = -Suby(1,1)
      Suby(2,1) = Suby(1,2)
      Subj(1) = 0.0D0
      Subj(2) = 0.0D0

      IF ( dabs(Om-2.D0*pi*P3(2)).LE.Om*Epsi ) THEN

         IF ( Om.NE.0.D0 ) THEN
            sq1 = sq12*P3(4)*dcos(P3(3))/P3(1)
            sq2 = sq12*P3(4)*dsin(P3(3))/P3(1)
         ELSE
            sq1 = P3(4)/P3(1)
            sq2 = 0.D0
         ENDIF

         Subj(1) = dcmplx(sq1,sq2)
         Subj(2) = -Subj(1)
      ENDIF
C
C      WRITE(6, 40)  ((I,J,SUBY(I,J),J=1,2),I=1,2)
C   40 FORMAT(2X,'LIBLIN EMF: '/
C     *       2X,'   SUBY(',I3,',',I3,')=',E13.6,',',E13.6/)

      END

      SUBROUTINE jdrive(Om,P1,L1,P2,L2,P3,L3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION Epsi , sq12
      INTEGER Kitu , L1 , L2 , L3 , Limerr
      COMMON /serv  / Epsi , Limerr , Kitu
      DOUBLE PRECISION P1 , P2 , P3 , pi , sq1 , sq2
      DOUBLE PRECISION Om
      DOUBLE COMPLEX Suby , Subj

      COMMON /subc  / Suby(15,15) , Subj(15)
      DIMENSION P1(L1) , P2(L2) , P3(L3)

C     PARAMETERS: COMMON-NON
C                INDIVIDUAL: G, F, FI, J
C

      sq12 = dsqrt(2.D0)
      pi = 4.D0*datan(1.D0)

      Suby(1,1) = dcmplx(P3(1),0.D0)
      Suby(2,2) = Suby(1,1)
      Suby(1,2) = -Suby(1,1)
      Suby(2,1) = Suby(1,2)
      Subj(1) = 0.0D0
      Subj(2) = 0.0D0

      IF ( dabs(Om-2.D0*pi*P3(2)).LE.Om*Epsi ) THEN

         IF ( Om.NE.0.D0 ) THEN
            sq1 = sq12*P3(4)*dcos(P3(3))
            sq2 = sq12*P3(4)*dsin(P3(3))
         ELSE
            sq1 = P3(4)
            sq2 = 0.D0
         ENDIF

         Subj(1) = dcmplx(sq1,sq2)
         Subj(2) = -Subj(1)
      ENDIF

      END
