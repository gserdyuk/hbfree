c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE LIBLIN(NAME,OM,NA1,LE1,NA2,LE2,NA3,LE3,NPOL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*4 NAMES(20),NAME(4)
      DOUBLE PRECISION PARAM
      COMMON/PARAMS/PARAM(1)
      DATA NAMES/'P2  ','E   ','LIB0','LL0 ','J   ','YTAB',           'S
     +TAB','LIB ','LANG','SHL ','KZ  ','XX  ',           'SMPL','MPL ','
     +DISC','INDS', 4*'XXXX' /
C
C
C
C     WRITE(6, 5) NAME,OM,NA1,LE1,NA2,LE2,NA3,LE3,NPOL
C   5 FORMAT(2X,'LIBLIN : N=',4A4,' OM=',E12.5,' NA1=',I5,' LE1=',I5/
C    *    11X,'NA2=',I5,' LE2=',I5,' NA3=',I5,' LE3=',I5,' NPOL=',I5)
C
C
      IF(NAME(1).NE.NAMES(3))GOTO 10
      IF(NAME(2).EQ.NAMES(4))CALL LINE1(OM,PARAM(NA1),LE1,PARAM(NA2),
     +                               LE2,PARAM(NA3),LE3)
   10 CONTINUE
      IF(NAME(1).EQ.NAMES(1))CALL RFFIPN(OM,PARAM(NA1),LE1,PARAM(NA2),
     +                                LE2,PARAM(NA3),LE3)
      IF(NAME(1).EQ.NAMES(2))CALL EMF   (OM,PARAM(NA1),LE1,PARAM(NA2),
     +                                LE2,PARAM(NA3),LE3)
      IF(NAME(1).EQ.NAMES(5))CALL JDRIVE(OM,PARAM(NA1),LE1,PARAM(NA2),
     +                                LE2,PARAM(NA3),LE3)
      IF(NAME(1).EQ.NAMES(6))CALL YTAB(OM,PARAM(NA1),LE1,PARAM(NA2),
     +                              LE2,PARAM(NA3),LE3,NPOL)
      IF(NAME(1).EQ.NAMES(7))CALL STAB(OM,PARAM(NA1),LE1,PARAM(NA2),
     +                              LE2,PARAM(NA3),LE3,NPOL)
      IF(NAME(1).EQ.NAMES(8).AND.NAME(2).EQ.NAMES(9)) CALL LANG(OM,
     +                               PARAM(NA1),LE1,PARAM(NA2),
     +                           LE2,PARAM(NA3),LE3,NPOL)
      IF(NAME(1).EQ.NAMES(8).AND.NAME(2).EQ.NAMES(13)) CALL SMPL(OM,
     +                                PARAM(NA1),LE1,PARAM(NA2),
     +                            LE2,PARAM(NA3),LE3,NPOL)
      IF(NAME(1).EQ.NAMES(8).AND.NAME(2).EQ.NAMES(14)) CALL MP(OM,
     +                              PARAM(NA1),LE1,PARAM(NA2),
     +                          LE2,PARAM(NA3),LE3,NPOL)
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

      IF(NAME(1).EQ.NAMES(16)) CALL INDSV(OM,PARAM(NA1),LE1,PARAM(NA2),
     +                                     LE2,PARAM(NA3),LE3)
C
      RETURN
C     DEBUG SUBTRACE
      END




      SUBROUTINE RFFIPN(OM,P1    ,L1,P2    ,L2,P3    ,L3)
C     CONVERT THE VALUE TO EPSI / REDUCE TO 0.1E-05 /
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /SERV/EPSI, LIMERR,KITU
      DOUBLE PRECISION P1,P2,P3,         SQ1,SQ2
      DOUBLE PRECISION   OM, PI
      DOUBLE COMPLEX SUBY       ,SUBJ

      COMMON/SUBC/SUBY(15,15),SUBJ(15)
      DIMENSION P1(L1),P2(L2),P3(L3)

      PI=4.D0*DATAN(1.D0)

      SUBY(1,1)=DCMPLX(1.D0/P3(1),0.0D0)
      SUBY(2,2)=SUBY(1,1)
      SUBY(1,2)=-SUBY(1,1)
      SUBY(2,1)=SUBY(1,2)
      SUBJ(1)=0.0D0
      SUBJ(2)=0.0D0

      NUMER=(L3-1)/3

      DO 10 I=1,NUMER
      IKA=3*(I-1)+2

      IF(DABS(OM-2.D0*PI*P3(IKA)).GT.OM*EPSI) GO TO 10

      IF(OM.NE.0.D0)GOTO 20
      SQ1=DSQRT(4.0D0*P3(IKA+2)/P3(1))
      SQ2=0.D0
      GOTO 30
   20 CONTINUE
      SQ1=DCOS(P3(IKA+1))*DSQRT(8.0D0*P3(IKA+2)/P3(1))
      SQ2=DSIN(P3(IKA+1))*DSQRT(8.0D0*P3(IKA+2)/P3(1))

C  ACCORDING TO THE DEFINITION, A POWER SOURCE WITH P POWER DELIVERS P POWER TO THE LOAD
C  UNDER THE CONDITION THAT ITS IMPEDANCE EQUALS THE LOAD IMPEDANCE
C  FROM THIS, WE GET   J = 2 * SQRT(P / R_internal)
C  NEXT, FOR FREQUENCIES F != 0, WE MULTIPLY BY SQRT(2) TO GET
C  THE INSTANTANEOUS VALUE
C
C   INSTANTANEOUS CURRENT VALUES MUST BE PRESENT IN ALL MODELS !!!
C
   30 CONTINUE

      SUBJ(1)=DCMPLX(SQ1,SQ2)
      SUBJ(2)=-SUBJ(1)
   10 CONTINUE
      RETURN
C     DEBUG                         SUBTRACE
      END

      SUBROUTINE EMF   (OM,P1    ,L1,P2    ,L2,P3    ,L3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/SERV/EPSI,LIMERR,KITU
      DOUBLE PRECISION P1,P2,P3,PI,SQ1,SQ2
      DOUBLE PRECISION   OM
      DOUBLE COMPLEX SUBY,SUBJ

      COMMON/SUBC/SUBY(15,15),SUBJ(15)
      DIMENSION P1(L1),P2(L2),P3(L3)

      SQ12=DSQRT(2.D0)
      PI=4.D0*DATAN(1.D0)

      SUBY(1,1)=DCMPLX(1.D0/P3(1),0.0D0)
      SUBY(2,2)=SUBY(1,1)
      SUBY(1,2)=-SUBY(1,1)
      SUBY(2,1)=SUBY(1,2)
      SUBJ(1)=0.0D0
      SUBJ(2)=0.0D0

      IF(DABS(OM-2.D0*PI*P3(2)).GT.OM*EPSI) GO TO 10

      IF(OM.NE.0.D0)GOTO 20
      SQ1=P3(4)/P3(1)
      SQ2=0.D0
      GOTO 30
   20 CONTINUE
      SQ1=SQ12*P3(4)*DCOS(P3(3))/P3(1)
      SQ2=SQ12*P3(4)*DSIN(P3(3))/P3(1)
   30 CONTINUE

      SUBJ(1)=DCMPLX(SQ1,SQ2)
      SUBJ(2)=-SUBJ(1)
   10 CONTINUE
C
C      WRITE(6, 40)  ((I,J,SUBY(I,J),J=1,2),I=1,2)
C   40 FORMAT(2X,'LIBLIN EMF: '/
C     *       2X,'   SUBY(',I3,',',I3,')=',E13.6,',',E13.6/)

      RETURN
      END

      SUBROUTINE JDRIVE  (OM,P1    ,L1,P2    ,L2,P3    ,L3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/SERV/EPSI,LIMERR,KITU
      DOUBLE PRECISION P1,P2,P3,PI,SQ1,SQ2
      DOUBLE PRECISION   OM
      DOUBLE COMPLEX SUBY,SUBJ

      COMMON/SUBC/SUBY(15,15),SUBJ(15)
      DIMENSION P1(L1),P2(L2),P3(L3)

C     PARAMETERS: COMMON-NON
C                INDIVIDUAL: G, F, FI, J
C

      SQ12=DSQRT(2.D0)
      PI=4.D0*DATAN(1.D0)

      SUBY(1,1)=DCMPLX(P3(1),0.D0)
      SUBY(2,2)=SUBY(1,1)
      SUBY(1,2)=-SUBY(1,1)
      SUBY(2,1)=SUBY(1,2)
      SUBJ(1)=0.0D0
      SUBJ(2)=0.0D0

      IF(DABS(OM-2.D0*PI*P3(2)).GT.OM*EPSI) GO TO 10

      IF(OM.NE.0.D0) GOTO 20
      SQ1=P3(4)
      SQ2=0.D0
      GOTO 30
   20 CONTINUE
      SQ1=SQ12*P3(4)*DCOS(P3(3))
      SQ2=SQ12*P3(4)*DSIN(P3(3))
   30 CONTINUE

      SUBJ(1)  =DCMPLX(SQ1,SQ2)
      SUBJ(2)  =-SUBJ(1)

   10 CONTINUE
      RETURN
      END
