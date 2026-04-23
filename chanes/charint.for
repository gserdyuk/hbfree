!*==C2I.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



c   these two function are for cast from Int*4 to char*4
      FUNCTION c2i(C)
      INTEGER*4 c2i , i
      CHARACTER*4 C , char
      EQUIVALENCE (i,char)
      char = C
      c2i = i
      END

      FUNCTION i2c(I)
      INTEGER*4 I , idint
      CHARACTER*4 char , i2c
      EQUIVALENCE (idint,char)
      idint = I
      i2c = char
      END
