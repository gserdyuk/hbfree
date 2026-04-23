c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



c   these two function are for cast from Int*4 to char*4
      function c2i(c)
      integer*4 c2i,i
      character*4 c,char
      equivalence (i,char)
      char=c
      c2i=i
      return
      end

      function i2c(i)
      integer*4 i,idint
      character*4 char,i2c
      equivalence (idint,char)
      idint=i
      i2c=char
      return
      end
