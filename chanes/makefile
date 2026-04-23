
CC	=	gcc
CPPC	=	c++
LIBR	=	ar
RANLIB	=	ranlib
PROG	=	../hbl

#F77	=	g77
#FFLAGS = -Wall -fautomatic \
#		-g -fdebug-kludge \
#		-fwritable-strings -malign-double \
#		-O
#-fbounds-check

ifeq ($(FC),intel)
#ifx (prep.: source /opt/intel/oneapi/setvars.sh)
	F77	=	ifx
	#FFLAGS = -O2
	FFLAGS = -warn all,noerrors,nointerfaces,noexternals,nointerfaces,nounused -O2
	#FFLAGS = -g -check uninit
	#LDLAGS = -g -check uninit
else
#gfortran
	F77	=	gfortran
	FFLAGS = -Wall -Wextra -Wuninitialized \
			-O2 \
			-malign-double \
			-fdec-char-conversions \
			-fallow-argument-mismatch \
			-fdefault-real-8 \
			-fdefault-double-8 \
			-freal-4-real-8
#			-std=legacy \
#		-g
endif

OBJ_DIR = obj


F_SRC =  \
	SPAGged/charint.f	\
	SPAGged/decode.f 	\
	SPAGged/detsyn.f 	\
	SPAGged/former.f 	\
	SPAGged/frequen.f 	\
	SPAGged/ftmas2.f 	\
	SPAGged/harm.f  	\
	SPAGged/initi.f 	\
	SPAGged/inkoor.f 	\
	SPAGged/koord.f  	\
	SPAGged/krdchk.f	\
	SPAGged/lena.f  	\
	SPAGged/len.f 	\
	SPAGged/lineq1.f 	\
	SPAGged/luback.f 	\
	SPAGged/lucan.f 	\
	SPAGged/lufrw.f 	\
	SPAGged/luslv.f 	\
	SPAGged/main.f  \
	SPAGged/otchm.f \
	SPAGged/packm.f \
	SPAGged/regul.f \
	SPAGged/restor.f	\
	SPAGged/sort.f  \
	SPAGged/sortuz.f	\
	SPAGged/stback.f	\
	SPAGged/stepfr.f	\
	SPAGged/sumdif.f	\
	SPAGged/suprim.f	\
	SPAGged/topo.f	\
	SPAGged/topoin.f	\
	SPAGged/zgradie.f	\
	SPAGged/ziny.f	\
	SPAGged/zlserch.f	\
	SPAGged/zmache.f	\
	SPAGged/znef.f	\
	SPAGged/zneinck.f	\
	SPAGged/znejac.f	\
	SPAGged/znewsol.f	\
	SPAGged/zstop0.f	\
	SPAGged/zstop.f 

LIB_DIR = ./libelem

LIB_SRC = \
	SPAGged/libelem/biptr.f \
	SPAGged/libelem/cbarr.f \
	SPAGged/libelem/cdiff.f \
	SPAGged/libelem/clin.f \
	SPAGged/libelem/cpoly.f \
	SPAGged/libelem/curt.f \
	SPAGged/libelem/cusd.f \
	SPAGged/libelem/discont.f \
	SPAGged/libelem/icujunc.f \
	SPAGged/libelem/icupoly.f \
	SPAGged/libelem/indsv.f \
	SPAGged/libelem/junc.f \
	SPAGged/libelem/lib0.f \
	SPAGged/libelem/liblin.f \
	SPAGged/libelem/libmod.f \
	SPAGged/libelem/lin.f \
	SPAGged/libelem/mdsch.f \
	SPAGged/libelem/mpl.f \
	SPAGged/libelem/poly5.f \
	SPAGged/libelem/shleif.f \
	SPAGged/libelem/stab.f \
	SPAGged/libelem/svlutl.f \
	SPAGged/libelem/test.f \
	SPAGged/libelem/ytab.f 

F_OBJ = \
	$(OBJ_DIR)/charint.o \
	$(OBJ_DIR)/decode.o \
	$(OBJ_DIR)/detsyn.o \
	$(OBJ_DIR)/former.o \
	$(OBJ_DIR)/frequen.o \
	$(OBJ_DIR)/ftmas2.o \
	$(OBJ_DIR)/harm.o \
	$(OBJ_DIR)/initi.o \
	$(OBJ_DIR)/inkoor.o \
	$(OBJ_DIR)/koord.o \
	$(OBJ_DIR)/krdchk.o \
	$(OBJ_DIR)/lena.o \
	$(OBJ_DIR)/len.o \
	$(OBJ_DIR)/lineq1.o \
	$(OBJ_DIR)/luback.o \
	$(OBJ_DIR)/lucan.o \
	$(OBJ_DIR)/lufrw.o \
	$(OBJ_DIR)/luslv.o \
	$(OBJ_DIR)/main.o \
	$(OBJ_DIR)/otchm.o \
	$(OBJ_DIR)/packm.o \
	$(OBJ_DIR)/regul.o \
	$(OBJ_DIR)/restor.o \
	$(OBJ_DIR)/sort.o \
	$(OBJ_DIR)/sortuz.o \
	$(OBJ_DIR)/stback.o \
	$(OBJ_DIR)/stepfr.o \
	$(OBJ_DIR)/sumdif.o \
	$(OBJ_DIR)/suprim.o \
	$(OBJ_DIR)/topo.o \
	$(OBJ_DIR)/topoin.o \
	$(OBJ_DIR)/zgradie.o \
	$(OBJ_DIR)/ziny.o \
	$(OBJ_DIR)/zlserch.o \
	$(OBJ_DIR)/zmache.o \
	$(OBJ_DIR)/znef.o \
	$(OBJ_DIR)/zneinck.o \
	$(OBJ_DIR)/znejac.o \
	$(OBJ_DIR)/znewsol.o \
	$(OBJ_DIR)/zstop0.o \
	$(OBJ_DIR)/zstop.o 
	
LIB_OBJ = \
	$(OBJ_DIR)/biptr.o \
	$(OBJ_DIR)/cbarr.o \
	$(OBJ_DIR)/cdiff.o \
	$(OBJ_DIR)/clin.o \
	$(OBJ_DIR)/cpoly.o \
	$(OBJ_DIR)/curt.o \
	$(OBJ_DIR)/cusd.o \
	$(OBJ_DIR)/discont.o \
	$(OBJ_DIR)/icujunc.o \
	$(OBJ_DIR)/icupoly.o \
	$(OBJ_DIR)/indsv.o \
	$(OBJ_DIR)/junc.o \
	$(OBJ_DIR)/lib0.o \
	$(OBJ_DIR)/liblin.o \
	$(OBJ_DIR)/libmod.o \
	$(OBJ_DIR)/lin.o \
	$(OBJ_DIR)/mdsch.o \
	$(OBJ_DIR)/mpl.o \
	$(OBJ_DIR)/poly5.o \
	$(OBJ_DIR)/shleif.o \
	$(OBJ_DIR)/stab.o \
	$(OBJ_DIR)/svlutl.o \
	$(OBJ_DIR)/test.o \
	$(OBJ_DIR)/ytab.o 

.all: $(PROG)

$(PROG): $(F_OBJ) $(LIB_OBJ)
		$(F77) $(LDLAGS) -o $@ $^
		
$(F_OBJ): $(OBJ_DIR)/%.o: SPAGged/%.f
		$(F77) -c $(FFLAGS) -o $@ $< 

$(LIB_OBJ): $(OBJ_DIR)/%.o: SPAGged/$(LIB_DIR)/%.f
		$(F77) -c $(FFLAGS) -o $@ $< 

touch:
	touch $(F_SRC) $(LIB_SRC) $(QUP_SRC)
		
clean:
	rm -f $(F_OBJ) $(LIB_OBJ)

#includes:
$(OBJ_DIR)/former.o:  circuit.i charint.i
$(OBJ_DIR)/frequen.o: circuit.i charint.i
$(OBJ_DIR)/len.o:     circuit.i charint.i
$(OBJ_DIR)/main.o:    circuit.i
$(OBJ_DIR)/sortuz.o:  circuit.i
$(OBJ_DIR)/suprim.o:  circuit.i
$(OBJ_DIR)/topo.o:    circuit.i
$(OBJ_DIR)/topoin.o:  circuit.i
$(OBJ_DIR)/znef.o:    circuit.i funcsize.i      
$(OBJ_DIR)/znejac.o:  circuit.i funcsize.i
$(OBJ_DIR)/inkoor.o:  funcsize.i
