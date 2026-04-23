
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
	#FFLAGS = -g -check uninit -O1
	#LDLAGS = -g -check uninit
else
#gfortran
	F77	=	gfortran
	FFLAGS = -Wall -Wextra -Wuninitialized \
			-O2 \
			-malign-double \
			-fdec-char-conversions \
			-fallow-argument-mismatch
#			-std=legacy \
#		-g
endif

OBJ_DIR = obj


F_SRC =  \
	charint.for	\
	decode.for 	\
	detsyn.for 	\
	former.for 	\
	frequen.for 	\
	ftmas2.for 	\
	harm.for  	\
	initi.for 	\
	inkoor.for 	\
	koord.for  	\
	krdchk.for	\
	lena.for  	\
	len.for 	\
	lineq1.for 	\
	luback.for 	\
	lucan.for 	\
	lufrw.for 	\
	luslv.for 	\
	main.for  \
	otchm.for \
	packm.for \
	regul.for \
	restor.for	\
	sort.for  \
	sortuz.for	\
	stback.for	\
	stepfr.for	\
	sumdif.for	\
	suprim.for	\
	topo.for	\
	topoin.for	\
	zgradie.for	\
	ziny.for	\
	zlserch.for	\
	zmache.for	\
	znef.for	\
	zneinck.for	\
	znejac.for	\
	znewsol.for	\
	zstop0.for	\
	zstop.for 

LIB_DIR = ./libelem

LIB_SRC = \
	libelem/biptr.for \
	libelem/cbarr.for \
	libelem/cdiff.for \
	libelem/clin.for \
	libelem/cpoly.for \
	libelem/curt.for \
	libelem/cusd.for \
	libelem/discont.for \
	libelem/icujunc.for \
	libelem/icupoly.for \
	libelem/indsv.for \
	libelem/junc.for \
	libelem/lib0.for \
	libelem/liblin.for \
	libelem/libmod.for \
	libelem/lin.for \
	libelem/mdsch.for \
	libelem/mpl.for \
	libelem/poly5.for \
	libelem/shleif.for \
	libelem/stab.for \
	libelem/svlutl.for \
	libelem/test.for \
	libelem/ytab.for 

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
		
$(F_OBJ): $(OBJ_DIR)/%.o: %.for
		$(F77) -c $(FFLAGS) -o $@ $< 

$(LIB_OBJ): $(OBJ_DIR)/%.o: libelem/%.for
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
