INSTDIR = /usr/local/bin
CFGDIR  = /etc

all: compile

compile:
	make -C chanes
	make -C spice2hbl

install:
	install s2h $(INSTDIR)
	install hbl $(INSTDIR)
	install scripts/hbfree $(INSTDIR)
	install spice2hbl/s2h.cfg $(CFGDIR)


uninstall:
	rm -f $(INSTDIR)/s2h $(INSTDIR)/hbl  $(INSTDIR)/hbfree
	rm -f $(CFGDIR)/s2h.cfg

