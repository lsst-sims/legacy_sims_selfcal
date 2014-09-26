CFLAGS=		-Wall -g 
LDFLAGS=	-lm 

TARGETS=        ampcoordtran.h _ccs2amp_lib.so

all:		$(TARGETS)

ampcoordtran.h:	fp_coordtran_h2c.perl lsst_fp_coordtran.h
		./fp_coordtran_h2c.perl lsst_fp_coordtran.h

ccs2amp_lib.o:	ccs2amp_lib.c
		gcc -fPIC -c ccs2amp_lib.c

ccs2amp_lib_wrap.o:	ccs2amp_lib.i
		swig -python ccs2amp_lib.i 
		gcc -I$$PYTHON_DIR/include/python2.7 -fPIC -c ccs2amp_lib_wrap.c

_ccs2amp_lib.so: ccs2amp_lib.o ccs2amp_lib_wrap.o
		gcc -shared ccs2amp_lib.o ccs2amp_lib_wrap.o -o _ccs2amp_lib.so

clean:		
		rm -f *.o ampcoordtran.h ccs2amp_lib_wrap.c _ccs2amp_lib.so ccs2amp_lib.py

