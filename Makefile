PREFIX = /home/dabergel/.local

LIBRARIES = libd_functions.a
HEADERS = d_gsloutput.h d_constants.h

d_gsloutput.o : d_gsloutput.h d_gsloutput.cpp
	g++ -c d_gsloutput.cpp

libd_functions.a : d_gsloutput.o
	ar -rs libd_functions.a d_gsloutput.o

install : libd_functions.a ${HEADERS}
	cp libd_functions.a ${PREFIX}/lib
	cp d_gsloutput.h ${PREFIX}/include
	cp d_constants.h ${PREFIX}/include

clean :
	rm ./d_gsloutput.o
	rm ./libd_functions.a
