PREFIX = /home/dabergel/.local

CFLAGS = -fopenmp

LIBRARIES = libd_functions.a
HEADERS = d_gsloutput.h d_constants.h

OBJECTS = d_functions.o d_gsloutput.o

d_functions.o : d_functions.h d_functions.cpp
	g++ ${CFLAGS} -c d_functions.cpp

d_gsloutput.o : d_gsloutput.h d_gsloutput.cpp
	g++ ${CFLAGS} -c d_gsloutput.cpp

libd_functions.a : ${OBJECTS}
	@for stub in ${OBJECTS} ; do \
		ar -rsv libd_functions.a $$stub ; \
	done

# Should add something here to check if destinations exist
install : libd_functions.a ${HEADERS}
	cp libd_functions.a ${PREFIX}/lib
	cp d_gsloutput.h ${PREFIX}/include
	cp d_functions.h ${PREFIX}/include
	cp d_constants.h ${PREFIX}/include

clean :
	@for stub in ${OBJECTS} ; do \
		if [ -e $$stub ] ; then rm -v $$stub ; fi ; \
	done
	@if [ -e libd_functions.a ] ; then rm -v ./libd_functions.a ; fi
	@if [ -e testing ] ; then rm -v ./testing ; fi

testing : testing.cpp ${OBJECTS}
	g++ -o testing ${CFLAGS} testing.cpp ${OBJECTS} -lgsl -lgslcblas
