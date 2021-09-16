# Build and install the SGP4 functions as a static library

SGP4_HEADER = SGP4.h

build:
	g++ -c -pedantic -Wall -W -fPIC -O2 -o SGP4.o SGP4.cpp
	ar ru libsgp4.a SGP4.o
	cp SGP4.h ../../include/sgp4.h
	cp libsgp4.a ../../lib/libsgp4.a
