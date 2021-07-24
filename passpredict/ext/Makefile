.PHONY: all build test clean debug

SOURCE_FILES = 	*.cpp \
				sgp4/SGP4.cpp \
				sofa/gmst82.c \
				sofa/utctai.c \
				sofa/taitt.c \
				sofa/obl80.c \
				sofa/nut80.c \
				sofa/eqeq94.c \
				sofa/anp.c \
				sofa/jd2cal.c \
				sofa/dat.c \
				sofa/cal2jd.c \
				sofa/anpm.c

OBJECT_FILES = 	observer.o \
				location.o \
				satellite.o \
				orbit.o \
				SGP4.o \
				gmst82.o \
				utctai.o \
				taitt.o \
				obl80.o \
				nut80.o \
				eqeq94.o \
				anp.o \
				jd2cal.o \
				dat.o \
				cal2jd.o \
				anpm.o

build:
	g++ -I. -Isgp4 -Isofa -O2 -std=c++11 -c main.cpp $(SOURCE_FILES)
	g++ main.o $(OBJECT_FILES) -lm -o main
	chmod +x main

debug:
	g++ -g -I. -Isgp4 -Isofa -c main.cpp $(SOURCE_FILES)
	g++ main.o $(OBJECT_FILES) -lm -o main
	chmod +x main
	gdb ./main

test:
	g++ -I. -Isgp4 -Isofa -Itests -w -std=c++11 -c tests/*.cpp
	g++ test_*.o $(OBJECT_FILES) -lm -o test_passpredict
	./test_passpredict

clean:
	rm -r *.o
