.PHONY: all build test clean debug

SOURCE_FILES = 	*.cpp \
				sgp4/SGP4.cpp \
				sofa/dtf2d.c \
				sofa/d2dtf.c \
				sofa/d2tf.c \
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
				sofa/anpm.c \
				ast2body/ast2Body.cpp \
				ast2body/astMath.cpp \
				ast2body/astTime.cpp

OBJECT_FILES = 	observer.o \
				location.o \
				satellite.o \
				passpredict.o \
				orbit.o \
				SGP4.o \
				dtf2d.o \
				d2dtf.o \
				d2tf.o \
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
				anpm.o \
				ast2Body.o \
				astMath.o \
				astTime.o

build:
	g++ -I. -Isgp4 -Isofa -Iast2body -O2 -Wall -std=c++11 -c main.cpp $(SOURCE_FILES)
	g++ main.o $(OBJECT_FILES) -lm -o main
	chmod +x main

debug:
	g++ -g3 -O0 -I. -Isgp4 -Isofa -Iast2body -c main.cpp $(SOURCE_FILES)
	g++ main.o $(OBJECT_FILES) -lm -o main
	# chmod +x main
	# gdb ./main

test:
	g++ -I. -Isgp4 -Isofa -Iast2body -Itests -w -Wall -std=c++11 -c tests/*.cpp
	g++ test_*.o $(OBJECT_FILES) -lm -o test_passpredict
	./test_passpredict

debugtest:
	g++ -g3 -O0 -I. -Isgp4 -Isofa -Iast2body -Itests -w -Wall -std=c++11 -c tests/*.cpp
	g++ test_*.o $(OBJECT_FILES) -lm -o test_passpredict

clean:
	rm -r *.o
