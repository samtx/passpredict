.PHONY: all build test clean debug


build:
	g++ -Isgp4 -Isofa -O2 -o main main.cpp sgp4/SGP4.cpp -lm
	chmod +x main

debug:
	g++ -g -Isgp4 -Isofa -o main main.cpp sgp4/SGP4.cpp -lm
	chmod +x main
	gdb ./main
