.PHONY: all build test clean debug


build:
	g++ -I../../include \
		-L../../lib \
		-O2 -o main main.cpp \
		-lsgp4 -lsofa_c -lm 
	chmod +x main

debug:
	g++ -g -I../../include \
		-L../../lib \
		-O2 -o main main.cpp \
		-lsgp4 -lsofa_c -lm
	chmod +x main