.PHONY: all build test clean


build:
	g++ -I../../include \
		-L../../lib -lsgp4 -lsofa_c -lm \
		-O2 -o main main.cpp
	chmod +x main