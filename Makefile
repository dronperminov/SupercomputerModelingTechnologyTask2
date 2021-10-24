COMPILER=g++
FLAGS=-O3 -Wall -pedantic -fopenmp

all: main

main: main.cpp
	$(COMPILER) $(FLAGS) main.cpp -o main