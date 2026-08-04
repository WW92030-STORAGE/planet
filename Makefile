# make > test.out

all:
	g++ main.cpp -o main -O3
	./main

ndebug:
	g++ -DNDEBUG main.cpp -o main -O3
	./main

