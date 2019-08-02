all: dfa

dfa:
	g++ -std=c++11 -Wall -fopenmp -O3 -march=native -o dfa dfa.cpp -g -msse2 -msse -maes aes.c
	cp dfa ../

clean:
	rm -f dfa
	rm -f ../dfa
	rm -f *.o *~

