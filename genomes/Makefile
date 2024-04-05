introns.o: introns.cpp
	g++ -c introns.cpp

run: introns.o
	g++ *.o -o introns

clean:
	rm *.out
	find . -perm +100 -type f -delete