Introns: Introns.o CodonTree.o
	g++ Introns.o CodonTree.o -o Introns

introns.o: Introns.cpp
	g++ -c Introns.cpp

CodonTree.o: CodonTree.h CodonTree.cpp
	g++ -c CodonTree.cpp

clean:
	rm *.o
	find . -perm +100 -type f -delete

run: Introns
	./Introns