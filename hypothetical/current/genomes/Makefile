Introns: Introns.o CodonTree.o Node.o
	g++ Introns.o CodonTree.o Node.o -o Introns

Introns.o: Introns.cpp
	g++ -c Introns.cpp

Node.o: Node.cpp Node.h
	g++ -c Node.cpp

CodonTree.o: CodonTree.h CodonTree.cpp
	g++ -c CodonTree.cpp

clean:
	rm *.o
	find . -perm +100 -type f -delete

run: Introns
	./Introns