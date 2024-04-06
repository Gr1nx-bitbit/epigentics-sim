introns: introns.o Node.o genes.o
	g++ introns.o Node.o genes.o

introns.o: introns.cpp
	g++ -c introns.cpp

Node.o: Node.cpp
	g++ -c Node.cpp

genes.o: genes.cpp
	g++ -c genes.cpp

clean:
	rm *.o
	find . -perm +100 -type f -delete