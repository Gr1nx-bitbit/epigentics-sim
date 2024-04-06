#ifndef AMINDOCODON_H
#define AMINOCODON_H

#include <iostream>
#include "AminoCodon.h"

class CodonTree {
    private:
        CodonTree* parent;
        CodonTree* u;
        bool nucU;
        CodonTree* c;
        bool nucC;
        CodonTree* a;
        bool nucA;
        CodonTree* g;
        bool nucG;
        //std::string aminoAcid;
        char* aminoAcid;
        bool setParent(CodonTree* parent);
        CodonTree* getParent();
        bool setU(CodonTree* u);
        CodonTree* getU();
        bool isU();
        bool setC(CodonTree* c);
        CodonTree* getC();
        bool isC();
        bool setA(CodonTree* a);
        CodonTree* getA();
        bool isA();
        bool setG(CodonTree* g);
        CodonTree* getG();
        bool isG();
        bool setAcid(char* acid);

    public:
        CodonTree();
        CodonTree(CodonTree* parent, int nucleotide);
        CodonTree(CodonTree* parent, int nucleotide, char* acid);
        ~CodonTree();
        void addAminoCodon(CodonTree* cursor, AminoCodon acPair, int index);
};

#endif