#ifndef AMINDOCODON_H
#define AMINOCODON_H

#include <iostream>
#include "AminoCodon.h"

class CodonTree {
    private:
        CodonTree* parent;
        CodonTree* u;
        CodonTree* c;
        CodonTree* a;
        CodonTree* g;
        //std::string aminoAcid;
        char* aminoAcid;
        bool setParent(CodonTree* parent);
        CodonTree* getParent();
        bool setU(CodonTree* u);
        CodonTree* getU();
        bool setC(CodonTree* c);
        CodonTree* getC();
        bool setA(CodonTree* a);
        CodonTree* getA();
        bool setG(CodonTree* g);
        CodonTree* getG();
        bool setAcid(char* acid);

    public:
        CodonTree();
        CodonTree(CodonTree* parent);
        CodonTree(CodonTree* parent, char* acid);
        ~CodonTree();
        void addAminoCodon(CodonTree* cursor, AminoCodon acPair, int index);
};

#endif