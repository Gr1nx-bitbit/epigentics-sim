#ifndef CODONTREE_H
#define CODONTREE_H

#include <string>
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
        AminoCodon acPair;
       
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
        bool setAcid(AminoCodon acPair);
        AminoCodon getAcid(void);
        CodonTree* getToTop(CodonTree* cursor);
        void cascadeDelete(CodonTree* cursor);

    public:
        CodonTree();
        CodonTree(CodonTree* parent, int nucleotide);
        CodonTree(CodonTree* parent, int nucleotide, AminoCodon acPair);
        ~CodonTree();
        void addAminoCodon(CodonTree* cursor, AminoCodon acPair, int index);
        void displayTree(CodonTree* head, CodonTree* cursor, char step, std::string path);
        AminoCodon getAminoCodon(std::string codon, CodonTree* cursor);
};

#endif