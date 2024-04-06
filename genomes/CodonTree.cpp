#include "CodonTree.h"

//I could make a macro function and inflate all
//the set functions instead of doing this... just
//saying

// #define CONDITION (AC, I, N) \
// ( \
//     ((AC).codon[(I) + 1] == (N)) \
// )

// #define IACTION ()

// #define EACTION

bool CodonTree::setParent(CodonTree* parent) {
    this->parent = parent;
    return true;
}

CodonTree* CodonTree::getParent() {
    return parent;
}

bool CodonTree::setU(CodonTree* u) {
    this->u = u;
    return true;
}

CodonTree* CodonTree::getU() {
    return u;
}

bool CodonTree::setC(CodonTree* c) {
    this->c = c;
    return true;
}

CodonTree* CodonTree::getC() {
    return c;
}

bool CodonTree::setA(CodonTree* a) {
    this->a = a;
    return true;
}

CodonTree* CodonTree::getA() {
    return a;
}

bool CodonTree::setG(CodonTree* g) {
    this->g = g;
    return true;
}

CodonTree* CodonTree::getG() {
    return g;
}

bool CodonTree::setAcid(char* acid) {
    aminoAcid = acid;
    return true;
}

CodonTree::CodonTree() {
    setParent(nullptr);
    setU(nullptr);
    setC(nullptr);
    setA(nullptr);
    setG(nullptr);
    setAcid('\0');
}

CodonTree::CodonTree(CodonTree* parent) {
    setParent(parent);
    setU(nullptr);
    setC(nullptr);
    setA(nullptr);
    setG(nullptr);
    setAcid('\0');
}

CodonTree::CodonTree(CodonTree* parent, char* acid) {
    setParent(parent);
    setU(nullptr);
    setC(nullptr);
    setA(nullptr);
    setG(nullptr);
    setAcid(acid);
}

// CodonTree::~CodonTree() {
//     CodonTree* cursor;
//     for (cursor = getParent()) {

//     }
// }

//sample invoke
//aminoCodon(head, acPair, -1);
void CodonTree::addAminoCodon(CodonTree* cursor, AminoCodon acPair, int index) {
    //We need a null head that just acts as a starting point so we can see what the
    //NEXT amino acid is
    //right now, I don't know if this first conditional sets the amino acid at
    //the nucleotide just before the leaf, or at the leaf nucleotide!
    if (acPair.codon[index + 1] == '\0') {
        cursor->setAcid(acPair.aminoAcid);
    } else if (acPair.codon[index + 1] == 'u') {
        if (cursor->getU() != nullptr) {
            addAminoCodon(cursor->getU(), acPair, (index + 1));
        } else {
            CodonTree* add = new CodonTree(cursor);
            cursor->setU(add);
            addAminoCodon(cursor->getU(), acPair, (index + 1));
        }
    } else if (acPair.codon[index + 1] == 'c') {
        if (cursor->getC() != nullptr) {
            addAminoCodon(cursor->getC(), acPair, (index + 1));
        } else {
            CodonTree* add = new CodonTree(cursor);
            cursor->setC(add);
            addAminoCodon(cursor->getC(), acPair, (index + 1));
        }
    } else if (acPair.codon[index + 1] == 'a') {
        if (cursor->getA() != nullptr) {
            addAminoCodon(cursor->getA(), acPair, (index + 1));
        } else {
            CodonTree* add = new CodonTree(cursor);
            cursor->setA(add);
            addAminoCodon(cursor->getA(), acPair, (index + 1));
        }
    } else if (acPair.codon[index + 1] == 'g') {
        if (cursor->getG() != nullptr) {
            addAminoCodon(cursor->getG(), acPair, (index + 1));
        } else {
            CodonTree* add = new CodonTree(cursor);
            cursor->setG(add);
            addAminoCodon(cursor->getG(), acPair, (index + 1));
        }
    }
}