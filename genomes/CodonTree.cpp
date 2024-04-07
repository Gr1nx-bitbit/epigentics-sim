#include "CodonTree.h"
#include <iostream>
using namespace std;
//#define DEBUG
#define DEBUGDISPLAY

//I could make a macro function and inflate all
//the set functions instead of doing this... just
//saying

// #define CONDITION (AC, I, N) \
// ( \
//     ((AC).codon[(I) + 1] == (N)) \
// )

// #define IACTION ()

// #define EACTION



/*
–––––––––––––––––––––
*                   *
*      PRIVATE      *
*                   *
–––––––––––––––––––––   
*/

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

bool CodonTree::isU() {
    return nucU;
}

bool CodonTree::setC(CodonTree* c) {
    this->c = c;
    return true;
}

CodonTree* CodonTree::getC() {
    return c;
}

bool CodonTree::isC() {
    return nucC;
}

bool CodonTree::setA(CodonTree* a) {
    this->a = a;
    return true;
}

CodonTree* CodonTree::getA() {
    return a;
}

bool CodonTree::isA() {
    return nucA;
}

bool CodonTree::setG(CodonTree* g) {
    this->g = g;
    return true;
}

CodonTree* CodonTree::getG() {
    return g;
}

bool CodonTree::isG() {
    return nucG;
}

bool CodonTree::setAcid(std::string acid) {
    aminoAcid = acid;
    return true;
}

string CodonTree::getAcid(void) {
    return aminoAcid;
}

/*
–––––––––––––––––––––
*                   *
*       PUBLIC      *
*                   *
–––––––––––––––––––––   
*/

CodonTree::CodonTree() {
    setParent(nullptr);
    setU(nullptr);
    nucU = false;
    setC(nullptr);
    nucC = false;
    setA(nullptr);
    nucA = false;
    setG(nullptr);
    nucG = false;
    setAcid("");
}

//parent is the parent of this node and int is the nucleotide
//to which this node will be. 1 is u, 2 is c, 3 is a, 4 is g.
CodonTree::CodonTree(CodonTree* parent, int nucleotide) {
    setParent(parent);
    setU(nullptr);
    setC(nullptr);
    setA(nullptr);
    setG(nullptr);
    setAcid("");
    if (nucleotide == 1) {
        nucU = true;
    } else if (nucleotide == 2) {
        nucC = true;
    } else if (nucleotide == 3) {
        nucA = true;
    } else if (nucleotide == 4) {
        nucG = true;
    }
}

CodonTree::CodonTree(CodonTree* parent, int nucleotide, std::string acid) {
    setParent(parent);
    setU(nullptr);
    setC(nullptr);
    setA(nullptr);
    setG(nullptr);
    setAcid(acid);
    if (nucleotide == 1) {
        nucU = true;
    } else if (nucleotide == 2) {
        nucC = true;
    } else if (nucleotide == 3) {
        nucA = true;
    } else if (nucleotide == 4) {
        nucG = true;
    }
}

//I'll figure out the destructor later. I'll have to use recursion to get through all the nodes!
// CodonTree::~CodonTree() {
//     CodonTree* cursor;
//     for (cursor = getParent(); cursor != nullptr; cursor = cursor->getParent());
//     CodonTree* traverse;
//     for (traverse = cursor; traverse->getU() != nullptr; traverse->getA()) {

//     }
// }

//sample invoke
//aminoCodon(head, acPair, -1);
void CodonTree::addAminoCodon(CodonTree* cursor, AminoCodon acPair, int index) {
    //We need a null head that just acts as a starting point so we can see what the
    //NEXT amino acid is
    //right now, I don't know if this first conditional sets the amino acid at
    //the nucleotide just before the leaf, or at the leaf nucleotide!
    #ifdef DEBUG
    cout << "The value at index " << (index + 1) << " is: " << acPair.codon[index + 1] << endl;
    #endif
    if (acPair.codon[index + 1] == '\0') {
        #ifdef DEBUG
        cout << "Setting amino acid to: " << acPair.aminoAcid << endl;
        #endif
        cursor->setAcid(acPair.aminoAcid);
    } else if (acPair.codon[index + 1] == 'u') {
        if (cursor->getU() != nullptr) {
            #ifdef DEBUG
            cout << "U pointer is defined. Calling add again with U pointer: " << cursor->getU() << endl;
            #endif
            addAminoCodon(cursor->getU(), acPair, (index + 1));
        } else {
            CodonTree* add = new CodonTree(cursor, 1);
            cursor->setU(add);
             #ifdef DEBUG
            cout << "U pointer is null. Allocating new node at: " << add << " and setting U pointer to it." << endl;
            #endif
            addAminoCodon(cursor->getU(), acPair, (index + 1));
        }
    } else if (acPair.codon[index + 1] == 'c') {
        if (cursor->getC() != nullptr) {
            #ifdef DEBUG
            cout << "C pointer is defined. Calling add again with C pointer: " << cursor->getC() << endl;
            #endif
            addAminoCodon(cursor->getC(), acPair, (index + 1));
        } else {
            CodonTree* add = new CodonTree(cursor, 2);
            cursor->setC(add);
             #ifdef DEBUG
            cout << "C pointer is null. Allocating new node at: " << add << " and setting C pointer to it." << endl;
            #endif
            addAminoCodon(cursor->getC(), acPair, (index + 1));
        }
    } else if (acPair.codon[index + 1] == 'a') {
        if (cursor->getA() != nullptr) {
            #ifdef DEBUG
            cout << "A pointer is defined. Calling add again with A pointer: " << cursor->getA() << endl;
            #endif
            addAminoCodon(cursor->getA(), acPair, (index + 1));
        } else {
            CodonTree* add = new CodonTree(cursor, 3);
            cursor->setA(add);
             #ifdef DEBUG
            cout << "A pointer is null. Allocating new node at: " << add << " and setting A pointer to it." << endl;
            #endif
            addAminoCodon(cursor->getA(), acPair, (index + 1));
        }
    } else if (acPair.codon[index + 1] == 'g') {
        if (cursor->getG() != nullptr) {
            #ifdef DEBUG
            cout << "G pointer is defined. Calling add again with G pointer: " << cursor->getG() << endl;
            #endif
            addAminoCodon(cursor->getG(), acPair, (index + 1));
        } else {
            CodonTree* add = new CodonTree(cursor, 4);
            cursor->setG(add);
             #ifdef DEBUG
            cout << "G pointer is null. Allocating new node at: " << add << " and setting G pointer to it." << endl;
            #endif
            addAminoCodon(cursor->getG(), acPair, (index + 1));
        }
    }
}

void CodonTree::displayTree(CodonTree* head, CodonTree* cursor, char step, std::string path) {
    path += step;
    //I need to add a path! I want to see where we end up at!
    if (cursor == nullptr) {
        return;
    } else {
            displayTree(head, cursor->getU(), 'u', path);
            displayTree(head, cursor->getC(), 'c', path);
            displayTree(head, cursor->getA(), 'a', path);
            displayTree(head, cursor->getG(), 'g', path);
            if (cursor->getAcid() != "") {
                cout << "Amino acid: " << cursor->getAcid() << endl;
                cout << "Path: " << path << endl;
            }
    }
}

string CodonTree::getAminoCodon(string codon, CodonTree* cursor) {
    for (int index = -1; true; index++) {
        if (codon[index + 1] == 'u') {
            if (cursor->getU() != nullptr) {
                cursor = cursor->getU();
            }
        } else if (codon[index + 1] == 'c') {
            if (cursor->getC() != nullptr) {
                cursor = cursor->getC();
            }
        } else if (codon[index + 1] == 'a') {
            if (cursor->getA() != nullptr) {
                cursor = cursor->getA();
            }
        } else if (codon[index + 1] == 'g') {
            if (cursor->getG() != nullptr) {
                cursor = cursor->getG();
            }
        } else if (codon[index + 1] == '\0') {
            return cursor->getAcid();
        }
    }
}