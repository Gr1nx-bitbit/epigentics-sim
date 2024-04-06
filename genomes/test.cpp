#include <iostream>
#include <string>
#include "CodonTree.h"
#include "AminoCodon.h"
using namespace std;

int main(void) {
    string codon = "aug";
    AminoCodon pair;
    pair.codon = codon;
    pair.aminoAcid = "Methionene";
    CodonTree* head = new CodonTree();
    head->addAminoCodon(head, pair, -1);
    return 0;
}