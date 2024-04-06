#include <iostream>
#include <string>
#include "CodonTree.h"
#include "AminoCodon.h"
using namespace std;

int main(void) {
    AminoCodon* pair = new AminoCodon;
    pair->codon = "aug";
    pair->aminoAcid = "Methionene";
    CodonTree* head = new CodonTree();
    head->addAminoCodon(head, *pair, -1);
    pair->codon = "auu";
    pair->aminoAcid = "Isoleucine";
    head->addAminoCodon(head, *pair, -1);
    return 0;
}