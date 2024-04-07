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
    char yes;
    head->displayTree(head, head, yes, "");
    return 0;
}

// void txtParser(void) {
//     string text;
//     vector<AminoCodon> output;
//     ifstream myFile("aminoCodon.txt");

//     while(getline(myFile, text, '\n')) {
//         AminoCodon tmp = parse(text);
//         output.push_back(tmp);
//     }

//     for (int i = 0; i < output.size(); i++) {
//         if (output[i].startCodon) { 
//             cout << "Start amino acid: " << output[i].aminoAcid << endl;
//         } else if (output[i].terminationCodon) {
//             cout << "Termination acid: " << output[i].aminoAcid << endl;
//         }
//     }
// }