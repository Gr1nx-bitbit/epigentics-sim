#include "CodonTree.h"
#include "AminoCodon.h"
#include <fstream>
#include <iostream>
//#include <vector>
using namespace std;

AminoCodon parse(string input);

int main(void) {
    CodonTree* head = new CodonTree();
    string text;
    ifstream myFile("aminoCodon.txt");

    while(getline(myFile, text, '\n')) {
        AminoCodon acPair = parse(text);
        head->addAminoCodon(head, acPair, -1);
    }

    char yes;
    head->displayTree(head, head, yes, "");
    return 0;
}

//this function takes in a line from the txt file 
//containing the codon and amino acid pairs and returns
//the pair in a AminoCodon
AminoCodon parse(string input) {
    bool colon = false;
    AminoCodon tmp;
    tmp.startCodon = false;
    tmp.terminationCodon = false;
    for (int index = 0; index < input.length(); index++) {
        if (input[index] == ':') {
            colon = true;            
        } else if (!colon && input[index] != ' ') {
            tmp.codon += input[index];
        } else {
            if (input[index] == '!') {
                tmp.startCodon = true;
                tmp.terminationCodon = false;
            } else if (input[index] == '/') {
                tmp.terminationCodon = true;
                tmp.startCodon = false;
            } else {
                tmp.aminoAcid += input[index];
            }
        }
    }

    return tmp;
}