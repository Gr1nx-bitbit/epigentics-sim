#include "../include/AminoCodon.h"
#include "../include/CodonTree.h"
#include <fstream>
#include <iostream>
//#include <vector>
using namespace std;

AminoCodon parse(string input);

int main(void) {
    CodonTree* head = new CodonTree();
    string aminoCodontext;
    ifstream myFile("aminoCodon.txt");

    while(getline(myFile, aminoCodontext, '\n')) {
        AminoCodon acPair = parse(aminoCodontext);
        head->addAminoCodon(head, acPair, -1);
    }

    char codon[3] = {'a', 'u', 'g'};
    string cool = codon;
    AminoCodon amino = head->getAminoCodon(cool, head);
    cout << amino.abbreviation << endl;
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