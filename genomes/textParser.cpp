#include <fstream>
#include <iostream>
#include <vector>
#include "AminoCodon.h"
using namespace std;

AminoCodon parse(string input);

int main(int argC, char* argV[]) {
    //instead of hardcode the path we can pass it in as a CLI argument
    //cout << argV[1] << endl;
    string text;
    vector<AminoCodon> output;
    ifstream myFile("aminoCodon.txt");

    while(getline(myFile, text, '\n')) {
        AminoCodon tmp = parse(text);
        output.push_back(tmp);
    }

    for (int i = 0; i < output.size(); i++) {
        if (output[i].startCodon) { 
            cout << "Start amino acid: " << output[i].aminoAcid << endl;
        } else if (output[i].terminationCodon) {
            cout << "Termination acid: " << output[i].aminoAcid << endl;
        }
    }

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