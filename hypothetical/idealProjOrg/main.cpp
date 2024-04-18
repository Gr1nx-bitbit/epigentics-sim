#include <iostream>
#include "AminoCodon.h"
#include "GenePosition.h"
#include "RnaAndExcision.h"
#include "Node.h"
#include "CodonTree.h"
#include "Introns.h"
#include <fstream>
#include <vector>
using namespace std;

//#define DEBUG


int main(void) {
    //takes a genetic sequence and returns the mRNA sequence with the introns excised and compares the original mRNA to 
    //the solely exon mRNA.
    string sequence = "cgagcacatcactaatataattggcgagggtccatctgataattgttgttgagttctgctaattgtcctcctcaatcgtccttcaagatcaattacgtccgatttctgtagagtcctgaatccaattaatcccgacttcggcggagattttcgtattgtttctgaaatcccttcaatgccttagatactagaaactgccttgatgcaatagtagtattgacagtatgtgccttgcaacttgtccttagaaccaacttgtaatggaacattcaaaccctgacgctacggcctactgaaagggaggaatatcgctgtaatgtttgtcccctactaagggaatgaccccgctccgcggttcccgtaaaccagtagacggcgcgtaccgaggtagtgctagatgtgacgcacaggcaccaacaactatataggttccattggatacgggtcagaccctcaatcagacatgatttggatccaccagtcaatctcggccgtgtatgacctactgtgtacactactggttagtaatttcatcagcaggagactcctatttcgtgtagcgtgtttgtttattgtcttatgcattcgtacaggcgatcatttactagaggcccttcacccagagttatggacgcaccgcgcgtatatcgagcgggatacacttgcggctgatggtcagcataatctttgggatgctatggctcgaatcgccttgcgaacagaatctgctgttttcaactgacatatcacctagaagagccccggtagctgtcccccgacgtgatgctcgaacttgatttcagtggggatcgattcttagtggaaggtcctactaatagcccacgtggacttgaacgagtcctatcatccatcaagagggggcgaatacctcacatagtgagcgtacggcatgttatacgaattgtaattgagctacgccccccgcgcaccgctatatatcttgccgcggagtgccaagcattcgtcccgtgaagtacccgcgttaaactatgtagggaagaaggccggcacccgtgtattcacgtcctcttcgaatcgcaaattgtaccttaattggttcattgaaaaacagttgaatgtgaaaatcgtaatgtgatgtactgttctcgctgacggtgatcgactgacccacggtaagatacacgattgattcttgccccggagacccttgtgttctccgccttctactatctcatgtagcacctgggcacgatggtggaccgcgctacagagaccaggggggcggtcctcgcgtgctggtcgtccgaccggcgaaagacctgaaagatttgtaggggttctaatatgtcaataactgatcgcgtgccaggcggtggaagctagtttaaggagcgtacgaccacgctagattactttgtagcgtaaattattatcgccaaccttattacccctagaccgtctagtagaacttagtttgaattccctggcgagtggtgtctagttttgttacatagtactcaacagtcgtcgacgttagggtcgcgtctcgaaagtttgttaccctttgacccgtcatcctagtataagatcttccagctagggtgtctgcgcgactcttgttgagatattactattccatgtgcttttgttagatcaagggctctaatccattgaggagattgaacgtggatttcattgagagagacctgtgttagtgtttgtgtagcacttaatgagggtgttgcggggcgcacaatccataaacacagataggtagctgcctatttgaccacggtgcgaactgtactggttcgccgcggaggtcccctagtgcagctcctgcgcgtcacctcgcaaacagattcctgcgtgaccacgacggtcttgttttacgccttatgacggtgatgcccccttgcaaatacgattaagggagagcggtgcggcccaagcctgcgagcgtgctcagagtatcagcacaacttgtgccccgccatggcgtgacgcccgcaggctttgctctctcgtcacgcgtcaagaacaaacaggcgtcgaacagttctttttccgcggtaagactgtcaaatagccatcgccatcaggcagaccaccatccgcgtcgtagaaataccagggatacgtttcttcttaacgtaattgtcgcttcgcggcaatagaactttgtgcgttcggacagtccgcgagacacgttttggaagtaggcgctatctgatttgtagaagagaagtacacaaccggccgcaggcccaacaggggtggcgtcgcgaagtcagtcgatcttcgccgcgaatcaaagtcccttaatagagaaaggagacttacgctgatttgtcgcgagcaagcagggaaatacgttcaatcaggttaggcatggagcttccgcggggttgccccccaccctccaaggcgggtgcggtttaaaaaaaccgtgcacggtataaaggtgcctgtaagcacacagtctgttctctacttgatcttcagctggaagtgaaagaaagtgaattctcagtcgtagggggtaaataggatcgcccggagcttatagtcgaggtgactcggtggcaacgatgatagatatcgacgtaatagtccctgtgcaatgatttgcttcctgggctggacgaatcactttgtccgatgagcttccgcaattaccgcgacgtctggcagtaagtgtacggtacggcgtgtgaaccaatcgtatgtaacatctgatatcgatcgttaatatcttatgccctggctgttcacctgcgagccagtacttcataactgccgataaatagcatcctgagtaatagagtattgtgatttgtgacttggtgaaaggtttccgatggggaagaaaagtttgaacgtacccgccacgtcaagcactgcagtctgccccggggtgcagagatgccacgaacgtcctattcgttctgaccgatgccaagtgaagaggtccctgcgaacctgagaggaggaaaatgaaaaagcgtgcttacaatt";
    string rnaSequence = dnaTOrna(sequence, sequence.length());
    Node* output = examine(rnaSequence, rnaSequence.length());
    rnaAndExcision preAmino = matureRNA(rnaSequence, rnaSequence.length(), output);
    //compare(rnaSequence, rnaSequence.length(), preAmino);

    //create a codonTree to store all the amino acids
    //so the problem is that the file isn't opening...
    CodonTree* head = new CodonTree();
    string aminoCodontext;
    ifstream myFile("aminoCodon.txt");

    //add all the amino acids and their path to the tree
    while(getline(myFile, aminoCodontext, '\n')) {
        AminoCodon acPair = parse(aminoCodontext);
        head->addAminoCodon(head, acPair, -1);
        cout << "hello!" << endl;
    }

    //this isn't working right now
    #ifdef DEBUG
    cout << "hello " << endl;
    char yes;
    head->displayTree(head, head, yes, "");
    cout << "bye" << endl;
    #endif
    vector<Node*> peptides = peptideSynthesis(preAmino.rna, preAmino.rna.length(), head);

    Node* cursor;
    int index = 0;

    #ifdef DEBUG
    cout << "Peptides size: " << peptides.size() << endl;
    

    
    for (cursor = peptides[index]; cursor;) {
        if (cursor->getAmino().startCodon) {
            cout << "Start amino: " << cursor->getAmino().aminoAcid << endl;
            cursor = cursor->getNext();
        } else if (cursor->getAmino().terminationCodon) {
            cout << "Termination codon: " << cursor->getAmino().aminoAcid << endl << endl;
            index++;
            cursor = peptides[index];
        } else {
            if (cursor->getAmino().aminoAcid != "")
            cout << "regular codon: " << cursor->getAmino().aminoAcid << endl;
            
            cursor = cursor->getNext();
        }
    }

    cout << "Peptides size: " << peptides.size() << endl;
    #endif

    //head->~CodonTree();
    //cout << peptides[4]->getAmino().aminoAcid << endl;
    for (int i = 0; i < peptides.size(); i++) {
        displayAbbreviations(peptides[i]);
        if (i != peptides.size() - 1) {
            cout << '\n';
        }
    }

    return 0;
}