#include <iostream>
#include "AminoCodon.h"
#include "CodonTree.h"
#include <fstream>
using namespace std;

struct genePosition {
    int startIndex;
    int endIndex;
};

struct rnaAndExcision {
    string* excised;
    genePosition* indecies;
    string rna;
};

class Node {
    private:
        genePosition gene;
        bool gen;
        string aminoAcid;
        Node* next;
    public:
        Node() {
            gene.startIndex = 0;
            gene.endIndex = 0;
            gen = false;
            aminoAcid = "";
            next = nullptr;
        }

        Node(genePosition geen) {
            gene = geen;
            gen = true;
            aminoAcid = "";
            next = nullptr;
        }

        Node(string acid) {
            gene.startIndex = 0;
            gene.endIndex = 0;
            gen = false;
            aminoAcid = acid;
            next = nullptr;
        }

        // Node(genePosition geen, Node* nxt) {
        //     gene = geen;
        //     isGene = true;
        //     next = nxt;
        // }

        genePosition getGene(void) {
            return gene;
        }

        void setGene(genePosition geen) {
            gene  = geen;
            gen = true;
        }

        bool isGene(void) {
            return gen;
        }

        string getAmino(void) {
            return aminoAcid;
        }

        void setAmino(string acid) {
            aminoAcid = acid;
        }

        Node* getNext(void) {
            return next;
        }

        void setNext(Node* nxt) {
            next = nxt;
        }
};

/*
––––––––––––––––––––––––––
*   FUNCTION SIGNATURES  *
––––––––––––––––––––––––––
*/

Node* examine(string sequence, int length);
rnaAndExcision matureRNA(string sequence, int length, Node* head);
void compare(string premature, int length, rnaAndExcision preAmino);
string dnaTOrna(string dna, int length);
AminoCodon parse(string input);
Node* peptideSynthesis(string matureRna, int length, CodonTree* head);

/*
––––––––––––––––––––––––––
*   FUNCTION SIGNATURES  *
––––––––––––––––––––––––––
*/

int main(void) {
    //takes a genetic sequence and returns the mRNA sequence with the introns excised and compares the original mRNA to 
    //the solely exon mRNA
    string sequence = "caggaatactggtcaagaggttagggcgcaccgtccgtcggacatgacacctgtagagctatgtagtgtacggtatcgatgcagggcagggcggtgcgta";
    string rnaSequence = dnaTOrna(sequence, sequence.length());
    Node* output = examine(rnaSequence, rnaSequence.length());
    rnaAndExcision preAmino = matureRNA(rnaSequence, rnaSequence.length(), output);
    compare(rnaSequence, rnaSequence.length(), preAmino);

    //create a codonTree to store all the amino acids
    CodonTree* head = new CodonTree();
    string aminoCodontext;
    ifstream myFile("aminoCodon.txt");

    //add all the amino acids and their path to the tree
    while(getline(myFile, aminoCodontext, '\n')) {
        AminoCodon acPair = parse(aminoCodontext);
        head->addAminoCodon(head, acPair, -1);
    }

    //this isn't working right now
    Node* peptides = peptideSynthesis(preAmino.rna, preAmino.rna.length(), head);
    Node* cursor;
    for (cursor = peptides; cursor; cursor = cursor->getNext()) {
        cout << cursor->getAmino() << endl;
    }

    return 0;
}

Node* examine(string sequence, int length) {
   bool start = false;
   Node* head = new Node();
   Node* cursor = head;
   genePosition genCur;
   genCur.startIndex = -1;
   genCur.endIndex = -1;

    for (int i = 0; i < length; i++)
    {
        if (!start && (sequence[i] == 'g') && (sequence[i+1] == 'u') && (i < (length - 1)) && (i != genCur.endIndex)) { //hey
            start = true;
            genCur.startIndex = i;
        } else if (start && (sequence[i] == 'a') && (sequence[i+1] == 'g') && (i < (length - 1)) && (i != genCur.startIndex)) {
            genCur.endIndex = i+1;
            if (!cursor->isGene()) {
                cursor->setGene(genCur);
            } else {
                Node* next = new Node(genCur);
                cursor->setNext(next);
                cursor = cursor->getNext();
            }
            start = false;
        }
    }

    return head;
}

rnaAndExcision matureRNA(string sequence, int length, Node* head) {
    string matureRNA;
    Node* cursor;
    int index = 0;
    for (cursor = head; cursor; cursor = cursor->getNext()) {
        index++;
    }
    rnaAndExcision returnRNA;
    returnRNA.excised = new string[index];
    returnRNA.indecies = new genePosition[index];
    string excisions[index];
    genePosition pairs[index];
    index = 0;
    for (cursor = head; cursor; cursor = cursor->getNext()) {
        pairs[index] = cursor->getGene();
        index++;
    }

    
    index = 0;
    for (int i = 0; i < length; i++) {
        if (i == pairs[index].startIndex) { 
            excisions[index] += sequence[i];
            returnRNA.indecies[index].startIndex = i;
            continue;
        }
        if ((i > pairs[index].startIndex) && (i < pairs[index].endIndex)) {
            excisions[index] += sequence[i];
            continue;
        }
        if (i == pairs[index].endIndex) {
            excisions[index] += sequence[i];
            returnRNA.indecies[index].endIndex = i;
            index++;
            continue;
        }
        matureRNA += sequence[i];
    }
    
    index = 0;
    for (string render : excisions) {
        returnRNA.excised[index] = render;
        index++;
    }

    returnRNA.rna = matureRNA;
    return returnRNA;
}

void compare(string premature, int length, rnaAndExcision preAmino) {
    cout << premature << endl;
    int index = 0;
    for (int i = 0; i < length; i++) {
        if (i == preAmino.indecies[index].startIndex) {
            cout << '^';
        } else if (i > preAmino.indecies[index].startIndex && i < preAmino.indecies[index].endIndex) {
            cout << '-';
        } else if (i == preAmino.indecies[index].endIndex) {
            cout << '^';
            index++;
        } else {
            cout << " ";
        }
        if (i == (length - 1)) {
            cout << endl;
        }
    }

    index = 0;
    int rnaIndex = 0;
    for (int i = 0; i < length; i++) {
        if (i == preAmino.indecies[index].startIndex) {
            cout << ' ';
        } else if (i > preAmino.indecies[index].startIndex && i < preAmino.indecies[index].endIndex) {
            cout << ' ';
        } else if (i == preAmino.indecies[index].endIndex) {
            cout << ' ';
            index++;
        } else {
            cout << preAmino.rna[rnaIndex];
            rnaIndex++;
            if (rnaIndex == (preAmino.rna.length())) {
                cout << endl;
            }
        }
    }
}

string dnaTOrna(string dna, int length) {
    string prematureRNA;
    for (int i = 0; i < length; i++) {
        if (dna[i] == 'a') {
            prematureRNA += 'u';
        } else if (dna[i] == 't') {
            prematureRNA += 'a';
        } else if (dna[i] == 'c') {
            prematureRNA += 'g';
        } else {
            prematureRNA += 'c';
        }
    }

    return prematureRNA;
}

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

//takes in a string matureRNA and converts it into a string of
//polypeptide chains. Takes codons mapped to their respective
//amino acid and replaces the codon with the amino acid
Node* peptideSynthesis(string matureRna, int length, CodonTree* head) {
    Node* node = new Node();
    Node* cursor = node;
    string codon = "";

    for (int i = 0; i < length; i++) {
        if ((i % 3 == 0) && (i != 0)) {
            if (cursor->getAmino() == "") {
                string amino = head->getAminoCodon(codon, head);
                cursor->setAmino(amino);
            } else {
                string amino = head->getAminoCodon(codon, head);
                Node* next = new Node(amino);
                cursor->setNext(next);
                cursor = cursor->getNext();
            }
            codon = matureRna[i];
        } else {
            codon += matureRna[i];
        }
    }

    return node;
}