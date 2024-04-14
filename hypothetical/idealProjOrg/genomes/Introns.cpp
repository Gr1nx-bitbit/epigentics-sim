#include <iostream>
#include <AminoCodon/AminoCodon.h>
#include <GenePosition/GenePosition.h>
#include <RnaAndExcision/RnaAndExcision.h>
#include <Node/Node.h>
#include <CodonTree/CodonTree.h>
#include <fstream>
#include <vector>
using namespace std;

//#define DEBUG


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
vector<Node*> peptideSynthesis(string matureRna, int length, CodonTree* head);
void displayAbbreviations(Node* root);

/*
––––––––––––––––––––––––––
*   FUNCTION SIGNATURES  *
––––––––––––––––––––––––––
*/

int main(void) {
    //takes a genetic sequence and returns the mRNA sequence with the introns excised and compares the original mRNA to 
    //the solely exon mRNA
    string sequence = "cgagcacatcactaatataattggcgagggtccatctgataattgttgttgagttctgctaattgtcctcctcaatcgtccttcaagatcaattacgtccgatttctgtagagtcctgaatccaattaatcccgacttcggcggagattttcgtattgtttctgaaatcccttcaatgccttagatactagaaactgccttgatgcaatagtagtattgacagtatgtgccttgcaacttgtccttagaaccaacttgtaatggaacattcaaaccctgacgctacggcctactgaaagggaggaatatcgctgtaatgtttgtcccctactaagggaatgaccccgctccgcggttcccgtaaaccagtagacggcgcgtaccgaggtagtgctagatgtgacgcacaggcaccaacaactatataggttccattggatacgggtcagaccctcaatcagacatgatttggatccaccagtcaatctcggccgtgtatgacctactgtgtacactactggttagtaatttcatcagcaggagactcctatttcgtgtagcgtgtttgtttattgtcttatgcattcgtacaggcgatcatttactagaggcccttcacccagagttatggacgcaccgcgcgtatatcgagcgggatacacttgcggctgatggtcagcataatctttgggatgctatggctcgaatcgccttgcgaacagaatctgctgttttcaactgacatatcacctagaagagccccggtagctgtcccccgacgtgatgctcgaacttgatttcagtggggatcgattcttagtggaaggtcctactaatagcccacgtggacttgaacgagtcctatcatccatcaagagggggcgaatacctcacatagtgagcgtacggcatgttatacgaattgtaattgagctacgccccccgcgcaccgctatatatcttgccgcggagtgccaagcattcgtcccgtgaagtacccgcgttaaactatgtagggaagaaggccggcacccgtgtattcacgtcctcttcgaatcgcaaattgtaccttaattggttcattgaaaaacagttgaatgtgaaaatcgtaatgtgatgtactgttctcgctgacggtgatcgactgacccacggtaagatacacgattgattcttgccccggagacccttgtgttctccgccttctactatctcatgtagcacctgggcacgatggtggaccgcgctacagagaccaggggggcggtcctcgcgtgctggtcgtccgaccggcgaaagacctgaaagatttgtaggggttctaatatgtcaataactgatcgcgtgccaggcggtggaagctagtttaaggagcgtacgaccacgctagattactttgtagcgtaaattattatcgccaaccttattacccctagaccgtctagtagaacttagtttgaattccctggcgagtggtgtctagttttgttacatagtactcaacagtcgtcgacgttagggtcgcgtctcgaaagtttgttaccctttgacccgtcatcctagtataagatcttccagctagggtgtctgcgcgactcttgttgagatattactattccatgtgcttttgttagatcaagggctctaatccattgaggagattgaacgtggatttcattgagagagacctgtgttagtgtttgtgtagcacttaatgagggtgttgcggggcgcacaatccataaacacagataggtagctgcctatttgaccacggtgcgaactgtactggttcgccgcggaggtcccctagtgcagctcctgcgcgtcacctcgcaaacagattcctgcgtgaccacgacggtcttgttttacgccttatgacggtgatgcccccttgcaaatacgattaagggagagcggtgcggcccaagcctgcgagcgtgctcagagtatcagcacaacttgtgccccgccatggcgtgacgcccgcaggctttgctctctcgtcacgcgtcaagaacaaacaggcgtcgaacagttctttttccgcggtaagactgtcaaatagccatcgccatcaggcagaccaccatccgcgtcgtagaaataccagggatacgtttcttcttaacgtaattgtcgcttcgcggcaatagaactttgtgcgttcggacagtccgcgagacacgttttggaagtaggcgctatctgatttgtagaagagaagtacacaaccggccgcaggcccaacaggggtggcgtcgcgaagtcagtcgatcttcgccgcgaatcaaagtcccttaatagagaaaggagacttacgctgatttgtcgcgagcaagcagggaaatacgttcaatcaggttaggcatggagcttccgcggggttgccccccaccctccaaggcgggtgcggtttaaaaaaaccgtgcacggtataaaggtgcctgtaagcacacagtctgttctctacttgatcttcagctggaagtgaaagaaagtgaattctcagtcgtagggggtaaataggatcgcccggagcttatagtcgaggtgactcggtggcaacgatgatagatatcgacgtaatagtccctgtgcaatgatttgcttcctgggctggacgaatcactttgtccgatgagcttccgcaattaccgcgacgtctggcagtaagtgtacggtacggcgtgtgaaccaatcgtatgtaacatctgatatcgatcgttaatatcttatgccctggctgttcacctgcgagccagtacttcataactgccgataaatagcatcctgagtaatagagtattgtgatttgtgacttggtgaaaggtttccgatggggaagaaaagtttgaacgtacccgccacgtcaagcactgcagtctgccccggggtgcagagatgccacgaacgtcctattcgttctgaccgatgccaagtgaagaggtccctgcgaacctgagaggaggaaaatgaaaaagcgtgcttacaatt";
    string rnaSequence = dnaTOrna(sequence, sequence.length());
    Node* output = examine(rnaSequence, rnaSequence.length());
    rnaAndExcision preAmino = matureRNA(rnaSequence, rnaSequence.length(), output);
    //compare(rnaSequence, rnaSequence.length(), preAmino);

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

Node* examine(string sequence, int length) {
   bool start = false;
   Node* head = new Node();
   Node* cursor = head;
   GenePosition genCur;
   genCur.startIndex = -1;
   genCur.endIndex = -1;

    for (int i = 0; i < length; i++)
    {
        if (!start && (sequence[i] == 'g') && (sequence[i+1] == 'u') && (i < (length - 1)) && (i != genCur.endIndex)) {
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
    returnRNA.indecies = new GenePosition[index];
    string excisions[index];
    GenePosition pairs[index];
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
    bool semiColon = false;
    AminoCodon tmp;
    tmp.startCodon = false;
    tmp.terminationCodon = false;
    for (int index = 0; index < input.length(); index++) {
        if (input[index] == ':') {
            colon = true;            
        } else if (input[index] == ';') {
            semiColon = true;
        } else if (!colon && input[index] != ' ') {
            tmp.codon += input[index];
        } else if (semiColon) {
            if (input[index] != ' ') {
                tmp.abbreviation += input[index];
            }
        } else {
            if (input[index] == '!') {
                tmp.startCodon = true;
                tmp.terminationCodon = false;
            } else if (input[index] == '/') {
                tmp.terminationCodon = true;
                tmp.startCodon = false;
            } else {
                if (input[index] != ' ') {
                    tmp.aminoAcid += input[index];
                }
            }
        }
    }

    return tmp;
}

//takes in a string matureRNA and converts it into a string of
//polypeptide chains. Takes codons mapped to their respective
//amino acid and replaces the codon with the amino acid
vector<Node*> peptideSynthesis(string matureRna, int length, CodonTree* head) {
    vector<Node*> sequences;
    Node* node = new Node();
    sequences.push_back(node);
    Node* cursor = sequences[0];
    string codon = "";
    bool start = false;
    int index = 0;

    for (int i = 0; i < length; i++) {

        if ((i % 3 == 0) && (i != 0)) {
            AminoCodon amino = head->getAminoCodon(codon, head);

            #ifdef DEBUG
            if (amino.startCodon) {
                cout << "Start acid: " << amino.aminoAcid << endl;
            } else if (amino.terminationCodon) {
                cout << "Terminating acid: " << amino.aminoAcid << endl;
            }
            #endif

            if (amino.startCodon && (cursor->getAmino().aminoAcid == "") && !start) {
                cursor->setAmino(amino);
                start = true;

                if ((length - i) < 3) {
                    sequences[index]->setEnd(cursor);
                }
            } else if (amino.startCodon && (cursor->getAmino().aminoAcid != "") && !start) {
                Node* next = new Node(amino);
                next->setParent(cursor);
                cursor->setNext(next);
                cursor = cursor->getNext();
                start = true;

                if ((length - i) < 3) {
                    sequences[index]->setEnd(cursor);
                }
            } else if (amino.terminationCodon && (cursor->getAmino().aminoAcid == "") && start) {
                cursor->setAmino(amino);
                sequences[index]->setEnd(cursor);

                Node* increment = new Node();
                sequences.push_back(increment);
                cursor = increment;
                start = false;
            } else if (amino.terminationCodon && (cursor->getAmino().aminoAcid != "") && start) {
                Node* next = new Node(amino);
                next->setParent(cursor);
                cursor->setNext(next);
                sequences[index]->setEnd(cursor->getNext());

                Node* increment = new Node();
                sequences.push_back(increment);
                cursor = increment;
                start = false;
            } else if (!amino.startCodon && !amino.terminationCodon && start && (cursor->getAmino().aminoAcid == "")) {
                cursor->setAmino(amino);

                //need to check if this is the last codon in the matureRNA
                if ((length - i) < 3) {
                    sequences[index]->setEnd(cursor);
                }
            } else if (!amino.startCodon && !amino.terminationCodon && start && (cursor->getAmino().aminoAcid != "")) {
                Node* next = new Node(amino);
                next->setParent(cursor);
                cursor->setNext(next);
                cursor = cursor->getNext();

                if ((length - i) < 3) {
                    sequences[index]->setEnd(cursor);
                }
            }
            codon = matureRna[i];
        } else {
            codon += matureRna[i];
        }
    }

    //only the last sequence has a chance at having no termination codon. 
    //Since peptide sequences can get into the tens of thousands, I'd rather 
    //add a last member variale to Node so we can instantly look at whether 
    //or not the sequence is terminated with an end codon and therefore free it
    //instead of looping through it twice and getting worse performance.
    //If the end node is not a termination codon then erase the sequence

    //this deletes the acPair but not the whole thing for some reason. 
    //I should probably just make a destructor for the node!
    //that results in a segmentation fault :(
    if (!sequences[index]->getEnd()->getAmino().terminationCodon) {
        for (cursor = sequences[index]->getEnd(); cursor;) {
            #ifdef DEBUG
            cout << "Deleting " << cursor->getAmino().aminoAcid << endl;
            cout << cursor->getAmino().aminoAcid << endl;
            #endif

            Node* last = cursor;
            cursor = cursor->getParent();
            delete last;
            last = nullptr;

            #ifdef DEBUG
            if (last == nullptr) {
                cout << "Deleted successfully" << endl << endl;
            } else {
                cout << "Proof: " << last->getAmino().aminoAcid << endl;
                cout << "not deleted for whatever reason" << endl << endl;
            }
            #endif
        }
        sequences.pop_back();
    }

    return sequences;
}

void displayAbbreviations(Node* root) {
    for (Node* cursor = root; cursor; cursor = cursor->getNext()) {
        cout << cursor->getAmino().abbreviation;
    }
    cout << '\n';
}