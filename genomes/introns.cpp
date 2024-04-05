#include <iostream>
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
        Node* next;
    public:
        Node() {
            gene.startIndex = 0;
            gene.endIndex = 0;
            gen = false;
            next = nullptr;
        }

        Node(genePosition geen) {
            gene = geen;
            gen = true;
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

        Node* getNext(void) {
            return next;
        }

        void setNext(Node* nxt) {
            next = nxt;
        }
};

Node* examine(string sequence, int length);
rnaAndExcision matureRNA(string sequence, int length, Node* head);
void compare(string premature, int length, rnaAndExcision preAmino);

int main(void) {
    string sequence = "acgatttcgatgtgccggtttattatatagccgcggcccaccttatagccgccggtataccaccgggcattggctacctcgcatggcaacgattctctca";
    Node* output = examine(sequence, sequence.length());
    rnaAndExcision preAmino = matureRNA(sequence, sequence.length(), output);
    cout << "intial sequence: " << sequence << endl;
    cout << "introns excised: " << preAmino.excised[0] << " and: " << preAmino.excised[1] << endl;
    cout << "new rna sequence: " << preAmino.rna << endl;
    //compare(sequence, sequence.length(), preAmino);
    return 0;
}

Node* examine(string sequence, int length) {
   bool start = false;
   Node* head = new Node();
   Node* cursor = head;
   genePosition genCur;

    for (int i = 0; i < length; i++)
    {
        if (!start && (sequence[i] == 'g') && (sequence[i+1] == 't') && (i < (length - 1))) {
            start = true;
            genCur.startIndex = i;
        } else if (start && (sequence[i] == 'a') && (sequence[i+1]) && (i < (length - 1))) {
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
        }
        if (i == (length - 1)) {
            cout << endl;
        }
    }
    cout << preAmino.rna << endl;
}