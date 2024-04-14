#include "Node.h"

Node::Node() {
    gene.startIndex = 0;
    gene.endIndex = 0;
    gen = false;
    acPair.aminoAcid = "";
    acPair.codon = "";
    acPair.startCodon = false;
    acPair.terminationCodon = false;
    next = nullptr;
    parent = nullptr;
    end = nullptr;
}

Node::Node(GenePosition geen) {
    gene = geen;
    gen = true;
    acPair.aminoAcid = "";
    acPair.codon = "";
    acPair.startCodon = false;
    acPair.terminationCodon = false;
    next = nullptr;
    parent = nullptr;
    end = nullptr;
}

Node::Node(AminoCodon acPair) {
    gene.startIndex = 0;
    gene.endIndex = 0;
    gen = false;
    this->acPair = acPair;
    next = nullptr;
    parent = nullptr;
    end = nullptr;
}