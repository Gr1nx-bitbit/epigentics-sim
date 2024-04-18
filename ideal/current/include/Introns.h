#ifndef INTRONS_H
#define INTRONS_H
#include <iostream>
#include <vector>
#include "Node.h"
#include "RnaAndExcision.h"
#include "AminoCodon.h"
#include "CodonTree.h"

Node* examine(std::string sequence, int length);
rnaAndExcision matureRNA(std::string sequence, int length, Node* head);
void compare(std::string premature, int length, rnaAndExcision preAmino);
std::string dnaTOrna(std::string dna, int length);
AminoCodon parse(std::string input);
std::vector<Node*> peptideSynthesis(std::string matureRna, int length, CodonTree* head);
void displayAbbreviations(Node* root);

#endif