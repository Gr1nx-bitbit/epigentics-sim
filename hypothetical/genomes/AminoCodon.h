#ifndef AMINOCODON_H
#define AMINOCODON_H

#include <string>

struct AminoCodon {
    std::string aminoAcid;
    std::string codon;
    bool startCodon;
    bool terminationCodon;
};

#endif