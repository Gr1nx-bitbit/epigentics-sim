#ifndef AMINOCODON_H
#define AMINOCODON_H

#include <string>

struct AminoCodon {
    std::string aminoAcid;
    std::string codon;
    std::string abbreviation;
    bool startCodon;
    bool terminationCodon;
};

#endif