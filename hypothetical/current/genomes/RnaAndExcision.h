#ifndef RNAANDEXCISION_H
#define RNAANDEXCISION_H

#include <string>
#include "GenePosition.h"

struct rnaAndExcision {
    std::string* excised;
    GenePosition* indecies;
    std::string rna;
};

#endif