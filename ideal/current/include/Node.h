#ifndef NODE_H
#define NODE_H

#include "AminoCodon.h"
#include "GenePosition.h"


class Node {
    private:
        GenePosition gene;
        bool gen;
        AminoCodon acPair;
        Node* next;
        Node* parent;
        Node* end;
    public:
        Node();
        Node(GenePosition geen);
        Node(AminoCodon acPair);

        // ~Node() {
        //     if (end && parent && next) {
        //         Node* p = getParent();
        //         p->setNext(getNext());
        //         p->setEnd(getEnd());
        //         parent = nullptr;
        //         next = nullptr;
        //         end = nullptr;
        //         delete this;
        //     } else if (!parent && next && end) {
        //         Node* n = getNext();
        //         n->setEnd(getEnd());
        //         next = nullptr;
        //         end = nullptr;
        //         delete this;
        //     } else if (parent && next && !end) {
        //         Node* p = getParent();
        //         p->setNext(getNext());
        //         delete this;
        //     }
        // }

        GenePosition getGene(void) 
        { return gene; }

        void setGene(GenePosition geen) 
        { gene  = geen; gen = true; }

        bool isGene(void) \
        { return gen; }

        AminoCodon getAmino(void)
        { return acPair; }

        void setAmino(AminoCodon acPair) 
        { this->acPair = acPair; }

        Node* getNext(void) 
        { return next; }

        void setNext(Node* nxt) 
        { next = nxt; }

        Node* getParent(void) 
        { return parent; }

        void setParent(Node* parent) 
        { this->parent = parent; }

        Node* getEnd(void) 
        { return end; }

        void setEnd(Node* end) 
        { this->end = end; }
};

#endif