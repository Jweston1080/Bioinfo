//
//  deBruijnByHash.cpp
//  k-assembler
//
//  Created by Joe Song on 11/24/15.
//  Copyright © 2015 Joe Song. All rights reserved.
//
//  Updated 3/19/2018

#include "k-assembler.hpp"

#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>

using namespace std;

struct DNAHasher
// Hash function used for DNA sequence
{
    std::size_t operator()(const string & seq) const
    {
        size_t val = 0;
        // TO DO: Write a DNA sequence hash function here
        // BEGIN your code here:
        for (const char &c : seq) {
            val <<= 2;
            switch (c) {
                case 'A': case 'a': val += 0; break;
                case 'C': case 'c': val += 1; break;
                case 'G': case 'g': val += 2; break;
                case 'T': case 't': val += 3; break;
                default: throw std::invalid_argument("Invalid DNA base");
            }
        // END your code above
        }
        return val;
    }
};


struct AlphabetHasher
// An example hash function used for the English alphabet
{
    std::size_t operator()(const string & seq) const
    {
        size_t val = 0;
        size_t max_width=20;
        for(size_t i=0; i<seq.size() && i<max_width; ++i) {
            val = val << 5;
            val += tolower(seq[i])-'a';
        }
        return val;
    }
};

// define the hash table class
// typedef unordered_multimap<string, size_t, AlphabetHasher> CSeqHash;

typedef unordered_multimap<string, size_t, DNAHasher> CSeqHash;

CSeqHash create_hash_table(const vector<string> & kmers)
// create one hash table by inserting both the prefix and suffix of each
//   k-mer. The prefix and suffix is the key. Associated with each element
//   in the hash table is the node id for that prefix or suffix in the
//   de Bruijn graph to be constructed.
{
    CSeqHash ht;
    size_t node_id=0; // the node id will be used in the de Bruijn graph
    for (auto i=0u; i<kmers.size(); ++i) {
        for(auto j=0u; j<2; ++j) { // j=0: prefix; j=1: suffix
            auto key = kmers[i].substr(j, kmers[i].length()-1);
            if (ht.find(key) == ht.end()) {
                ht.insert(make_pair(key, node_id ++));
            }
        }
    }
    return ht;
}

void create_deBruijn_graph_by_hashing(const vector<string> & kmers, DiGraph & g)
// create a deBruijn graph by inserting all k-mers into the graph by hashing
{
    // TO DO:
    
    // BEGIN your code below:
    
    // create one hash table for both the k-1 prefix and suffix of
    //   each k-mer
    CSeqHash ht = create_hash_table(kmers);
    g.m_nodes.resize(ht.size());
    // initialize an empty node vector for graph g
    // for each k-mer
    for (const auto &kmer : kmers)
    {
        // find the prefix node id from_id from the hash table
        string prefix = kmer.substr(0, kmer.length() - 1);
        string suffix = kmer.substr(1, kmer.length() - 1);
        // update node from_id's label to prefix if necessary
        // find the suffix node id to_id from the hash table
        // update node to_id's label to suffix if necessary
        size_t from_id = ht.find(prefix)->second;
        size_t to_id = ht.find(suffix)->second;
        // create a new edge (from_id, to_id) by inserting node
        //   to_id into the adjaceny list of node from_id
        if (g.m_nodes[from_id].m_label.empty()) {
            g.m_nodes[from_id].m_label = prefix;
        }
        if (g.m_nodes[to_id].m_label.empty()) {
            g.m_nodes[to_id].m_label = suffix;
        }
        // update the number of incoming edges of node to_id
        g.m_nodes[from_id].m_outgoing.push_back(to_id);
        // transfer the nodes from the vector to the graph
        g.m_nodes[to_id].m_num_of_incoming++;    

   } // end for loop

    // END your code above

}
