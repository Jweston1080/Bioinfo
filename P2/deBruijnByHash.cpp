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
//task 1 I NEED TO DO
{
    std::size_t operator()(const string & seq) const
    {
        size_t val = 0;
        // TO DO: Write a DNA sequence hash function here
        // BEGIN your code here:
        
        // get length of kmer
        size_t k = seq.size();
        /**
         * We are doing this to handle odd and even kmers,
         * even obviously works since we are diviing by two.
         * Since int rounds down, diving by 2 then subtracting result
         * from k will work for odd kmers as well
        */

        // Get the prefix and suffix of the k-mer
        size_t prefix_len = k / 2;   
        size_t suffix_len = k - prefix_len;
        string prefix = seq.substr(0, prefix_len);
        string suffix = seq.substr(prefix_len);
        
        // Concatenate the prefix and suffix and hash the resulting string
        string prefix_suffix = prefix + suffix;
        for (size_t i = 0; i < prefix_suffix.size(); ++i) {
            val = (val << 2) | (prefix_suffix[i] & 3); // combine two bits of each nucleotide to a byte
        }
        
        // END your code above
        return val;
    }
};


struct AlphabetHasher
// An example hash function used for the English alphabet
//task 1 helpful example
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
// related to my task 1 but nothing needs to be done here, DNAHasher is
// implicitly called
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
    /*
    // BEGIN your code below:
    
    // create one hash table for both the k-1 prefix and suffix of
    //   each k-mer
    CSeqHash ht = create_hash_table(kmers); //My task 1 is required for this
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
    */


    // create one hash table for both the k-1 prefix and suffix of each k-mer
    CSeqHash ht = create_hash_table(kmers);
    // initialize an empty node vector for graph g
    vector<Node> nodes(ht.size());
    for (auto itr = ht.begin(); itr != ht.end(); ++itr) {
        nodes[itr->second].m_label = itr->first;
    }
    // for each k-mer
    for (auto i = 0u; i < kmers.size(); ++i) {
        // find the prefix node id from_id from the hash table
        auto prefix_key = kmers[i].substr(0, kmers[i].length() - 1);
        auto itr = ht.find(prefix_key);
        auto from_id = itr->second;
        // find the suffix node id to_id from the hash table
        auto suffix_key = kmers[i].substr(1, kmers[i].length() - 1);
        itr = ht.find(suffix_key);
        auto to_id = itr->second;
        // create a new edge (from_id, to_id) by inserting node to_id into the adjacency list of node from_id
        nodes[from_id].m_outgoing.push_back(to_id);
        // update the number of incoming edges of node to_id
        nodes[to_id].m_num_of_incoming++;
    }
    // transfer the nodes from the vector to the graph
    g.m_nodes.swap(nodes);


   } // end for loop

    // END your code above


