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
      // DNA sequence hash function implementation
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
        for (size_t i = 0; i < seq.size(); ++i) {
            val = (val << 2) | (seq[i] & 3); // combine two bits of each nucleotide to a byte
        
        // END your code above
        }
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

void create_deBruijn_graph_by_hashing(const vector<string> &kmers, DiGraph &g)
{
    CSeqHash node_map; // hash table for prefix and suffix of each k-mer
    vector<Node> nodes; // node vector for graph g
    std::cout << "Entered create_deBruijn_graph_by_hashing()";

    // iterate over each k-mer and create corresponding nodes and edges
    for (const auto &kmer : kmers)
    {
        string prefix = kmer.substr(0, kmer.length() - 1);
        string suffix = kmer.substr(1);

        // check if prefix node exists
        auto range = node_map.equal_range(prefix);
        bool prefix_node_exists = false;
        for (auto it = range.first; it != range.second; ++it)
        {
            size_t node_id = it->second;
            if (nodes[node_id].m_label == prefix)
            {
                prefix_node_exists = true;
                break;
            }
        }

        size_t from_id;
        if (prefix_node_exists)
        {
            from_id = range.first->second;
        }
        else
        {
            from_id = nodes.size();
            nodes.emplace_back(Node{prefix});
            node_map.emplace(prefix, from_id);
        }

        // check if suffix node exists
        range = node_map.equal_range(suffix);
        bool suffix_node_exists = false;
        size_t to_id;
        for (auto it = range.first; it != range.second; ++it)
        {
            size_t node_id = it->second;
            if (nodes[node_id].m_label == suffix)
            {
                suffix_node_exists = true;
                to_id = node_id;
                break;
            }
        }

        if (!suffix_node_exists)
        {
            to_id = nodes.size();
            nodes.emplace_back(Node{suffix});
            node_map.emplace(suffix, to_id);
        }

        // add edge
        nodes[from_id].m_outgoing.push_back(to_id);
        ++nodes[to_id].m_num_of_incoming;
    }

    // assign nodes to graph
    g.m_nodes = std::move(nodes);
    std::cout << "Entered create_deBruijn_graph_by_hashing()";

}

