//
//  EulerPath.cpp
//  k-assembler
//
//  Created by Joe Song on 11/23/15.
//  Copyright © 2015-2021 Joe Song. All rights reserved.
//
//  Updated 3/19/2018
//  Updated 9/19/2021

#include "k-assembler.hpp"
#include <stack>
using namespace std;

size_t source(const DiGraph & g)
// Find a source node from g: the node has one more
// outgoing edge than incoming edge. When such a node
// does not exist, a value of the total number of nodes
// is returned.
{
    size_t i;
    
    for (i=0; i<g.m_nodes.size(); ++i) {
        if (g.m_nodes[i].m_outgoing.size()
            == g.m_nodes[i].m_num_of_incoming + 1) {
            break;
        }
    }
    return i;
}

size_t sink(const DiGraph & g)
// Find a source node from g: the node has one more
// outgoing edge than incoming edge. When such a node
// does not exist, a value of the total number of nodes
// is returned.
{
    size_t i;
    for (i=0; i<g.m_nodes.size(); ++i) {
        if (g.m_nodes[i].m_outgoing.size() + 1
            == g.m_nodes[i].m_num_of_incoming) {
            break;
        }
    }
    return i;
}


list<size_t> find_Eulerian_cycle(DiGraph & g)
{
    list<size_t> cycle; // main cycle
    size_t src = source(g);    // find the source node
    size_t dest = sink(g);     // find the sink node

    // Check if the graph has an Eulerian cycle
    if (!has_Eulerian_path(g)) {
        throw "The graph does not have an Eulerian cycle!";
    }

    stack<size_t> nodeStack;

    // Start from an arbitrary vertex
    nodeStack.push(0);

    while (!nodeStack.empty()) {
        size_t curNode = nodeStack.top();

        if (!g.m_nodes[curNode].m_outgoing.empty()) {

            // Traverse the outgoing edges of the current node
            size_t neighbor = g.m_nodes[curNode].m_outgoing.front(); // Only process one adjacent node at a time
            g.m_nodes[curNode].m_outgoing.pop_front(); // Delete that node from the adjacency list
            nodeStack.push(neighbor);
        }
        else {
            // If we have visited this node before, it means we have completed a cycle
            // Add this node to the cycle and remove it from the stack
            cycle.push_back(curNode);
            nodeStack.pop();
        }
    }

    // Reverse the order of nodes in the cycle to get the correct orientation
    cycle.reverse();

    return cycle;
}

list<size_t> find_Eulerian_path(DiGraph & g)
// find an Eulerian path from graph g, assuming g has such a path
{
    list <size_t> path, cycle;
    
    size_t src = source(g);    // find the source node
    size_t dest = sink(g);     // find the sink node
    
    // In the special case graph with only cycles
    //   and no source nor sink, we choose node 0 as the
    //   start and end of the Eulerian path:
    src = src >= g.m_nodes.size() ? 0 : src;
    dest = dest >= g.m_nodes.size() ? 0 : dest;
    
    vector<Node> & nodes = g.m_nodes;
    
    // add an edge from the sink node to the source node
    nodes[dest].m_outgoing.push_back(src);
    
    // increase the incoming degree of the source node by one
    nodes[src].m_num_of_incoming ++;
    
    cycle = find_Eulerian_cycle(g);
    
    list<size_t>::iterator pos_src, pos_dest;
    
    for (pos_dest = cycle.begin(); pos_dest != cycle.end(); ++ pos_dest) {
        pos_src = pos_dest;
        pos_src ++;
        if(pos_src != cycle.end()) {
            if (*pos_src == src && *pos_dest == dest) {
                break;
            }
        } else {
            break;
        }
    }
    
    if (pos_src != cycle.end() && pos_dest != cycle.end()) {
        
        /*
        // remove the last element on the cycle which is the same
        //   with the first element on the cycle:
        cycle.pop_back();
        
        // Formulate a path from the cycle such that src is
        //   the first and dest is the last on the path:
        path.splice(path.end(), cycle, pos_src, cycle.end());
        path.splice(path.end(), cycle, cycle.begin(), ++ pos_dest);
         */
        
        auto pos = pos_src;
        do {
            path.push_back(*pos);
            pos ++;
            if (pos == cycle.end()) {
                pos = cycle.begin();
                pos ++; // skip the first node in the cycle
                        // which is the same with the last
                        // last node in the cycle.
            }
        } while (pos != pos_src);
        
    } else {
        throw "Searching for Eulerian path has failed!";
    }
    
    // return the path
    return path;
}

bool has_Eulerian_path(const DiGraph & g)
// determine if graph g has an Eulerian path. This path could
//   be a cycle in special cases.
{
    bool exist = true;
    
    size_t numSources=0;
    size_t numSinks=0;
    
    for (auto & node : g.m_nodes) {
        size_t out = node.m_outgoing.size();
        size_t in = node.m_num_of_incoming;
        if(out == in) { // check for intermediate balanced node
            continue;
        } else if(out == in + 1) { // check for source node
            numSources ++;
            if (numSources > 1) {
                exist = false;
                break;
            }
        } else if(out + 1 == in) { // check for sink node
            numSinks ++;
            if (numSinks > 1) {
                exist = false;
                break;
            }
        } else {
            exist = false;
            break;
        }
    }
    
    return exist;
}


