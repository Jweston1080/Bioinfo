//
//  main.cpp
//  k-assembler
//
//  Created by Joe Song on 11/19/15.
//  Copyright © 2015 Joe Song. All rights reserved.
//
//  Updated 3/19/2018

#include <string>
using namespace std;

extern void test_seq_assembly();

extern void test_1(const string & method);

int main(int argc, const char * argv[])
{
    test_seq_assembly();

    return 0;
}
