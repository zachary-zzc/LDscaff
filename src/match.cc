/* Author: zicheng*/

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <getopt.h>
#include <cstdlib>
#include <string.h>
#include <time.h>

#include "graph.h"

#define INF 65535

using namespace std;

const static char* myopts = "s:n:c:m:l:t:r:h";

static char* sample_list = NULL;
static char* log_file = (char*)"edge.log";
static char* tmp_folder = (char*)"/tmp";
static int nsample = 0;
static int nscaf = 0;
static int nmarker = 100;
static float minor = 0.2;

static void Usage() {
    fprintf(stderr, 
            "Program: matching (Link scaffolds by LD)\n"
            "Version: 1.0\n"
            "Contact: Zicheng ZHAO <zachazhao2-c@my.cityu.edu.hk>\n"
            "\n"
            "Usage: matching -s <sample.list> -n <nsample> -c <nscaf> [options]\n"
            "\n"
            "Required options:\n"
            "         -s FILE   Sample list file. This file should contain one filename per line.\n"
            "                   Each file in the list should contain genotype data with one scaffold\n"
            "                   per line and space-separated marker values.\n"
            "                   Example format of a sample file:\n"
            "                   0 1 2 0 2 1 ... (for scaffold 1)\n"
            "                   1 0 1 2 0 0 ... (for scaffold 2)\n"
            "                   ...\n"
            "\n"
            "         -n INT    Number of samples in the sample list. This should match the\n"
            "                   number of files listed in the sample list file.\n"
            "\n"
            "         -c INT    Number of scaffolds in each sample file. Each sample file\n"
            "                   should have the same number of scaffolds, with one scaffold per line.\n"
            "\n"
            "Other options:\n"
            "         -t STR    Temp directory [default: /tmp]\n"
            "         -m INT    Number of markers to use from each scaffold [default: 100]\n"
            "                   The program will use up to this many markers from the beginning\n"
            "                   of each scaffold line.\n"
            "         -l STR    Log file for edge weights [default: edge.log]\n"
            "         -r FLOAT  Minor allele cutoff [default: 0.2]\n"
            "         -h        Print this help message\n");
    return;
}

static int parse_opt(int argc, char** argv) {
    int ch;
    
    while ((ch = getopt(argc, argv, myopts)) != -1) {
        switch(ch) {
            case 's':
                sample_list = strdup(optarg);
                break;
            case 'n':
                nsample = atoi(optarg);
                break;
            case 'c':
                nscaf = atoi(optarg);
                break;
            case 't':
                tmp_folder = strdup(optarg);
                break;
            case 'm':
                nmarker = atoi(optarg);
                break;
            case 'l':
                log_file = strdup(optarg);
                break;
            case 'r':
                minor = atof(optarg);
                break;
            case 'h':
                Usage();
                return 1;
            case '?':
                fprintf(stderr, "Unrecognized option - %c\n", optopt);
                return 1;
            default:
                return 1;
        }
    }
    
    if (sample_list == NULL) {
        fprintf(stderr, "Error: Need sample list file (-s)\n");
        Usage();
        return 1;
    }
    
    if (nsample <= 0) {
        fprintf(stderr, "Error: Need number of samples (-n)\n");
        Usage();
        return 1;
    }
    
    if (nscaf <= 0) {
        fprintf(stderr, "Error: Need number of scaffolds (-c)\n");
        Usage();
        return 1;
    }
    
    return 0;
}

int main(int argc, char** argv) {
    if (parse_opt(argc, argv) != 0) {
        return 1;
    }

    // Verify sample list file exists
    ifstream test_file(sample_list);
    if (!test_file.is_open()) {
        fprintf(stderr, "Error: Cannot open sample list file: %s\n", sample_list);
        return 1;
    }
    test_file.close();

    clock_t start, finish;
    double totalTime;
    start = clock();

    cout << "=====================================================" << endl;
    cout << "Starting matching with parameters:" << endl;
    cout << "Sample list file:     " << sample_list << endl;
    cout << "Number of samples:    " << nsample << endl;
    cout << "Number of scaffolds:  " << nscaf << endl;
    cout << "Number of markers:    " << nmarker << endl;
    cout << "Log file:             " << log_file << endl;
    cout << "Temp folder:          " << tmp_folder << endl;
    cout << "Minor allele cutoff:  " << minor << endl;
    cout << "=====================================================" << endl;

    // Call the processing function from graph.cc
    int result = processMatching(sample_list, nsample, nscaf, nmarker, log_file);

    finish = clock();
    totalTime = (double)(finish - start) / CLOCKS_PER_SEC;
    cout << "=====================================================" << endl;
    cout << "[main] Finish, total time: " << totalTime << " seconds" << endl;
    cout << "=====================================================" << endl;

    return result;
}
