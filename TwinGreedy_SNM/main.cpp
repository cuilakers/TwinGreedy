#include "Competitor.h"
#include "TwinGreedyFast.h"
#include "Graph.h"
#include "TwinGreedy.h"


// To generate random graph with different nodes, please refer
// https://github.com/fahrbach/icml-2019-non-monotone-submodular-maximization
int main(int argc, char* argv[]){
    //srand(time(0));
    string fn = argv[1];
    int low = stoi(argv[2]);
    int high = stoi(argv[3]);
    int inc = stoi(argv[4]);
    int h = stoi(argv[5]);
    printf("h = %d\n", h);
    Graph* g_ptr = new Graph(fn, h);
    printf("graph %s\n", fn.c_str());
    double eps = 0.02;
    Competitor comp(g_ptr);
    for(int k = low; k <= high; k+=inc){
        printf("k = %d\n", k);
            
        TwinGreedyFast twin_fast(g_ptr, eps, k);
        auto s1 = twin_fast.alg();
        printf("twin-fast cut value = %d\n", (int)s1.first);
        printf("seed size = %d\n", s1.second.size());

        auto s2 = comp.fantom(k);
        printf("fantom cut value = %d\n", (int)s2.first);
        printf("seed size = %d\n", s2.second.size());

        auto s3 = comp.SampleGreedy(k);
        printf("samplegreedy cut value = %d\n", (int)s3.first);
        printf("seed size = %d\n", s3.second.size());

        auto s4 = comp.RRG(k);
        printf("RRG cut value = %d\n", (int)s4.first);
        printf("seed size = %d\n", s4.second.size());

        auto s5 = comp.random(k);
        printf("random cut value = %d\n", (int)s5.first);
        printf("seed size = %d\n", s5.second.size());

        auto s6 = comp.random_fantom(k);
        printf("random fantom cut value = %d\n", (int)s6.first);
        printf("seed size = %d\n", s6.second.size());

        TwinGreedy twin(g_ptr, k);
        auto s7 = twin.alg();
        printf("twin cut value = %d\n", (int)s7.first);
        printf("seed size = %d\n", s7.second.size());
    }
    delete g_ptr;
    g_ptr = NULL;
}