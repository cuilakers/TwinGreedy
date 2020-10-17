#ifndef GRAPH_H
#define GRAPH_H
#include <vector>
#include <cfloat>
#include <unordered_set>
#include <unordered_map>
#include "pair_hash.h"
using namespace std;
using NodeSet = unordered_set<int>;
using Solution = pair<double, NodeSet>;
using Edge = pair<int, double>;


class Graph{
public:
    int n;
    int h;
    unordered_set<pair<int, int>, pair_hash> inserted_edge;
    vector<vector<Edge>> adj;
    unordered_set<int> nodes;
    unordered_map<int, int> kind;
    Graph(string filename, int h);
};

Graph::Graph(string filename, int h){
    this->h = h;
    freopen(filename.c_str(), "r", stdin);
    scanf("%d", &n);
    printf("n = %d\n", n);
    adj.resize(n);
    int u, v;
    double w;
    int m = 0;
    while(scanf("%d %d %lf", &u, &v, &w) != EOF){
        if(inserted_edge.count({v, u})) continue;
        nodes.insert(u);
        nodes.insert(v);
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
        m += 2;
        //printf("edge %d %d weight %f\n", u, v, weight);
        inserted_edge.insert({u, v});
    }
    printf("nodes size = %d\n", nodes.size());
    printf("edge size = %d\n", m);
    for(auto u : nodes){
        kind[u] = rand()%h;
    }
}

#endif
