#ifndef TWINGREEDY_H
#define TWINGREEDY_H
#include "Graph.h"
#include <math.h>
#include <unordered_map>

class TwinGreedy{
    Graph* G_ptr;
    int k;
    int h;
    NodeSet s[2];
    bool search[2];
    NodeSet M[2];
    vector<int> cnt[2]; 
    unordered_map<int, int> memo[2];
    bool update[2];
    int argmax(pair<double, int>& res);
public:
    TwinGreedy(Graph* g_ptr, int k){
        this->k = k;
        G_ptr = g_ptr;
        h = G_ptr->h;
        search[0] = search[1] = true;
        query_time = 0;
        cnt[0] = vector<int>(h, 0);
        cnt[1] = vector<int>(h, 0);
        M[0] = g_ptr->nodes;
        M[1] = g_ptr->nodes;
    }
    bool check(Solution& s);
    Solution alg();
    int query_time;
};

bool TwinGreedy::check(Solution& s){
    double f_val = s.first;
    NodeSet V = s.second;
    double f = 0.0;
    for(auto u : V){
        for(auto e : G_ptr->adj[u]){
            auto v = e.first;
            auto w = e.second;
            if(!V.count(v)) f += w;
        }
    }
    printf("check value = %f\n", f);
    printf("solution value = %f\n", f_val);
    return f == f_val;
}

int TwinGreedy::argmax(pair<double, int>& res){
    pair<double, int> item1, item2;
    if(search[0]){
        if(!update[0]){
            double max_margin = DBL_MIN;
            int max_node = 0;
            for(auto& item : memo[0]){
                if(cnt[0][G_ptr->kind[item.first]] == k) continue;
                if(max_margin < item.second){
                    max_margin = item.second;
                    max_node = item.first;
                }
            }
            item1 = {max_margin, max_node};
        }
        else {
            double max_margin = DBL_MIN;
            int max_node = 0;
            for(auto u : M[0]){
                if(cnt[0][G_ptr->kind[u]] == k){
                    continue;
                }
                ++query_time;
                double margin = 0.0;
                for(auto e : G_ptr->adj[u]){
                    auto v = e.first;
                    auto w = e.second;
                    if(s[0].count(v)) margin -= w;
                    else margin += w;
                }
                memo[0][u] = margin;
                if(margin > max_margin){
                    max_margin = margin;
                    max_node = u;
                }
            }
            item1 = {max_margin, max_node};
        }
    }
    if(search[1]){
        if(!update[1]){
            double max_margin = DBL_MIN;
            int max_node = 0;
            for(auto& item : memo[1]){
                if(cnt[1][G_ptr->kind[item.first]] == k) continue;
                if(max_margin < item.second){
                    max_margin = item.second;
                    max_node = item.first;
                }
            }
            item2 = {max_margin, max_node};
        }
        else {
            double max_margin = DBL_MIN;
            int max_node = 0;
            for(auto u : M[1]){
                if(cnt[1][G_ptr->kind[u]] == k){
                    continue;
                }
                ++query_time;
                double margin = 0.0;
                for(auto e : G_ptr->adj[u]){
                    auto v = e.first;
                    auto w = e.second;
                    if(s[1].count(v)) margin -= w;
                    else margin += w;
                }
                memo[1][u] = margin;
                if(margin > max_margin){
                    max_margin = margin;
                    max_node = u;
                }
            }
            item2 = {max_margin, max_node};
        }
    }
    if(search[0] && search[1]){
        if(item1.first >= item2.first){
            res = item1; return 0;
        }
        else {
            res = item2; return 1;
        }
    }
    if(search[0]){
        res = item1; return 0;
    }
    if(search[1]) {
        res = item2; return 1;
    }
}

Solution TwinGreedy::alg(){
    clock_t start, end;
    start = clock();
    update[0] = update[1] = true;
    int insert_cnt = 0;
    while(search[0] || search[1]){
        ++insert_cnt;
        pair<double, int> res;
        int i = argmax(res);
        double margin = res.first;
        int u = res.second;
        if(margin <= 0) {
            //printf("margin < 0\n");
            break;
        }
        if(s[i].count(u)){
            //printf("repeat node\n");
            break;
        }
        s[i].insert(u);
        M[0].erase(u);
        M[1].erase(u);
        update[i] = true;
        update[(i+1)%2] = false;
        memo[0].erase(u);
        memo[1].erase(u);
        int t = G_ptr->kind[u];
        cnt[i][t]++;
        if(s[i].size() == h * k) search[i] = false;
    }
    double f_val[2];
    f_val[0] = f_val[1] = 0;
    for(int i = 0; i < 2; i++){
        for(auto u : s[i]){
            for(auto e : G_ptr->adj[u]){
                auto v = e.first;
                auto w = e.second;
                if(!s[i].count(v)) f_val[i]+=w;
            }
        }
    }
    Solution res;
    if(f_val[0] >= f_val[1]) res = {f_val[0], s[0]};
    else res = {f_val[1], s[1]};
    end = clock();
    printf("twin greedy time = %fs\n", double(end-start)/CLOCKS_PER_SEC);
    printf("twin query time = %d\n", query_time);
    return res;
}
#endif