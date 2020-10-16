#ifndef TWINGREEDYFAST_H
#define TWINGREEDYFAST_H
#include "Graph.h"
#include <math.h>
#include <unordered_map>
#include <math.h>

class TwinGreedyFast{
    Graph* G_ptr;
    int v_star;
    double f_v_star;
    vector<double> singleton_cut;
    void cal_singleton_cut();
    double epsilon;
    int k;
    int h;
    double tau_max, tau_min, tau;
    unordered_set<int> s[2];
    bool search[2];
    unordered_set<int> V; // candidate nodes
    double f_margin[2];
    vector<int> cnt[2]; 
    int argmax(int u);
public:
    TwinGreedyFast(Graph* g_ptr, double eps, int k){
        this->k = k;
        G_ptr = g_ptr;
        epsilon = eps;
        h = G_ptr->h;
        cal_singleton_cut();
        search[0] = search[1] = true;
        query_time = 0;
        cnt[0] = vector<int>(h, 0);
        cnt[1] = vector<int>(h, 0);
    }
    bool check(Solution& s);
    Solution alg();
    unsigned long long query_time;
};

bool TwinGreedyFast::check(Solution& s){
    double f_val = s.first;
    NodeSet V = s.second;
    double f = 0;
    for(auto u : V){
        for(auto e : G_ptr->adj[u]){
            int v = e.first;
            double w = e.second;
            if(!V.count(v)) f += w;
        }
    }
    printf("check value = %f\n", f);
    printf("solution value = %f\n", f_val);
    return f == f_val;
}

void TwinGreedyFast::cal_singleton_cut(){
    f_v_star = 0;
    for(auto u : G_ptr->nodes){
        double sum = 0.0;
        for(auto e : G_ptr->adj[u]){
            sum += e.second;
        }
        if(f_v_star < sum){
            f_v_star = sum;
            v_star = u;
        }
    }
    tau_max = f_v_star;
    //tau_min = epsilon * f_v_star / G_ptr->n;
    //tau_min = epsilon * f_v_star / G_ptr->n / h;
    tau_min = epsilon * f_v_star / k / h;
    tau = tau_max;
    V = G_ptr->nodes;
}

int TwinGreedyFast::argmax(int u){
    f_margin[0] = f_margin[1] = 0;
    if(search[0] && search[1]){
        for(auto e : G_ptr->adj[u]){
            int v = e.first;
            double w = e.second;
            if(s[0].count(v)) f_margin[0] -= w;
            else f_margin[0] += w;
            if(s[1].count(v)) f_margin[1] -= w;
            else f_margin[1] += w;
        }
        query_time += 2;
        if(f_margin[0] >= f_margin[1]) return 0;
        else return 1;
    }
    if(search[0]) {
        for(auto e : G_ptr->adj[u]){
            int v = e.first;
            double w = e.second;
            if(s[0].count(v)) f_margin[0] -= w;
            else f_margin[0] += w;
        }
        query_time += 1;
        return 0;
    }
    if(search[1]) {
        for(auto e : G_ptr->adj[u]){
            int v = e.first;
            double w = e.second;
            if(s[1].count(v)) f_margin[1] -= w;
            else f_margin[1] += w;
        }
        query_time += 1;
        return 1;
    }
}

Solution TwinGreedyFast::alg(){
    clock_t start, end;
    start = clock();
    double z;
    int for_cnt = 0;
    do{
        //printf("tau_max = %f\t", tau_max);
        //printf("tau_min = %f\t", tau_min);
        //printf("current tau = %f\t", tau);
        //printf("tau_max/tau_min = %f\t", tau_max/tau_min);
        //printf("s_0 size = %d\t", s[0].size());
        //printf("s_1 size = %d\t", s[1].size());
        //printf("query_time = %d\t", query_time);
        //printf("\n");
        for_cnt++;
        NodeSet M = V;
        //z = 0.0;
        while(!M.empty()){
            //int u = *next(M.begin(), rand() % M.size()); //too slow
            int u = *M.begin();
            int t = G_ptr->kind[u];
            
            if(cnt[0][t] == k && cnt[1][t] == k){
                M.erase(u); continue;
            }

            double delta[2] = {DBL_MIN, DBL_MIN};
            for(int i = 0; i < 2; i++){
                if(cnt[i][t] < k){
                    delta[i] = 0.0;
                    for(auto e : G_ptr->adj[u]){
                        int v = e.first;
                        double w = e.second;
                        if(s[i].count(v)) delta[i] -= w;
                        else delta[i] += w;
                    }
                    query_time++;
                }
            }
            //if(delta[0] == DBL_MIN && delta[1] == DBL_MIN) break;
            //if(delta[0] == DBL_MIN || delta[1] == DBL_MIN) M.erase(u);
            int i = delta[0] > delta[1] ? 0 : 1;
            if(delta[i] >= tau){
                s[i].insert(u);
                V.erase(u);
                cnt[i][t]++;
                //if(s[i].size() == h * k) search[i] = false;
            }
            //if(z < delta[i]) z = delta[i];
            M.erase(u);
        }
        tau /= (1+epsilon);
    } while(!V.empty() && tau >= tau_min/(1+epsilon) && (s[0].size() < k * h || s[1].size() < k * h));
    end = clock();
    //printf("tau_max = %f\t", tau_max);
    //printf("tau_min = %f\t", tau_min);
    //printf("current tau = %f\t", tau);
    //printf("tau_max/tau_min = %f\t", tau_max/tau_min);
    //printf("s_0 size = %d\t", s[0].size());
    //printf("s_1 size = %d\t", s[1].size());
    //printf("query_time = %d\t", query_time);
    //printf("\n");
    //printf("s_0 size = %d\n", s[0].size());
    //printf("s_1 size = %d\n", s[1].size());
    //printf("for count = %d\n", for_cnt);
    printf("twin-fast greedy time = %fs\n", double(end-start)/CLOCKS_PER_SEC);
    printf("twin-fast query time = %d\n", query_time);
    double f_val[2];
    f_val[0] = f_val[1] = 0.0;
    for(int i = 0; i < 2; i++){
        for(auto u : s[i]){
            for(auto e : G_ptr->adj[u]){
                int v = e.first;
                double w = e.second;
                if(!s[i].count(v)) f_val[i] += w;
            }
        }
    }
    if(f_val[0] >= f_val[1]) return {f_val[0], s[0]};
    else return {f_val[1], s[1]};
}
#endif