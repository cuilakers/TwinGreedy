#ifndef COMPETITOR_H
#define COMPETITOR_H
#include "Graph.h"
#include <vector>
#include <queue>
#include <unordered_map>
#include <math.h>
#include <algorithm>
using namespace std;

class Competitor{
    Graph* G_ptr;
    Solution USM(NodeSet M);
    Solution random_USM(NodeSet M);
    Solution greedy(NodeSet M, int k);
    int h;
    vector<int> cnt;
    vector<int> cnt_M;
    static bool cmp(const pair<double, int>& a, const pair<double, int>&b){
        return a.first > b.first;
    }
public:
    Competitor(Graph* g_ptr){
        G_ptr = g_ptr;
        h = G_ptr->h;
    }
    Solution SampleGreedy(int k);
    Solution fantom(int k);
    Solution random_fantom(int k);
    Solution RRG(int k);
    Solution random(int k);
    bool check(Solution& s);
    int query_time;
};

bool Competitor::check(Solution& s){
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

Solution Competitor::greedy(NodeSet M, int k){
    cnt = vector<int>(h, 0);
    NodeSet S;
    int res = h;
    while(res > 0 && !M.empty()){
        double max_margin = 0.0;
        int max_node = -1;
        for(auto u : M){
            if(cnt[G_ptr->kind[u]] == k) {
                continue;
            }
            query_time ++;
            double margin = 0.0;
            for(auto e : G_ptr->adj[u]){
                int v = e.first;
                double w = e.second;
                if(!S.count(v)) margin += w;
                else margin -= w;
            }
            //printf("node %d margin = %f\n", u, margin);
            if(margin > max_margin){
                max_margin = margin;
                max_node = u;
            }
        }
        if(max_node == -1){
            //printf("max margin = %f\n", max_margin);
            break;
        }
        S.insert(max_node);
        M.erase(max_node);
        int t = G_ptr->kind[max_node];
        cnt[t]++;
        if(cnt[t] == k) {
            res--;
            NodeSet M_aux = M;
            for(auto u : M_aux){
                if(G_ptr->kind[u] == t) M.erase(u); 
            }
        }
        for(int i = 0; i < h; i++){
            if(cnt[i] > k){
                printf("type %d cnt = %d\n", i, cnt[i]);
            }
        }
    }
    double f = 0.0;
    for(auto u : S){
        for(auto e : G_ptr->adj[u]){
            int v = e.first;
            double w = e.second;
            if(!S.count(v)) f += w;
        }
    }
    return {f, S};
}

Solution Competitor::USM(NodeSet M){
    NodeSet x;
    NodeSet y = M;
    double f_val = 0;
    for(auto u : M){
        double a = 0.0, b = 0.0;
        for(auto e : G_ptr->adj[u]){
            int v = e.first;
            double w = e.second;
            if(x.count(v)) a -= w;
            else if(M.count(v)) a += w;
            if(y.count(v)) b += w;
            else if(M.count(v)) b -= w;
        }
        if(a >= b) {
            x.insert(u); f_val += a;
        }
        else y.erase(u);
    }
    // Solution s{f_val, x};
    // if(check(s)) printf("USM correct!\n");
    if(x == M) printf("select all node\n");
    return {f_val, x};
}

Solution Competitor::random_USM(NodeSet M){
    const double eps = 1e-8;
    NodeSet x;
    NodeSet y = M;
    double f_val = 0;
    for(auto u : M){
        double a = 0.0, b = 0.0;
        for(auto e : G_ptr->adj[u]){
            int v = e.first;
            double w = e.second;
            if(x.count(v)) a -= w;
            else if(M.count(v)) a += w;
            if(y.count(v)) b += w;
            else if(M.count(v)) b -= w;
        }
        double ap = max(a, 0.0);
        double bp = max(b, 0.0);
        if(ap < eps && bp < eps){
            x.insert(u); f_val += a;
            //printf("a' and b' = 0\n");
        }
        else{
            double p1 = (double)rand()/(double)RAND_MAX;
            double p2 = (double)ap/(ap+bp);
            if(p2 < p1){
                x.insert(u); f_val += a;
            }
            else {
                y.erase(u);
            }
        }
    }
    // Solution s{f_val, x};
    // if(check(s)) printf("USM correct!\n");
    if(x == M) printf("select all node\n");
    return {f_val, x};
}

Solution Competitor::fantom(int k){
    query_time = 0;
    clock_t start, end;
    start = clock();
    vector<Solution> A;
    NodeSet N = G_ptr->nodes;
    for(int i = 0; i < 2; i++){
        auto s = greedy(N, k);
        //printf("fantom %d: ", i);
        auto s_p = USM(s.second);
        A.push_back(s); A.push_back(s_p);
        //printf("before %f, after %f\n", s.first, s_p.first);
        for(auto u : s.second){
            N.erase(u);
        }
    }
    end = clock();
    printf("fantom time = %fs\n", double(end-start)/CLOCKS_PER_SEC);
    printf("fantom query time = %d\n", query_time);
    Solution res = {0.0, {}};
    for(int i = 0; i < A.size(); i++){
        if(res.first < A[i].first) res = A[i];
    }
    return res;
}

Solution Competitor::random_fantom(int k){
    query_time = 0;
    clock_t start, end;
    start = clock();
    vector<Solution> A;
    NodeSet N = G_ptr->nodes;
    for(int i = 0; i < 2; i++){
        auto s = greedy(N, k);
        auto s_p = random_USM(s.second);
        A.push_back(s); A.push_back(s_p);
        for(auto u : s.second){
            N.erase(u);
        }
    }
    end = clock();
    printf("random fantom time = %fs\n", double(end-start)/CLOCKS_PER_SEC);
    printf("random fantom query time = %d\n", query_time);
    Solution res = {0.0, {}};
    for(int i = 0; i < A.size(); i++){
        if(res.first < A[i].first) res = A[i];
    }
    return res;
}

Solution Competitor::SampleGreedy(int k){
    query_time = 0;
    clock_t start, end;
    start = clock();
    NodeSet N = G_ptr->nodes;
    for(auto u : G_ptr->nodes){
        double p = (double)rand()/(double)RAND_MAX;
        if(p < 0.5) N.erase(u); 
    }
    auto res = greedy(N, k);
    end = clock();
    printf("samplegreedy time = %fs\n", (double)(end-start)/CLOCKS_PER_SEC);
    printf("samplegreedy query time = %d\n", query_time);
    return res;
}

Solution Competitor::random(int k){
    cnt = vector<int>(h, 0);
    NodeSet S;
    for(auto u : G_ptr->nodes){
        S.insert(u);
        if(S.size() == h * k) break;
    }
    double f_val = 0;
    for(auto u : S){
        for(auto e : G_ptr->adj[u]){
            int v = e.first;
            double w = e.second;
            if(!S.count(v)) f_val += w;
        }
    }
    return {f_val, S};
}

Solution Competitor::RRG(int k){
    query_time = 0;
    cnt = vector<int>(h, 0);
    cnt_M = vector<int>(h, 0);
    int h = G_ptr->h;
    int r = h * k;
    NodeSet S_i;
    NodeSet M_res = G_ptr->nodes;
    vector<pair<double, int>> Q;
    bool flag = true; // if true, calculate Q from scratch
    clock_t start, end;
    start = clock();
    for(int i = 0; i < r; i++){
        if(flag){
            Q.clear();
            for(auto v : M_res){
                int t = G_ptr->kind[v];
                if(cnt[t] < k){
                    double margin = 0.0;
                    for(auto e : G_ptr->adj[v]){
                        int u = e.first;
                        double w = e.second;
                        if(!S_i.count(u)) margin += w;
                        else margin -= w;
                    }
                    Q.push_back({margin, v});
                    query_time++;
                }
            }
            sort(Q.begin(), Q.end(), cmp);
        }
        vector<int> M_i;
        cnt_M = vector<int>(h, 0);
        int res = h;
        for(auto& item : Q){
            double margin = item.first;
            int u = item.second;
            int t = G_ptr->kind[u];
            if(margin <= 1e-10) break;
            if(cnt_M[t] < k){
                cnt_M[t]++;
                M_i.push_back(u);
                if(cnt_M[t] == k) res--;
                if(res == 0) break;
            }
        }
        while(M_i.size() < r - i) M_i.push_back(-1);
        //random_shuffle(M_i.begin(), M_i.end());
        //int u = M_i[0];
        int u = M_i[rand()%M_i.size()];
        if(u != -1) {
            S_i.insert(u);
            M_res.erase(u);
            cnt[G_ptr->kind[u]]++;
            flag = true;
        }
        else flag = false;
    }
    end = clock();
    printf("RRG time = %fs\n", double(end-start)/CLOCKS_PER_SEC);
    printf("RRG query time = %d\n", query_time);
    double f_val = 0.0;
    for(auto u : S_i){
        for(auto e : G_ptr->adj[u]){
            int v = e.first;
            double w = e.second;
            if(!S_i.count(v)) f_val += w;
        }
    }
    return {f_val, S_i};
}
#endif