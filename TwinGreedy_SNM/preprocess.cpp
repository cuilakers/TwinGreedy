#include <stdio.h>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

int main(){
    int n;
    string filename = "old-er3k.el";
    freopen(filename.c_str(), "r", stdin);
    scanf("%d", &n);
    printf("n = %d\n", n);
    unordered_map<int, int> map_id;
    vector<vector<int>> adj(n);
    int u, v;
    int cnt = 0;
    while(scanf("%d %d", &u, &v) != EOF){
        if(!map_id.count(u)) map_id[u] = cnt++;
        if(!map_id.count(v)) map_id[v] = cnt++;
        adj[map_id[u]].push_back(map_id[v]);
    }
    printf("cnt = %d\n", cnt);
    FILE* fp;
    fp = fopen("old-er3k-post.el", "w+");
    fprintf(fp, "%d\n", n);
    for(int u = 0; u < n; u++){
        for(auto v : adj[u]){
            double w = (double)rand()/RAND_MAX;
            fprintf(fp, "%d %d %f\n", u, v, w);
        }
    }
    fclose(fp);
    return 0;
}