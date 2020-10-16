#ifndef _TIMGRAPH_C_
#define _TIMGRAPH_C_

#include "TimGraph.h"

namespace _Cide{
        
    TimGraph::TimGraph(advertiser *adv,  float eps, int nrNodes, int NrEdges)
    {

        this->epsilon = eps;
        this->adv = adv;//TIM图对应广告商
        this->n = nrNodes;
        this->m = NrEdges;


        this->seedSet.clear();//std::vector<int> seedSet; // keep like this to be able to update in this order
        this->usersExamined.clear();
        this->criterQueue.clear();
        this->seedMgRevs.clear();

        visit_mark = std::vector<int>(n,0);//访问过的结点号，size是n
        visit = std::vector<bool>(n,false);//访问过的标记，size是n
        hyper_degree = std::vector<int>(n,0);//超图的度
        //nodeCTRs = std::vector<float>(n,0);
        
        for(int i = 0; i < 12; i++)
            sfmt_init_gen_rand(&sfmtSeed, i+1234);

        currentPayment=0;
        currentRevenue=0;
        currentSeedCosts=0;
    }
    
    TimGraph::~TimGraph(void) {

    }
    void TimGraph::RRset_generation()
    {

        for(int i = 0; i < theta; i++)//构造初始采样超图，超边数满足阈值theta
        {
            hyperGT_adv.push_back(std::vector<int>());
            isCovered.push_back(false);
            BuildHypergraphNode(sfmt_genrand_uint32(&sfmtSeed)%n, i, true, hyperGT_adv);//随机起始结点之前已经生成好了，每个广告商在一个RR集上的起始结点是相同的
        }

        //hyperGT是采样算法sample生成的超图
        hyperG_adv.clear();//其实应该不需要，先留着吧
        for(int i = 0; i < n; i++)
            hyperG_adv.push_back(std::vector<int>());

        for(int i = 0; i < theta; i++)
        {
            for(int j = 0; j < (int) hyperGT_adv[i].size(); j++) //遍历超图中的每个超边的每个结点
            {
                int t = hyperGT_adv[i][j];
                hyperG_adv[t].push_back(i);//保存每个结点出现在了哪些超边中，便于之后算结点的覆盖
            }
        }
        for(int i = 0; i < n; i++)
        {
            hyper_degree[i] = hyperG_adv[i].size();//每个结点出现在不同超边中的次数就是结点的度
        }

        theta_old=theta;
    }


    void TimGraph::opim_assign_best_node(int node)
    {
        //is_selected[node] = false;//该结点已经被选，则不能再被选，置为false
        hyper_degree[node] = 0;//结点的度置为0

        seedSet.push_back(node);//选为该广告商的结点

        for(int j = 0; j < (int) hyperG_adv[node].size(); j++)//遍历被选结点出现过的超边
        { //for each RR set t covered by this node
            int t = hyperG_adv[node][j]; // rr set t，超边t
            if(!isCovered[t])//如果这条超边未被结点覆盖过
            { //if the RR set was not covered before - since we do not delete covered from hyperG_adv except new RR set genration
                isCovered[t]=true; // bestnode covers it,则现在标记为被覆盖了
                for(int z = 0; z < (int) hyperGT_adv[t].size(); z++)//遍历该超边涉及的结点
                { //decrease the degree of parents who are associated with this RR set
                    int item = hyperGT_adv[t][z];//获取结点标号
                    hyper_degree[item]--;//因为已选的种子结点已经覆盖了这些超边，因此在在要选择新的种子结点就不能把这些超边作为度的贡献，从而实现了边际增益的功能
                    //usersExamined.insert(item);//保存新被选的种子结点覆盖的超边会影响到的结点，这些结点的信息需要更新
                    if(hyper_degree[item] < 0)//已经已经小于0了则直接置为0
                    {
                        hyper_degree[item] = 0;
                    }
                }
            }
        }
    }
    void TimGraph::opim_help_cal_f(int node)
    {
        hyper_degree[node] = 0;//结点的度置为0
        for(int j = 0; j < (int) hyperG_adv[node].size(); j++)//遍历被选结点出现过的超边
        { //for each RR set t covered by this node
            int t = hyperG_adv[node][j]; // rr set t，超边t
            if(!isCovered[t])//如果这条超边未被结点覆盖过
            { //if the RR set was not covered before - since we do not delete covered from hyperG_adv except new RR set genration
                isCovered[t]=true; // bestnode covers it,则现在标记为被覆盖了
                for(int z = 0; z < (int) hyperGT_adv[t].size(); z++)//遍历该超边涉及的结点
                { //decrease the degree of parents who are associated with this RR set
                    int item = hyperGT_adv[t][z];//获取结点标号
                    hyper_degree[item]--;//因为已选的种子结点已经覆盖了这些超边，因此在在要选择新的种子结点就不能把这些超边作为度的贡献，从而实现了边际增益的功能

                    //usersExamined.insert(item);//保存新被选的种子结点覆盖的超边会影响到的结点，这些结点的信息需要更新

                    if(hyper_degree[item] < 0)//已经已经小于0了则直接置为0
                    {
                        hyper_degree[item] = 0;
                    }
                }
            }
        }
    }
    void TimGraph::opim_reverse_assign_best_node(int node,int degree)
    {
        is_selected[node] = true;//该结点已经被选，则不能再被选，置为false
        hyper_degree[node] = degree;//结点的度置为degree
        //seedSet.push_back(node);//选为该广告商的结点

        for(int j = 0; j < (int) hyperG_adv[node].size(); j++)//遍历被选结点出现过的超边
        { //for each RR set t covered by this node
            int t = hyperG_adv[node][j]; // rr set t，超边t
            if(isCovered[t])//如果这条超边被结点覆盖过
            { //if the RR set was not covered before - since we do not delete covered from hyperG_adv except new RR set genration
                isCovered[t]= false; // bestnode covers it,则现在标记为被覆盖了
                for(int z = 0; z < (int) hyperGT_adv[t].size(); z++)//遍历该超边涉及的结点
                { //decrease the degree of parents who are associated with this RR set
                    int item = hyperGT_adv[t][z];//获取结点标号
                    hyper_degree[item]++;//因为已选的种子结点已经覆盖了这些超边，因此在在要选择新的种子结点就不能把这些超边作为度的贡献，从而实现了边际增益的功能
                    //usersExamined.insert(item);//保存新被选的种子结点覆盖的超边会影响到的结点，这些结点的信息需要更新
                    if(hyper_degree[item] < 0)//已经已经小于0了则直接置为0
                    {
                        hyper_degree[item] = 0;
                    }
                }
            }
        }
    }
    int TimGraph::BuildHypergraphNode(int uStart, int hyperiiid, bool addHyperEdge, std::vector< std::vector<int> > &hyperGT)
    {
        int n_visit_edge=1;
        
        if(addHyperEdge)
        {
            // 			ASSERT((int)hyperGT.size() > hyperiiid);
            hyperGT[hyperiiid].push_back(uStart);
        }
        
        int n_visit_mark=0;
        
        q.clear();//双向队列
        q.push_back(uStart);
        // 		ASSERT(n_visit_mark < n);

        visit_mark[n_visit_mark++] = uStart;//图中uStart这个结点已访问，记录下来该标号并标记为true（已访问）
        visit[uStart] = true;
        //随机选择一个结点uStart作为采样算法的启动结点，加入到队列中
        while(!q.empty())
        {
            int expand = q.front();//队列不为空，说明还有可采样的结点，则取队列的首结点开始采样
            q.pop_front();
            
            int i = expand;
            for(int j = 0; j < (int)graphT[i].size(); j++)
            {
                //int u=expand;
                int v = graphT[i][j]; //parent of u in the original graph G，graphT存储了所有与结点u相关的边（v，u）的源结点v，也即graphT[i][j]
                n_visit_edge++;//indegree++，也即宽度++
                float randFloat= float(sfmt_genrand_uint32(&sfmtSeed))/ float(RAND_MAX)/2;
                if(randFloat > adv->probT[i][j])//产生随机数，大于边上的概率则continue，否则则访问该结点，应该是个IC模型
                    continue;
                if(visit[v])//如果已经访问过该结点则continue
                    continue;
                if(!visit[v])
                {
                    // 					ASSERT(n_visit_mark < n);
                    visit_mark[n_visit_mark++]=v;//标记为已访问
                    visit[v]=true;
                }
                q.push_back(v);
                
                if(addHyperEdge)
                {
                    // 					ASSERT((int)hyperGT.size() > hyperiiid);
                    hyperGT[hyperiiid].push_back(v);//超边中添加该结点
                }
            }
        }
        
        for(int i = 0; i < n_visit_mark; i++)
            visit[visit_mark[i]]=false;//重新置为false，便于下次使用
        
        return n_visit_edge;
        // returns number of edges considered, i.e., width of the RR set created from uStart，即为KPT算法返回宽度w（R）


    }

    
    
}

#endif

