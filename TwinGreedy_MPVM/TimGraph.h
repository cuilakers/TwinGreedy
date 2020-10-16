#ifndef _TIMGRAPH_H_
#define _TIMGRAPH_H_

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <set>
#include <algorithm>
#include <vector>
#include <deque>
#include <utility>
#include "sfmt/SFMT.h"
#include "advertiser.h"
#include "utils.h"
#include <unordered_set>
#include <map>

#define IF_TRACE(args) ;

namespace _Cide {
	
	
	class TimGraph {
		
	public:

        std::vector< std::vector<int> > hyperG_temp;//保存每个结点出现在了哪些超边中，便于之后算结点的覆盖
        std::vector< std::vector<int> > hyperGT_temp;//保存每个超边所涉及的结点
        std::set<int> seedSetTemp; // temporary seed set created for Kpt estimation purposes
        std::vector<int> nodeaval;

		sfmt_t sfmtSeed;

		float epsilon;
        float lambda;
		int n, m;

        TimGraph(advertiser *adv,  float eps, int nrNodes, int nrEdges);
        ~TimGraph(void);

        void RRset_generation();

        void opim_assign_best_node(int node);
        void opim_reverse_assign_best_node(int node,int degree);

        advertiser *adv;//每个广告商对应一个图，该图最初是一个完整的社交图


        //使用前需要拷贝
        std::vector< std::vector<int> > hyperG_adv;//保存每个结点出现在了哪些超边中，便于之后算结点的覆盖
        std::vector< std::vector<int> > hyperGT_adv;//保存每个超边所涉及的结点，即超边（RR集）本身
        std::vector<int> hyper_degree; // 0 to n, her node icin var - bu RR set yaratilirken doldurulmali, hyperG_adv'den dolduruluyor (hyperG_adv[i].size())
        std::vector<bool> isCovered; // her RR_id icin - bu da RR set yapilirken de

        vector<bool> node_aval;

        int64 theta;//RR集的大小
        int64 theta_old;

        //使用前需要清零
        multimap<float,int> criterQueue; // the float keys differ wrt cost-sensitive or not
        set<int> usersExamined;//被影响的结点
        std::vector<int> seedSet; // 每个广告商的种子集
        std::vector<float> seedMgRevs;
        float currentPayment;
        float currentRevenue;
        float currentSeedCosts;


        //建立超边时使用，之后不再使用，因为无需拷贝无需清零
        std::vector<bool> visit;
        std::vector<int> visit_mark;
        std::deque<int> q;//双向队列

		int BuildHypergraphNode(int uStart, int hyperiiid, bool addHyperEdge, std::vector< std::vector<int> > &hyperGT);


        void opim_help_cal_f(int node);
    };
	
}


#endif
