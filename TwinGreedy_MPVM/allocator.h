#ifndef ALLOCATOR_H
#define ALLOCATOR_H


#include <ctime>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include "anyoption.h"
#include "utils.h"
#include "advertiser.h"
#include "TimGraph.h"
#include <cmath>
#include <algorithm>
#include <list>
#include <unordered_set>
#include "random"
#include "time.h"

namespace _Cide{
    
    typedef std::vector<advertiser*> advertiserList;//把advertiserList定义为vector<advertiser*>数组的别名
    typedef std::vector<TimGraph*> TimGraphList;

    class allocator {
        
    public:
        
        time_t startTime;

        float alpha;
        AnyOption *opt;//option对象
        int n, m, nrTopics, nrCompanies;

        string delim, probGraphFile;
        float epsilon;
        float theta_0;
        int b = -1;


        advertiserList *advList;
        TimGraphList *timList1,*timList2,*timListTemp1,*timListTemp2;
        TimGraph *timTemp,*timTempForCalVal;

        void TwinGreedyFastFinal();
        void ResidualRandom();
        void TwinGreedy();
        void SampleGreedyFinal();
        void greedy(list<pair<int,int>> M);
        void fantom();
        float pairgreedy(TimGraphList *timlist);

        bool first_u1;
        bool first_u2;

        void copy_tim(TimGraph *src,TimGraph *dest);
        void copy_tim_fantom(TimGraph *src, TimGraph *dest);

        TimGraphList *timList;

        allocator(AnyOption* opt);
        ~allocator();

        int windowSize; 

        void readTICGraph();
        void readItemDistsFile();
        void readIncentiveCosts();

        float USM(vector<pair<int, int>> S);

        double cal_f_value_for_usm(int H_num, int u);
    };

}


#endif
