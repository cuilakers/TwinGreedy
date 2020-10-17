#include "allocator.h"
#include "anyoption.h"
#include "memoryusage.h"
namespace _Cide{

    allocator::allocator(AnyOption* opt1)
    {
        opt = opt1;
        delim = " \t";//分隔符是制表符TAB

        //读取算法的参数
        n = strToInt(opt->getValue("n"));//社交网络的结点数，即用户数
        m = strToInt(opt->getValue("m"));//社交网络的边数f
        nrTopics = strToInt(opt->getValue("nrTopics"));
        nrCompanies = strToInt(opt->getValue("nrCompanies"));
        alpha = strToFloat(opt->getValue("alpha"));

        epsilon = strToFloat(opt->getValue("epsilon"));
        theta_0 = strToFloat(opt->getValue("theta_0"));
        lambda = strToFloat(opt->getValue("lambda"));
        max_node=strToInt(opt->getValue("max_node"));

        for(int i = 0; i < n; i++) // keeps p^z_{uv} for each topic z
            graphT.push_back(vector<int>());//全局变量，给图中全部n个结点各一个vector数组，用于保存入边的源结点

        advList = new advertiserList();

        timList1 = new TimGraphList();
        timList2 = new TimGraphList();
        timListTemp1=new TimGraphList();
        timListTemp2=new TimGraphList();

        // create advertiser objects and initialize their TIC prob-vector probT, aligned with graphT
        for(int i = 0; i < nrCompanies; i++) //遍历所有广告商
        {
            advertiser *aa = new advertiser(i,nrTopics);//初始化广告商的对象

            advList->push_back(aa);//advList存储所有广告商的对象
            //构建每个广告商对应的timgraph
            _Cide::TimGraph *tim1 = new _Cide::TimGraph(aa, epsilon, n, m);
            _Cide::TimGraph *tim2 = new _Cide::TimGraph(aa, epsilon, n, m);
            _Cide::TimGraph *timTemp1 = new _Cide::TimGraph(aa, epsilon, n, m);
            _Cide::TimGraph *timTemp2 = new _Cide::TimGraph(aa, epsilon, n, m);


            timTempForCalVal = new _Cide::TimGraph(aa, epsilon, n, m);

            timList1->push_back(tim1);//timList存储所有广告商的timgraph
            timList2->push_back(tim2);
            timListTemp1->push_back(timTemp1);
            timListTemp2->push_back(timTemp2);

            aa->maxCost = 0;
            aa->minCost = (float) n ;
            aa->seedUserCosts.resize(n);//最多n个结点（用户），所以重设为n
            for (int i = 0; i < n; i++)
                aa->probT.push_back(std::vector< float>());

        }

        first_u1=true;
        first_u2=true;

        readItemDistsFile();
        readTICGraph();

        is_selected=std::vector<bool>(n,true);//初始化为true,表示该结点可选
        readIncentiveCosts(); // reading incentive costs as matrix nodes X ads
        B=B*nrCompanies;
        cout<<max_node<<endl;
        cout<<"h: "<<nrCompanies<<endl;

        TwinGreedyFastFinal();
        //TwinGreedy();
    }
    bool cmp(pair<pair<int, int>, float> a, pair<pair<int, int>, float> b)
    {
        return a.second > b.second;
    }
    
    void allocator::TwinGreedyFastFinal()
    {
        advertiser *adv;
        TimGraph *tim1,*tim2;
        sfmt_t sfmt_seed;
        //随机种子
        for(int i = 0; i < 12; i++)//随机种子
            sfmt_init_gen_rand(&sfmt_seed, i+1234);
        //sfmt_init_gen_rand(&sfmt_seed,time(NULL));

        //generate rrset,then copy
        tim1 = timList1->at(0);
        tim1->theta_old=theta_0;
        tim1->theta=theta_0;
        tim1->RRset_generation();

        TimGraph *tim=timList1->at(0);
        //S_1
        for(int i = 1; i < nrCompanies; i++)//h_1~h_n copy rrset from h_0
        {
            //generate 每个广告商的RR集
            tim1 = timList1->at(i);
            copy_tim(tim,tim1);
        }

        //S_2
        //*
        for(int i = 0; i < nrCompanies; i++)//h_0~h_n copy rrset from h_0
        {
            //generate 每个广告商的RR集
            tim2 = timList2->at(i);
            copy_tim(tim,tim2);
        }//*/

        float v_star=0.0;
        adv = advList->at(0);
        tim = timList1->at(0);

        for(int i=0;i<n;i++)
        {
            float temp=(( float) n * (( float) tim->hyper_degree[i] / tim->theta)) * adv->cpe-  lambda*adv->seedUserCosts[i]+lambda*B;
            if(temp>v_star) v_star=temp;
        }
        float tau_max=v_star;
        float tau_min=epsilon*v_star/(float)max_node;
       // float tau_min=epsilon*v_star/(float)n;
        //float tau_min=epsilon*v_star/((float)n*(float)nrCompanies);
        float tau=tau_max;
        cout<<"RRset Size: "<<theta_0<<endl;
        //cout<<"Budget: "<<B<<endl;

        //list<pair<int,int>> N;
        vector<pair<int,int>> M;
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<nrCompanies;j++)
            {
                //N.push_back(pair<int, int>(i,j));
                M.push_back(pair<int, int>(i,j));
            }
        }

        //int random_num;
        pair<int,int> u;
        float delta1,delta2;
        vector<bool> S1_selected(n, true);//初始化为false,表示该结点可选
        vector<bool> S2_selected(n, true);//初始化为false,表示该结点可选


        int now_index=0;

        vector<pair<int, int>> N_temp;
        float Z_max=-999999999;

        long long int mariginal_count=0;
        long long int martoid_count=0;
        int sum_node1=0;
        int sum_node2=0;
        float max_delta;
        bool S1full=false;
        bool S2full=false;
	time(&startTime);

        bool MorN=true;
        while(1)
        {
            if(MorN) {
                now_index = 0;
                N_temp.clear();
                Z_max = -999999999;
                while (now_index < M.size()) {
                    //get first element;
                    //*
                    u.first = M[now_index].first;
                    u.second = M[now_index].second;
                    // */

                    delta1 = delta2 = -1;
                    martoid_count++;
                    if (sum_node1 < max_node) {
                        if (S1_selected[u.first]) {
                            tim = timList1->at(u.second);
                            adv = advList->at(u.second);
                            delta1 = ((float) n * ((float) tim->hyper_degree[u.first] / tim->theta)) * adv->cpe -
                                     lambda * adv->seedUserCosts[u.first];
                            mariginal_count++;
                            if (first_u1) {
                                delta1 += lambda * B;
                                //first_u1= false;
                            }
                        }
                    } else {
                        S1full = true;
                    }
                    martoid_count++;
                    if (sum_node2 < max_node) {
                        if (S2_selected[u.first]) {
                            tim = timList2->at(u.second);
                            adv = advList->at(u.second);
                            delta2 = ((float) n * ((float) tim->hyper_degree[u.first] / tim->theta)) * adv->cpe -
                                     lambda * adv->seedUserCosts[u.first];
                            mariginal_count++;
                            if (first_u2) {
                                delta2 += lambda * B;
                                //first_u2= false;
                            }
                        }
                    } else {
                        S2full = true;
                    }
                    if (S1full && S2full) break;
                    int S_num;
                    if (delta1 > delta2) {
                        max_delta = delta1;
                        S_num = 1;
                    } else {
                        max_delta = delta2;
                        S_num = 2;
                    }

                    if (max_delta >= tau)
                    {
                        if (S_num == 1) {
                            S1_selected[u.first] = false;
                            tim = timList1->at(u.second);
                            adv = advList->at(u.second);

                            tim->opim_assign_best_node(u.first);
                            sum_node1++;

                            tim->currentRevenue += delta1;

                            first_u1 = false;

                        }
                        if (S_num == 2)
                        {
                            S2_selected[u.first] = false;
                            tim = timList2->at(u.second);
                            adv = advList->at(u.second);

                            tim->opim_assign_best_node(u.first);
                            sum_node2++;

                            tim->currentRevenue += delta2;

                            first_u2 = false;

                        }

                    } else if(max_delta>0){
                        N_temp.push_back(u);
                    }

                    now_index++;//M = M-{u}
                    if (max_delta > Z_max) Z_max = max_delta;

                }
                if (Z_max < 0) {
                    cout << "quit for Z_max!" << endl;
                    break;
                }
                if (N_temp.empty()) {
                    cout << "quit for N.empty!" << endl;
                    break;
                }
                if (tau < (tau_min / (1 + epsilon))) {
                    cout << "quit for tau!" << endl;
                    break;
                }
                if (S1full && S2full) {
                    cout << "quit for max_node!" << endl;
                    break;
                }
                tau = tau / (1.0 + epsilon);
                MorN=false;
            }
            else{
                now_index = 0;
                M.clear();
                Z_max = -999999999;
                while (now_index < N_temp.size()) {
                    //get first element;
                    //*
                    u.first = N_temp[now_index].first;
                    u.second = N_temp[now_index].second;
                    // */

                    delta1 = delta2 = -1;
                    martoid_count++;
                    if (sum_node1 < max_node) {
                        if (S1_selected[u.first]) {
                            tim = timList1->at(u.second);
                            adv = advList->at(u.second);
                            delta1 = ((float) n * ((float) tim->hyper_degree[u.first] / tim->theta)) * adv->cpe -
                                     lambda * adv->seedUserCosts[u.first];
                            mariginal_count++;
                            if (first_u1) {
                                delta1 += lambda * B;
                                //first_u1= false;
                            }
                        }
                    } else {
                        S1full = true;
                    }
                    martoid_count++;
                    if (sum_node2 < max_node) {
                        if (S2_selected[u.first]) {
                            tim = timList2->at(u.second);
                            adv = advList->at(u.second);
                            delta2 = ((float) n * ((float) tim->hyper_degree[u.first] / tim->theta)) * adv->cpe -
                                     lambda * adv->seedUserCosts[u.first];
                            mariginal_count++;
                            if (first_u2) {
                                delta2 += lambda * B;
                                //first_u2= false;
                            }
                        }
                    } else {
                        S2full = true;
                    }
                    if (S1full && S2full) break;
                    int S_num;
                    if (delta1 > delta2) {
                        max_delta = delta1;
                        S_num = 1;
                    } else {
                        max_delta = delta2;
                        S_num = 2;
                    }

                    if (max_delta >= tau) {
                        if (S_num == 1) {
                            S1_selected[u.first] = false;
                            tim = timList1->at(u.second);
                            adv = advList->at(u.second);

                            tim->opim_assign_best_node(u.first);
                            sum_node1++;

                            tim->currentRevenue += delta1;

                            first_u1 = false;

                        }
                        if (S_num == 2) {
                            S2_selected[u.first] = false;
                            tim = timList2->at(u.second);
                            adv = advList->at(u.second);

                            tim->opim_assign_best_node(u.first);
                            sum_node2++;

                            tim->currentRevenue += delta2;

                            first_u2 = false;

                        }

                    } else if(max_delta>0) {
                        M.push_back(u);
                    }

                    now_index++;
                    if (max_delta > Z_max) Z_max = max_delta;

                }
                if (Z_max < 0) {
                    cout << "quit for Z_max!" << endl;
                    break;
                }
                if (M.empty()) {
                    cout << "quit for N.empty!" << endl;
                    break;
                }
                if (tau < (tau_min / (1 + epsilon))) {
                    cout << "quit for tau!" << endl;
                    break;
                }
                if (S1full && S2full) {
                    cout << "quit for max_node!" << endl;
                    break;
                }
                tau = tau / (1.0 + epsilon);
                MorN=true;
            }
        }

        float totalDuration = getRunningTime(startTime); // in seconds
        float totalMemory = disp_mem_usage(); // in MB
        cout << "总时间开销: " << totalDuration << " s 总内存开销: " << totalMemory <<" MB "<< endl;

        // write results to master and adv specific files
        float total_revenue1 = 0.0;
        float total_revenue2 = 0.0;
        int total_seedsize1=0;
        int total_seedsize2=0;

        cout << "writing the results to output files.." << endl;

        for(int i = 0; i < nrCompanies; i++) //遍历每个广告商，输出到对应的txt文件中
        {
            adv = advList->at(i);
            tim1 = timList1->at(i);
            tim2 = timList2->at(i);

            total_revenue1 += tim1->currentRevenue;//累计整个分配的收益等
            total_revenue2 += tim2->currentRevenue;//累计整个分配的收益等
            total_seedsize1+=tim1->seedSet.size();
            total_seedsize2+=tim2->seedSet.size();
        }
        float total1 = 0.0;
        float total2 = 0.0;
        total1=total_revenue1;
        total2=total_revenue2;
        cout<<"S1:"<<endl;
        cout<<"f(S1): "<<total1<<endl;
        cout<<"Seed Size: "<<total_seedsize1<<endl;
        cout<<"S2:"<<endl;
        cout<<"f(S2): "<<total2<<endl;
        cout<<"Seed Size: "<<total_seedsize2<<endl;

        if(total1>total2)
        {
            cout<<"S*:"<<endl;
            cout<<"f(S*): "<<total1<<endl;
            cout<<"Seed Size: "<<total_seedsize1<<endl;
        }
        else{
            cout<<"S*:"<<endl;
            cout<<"f(S*): "<<total2<<endl;
            cout<<"Seed Size: "<<total_seedsize2<<endl;
        }
        cout<<"mariginal times: "<<mariginal_count<<endl;
        cout<<"martoid times: "<<martoid_count<<endl;
        cout<<"max node: "<<max_node<<endl;
        cout<<"eps: "<<epsilon<<endl;
    }

    void allocator::TwinGreedy()
    {
        advertiser *adv;
        TimGraph *tim1,*tim2;
        sfmt_t sfmt_seed;
        //随机种子
        for(int i = 0; i < 12; i++)//随机种子
            sfmt_init_gen_rand(&sfmt_seed, i+1234);
        //sfmt_init_gen_rand(&sfmt_seed,time(NULL));

        //generate rrset,then copy
        tim1 = timList1->at(0);
        tim1->theta_old=theta_0;
        tim1->theta=theta_0;
        tim1->RRset_generation();

        cout<<"RRset Size: "<<theta_0<<endl;
        TimGraph *tim=timList1->at(0);
        //S_1
        for(int i = 1; i < nrCompanies; i++)//h_1~h_n copy rrset from h_0
        {
            //generate 每个广告商的RR集
            tim1 = timList1->at(i);
            copy_tim(tim,tim1);
        }

        //S_2
        //*
        for(int i = 0; i < nrCompanies; i++)//h_0~h_n copy rrset from h_0
        {
            //generate 每个广告商的RR集
            tim2 = timList2->at(i);
            copy_tim(tim,tim2);
        }//*/

        for(int i = 0; i < nrCompanies; i++)//h_1~h_n copy rrset from h_0
        {
            //generate 每个广告商的RR集
            tim1 = timList1->at(i);
            tim2 = timList2->at(i);

            tim1->node_aval.resize(n,true);
            tim2->node_aval.resize(n,true);
        }
        int aval_num1=nrCompanies*n;
        int aval_num2=nrCompanies*n;
        //int random_num;

        long long int mariginal_count=0;
        long long int martoid_count=0;
        int sum_node1=0;
        int sum_node2=0;

        float max_revenue1=-999999999;
        pair<int,int> argmax_node1;
        bool need_argmax1=true;
        float max_revenue2=-999999999;
        pair<int,int> argmax_node2;
        bool need_argmax2=true;

        bool S1full=false;
        bool S2full=false;
        startTime = clock();
        while((aval_num1>0)||(aval_num2>0))
        {
            if(sum_node1<max_node)
            {
                if (need_argmax1)
                {
                    max_revenue1 = -999999999;
                    if(aval_num1>0)
                    {
                        for (int i = 0; i < nrCompanies; i++)//h_1~h_n copy rrset from h_0
                        {
                            //generate 每个广告商的RR集
                            tim = timList1->at(i);
                            adv = advList->at(i);
                            for (int j = 0; j < n; j++) {
                                martoid_count++;
                                if (tim->node_aval[j]) {
                                    mariginal_count++;

                                    float temp = ((float) n * ((float) tim->hyper_degree[j] / tim->theta)) * adv->cpe -
                                                 lambda * adv->seedUserCosts[j];
                                    if (first_u1) {
                                        temp += lambda * B;
                                    }
                                    if (temp > max_revenue1) {
                                        max_revenue1 = temp;
                                        argmax_node1.first = j;
                                        argmax_node1.second = i;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                S1full= true;
                max_revenue1 = -999999999;
            }
            if(sum_node2<max_node)
            {
                if (need_argmax2)
                {
                    max_revenue2 = -999999999;
                    if(aval_num1>0)
                    {
                        for (int i = 0; i < nrCompanies; i++)//h_1~h_n copy rrset from h_0
                        {
                            //generate 每个广告商的RR集
                            tim = timList2->at(i);
                            adv = advList->at(i);
                            for (int j = 0; j < n; j++) {
                                martoid_count++;
                                if (tim->node_aval[j]) {
                                    mariginal_count++;

                                    float temp = ((float) n * ((float) tim->hyper_degree[j] / tim->theta)) * adv->cpe -
                                                 lambda * adv->seedUserCosts[j];
                                    if (first_u2) {
                                        temp += lambda * B;
                                    }
                                    if (temp > max_revenue2) {
                                        max_revenue2 = temp;
                                        argmax_node2.first = j;
                                        argmax_node2.second = i;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                S2full= true;
                max_revenue2 = -999999999;
            }
            if(S1full&&S2full) break;
            if((max_revenue1<=0)&&(max_revenue2<=0))
            {
                cout<<"quit for argmax<=0!"<<endl;
                break;
            }

            if(max_revenue1>max_revenue2)
            {
                tim=timList1->at(argmax_node1.second);
                tim->opim_assign_best_node(argmax_node1.first);
                sum_node1++;
                first_u1=false;
                tim->currentRevenue+=max_revenue1;
                //cout<<"S1 choose ("<<argmax_node1.first<<" , "<<argmax_node1.second<<")"<<" currentRev"<<tim->currentRevenue<<endl;

                for(int i=0;i<nrCompanies;i++)
                {
                    tim=timList1->at(i);
                    tim->node_aval[argmax_node1.first] = false;
                    aval_num1--;
			
			
                tim=timList2->at(i);
                tim->node_aval[argmax_node1.first] = false;
                aval_num2--;
                }


                need_argmax1=true;
                need_argmax2=false;
            }
            else{
                tim=timList2->at(argmax_node2.second);
                tim->opim_assign_best_node(argmax_node2.first);
                sum_node2++;
                first_u2=false;
                tim->currentRevenue+=max_revenue2;
                //cout<<"S2 choose ("<<argmax_node2.first<<" , "<<argmax_node2.second<<")"<<" currentRev"<<tim->currentRevenue<<endl;


                for(int i=0;i<nrCompanies;i++)
                {
                    tim=timList2->at(i);
                    tim->node_aval[argmax_node2.first] = false;
                    aval_num2--;
			
		tim=timList1->at(i);
                tim->node_aval[argmax_node2.first] = false;
                aval_num1--;

                }


                need_argmax2=true;
                need_argmax1=false;
            }

        }

        float totalDuration = getRunningTime(startTime); // in seconds
        float totalMemory = disp_mem_usage(); // in MB
        cout << "总时间开销: " << totalDuration << " s 总内存开销: " << totalMemory <<" MB "<< endl;

        // write results to master and adv specific files
        float total_revenue1 = 0.0;
        float total_revenue2 = 0.0;
        int total_seedsize1=0;
        int total_seedsize2=0;

        cout << "writing the results to output files.." << endl;

        for(int i = 0; i < nrCompanies; i++) //遍历每个广告商，输出到对应的txt文件中
        {
            adv = advList->at(i);
            tim1 = timList1->at(i);
            tim2 = timList2->at(i);

            total_revenue1 += tim1->currentRevenue;//累计整个分配的收益等
            total_revenue2 += tim2->currentRevenue;//累计整个分配的收益等
            total_seedsize1+=tim1->seedSet.size();
            total_seedsize2+=tim2->seedSet.size();
        }
        float total1 = 0.0;
        float total2 = 0.0;
        total1=total_revenue1;
        total2=total_revenue2;
        cout<<"S1:"<<endl;
        cout<<"f(S1): "<<total1<<endl;
        cout<<"Seed Size: "<<total_seedsize1<<endl;
        cout<<"S2:"<<endl;
        cout<<"f(S2): "<<total2<<endl;
        cout<<"Seed Size: "<<total_seedsize2<<endl;

        if(total1>total2)
        {
            cout<<"S*:"<<endl;
            cout<<"f(S*): "<<total1<<endl;
            cout<<"Seed Size: "<<total_seedsize1<<endl;
        }
        else{
            cout<<"S*:"<<endl;
            cout<<"f(S*): "<<total2<<endl;
            cout<<"Seed Size: "<<total_seedsize2<<endl;
        }
        cout<<"mariginal times: "<<mariginal_count<<endl;
        cout<<"martoid times: "<<martoid_count<<endl;
        cout<<"max node: "<<max_node<<endl;
    }
    
    void allocator::copy_tim(TimGraph *src, TimGraph *dest)
    {
        dest->theta=src->theta;
        //dest->theta_old=src->theta;
        //清空原vector数组并复制
        dest->hyperGT_adv.assign(src->hyperGT_adv.begin(),src->hyperGT_adv.end());
        dest->isCovered.assign(src->isCovered.begin(),src->isCovered.end());
        dest->hyperG_adv.assign(src->hyperG_adv.begin(),src->hyperG_adv.end());
        dest->hyper_degree.assign(src->hyper_degree.begin(),src->hyper_degree.end());
    }

    void allocator::readTICGraph() {
        
        string probGraphFile = opt->getValue("probGraphFile");//社交网络图，每行是一条边，和在对应主题上生效的概率
        cout << "Reading file " << probGraphFile << endl;
        ifstream myfile (probGraphFile.c_str(), ios::in);
        
        float *dists;
        float p;
        advertiser *advTemp;
        
        int nrEdges = 0;
        set<int> nodes; // kontrol amacli simdi
        
        if (myfile.is_open()) {
            while (! myfile.eof() )	{
                std::string line;
                getline (myfile,line);
                if (line.empty()) continue;
                nrEdges++;
                
                std::string::size_type pos = line.find_first_of(delim);
                int prevpos = 0;
                
                //first user
                string str = line.substr(prevpos, pos-prevpos);
                int u1 = strToInt(str);
                
                //second user
                prevpos = line.find_first_not_of(delim, pos);
                pos = line.find_first_of(delim, prevpos);
                int u2 = strToInt(line.substr(prevpos, pos-prevpos));
                
                if (u1 == u2)//如果边是结点自身到自身，说明有误，直接跳过
                    continue;
                
                graphT[u2].push_back(u1); //insert to the transposed graph，为每个目标结点保存入边的源结点
                
                // kontrol amacli
                nodes.insert(u1);//nodes是一个set集合，因此用来统计所有结点
                nodes.insert(u2);
                
                prevpos = line.find_first_not_of(delim, pos);
                
                str = line.substr(prevpos);
                dists = new float[nrTopics];
                stringTokenizer(str, dists, nrTopics, delim);//把在所有主题上的概率保存下来
                
                for(int i = 0; i < nrCompanies; i++)
                {
                    advTemp = advList->at(i);//临时广告商对象，遍历所有广告商
                    p = 0.0;
                    for(int j = 0; j < nrTopics; j++)
                        p += (dists[j] * advTemp->gamma[j]);//遍历所有主题，边的概率乘以广告商在主题上的分布概率才是真正的影响概率p^i_{uv}
                    //累加起来得到的p就是这只广告在这个边上可能会产生影响力的概率（没错，一只广告所有主题上概率的累加才是合理的）。
                    advTemp->probT[u2].push_back(p);//保存一个广告商，在选择一个用户激活宣传后，该用户通过不同的边产生影响力的概率
                }
            }
            
            myfile.close();
        }
        
        else
            cout << "Can't open friendship graph file " << probGraphFile << endl;
        
        cout << "Built transposed adj matrix from file " << endl;
        cout << "number of nodes " << nodes.size() << endl;
        cout << "number of edges " << nrEdges << endl;
        
    }
    
    void allocator::readItemDistsFile() {
        
        cout << "reading item distributions file " << endl;
        string itemDistsFile = opt->getValue("itemDistsFile");
        ifstream myfile(itemDistsFile.c_str(), ios::in);
        float *tokens;
        
        int advIndex = 0;
        
        if(myfile.is_open()) {
            while(!myfile.eof()) {
                std::string line;
                getline(myfile, line);//逐行读取
                if(line.empty())
                    continue;
                tokens = new  float[nrTopics];//主题数大小的float数组
                stringTokenizer(line, tokens, nrTopics, delim);//读取每行的广告商在不同主题上的概率分布
                advList->at(advIndex++)->setItemDist(tokens, nrTopics);
                if(advIndex >= nrCompanies)//超过了公司数，文件中剩余的概率分布就不需要了
                    break;
            }
            
            myfile.close();
        }
        
        else {
            cout << "problem opening the item distributions file, exiting... " << itemDistsFile <<  endl;
            exit(1);
        }
        
    }

    // this reads a cost file in the form of a cost matrix in the form of nodes X ads -- might be different for the scalability version
    void allocator::readIncentiveCosts() {
        //每行表示一个用户被不同广告商收买所需的费用
        //该费用通过不同的模型由影响力估得
        string readIncentiveCostsFile = opt->getValue("incentiveCostsFile");
        ifstream myfile(readIncentiveCostsFile.c_str(), ios::in);
        
        int lineIndex = 0;
        float *tokens;
        float tempCostToken = 0.0;


        if(myfile.is_open()) {
            while(!myfile.eof())
            {
                std::string line;
                getline(myfile, line);
                if(line.empty())
                    continue;
                tokens = new float[nrCompanies];
                stringTokenizer(line, tokens, nrCompanies, delim);//按广告商数目切割每行
/*
                for (int i = 0; i < nrCompanies; i++)
                {
                    if(tokens[i] < 1.0)
                        tokens[i]= 1.1;
                }*/
                for(int i = 0; i < nrCompanies; i++)
                {
                    //采用不同的激励模型，alpha是控制大小的系数
                    //线性
                    if(string(opt->getValue("costFunctionType")).compare("l") == 0) { // linear
                        tempCostToken = tokens[i] * alpha;
                    }
                    //所有用户相同的，constan
                    else if(string(opt->getValue("costFunctionType")).compare("u") == 0) { // uniform -- reads uniform input
                        tempCostToken = tokens[i] * alpha;
                    }
                    //二次函数
                    else if(string(opt->getValue("costFunctionType")).compare("q") == 0) { // quadratic
                        tempCostToken = tokens[i] * tokens[i] * alpha;
                    }
                    //次线性，即Log
                    else if(string(opt->getValue("costFunctionType")).compare("s") == 0) { // sublinear
                        tempCostToken = log(tokens[i]) * alpha;
                        if(tempCostToken==0) tempCostToken=log(tokens[i]+0.1) * alpha;
                    }
                    //随机的
                    else if(string(opt->getValue("costFunctionType")).compare("r") == 0) { // random
                        
                        tempCostToken = tokens[i] * alpha; // daha sonra shuffle edilecek bu
                    }

                    //依次存储到list数据结构中，每次循环对应存储到当前广告商所对应于每个用户标号的cost
                    advList->at(i)->seedUserCosts[lineIndex] = tempCostToken;
                    //保存对于一个广告商来说，最贵的代言人和最便宜的代言人
                    if(advList->at(i)->maxCost < tempCostToken)
                        advList->at(i)->maxCost = tempCostToken;
                    
                    if(advList->at(i)->minCost > tempCostToken)
                        advList->at(i)->minCost = tempCostToken;

                }
                B+= tempCostToken;
                lineIndex++;
            }
            
            myfile.close();
            
            //            cout << " coost - kontrol : " << advList->at(0)->seedUserCosts[0] << endl;
            // for random cost function shuffle the vectors here
            //如果选择激励费用是随机的，则之前正常读入cost，然后在这里打乱
            if(string(opt->getValue("costFunctionType")).compare("r") == 0)
            { // random
                for (int i =0 ; i < nrCompanies; i++)
                {
                    //cout << "adv " << i << " cost kontrol " << advList->at(i)->maxCost << endl;
                    //                    cout << "randomness kontrol b4 : " << advList->at(i)->seedUserCosts[0] << endl;
			        std::srand(std::time(0));
                    std::random_shuffle(advList->at(i)->seedUserCosts.begin(), advList->at(i)->seedUserCosts.end());
                    //                    cout << "randomness kontrol after : " << advList->at(i)->seedUserCosts[0] << endl;
                }
            }
            //输出对于每个广告商来说，最贵的代言人和最便宜的代言人
            //*
            cout << "cost control " << endl;
            for (int i = 0; i < nrCompanies; i++)
            {
                cout << "max for adv " << i << " " << advList->at(i)->maxCost << endl;
                cout << "min for adv " << i << " " << advList->at(i)->minCost << endl;
            }//*/

        }
        
        else {
            cout << "problem opening the incentive costs file, exiting... " << readIncentiveCostsFile <<  endl;
            exit(1);
        }
        
    }
    
    allocator::~allocator()
    {
        //cout << "destructor called " << endl;
    }
}
