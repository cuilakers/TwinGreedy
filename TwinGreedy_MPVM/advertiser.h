#ifndef _ADVERTISER_H
#define _ADVERTISER_H

namespace _Cide{

    class advertiser{
        
    public:
        
        std::vector< std::vector< float> > probT;//保存一个广告商，在选择一个用户激活宣传后，该用户通过不同的边产生影响力的概率（也即一个结点不同边上的概率）
        int advertiserID;//广告商的ID
        float *gamma;//广告商在不同主题上的分布概率
        std::vector<float> seedUserCosts; // advertiser-specific 每个广告商所选的每个种子的花费？
        float maxCost; // to be used in the timgraph object
        float minCost; // to be used in the timgraph object

        float cpe;
        
        advertiser(int id, int nrTopics)
        {
            this->advertiserID = id;
            this->gamma = new  float[nrTopics];
            this->cpe=1;
        }
        
        ~advertiser() {}
        
        void setItemDist( float *temp, int size) {
            for(int i = 0; i < size; i++)
                this->gamma[i] = temp[i];//设置每个广告的概率分布，保存到数组结构中
        }

    };


}


#endif
