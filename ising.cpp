#include "ising.h"

/************************************************************/
/*                         构造函数                          */
/************************************************************/
chain::chain(int L)
{
	len = L;

	//初始化长链状态
	for(int i = 0; i != L; ++i) {
		if(get_random()<0.5)
			state.push_back(1);
		else
			state.push_back(-1);
	}

	//初始化无向图链接矩阵
	table=new int*[L];
	for(int i=0;i<L;i++)
		table[i]=new int[L];

	for(int i=0;i<L;i++) {
		for(int j=0;j<L;j++)
			table[i][j]=MAXMIN;
	}

	for(int i=0; i<L; ++i) {
		for(int j=0; j<L; ++j) {
			if( abs(i-j)==1 || abs(i-j)==L-1 ) {
				table[i][j] = 1;
				table[j][i] = 1;
			}

		}
	}

	//初始化储存连接点的vector
	vector<int> a;
	a.push_back(1);a.push_back(L-1);
	link_point.push_back(a);
	for(int i=1; i != L-1; ++i) {
		vector<int> a;
		a.push_back(i-1);a.push_back(i+1);
		link_point.push_back(a);
	}
	vector<int> b;
	b.push_back(0);b.push_back(L-2);
	link_point.push_back(b);	
}

/****************************************************************************/
/*                               搜索最短路径                                 */
/****************************************************************************/
vector<int> chain::path(int k)
{  
	
	int *dist,*prev;
	dist=new int[len];
	prev=new int[len];


    int i,j;
    bool s[MAXMIN];   //maxint是个非常大的数
	int count=1;

    for (i=0;i<len;++i) {
        dist[i] = table[k][i];
        s[i] = false;
        if (dist[i] == MAXMIN) prev[i] = 0;   //将该点的前一个点赋为0,应为它不与v点直接相连
        else prev[i] = k;   
    }

    dist[k] = 0; s[k] = true;      //与prim不同的是初始时从源点出发 
    for (i=0;i<len-1;++i) {
        int temp = MAXMIN;
        int u = k;
        for (j=0;j<len;++j) {
            if ((!s[j])&&(dist[j]<temp)) {         //先找和V点距离最短并且没有被访问过的点u {
                u = j;                           
                temp = dist[j];
            }
        }

        s[u] = true;                              //找到之后访问
		
        for (j=0;j<len;++j) {
            if ((!s[j])&&table[u][j]<MAXMIN) {            //再找与u点相邻的点
                int newdist = dist[u] + table[u][j];     
                if (newdist<dist[j]) {
                    dist[j] = newdist;           //如果存在点与u点的距离加上u点到原点的距离比原dist[j]要小，那么newidist为从源点到该点的最短特殊路径
                    prev[j] = u;       //将j的前驱记为u
                }
            }
        }
    }

    vector<int> shortpath;
	for(j=0;j<len;j++) {
		int num;
		num=j;
		count=1;
		if(j==k) {
			shortpath.push_back(0);
			continue;
		}
		else {
			if(prev[j]==k) {
				 shortpath.push_back(1);
				 continue;
			 }
			 else {
				while(prev[num]!=k) {
					count++;
					num=prev[num];
				}
				shortpath.push_back(count);
			 }
		 }
	 }

	 return shortpath;
			 
}
/**************************************************************************************/
/*                                  建立长程连接                                        */
/**************************************************************************************/
void chain::long_link()
{
	vector<int> dis_vector;
	int j; //与点i相连接的点
	for(int i = 0; i != len; ++i) {

		if(get_random() < probability) {
			//找出每个点和i点的拓扑距离
			dis_vector = path(i);
			//找出最远距离是多少
			int max_distance = vector_max(dis_vector);
			//生成累计概率矩阵
			vector<double> prob = prob_vector(2, max_distance);
			//找出到底和哪个距离的点连接
			double p = get_random();			
			int count = 0;
			for(vector<double>::iterator iter=prob.begin(); iter!=prob.end(); iter++) {
				if( p < *iter)
					break;
				count++;	
			}
			int dist = max_distance - count;
			//找出具有该dist的点d的序号
			vector<int> index_vector = index_find(dis_vector, dist);
			//随机选出一个点与第i点连接，因为可能同时有几个点与点i具有相同距离
			if(index_vector.empty())
				cerr<<"Wrong"<<endl;
			else {
				int size_index = index_vector.size();
				int location = floor(size_index * get_random());

				j = index_vector[location];
				//修改储存连接点的二维vector
				link_point[i].push_back(j);
				//修改无向图的连接矩阵
				table[i][j] = 1;
				table[j][i] = 1;
			}

		}

	}
}
/*************************************************************************************/
/*                             随机选取一个点进行翻转操作                                */
/*************************************************************************************/
void chain::manipulate(double temp)
{
	int k = floor(get_random() * len);
	//计算相邻点磁矩之和
	int sum = 0;
	for(vector<int>::iterator iter=link_point[k].begin(); iter!=link_point[k].end(); iter++) {
		sum += state[*iter];
	}
	//计算E0
	int E0 = -state[k] * sum;
	//计算E1
	int reverse = -state[k];
	int E1 = -reverse * sum;

	//判断是否翻转
	if(E1 < E0) {
		state[k] = reverse;
	}
	else {
		double line = exp(-(E1 - E0)/(kb * temp)); //是否接受新状态的概率
		if(get_random() < line)
			state[k] = reverse;
	}
}

/*************************************************************************************/
/*                             一个蒙特卡洛步数                                        */
/*************************************************************************************/
void chain::mcs(double temp)
{
	for(int i=0; i < len; ++i) {
		manipulate(temp);
	}
}

/*************************************************************************************/
/*                           得到体系当前磁化强度                                       */
/*************************************************************************************/

double chain::get_magnetization()
{
	double s;
	for(vector<int>::iterator iter=state.begin(); iter!=state.end(); iter++) 
		s += *iter;
	return fabs(s)/len;
}

/***************************************************************************************/
/*                                辅助函数                                              */
/***************************************************************************************/
vector<int> index_find(vector<int> dis_vector, int dist)
{
	vector<int> result;
	int count = 0;
	for(vector<int>::iterator iter=dis_vector.begin(); iter!=dis_vector.end(); iter++) {
		if( *iter == dist)
			result.push_back(count);
		count++;	
	}
	return result;
}

//生成累计概率向量
vector<double> prob_vector(const int min, const int max)
{
	double s = 0;
	vector<double> prob;
	for(int l = max; l >= min; --l) {
		s +=  pow(l,-m);
		prob.push_back(s);
	}
	for(vector<double>::iterator iter=prob.begin(); iter!=prob.end(); iter++) {
		*iter /= s;
	}
	return prob;
}


//生成0-1上随机数
double get_random()
{
	return rand()/(double)RAND_MAX;
}

int vector_max(vector<int> mmm)
{
	int max = mmm[0];
	for(vector<int>::iterator iter=mmm.begin(); iter!=mmm.end(); iter++) {
		if( *iter > max)
			max = *iter;
	}
	return max;
}


