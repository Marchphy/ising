#ifndef ISING_H
#define ISING_H


#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#define MAXMIN 9999

#define LENGTH 501
#define m 0.2
#define probability 1
#define kb 1
#define step 10000


using namespace std;

class chain {
	private:
		vector<int> state; //描述系统当前状态

		vector<vector<int> > link_point;

		int len; //长链的晶格数目
		int **table; //储存无向图的链接矩阵

	public:
		chain(int L); //有参数构造函数，L为格点数目，随机分配+1或者-1给各个格点

		vector<int> path(int k); //返回一个vector，包含k点到所有点的最短拓扑距离

		void long_link(); //建立长程连接

		void manipulate(double temp); //随机选取一个点进行操作

		void mcs(double temp); //一个蒙特卡洛不

		double get_magnetization(); //得到磁化强度

		vector<vector<int> > aaaaaa()
		{
			return link_point;
		}

		int** bbbbbbb()
		{
			return table;
		}

};





//生成0-1上随机数
double get_random();

//返回一个vector中具有特定值的序号
vector<int> index_find(vector<int> dis_vector, int dist);

//生成概率分布vector
vector<double> prob_vector(const int min,const int max);

//返回vector<int> 中的最大值
int vector_max(vector<int> mmm);



#endif