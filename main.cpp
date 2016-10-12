#include "ising.h"

int main()
{
	srand((unsigned)time(NULL));
	chain march(LENGTH);

	march.long_link();


	ofstream myfile("data.txt",ios::out);


	for(double temp = 0.1; temp < 0.5 ;temp += 0.01 ) {

		double m_sum = 0.0;
		double m_2_sum = 0.0;
		//double m_4_sum = 0.0;
		for(int i=0; i != step; ++i) {
			march.mcs(temp);
			if(i >= step/2) {
				m_sum += march.get_magnetization();
			//cout<<march.get_magnetization()<<endl;
				m_2_sum += pow(march.get_magnetization(),2);
			//m_4_sum += pow(march.get_magnetization(),4);
				cout<<i<<endl;
			}
	}

	//cout<<"\t"<<m_2_sum<<"\t"<<m_4_sum<<endl;

	//计算磁化率
	double susceptibility = (m_2_sum/(step/2) - (m_sum/(step/2)) * (m_sum/(step/2))) / (LENGTH * temp);
	//计算累积量
	//double Binder = 1 - (m_4_sum/(step/2)) / (3 * (m_2_sum/(step/2)) * (m_2_sum/(step/2)));

	myfile<<temp<<"\t"<<susceptibility<<endl;
	//myfile<<temp<<"\t"<<Binder<<endl;

}
	myfile.close();

	return 0;
}
