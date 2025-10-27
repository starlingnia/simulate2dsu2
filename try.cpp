#include"generate.cpp"
#include <iostream>
#include <fstream>
using namespace std;
int main() {
  cout<< "input steps:" << endl;
  int steps;
  cin>> steps ;
  double movestep;
  movestep = 0.544;
   std::ofstream file("results.txt");
   file << "step" << "\t" << "integral_value\n";
  double wloop;

dof lati = randomdof(dis(gen));
for (int i = 0; i < steps; ++i) {
        // 生成一个随机的提议状态 (proposal)
        dof lwm = randomdof(movestep*dis(gen));

       dof ln = updata(lati, lwm);

        double alpha = target_distribution(ln) / target_distribution(lati);
        
        // 接受或拒绝提议状态
        if (dis(gen) < alpha) {
            lati = ln; // 接受提议
        }
        
        // 累加函数 g(x) = x^2 的值
   wloop = real((lati(4,4,0).conj() * lati(3,4,1).conj()*lati(3,3,0).conj() * lati(4,3,1).conj()).trac());
    file << i << "\t" << wloop<< "\t" << target_distribution(lati) << "\n";
       
    }

    file.close();

return 0;
}




