#include"generate.cpp"
#include <iostream>
using namespace std;
int main() {
    tdumatrix u;
   u = generate_tdu_matrix(0.24);
   cout<< "Generated TDU Matrix:" << endl;
   cout<< u(0,0) << " " << u(0,1) << endl;
   cout<< u(1,0) << " " << u(1,1) << endl;
    cout<< "Determinant: " << u.determinant() << endl;
   cout<< "Determinant: " << u.determinant()*conj(u.determinant()) << endl;

  cout<< "input steps:" << endl;
  int steps;
  cin>> steps ;
  double movestep;
  movestep = 0.334;
    double integral_sum = 0.0; // 用于累加 g(x) 的值
    std::vector<double> results; // 存储每一步的平均值

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
integral_sum = integral_sum + real(conj(lati(0,0,0).determinant()) *
                                 conj(lati(3,0,1).determinant())*lati(3,3,0).determinant() * lati(0,3,1).determinant());;

        results.push_back(integral_sum / (i + 1.0));
    }
    cout<< results[steps-1]<< endl;










return 0;
}




