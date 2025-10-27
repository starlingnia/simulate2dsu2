#include <complex>
#include <vector>
#include <cmath>
#include <random>
using dcomplex = std::complex<double>;
   static std::random_device rd;
   static std::mt19937 gen(rd());
   static std::uniform_real_distribution<> dis(0.0, 1.0);


typedef struct randommatrix {
  virtual std::complex<double> determinant() = 0;
  int N;
  std::vector<dcomplex> a;
  randommatrix( int n) : N(n), a(n * n, dcomplex(0.0, 0.0)) {}
  dcomplex& operator()(int i, int j){ return a[i*N + j]; }
  dcomplex trac() {
    dcomplex tr = dcomplex(0.0, 0.0);
    for (int i = 0; i < N; ++i) {
      tr += a[i * N + i];
    }
    return tr;
  }
} randommatrix;

typedef struct tdumatrix: public randommatrix {
   tdumatrix() : randommatrix(2) {}
   std::complex<double> determinant() override {
    return a[0]*a[3] - a[1]*a[2];
  }
    tdumatrix conj() {
    tdumatrix result;
    result(0,0) = std::conj(a[0]);
    result(0,1) = std::conj(a[2]);
    result(1,0) = std::conj(a[1]);
    result(1,1) = std::conj(a[3]);
    return result; }
    tdumatrix operator*(tdumatrix B) {
      tdumatrix C;
     C(0,0) = a[0]*B(0,0) + a[1]*B(1,0);
     C(0,1) = a[0]*B(0,1) + a[1]*B(1,1);
    C(1,0) = a[2]*B(0,0) + a[3]*B(1,0);
     C(1,1) = a[2]*B(0,1) + a[3]*B(1,1);
      return C;
    } 

} td;

tdumatrix generate_tdu_matrix(double epsilon) {
  tdumatrix M;
  double theta1 = 2.0 * M_PI * dis(gen);
  double theta2 = 2.0 * M_PI * dis(gen);
  double thetaf = 2.0 * M_PI * dis(gen);
  double r1 = std::sqrt(-std::log(dis(gen))) * epsilon;
  double r2 = std::sqrt(1.0 - r1 * r1);

  M(0, 0) = dcomplex(std::cos(theta1), std::sin(theta1)) * (r1);
  M(0, 1) = dcomplex(std::cos(theta2), std::sin(theta2)) * r2;
  M(1, 0) = dcomplex(std::cos(theta2+thetaf), -std::sin(theta2+thetaf)) * r2;
  M(1, 1) = dcomplex(-std::cos(theta1+thetaf), std::sin(theta1+thetaf)) * (r1);

  return M;
}

tdumatrix multi(tdumatrix A, tdumatrix B) {
  tdumatrix C;
  C(0, 0) = A(0, 0) * B(0, 0) + A(0, 1) * B(1, 0);
  C(0, 1) = A(0, 0) * B(0, 1) + A(0, 1) * B(1, 1);
  C(1, 0) = A(1, 0) * B(0, 0) + A(1, 1) * B(1, 0);
  C(1, 1) = A(1, 0) * B(0, 1) + A(1, 1) * B(1, 1);
  return C;
}
typedef struct  latticesitedof
{ 
  int M = 6;
  int N = 6;  
  std::vector <tdumatrix > site;
  std::vector <tdumatrix> ssite;
  // 添加构造函数来初始化向量的大小
  latticesitedof() : site(M * N), ssite(M * N) {} // 初始化为 16 个 tdumatrix 的向量
  
 tdumatrix& operator()(int i, int j, int mu){ 
  if (mu ==0)
  return site[i*M + j];
  else
  return ssite[i*N + j]; }
}dof;
// 移除 randomdof 中的 push_back，直接使用索引赋值
dof randomdof(double move) {
  dof lat; // 调用新的构造函数，site 和 ssite 已有 16 个元素
  for (int i = 0; i < 36; ++i) {
    // 使用索引而不是 push_back
    lat.site[i] = generate_tdu_matrix(move * dis(gen));
    lat.ssite[i] = generate_tdu_matrix(move * dis(gen));
  }
  return lat;
}

dof updata(dof lat, dof ran) {
  dof latnew = lat;
  int size = lat.M * lat.N; // 16
  
  for (int i = 0; i < size; ++i) {
    // site
    tdumatrix delta_site = ran.site[i]; // ran.site 只包含 16 个
    latnew.site[i] = multi(lat.site[i], delta_site);
    
    // ssite - 也要更新
    tdumatrix delta_ssite = ran.ssite[i];
    latnew.ssite[i] = multi(lat.ssite[i], delta_ssite);

  }return latnew;
}

double target_distribution(dof x){
  dcomplex dete = dcomplex(0.0, 0.0);
  int M = x.M, N = x.N; // 4x4 
  for (int i = 0; i < M; ++i) { // i for row/site
    for (int j = 0; j < N; ++j) { // j for column/site
     // 确保索引 i, j 分别对应 M 和 N 的模
     // x((i+1)%M, (j+1)%N, 0)
     // x(i, (j+1)%N, 1)
     // x(i, j, 0)
     // x((i+1)%M, j, 1)
     // 假设 M 是行（垂直），N 是列（水平）

     dete = dete + multi(multi(x((i+1)%M,(j+1)%N,0).conj(), x(i,(j+1)%N,1).conj()), multi(x(i,j,0), x((i+1)%M,j,1))).trac();
   
   
    }
  }
  return exp(-0.5 * real(dete)) / sqrt(2 * M_PI);
}