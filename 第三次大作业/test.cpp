#include<algorithm>
#include<iostream>
#include<fstream>
#include<iostream>
#include<complex>
using namespace std;
using com = complex<double>;
//x轴离散成2*N0个格子
constexpr int N0 = 20000;
//x间隔为delta_x
constexpr double delta_x = 0.1;
//
constexpr double E = 2.71828182845904523536028747;
constexpr double PI = 3.14159265358979323846264338;
//归一化
template<int n>
void normalize(com* x) {
	double norm=0;
	for (int i = 0; i < n; i++) {
		norm += pow(abs(x[i]),2)*delta_x;
	}
	norm = sqrt(norm);
	for (int i = 0; i < n; i++) {
		x[i] = x[i] / norm;
	}
}
//特征矢+特征值
struct eigen {
	com* vec;
	double val;
};
//打印向量
template<int n>
void print_sol(com x[n]) {
	for (int i = 0; i < n; i++) { cout << i + 1 << "   " << x[i] << endl; }
}

//解对称三对角线性方程组
template<int n>
com* solve_eq(com c_[n][2], com b_[n]) {
	
	//complex c[n][2];
	com(*c)[2] = new com[n][2];
	com* b = new com[n];
	for (int i = 0; i < n; i++) {
		c[i][0] = c_[i][0];
		c[i][1] = c_[i][1];
		b[i] = b_[i];
	}
	//计算L矩阵元
	for (int i = 1; i <= n; i++) {
		int r = max(1, i - 1);
		for (int j = r; j <= i; j++) {
			for (int k = r; k <= j - 1; k++) {
				//由于数组下标从0开始,故需要-1
				c[i - 1][j - i + 2 - 1] = c[i - 1][j - i + 2 - 1] - c[i - 1][k - i + 2 - 1] * c[j - 1][k - j + 2 - 1] / c[k - 1][2 - 1];
			}
		}
	}
	//计算b'
	for (int i = 1; i <= n; i++) {
		int r = max(1, i - 1);
		for (int j = r; j <= i - 1; j++) {
			b[i - 1] = b[i - 1] - c[i - 1][j - i + 2 - 1] * b[j - 1] / c[j - 1][2 - 1];
		}
	}
	//上三角回代
	for (int i = n; i >= 1; i--) {
		int t = min(n, i + 1);
		for (int j = i + 1; j <= t; j++) {
			b[i - 1] = b[i - 1] - c[j - 1][i - j + 2 - 1] * b[j - 1];
		}
		b[i - 1] = b[i - 1] / c[i - 1][2 - 1];
	}
	delete[] c;
	return b;
}

//反幂法求本征问题
template<int N>
eigen inver_exp(com C[2*N-1][2], com u_[2*N-1]) {
	//初始化u,lambda
	constexpr int n = 2 * N - 1;
	com* u=new com[n];
	for (int i = 0; i < n; i++) {
		u[i] = u_[i];
	}
	com lambda = 0;
	com* v=new com[n];
	for (int k = 0; k < 10000; k++) {
		delete[] v;
		v = solve_eq<n>(C, u);
		//求v绝对值最大的分量max
		com max = 0;
		for (int i = 0; i < n; i++) {
			if (abs(max) < abs(v[i])) {
				max = v[i];
			}
		}
		//储存前一次本征值
		com lambda_before = lambda;
		lambda = 1. / max;
		//u=v*lambda
		for (int i = 0; i < n; i++) {
			u[i] = v[i]*lambda;
		}
		//检查精度是否合适
		if (abs(lambda_before - lambda)/abs(lambda) < 1e-10) {
			break;
		}
	}
	delete[] v;
	return eigen{ u,abs(lambda) };
}


//乘高斯吸收函数
template<int n>
void absorb(com A[n]) {
	double x0 = 0.75*N0*delta_x;
	double x = -N0 * delta_x;
	for (int i = 0; i < n; i++) {
		x += delta_x;
		if (abs(x) > x0) {
			A[i] = A[i] * pow(E, -pow((abs(x) - x0) / 0.2, 2));
		}

	}
}

//波函数传播delta_t
//其中I(10^16为单位),omega(原子单位),N是激光场参数
template<int n>
void Udelta_t(com psai[n],double t,double I,double omega,int N,double delta_t) {
	//构造哈密顿量H,装在C中
	com(*C)[2] = new com[n][2];
	//constexpr int n = 2 * N0 - 1;
	C[0][0] = 0;
	double x = -N0 * delta_x;//x是坐标
	for (int i = 1; i < n; i++) {//构造非对角元
		C[i][0] = -0.5 / pow(delta_x, 2);
	}
	//带入原子单位换算公式求sqrt(I)
	const double sqrtI=sqrt(I/3.5094448314);
	for (int i = 0; i < n; i++) {//构造对角元
		x += delta_x;
		C[i][1] = x*sqrtI*pow(sin(omega*t/2./N),2)*sin(omega*t)
			+1 / pow(delta_x, 2) - 1 / sqrt(2 + pow(x, 2)) ;
	}
	for (int i = 0; i < n; i++) {//*0.5*i*delta_t
		com i0{ 0.,1. };
		C[i][1] = C[i][1] * (-0.5*delta_t*i0);
		C[i][0] = C[i][0] * (-0.5*delta_t*i0);
	}


	//计算psai*(1-0.5*i*H*delta_t)
	com *psai_temp=new com[n];//临时储存psai
	for (int i = 0; i < n; i++) {
		psai_temp[i] = psai[i];
	}
	for (int k = 1; k < n ; k++) {//开始计算
		psai[k] = psai_temp[k] * (1. - C[k][1]) - psai_temp[k - 1] * C[k][0] - psai_temp[k + 1] * C[k + 1][0];
	}
	psai[0]= psai_temp[0] * (1. - C[0][1]) - psai_temp[1] * C[1][0];//计算端点值
	psai[n-1]= psai_temp[n-1] * (1. - C[n-1][1]) - psai_temp[n-2] * C[n-1][0];
	//计算psai*(1+0.5*i*H*delta_t)^-1
	for (int i = 0; i < n; i++) {//临时储存psai
		psai_temp[i] = psai[i];
	}
	for (int i = 0; i < n; i++) {//构造系数矩阵(1+0.5*i*H*delta_t)
		C[i][1] = 1. + C[i][1];
	}
	delete[] psai;
	psai = solve_eq<n>(C, psai_temp);//开始计算
	delete[] C;
	delete[] psai_temp;
}
//生成t=0的波函数
com* generate_t0() {
	//构造三对角型的系数矩阵A,装在C中
	com (*C)[2]=new com[2*N0-1][2];
	constexpr int n = 2 * N0 - 1;
	C[0][0] = 0;
	double x = -N0 * delta_x;
	for (int i = 1; i < n; i++) {//构造非对角元
		C[i][0] = -0.5 / pow(delta_x, 2);
	}
	for (int i = 0; i < n; i++) {//构造对角元
		x += delta_x;
		C[i][1] = 1 / pow(delta_x, 2) - 1 / sqrt(2 + pow(x, 2)) + 0.48;
	}
	//u 为初始向量
	com* u = new com[n];
	for (int i = 0; i < n; i++) {
		u[i] = 1 / sqrt(2 * N0*delta_x);
	}

	eigen psai0 = inver_exp<N0>(C, u);
	normalize<n>(psai0.vec);
	delete[] u,C;
	return psai0.vec;
}
//生成tx时刻波函数
com* generate_tf(double tx, double I, double omega, int N) {
	//t间隔delta_t
	const double delta_t = 0.5;//!!!!
	//获得基态波函数
	com* psai = generate_t0();
	for (double t = 0; t < tx + delta_t; t += delta_t) {
		Udelta_t<2*N0-1>(psai, t,I,omega,N,delta_t);//波函数传播
		absorb<2 * N0 - 1>(psai);//乘吸收函数
		if(1){
			cout << psai[N0] <<"  "<<norm(psai[N0])<< endl;//!!!!!!!
		}
	}
	return psai;
}





//求末态电离波函数(针对第(2)问)
com* generate_psaif() {
	//初始化激光数据
	double I = 1.;
	double omega = 1.;
	int N = 18;
	double tf = 2 * N*PI / omega;
	//末态波函数
	com* psaif = generate_tf(tf, I, omega, N);
	//初态波函数
	com* psai0 = generate_t0();
	//基态概率幅p
	com p=0;
	for (int i = 0; i < 2 * N0 - 1; i++) {
		p += conj(psai0[i])*psaif[i]*delta_x;
	}
	//减去基态部分
	for (int i = 0; i < 2 * N0 - 1; i++) {
		psaif[i] = psaif[i] - p * psai0[i];
	}
	delete[] psai0;
	return psaif;
}














//第1问测试函数
void test1() {
	fstream os;
	os.open("temp1.txt");
	//构造三对角型的系数矩阵A,装在C中
	static com C[2 * N0 - 1][2];
	constexpr int n = 2 * N0 - 1;
	C[0][0] = 0;
	double x = -N0 * delta_x;
	for (int i = 1; i < n; i++) {//构造非对角元
		C[i][0] = -0.5 / pow(delta_x, 2);
	}
	for (int i = 0; i < n; i++) {//构造对角元
		x += delta_x;
		C[i][1] = 1 / pow(delta_x, 2) - 1 / sqrt(2 + pow(x, 2))+0.48;
	}
	//u 为初始向量
	com* u = new com[n];
	for (int i = 0; i < n; i++) {
		u[i] = 1 / sqrt(2 * N0*delta_x);
	}

	eigen psai0 = inver_exp<N0>(C, u);
	normalize<n>(psai0.vec);
	x = -N0 * delta_x;
	for (int i = 0; i < n; i++) {
		x += delta_x;
		os << x << "  " << pow(abs(psai0.vec[i]),2)<<endl;
	}
	cout << psai0.val << endl;
	delete[] u;
}
//

void testexp() {
	constexpr int n = 10000000;
	com* b = new com[n];
	b[0] = 60;
	b[n - 1] = 60;
	for (int i = 1; i < n - 1; i++) {
		b[i] = 120;
	}
	static com C[n][2];
	C[0][0] = 0;
	C[0][1] = 5;
	C[n - 1][1] = 5;
	for (int i = 1; i <= n - 2; i++) {
		C[i][1] = 6;
	}
	for (int i = 1; i <= n - 1; i++) {
		C[i][0] = 4;
	}


	com* x = solve_eq<n>(C, b);
	print_sol<n>(x);
}
int main() {
	generate_psaif();
	system("pause");
}