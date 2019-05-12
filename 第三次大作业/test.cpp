#include<algorithm>
#include<iostream>
#include<fstream>
#include<iostream>
#include<complex>
#include"fft.h"
#include"const.h"
using namespace std;
using com = complex<double>;
//x轴离散成2*N0个格子
//constexpr int N0 = 16384;
//x间隔为delta_x
//constexpr double delta_x = 0.1;
//归一化
template<int n>
void normalize(com* x,double delta_x) {
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
void absorb(com A[n],double delta_x) {
	int N0 = (n + 1) / 2;
	double x0 = 0.75*N0*delta_x;
	double x = -N0 * delta_x;
	for (int i = 0; i < n; i++) {
		x += delta_x;
		if (abs(x) > x0) {
			com temp = pow(E, -pow((abs(x) - x0) / 0.2, 2));

			com temp2 = A[i];


			A[i] = temp2 * temp;
		}

	}
}


//波函数传播delta_t
//其中I(10^16为单位),omega(原子单位),N是激光场参数
template<int N0>
void Udelta_t(com psai[2*N0-1], double t, double I, double omega, int N, double delta_t,double delta_x) {
	constexpr int n = 2 * N0 - 1;
	//构造哈密顿量H,装在C中
	com(*C)[2] = new com[n][2];
	//constexpr int n = 2 * N0 - 1;
	C[0][0] = 0;
	double x = -N0 * delta_x;//x是坐标
	for (int i = 1; i < n; i++) {//构造非对角元
		C[i][0] = -0.5 / pow(delta_x, 2);
	}
	//带入原子单位换算公式求sqrt(I)
	const double sqrtI = sqrt(I / 3.5094448314);
	for (int i = 0; i < n; i++) {//构造对角元
		x += delta_x;
		C[i][1] = x * sqrtI*pow(sin(omega*t / 2. / N), 2)*sin(omega*t)
			+ 1 / pow(delta_x, 2) - 1 / sqrt(2 + pow(x, 2));
	}
	for (int i = 0; i < n; i++) {//*0.5*i*delta_t
		com i0{ 0.,1. };
		C[i][1] = C[i][1] * (-0.5*delta_t*i0);
		C[i][0] = C[i][0] * (-0.5*delta_t*i0);
	}


	//计算psai*(1-0.5*i*H*delta_t)
	com *psai_temp = new com[n];//临时储存psai
	for (int i = 0; i < n; i++) {
		psai_temp[i] = psai[i];
	}
	for (int k = 1; k < n; k++) {//开始计算
		psai[k] = psai_temp[k] * (1. - C[k][1]) - psai_temp[k - 1] * C[k][0] - psai_temp[k + 1] * C[k + 1][0];
	}
	psai[0] = psai_temp[0] * (1. - C[0][1]) - psai_temp[1] * C[1][0];//计算端点值
	psai[n - 1] = psai_temp[n - 1] * (1. - C[n - 1][1]) - psai_temp[n - 2] * C[n - 1][0];
	//计算psai*(1+0.5*i*H*delta_t)^-1
	for (int i = 0; i < n; i++) {//临时储存psai
		psai_temp[i] = psai[i];
	}
	for (int i = 0; i < n; i++) {//构造系数矩阵(1+0.5*i*H*delta_t)
		C[i][1] = 1. + C[i][1];
	}

	/*
	delete[] psai;
	psai = nullptr;
	psai = solve_eq<n>(C, psai_temp);//开始计算
	*/
	com* temp_solve = solve_eq<n>(C, psai_temp);
	for (int i = 0; i < n; i++) {
		psai[i] = temp_solve[i];
	}
	delete[] temp_solve;
	
	delete[] C;
	delete[] psai_temp;
}


//生成t=0的波函数
template<int N0>
com* generate_t0(double delta_x) {
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
	normalize<n>(psai0.vec,delta_x);
	delete[] u,C;
	return psai0.vec;
}

//求末态电离波函数(针对第(2)问)并且输出每一时刻的布居数
com* generate_psaif() {

	//x轴离散成2*N0个格子
	constexpr int N0 = 16384;
	//x间隔为delta_x
	constexpr double delta_x = 0.1;
	fstream os;
	os.open( "temp2_Pt.txt" );
	//初始化激光数据
	double I = 1.;
	double omega = 1.;
	int N = 18;
	double tf = 2 * N*PI / omega;
	//t迭代步长delta_t
	const double delta_t = 0.05;
	//初态波函数
	com* psai0 = generate_t0<N0>(delta_x);
	//求末态波函数,并输出每一时刻的布居数
	//!!!!!!!!com* psaif = generate_tf(tf, I, omega, N);
	com* psaif = generate_t0<N0>(delta_x);
	for (double t = 0; t < tf + delta_t; t += delta_t) {
		Udelta_t<N0>(psaif, t, I, omega, N, delta_t,delta_x);//波函数传播
		absorb<2 * N0 - 1>(psaif,delta_x);//乘吸收函数
		//输出时间&布居数
		com p = 0;
		for (int i = 0; i < 2 * N0 - 1; i++) {
			p += conj(psai0[i])*psaif[i] * delta_x;
		}
		os <<t<<"\t"<< abs(p)*abs(p)<<endl;
		cout << t << endl;//!!!!!!!!!
	}
	//基态概率幅p
	com p=0;
	for (int i = 0; i < 2 * N0 - 1; i++) {
		p += conj(psai0[i])*psaif[i]*delta_x;
	}
	cout << abs(p)*abs(p);//!!!!!!!!
	//减去基态部分
	for (int i = 0; i < 2 * N0 - 1; i++) {
		psaif[i] = psaif[i] - p * psai0[i];
	}
	delete[] psai0;
	return psaif;
}

//求动量谱(针对第(2)问)
void momentum_distribute() {
	//x轴离散成2*N0个格子
	constexpr int N0 = 16384;
	//x间隔为delta_x
	constexpr double delta_x = 0.1;
	fstream osk;
	osk.open("temp2_k.txt");
	//输出每一时刻的布居数
	//并且输出末态电离函数
	com* psaif = generate_psaif();
	//k的步长
	const double delta_k = 0.05;
	//输出动量谱
	com pk = 0.;
	double x = -N0 * delta_x;
	for (double k = -2.5; k < 2.5; k += delta_k) {
		//动量概率幅
		pk = 0.;
		x = -N0 * delta_x;
		for (int i = 0; i < 2 * N0 - 1; i++) {
			x += delta_x;
			pk += psaif[i] * pow(E, -com{ 0.,1. }*k*x)*delta_x;
		}
		osk << k << "\t" << abs(pk)*abs(pk) << endl;
	}
	delete[] psaif;
}
//第2问测试函数
void test2() {
	momentum_distribute();
}

//第1问测试函数
void test1() {
	//x轴离散成2*N0个格子
	constexpr int N0 = 16384;
	//x间隔为delta_x
	constexpr double delta_x = 0.1;
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
	normalize<n>(psai0.vec,delta_x);
	x = -N0 * delta_x;
	for (int i = 0; i < n; i++) {
		x += delta_x;
		os << x << "  " << pow(abs(psai0.vec[i]),2)<<endl;
	}
	cout << psai0.val << endl;
	delete[] u;
}

//求电偶极矩和t的函数
template<int N0>
com* generate_dt(int Ndata,double delta_x) {
	//初始化激光数据
	double I = 0.02;
	double omega300 = 45.5633525316/300;
	int N = 48;
	double tf = 2 * N*PI / omega300;
	//对应的delta_t
	double delta_t = tf / Ndata;
	com* dt = new com[Ndata];	
	com* psai = generate_t0<N0>(delta_x);
	double t = 0;
	for (int i = 0; i < Ndata;i++) {
		Udelta_t<N0>(psai, t, I, omega300, N, delta_t,delta_x);//波函数传播
		absorb<2 * N0 - 1>(psai,delta_x);//乘吸收函数
		//求t时刻dt
		double x = -N0 * delta_x;
		for (int j = 0; j < 2 * N0 - 1; j++) {
			x += delta_x;
			dt[i] += conj(psai[j])*psai[j]*x * delta_x;
		}
		t += delta_t;
		if (i % 100 == 0) {//!!!!!!!!!!!!!
			cout << i << endl;
		}
	}
	//!!!!!!!!!!!
	fstream otemp;
	otemp.open("temp.txt");
	double t_temp = 0;
	for (int i = 0; i < Ndata; i++) {
		otemp <<t_temp<<"\t"<< abs(dt[i]) << endl;
		t_temp += delta_t;
	}
	delete[] psai;
	return dt;
}

//求加速度和t的函数
template<int N0>
com* generate_at(int Ndata,double delta_x,double I,double omega300,int N) {
	constexpr int n = 2 * N0 - 1;
	double tf = 2 * N*PI / omega300;
	double sqrtI = sqrt(I / 3.5094448314);
	//对应的delta_t
	double delta_t = tf / Ndata;
	com* at = new com[Ndata];
	com* psai = generate_t0<N0>(delta_x);
	double t = 0;
	for (int i = 0; i < Ndata; i++) {
		Udelta_t<N0>(psai, t, I, omega300, N, delta_t, delta_x);//波函数传播
		absorb<n>(psai,delta_x);//乘吸收函数
		//求t时刻at
		double x = -N0 * delta_x;
		for (int j = 0; j < 2 * N0 - 1; j++) {
			x += delta_x;
			//受力算符 -dV/dx+E
			com F = -x / pow(x*x + 2, 1.5) + sqrtI * pow(sin(omega300*t / 2. / N), 2)*sin(omega300*t);
			at[i] += conj(psai[j])*psai[j] *F* delta_x;
		}
		t += delta_t;
		if (i % 100 == 0) {//!!!!!!!!!!!!!
			cout << i << endl;
		}
	}
	delete[] psai;
	return at;
}

//利用电偶极矩输出辐射谱
void A_by_d() {
	//x轴离散成2*N0个格子
	constexpr int N0 = 16384;
	//x间隔为delta_x
	constexpr double delta_x = 0.1;
	fstream os;
	os.open("temp3_d.txt");
	//取4096个数据点
	const int Ndata = 4096;
	//激光数据
	int N = 48;
	double omega300 = 45.5633525316 / 300;
	double tf = 2 * N*PI / omega300;
	//求辐射谱
	com* dt = generate_dt<N0>(Ndata,delta_x);
	com* A=new com[Ndata];
	fft<Ndata>(dt, A);//做fft
	//乘上系数并输出
	double omega = 0;
	double delta_omega = 2 * PI / tf;
	for (int i = 0; i < Ndata; i++) {
		A[i]=-1 / sqrt(2 * PI)*pow(omega, 2)*A[i];
		//输出omega/omega300 , log10A
		os << omega/omega300 <<"\t"<< log10(norm(A[i]))<<endl;
		omega += delta_omega;
	}

	delete[] dt,A;
}

//利用加速度输出辐射谱
void A_by_a() {
	//初始化激光数据
	double I = 0.02;
	double omega300 = 45.5633525316 / 300;
	int N = 48;
	double tf = 2 * N*PI / omega300;
	//x轴离散成2*N0个格子
	constexpr int N0 = 16384;
	//x间隔为delta_x
	constexpr double delta_x = 0.1;
	fstream os;
	os.open("temp3_a.txt");//!!!!!!!!!!!
	//取4096个数据点
	const int Ndata = 4096;
	//求辐射谱
	com* at = generate_at<N0>(Ndata,delta_x,I,omega300,N);
	com* A = new com[Ndata];
	fft<Ndata>(at, A);//做fft
	//乘上系数并输出
	double omega = 0;
	double delta_omega = 2 * PI / tf;
	for (int i = 0; i < Ndata; i++) {
		A[i] = 1 / sqrt(2 * PI)*A[i];
		//输出omega/omega300 , log10A
		os << omega / omega300 << "\t" << log10(norm(A[i])) << endl;
		omega += delta_omega;
	}

	delete[] at, A;
}

//第3问测试函数
void test3() {
	//运行A_by_d();或A_by_a();
	//A_by_d();
	A_by_a();
	//不要同时运行这两个函数,否则会出错
}
//第四问测试函数
void test4() {
	ofstream os4;
	os4.open( "temp4.txt" );//!!!!!!!!!!!!!!!!!!!!
	//初始化激光数据
	double I = 0.01;
	double omega400 = 45.5633525316 / 400;
	int N = 4;
	double tf = 2 * N*PI / omega400;
	//迭代数据
	constexpr int N0 = 16384;//取2N0-1个坐标点
	double delta_x = 0.1;
	int Ndata = 400;//取800个时间点!!!!!!!!!!!!!!!!!!!!!
	double delta_t = tf / Ndata;
	//输出的omega步长
	double delta_omega = omega400 / 20;
	//生成加速度
	com* at = generate_at<N0>(Ndata, delta_x, I, omega400, N);
	//输出功率谱
	for (double omega = 0; omega < 20 * omega400; omega += delta_omega) {
		for (double t0 = 0; t0 <= tf; t0 += delta_t) {
			com A = 0;
			double t = 0;
			for (int i=0; i<Ndata;i++, t += delta_t) {
				A += at[i] * pow(E, -com{ 0.,1. }*omega*t - pow((t - t0), 2) / 2 / 15 / 15)*delta_t;
			}
			os4 << t0 << "\t" << omega / omega400 << "\t" << log10(norm(A))<<endl;
		}
	}
}

//第四问测试函数
void test4_usefft() {
	ofstream os4;
	os4.open("temp.txt");//!!!!!!!!!!!!!!!!!!!!
	//初始化激光数据
	double I = 0.01;
	double omega400 = 45.5633525316 / 400;
	int N = 4;
	double tf = 2 * N*PI / omega400;
	//迭代数据
	constexpr int N0 = 16384;//取2N0-1个坐标点
	double delta_x = 0.1;
	constexpr int Ndata = 512;//取800个时间点!!!!!!!!!!!!!!!!!!!!!
	double delta_t = tf / Ndata;
	//输出的omega步长
	double delta_omega = omega400 / 20;
	//生成加速度
	com* at = generate_at<N0>(Ndata, delta_x, I, omega400, N);
	//输出功率谱
	/*
	for (double omega = 0; omega < 20 * omega400; omega += delta_omega) {
		

		for (double t0 = 0; t0 <= tf; t0 += delta_t) {
			com A = 0;
			double t = 0;
			for (int i = 0; i < Ndata; i++, t += delta_t) {
				A += at[i] * pow(E, -com{ 0.,1. }*omega*t - pow((t - t0), 2) / 2 / 15 / 15)*delta_t;
			}
			os4 << t0 << "\t" << omega / omega400 << "\t" << log10(norm(A)) << endl;
		}
	}*/
	for (int t0 = 0; t0 <= tf; t0 += tf / 20) {
		com* a_temp = new com[Ndata];
		com* result = new com[Ndata];
		double t = 0;
		for (int i = 0; i < Ndata; i++) {
			a_temp[i] = at[i] * pow(E, -pow((t - t0), 2) / 2. / 15. / 15.);
			t += delta_t;
		}
		fft<Ndata>(a_temp, result);
		int k = 0;
		for (double omega = 0; omega < 20 * omega400; omega += 2 * PI / tf) {
			os4 << t0 << "\t" << omega / omega400 << "\t" << log10(norm(result[k])) << endl;
			k++;
		}
		delete[] a_temp, result;
		cout << t0 / tf * 20 << endl;//!!!!!!!!!
	}
	delete[] at;
}



int main() {
	//testbug3();
	test4_usefft();
	system("pause");
	;
	;
	;
	;//测试3
}


