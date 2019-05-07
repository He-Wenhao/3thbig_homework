#include<algorithm>
#include<iostream>
using namespace std;
constexpr int N0 = 10000;
constexpr double delta_x = 0.01;
template<int n>
double* solve_eq(double c_[n][2], double b_[n]) {
	
	double c[n][2];
	double* b = new double[n];
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
	return b;
}

//反幂法求本征问题
template<int N>
double* inver_exp(double C[2*N-1][2]) {
	//初始化u,lambda
	constexpr int n = 2 * N - 1;
	double u[n];
	for (int i = 0; i < n; i++) {
		u[i] = 1 / sqrt(2 * N*delta_x);
		
	}
	double lambda = 0;
	double* v=new double[n];
	for (int k = 0; k < 10000; k++) {
		delete[] v;
		v = solve_eq<n>(C, u);
		//求v绝对值最大的分量max
		double max = 0;
		for (int i = 0; i < n; i++) {
			if (abs(max) < abs(v[i])) {
				max = v[i];
			}
		}
		//储存前一次本征值
		double lambda_before = lambda;
		lambda = 1 / max;
		//u=v*lambda
		for (int i = 0; i < n; i++) {
			u[i] = v[i]*lambda;
		}
		//检查精度是否合适
		if (abs(lambda_before - lambda)/abs(lambda) < 1e-10) {
			break;
		}
	}
	cout << lambda;
	return u;
}


void test1() {

	//构造三对角型的系数矩阵A,装在C中
	double C[2 * N0 - 1][2];
	constexpr int n = 2 * N0 - 1;
	C[0][0] = 0;
	double x = -N0 * delta_x;
	for (int i = 1; i < n; i++) {
		C[i][0] = -0.5 / pow(delta_x, 2);
	}
	for (int i = 0; i < n; i++) {
		x += delta_x;
		C[i][1] = 1 / pow(delta_x, 2) - 1 / sqrt(2 + pow(x, 2))+0.48;
	}

	inver_exp<N0>(C);
}
//打印向量
template<int n>
void print_sol(double x[n]) {
	for (int i = 0; i < n; i++) { cout << i + 1 << "   " << x[i] << endl; }
}


void testexp() {
	constexpr int n = 55;
	double* b = new double[n];
	b[0] = 60;
	b[n - 1] = 60;
	for (int i = 1; i < n - 1; i++) {
		b[i] = 120;
	}
	double C[n][2];
	C[0][0] = 0;
	C[0][1] = 5;
	C[n - 1][1] = 5;
	for (int i = 1; i <= n - 2; i++) {
		C[i][1] = 6;
	}
	for (int i = 1; i <= n - 1; i++) {
		C[i][0] = 4;
	}


	double* x = solve_eq<n>(C, b);
	print_sol<n>(x);
}
int main() {
	test1();
	system("pause");
}