#pragma once
//一些辅助函数
#include<complex>
using namespace std;
//归一化
template<int n>
void normalize(complex<double>* x, double delta_x) {
	double norm = 0;
	for (int i = 0; i < n; i++) {
		norm += pow(abs(x[i]), 2)*delta_x;
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
//用来装打印向量
template<int n>
void print_sol(complex<double> x[n]) {
	for (int i = 0; i < n; i++) { cout << i + 1 << "   " << x[i] << endl; }
}
