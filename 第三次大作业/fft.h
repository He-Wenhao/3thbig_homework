#pragma once
#include<complex>
#include"const.h"
using com = std::complex<double>;
//FFT

//递归
template<int n>//n必须为2的幂次
void fft_origin(com x[n], com y[n], com omega) {
	if (n == 1) {
		y[0] = x[0];//底部循环
	}
	else {
		//检查n是否是2的幂次
		if (n % 2 == 1) {
			throw std::runtime_error{ "not 2^i" };
		}
		com* p = new com[n / 2];
		com* s = new com[n / 2];
		for (int k = 0; k < n / 2; k++) {
			p[k] = x[2 * k];
			s[k] = x[2 * k + 1];
		}//分成奇偶子列
		com* q = new com[n / 2];
		com* t = new com[n / 2];
		//递归调用fft
		constexpr int k0 = n / 2;
		fft_origin<k0>(p, q, pow(omega, 2.));
		fft_origin<k0>(s, t, pow(omega, 2.));
		for (int k = 0; k < n; k++) {
			y[k] = q[k % (n / 2)] + pow(omega, k)*t[k % (n / 2)];
		}//组合
		delete[] p, s, q, t;
	}
}

//底层调用
template<>
void fft_origin<1>(com x[1], com y[1], com omega) {
	y[0] = x[0];
}


//真正调用的fft
template<int n>//n必须为2的幂次
void fft(com x[n], com y[n]) {
	com omega = pow(E, -com{ 0,1 }*2. * PI / double(n));
	fft_origin<n>(x, y, omega);
}


//真正调用的逆变换
template<int n>//n必须为2的幂次
void fft_inv(com x[n], com y[n]) {
	com omega = pow(E, +com{ 0,1 }*2. * PI / double(n));
	fft_origin<n>(x, y, omega);
	for (int i = 0; i < n; i++) {
		y[i] = y[i] / double(n);
	}
}
