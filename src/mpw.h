#ifndef MPW_H
#define MPW_H

#include <iostream>
#include <string>
#include <cmath>
#include <boost/variant.hpp>
#include <tuple>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include "mpfr.h"
#include "mpc.h"

namespace mpw_ns
{
	class mpw_defs
	{
		friend class mpcw;

		private:
			static mpfr_prec_t mpw_b_prec;
			static size_t mpw_b_prec_overshoot;
			static size_t mpw_d_prec;
			static mpfr_rnd_t mpfr_rnd_type;
			static mpc_rnd_t mpcw_rnd_type;
			static size_t str_length;
			static size_t mpw_max_iters;
			static size_t mpw_d_prec_buffer;
			static int mpw_tol_exp;
			static bool mpw_use_Householder_prec;
			mpw_defs() {};

		public:
			static void set_d_prec(const size_t new_prec);
			static void set_rnd_type(const mpfr_rnd_t new_rnd);
			static void set_mpcw_str_length(const size_t length);
			static void set_max_iters(const size_t new_iters);
			static void set_d_prec_buffer(const size_t new_buffer);
			static void set_tol_exp(const int new_tol);
			static void set_use_Householder_prec(const bool new_prec);
			static mpfr_prec_t b_prec();
			static size_t d_prec();
			static size_t max_iters();
			static size_t d_prec_buffer();
			static mpfr_rnd_t rnd_type();
			static size_t mpcw_str_length();
			static int tol_exp();
			static bool use_Householder_prec();
			static std::string mpfr_to_str(const mpfr_t val, const size_t length);
	};


	class mpcw
	{
		private:
			mpc_t mpc_l;

		public:
			mpcw();
			mpcw(const int re_in);
			mpcw(const int re_in, const size_t new_d_prec);
			mpcw(const int re_in, const int im_in);
			mpcw(const int re_in, const int im_in, const size_t new_d_prec);
			mpcw(const double re_in);
			mpcw(const mpfr_t re_in);
			mpcw(const double re_in, const double im_in);
			mpcw(const mpfr_t re_in, const mpfr_t im_in);
			mpcw(const mpcw& in);
			mpcw(const mpcw& in, const size_t new_d_prec);
			mpcw(const mpcw& re_in, const mpcw& im_in);
			~mpcw();
			std::string str(const size_t length) const;
			bool isnan() const;
			bool is_real() const;
			bool is_imag() const;
			void operator=(const int other);
			void operator=(const double other);
			void operator=(const mpcw& other);
			mpcw operator+(const mpcw& other) const;
			mpcw operator+(const double other) const;
			mpcw operator+(const int other) const;
			mpcw& operator+=(const mpcw& other);
			mpcw& operator+=(const double other);
			mpcw& operator+=(const int other);
			mpcw operator-() const;
			mpcw operator-(const mpcw& other) const;
			mpcw operator-(const double other) const;
			mpcw operator-(const int other) const;
			mpcw& operator-=(const int other);
			mpcw& operator-=(const double other);
			mpcw& operator-=(const mpcw& other);
			mpcw operator*(const mpcw& other) const;
			mpcw operator*(const double other) const;
			mpcw operator*(const int other) const;
			mpcw& operator*=(const mpcw& other);
			mpcw& operator*=(const double other);
			mpcw& operator*=(const int other);
			mpcw operator/(const mpcw& other) const;
			mpcw operator/(const double other) const;
			mpcw operator/(const int other) const;
			mpcw& operator/=(const mpcw& other);
			mpcw& operator/=(const double other);
			mpcw& operator/=(const int other);
			mpcw operator^(const mpcw& other) const;
			mpcw operator^(const double other) const;
			mpcw operator^(const int other) const;
			bool operator==(const mpcw& other) const;
			bool operator==(const double other) const;
			bool operator==(const int other) const;
			bool operator>(const mpcw& other) const;
			bool operator>(const double other) const;
			bool operator>(const int other) const;
			bool operator>=(const mpcw& other) const;
			bool operator>=(const double other) const;
			bool operator>=(const int other) const;
			bool operator<(const mpcw& other) const;
			bool operator<(const double other) const;
			bool operator<(const int other) const;
			bool operator<=(const mpcw& other) const;
			bool operator<=(const double other) const;
			bool operator<=(const int other) const;
			bool operator!=(const mpcw& other) const;
			bool operator!=(const double other) const;
			bool operator!=(const int other) const;
			mpcw abs() const;
			mpcw real() const;
			mpcw imag() const;
			mpcw conj() const;
			mpcw angle() const;
			size_t length();
			mpcw T() const;
			mpcw H() const;
			mpcw sqrt() const;
			mpcw exp() const;
			mpcw sin() const;
			mpcw cos() const;
			mpcw tan() const;
			mpcw log() const;
			mpcw log10() const;
			mpcw sinh() const;
			mpcw cosh() const;
			mpcw tanh() const;
			mpcw asin() const;
			mpcw acos() const;
			mpcw atan() const;
			mpcw asinh() const;
			mpcw acosh() const;
			mpcw atanh() const;
	};


	class mpcm
	{
		private:
			mpcw* data;
			size_t n_rows;
			size_t n_cols;
			size_t _index(const size_t r, const size_t c) const;
			size_t _max_index() const;
			void _initialize_pivot(uint& r_p, uint& c_p, mpcm& max_vals, uint* max_indexes);
			void _update_pivot(uint& r_p, uint& c_p, mpcm& max_vals, uint* max_indexes);

		public:
			mpcm();
			mpcm(const size_t n_rows, const size_t n_cols);
			mpcm(const size_t n_rows, const size_t n_cols, const int** ext_data);
			mpcm(const size_t n_rows, const size_t n_cols, const double** ext_data);
			mpcm(const size_t n_rows, const size_t n_cols, const mpcw** ext_data);
			mpcm(const mpcm& other);
			mpcm(const mpcm& other, const size_t new_d_prec);
			~mpcm();
			std::string str(const size_t length) const;
			size_t N_rows() const;
			size_t N_cols() const;
			const mpcw& operator()(const size_t r, const size_t c) const;
			mpcw& operator()(const size_t r, const size_t c);
			bool is_Symmetric() const;
			bool is_Hermitian() const;
			bool is_Square() const;
			bool is_row_vector() const;
			bool is_col_vector() const;
			mpcm T() const;
			mpcm H() const;
			mpcm inv() const;
			mpcw max_diff(const mpcm& other) const;
			std::tuple<mpcm, mpcm> Householder_double_sided() const;
			std::tuple<mpcm, mpcm> Householder_double_sided(const size_t prec_buffer) const;
			std::tuple<mpcm, mpcm> Householder_single_sided() const;
			std::tuple<mpcm, mpcm> Householder_single_sided(const size_t prec_buffer) const;
			void Jacobi_complex_rotator(mpcm& J, const size_t r_p, const size_t c_p);
			void Jacobi_real_rotator(mpcm& J, const size_t r_p, const size_t c_p);
			void Jacobi_rotator(mpcm& J, const size_t r_p, const size_t c_p, const mpcw& x1, const mpcw& x2, const mpcw& x3, const mpcw& x4);
			std::tuple<mpcm, mpcm> Hermitian_eig() const;
			std::tuple<mpcm, mpcm> Hermitian_eig(const size_t max_iter, const mpcw& tol, const size_t prec_buffer, const bool use_Householder_prec) const;
			std::tuple<mpcm, mpcm> Takagi_fact(const mpcm& M, const size_t max_iter, const size_t prec_buffer) const;
			std::tuple<mpcm, mpcm, mpcm> SVD(const mpcm& M, const size_t max_iter, const size_t prec_buffer) const;
			std::tuple<mpcm, mpcm> QR(const mpcm& M, const size_t max_iter, const size_t prec_buffer) const;
			std::tuple<mpcm, mpcm> LQ(const mpcm& M, const size_t max_iter, const size_t prec_buffer) const;
			size_t length() const;
			mpcm diag() const;
			mpcw scalar() const;
			mpcm get_row(size_t row_index) const;
			void set_row(size_t row_index, const mpcm& new_row);
			mpcm get_column(size_t row_index) const;
			void set_column(size_t row_index, const mpcm& new_row);
			mpcm normalized() const;
			void operator=(const mpcm& other);
			mpcm operator+(const mpcm& other) const;
			mpcm operator+(const int other) const;
			mpcm operator+(const double other) const;
			mpcm operator+(const mpcw& other) const;
			mpcm operator-() const;
			mpcm operator-(const mpcm& other) const;
			mpcm operator-(const int other) const;
			mpcm operator-(const double other) const;
			mpcm operator-(const mpcw& other) const;
			mpcm operator*(const mpcm& other) const;
			mpcm operator*(const int other) const;
			mpcm operator*(const double other) const;
			mpcm operator*(const mpcw& other) const;
			mpcm operator/(const mpcm& other) const;
			mpcm operator/(const int other) const;
			mpcm operator/(const double other) const;
			mpcm operator/(const mpcw& other) const;
			mpcm abs() const;
			mpcm real() const;
			mpcm imag() const;
			mpcm conj() const;
			mpcm angle() const;
			mpcm sqrt() const;
			mpcm exp() const;
			mpcm sin() const;
			mpcm cos() const;
			mpcm tan() const;
			mpcm log() const;
			mpcm log10() const;
			mpcm sinh() const;
			mpcm cosh() const;
			mpcm tanh() const;
			mpcm asin() const;
			mpcm acos() const;
			mpcm atan() const;
			mpcm asinh() const;
			mpcm acosh() const;
			mpcm atanh() const;
	};


	mpcw operator+(const int in1, const mpcw& in2);
	mpcw operator+(const double in1, const mpcw& in2);
	mpcw operator-(const int in1, const mpcw& in2);
	mpcw operator-(const double in1, const mpcw& in2);
	mpcw operator*(const int in1, const mpcw& in2);
	mpcw operator*(const double in1, const mpcw& in2);
	mpcw operator/(const int in1, const mpcw& in2);
	mpcw operator/(const double in1, const mpcw& in2);
	mpcm operator+(const int in1, const mpcm& in2);
	mpcm operator+(const double in1, const mpcm& in2);
	mpcm operator+(const mpcw& in1, const mpcm& in2);
	mpcm operator-(const int in1, const mpcm& in2);
	mpcm operator-(const double in1, const mpcm& in2);
	mpcm operator-(const mpcw& in1, const mpcm& in2);
	mpcm operator*(const int in1, const mpcm& in2);
	mpcm operator*(const double in1, const mpcm &in2);
	mpcm operator*(const mpcw& in1, const mpcm &in2);
	mpcm operator/(const int in1, const mpcm& in2);
	mpcm operator/(const double in1, const mpcm &in2);
	mpcm operator/(const mpcw& in1, const mpcm &in2);

	template <typename _T> inline std::string str(const _T& in, const size_t length ){return in.str(length);};
	template <typename _T> inline std::string str(const _T& in){return in.str(mpw_defs::mpcw_str_length());};
	inline mpcw abs(const mpcw &in){return in.abs();};
	inline mpcm abs(const mpcm& in){return in.abs();};
	inline mpcw real(const mpcw& in){return in.real();};
	inline mpcm real(const mpcm& in){return in.real();};
	inline mpcw imag(const mpcw& in){return in.imag();};
	inline mpcm imag(const mpcm& in){return in.imag();};
	inline mpcw conj(const mpcw &in){return in.conj();};
	inline mpcm conj(const mpcm& in){return in.conj();};
	inline mpcw angle(const mpcw &in){return in.angle();};
	inline mpcm angle(const mpcm& in){return in.angle();};
	inline mpcw scalar(const mpcm& in){return in.scalar();};
	inline mpcw max_diff(const mpcm& in1, const mpcm& in2){return in1.max_diff(in2);};
	template <typename T> inline T conj(const T &in){return in.conj();};
	template <typename T> inline T sqrt(const T &in){return in.sqrt();};
	template <typename T> inline T exp(const T &in){return in.exp();};
	template <typename T> inline T sin(const T &in){return in.sin();};
	template <typename T> inline T cos(const T &in){return in.cos();};
	template <typename T> inline T tan(const T &in){return in.tan();};
	template <typename T> inline T log(const T &in){return in.log();};
	template <typename T> inline T log10(const T &in){return in.log10();};
	template <typename T> inline T sinh(const T &in){return in.sinh();};
	template <typename T> inline T cosh(const T &in){return in.cosh();};
	template <typename T> inline T tanh(const T &in){return in.tanh();};
	template <typename T> inline T asin(const T &in){return in.asin();};
	template <typename T> inline T acos(const T &in){return in.acos();};
	template <typename T> inline T atan(const T &in){return in.atan();};
	template <typename T> inline T asinh(const T &in){return in.asinh();};
	template <typename T> inline T acosh(const T &in){return in.acosh();};
	template <typename T> inline T atanh(const T &in){return in.atanh();};

	mpcm zeros(const size_t n_rows, const size_t n_cols);
	mpcm ones(const size_t n_rows, const size_t n_cols);
	mpcm I(const size_t n_rows, const size_t n_cols);
	mpcm I(const size_t n_rows, const size_t n_cols, const size_t new_d_prec);

	mpcw rand();
	mpcw crand();
	mpcm rand(const size_t n_rows, const size_t n_cols);
	mpcm crand(const size_t n_rows, const size_t n_cols);
	mpcm int_test(const size_t n_rows, const size_t n_cols);
}

#endif
