#ifndef MPW_H
#define MPW_H

#include <iostream>
#include <string>
#include <cmath>
#include <boost/variant.hpp>
#include <tuple>
#include "mpfr.h"
#include "mpc.h"

namespace mpw_ns
{
	class mpw_defs
	{
		friend class mprw;
		friend class mpcw;

		private:
			static mpfr_prec_t mpw_b_prec;
			static size_t mpw_d_prec;
			static mpfr_rnd_t mprw_rnd_type;
			static mpc_rnd_t mpcw_rnd_type;
			mpw_defs() {};

		public:
			static void set_b_prec(const mpfr_prec_t new_prec);
			static void set_d_prec(const size_t new_prec);
			static void set_rnd_type(const mpfr_rnd_t new_rnd);
			mpfr_prec_t b_prec() const;
			size_t d_prec() const;
			mpfr_rnd_t rnd_type();
	};
	mpfr_prec_t mpw_defs::mpw_b_prec = 56;
	mpfr_rnd_t mpw_defs::mprw_rnd_type = MPFR_RNDN;
	mpc_rnd_t mpw_defs::mpcw_rnd_type = MPC_RNDNN;

	class mprw
	{
		friend class mpcw;

		private:
			mpfr_t mpr_l;

		public:
			mprw();
			mprw(const int in);
			mprw(const double in);
			mprw(const mprw &in);
			~mprw();
			std::string str(const size_t length = 10) const;
			bool isnan();
			bool operator>(const mprw &other) const;
			bool operator<(const mprw &other) const;
			bool operator>=(const mprw &other) const;
			bool operator<=(const mprw &other) const;
			bool operator==(const mprw &other) const;
			mprw operator+(const mprw &other) const;
			mprw operator-(const mprw &other) const;
			mprw operator-() const;
			mprw operator*(const mprw &other) const;
			mprw operator/(const mprw &other) const;
			mprw abs() const;
			mprw real() const;
			mprw imag() const;
			mprw conj() const;
			mprw sqrt() const;
			mprw exp() const;
			mprw sin() const;
			mprw cos() const;
			mprw tan() const;
			mprw log() const;
			mprw log10() const;
			mprw sinh() const;
			mprw cosh() const;
			mprw tanh() const;
			mprw asin() const;
			mprw acos() const;
			mprw atan() const;
			mprw asinh() const;
			mprw acosh() const;
			mprw atanh() const;

	};

	class mpcw
	{
		friend class mprw;

		private:
			mpc_t mpc_l;

		public:
			mpcw();
			mpcw(const mprw &re_in);
			mpcw(const int re_in);
			mpcw(const int re_in, const int im_in);
			mpcw(const double re_in);
			mpcw(const double re_in, const double im_in);
			mpcw(const mprw &re_in, const mprw &im_in);
			mpcw(const mpcw &in);
			~mpcw();
			std::string str(const size_t length) const;
			bool isnan();
			mpcw operator+(const mpcw &other) const;
			mpcw operator-(const mpcw &other) const;
			mpcw operator-() const;
			mpcw operator*(const mpcw &other) const;
			mpcw operator/(const mpcw &other) const;
			bool operator==(const mpcw &other) const;
			mprw abs() const;
			mprw real() const;
			mprw imag() const;
			mpcw conj() const;
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

	template <typename _T> class mpbm
	{
		protected:
			_T *data;
			size_t n_rows;
			size_t n_cols;
			size_t _index(const size_t r, const size_t c) const;

		public:
			mpbm();
			mpbm(const size_t rows, const size_t cols);
			mpbm(const size_t rows, const size_t cols, const int **ext_data);
			mpbm(const size_t rows, const size_t cols, const double **ext_data);
			mpbm(const size_t n_rows, const size_t n_cols);
			mpbm(const size_t n_rows, const size_t n_cols, const int **ext_data);
			mpbm(const size_t n_rows, const size_t n_cols, const double **ext_data);
			mpbm(const size_t rows, const size_t cols, const _T **ext_data);
			mpbm(const mpbm &other);
			~mpbm();
			std::string& str(const size_t length) const;
			size_t N_rows() const;
			size_t N_cols() const;
			const _T &operator()(const size_t r, const size_t c) const;
			_T &operator()(const size_t r, const size_t c);
			bool is_Symmetric() const;
			bool is_Hermitian() const;
			bool is_Square() const;
			bool is_row_vector() const;
			bool is_col_vector() const;
	};

	class mprm: public mpbm<mprw>
	{
		public:
			mprm();
			mprm(const size_t n_rows, const size_t n_cols);
			mprm(const size_t n_rows, const size_t n_cols, const int **ext_data);
			mprm(const size_t n_rows, const size_t n_cols, const double **ext_data);
			mprm(const size_t n_rows, const size_t n_cols, const mprw **ext_data);
			mprm operator+(const mprm& other) const;
			mprm operator+(const int other) const;
			mprm operator+(const double other) const;
			mprm operator+(const mprw& other) const;
			mprm operator-() const;
			mprm operator-(const mprm& other) const;
			mprm operator-(const int other) const;
			mprm operator-(const double other) const;
			mprm operator-(const mprw& other) const;
			mprm operator*(const mprm& other) const;
			mprm operator*(const int other) const;
			mprm operator*(const double other) const;
			mprm operator*(const mprw& other) const;
			mprm operator/(const mprm& other) const;
			mprm operator/(const int other) const;
			mprm operator/(const double other) const;
			mprm operator/(const mprw& other) const;
			mpbm T() const;
			mpbm H() const;
			mpbm conj() const;
			mpbm sqrt() const;
			mpbm exp() const;
			mpbm sin() const;
			mpbm cos() const;
			mpbm tan() const;
			mpbm log() const;
			mpbm log10() const;
			mpbm sinh() const;
			mpbm cosh() const;
			mpbm tanh() const;
			mpbm asin() const;
			mpbm acos() const;
			mpbm atan() const;
			mpbm asinh() const;
			mpbm acosh() const;
			mpbm atanh() const;
			std::tuple<mpbm<_T>, mpbm<_T>> Hessenberg();
			std::tuple<mpbm<_T>, mpbm<_T>> eig();
	};

	class mpcm: public mpbm<mpcw>
	{
		public:
			mpcm();
			mpcm(const size_t n_rows, const size_t n_cols);
			mpcm(const size_t n_rows, const size_t n_cols, const int **ext_data);
			mpcm(const size_t n_rows, const size_t n_cols, const double **ext_data);
			mpcm(const size_t n_rows, const size_t n_cols, const mprw **ext_data);
			mpcm(const size_t n_rows, const size_t n_cols, const mpcw **ext_data);
			mpcm operator+(const mpcm& other) const;
			mpcm operator+(const mprm& other) const;
			mpcm operator+(const int other) const;
			mpcm operator+(const double other) const;
			mpcm operator+(const mprw& other) const;
			mpcm operator-() const;
			mpcm operator-(const mpcm& other) const;
			mpcm operator-(const mprm& other) const;
			mpcm operator-(const int other) const;
			mpcm operator-(const double other) const;
			mpcm operator-(const mprw& other) const;
			mpcm operator*(const mprm& other) const;
			mpcm operator*(const int other) const;
			mpcm operator*(const double other) const;
			mpcm operator*(const mprw& other) const;
			mpcm operator/(const mprm& other) const;
			mpcm operator/(const int other) const;
			mpcm operator/(const double other) const;
			mpcm operator/(const mprw& other) const;
	};

	mprw operator+(const int in1, const mprw &in2);
	mprw operator+(const double in1, const mprw &in2);
	mprw operator-(const int in1, const mprw &in2);
	mprw operator-(const double in1, const mprw &in2);
	mprw operator*(const int in1, const mprw &in2);
	mprw operator*(const double in1, const mprw &in2);
	mprw operator/(const int in1, const mprw &in2);
	mprw operator/(const double in1, const mprw &in2);

	mpcw operator+(const int in1, const mpcw &in2);
	mpcw operator+(const double in1, const mpcw &in2);
	mpcw operator+(const mprw &in1, const mpcw &in2);
	mpcw operator-(const int in1, const mpcw &in2);
	mpcw operator-(const double in1, const mpcw &in2);
	mpcw operator-(const mprw &in1, const mpcw &in2);
	mpcw operator*(const int in1, const mpcw &in2);
	mpcw operator*(const double in1, const mpcw &in2);
	mpcw operator*(const mprw &in1, const mpcw &in2);
	mpcw operator/(const int in1, const mpcw &in2);
	mpcw operator/(const double in1, const mpcw &in2);
	mpcw operator/(const mprw &in1, const mpcw &in2);

	template <typename _T> std::string str(const _T& in, const size_t length);
	inline mprw abs(const mprw& in){return in.abs();};
	inline mprw abs(const mpcw &in){return in.abs();};
	inline mprm abs(const mprm& in){return in.abs();};
	mprm abs(const mpcm& in);
	inline mprw real(const mprw& in){return in.real();};
	inline mprw real(const mpcw& in){return in.real();};
	mprm real(const mprm& in);
	mprm real(const mpcm& in);
	inline mprw imag(const mprw& in){return in.imag();};
	inline mprw imag(const mpcw& in){return in.imag();};
	mprm imag(const mprm& in);
	mprm imag(const mpcm& in);
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

	mpbm<mprw> zeros(const size_t n_rows, const size_t n_cols);
	mpbm<mprw> ones(const size_t n_rows, const size_t n_cols);
	mpbm<mprw> I(const size_t n_rows, const size_t n_cols);

	mprw rand();
	mpcw crand();
	mprm rand(const size_t n_rows, const size_t n_cols);
	mpcm crand(const size_t n_rows, const size_t n_cols);

}

#endif
