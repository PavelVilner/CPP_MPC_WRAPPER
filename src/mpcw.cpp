//============================================================================
// Name        : MPCW.cpp
// Author      : Pavel Vilner
// Version     :
// Copyright   : Your copyright notice
// Description : Implementation of an MPC wrapper
//============================================================================

#include "mpw.h"

mpw_ns::mpcw::mpcw()
{
	mpc_init2(this->mpc_l, mpw_ns::mpw_defs::mpw_b_prec);
}

mpw_ns::mpcw::mpcw(const mpfr_t re_in)
{
	mpc_init2(this->mpc_l, mpw_ns::mpw_defs::mpw_b_prec);
	mpc_set_fr(this->mpc_l, re_in, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::mpcw(const int re_in)
{
	mpc_init2(this->mpc_l, mpw_ns::mpw_defs::mpw_b_prec);
	mpc_set_si(this->mpc_l, re_in, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::mpcw(const int re_in, const size_t new_d_prec)
{
	mpc_init2(this->mpc_l, mpfr_prec_t(ceil((new_d_prec+1) * 3.322) + 5));
	mpc_set_si(this->mpc_l, re_in, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::mpcw(const int re_in, const int im_in)
{
	mpc_init2(this->mpc_l, mpw_ns::mpw_defs::mpw_b_prec);
	mpc_set_si_si(this->mpc_l, re_in, im_in, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::mpcw(const int re_in, const int im_in, const size_t new_d_prec)
{
	mpc_init2(this->mpc_l, mpfr_prec_t(ceil((new_d_prec+1) * 3.322) + 5));
	mpc_set_si_si(this->mpc_l, re_in, im_in, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::mpcw(const double re_in)
{
	mpc_init2(this->mpc_l, mpw_ns::mpw_defs::mpw_b_prec);
	mpc_set_d(this->mpc_l, re_in, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::mpcw(const double re_in, const double im_in)
{
	mpc_init2(this->mpc_l, mpw_ns::mpw_defs::mpw_b_prec);
	mpc_set_d_d(this->mpc_l, re_in, im_in, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::mpcw(const mpfr_t re_in, const mpfr_t im_in)
{
	mpc_init2(this->mpc_l, mpw_ns::mpw_defs::mpw_b_prec);
	mpc_set_fr_fr(this->mpc_l, re_in, im_in, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::mpcw(const mpcw &in)
{
	mpc_init2(this->mpc_l, mpw_ns::mpw_defs::mpw_b_prec);
	mpc_set(this->mpc_l, in.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::mpcw(const mpcw &in, const size_t new_d_prec)
{
	mpc_init2(this->mpc_l, mpfr_prec_t(ceil((new_d_prec+1) * 3.322) + 5));
	mpc_set(this->mpc_l, in.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::~mpcw()
{
	mpc_clear(this->mpc_l);
}

std::string mpw_ns::mpcw::str(const size_t length = mpw_ns::mpw_defs::str_length) const
{
	std::string filler = std::string("");
	if (std::floor(length/2.0) == length/2.0)
	{
		filler = std::string(" ");
	}
	mpfr_t aux_re;
	mpfr_init(aux_re);
	mpfr_t aux_im;
	mpfr_init(aux_im);
	mpfr_t aux_re_abs;
	mpfr_init(aux_re_abs);
	mpfr_t aux_im_abs;
	mpfr_init(aux_im_abs);
	mpc_real(aux_re, this->mpc_l, mpw_ns::mpw_defs::mpfr_rnd_type);
	mpc_imag(aux_im, this->mpc_l, mpw_ns::mpw_defs::mpfr_rnd_type);
	mpfr_abs(aux_re_abs, aux_re, mpw_ns::mpw_defs::mpfr_rnd_type);
	mpfr_abs(aux_im_abs, aux_im, mpw_ns::mpw_defs::mpfr_rnd_type);
	std::string sign_re = std::string(" ");
	std::string sign_im = std::string(" + j");
	if (mpfr_cmp_si(aux_re, 0) < 0)
	{
		sign_re = std::string("-");
	}
	if (mpfr_cmp_si(aux_im, 0) < 0)
	{
		sign_im = std::string(" - j");
	}
	std::string abs_re_s = mpw_ns::mpw_defs::mpfr_to_str(aux_re_abs, size_t((length-5)/2));
	std::string abs_im_s = mpw_ns::mpw_defs::mpfr_to_str(aux_im_abs, size_t((length-5)/2));
	mpfr_clear(aux_re);
	mpfr_clear(aux_im);
	mpfr_clear(aux_re_abs);
	mpfr_clear(aux_im_abs);
	return sign_re + abs_re_s + sign_im + abs_im_s + filler;
}

bool mpw_ns::mpcw::isnan() const
{
	return this->real().isnan() || this->imag().isnan();
}

bool mpw_ns::mpcw::is_real() const
{
	return this->imag()== mpw_ns::mpcw(0);
}

bool mpw_ns::mpcw::is_imag() const
{
	return this->real()== mpw_ns::mpcw(0);
}

void mpw_ns::mpcw::operator=(const int other)
{
	mpc_set_si(this->mpc_l, other, mpw_ns::mpw_defs::mpcw_rnd_type);
}

void mpw_ns::mpcw::operator=(const double other)
{
	mpc_set_d(this->mpc_l, other, mpw_ns::mpw_defs::mpcw_rnd_type);
}

void mpw_ns::mpcw::operator=(const mpw_ns::mpcw& other)
{
	mpc_set(this->mpc_l, other.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw mpw_ns::mpcw::operator+(const int other) const
{
	return *this + mpw_ns::mpcw(other);
}

mpw_ns::mpcw mpw_ns::mpcw::operator+(const double other) const
{
	return *this + mpw_ns::mpcw(other);
}

mpw_ns::mpcw mpw_ns::mpcw::operator+(const mpcw& other) const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_add(res.mpc_l, this->mpc_l, other.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw& mpw_ns::mpcw::operator+=(const mpcw& other)
{
	mpc_add(this->mpc_l, this->mpc_l, other.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return *this;
}

mpw_ns::mpcw& mpw_ns::mpcw::operator+=(const int other)
{
	return *this += mpw_ns::mpcw(other);
}

mpw_ns::mpcw& mpw_ns::mpcw::operator+=(const double other)
{
	return *this += mpw_ns::mpcw(other);
}

mpw_ns::mpcw mpw_ns::mpcw::operator-(const int other) const
{
	return *this - mpw_ns::mpcw(other);
}

mpw_ns::mpcw mpw_ns::mpcw::operator-(const double other) const
{
	return *this - mpw_ns::mpcw(other);
}

mpw_ns::mpcw mpw_ns::mpcw::operator-(const mpcw& other) const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_sub(res.mpc_l, this->mpc_l, other.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw& mpw_ns::mpcw::operator-=(const int other)
{
	return *this -= mpw_ns::mpcw(other);
}

mpw_ns::mpcw& mpw_ns::mpcw::operator-=(const double other)
{
	return *this -= mpw_ns::mpcw(other);
}

mpw_ns::mpcw& mpw_ns::mpcw::operator-=(const mpcw& other)
{
	mpc_sub(this->mpc_l, this->mpc_l, other.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return *this;
}

mpw_ns::mpcw mpw_ns::mpcw::operator-() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_sub(res.mpc_l, mpw_ns::mpcw(0,0).mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::operator*(const int other) const
{
	return *this * mpw_ns::mpcw(other);
}

mpw_ns::mpcw mpw_ns::mpcw::operator*(const double other) const
{
	return *this * mpw_ns::mpcw(other);
}

mpw_ns::mpcw mpw_ns::mpcw::operator*(const mpcw& other) const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_mul(res.mpc_l, this->mpc_l, other.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw& mpw_ns::mpcw::operator*=(const int other)
{
	return *this *= mpw_ns::mpcw(other);
}

mpw_ns::mpcw& mpw_ns::mpcw::operator*=(const double other)
{
	return *this *= mpw_ns::mpcw(other);
}

mpw_ns::mpcw& mpw_ns::mpcw::operator*=(const mpcw& other)
{
	mpc_mul(this->mpc_l, this->mpc_l, other.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return *this;
}

mpw_ns::mpcw mpw_ns::mpcw::operator/(const int other) const
{
	return *this / mpw_ns::mpcw(other);
}

mpw_ns::mpcw mpw_ns::mpcw::operator/(const double other) const
{
	return *this / mpw_ns::mpcw(other);
}

mpw_ns::mpcw mpw_ns::mpcw::operator/(const mpcw& other) const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_div(res.mpc_l, this->mpc_l, other.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw& mpw_ns::mpcw::operator/=(const int other)
{
	return *this /= mpw_ns::mpcw(other);
}

mpw_ns::mpcw& mpw_ns::mpcw::operator/=(const double other)
{
	return *this /= mpw_ns::mpcw(other);
}

mpw_ns::mpcw& mpw_ns::mpcw::operator/=(const mpcw& other)
{
	mpc_div(this->mpc_l, this->mpc_l, other.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return *this;
}

mpw_ns::mpcw mpw_ns::mpcw::operator^(const mpcw& other) const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_pow(res.mpc_l, this->mpc_l, other.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::operator^(const double other) const
{
	return *this ^ mpw_ns::mpcw(other);
}

mpw_ns::mpcw mpw_ns::mpcw::operator^(const int other) const
{
	return *this ^ mpw_ns::mpcw(other);
}

bool mpw_ns::mpcw::operator==(const int other) const
{
	return *this == mpw_ns::mpcw(other);
}

bool mpw_ns::mpcw::operator==(const double other) const
{
	return *this == mpw_ns::mpcw(other);
}

bool mpw_ns::mpcw::operator==(const mpcw &other) const
{
	return mpc_cmp(this->mpc_l, other.mpc_l) == 0;
}

bool mpw_ns::mpcw::operator>(const mpcw& other) const
{
	if(mpw_ns::imag(*this) != 0 || mpw_ns::imag(other) != 0)
	{
		throw std::invalid_argument("Can't compare non-real numbers!");
	}
	return MPC_INEX_RE(mpc_cmp(this->mpc_l, other.mpc_l)) > 0;
}

bool mpw_ns::mpcw::operator>(const double other) const
{
	if(mpw_ns::imag(*this) != 0)
	{
		throw std::invalid_argument("Can't compare non-real numbers!");
	}
	return *this > mpw_ns::mpcw(other);
}

bool mpw_ns::mpcw::operator>(const int other) const
{
	if(mpw_ns::imag(*this) != 0)
	{
		throw std::invalid_argument("Can't compare non-real numbers!");
	}
	return *this > mpw_ns::mpcw(other);
}

bool mpw_ns::mpcw::operator>=(const mpcw& other) const
{
	if(mpw_ns::imag(*this) != 0 || mpw_ns::imag(other) != 0)
	{
		throw std::invalid_argument("Can't compare non-real numbers!");
	}
	return MPC_INEX_RE(mpc_cmp(this->mpc_l, other.mpc_l)) >= 0;
}

bool mpw_ns::mpcw::operator>=(const double other) const
{
	if(mpw_ns::imag(*this) != 0 || mpw_ns::imag(other) != 0)
	{
		throw std::invalid_argument("Can't compare non-real numbers!");
	}
	return *this >= mpw_ns::mpcw(other);
}

bool mpw_ns::mpcw::operator>=(const int other) const
{
	if(mpw_ns::imag(*this) != 0 || mpw_ns::imag(other) != 0)
	{
		throw std::invalid_argument("Can't compare non-real numbers!");
	}
	return *this >= mpw_ns::mpcw(other);
}

bool mpw_ns::mpcw::operator<(const mpcw& other) const
{
	if(mpw_ns::imag(*this) != 0 || mpw_ns::imag(other) != 0)
	{
		throw std::invalid_argument("Can't compare non-real numbers!");
	}
	return MPC_INEX_RE(mpc_cmp(this->mpc_l, other.mpc_l)) < 0;
}

bool mpw_ns::mpcw::operator<(const double other) const
{
	if(mpw_ns::imag(*this) != 0 || mpw_ns::imag(other) != 0)
	{
		throw std::invalid_argument("Can't compare non-real numbers!");
	}
	return *this < mpw_ns::mpcw(other);
}

bool mpw_ns::mpcw::operator<(const int other) const
{
	if(mpw_ns::imag(*this) != 0 || mpw_ns::imag(other) != 0)
	{
		throw std::invalid_argument("Can't compare non-real numbers!");
	}
	return *this < mpw_ns::mpcw(other);
}

bool mpw_ns::mpcw::operator<=(const mpcw& other) const
{
	if(mpw_ns::imag(*this) != 0 || mpw_ns::imag(other) != 0)
	{
		throw std::invalid_argument("Can't compare non-real numbers!");
	}
	return MPC_INEX_RE(mpc_cmp(this->mpc_l, other.mpc_l)) <= 0;
}

bool mpw_ns::mpcw::operator<=(const double other) const
{
	if(mpw_ns::imag(*this) != 0 || mpw_ns::imag(other) != 0)
	{
		throw std::invalid_argument("Can't compare non-real numbers!");
	}
	return *this <= mpw_ns::mpcw(other);
}

bool mpw_ns::mpcw::operator<=(const int other) const
{
	if(mpw_ns::imag(*this) != 0 || mpw_ns::imag(other) != 0)
	{
		throw std::invalid_argument("Can't compare non-real numbers!");
	}
	return *this <= mpw_ns::mpcw(other);
}

bool mpw_ns::mpcw::operator!=(const mpcw& other) const
{
	return mpc_cmp(this->mpc_l, other.mpc_l) != 0;
}

bool mpw_ns::mpcw::operator!=(const double other) const
{
	return mpc_cmp(this->mpc_l, mpcw(other).mpc_l) != 0;
}

bool mpw_ns::mpcw::operator!=(const int other) const
{
	return mpc_cmp_si(this->mpc_l, other) != 0;
}

mpw_ns::mpcw mpw_ns::mpcw::abs() const
{
	mpfr_t aux;
	mpfr_init(aux);
	mpc_abs(aux, this->mpc_l, mpw_ns::mpw_defs::mpfr_rnd_type);
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_set_fr(res.mpc_l, aux, mpw_ns::mpw_defs::mpfr_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::real() const
{
	mpfr_t aux;
	mpfr_init(aux);
	mpc_real(aux, this->mpc_l, mpw_ns::mpw_defs::mpfr_rnd_type);
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_set_fr(res.mpc_l, aux, mpw_ns::mpw_defs::mpfr_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::imag() const
{
	mpfr_t aux;
	mpfr_init(aux);
	mpc_imag(aux, this->mpc_l, mpw_ns::mpw_defs::mpfr_rnd_type);
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_set_fr(res.mpc_l, aux, mpw_ns::mpw_defs::mpfr_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::conj() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_conj(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::angle() const
{
	mpfr_t aux;
	mpfr_init2(aux, mpw_ns::mpw_defs::mpw_b_prec);
	mpc_arg(aux, this->mpc_l, mpw_ns::mpw_defs::mpfr_rnd_type);
	mpw_ns::mpcw res = mpw_ns::mpcw(aux);
	return res;
}

size_t mpw_ns::mpcw::length()
{
	return 1;
}

mpw_ns::mpcw mpw_ns::mpcw::H() const
{
	return this->conj();
}

mpw_ns::mpcw mpw_ns::mpcw::sqrt() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_sqrt(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::exp() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_exp(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::sin() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_sin(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::cos() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_cos(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::tan() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_tan(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::log() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_log(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::log10() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_log10(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::sinh() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_sinh(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::cosh() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_cosh(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::tanh() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_tanh(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::asin() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_asin(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::acos() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_acos(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::atan() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_atan(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::asinh() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_asinh(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::acosh() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_acosh(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::atanh() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_atanh(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::operator+(const int in1, const mpw_ns::mpcw& in2)
{
	return mpw_ns::mpcw(in1) + in2;
}

mpw_ns::mpcw mpw_ns::operator+(const double in1, const mpw_ns::mpcw& in2)
{
	return mpw_ns::mpcw(in1) + in2;
}

mpw_ns::mpcw mpw_ns::operator-(const int in1, const mpw_ns::mpcw& in2)
{
	return mpw_ns::mpcw(in1) - in2;
}

mpw_ns::mpcw mpw_ns::operator-(const double in1, const mpw_ns::mpcw& in2)
{
	return mpw_ns::mpcw(in1) - in2;
}

mpw_ns::mpcw mpw_ns::operator*(const int in1, const mpw_ns::mpcw& in2)
{
	return mpw_ns::mpcw(in1) * in2;
}

mpw_ns::mpcw mpw_ns::operator*(const double in1, const mpw_ns::mpcw& in2)
{
	return mpw_ns::mpcw(in1) * in2;
}

mpw_ns::mpcw mpw_ns::operator/(const int in1, const mpw_ns::mpcw& in2)
{
	return mpw_ns::mpcw(in1) / in2;
}

mpw_ns::mpcw mpw_ns::operator/(const double in1, const mpw_ns::mpcw &in2)
{
	return mpw_ns::mpcw(in1) / in2;
}

mpw_ns::mpcw mpw_ns::rand()
{
	mpw_ns::mpcw res = mpw_ns::mpcw(double(std::rand()) / double(std::numeric_limits<int>::max()));
	return res;
}

mpw_ns::mpcw mpw_ns::crand()
{
	mpw_ns::mpcw res = mpw_ns::mpcw(double(std::rand()) / double(std::numeric_limits<int>::max()), double(std::rand()) / double(std::numeric_limits<int>::max()));
	return res;
}


