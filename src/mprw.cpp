//============================================================================
// Name        : MPRW.cpp
// Author      : Pavel Vilner
// Version     :
// Copyright   : Your copyright notice
// Description : Implementation of an MPFR wrapper
//============================================================================

#include "mpw.h"

mpw_ns::mprw::mprw()
{
	mpfr_init2(this->mpr_l, mpw_ns::mpw_defs::mpw_b_prec);
}

mpw_ns::mprw::mprw(const int in)
{
	mpfr_init2(this->mpr_l, mpw_ns::mpw_defs::mpw_b_prec);
	mpfr_set_si(this->mpr_l, in, mpw_ns::mpw_defs::mprw_rnd_type);
}

mpw_ns::mprw::mprw(const double in)
{
	mpfr_init2(this->mpr_l, mpw_ns::mpw_defs::mpw_b_prec);
	mpfr_set_d(this->mpr_l, in, mpw_ns::mpw_defs::mprw_rnd_type);
}

mpw_ns::mprw::mprw(const mprw &in)
{
	mpfr_init2(this->mpr_l, mpw_ns::mpw_defs::mpw_b_prec);
	mpfr_set(this->mpr_l, in.mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
}

mpw_ns::mprw::~mprw()
{
	mpfr_clear(this->mpr_l);
}

std::string mpw_ns::mprw::str(const size_t length) const
{
	mpfr_exp_t exp_l;
	char* qqq = mpfr_get_str(NULL, &exp_l, 10, length-3, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	std::string exp_str = std::to_string(exp_l - 1);
	mpfr_free_str(qqq);
	qqq = mpfr_get_str(NULL, &exp_l, 10, length-3-exp_str.length(), this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	std::string str_ret = std::string(qqq);
	if (*this >= 0)
	{
		str_ret.insert(1, 1, '.');
		str_ret.insert(0, 1, ' ');
	}
	else
	{
		str_ret.insert(2, 1, '.');
	}
	str_ret.append(1, 'e') += std::to_string(exp_l-1);
	mpfr_free_str(qqq);
	return str_ret;
}

bool mpw_ns::mprw::isnan()
{
	return (mpfr_nan_p(this->mpr_l) != 0);
}

mpw_ns::mprw mpw_ns::mprw::operator+(const mprw &other) const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_add(res.mpr_l, this->mpr_l, other.mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::operator-(const mprw &other) const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_sub(res.mpr_l, this->mpr_l, other.mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::operator-() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_sub(res.mpr_l, mpw_ns::mprw(0).mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::operator*(const mprw &other) const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_mul(res.mpr_l, this->mpr_l, other.mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::operator/(const mprw &other) const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_div(res.mpr_l, this->mpr_l, other.mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

bool mpw_ns::mprw::operator>(const mpw_ns::mprw &other) const
{
	return (mpfr_cmp(this->mpr_l, other.mpr_l) > 0);
}

bool mpw_ns::mprw::operator<(const mpw_ns::mprw &other) const
{
	return (mpfr_cmp(this->mpr_l, other.mpr_l) < 0);
}

bool mpw_ns::mprw::operator>=(const mpw_ns::mprw &other) const
{
	return (mpfr_cmp(this->mpr_l, other.mpr_l) >= 0);
}

bool mpw_ns::mprw::operator<=(const mpw_ns::mprw &other) const
{
	return (mpfr_cmp(this->mpr_l, other.mpr_l) <= 0);
}

bool mpw_ns::mprw::operator==(const mpw_ns::mprw &other) const
{
	return (mpfr_cmp(this->mpr_l, other.mpr_l) == 0);
}

mpw_ns::mprw mpw_ns::mprw::abs() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_abs(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}
mpw_ns::mprw mpw_ns::mprw::real() const
{
	return mpw_ns::mprw(*this);
}

mpw_ns::mprw mpw_ns::mprw::imag() const
{
	return mpw_ns::mprw(0);
}

mpw_ns::mprw mpw_ns::mprw::conj() const
{
	return mpw_ns::mprw(*this);
}

mpw_ns::mprw mpw_ns::mprw::exp() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_exp(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::sin() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_sin(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::cos() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_cos(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::tan() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_tan(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::log() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_log(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::log10() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_log10(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::sinh() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_sinh(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::cosh() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_cosh(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::tanh() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_tanh(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::asin() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_asin(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::acos() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_acos(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::atan() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_atan(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::asinh() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_asinh(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::acosh() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_acosh(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mprw::atanh() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpfr_atanh(res.mpr_l, this->mpr_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::operator+(const int in1, const mpw_ns::mprw &in2)
{
	return mpw_ns::mprw(in1) + in2;
}

mpw_ns::mprw mpw_ns::operator+(const double in1, const mpw_ns::mprw &in2)
{
	return mpw_ns::mprw(in1) + in2;
}

mpw_ns::mprw operator-(const int in1, const mpw_ns::mprw &in2)
{
	return mpw_ns::mprw(in1) - in2;
}

mpw_ns::mprw mpw_ns::operator-(const double in1, const mpw_ns::mprw &in2)
{
	return mpw_ns::mprw(in1) - in2;
}

mpw_ns::mprw mpw_ns::operator*(const int in1, const mpw_ns::mprw &in2)
{
	return mpw_ns::mprw(in1) * in2;
}
mpw_ns::mprw mpw_ns::operator*(const double in1, const mpw_ns::mprw &in2)
{
	return mpw_ns::mprw(in1) * in2;
}
mpw_ns::mprw mpw_ns::operator/(const int in1, const mpw_ns::mprw &in2)
{
	return mpw_ns::mprw(in1) / in2;
}

mpw_ns::mprw mpw_ns::operator/(const double in1, const mpw_ns::mprw &in2)
{
	return mpw_ns::mprw(in1) / in2;
}

