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
	mpc_init2(this->mpc_l, this->mpw_ns::mpw_defs::mpw_b_prec);
}

mpw_ns::mpcw::mpcw(const mprw &re_in)
{
	mpc_init2(this->mpc_l, this->mpw_ns::mpw_defs::mpw_b_prec);
	mpc_set_fr(this->mpc_l, re_in.mpr_l, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::mpcw(const int re_in)
{
	mpc_init2(this->mpc_l, this->mpw_ns::mpw_defs::mpw_b_prec);
	mpc_set_si(this->mpc_l, re_in, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::mpcw(const int re_in, const int im_in)
{
	mpc_init2(this->mpc_l, this->mpw_ns::mpw_defs::mpw_b_prec);
	mpc_set_si_si(this->mpc_l, re_in, im_in, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::mpcw(const double re_in)
{
	mpc_init2(this->mpc_l, this->mpw_ns::mpw_defs::mpw_b_prec);
	mpc_set_d(this->mpc_l, re_in, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::mpcw(const double re_in, const double im_in)
{
	mpc_init2(this->mpc_l, this->mpw_ns::mpw_defs::mpw_b_prec);
	mpc_set_d_d(this->mpc_l, re_in, im_in, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::mpcw(const mprw &re_in, const mprw &im_in)
{
	mpc_init2(this->mpc_l, this->mpw_ns::mpw_defs::mpw_b_prec);
	mpc_set_fr_fr(this->mpc_l, re_in.mpr_l, im_in.mpr_l, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::mpcw(const mpcw &in)
{
	mpc_init2(this->mpc_l, this->mpw_ns::mpw_defs::mpw_b_prec);
	mpc_set(this->mpc_l, in.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
}

mpw_ns::mpcw::~mpcw()
{
	mpc_clear(this->mpc_l);
}

std::string mpw_ns::mpcw::str(const size_t length = 15) const
{
	std::string filler = std::string("");
	if (std::floor(length/2.0) == length/2.0)
	{
		filler = std::string(" ");
	}
	mpw_ns::mprw re_p = this->real();
	mpw_ns::mprw im_p = this->imag();
	std::string sign_re = std::string(" ");
	std::string sign_im = std::string(" + j");
	std::string abs_re_s = this->real().abs().str(int((length-3)/2)).erase(0,1);
	std::string abs_im_s = this->imag().abs().str(int((length-3)/2)).erase(0,1);
	if (re_p < 0)
	{
		sign_re = std::string("-");
	}
	if (im_p < 0)
	{
		sign_im = std::string(" - j");
	}
	return sign_re + abs_re_s + sign_im + abs_im_s + filler;
}

bool mpw_ns::mpcw::isnan()
{
	return this->real().isnan() || this->imag().isnan();
}

mpw_ns::mpcw mpw_ns::mpcw::operator+(const mpcw &other) const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_add(res.mpc_l, this->mpc_l, other.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::operator-(const mpcw &other) const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_sub(res.mpc_l, this->mpc_l, other.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::operator-() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_sub(res.mpc_l, mpw_ns::mpcw(0,0).mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::operator*(const mpcw &other) const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_mul(res.mpc_l, this->mpc_l, other.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::operator/(const mpcw &other) const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_div(res.mpc_l, this->mpc_l, other.mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
}

bool mpw_ns::mpcw::operator==(const mpcw &other) const
{
	return mpc_cmp(this->mpc_l, other.mpc_l) == 0;
}

mpw_ns::mprw mpw_ns::mpcw::abs() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpc_abs(res.mpr_l, this->mpc_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mpcw::real() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpc_real(res.mpr_l, this->mpc_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mprw mpw_ns::mpcw::imag() const
{
	mpw_ns::mprw res = mpw_ns::mprw();
	mpc_imag(res.mpr_l, this->mpc_l, mpw_ns::mpw_defs::mprw_rnd_type);
	return res;
}

mpw_ns::mpcw mpw_ns::mpcw::conj() const
{
	mpw_ns::mpcw res = mpw_ns::mpcw();
	mpc_conj(res.mpc_l, this->mpc_l, mpw_ns::mpw_defs::mpcw_rnd_type);
	return res;
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

mpw_ns::mpcw mpw_ns::operator+(const int in1, const mpw_ns::mpcw &in2)
{
	return mpw_ns::mpcw(in1) + in1;
}

mpw_ns::mpcw mpw_ns::operator+(const double in1, const mpw_ns::mpcw &in2)
{
	return mpw_ns::mpcw(in1) + in1;
}

mpw_ns::mpcw mpw_ns::operator+(const mpw_ns::mprw &in1, const mpw_ns::mpcw &in2)
{
	return mpw_ns::mpcw(in1) + in1;
}

mpw_ns::mpcw mpw_ns::operator-(const int in1, const mpw_ns::mpcw &in2)
{
	return mpw_ns::mpcw(in1) - in1;
}

mpw_ns::mpcw mpw_ns::operator-(const double in1, const mpw_ns::mpcw &in2)
{
	return mpw_ns::mpcw(in1) - in1;
}

mpw_ns::mpcw mpw_ns::operator-(const mpw_ns::mprw &in1, const mpw_ns::mpcw &in2)
{
	return mpw_ns::mpcw(in1) - in1;
}

mpw_ns::mpcw mpw_ns::operator*(const int in1, const mpw_ns::mpcw &in2)
{
	return mpw_ns::mpcw(in1) * in1;
}

mpw_ns::mpcw mpw_ns::operator*(const double in1, const mpw_ns::mpcw &in2)
{
	return mpw_ns::mpcw(in1) * in1;
}

mpw_ns::mpcw mpw_ns::operator*(const mpw_ns::mprw &in1, const mpw_ns::mpcw &in2)
{
	return mpw_ns::mpcw(in1) * in1;
}

mpw_ns::mpcw mpw_ns::operator/(const int in1, const mpw_ns::mpcw &in2)
{
	return mpw_ns::mpcw(in1) / in2;
}

mpw_ns::mpcw mpw_ns::operator/(const double in1, const mpw_ns::mpcw &in2)
{
	return mpw_ns::mpcw(in1) / in2;
}

mpw_ns::mpcw mpw_ns::operator/(const mpw_ns::mprw &in1, const mpw_ns::mpcw &in2)
{
	return mpw_ns::mpcw(in1) / in2;
}


