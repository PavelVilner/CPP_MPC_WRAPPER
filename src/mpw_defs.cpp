
#include "mpw.h"

size_t mpw_ns::mpw_defs::mpw_b_prec_overshoot = 5;
size_t mpw_ns::mpw_defs::mpw_d_prec = 16;
mpfr_prec_t mpw_ns::mpw_defs::mpw_b_prec = mpfr_prec_t(ceil((mpw_ns::mpw_defs::mpw_d_prec+1) * 3.322) + mpw_ns::mpw_defs::mpw_b_prec_overshoot);
mpfr_rnd_t mpw_ns::mpw_defs::mpfr_rnd_type = MPFR_RNDN;
mpc_rnd_t mpw_ns::mpw_defs::mpcw_rnd_type = MPC_RNDNN;
size_t mpw_ns::mpw_defs::str_length = 15;
size_t mpw_ns::mpw_defs::mpw_d_prec_buffer = 10;
size_t mpw_ns::mpw_defs::mpw_max_iters = 100;

void mpw_ns::mpw_defs::set_d_prec(const size_t new_prec)
{
	mpw_ns::mpw_defs::mpw_b_prec = mpfr_prec_t(ceil((new_prec+1) * 3.322) + 5);
	mpw_ns::mpw_defs::mpw_d_prec = new_prec;
}

void mpw_ns::mpw_defs::set_rnd_type(const mpfr_rnd_t new_rnd)
{
	mpw_ns::mpw_defs::mpfr_rnd_type = new_rnd;
	if (new_rnd == MPFR_RNDN)
	{
		mpw_ns::mpw_defs::mpcw_rnd_type = MPC_RNDNN;
	}
	else if (new_rnd == MPFR_RNDZ)
	{
		mpw_ns::mpw_defs::mpcw_rnd_type = MPC_RNDZZ;
	}
	else if (new_rnd == MPFR_RNDU)
	{
		mpw_ns::mpw_defs::mpcw_rnd_type = MPC_RNDUU;
	}
	else if (new_rnd == MPFR_RNDD)
	{
		mpw_ns::mpw_defs::mpcw_rnd_type = MPC_RNDDD;
	}
}

void mpw_ns::mpw_defs::set_mpcw_str_length(const size_t length)
{
	mpw_ns::mpw_defs::str_length = length;
}
mpfr_prec_t mpw_ns::mpw_defs::b_prec()
{
	return mpw_ns::mpw_defs::mpw_b_prec;
}

size_t mpw_ns::mpw_defs::d_prec()
{
	return mpw_ns::mpw_defs::mpw_d_prec;
}

mpfr_rnd_t mpw_ns::mpw_defs::rnd_type()
{
	return mpw_ns::mpw_defs::mpfr_rnd_type;
}

size_t mpw_ns::mpw_defs::mpcw_str_length()
{
	return mpw_ns::mpw_defs::str_length;
}

std::string mpw_ns::mpw_defs::mpfr_to_str(const mpfr_t val, const size_t length)
{
	char* aux_char = NULL;
	mpfr_exp_t aux_exp;
	aux_char = mpfr_get_str(aux_char, &aux_exp, 10, length, val, mpw_ns::mpw_defs::mpfr_rnd_type);
	std::string exp_string = std::string("e") + std::to_string(int(aux_exp)-1);
	aux_char = mpfr_get_str(aux_char, &aux_exp, 10, length - 1 - exp_string.length(), val, mpw_ns::mpw_defs::mpfr_rnd_type);
	std::string val_string = std::string(aux_char);
	val_string.insert(1, 1, '.');
	mpfr_free_str(aux_char);
	return val_string + exp_string;
}
