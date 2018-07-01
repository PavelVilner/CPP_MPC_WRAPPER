
#include "mpw.h"
#include <stdexcept>

size_t mpw_ns::mpcm::_index(const size_t r, const size_t c) const
{
	return this->n_rows*c + r;
}

size_t mpw_ns::mpcm::_max_index() const
{
	return this->n_rows * this->n_cols;
}

void mpw_ns::mpcm::_initialize_pivot(uint& r_p, uint& c_p, mpw_ns::mpcm& max_vals, uint* max_indexes, const size_t new_d_prec)
{
	max_vals.data[0] = mpw_ns::mpcw(0, new_d_prec); // the zero element is unused, so it's utilized to hold the total maximum
	max_indexes[0] = 0;
	r_p = 0;
	c_p = 0;
	for (uint r = 1; r < this->n_rows; r++)
	{
		max_vals.data[r] = mpw_ns::mpcw(0, new_d_prec);
		for (uint c = 0; c < r; c++)
		{
//			std::cout << mpw_ns::str(mpw_ns::abs(this->data[this->_index(r,c)])) << std::endl;
//			std::cout << mpw_ns::str(max_vals.data[r]) << std::endl;
//			std::cout << (mpw_ns::abs(this->data[this->_index(r,c)]) > max_vals.data[r]) << std::endl;
			if (mpw_ns::abs(this->data[this->_index(r,c)]) > max_vals.data[r])
			{
				max_vals.data[r] = mpw_ns::abs(this->data[this->_index(r,c)]);
				max_indexes[r] = c;
			}
			if (mpw_ns::abs(this->data[this->_index(r,c)]) > max_vals.data[0])
			{
				max_vals.data[0] = mpw_ns::abs(this->data[this->_index(r,c)]);
				r_p = r;
				c_p = c;
			}
		}
	}
}

void mpw_ns::mpcm::_update_pivot(uint& r_p, uint& c_p, mpw_ns::mpcm& max_vals, uint* max_indexes, const size_t new_d_prec)
{
	max_vals.data[0] = mpw_ns::mpcw(0, new_d_prec); // the zero element is unused, so it's utilized to hold the total maximum
	max_indexes[0] = 0;
	uint new_r_p = r_p;
	uint new_c_p = c_p;
	max_vals.data[c_p] = mpw_ns::mpcw(0, new_d_prec);
	for (uint q = 0; q < c_p; q++)
	{
//		std::cout << mpw_ns::str(mpw_ns::abs(this->data[this->_index(c_p,q)])) << std::endl;
//		std::cout << mpw_ns::str(max_vals.data[c_p]) << std::endl;
//		bool qqq = mpw_ns::abs(this->data[this->_index(c_p,q)]) > max_vals.data[c_p];
		if (mpw_ns::abs(this->data[this->_index(c_p,q)]) > max_vals.data[c_p])
		{
			max_vals.data[c_p] = mpw_ns::abs(this->data[this->_index(c_p,q)]);
			max_indexes[c_p] = q;
		}
		if (mpw_ns::abs(this->data[this->_index(c_p,q)]) > max_vals.data[0])
		{
			max_vals.data[0] = mpw_ns::abs(this->data[this->_index(c_p,q)]);
			new_r_p = c_p;
			new_c_p = q;
		}
	}
	max_vals.data[r_p] = mpw_ns::mpcw(0, new_d_prec);
	for (uint q = 0; q < r_p; q++)
	{
		if (mpw_ns::abs(this->data[this->_index(r_p,q)]) > max_vals.data[r_p])
		{
			max_vals.data[r_p] = mpw_ns::abs(this->data[this->_index(r_p,q)]);
			max_indexes[r_p] = q;
		}
		if (mpw_ns::abs(this->data[this->_index(r_p,q)]) > max_vals.data[0])
		{
			max_vals.data[0] = mpw_ns::abs(this->data[this->_index(r_p,q)]);
			new_r_p = r_p;
			new_c_p = q;
		}
	}
	for (uint q = c_p+1; q < this->n_rows; q++)
	{
		if (mpw_ns::abs(this->data[this->_index(q,c_p)]) > max_vals.data[q])
		{
			max_vals.data[q] = mpw_ns::abs(this->data[this->_index(q, c_p)]);
			max_indexes[q] = c_p;
		}
		if (mpw_ns::abs(this->data[this->_index(q,c_p)]) > max_vals.data[0])
		{
			max_vals.data[0] = mpw_ns::abs(this->data[this->_index(q, c_p)]);
			new_r_p = q;
			new_c_p = c_p;
		}
	}
	for (uint q = r_p+1; q < this->n_rows; q++)
	{
		if (mpw_ns::abs(this->data[this->_index(q,r_p)]) > max_vals.data[q])
		{
			max_vals.data[q] = mpw_ns::abs(this->data[this->_index(q, r_p)]);
			max_indexes[q] = r_p;
		}
		if (mpw_ns::abs(this->data[this->_index(q,r_p)]) > max_vals.data[0])
		{
			max_vals.data[0] = mpw_ns::abs(this->data[this->_index(q, r_p)]);
			new_r_p = q;
			new_c_p = r_p;
		}
	}
	r_p = new_r_p;
	c_p = new_c_p;
}

mpw_ns::mpcm::mpcm()
{
	this->n_rows = 0;
	this->n_cols = 0;
	this->data = NULL;
}

mpw_ns::mpcm::mpcm(const size_t rows, const size_t cols)
{
	this->n_rows = rows;
	this->n_cols = cols;
	this->data = new mpw_ns::mpcw[rows*cols];
}

mpw_ns::mpcm::mpcm(const size_t rows, const size_t cols, const int** ext_data)
{
	this->n_rows = rows;
	this->n_cols = cols;
	this->data = new mpw_ns::mpcw[rows*cols];
	for (uint c = 0; c < cols; c++)
	{
		for (uint r = 0; r < rows; r++)
		{
			this->data[this->_index(r,c)] = mpw_ns::mpcw(ext_data[c][r]);
		}
	}
}

mpw_ns::mpcm::mpcm(const size_t rows, const size_t cols, const double** ext_data)
{
	this->n_rows = rows;
	this->n_cols = cols;
	this->data = new mpw_ns::mpcw[rows*cols];
	for (uint c = 0; c < cols; c++)
	{
		for (uint r = 0; r < rows; r++)
		{
			this->data[this->_index(r,c)] = mpw_ns::mpcw(ext_data[c][r]);
		}
	}
}

mpw_ns::mpcm::mpcm(const size_t rows, const size_t cols, const mpw_ns::mpcw** ext_data)
{
	this->n_rows = rows;
	this->n_cols = cols;
	this->data = new mpw_ns::mpcw[rows*cols];
	for (uint c = 0; c < cols; c++)
	{
		for (uint r = 0; r < rows; r++)
		{
			this->data[this->_index(r,c)] = mpw_ns::mpcw(ext_data[c][r]);
		}
	}
}

mpw_ns::mpcm::mpcm(const mpw_ns::mpcm& other)
{
	this->n_rows = other.n_rows;
	this->n_cols = other.n_cols;
	this->data = new mpw_ns::mpcw[this->n_rows*this->n_cols];
	for (uint _i = 0; _i < this->n_rows*this->n_cols; _i++)
	{
		this->data[_i] = other.data[_i];
	}
}

mpw_ns::mpcm::mpcm(const mpw_ns::mpcm& other, const size_t new_d_prec)
{
	this->n_rows = other.n_rows;
	this->n_cols = other.n_cols;
	this->data = new mpw_ns::mpcw[this->n_rows*this->n_cols];
	for (uint _i = 0; _i < this->n_rows*this->n_cols; _i++)
	{
		this->data[_i] = mpw_ns::mpcw(other.data[_i], new_d_prec);
	}
}

mpw_ns::mpcm::~mpcm()
{
	delete[] this->data;
}

std::string mpw_ns::mpcm::str(const size_t length) const
{
	std::string res = std::string("");
	for (uint r = 0; r < this->n_rows; r++)
	{
		res.append((length+2)*this->n_cols + this->n_cols + 1, '-');
		res += "\n";
		for (uint c = 0; c < this->n_cols; c++)
		{
			res += "| ";
			res += this->data[this->_index(r,c)].str(length);
			res += " ";
		}
		res += "|";
		res += "\n";
	}
	res.append((length+2)*this->n_cols + this->n_cols + 1, '-');
	res += "\n";
	return res;
}

mpw_ns::mpcw mpw_ns::mpcm::max_diff(const mpcm& other) const
{
	if (this->n_rows != other.n_rows || this->n_cols != other.n_cols)
	{
		throw std::invalid_argument("Matrix dimensions must agree!");
	}
	mpw_ns::mpcw res = mpw_ns::mpcw(0);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		if (mpw_ns::abs(this->data[_i] - other.data[_i]) > res)
		{
			res = mpw_ns::abs(this->data[_i] - other.data[_i]);
		}
	}
	return res;
}

const mpw_ns::mpcw& mpw_ns::mpcm::operator()(const size_t r, const size_t c) const
{
	return this->data[this->_index(r,c)];
}

mpw_ns::mpcw& mpw_ns::mpcm::operator()(const size_t r, const size_t c)
{
	return this->data[this->_index(r,c)];
}

size_t mpw_ns::mpcm::N_rows() const
{
	return this->n_rows;
}

size_t mpw_ns::mpcm::N_cols() const
{
	return this->n_cols;
}

bool mpw_ns::mpcm::is_Symmetric() const
{
	if (this->is_Square())
	{
		for (uint r = 0; r < this->n_rows; r++)
		{
			for (uint c = r; c < this->n_cols; c++)
			{
				if (this->data[this->_index(r,c)] != this->data[this->_index(c,r)])
				{
					return false;
				}
			}
		}
		return true;
	}
	else
	{
		return false;
	}
}

bool mpw_ns::mpcm::is_Hermitian() const
{
	if (this->is_Square())
	{
		for (uint r = 0; r < this->n_rows; r++)
		{
			for (uint c = r; c < this->n_cols; c++)
			{
				if (this->data[this->_index(r,c)] != mpw_ns::conj(this->data[this->_index(c,r)]))
				{
					return false;
				}
			}
		}
		return true;
	}
	else
	{
		return false;
	}
}

bool mpw_ns::mpcm::is_Square() const
{
	return this->n_rows == this->n_cols;
}

bool mpw_ns::mpcm::is_row_vector() const
{
	return this->n_rows == 1;
}

bool mpw_ns::mpcm::is_col_vector() const
{
	return this->n_cols == 1;
}

mpw_ns::mpcm mpw_ns::mpcm::T() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_cols, this->n_rows);
	for (uint r = 0; r < res.n_rows; r++)
	{
		for (uint c = 0; c < res.n_cols; c++)
		{
			res.data[res._index(r,c)] = this->data[this->_index(c,r)];
		}
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::H() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_cols, this->n_rows);
	for (uint r = 0; r < res.n_rows; r++)
	{
		for (uint c = 0; c < res.n_cols; c++)
		{
			res.data[res._index(r,c)] = mpw_ns::conj(this->data[this->_index(c,r)]);
		}
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::normalized() const
{
	if (!this->is_col_vector() && !this->is_row_vector())
	{
		throw std::invalid_argument("Only row or column vectors can be normalized!");
	}
	mpw_ns::mpcm res = *this;
	if (this->is_col_vector())
	{
		res = res / mpw_ns::sqrt((res.H()*res).scalar());
	}
	else
	{
		res = res / mpw_ns::sqrt((res*res.H()).scalar());
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::inv() const
{
	if (!this->is_Square())
	{
		throw std::invalid_argument("Matrix must be square!");
	}
	return mpw_ns::I(this->n_rows, this->n_cols) / *this;
}

std::tuple<mpw_ns::mpcm, mpw_ns::mpcm> mpw_ns::mpcm::Householder_double_sided() const
{
	if (!this->is_Square())
	{
		throw std::invalid_argument("Matrix must be square!");
	}
	if (this->n_rows <= 2)
	{
		throw std::invalid_argument("Matrix must have more than 2 rows / columns!");
	}
	mpw_ns::mpcw i = mpw_ns::mpcw(0,1);
	mpw_ns::mpcm res = *this;
	mpw_ns::mpcm P = mpw_ns::I(this->n_rows, this->n_cols);
	for (uint c = 0; c < res.n_cols - 2; c++)
	{
			mpw_ns::mpcm loc_col = res.get_column(c);
			mpw_ns::mpcm tar_col = mpw_ns::zeros(res.n_rows,1);
			for (uint _i = 0; _i <= c; _i++)
			{
				tar_col.data[_i] = loc_col.data[_i];
			}
			mpw_ns::mpcw aux = mpw_ns::mpcw(0);
			for (uint _i = c+1; _i < loc_col.n_rows; _i++)
			{
				aux += mpw_ns::conj(loc_col.data[_i]) * loc_col.data[_i];
			}
			tar_col.data[c+1] = mpw_ns::sqrt(aux) * mpw_ns::exp(i*mpw_ns::angle(loc_col.data[c+1]));
			mpw_ns::mpcm u = (loc_col - tar_col).normalized();
			res = res - 2*u*(u.H()*res) - 2*(res*u)*u.H() + 4*(u*((u.H()*(res*u))*u.H()));
			P = P - 2*u*(u.H()*P);
	}
	return std::make_tuple(res, P.H());
}

std::tuple<mpw_ns::mpcm, mpw_ns::mpcm> mpw_ns::mpcm::Householder_single_sided() const
{
	throw std::invalid_argument("Not realized yet!");
}

void mpw_ns::mpcm::Jacobi_complex_rotator(mpw_ns::mpcm& J, const size_t r_p, const size_t c_p)
{
	mpw_ns::mpcw i = mpw_ns::mpcw(0,1);
	mpw_ns::mpcw theta = - mpw_ns::angle(this->data[this->_index(r_p, c_p)]) / 2;
	mpw_ns::mpcw x1 = mpw_ns::mpcw(1)/mpw_ns::sqrt(mpw_ns::mpcw(2)) * mpw_ns::exp(-i*theta);
	mpw_ns::mpcw x2 = mpw_ns::mpcw(1)/mpw_ns::sqrt(mpw_ns::mpcw(2)) * mpw_ns::exp(i*theta);
	mpw_ns::mpcw x3 = mpw_ns::mpcw(1)/mpw_ns::sqrt(mpw_ns::mpcw(2)) * mpw_ns::exp(-i*theta);
	mpw_ns::mpcw x4 = -mpw_ns::mpcw(1)/mpw_ns::sqrt(mpw_ns::mpcw(2)) * mpw_ns::exp(i*theta);
	this->Jacobi_rotator(J, r_p, c_p, x1, x2, x3, x4);
}

void mpw_ns::mpcm::Jacobi_real_rotator(mpw_ns::mpcm& J, const size_t r_p, const size_t c_p)
{
	mpw_ns::mpcw aux = 2*this->data[this->_index(c_p, r_p)] / (this->data[this->_index(r_p, r_p)] - this->data[this->_index(c_p, c_p)]);
	mpw_ns::mpcw theta = mpw_ns::atan(aux) / 2;
	mpw_ns::mpcw x1 = mpw_ns::cos(theta);
	mpw_ns::mpcw x2 = mpw_ns::cos(theta);
	mpw_ns::mpcw x3 = mpw_ns::sin(theta);
	mpw_ns::mpcw x4 = -mpw_ns::sin(theta);
	this->Jacobi_rotator(J, r_p, c_p, x1, x2, x3, x4);
}

void mpw_ns::mpcm::Jacobi_rotator(mpw_ns::mpcm& J, const size_t r_p, const size_t c_p,
							 	  const mpw_ns::mpcw& x1, const mpw_ns::mpcw& x2, const mpw_ns::mpcw& x3, const mpw_ns::mpcw& x4)
{
	size_t N = this->n_rows;
	mpw_ns::mpcm M_c_row = this->get_row(c_p);
	mpw_ns::mpcm M_r_row = this->get_row(r_p);
	mpw_ns::mpcm M_c_col = this->get_column(c_p);
	mpw_ns::mpcm M_r_col = this->get_column(r_p);
	mpw_ns::mpcm J_c_row = J.get_row(c_p);
	mpw_ns::mpcm J_r_row = J.get_row(r_p);
	// c-row
	for (uint q = 0; q <= c_p; q++)
	{
		if (q == c_p)
		{
			this->data[this->_index(c_p,q)] = (x1*M_c_row.data[c_p] + x4*M_r_row.data[c_p])*mpw_ns::conj(x1) +
									  (x1*M_c_row.data[r_p] + x4*M_r_row.data[r_p])*mpw_ns::conj(x4);
		}
		else
		{
			this->data[this->_index(c_p,q)] = x1*M_c_row.data[q] + x4*M_r_row.data[q];
		}
		this->data[this->_index(q, c_p)] = mpw_ns::conj(this->data[this->_index(c_p, q)]);
	}
	// r-row
	for (uint q = 0; q <= r_p; q++)
	{
		if (q == r_p)
		{
			this->data[this->_index(r_p,q)] = (x3*M_c_row.data[c_p] + x2*M_r_row.data[c_p])*mpw_ns::conj(x3) +
										   	  (x3*M_c_row.data[r_p] + x2*M_r_row.data[r_p])*mpw_ns::conj(x2);
		}
		else if (q == c_p)
		{
			this->data[this->_index(r_p,q)] = (x3*M_c_row.data[c_p] + x2*M_r_row.data[c_p])*mpw_ns::conj(x1) +
											  (x3*M_c_row.data[r_p] + x2*M_r_row.data[r_p])*mpw_ns::conj(x4);
		}
		else
		{
			this->data[this->_index(r_p, q)] = x3*M_c_row.data[q] + x2*M_r_row.data[q];
		}
		this->data[this->_index(q, r_p)] = mpw_ns::conj(this->data[this->_index(r_p, q)]);
	}
	// c-col
	for (uint q = c_p+1; q < N; q++)
	{
		if (q != r_p)
		{
			this->data[this->_index(q,c_p)] = mpw_ns::conj(x1)*M_c_col.data[q] + mpw_ns::conj(x4)*M_r_col.data[q];
			this->data[this->_index(c_p, q)] = mpw_ns::conj(this->data[this->_index(q, c_p)]);
		}
	}
	// r-col
	for (uint q = r_p+1; q < N; q++)
	{
		this->data[this->_index(q, r_p)] = mpw_ns::conj(x3)*M_c_col.data[q] + mpw_ns::conj(x2)*M_r_col.data[q];
		this->data[this->_index(r_p, q)] = mpw_ns::conj(this->data[this->_index(q, r_p)]);
	}
	// J-operations, are applied only from one side
	for (uint q = 0; q < N; q++)
	{
		J.data[J._index(c_p,q)] = x1*J_c_row.data[q] + x4*J_r_row.data[q];
		J.data[J._index(r_p,q)] = x3*J_c_row.data[q] + x2*J_r_row.data[q];
	}

}

std::tuple<mpw_ns::mpcm, mpw_ns::mpcm> mpw_ns::mpcm::Hermitian_eig(const size_t max_iter, const mpw_ns::mpcw& tol,
																   const size_t prec_buffer) const
{
	if (!this->is_Square())
	{
		throw std::invalid_argument("The matrix must be square!");
	}
	if (!this->is_Hermitian())
	{
		throw std::invalid_argument("The matrix must be Hermitian!");
	}
	mpw_ns::mpcm M_loc = mpw_ns::mpcm(*this, mpw_ns::mpw_defs::d_prec() + prec_buffer);
	mpw_ns::mpcm J = mpw_ns::I(this->n_rows, this->n_cols, mpw_ns::mpw_defs::d_prec() + prec_buffer);
	uint current_iter = 0;
	uint r_p = 0;
	uint c_p = 0;
	mpw_ns::mpcm max_vals = mpw_ns::zeros(this->n_rows, 1);
	uint* max_indexes = new uint[this->n_rows];
	M_loc._initialize_pivot(r_p, c_p, max_vals, max_indexes, mpw_ns::mpw_defs::d_prec() + prec_buffer);
	bool break_flag = false;
	while (current_iter < max_iter && !break_flag)
	{
		current_iter++;
//		if (current_iter % 100 == 0)
//		{
//			std::cout << "Iter " << current_iter << std::endl;
//		}
		mpw_ns::mpcw max_loc = max_vals.data[0];
		for (uint _i = 0; _i < M_loc.n_rows; _i++)
		{
			if (max_loc < max_vals.data[_i])
			{
				max_loc = max_vals.data[_i];
			}
		}
		if (max_loc < tol)
		{
			break_flag = true;
		}
		M_loc.Jacobi_complex_rotator(J, r_p, c_p);
		M_loc.Jacobi_real_rotator(J, r_p, c_p);
		M_loc._update_pivot(r_p, c_p, max_vals, max_indexes, mpw_ns::mpw_defs::d_prec() + prec_buffer);
	}
	delete[] max_indexes;
	if (!break_flag)
	{
		throw std::runtime_error("Jacobi Hermitian eigenvalues algorithm failed to converge to the desired precision!");
	}
	return std::make_tuple(M_loc, J.H());
}

size_t mpw_ns::mpcm::length() const
{
	if (this->n_rows >= this->n_cols)
	{
		return this->n_rows;
	}
	else
	{
		return this->n_cols;
	}
}

mpw_ns::mpcm mpw_ns::mpcm::diag() const
{
	if (this->n_rows == 1 && this->n_cols == 1)
	{
		return *this;
	}
	else if (this->n_rows == 1 || this->n_cols == 1)
	{
		mpw_ns::mpcm res = mpw_ns::zeros(this->length(), this->length());
		for (uint _i = 0; _i < this->length(); _i++)
		{
			res.data[res._index(_i, _i)] = this->data[_i];
		}
		return res;
	}
	else if (this->is_Square())
	{
		mpw_ns::mpcm res = mpw_ns::zeros(this->n_rows, 1);
		for (uint _i = 0; _i < this->n_rows; _i++)
		{
			res.data[_i] = this->data[this->_index(_i, _i)];
		}
		return res;
	}
	else
	{
		throw std::invalid_argument("Matrix must be square!");
	}
}

mpw_ns::mpcw mpw_ns::mpcm::scalar() const
{
	return this->data[0];
}

mpw_ns::mpcm mpw_ns::mpcm::get_row(size_t row_index) const
{
	if (row_index >= this->n_rows)
	{
		throw std::invalid_argument("Row index is larger than the matrix max row index!");
	}
	mpw_ns::mpcm res = mpw_ns::mpcm(1, this->n_cols);
	for (uint c = 0; c < this->n_cols; c++)
	{
		res.data[c] = this->data[this->_index(row_index, c)];
	}
	return res;
}

void mpw_ns::mpcm::set_row(size_t row_index, const mpw_ns::mpcm& new_row)
{
	if (row_index >= this->n_rows)
	{
		throw std::invalid_argument("Row index is larger than the matrix max row index!");
	}
	if (new_row.n_cols != this->n_cols || new_row.n_rows != 1)
	{
		throw std::invalid_argument("Matrix dimensions must agree!");
	}
	for (uint c = 0; c < this->n_cols; c++)
	{
		this->data[this->_index(row_index, c)] = new_row.data[c];
	}
}

mpw_ns::mpcm mpw_ns::mpcm::get_column(size_t col_index) const
{
	if (col_index >= this->n_cols)
	{
		throw std::invalid_argument("Column index is larger than the matrix max column index!");
	}
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, 1);
	for (uint r = 0; r < this->n_cols; r++)
	{
		res.data[r] = this->data[this->_index(r, col_index)];
	}
	return res;
}

void mpw_ns::mpcm::set_column(size_t col_index, const mpw_ns::mpcm& new_col)
{
	if (col_index >= this->n_cols)
	{
		throw std::invalid_argument("Column index is larger than the matrix max column index!");
	}
	if (new_col.n_rows != this->n_rows || new_col.n_cols != 1)
	{
		throw std::invalid_argument("Matrix dimensions must agree!");
	}
	for (uint r = 0; r < this->n_rows; r++)
	{
		this->data[this->_index(r, col_index)] = new_col.data[r];
	}
}

mpw_ns::mpcw mpw_ns::mpcw::T() const
{
	return *this;
}

mpw_ns::mpcm mpw_ns::mpcm::conj() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i ++)
	{
		res.data[_i] = mpw_ns::conj(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::real() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i ++)
	{
		res.data[_i] = mpw_ns::real(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::imag() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i ++)
	{
		res.data[_i] = mpw_ns::imag(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::sqrt() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = mpw_ns::sqrt(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::exp() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = mpw_ns::exp(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::sin() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = mpw_ns::sin(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm:: cos() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = mpw_ns::cos(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm:: tan() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = mpw_ns::tan(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::log() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = mpw_ns::log(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm:: log10() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = mpw_ns::log10(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm:: sinh() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = mpw_ns::sinh(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm:: cosh() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = mpw_ns::cosh(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm:: tanh() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = mpw_ns::tanh(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm:: asin() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = mpw_ns::asin(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm:: acos() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = mpw_ns::acos(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm:: atan() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = mpw_ns::atan(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm:: asinh() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = mpw_ns::asinh(this->data[_i]);
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm:: acosh() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = mpw_ns::acosh(this->data[_i]);
	}
	return res;
}

void mpw_ns::mpcm::operator=(const mpw_ns::mpcm& other)
{
	delete[] this->data;
	this->n_rows = other.n_rows;
	this->n_cols = other.n_cols;
	this->data = new mpw_ns::mpcw[this->n_rows * this->n_cols];
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		this->data[_i] = other.data[_i];
	}
}

mpw_ns::mpcm mpw_ns::mpcm::operator+(const mpw_ns::mpcm& other) const
{
	if ((this->n_rows != other.n_rows) || (this->n_cols != other.n_cols))
	{
		throw std::invalid_argument("Matrix dimensions must agree!");
	}
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = this->data[_i] + other.data[_i];
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator+(const int other) const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = this->data[_i] + other;
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator+(const double other) const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = this->data[_i] + other;
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator+(const mpw_ns::mpcw& other) const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = this->data[_i] + other;
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator-() const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = -this->data[_i];
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator-(const mpw_ns::mpcm& other) const
{
	if ((this->n_rows != other.n_rows) || (this->n_cols != other.n_cols))
	{
		throw std::invalid_argument("Matrix dimensions must agree!");
	}
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = this->data[_i] - other.data[_i];
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator-(const int other) const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = this->data[_i] - other;
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator-(const double other) const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = this->data[_i] - other;
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator-(const mpw_ns::mpcw& other) const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = this->data[_i] - other;
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator*(const mpw_ns::mpcm& other) const
{
	if (this->n_cols != other.n_rows)
	{
		throw std::invalid_argument("Matrix dimensions must agree!");
	}
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, other.n_cols);
	if (this->n_rows == 1 && other.n_cols == 1)
	{
		res.data[0] = mpw_ns::mpcw(0);
		for (uint _i = 0; _i < this->n_cols; _i++)
		{
			res.data[0] += this->data[_i] * other.data[_i];
		}
		return res;
	}
	if (this->n_rows == 1)
	{
		for (uint c = 0; c < res.n_cols; c++)
		{
			res.data[c] = mpw_ns::mpcw(0);
			for (uint _i = 0; _i < other.n_rows; _i++)
			{
				res.data[c] += this->data[_i] * other.data[other._index(_i,c)];
			}
		}
		return res;
	}
	if (other.n_cols == 1)
	{
		for (uint r = 0; r < res.n_rows; r++)
		{
			res.data[r] = mpw_ns::mpcw(0);
			for (uint _i = 0; _i < other.n_rows; _i++)
			{
				res.data[r] += this->data[this->_index(r,_i)] * other.data[_i];
			}
		}
		return res;
	}
	for (uint r = 0; r < this->n_rows; r++)
	{
		for (uint c = 0; c < other.n_cols; c++)
		{
			res.data[res._index(r,c)] = mpw_ns::mpcw(0,0);
			for (uint _i = 0; _i < this->n_cols; _i++)
			{
				res.data[res._index(r,c)] += this->data[this->_index(r,_i)] * other.data[other._index(_i,c)];
			}
		}
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator*(const int other) const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = this->data[_i] * other;
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator*(const double other) const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = this->data[_i] * other;
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator*(const mpw_ns::mpcw& other) const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i  = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = this->data[_i] * other;
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator/(const mpw_ns::mpcm& other) const
{
	// REALIZES GRAM-SCHMIDT ALGORITHM
	if (!this->is_Square() || !other.is_Square() || this->n_rows != other.n_rows)
	{
		throw std::invalid_argument("Matrix dimensions must agree!");
	}
	mpw_ns::mpcm res = *this;
	mpw_ns::mpcm op = other;
	// lower-triangularization
	for (uint r = 0; r < op.n_rows-1; r++)
	{
		// resolve a possibility of the leading element being zero
		uint r_init = r;
		while (op.data[op._index(r_init, r)] == 0 && r_init < op.n_rows)
		{
			r_init++;
		}
		if (r_init >= op.n_rows)
		{
			throw std::invalid_argument("The denominator matrix appears to be non-invertible!");
		}
		// normalize the leading element of the current row to 1
		res.set_row(r_init, res.get_row(r_init) / op.data[op._index(r_init,r)]);
		op.set_row(r_init, op.get_row(r_init) / op.data[op._index(r_init,r)]);
		if (r_init > r)
		{
			res.set_row(r, res.get_row(r) + res.get_row(r_init));
			op.set_row(r, op.get_row(r) + op.get_row(r_init));
		}
		// go over the rows below the current
		for (uint r_loc = r+1; r_loc < other.n_rows; r_loc++)
		{
			if (op.data[op._index(r_loc,r)] == 0)
			{
				throw std::invalid_argument("The denominator matrix appears to be non-invertible!");
			}
			res.set_row(r_loc, res.get_row(r_loc) - res.get_row(r)*op.data[op._index(r_loc,r)]);
			op.set_row(r_loc, op.get_row(r_loc) - op.get_row(r)*op.data[op._index(r_loc,r)]);
		}
	}
	if (op.data[op._index(op.n_rows-1, op.n_cols-1)] == 0)
	{
		throw std::invalid_argument("The denominator matrix appears to be non-invertible!");
	}
	res.set_row(res.n_rows-1, res.get_row(res.n_rows-1) / op.data[op._index(op.n_rows-1, op.n_cols-1)]);
	op.set_row(op.n_rows-1, op.get_row(op.n_rows-1) / op.data[op._index(op.n_rows-1, op.n_cols-1)]);
	// upper-triangularization
	for (int r = op.n_rows-1; r > 0; r--)
	{
		for (int r_loc = r - 1; r_loc >= 0; r_loc--)
		{
			res.set_row(r_loc, res.get_row(r_loc) - res.get_row(r)*op.data[op._index(r_loc,r)]);
			op.set_row(r_loc, op.get_row(r_loc) - op.get_row(r)*op.data[op._index(r_loc,r)]);
		}
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator/(const int other) const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i  = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = this->data[_i] / other;
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator/(const double other) const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = this->data[_i] / other;
	}
	return res;
}

mpw_ns::mpcm mpw_ns::mpcm::operator/(const mpcw& other) const
{
	mpw_ns::mpcm res = mpw_ns::mpcm(this->n_rows, this->n_cols);
	for (uint _i = 0; _i < this->_max_index(); _i++)
	{
		res.data[_i] = this->data[_i] / other;
	}
	return res;
}

mpw_ns::mpcm mpw_ns::zeros(const size_t n_rows, const size_t n_cols)
{
	mpw_ns::mpcm res = mpw_ns::mpcm(n_rows, n_cols);
	for (uint r = 0; r < n_rows; r++)
	{
		for (uint c = 0; c < n_cols; c++)
		{
			res(r,c) = mpw_ns::mpcw(0,0);
		}
	}
	return res;
}

mpw_ns::mpcm mpw_ns::ones(const size_t n_rows, const size_t n_cols)
{
	mpw_ns::mpcm res = mpw_ns::mpcm(n_rows, n_cols);
	for (uint r = 0; r < n_rows; r++)
	{
		for (uint c = 0; c < n_cols; c++)
		{
			res(r,c) = mpw_ns::mpcw(1,0);
		}
	}
	return res;
}

mpw_ns::mpcm mpw_ns::I(const size_t n_rows, const size_t n_cols)
{
	mpw_ns::mpcm res = mpw_ns::mpcm(n_rows, n_cols);
	for (uint r = 0; r < n_rows; r++)
	{
		for (uint c = 0; c < n_cols; c++)
		{
			if (r == c)
			{
				res(r,c) = mpw_ns::mpcw(1,0);
			}
			else
			{
				res(r,c) = mpw_ns::mpcw(0,0);
			}
		}
	}
	return res;
}

mpw_ns::mpcm mpw_ns::I(const size_t n_rows, const size_t n_cols, const size_t new_d_prec)
{
	mpw_ns::mpcm res = mpw_ns::mpcm(n_rows, n_cols);
	for (uint r = 0; r < n_rows; r++)
	{
		for (uint c = 0; c < n_cols; c++)
		{
			if (r == c)
			{
				res(r,c) = mpw_ns::mpcw(1,new_d_prec);
			}
			else
			{
				res(r,c) = mpw_ns::mpcw(0,new_d_prec);
			}
		}
	}
	return res;
}

mpw_ns::mpcm mpw_ns::rand(const size_t n_rows, const size_t n_cols)
{
	mpw_ns::mpcm res = mpw_ns::mpcm(n_rows, n_cols);
	std::srand(std::time(nullptr));
	for (uint r = 0; r < n_rows; r++)
	{
		for (uint c = 0; c < n_cols; c++)
		{
			res(r,c) = mpw_ns::rand();
		}
	}
	return res;
}

mpw_ns::mpcm mpw_ns::crand(const size_t n_rows, const size_t n_cols)
{
	mpw_ns::mpcm res = mpw_ns::mpcm(n_rows, n_cols);
	std::srand(std::time(nullptr));
	for (uint r = 0; r < n_rows; r++)
	{
		for (uint c = 0; c < n_cols; c++)
		{
			res(r,c) = mpw_ns::crand();
		}
	}
	return res;
}

mpw_ns::mpcm mpw_ns::operator+(const int in1, const mpw_ns::mpcm &in2)
{
	return in2 + in1;
}

mpw_ns::mpcm mpw_ns::operator+(const double in1, const mpw_ns::mpcm &in2)
{
	return in2 + in1;
}

mpw_ns::mpcm mpw_ns::operator+(const mpw_ns::mpcw& in1, const mpw_ns::mpcm &in2)
{
	return in2 + in1;
}

mpw_ns::mpcm mpw_ns::operator-(const int in1, const mpw_ns::mpcm &in2)
{
	return -in2 + in1;
}

mpw_ns::mpcm mpw_ns::operator-(const double in1, const mpw_ns::mpcm &in2)
{
	return -in2 + in1;
}

mpw_ns::mpcm mpw_ns::operator-(const mpw_ns::mpcw& in1, const mpw_ns::mpcm &in2)
{
	return -in2 + in1;
}

mpw_ns::mpcm mpw_ns::operator*(const int in1, const mpw_ns::mpcm& in2)
{
	return in2 * in1;
}

mpw_ns::mpcm mpw_ns::operator*(const double in1, const mpw_ns::mpcm &in2)
{
	return in2 * in1;
}

mpw_ns::mpcm mpw_ns::operator*(const mpw_ns::mpcw& in1, const mpw_ns::mpcm &in2)
{
	return in2 * in1;
}

mpw_ns::mpcm mpw_ns::operator/(const int in1, const mpw_ns::mpcm& in2)
{
	return in1 * (in2.inv());
}

mpw_ns::mpcm mpw_ns::operator/(const double in1, const mpw_ns::mpcm& in2)
{
	return in1 * (in2.inv());
}

mpw_ns::mpcm mpw_ns::operator/(const mpw_ns::mpcw& in1, const mpw_ns::mpcm& in2)
{
	return in1 * in2.inv();
}

mpw_ns::mpcm mpw_ns::int_test(const size_t n_rows, const size_t n_cols)
{
	mpw_ns::mpcm res = mpw_ns::mpcm(n_rows, n_cols);
	int aux = 1;
	for (uint r = 0; r < n_rows; r++)
	{
		for (uint c = 0; c < n_cols; c++)
		{
			res(r,c) = mpw_ns::mpcw(aux);
			aux++;
		}
	}
	return res;
}


