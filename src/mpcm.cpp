
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
		res.data[res._index(0,c)] = this->data[this->_index(row_index, c)];
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
		this->data[this->_index(row_index, c)] = new_row.data[new_row._index(0, c)];
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


