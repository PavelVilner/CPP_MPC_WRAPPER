
#include "mpw.h"

template<typename _T> size_t mpw_ns::mpbm<_T>::_index(const size_t r, const size_t c) const
{
	return this->n_rows*c + r;
}

template<typename _T> mpw_ns::mpbm<_T>::mpbm()
{
	this->n_rows = 0;
	this->n_cols = 0;
	this->data = NULL;
}

template<typename _T> mpw_ns::mpbm<_T>::mpbm(const size_t rows, const size_t cols)
{
	this->n_rows = rows;
	this->n_cols = cols;
	this->data = new _T[rows*cols];
}

template<typename _T> mpw_ns::mpbm<_T>::mpbm(const size_t rows, const size_t cols, const int **ext_data);
template<typename _T> mpw_ns::mpbm<_T>::mpbm(const size_t rows, const size_t cols, const double **ext_data);
template<typename _T> mpw_ns::mpbm<_T>::mpbm(const size_t rows, const size_t cols, const _T **ext_data);
template<typename _T> mpw_ns::mpbm<_T>::mpbm(const mpw_ns::mpbm &other);

template<typename _T> mpw_ns::mpbm<_T>::~mpbm()
{
	delete[] this->data;
}

template<typename _T> std::string& mpw_ns::mpbm<_T>::str(const size_t length) const
{
	std::string res = std::string("");
	for (int r = 0; r < this->n_rows; r++)
	{
		res.append("-", length*(this->n_cols+2)) += "\n";
		for (int c = 0; c < this->n_cols; c++)
		{
			res += "| ";
			res += this->data[this->_index(r,c)].str();
			res += " ";
		}
		res += "|";
		res += "\n";
	}
	res.append("-", length*(this->n_cols+2)) += "\n";
	return res;
}

template <typename _T> const _T& mpw_ns::mpbm<_T>::operator()(const size_t r, const size_t c) const
{
	return this->data[this->_index(r,c)];
}

template <typename _T> _T& mpw_ns::mpbm<_T>::operator()(const size_t r, const size_t c)
{
	return this->data[this->_index(r,c)];
}

template<typename _T> size_t mpw_ns::mpbm<_T>::N_rows() const
{
	return this->n_rows;
}

template<typename _T> size_t mpw_ns::mpbm<_T>::N_cols() const
{
	return this->n_cols;
}

template<typename _T> bool mpw_ns::mpbm<_T>::is_Symmetric() const
{
	if (this->is_Square())
	{
		for (int r = 0; r < this->n_rows; r++)
		{
			for (int c = r; c < this->n_cols; c++)
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

template<typename _T> bool mpw_ns::mpbm<_T>::is_Hermitian() const
{
	if (this->is_Square())
	{
		for (int r = 0; r < this->n_rows; r++)
		{
			for (int c = r; c < this->n_cols; c++)
			{
				if (this->data[this->_index(r,c)] != conj(this->data[this->_index(c,r)]))
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

template<typename _T> bool mpw_ns::mpbm<_T>::is_Square() const
{
	return this->n_rows == this->n_cols;
}

template<typename _T> bool mpw_ns::mpbm<_T>::is_row_vector() const
{
	return this->n_rows == 1;
}

template<typename _T>bool mpw_ns::mpbm<_T>::is_col_vector() const
{
	return this->n_cols == 1;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>::T() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = this->data[this->_index(c,r)];
			res.data[res._index(c,r)] = this->data[this->_index(r,c)];
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>::H() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = conj(this->data[this->_index(c,r)]);
			res.data[res._index(c,r)] = conj(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>::conj() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = conj(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>::sqrt() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = sqrt(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>:: exp() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = exp(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>:: sin() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = sin(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>:: cos() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = cos(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>:: tan() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = tan(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>:: log() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = log(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>:: log10() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = log10(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>:: sinh() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = sinh(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>:: cosh() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = cosh(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>:: tanh() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = cos(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>:: asin() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = asin(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>:: acos() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = acos(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>:: atan() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = atan(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>:: asinh() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = asinh(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>:: acosh() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = acosh(this->data[this->_index(r,c)]);
		}
	}
	return res;
}

template <typename _T> mpw_ns::mpbm<_T> mpw_ns::mpbm<_T>:: atanh() const
{
	mpw_ns::mpbm<_T> res = mpw_ns::mpbm<_T>(this->n_cols, this->n_rows);
	for (int r = 0; r < this->n_rows; r++)
	{
		for (int c = r; c < this->n_cols; c++)
		{
			res.data[res._index(r,c)] = atanh(this->data[this->_index(r,c)]);
		}
	}
	return res;
}
