/*
 * Panagiotis Skrimponis
 */

#pragma once
#include "cs_headerdef.h"

#define A (*this)

template <class T>
class Matrix
{
public:
	Matrix(std::string, int, int);
	~Matrix();
	inline int size() const;
	inline int nodeID() const;
	inline int branchID() const;
	inline T &operator() (int, int);
	void LU();
	void CG(T *, T *, double);
	void BiCG(double *, double *, double);
	void BiCG(std::complex<double> *, std::complex<double> *, double);
	void Cholesky();
	void solve(T *, T *);
	void solveSPD(T *, T *);
private:
	std::string _label;
	int _size, _nodeID, _branchID;
	T *_matrix;
	int *_P;
};

template <class T>
inline int Matrix<T>::size() const
{
	return _size;
}

template <class T>
inline int Matrix<T>::nodeID() const
{
	return _nodeID;
}

template <class T>
inline int Matrix<T>::branchID() const
{
	return _branchID;
}

template <class T>
inline T &Matrix<T>::operator()(int i, int j)
{
	return _matrix[i * _size + j];
}

template <class T>
Matrix<T>::Matrix(std::string label, int nodeID, int branchID) : _label(label), _nodeID(nodeID), _branchID(branchID)
{
	_size = nodeID + branchID;
	_matrix = new T[_size * _size];
	_P = new int[_size];
	for(auto i = 0; i < _size; ++i) _P[i] = i;
}

template <class T>
Matrix<T>::~Matrix()
{
	delete [] _P, _matrix;
}

template <class T>
void Matrix<T>::LU()
{
	_label = _label+" (LU)";
	int i=0, j=0, k=0, max_pos=1, temp=0;
	double max_val = abs(A(1, 1));
	
	for(i = 2; i < _size; ++i)
	{
		if(abs(A(i, 1)) > max_val)
		{
			max_val = abs(A(i, 1));
			max_pos = i;
		}
	}

	for(k = 1; k < _size; ++k)
	{
		temp = _P[k];
		_P[k] = max_pos;
		_P[max_pos] = temp;

		max_val = 0;
		max_pos = k+1;

		for(i = k + 1; i < _size; ++i)
		{	
			A(_P[i], k)   /= A(_P[k], k);
			A(_P[i], k+1) -= A(_P[i], k) * A(_P[k], k+1);

			if (abs(A(_P[i], k+1)) > max_val) {
				max_val = abs(A(_P[i], k+1));
				max_pos = i;
			}

			for(j = k + 2; j < _size; ++j)
				A(_P[i], j) -= A(_P[i], k) * A(_P[k], j);
		}
	}
}

template <class T>
void Matrix<T>::solve(T *x, T *b)
{
	unsigned int k, j;
	T *y = new T[_size];
	// Forward Solve
	for(k = 1; k < _size; ++k)
	{  
		y[k] = b[_P[k]];
		for(j = 1; j < k; ++j)
			y[k] -= A(_P[k], j) * y[j];
	}
	// Backward Solve
	for(k = _size - 1; k > 0; --k)
	{
		for(j = k + 1; j < _size; ++j)
			y[k] -= A(_P[k], j) * x[j];
		x[k] = y[k] / A(_P[k], k);
	}
	delete y;
}

template <class T>
void Matrix<T>::Cholesky()
{
	_label = _label + " (Cholesky)";
	unsigned int i, j, k;
	T sum;

	for(k = 1; k < _size; ++k)
	{
		sum = 0.0;
		for(j = 1; j < k; ++j)
			sum += A(k, j) * A(k, j);
		sum = A(k, k) - sum;

		A(k, k) = sqrt(sum);
		for(i = k + 1; i < _size; ++i)
		{
			sum = 0.0;
			for(j = 1; j < k; ++j)
				sum += A(i, j) * A(k, j);
			sum = A(k, i) - sum;

			A(i, k) = (sum / A(k, k));
		}
	}
}

template <class T>
void Matrix<T>::solveSPD(T *x, T *b)
{
	unsigned int k, j;
	T *y = new T[_size];
	// Forward Solve
	for(k = 1; k < _size; k++)
	{
		y[k] = b[k];
		for(j = 1; j < k; j++)
			y[k] -= A(k, j) * y[j];
		y[k] /= A(k, k);
	}
	// Backward Solve
	for(k = _size - 1; k > 0; k--)
	{
		for(j = k + 1; j < _size; j++)
			y[k] -= A(j, k) * x[j];
		x[k] = y[k] / A(k, k);
	}
	delete y;
}

template <class T>
void Matrix<T>::CG(T *x, T *b, double ITOL)
{
	unsigned int i, j, iter=0;
	double r_norm, b_norm;
	T rho, rho1, beta, alpha;
	T * p = new T[_size];
	T * q = new T[_size];
	T * r = new T[_size];
	T * z = new T[_size];		
	for(i=1;i<_size;++i)
	{
		r[i] = b[i];
		for(j=1;j<_size;++j)
			r[i] -= A(i,j) * x[j];
	}

	r_norm = 0.0;
	for(i = 1; i < _size; ++i) r_norm += r[i] * r[i];
	r_norm = sqrt(r_norm);

	b_norm = 0.0;
	for(i = 1; i < _size; ++i) b_norm += b[i] * b[i];
	b_norm = sqrt(b_norm);
	if(b_norm == 0) b_norm = 1;
	
	while((r_norm/b_norm > ITOL) && (iter < _size))
	{
		iter++;

		for(i = 1; i < _size; ++i) z[i] = (A(i,i) != 0) ? r[i] / A(i,i) : r[i];

		rho = 0.0;
		for(i = 1; i < _size; ++i) rho += r[i] * z[i];

		if(iter == 1)
			for(i = 1; i < _size; ++i) p[i] = z[i];
		else
		{
			beta = rho / rho1;
			for(i = 1; i < _size; ++i) p[i] = z[i] + beta * p[i];
		}

		rho1 = rho;

		for(i = 1; i < _size; ++i)
		{
			q[i] = 0.0;
			for(j = 1; j < _size; ++j)
				q[i] += A(i,j) * p[j];
		}

		alpha = 0.0;
		for(i = 1; i < _size; ++i) alpha += p[i] * q[i];
		if(alpha == 0)
		{
			std::cout << "\033[1;31mERROR: CG Failed\033[0m \n";
			return;
		}
		alpha = rho / alpha;

		for(i = 1; i < _size; ++i) x[i] += alpha * p[i];
		for(i = 1; i < _size; ++i) r[i] -= alpha * q[i];
		r_norm = 0.0;
		for(i = 1; i < _size; ++i) r_norm += r[i] * r[i];
		r_norm = sqrt(r_norm);
	}
	delete p, q, r, z;
}


template <class T>
void Matrix<T>::BiCG(double *x, double *b, double ITOL)
{
	unsigned int i, j, iter=0;
	double rho, rho1, beta, alpha, r_norm, b_norm, omega;
	double * p = new double[_size];
	double * q = new double[_size];
	double * r = new double[_size];
	double * z = new double[_size];
	double * p_bar = new double[_size];
	double * q_bar = new double[_size];
	double * r_bar = new double[_size];
	double * z_bar = new double[_size];
	
	for(i = 1; i < _size; ++i)
	{
		r[i] = b[i];
		for(j = 1; j < _size; ++j)
			r[i] -= A(i,j) * x[j];
	}
	
	for(i = 1; i < _size; ++i) r_bar[i] = r[i];
		

	r_norm = 0.0;
	for(i = 1; i < _size; ++i) r_norm += r[i] * r[i];
	r_norm = sqrt(r_norm);

	b_norm = 0.0;
	for(i = 1; i < _size; ++i) b_norm += b[i] * b[i];
	b_norm = sqrt(b_norm);
	if(b_norm == 0) b_norm = 1;

	while((r_norm/b_norm > ITOL) && (iter < _size))
	{
		iter++;

		for(i = 1;i < _size; ++i) z[i] = (A(i,i) != 0) ? (double) r[i] / A(i,i) : r[i];
		for(i = 1;i < _size; ++i) z_bar[i] = (A(i,i) != 0) ? (double) r_bar[i] / A(i,i) : r_bar[i];

		rho = 0.0;
		for(i = 1; i < _size; ++i) rho += r_bar[i] * z[i];
		if(rho == 0)
		{
			std::cout << "\033[1;31mERROR: BiCG Failed\033[0m \n";
			return;
		}
		if(iter == 1)
		{
			for(i = 1;i < _size; ++i) p[i] = z[i];
			for(i = 1;i < _size; ++i) p_bar[i] = z_bar[i];
		}
		else
		{
			beta = rho / rho1;
			for(i = 1; i < _size; ++i) p[i] = z[i] + beta * p[i];
			for(i = 1; i < _size; ++i) p_bar[i] = z_bar[i] + beta * p_bar[i];
		}

		rho1 = rho;

		for(i = 1; i < _size; ++i)
		{
			q[i] = 0.0;
			for(j = 1; j < _size; ++j)
				q[i] += A(i,j) * p[j];
		}

		for(i = 1; i < _size; ++i)
		{
			q_bar[i] = 0.0;
			for(j = 1;j < _size; ++j)
				q_bar[i] += A(j,i) * p_bar[j];
		}

		alpha = 0.0;
		for(i = 1; i < _size; ++i) alpha += p_bar[i] * q[i];
		if(alpha == 0)
		{
			std::cout << "\033[1;31mERROR: BiCG Failed\033[0m \n";
			return;
		}
		alpha = (double) rho / alpha;

		for(i = 1; i < _size; ++i) x[i] += alpha * p[i];
		for(i = 1; i < _size; ++i) r[i] -= alpha * q[i];
		for(i = 1; i < _size; ++i) r_bar[i] -= alpha * q_bar[i];

		r_norm = 0.0;
		for(i = 1; i < _size;++i) r_norm += r[i] * r[i];
		r_norm = sqrt(r_norm);
	}


	delete p, q, r, z, p_bar, q_bar, r_bar, z_bar;
}

template <class T>
void Matrix<T>::BiCG(std::complex<double> *x, std::complex<double> *b, double ITOL)
{
	unsigned int i, j, iter=0;
	double r_norm, b_norm;
	std::complex<double> rho, rho1, beta, alpha;
	std::complex<double> * r = new std::complex<double>[_size];
	std::complex<double> * r_bar = new std::complex<double>[_size];
	std::complex<double> * z = new std::complex<double>[_size];
	std::complex<double> * z_bar = new std::complex<double>[_size];
	std::complex<double> * p = new std::complex<double>[_size];
	std::complex<double> * p_bar = new std::complex<double>[_size];
	std::complex<double> * q = new std::complex<double>[_size];
	std::complex<double> * q_bar = new std::complex<double>[_size];
	
	for(i = 1; i < _size; ++i)
	{
		r[i] = b[i];
		for(j = 1; j < _size; ++j)
			r[i] -= A(i,j) * x[j];
	}
	
	for(i = 1; i < _size; ++i) r_bar[i] = r[i];
		

	r_norm = 0.0;
	for(i = 1; i < _size; ++i) r_norm += abs(r[i] * r[i]);
	r_norm = sqrt(r_norm);

	b_norm = 0.0;
	for(i = 1; i < _size; ++i) b_norm += abs(b[i] * b[i]);
	b_norm = sqrt(b_norm);
	if(b_norm == 0) b_norm = 1;

	while((r_norm/b_norm > ITOL) && (iter < _size))
	{
		iter++;

		for(i = 1;i < _size; ++i) z[i] = (A(i,i) != std::complex<double>(0.0, 0.0)) ? r[i] / A(i,i) : r[i];
		for(i = 1;i < _size; ++i) z_bar[i] = (A(i,i) != std::complex<double>(0.0, 0.0)) ? r_bar[i] / conj(A(i,i)) : r_bar[i];

		rho = std::complex<double>(0.0, 0.0);
		for(i = 1; i < _size; ++i) rho += conj(r_bar[i]) * z[i];
		if(rho == std::complex<double>(0.0, 0.0))
		{
			std::cout << "\033[1;31mERROR: BiCG Failed\033[0m \n";
			return;
		}
		if(iter == 1)
		{
			for(i = 1;i < _size; ++i) p[i] = z[i];
			for(i = 1;i < _size; ++i) p_bar[i] = z_bar[i];
		}
		else
		{
			beta = rho / rho1;
			for(i = 1; i < _size; ++i) p[i] = z[i] + beta * p[i];
			for(i = 1; i < _size; ++i) p_bar[i] = z_bar[i] + conj(beta) * p_bar[i];
		}

		rho1 = rho;

		for(i = 1; i < _size; ++i)
		{
			q[i] = std::complex<double>(0.0, 0.0);
			for(j = 1; j < _size; ++j)
				q[i] += A(i,j) * p[j];
		}

		for(i = 1; i < _size; ++i)
		{
			q_bar[i] = std::complex<double>(0.0, 0.0);
			for(j = 1;j < _size; ++j)
				q_bar[i] += conj(A(j,i)) * p_bar[j];
		}

		alpha = std::complex<double>(0.0, 0.0);
		for(i = 1; i < _size; ++i) alpha += conj(p_bar[i]) * q[i];
		if(alpha == std::complex<double>(0.0, 0.0))
		{
			std::cout << "\033[1;31mERROR: BiCG Failed\033[0m \n";
			return;
		}
		alpha = rho / alpha;

		for(i = 1; i < _size; ++i) x[i] += alpha * p[i];
		for(i = 1; i < _size; ++i) r[i] -= alpha * q[i];
		for(i = 1; i < _size; ++i) r_bar[i] -= alpha * q_bar[i];

		r_norm = 0.0;
		for(i = 1; i < _size; ++i) r_norm += abs(r[i] * r[i]);
		r_norm = sqrt(r_norm);
	}
	delete p, q, r, z, p_bar, q_bar, r_bar, z_bar;
}