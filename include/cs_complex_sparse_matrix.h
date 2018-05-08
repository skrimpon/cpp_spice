/*
 * Panagiotis Skrimponis
 */

#pragma once
#include "cs.h"
#include "cs_headerdef.h"
#define Ci(k) _matrix->i[k]
#define Cp(k) _matrix->p[k]
#define Cx(k) _matrix->x[k]

template <class T>
class ComplexSparseMatrix
{
public:
	ComplexSparseMatrix(std::string, int, int, int);
	ComplexSparseMatrix(std::string, int, int, cs_ci *);
	~ComplexSparseMatrix();
	void LU();
	void BiCG(std::complex<double> *, std::complex<double> *, double);
	void solve(T *,T *);
	void compressMatrix();
	inline void operator()(int idx, int i, int j, T x);
	inline T operator()(int idx);
	inline int size() const;
	inline int nodeID() const;
	inline int branchID() const;
	inline int nonZero() const;
	inline cs_ci * matrix() const;
private:
	std::string _label;
	int _size, _nodeID, _branchID, _nonZero;
	T *_M;
	cs_ci *_matrix, *_A;
	cs_cis * _S;
	cs_cin * _N;
};

template <class T> 
inline T ComplexSparseMatrix<T>::operator()(int idx)
{
	return _matrix->x[idx];
}

template <class T> 
inline void ComplexSparseMatrix<T>::operator()(int idx, int i, int j, T x)
{
	_A->i[idx] = i; _A->p[idx] = j; _A->x[idx] = x; 
}

template <class T> 
inline int ComplexSparseMatrix<T>::size() const
{
	return _size;
}

template <class T> 
inline int ComplexSparseMatrix<T>::nodeID() const
{
	return _nodeID;
}

template <class T> 
inline int ComplexSparseMatrix<T>::branchID() const
{
	return _branchID;
}

template <class T> 
inline int ComplexSparseMatrix<T>::nonZero() const
{
	return _nonZero;
}

template <class T> 
inline cs_ci * ComplexSparseMatrix<T>::matrix() const
{
	return _matrix;
}

template <class T>
ComplexSparseMatrix<T>::ComplexSparseMatrix(std::string label, int nodeID, int branchID, int nonZero) :
_label(label), _nodeID(nodeID), _branchID(branchID), _nonZero(nonZero)
{
	_size = nodeID + branchID - 1;
	_A = cs_ci_spalloc(_size, _size, _nonZero, 1, 1);
	_A->nz = _nonZero;
	_M = new T[_size];
}

template <class T>
ComplexSparseMatrix<T>::ComplexSparseMatrix(std::string label, int nodeID, int branchID, cs_ci * matrix) :
_label(label), _nodeID(nodeID), _branchID(branchID), _matrix(matrix)
{
	_size = nodeID + branchID - 1;
	_M = new T[_size];
	for(auto k = 0; k < _size; ++k)
	{
		_M[k] = 0.0;
		for(auto j = Cp(k); (j < Cp(k+1)); ++j)
			if(Ci(j) == k)
				_M[k] = Cx(j);
	}
}

template <class T>
ComplexSparseMatrix<T>::~ComplexSparseMatrix()
{   
	delete [] _M;
}

template <class T>
void ComplexSparseMatrix<T>::compressMatrix()
{
	_matrix = cs_ci_compress(_A);
	cs_ci_spfree(_A);
	cs_ci_dupl(_matrix);
	// Get the Diag of C Matrix
	for(auto k = 0; k < _size; ++k)
	{
		_M[k] = 0.0;
		for(auto j = Cp(k); (j < Cp(k+1)); ++j)
			if(Ci(j) == k)
				_M[k] = Cx(j);
	}
}

template <class T>
void ComplexSparseMatrix<T>::LU()
{	
	_S = cs_ci_sqr(2, _matrix, 0);
	_N = cs_ci_lu(_matrix, _S, 1);
}

template <class T>
void ComplexSparseMatrix<T>::solve(T *x, T *b)
{
	T *y = new T[_size];
	for(auto i=0; i < _size; ++i) y[i] = 0;
	cs_ci_ipvec(_N->pinv, b, y, _size);
	cs_ci_lsolve(_N->L, y);
	cs_ci_usolve(_N->U, y);
	cs_ci_ipvec(_S->q, y, x, _size);
	delete y;
}

template <class T>
void ComplexSparseMatrix<T>::BiCG(std::complex<double> *x, std::complex<double> *b, double ITOL)
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

	for(j = 0; j < _size; j++)
		r[j] = b[j];
	for(j = 0; j < _size; j++)
	{		
		for(i = Cp(j) ; i < Cp(j+1); ++i)
			r[Ci(i)] -= Cx(i) * x[j];
	}
	for(i = 0; i < _size; ++i) r_bar[i] = r[i];				

	r_norm = 0.0;
	for(i = 0; i < _size; ++i) r_norm += abs(r[i] * r[i]);
	r_norm = sqrt(r_norm);

	b_norm = 0.0;
	for(i = 0; i < _size; ++i) b_norm += abs(b[i] * b[i]);
	b_norm = sqrt(b_norm);

	if(b_norm == 0.0) b_norm = 1;

	while((r_norm/b_norm > ITOL) && (iter < _size))
	{
		iter++;

		for(i = 0; i < _size; ++i) z[i] = (_M[i] != std::complex<double>(0.0, 0.0)) ? r[i] / _M[i] : r[i];
		for(i = 0; i < _size; ++i) z_bar[i] = (_M[i] != std::complex<double>(0.0, 0.0)) ? r_bar[i] / std::conj(_M[i]): r_bar[i];

		rho = std::complex<double>(0.0, 0.0);
		for(i = 0; i < _size; ++i) rho += conj(r_bar[i]) * z[i];
		if(rho == std::complex<double>(0.0, 0.0))
		{
			std::cout << "\033[1;31mERROR: BiCG Failed\033[0m \n";
			return;
		}
		if(iter == 1)
		{
			for(i = 0; i < _size; ++i) p[i] = z[i];
			for(i = 0; i < _size; ++i) p_bar[i] = z_bar[i];
		}
		else
		{
			beta = rho / rho1;
			for(i = 0; i < _size; ++i) p[i] = z[i] + beta * p[i];
			for(i = 0; i < _size; ++i) p_bar[i] = z_bar[i] + std::conj(beta) * p_bar[i];
		}

		rho1 = rho;

		for(j = 0; j < _size; ++j)	q[j] = 0.0;	
		for(j = 0; j < _size; ++j)
		{			
			for(i = Cp(j); i < Cp(j+1); i++)
				q[Ci(i)] += Cx(i) * p[j];
		}

		for(j = 0; j < _size; ++j)
		{			
			q_bar[j] = 0.0;	
			for(i = Cp(j); i < Cp(j+1); ++i)
				q_bar[j] += std::conj(Cx(i)) * p_bar[Ci(i)];
		}
		
		alpha = 0.0;
		for(i = 0; i < _size; ++i) alpha += std::conj(p_bar[i]) * q[i];
		if(alpha == 0.0)
		{
			std::cout << "\033[1;31mERROR: BiCG Failed\033[0m \n";
			return;
		}
		alpha = rho / alpha;

		for(i = 0; i < _size; ++i) x[i] += alpha * p[i];
		for(i = 0; i < _size; ++i) r[i] -= alpha * q[i];
		for(i = 0; i < _size; ++i) r_bar[i] -= std::conj(alpha) * q_bar[i];

		r_norm = 0.0;
		for(i = 0; i < _size; ++i) r_norm += abs(r[i] * r[i]);
		r_norm = sqrt(r_norm);
	}
	delete p, q, r, z, p_bar, q_bar, r_bar, z_bar;
}