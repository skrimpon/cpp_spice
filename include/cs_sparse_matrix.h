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
class SparseMatrix
{
public:
	SparseMatrix(std::string, int, int, int);
	SparseMatrix(std::string, int, int, cs *);
	~SparseMatrix();
	void LU();
	void CG(T *, T *, double);
	void BiCG(double *, double *, double);
	void Cholesky();
	void solve(T *,T *);
	void solveSPD(T *,T *);
	void compressMatrix();
	inline void operator()(int idx, int i, int j, T x);
	inline int size() const;
	inline int nodeID() const;
	inline int branchID() const;
	inline int nonZero() const;
	inline cs_di * matrix() const;
private:
	std::string _label;
	int _size, _nodeID, _branchID, _nonZero;
	T *_M;
	cs_di *_matrix, *_A;
	cs_dis * _S;
	cs_din * _N;
};

template <class T> 
inline void SparseMatrix<T>::operator()(int idx, int i, int j, T x)
{
	_A->i[idx] = i; _A->p[idx] = j; _A->x[idx] = x; 
}

template <class T> 
inline int SparseMatrix<T>::size() const
{
	return _size;
}

template <class T> 
inline int SparseMatrix<T>::nodeID() const
{
	return _nodeID;
}

template <class T> 
inline int SparseMatrix<T>::branchID() const
{
	return _branchID;
}

template <class T> 
inline int SparseMatrix<T>::nonZero() const
{
	return _nonZero;
}

template <class T> 
inline cs_di * SparseMatrix<T>::matrix() const
{
	return _matrix;
}

template <class T>
SparseMatrix<T>::SparseMatrix(std::string label, int nodeID, int branchID, int nonZero) :
_label(label), _nodeID(nodeID), _branchID(branchID), _nonZero(nonZero)
{
	_size = nodeID + branchID - 1;
	_A = cs_di_spalloc(_size, _size, _nonZero, 1, 1);
	_A->nz = _nonZero;
	_M = new T[_size];
}

template <class T>
SparseMatrix<T>::SparseMatrix(std::string label, int nodeID, int branchID, cs_di * matrix) :
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
SparseMatrix<T>::~SparseMatrix()
{   
	delete [] _M;
}

template <class T>
void SparseMatrix<T>::compressMatrix()
{
	_matrix = cs_di_compress(_A);
	cs_di_spfree(_A);
	cs_di_dupl(_matrix);
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
void SparseMatrix<T>::LU()
{	
	_S = cs_di_sqr(2, _matrix, 0);
	_N = cs_di_lu(_matrix, _S, 1);
}

template <class T>
void SparseMatrix<T>::Cholesky()
{
	_S = cs_di_schol(1, _matrix);
	_N = cs_di_chol(_matrix, _S);
}

template <class T>
void SparseMatrix<T>::solve(T *x, T *b)
{
	T *y = new T[_size];
	for(auto i=0; i < _size; ++i) y[i] = 0;
	cs_di_ipvec(_N->pinv, b, y, _size);
	cs_di_lsolve(_N->L, y);
	cs_di_usolve(_N->U, y);
	cs_di_ipvec(_S->q, y, x, _size);
	delete y;
}

template <class T>
void SparseMatrix<T>::solveSPD(T *x, T *b)
{
	T *y = new T[_size];
	for(int i=0; i < _size; ++i) y[i] = 0;
	cs_di_ipvec(_S->pinv, b, y, _size);
	cs_di_lsolve(_N->L, y);
	cs_di_ltsolve(_N->L, y);
	cs_di_pvec(_S->pinv, y, x, _size);
	delete y;
}

template <class T>
void SparseMatrix<T>::CG(T *x, T *b, double ITOL)
{
	unsigned int i, j, k, iter=0;
	double r_norm, b_norm;
	T rho, rho1, beta, alpha;
	T * p = new T[_size];
	T * q = new T[_size];
	T * r = new T[_size];
	T * z = new T[_size];	
	for(j = 0; j < _size; ++j) r[j] = b[j];	
	for(j = 0; j < _size; ++j)
	{		
		for(i = Cp(j) ; i < Cp(j+1); ++i)
			r[Ci(i)] -= Cx(i) * x[j];
	}

	r_norm = 0.0;
	for(i = 0; i < _size; ++i) r_norm += r[i] * r[i];
	r_norm = sqrt(r_norm);

	b_norm = 0.0;
	for(i = 0; i < _size; ++i) b_norm += b[i] * b[i];
	b_norm = sqrt(b_norm);
	if(b_norm == 0) b_norm = 1;

	while((r_norm/b_norm > ITOL) && (iter < _size))
	{
		iter++;

		for(i = 0; i < _size; ++i) z[i] = (_M[i] != 0) ? (double) r[i] / _M[i] : r[i];

		rho = 0.0;
		for (i = 0; i < _size; ++i) rho += r[i] * z[i];
		
		if(iter == 1)
			for(i = 0; i < _size; ++i) p[i] = z[i];
		else
		{
			beta = rho / rho1;
			for(i = 0; i < _size; ++i) p[i] = z[i] + beta * p[i];
		}

		rho1 = rho;

		for(j = 0; j < _size; ++j)	q[j] = 0.0;	
		for(j = 0; j < _size; ++j)
		{			
			for(i = Cp(j); i < Cp(j+1); ++i)
				q[Ci(i)] += Cx(i) * p[j];
		}

		alpha = 0.0;
		for(i = 0; i < _size; i++) alpha += p[i] * q[i];

		if(alpha == 0)
		{
			std::cout << "\033[1;31mERROR: CG Failed\033[0m \n";
			return;
		}

		alpha = rho / alpha;

		for(i = 0; i < _size; ++i)x[i] += alpha * p[i];
		for(i = 0; i < _size; ++i) r[i] -= alpha * q[i];
		r_norm = 0.0;
		for(i = 0; i < _size; ++i) r_norm += r[i] * r[i];
		r_norm = sqrt(r_norm);
	}

	delete p, q, r, z;
}

template <class T>
void SparseMatrix<T>::BiCG(double *x, double *b, double ITOL)
{
	unsigned int i, j, k, iter=0;
	double rho, rho1, beta, alpha, r_norm, b_norm;
	double * p = new double[_size];
	double * q = new double[_size];
	double * r = new double[_size];
	double * z = new double[_size];
	double * p_bar = new double[_size];
	double * q_bar = new double[_size];
	double * r_bar = new double[_size];
	double * z_bar = new double[_size];

	for(j = 0; j < _size; j++) r[j] = b[j];	
	for(j = 0; j < _size; j++)
	{		
		for(i = Cp(j) ; i < Cp(j+1); ++i)
			r[Ci(i)] -= Cx(i) * x[j];
	}


	for(i = 0; i < _size; ++i) r_bar[i] = r[i];				

	r_norm = 0.0;
	for(i = 0; i < _size; ++i) r_norm += r[i] * r[i];
	r_norm = sqrt(r_norm);

	b_norm = 0.0;
	for(i = 0; i < _size; ++i) b_norm += b[i] * b[i];
	b_norm = sqrt(b_norm);
	if(b_norm == 0) b_norm = 1;

	while((r_norm/b_norm > ITOL) && (iter < _size))
	{
		iter++;

		for(i = 0; i < _size; ++i) z[i] = (_M[i] != 0) ? (double) r[i] / _M[i] : r[i];
		for(i = 0; i < _size; ++i) z_bar[i] = (_M[i] != 0) ? (double) r_bar[i] / _M[i]: r_bar[i];

		rho = 0.0;
		for(i = 0; i < _size; ++i) rho += r_bar[i] * z[i];
		if(rho == 0)
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
			for(i = 0; i < _size; ++i) p_bar[i] = z_bar[i] + beta * p_bar[i];
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
				q_bar[j] += Cx(i) * p_bar[Ci(i)];
		}
		
		alpha = 0.0;
		for(i = 0; i < _size; ++i) alpha += p_bar[i] * q[i];
		if(alpha == 0)
		{
			std::cout << "\033[1;31mERROR: BiCG Failed\033[0m \n";
			return;
		}
		alpha = (double) rho / alpha;

		for(i = 0; i < _size; ++i) x[i] += alpha * p[i];
		for(i = 0; i < _size; ++i) r[i] -= alpha * q[i];
		for(i = 0; i < _size; ++i) r_bar[i] -= alpha * q_bar[i];

		r_norm = 0.0;
		for(i = 0; i < _size; ++i) r_norm += r[i] * r[i];
		r_norm = sqrt(r_norm);
	}
	delete p, q, r, z, p_bar, q_bar, r_bar, z_bar;
}