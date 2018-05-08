/*
 * Panagiotis Skrimponis
 */

#pragma once
#include "cs_headerdef.h"
/*
 * TransientComponent Class
 */
class TransientComponent
{
public:
	TransientComponent(char, std::string, std::stringstream&, double, int, int, int); // Class Constructor.
	~TransientComponent();
	inline char type() const;
	inline int plus() const;
	inline int minus() const;
	inline int branch() const;
	inline double dc_value() const;
	double value(double) const;
private:
	char _component_type;
	int _component_plus, _component_minus, _component_branch;
	std::string _name, _type;
	double _dc_value, _i1, _i2, _td1, _tc1, _td2, _tc2, _td, _tr, _tf, _pw, _per, _ia, _fr, _df, _ph;
	std::vector < double > _pwl_time, _pwl_value;
};

inline char TransientComponent::type() const
{
	return _component_type;
}

inline int TransientComponent::plus() const
{
	return _component_plus;
}

inline int TransientComponent::minus() const
{
	return _component_minus;
}

inline int TransientComponent::branch() const
{
	return _component_branch;
}

inline double TransientComponent::dc_value() const
{
	return _dc_value;
}
/*
 * TransientAnalysis Class
 */
class TransientAnalysis
{
public:
	TransientAnalysis(double, double); // Class Constructor.
	~TransientAnalysis();
	inline void add_node(std::string);
	template <class T>
	void analyse(T&, T&, double *, double *, double, bool, bool, bool, std::vector<TransientComponent>&);
	template <class T>
	void sp_analyse(T&, T&, double *, double *, double, bool, bool, bool, std::vector<TransientComponent>&);
	void update(std::unordered_map < std::string, int >&);
private:
	int * _nodes; 
	std::vector < std::string > _node_vector;
	double _step, _stop;
};
inline void TransientAnalysis::add_node(std::string node)
{
	_node_vector.push_back(node);
}



template <class T>
void TransientAnalysis::analyse(T& C, T& G, double * x, double * b, double itol, bool iter, bool spd, bool method, std::vector<TransientComponent>& TransientComponentVector)
{
	unsigned int i, size = C.size();
	double * e = new double[size]; for (auto j = 0; j < size; ++j) e[j]=0;
	std::ofstream * oFile_v = new std::ofstream[_node_vector.size()];
	double sum;
	
	for(auto &it : TransientComponentVector)
	{
		if (it.type() == 'V')
			b[it.branch() + C.nodeID()] = it.value(0);
		else
		{
			b[it.plus ()] += it.dc_value() - it.value(0);
			b[it.minus()] -= it.dc_value() - it.value(0);
		}
	}
	
	for(i = 0; i < _node_vector.size(); ++i)
	{
		oFile_v[i].open("Node_" + _node_vector[i] + "_Transient_" + std::to_string(_step) + "_" + std::to_string(_stop) + ".txt");
	}

	if(method)
	{	
		T BE("Backward Euler Matrix", C.nodeID(), C.branchID());
		for (auto i = 1; i < size; ++i)
		{
			for (auto j = 1; j < size; ++j)
			{
				BE(i, j) = C(i, j) + (G(i, j) / _step);
			}
		}

		if(iter)
		{
			spd ? C.CG(x, b, itol) : C.BiCG(x, b, itol);
		}
		else
		{
			spd ? BE.Cholesky() : BE.LU();
			spd ? C.Cholesky() : C.LU();
			spd ? C.solveSPD(x, b) : C.solve(x, b);
		}
		for(auto i = 0; i < _node_vector.size(); ++i)
			oFile_v[i] << 0 << " " << x[_nodes[i]] << std::endl;
		for (double s = _step ; s <= _stop; s += _step)
		{
			for(auto &it : TransientComponentVector)
			{
				if (it.type() == 'V')
					b[it.branch() + C.nodeID()] = it.value(s);
				else
				{
					b[it.plus ()] += it.value(s - _step) - it.value(s);
					b[it.minus()] -= it.value(s - _step) - it.value(s);
				}
			}
			for (auto i = 1; i < size; ++i)
			{
				sum = 0.0;
				for (auto j = 1; j < size; ++j)
				{
					sum += x[j] * G(i,j);
				}
				e[i] = b[i] + (sum/_step);
			}
			if(iter)
			{
				spd ? BE.CG(x, e, itol) : BE.BiCG(x, e, itol);
			}
			else
			{
				for (auto j = 0; j < size; ++j) x[j]=0;
				spd ? BE.solveSPD(x, e) : BE.solve(x, e);
			}
			for(auto i = 0; i < _node_vector.size(); ++i)
				oFile_v[i] << s << " " << x[_nodes[i]] << std::endl;
		}
	}
	else
	{
		T TR("Trapezoidal Matrix", C.nodeID(), C.branchID());
		double * b_prev = new double[size];
		for (auto i = 0; i < size; ++i)
		{
			for (auto j = 0; j < size; ++j)
			{
				TR(i, j) = C(i, j) + (2 * G(i, j) / _step);
			}
		}
		T DC("MNA DC Matrix", C.nodeID(), C.branchID());
		for (auto i = 0; i < size; ++i)
			for (auto j = 0; j < size; ++j)
				DC(i, j) = C(i, j);

		if(iter)
		{
			spd ? C.CG(x, b, itol) : C.BiCG(x, b, itol);
		}
		else
		{
			spd ? TR.Cholesky() : TR.LU();
			spd ? C.Cholesky() : C.LU();
			spd ? C.solveSPD(x, b) : C.solve(x, b);
		}
		for (auto i = 0; i < size; ++i) b_prev[i] = b[i];
		for (auto i = 0; i < _node_vector.size(); ++i)
				oFile_v[i] << 0 << " " << x[_nodes[i]] << std::endl;
		for (double s = _step ; s <= _stop; s += _step)
		{
			for(auto &it : TransientComponentVector)
			{
				if (it.type() == 'V')
					b[it.branch() + C.nodeID()] = it.value(s);
				else
				{
					b[it.plus ()] += it.value(s - _step) - it.value(s);
					b[it.minus()] -= it.value(s - _step) - it.value(s);
				}
			}

			for (auto i = 1; i < size; ++i)
			{
				sum = 0.0;
				for (auto j = 1; j < size; ++j)
				{
					sum += x[j] * (DC(i,j) - (2.0 * G(i,j)/_step));
				}
				e[i] = b[i] + b_prev[i] - sum;
			}
			if(iter)
			{
				spd ? TR.CG(x, e, itol) : TR.BiCG(x, e, itol);
			}
			else
			{
				for (auto j = 0; j < size; ++j) x[j]=0;
				spd ? TR.solveSPD(x, e) : TR.solve(x, e);
			}
			for (auto i = 0; i < size; ++i) b_prev[i] = b[i];
			for (auto i = 0; i < _node_vector.size(); ++i)
				oFile_v[i] << s << " " << x[_nodes[i]] << std::endl;
		}
		delete [] b_prev;
	}
	delete [] e, oFile_v;
}

#define Ci(k) G.matrix()->i[k]
#define Cp(k) G.matrix()->p[k]
#define Cx(k) G.matrix()->x[k]

#define Gi(k) BE.matrix()->i[k]
#define Gp(k) BE.matrix()->p[k]
#define Gx(k) BE.matrix()->x[k]

template <class T>
void TransientAnalysis::sp_analyse(T& C, T& G, double * x, double * b, double itol, bool iter, bool spd, bool method, std::vector<TransientComponent>& TransientComponentVector)
{
	unsigned int i, j, size = C.size();
	double * e = new double[size]; for (auto j = 0; j < size; ++j) e[j]=0;
	std::ofstream * oFile_v = new std::ofstream[_node_vector.size()];
	double sum;
	
	for(auto &it : TransientComponentVector)
	{
		if (it.type() == 'V')
			if((it.branch() + C.nodeID()) != 0)
			b[it.branch() + C.nodeID() - 1] = it.value(0);
		else
		{
			if(it.plus () != 0)
				b[it.plus () - 1] += it.dc_value() - it.value(0);
			if(it.minus () != 0)
				b[it.minus() - 1] -= it.dc_value() - it.value(0);
		}
	}
	
	for(i = 0; i < _node_vector.size(); ++i)
	{
		oFile_v[i].open("Node_" + _node_vector[i] + "_Transient_" + std::to_string(_step) + "_" + std::to_string(_stop) + ".txt");
	}

	if(method)
	{	
		T BE("Backward Euler Matrix", C.nodeID(), C.branchID(), cs_add(C.matrix(), G.matrix(), 1.0, 1.0/_step));

		if(iter)
		{
			spd ? C.CG(x, b, itol) : C.BiCG(x, b, itol);
		}
		else
		{
			spd ? BE.Cholesky() : BE.LU();
			spd ? C.Cholesky() : C.LU();
			spd ? C.solveSPD(x, b) : C.solve(x, b);
		}
		
		for(auto i = 0; i < _node_vector.size(); ++i)
			oFile_v[i] << 0 << " " << x[_nodes[i]-1] << std::endl;

		for (double s = _step ; s <= _stop; s += _step)
		{
			for(auto &it : TransientComponentVector)
			{
				if (it.type() == 'V')
					if((it.branch() + C.nodeID()) != 0)
						b[it.branch() + C.nodeID() - 1] = it.value(s);
				else
				{
					if(it.plus () != 0)
						b[it.plus () - 1] += it.value(s - _step) - it.value(s);				
					if(it.minus () != 0)
						b[it.minus() - 1] -= it.value(s - _step) - it.value(s);
				}
			}

			for(j = 0; j < size; ++j)
			{			
				sum = 0.0;	
				for(i = Cp(j); i < Cp(j+1); ++i)
					sum += Cx(i) * x[Ci(i)];
				e[j] = b[j] + (sum/_step);
			}

			if(iter)
			{
				spd ? BE.CG(x, e, itol) : BE.BiCG(x, e, itol);
			}
			else
			{
				for (auto j = 0; j < size; ++j) x[j]=0;
				spd ? BE.solveSPD(x, e) : BE.solve(x, e);
			}

			for(auto i = 0; i < _node_vector.size(); ++i)
				oFile_v[i] << s << " " << x[(_nodes[i]-1)] << std::endl;
		}
	}
	else
	{

		T TR("Trapezoidal Matrix", C.nodeID(), C.branchID(), cs_add(C.matrix(), G.matrix(), 1.0, 2.0/_step));
		T BE("Temp Matrix", C.nodeID(), C.branchID(), cs_add(C.matrix(), G.matrix(), 1.0, -2.0/_step));

		double * b_prev = new double[size];

		if(iter)
		{
			spd ? C.CG(x, b, itol) : C.BiCG(x, b, itol);
		}
		else
		{
			spd ? TR.Cholesky() : TR.LU();
			spd ? C.Cholesky() : C.LU();
			spd ? C.solveSPD(x, b) : C.solve(x, b);
		}
		for (auto i = 0; i < size; ++i) b_prev[i] = b[i];
		for (auto i = 0; i < _node_vector.size(); ++i)
				oFile_v[i] << 0 << " " << x[_nodes[i]-1] << std::endl;
		for (double s = _step ; s <= _stop; s += _step)
		{

			for(auto &it : TransientComponentVector)
			{
				if (it.type() == 'V')
					if((it.branch() + C.nodeID()) != 0)
						b[it.branch() + C.nodeID() - 1] = it.value(s);
				else
				{
					if(it.plus () != 0)
						b[it.plus () - 1] += it.value(s - _step) - it.value(s);				
					if(it.minus () != 0)
						b[it.minus() - 1] -= it.value(s - _step) - it.value(s);
				}
			}

			for(j = 0; j < size; ++j)
			{			
				sum = 0.0;	
				for(i = Gp(j); i < Gp(j+1); ++i)
					sum += Gx(i) * x[Gi(i)];
				e[j] = b[j] + b_prev[j] - sum;
			}

			if(iter)
			{
				spd ? TR.CG(x, e, itol) : TR.BiCG(x, e, itol);
			}
			else
			{
				for (auto j = 0; j < size; ++j) x[j]=0;
				spd ? TR.solveSPD(x, e) : TR.solve(x, e);
			}
			for (auto i = 0; i < size; ++i) b_prev[i] = b[i];
			for (auto i = 0; i < _node_vector.size(); ++i)
				oFile_v[i] << s << " " << x[_nodes[i]-1] << std::endl;
		}
		delete [] b_prev;
	}
	delete [] e, oFile_v;
}