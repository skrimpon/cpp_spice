/*
 * Panagiotis Skrimponis
 */

#pragma once
#include "cs_headerdef.h"
/*
 * ACComponent Class
 */
class ACComponent
{
public:
	ACComponent(char, std::string, double, double, double, int, int, int); // Class Constructor.
	~ACComponent();
	inline char type() const;
	inline int plus() const;
	inline int minus() const;
	inline int branch() const;
	inline std::complex<double> value() const;
private:
	char _component_type;
	int _component_plus, _component_minus, _component_branch;
	std::string _name, _type;
	double _dc_value;
	std::complex <double> _component_complex;
};

inline char ACComponent::type() const
{
	return _component_type;
}

inline int ACComponent::plus() const
{
	return _component_plus;
}

inline int ACComponent::minus() const
{
	return _component_minus;
}

inline int ACComponent::branch() const
{
	return _component_branch;
}
inline std::complex<double> ACComponent::value() const
{
	return _component_complex;
}
/*
 * ACSweepAnalysis Class
 */
class ACSweepAnalysis
{
public:
	ACSweepAnalysis(std::string, double, double, double);
	~ACSweepAnalysis();
	inline void add_node(std::string);
	void update(std::unordered_map < std::string, int >&);
	template <class T>
	void analyse(T&, std::complex<double> *, double, bool, bool, std::vector<ACComponent>&);
	template <class T>
	void sp_analyse(T&, std::complex<double> *, double, bool, bool, std::vector<ACComponent>&);
private:
	int * _nodes; 
	std::vector < std::string > _node_vector;
	std::string _sweep;
	double _start, _step, _stop;
};

inline void ACSweepAnalysis::add_node(std::string node)
{
	_node_vector.push_back(node);
}

template <class T>
void ACSweepAnalysis::analyse(T& MNA_AC, std::complex<double> * MNA_b, double itol, bool iter, bool spd, std::vector<ACComponent>& ACComponentVector)
{	
	unsigned int i, j, k, size = MNA_AC.size(), max_step = _step;
	std::ofstream * oFile_v = new std::ofstream[ _node_vector.size() ];
	std::complex<double> * MNA_x = new std::complex<double>[size];
	for(i=0; i < MNA_AC.size(); i++) MNA_x[i] = std::complex<double>(0.0, 0.0);

	if (_sweep == "lin")
		_step = (_stop - _start)/(_step - 1.0);
	else {
		_step = log(10)/_step;
		std::cout << "LOG\n";
	}

	for(i = 0; i < _node_vector.size(); ++i)
		oFile_v[i].open("Node_" + _node_vector[i] + "_AC_Sweep_" + std::to_string(_step) + "_" + std::to_string(_stop) + ".txt");

	for (double s = _start, k = 0; k < max_step; k++)
	{
		T ABS("MNA ABS Matrix", MNA_AC.nodeID(), MNA_AC.branchID());

		for(i=1; i < MNA_AC.size(); i++)
			for(j=1; j < MNA_AC.size(); j++)
				ABS(i,j) = std::complex<double>(MNA_AC(i,j).real(), s * 2 * M_PI * MNA_AC(i,j).imag());
		
		if(iter)
	    {
			ABS.BiCG(MNA_x, MNA_b, itol);
			for(i = 0; i < _node_vector.size(); ++i)
				oFile_v[i] << s << " " << MNA_x[_nodes[i]] << std::endl;
	    }
	    else
	    {
			for(i = 0; i < ABS.size(); ++i) MNA_x[i] = std::complex<double>(0.0, 0.0);
			ABS.LU(); ABS.solve(MNA_x, MNA_b);
			for(i = 0; i < _node_vector.size(); ++i)
				oFile_v[i] << s << " " << MNA_x[_nodes[i]] << std::endl;
		}
		if (_sweep == "lin")
			s += _step;
		else
			s *= exp(_step);
	}
	for(i = 0; i < _node_vector.size(); ++i)
		oFile_v[i].close();
}

#define Ci(k) MNA_AC.matrix()->i[k]
#define Cp(k) MNA_AC.matrix()->p[k]
#define Cx(k) MNA_AC.matrix()->x[k]

template <class T>
void ACSweepAnalysis::sp_analyse(T& MNA_AC, std::complex<double> * MNA_b, double itol, bool iter, bool spd, std::vector<ACComponent>& ACComponentVector)
{	
	unsigned int i, j, k, size = MNA_AC.size(), max_step = _step;

	std::ofstream * oFile_v = new std::ofstream[ _node_vector.size() ];

	std::complex<double> * MNA_x = new std::complex<double>[size];

	for(i=0; i < MNA_AC.size(); i++) MNA_x[i] = std::complex<double>(0.0, 0.0);

	_step = (_stop - _start)/(_step - 1.0);

	for(i = 0; i < _node_vector.size(); ++i)
	{
		oFile_v[i].open("Node_" + _node_vector[i] + "_AC_Sweep_" + std::to_string(_step) + "_" + std::to_string(_stop) + ".txt");
	}

	cs_ci * csi_miami = cs_ci_add(MNA_AC.matrix(), MNA_AC.matrix(), 0.0, 1.0);

	T ABS("MNA ABS Matrix", MNA_AC.nodeID(), MNA_AC.branchID(), csi_miami);

	for (double s = _start, k = 0; k < max_step; k++, s += _step)
	{
		for(j = 0; j < MNA_AC.size(); ++j)
			for(i = Cp(j); i < Cp(j+1); ++i)
			{
				ABS.matrix()->i[i] = Ci(i);
				ABS.matrix()->p[i] = Cp(i);
				ABS.matrix()->x[i] = std::complex<double>(Cx(i).real(), s * 2 * M_PI * Cx(i).imag());
			}

		if(iter)
	    {
			for(i = 0; i < ABS.size(); ++i) MNA_x[i] = std::complex<double>(0.0, 0.0);
			ABS.BiCG(MNA_x, MNA_b, itol);
			for(i = 0; i < _node_vector.size(); ++i)
				oFile_v[i] << s << " " << MNA_x[_nodes[i]-1] << std::endl;
	    }
	    else
	    {
			for(i = 0; i < ABS.size(); ++i) MNA_x[i] = std::complex<double>(0.0, 0.0);
			ABS.LU();
			ABS.solve(MNA_x, MNA_b);
			for(i = 0; i < _node_vector.size(); ++i)
				oFile_v[i] << s << " " << MNA_x[_nodes[i]-1] << std::endl;
		}
	}
	for(i = 0; i < _node_vector.size(); ++i)
		oFile_v[i].close();
}