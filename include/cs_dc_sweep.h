/*
 * Panagiotis Skrimponis
 */

#pragma once
#include "cs_headerdef.h"

class DCSweepAnalysis
{
public:
	DCSweepAnalysis(char, std::string, double, double, double);
	~DCSweepAnalysis();
	template <class T>
	void analyse(T &, double *, double, bool, bool);
	template <class T>
	void sp_analyse(T &, double *, double, bool, bool);
	inline std::string name() const;
	void update(int, int, int, double, std::unordered_map < std::string, int >&);
	inline void add_node(std::string);
private:
	char _type;
	std::string _name;
	int _plus, _minus, _branch;
	double _start, _step, _stop, _value;
	std::vector < std::string > _node_vector;
	int * _nodes; 
};

inline std::string DCSweepAnalysis::name() const
{
	return _name;
}

inline void DCSweepAnalysis::add_node(std::string node)
{
	_node_vector.push_back(node);
}

template <class T>
void DCSweepAnalysis::analyse(T& MNA_DC, double * MNA_b, double itol, bool iter, bool spd)
{
	double * MNA_x = new double[MNA_DC.size()];
	
	unsigned int i, k;
	double temp;

	std::ofstream * oFile_v = new std::ofstream[ _node_vector.size() ];
	for(i = 0; i < _node_vector.size(); ++i)
	{
		oFile_v[i].open("Node_" + _node_vector[i] + "_DC_Sweep_" + _name + "_" + std::to_string(_start) + "_" + std::to_string(_stop) + "_" + std::to_string(_step) + ".txt");
	}

	switch ( _type ) {
	case('V') :
		k = _branch + MNA_DC.nodeID();
		temp = MNA_b[k];
		do {
			MNA_b[k] = _start;
		    if(iter)
		    {
				spd ? MNA_DC.CG(MNA_x, MNA_b, itol) : MNA_DC.BiCG(MNA_x, MNA_b, itol);
				for(i = 0; i < _node_vector.size(); ++i)
					oFile_v[i] << _start << " " << MNA_x[_nodes[i]] << std::endl;
		    }
		    else
		    {
				for(i = 0; i < MNA_DC.size(); ++i) MNA_x[i] = 0;
				spd ? MNA_DC.solveSPD(MNA_x, MNA_b) : MNA_DC.solve(MNA_x, MNA_b);
				for(i = 0; i < _node_vector.size(); ++i)
					oFile_v[i] << _start << " " << MNA_x[_nodes[i]] << std::endl;
			}
			_start += _step;
		} while( _stop >= _start );
		MNA_b[k] = temp;
		break;
	case('I') :
		double temp_p = MNA_b[_plus];
		double temp_m = MNA_b[_minus];
		MNA_b[_plus] += _value - _start;
		MNA_b[_minus] -= _value - _start;
		do {
		    if(iter)
		    {
				spd ? MNA_DC.CG(MNA_x, MNA_b, itol) : MNA_DC.BiCG(MNA_x, MNA_b, itol);
				for(i = 0; i < _node_vector.size(); ++i)
					oFile_v[i] << _start << " " << MNA_x[_nodes[i]] << std::endl;
		    }
		    else
		    {
				for(i = 0; i < MNA_DC.size(); ++i) MNA_x[i] = 0;
				spd ? MNA_DC.solveSPD(MNA_x, MNA_b) : MNA_DC.solve(MNA_x, MNA_b);
				for(i = 0; i < _node_vector.size(); ++i)
					oFile_v[i] << _start << " " << MNA_x[_nodes[i]] << std::endl;
			}
			_start += _step;
			MNA_b[_plus] -= _step;
			MNA_b[_minus] += _step;
		} while( _stop >= _start );
		MNA_b[_plus] = temp_p;
		MNA_b[_minus] = temp_m;
	}

	for(i = 0; i < _node_vector.size(); ++i)
		oFile_v[i].close();
	delete [] oFile_v, MNA_x;
}

template <class T>
void DCSweepAnalysis::sp_analyse(T& MNA_DC, double * MNA_b, double itol, bool iter, bool spd)
{
	double * MNA_x = new double[MNA_DC.size()];
	
	unsigned int i, k;
	double temp;

	std::ofstream * oFile_v = new std::ofstream[ _node_vector.size() ];
	for(i = 0; i < _node_vector.size(); ++i)
	{
		oFile_v[i].open("Node_" + _node_vector[i] + "_DC_Sweep_" + _name + "_" + std::to_string(_start) + "_" + std::to_string(_stop) + "_" + std::to_string(_step) + ".txt");
	}

	switch ( _type ) {
	case('V') :
		k = _branch + MNA_DC.nodeID() - 1;
		temp = MNA_b[k];
		do {
			MNA_b[k] = _start;
		    if(iter)
		    {
				spd ? MNA_DC.CG(MNA_x, MNA_b, itol) : MNA_DC.BiCG(MNA_x, MNA_b, itol);
				for(i = 0; i < _node_vector.size(); ++i)
					oFile_v[i] << _start << " " << MNA_x[_nodes[i]-1] << std::endl;
		    }
		    else
		    {
				for(i = 0; i < MNA_DC.size(); ++i) MNA_x[i] = 0;
				spd ? MNA_DC.solveSPD(MNA_x, MNA_b) : MNA_DC.solve(MNA_x, MNA_b);
				for(i = 0; i < _node_vector.size(); ++i)
					oFile_v[i] << _start << " " << MNA_x[_nodes[i]-1] << std::endl;
			}
			_start += _step;
		} while( _stop >= _start );
		MNA_b[k] = temp;
		break;
	case('I') :
		double temp_p, temp_m;
		if(_plus != 0)
		{
			temp_p = MNA_b[_plus-1];
			MNA_b[_plus] += _value - _start;
		}
		if(_minus != 0)
		{
			temp_m = MNA_b[_minus-1];
			MNA_b[_minus] -= _value - _start;
		}
		do {
		    if(iter)
		    {
				spd ? MNA_DC.CG(MNA_x, MNA_b, itol) : MNA_DC.BiCG(MNA_x, MNA_b, itol);
				for(i = 0; i < _node_vector.size(); ++i)
					oFile_v[i] << _start << " " << MNA_x[_nodes[i]] << std::endl;
		    }
		    else
		    {
				for(i = 0; i < MNA_DC.size(); ++i) MNA_x[i] = 0;
				spd ? MNA_DC.solveSPD(MNA_x, MNA_b) : MNA_DC.solve(MNA_x, MNA_b);
				for(i = 0; i < _node_vector.size(); ++i)
					oFile_v[i] << _start << " " << MNA_x[_nodes[i]] << std::endl;
			}
			_start += _step;
			if(_plus != 0)
				MNA_b[_plus-1] -= _step;
			if(_minus != 0)
				MNA_b[_minus-1] += _step;
		} while( _stop >= _start );
		if(_plus != 0)
			MNA_b[_plus-1] = temp_p;
		if(_minus != 0)
			MNA_b[_minus-1] = temp_m;
	}

	for(i = 0; i < _node_vector.size(); ++i)
		oFile_v[i].close();
	delete [] oFile_v, MNA_x;
}