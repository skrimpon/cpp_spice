/*
 * Panagiotis Skrimponis
 */

#include "cs_dc_sweep.h"

DCSweepAnalysis::DCSweepAnalysis(char type, std::string name, double start, double step, double stop):
_type(type), _name(name), _start(start), _step(step), _stop(stop)
{

}

DCSweepAnalysis::~DCSweepAnalysis()
{
	delete [] _nodes;
}

void DCSweepAnalysis::update(int plus, int minus, int branch, double value, std::unordered_map < std::string, int >&nodeHash)
{
	_plus = plus; _minus = minus; _branch = branch; _value = value;
	unsigned int i = 0;
	_nodes = new int[_node_vector.size()];
	for(auto& it : _node_vector)
	{
		_nodes[i++] = nodeHash[it];
	}
}
