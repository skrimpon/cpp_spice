/*
 * Panagiotis Skrimponis
 */

#include "cs_ac_sweep.h"
/*
 * ACComponent Class
 */
ACComponent::ACComponent(char type, std::string name, double dc_value, double mag, double phase, int plus, int minus, int branch) :
 _name(name), _dc_value(dc_value), _component_plus(plus), _component_minus(minus), _component_branch(branch), _component_type(type)
{
	_component_complex = std::polar(mag, phase * 2 * M_PI / 360.0);
}

ACComponent::~ACComponent()
{

}
/*
 * ACSweepAnalysis Class
 */
 ACSweepAnalysis::ACSweepAnalysis(std::string sweep, double start, double step, double stop): 
_sweep(sweep), _start(start), _step(step), _stop(stop)
{

}

ACSweepAnalysis::~ACSweepAnalysis()
{
	delete [] _nodes;
}

void ACSweepAnalysis::update(std::unordered_map < std::string, int >& nodeHash)
{
	int i = 0;
	_nodes = new int[_node_vector.size()];
	for(auto& it : _node_vector)
	{
		_nodes[i++] = nodeHash[it];
	}
}