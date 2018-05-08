/*
 * Panagiotis Skrimponis
 */

#include "cs_transient.h"
/*
 * TransientComponent Class
 */
TransientComponent::TransientComponent(char type, std::string name, std::stringstream& iLine, double dc_value, int plus, int minus, int branch) :
 _name(name), _dc_value(dc_value), _component_plus(plus), _component_minus(minus), _component_branch(branch), _component_type(type)
{
	double t, v;
	iLine >> _type;
	if (_type == "EXP") iLine >> _i1 >> _i2 >> _td1 >> _tc1 >> _td2 >> _tc2;
	else if (_type == "PULSE") iLine >> _i1 >> _i2 >> _td >> _tr >> _tf >> _pw >> _per;
	else if (_type == "SIN") iLine >> _i1 >> _ia >> _fr >> _td >> _df >> _ph;
	else // PWL
	{
		iLine >> t >> v;
		while (!iLine.eof())
		{
			_pwl_time.push_back(t);
			_pwl_value.push_back(v);
			iLine >> t >> v;
		}
	}
}

TransientComponent::~TransientComponent()
{

}

double TransientComponent::value(double t) const
{
	if (_type == "EXP")
	{
		if ( t < _td1 )							return _i1;
		else if ( t < _td2 )					return (_i1 + (_i2 - _i1) * (1 - exp((_td1-t)/_tc1)));
		else 									return (_i1 + ((_i2 - _i1) * (exp((_td2-t)/_tc2))) - exp((_td1-t)/_tc1));	
	}
	else if (_type == "SIN")
	{
		if( t <= _td )							return (_i1 + _ia * sin((2 * M_PI * _ph)/360));
		else									return (_i1 + _ia * sin((2 * M_PI * _fr * (t - _td)) + ((2 * M_PI *_ph)/360)) * exp((_td-t)/_df));
	}
	else if (_type == "PULSE")
	{
		t = fmod(t, _per);
		if ( t < _td )							return _i1;
		else if ( t < (_td + _tr) )				return _i1 + (((_i2 - _i1)/ _tr) * (t - _td));
		else if ( t < (_td + _tr + _pw) )		return _i2;
		else if ( t < (_td + _tr + _pw + _tf) )	return _i2 + (((_i1 - _i2)/ _tf) * (t - (_td + _tr + _pw)));
		else									return _i1;
	}
	else // PWL
	{
		if( _pwl_time[0] > t )					return _pwl_value[0];
		for(auto i = 1; i < _pwl_time.size(); ++i)
			if(_pwl_time[i] > t)				return ((_pwl_value[i] - _pwl_value[i-1])/(_pwl_time[i]-_pwl_time[i-1])*(t-_pwl_time[i-1])+_pwl_value[i-1]);
												return _pwl_value[_pwl_time.size()-1];
	}
}
/*
 * TransientAnalysis Class
 */
TransientAnalysis::TransientAnalysis(double step, double stop): _step(step), _stop(stop)
{

}

TransientAnalysis::~TransientAnalysis()
{
	delete [] _nodes;
}

void TransientAnalysis::update(std::unordered_map < std::string, int >& nodeHash)
{
	unsigned int i = 0;
	_nodes = new int[_node_vector.size()];
	for(auto& it : _node_vector)
	{
		_nodes[i++] = nodeHash[it];
	}
}