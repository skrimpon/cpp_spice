/*
 * Panagiotis Skrimponis
 */

#pragma once
#include "cs_headerdef.h"

class Component
{
public:
	Component(char, std::string, int, int, double);	// Class Constructor.
	~Component();
	inline void set_branch(int);					// Set Component branch.
	inline char type() const;						// Query Component type.
	inline std::string name() const;				// Query Component name.
	inline int plus() const;						// Query Component <+> node.
	inline int minus() const;						// Query Component <-> node.
	inline int branch() const;						// Query Component branch.
	inline double value() const;					// Query Component value.
private:
	char _type;
	std::string _name;
	int _branch, _plus_node, _minus_node;
	double _value;
};

inline void Component::set_branch(int branch) {
	_branch = branch;
}

inline char Component::type() const {
	return _type;
}

inline std::string Component::name() const {
	return _name;
}

inline int Component::plus() const {
	return _plus_node;
}

inline int Component::minus() const {
	return _minus_node;
}

inline double Component::value() const {
	return _value;
}

inline int Component::branch() const {
	return _branch;
}