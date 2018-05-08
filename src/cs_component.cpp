/*
 * Panagiotis Skrimponis
 */

#include "cs_component.h"
// Class constructor.
Component::Component( char type, std::string name, int plus_node, int minus_node, double value) : 
_type(type), _name(name), _plus_node(plus_node), _minus_node(minus_node), _value(value) {}
// Class destructor.
Component::~Component() {}