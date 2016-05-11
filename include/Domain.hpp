#ifndef _DOMAIN_H
#define _DOMAIN_H

//#include "RegularMesh.hpp"

class Domain{
public:

	virtual unsigned int array_node_offset(int proc) const = 0;
	virtual unsigned int array_element_offset(int proc) const = 0;
	virtual unsigned int nnodes(int proc) const = 0;
	virtual unsigned int nelements(int proc) const = 0;
	virtual unsigned int nnodes() const = 0;
	virtual unsigned int nelements() const = 0;

};

#endif