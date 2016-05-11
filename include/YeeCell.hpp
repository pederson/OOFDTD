#ifndef _YEECELL_H
#define _YEECELL_H


#include <map>
#include "RegularMesh.hpp"
#include "FDTDCommon.hpp"

using namespace std;

class YeeCell{
public:

	YeeCell();

	// update functions for the Yee algorithm
	virtual void update_Dx(unsigned int nstep, double dt, double dx);
	virtual void update_Dy(unsigned int nstep, double dt, double dx);
	virtual void update_Dz(unsigned int nstep, double dt, double dx);

	virtual void update_Bx(unsigned int nstep, double dt, double dx);
	virtual void update_By(unsigned int nstep, double dt, double dx);
	virtual void update_Bz(unsigned int nstep, double dt, double dx);

	// update using permittivity and permeability
	virtual void update_Ex(unsigned int nstep, double dt, double dx);
	virtual void update_Ey(unsigned int nstep, double dt, double dx);
	virtual void update_Ez(unsigned int nstep, double dt, double dx);
	virtual void update_Hx(unsigned int nstep, double dt, double dx);
	virtual void update_Hy(unsigned int nstep, double dt, double dx);
	virtual void update_Hz(unsigned int nstep, double dt, double dx);

	// accessors
	double Ex() const {return E_x;};
	double Ey() const {return E_y;};
	double Ez() const {return E_z;};
	double Hx() const {return H_x;};
	double Hy() const {return H_y;};
	double Hz() const {return H_z;};

	void set_neighbor(BoundaryLocation bl, YeeCell * neighbor);
	void set_properties(double permittiv, double permeabil);

	virtual void print_summary(ostream & os=std::cout) const;
	unsigned int num_neighbors() const;
	map<BoundaryLocation, YeeCell *> get_neighbors();

protected:
	// neighboring cells
	// (can throw error if there is no neighbor!!)
	map<BoundaryLocation, YeeCell *> neighbors;

	double E_x, E_y, E_z;
	double H_x, H_y, H_z;
	double D_x, D_y, D_z;
	double B_x, B_y, B_z;

	double permittivity, permeability;

};

#endif