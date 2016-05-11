#ifndef _GHOSTCELL_H
#define _GHOSTCELL_H

#include "YeeCell.hpp"
#include "mpitools.hpp"


class GhostCell : public YeeCell{
public:
	GhostCell(double & buff);

	// updaters for boundary cell should do nothing
	void update_Dx(unsigned int nstep, double dt, double dx);
	void update_Dy(unsigned int nstep, double dt, double dx);
	void update_Dz(unsigned int nstep, double dt, double dx);
	void update_Bx(unsigned int nstep, double dt, double dx);
	void update_By(unsigned int nstep, double dt, double dx);
	void update_Bz(unsigned int nstep, double dt, double dx);

	void fetch_data(Field f);

private:

	double * buffer_out;
};

#endif