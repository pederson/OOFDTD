#ifndef _PECCELL_H
#define _PECCELL_H

#include "YeeCell.hpp"

class PECCell : public YeeCell{
public:

	PECCell();

	// updaters for boundary cell should do nothing
	void update_Dx(unsigned int nstep, double dt, double dx);
	void update_Dy(unsigned int nstep, double dt, double dx);
	void update_Dz(unsigned int nstep, double dt, double dx);
	void update_Bx(unsigned int nstep, double dt, double dx);
	void update_By(unsigned int nstep, double dt, double dx);
	void update_Bz(unsigned int nstep, double dt, double dx);

	void print_summary() const;

private:

};

#endif