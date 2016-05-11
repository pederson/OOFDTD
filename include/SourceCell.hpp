#ifndef _SOURCECELL_H
#define _SOURCECELL_H

#include "YeeCell.hpp"
#include "SignalGenerator.hpp"


class SourceCell : public YeeCell{
public:
	SourceCell(SignalGenerator * sg, Field f);

	// updaters for source cell should do nothing
	void update_Dx(unsigned int nstep, double dt, double dx);
	void update_Dy(unsigned int nstep, double dt, double dx);
	void update_Dz(unsigned int nstep, double dt, double dx);
	void update_Bx(unsigned int nstep, double dt, double dx);
	void update_By(unsigned int nstep, double dt, double dx);
	void update_Bz(unsigned int nstep, double dt, double dx);

	void print_summary(ostream & os=std::cout) const;

private:

	SignalGenerator * signal;
	Field fsrc;
};

#endif