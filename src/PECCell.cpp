#include "PECCell.hpp"


PECCell::PECCell(){

}

// updaters for boundary cell should do nothing
void PECCell::update_Dx(unsigned int nstep, double dt, double dx){};

void PECCell::update_Dy(unsigned int nstep, double dt, double dx){};

void PECCell::update_Dz(unsigned int nstep, double dt, double dx){};

void PECCell::update_Bx(unsigned int nstep, double dt, double dx){};

void PECCell::update_By(unsigned int nstep, double dt, double dx){};

void PECCell::update_Bz(unsigned int nstep, double dt, double dx){};

void PECCell::print_summary() const{
	cout << "PECCell " ;
}
