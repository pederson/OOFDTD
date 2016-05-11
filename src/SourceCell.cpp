#include "SourceCell.hpp"

using namespace std;

SourceCell::SourceCell(SignalGenerator * sg, Field f){
	signal = sg;
	fsrc = f;
}

// updaters for source cell should do nothing
void SourceCell::update_Dx(unsigned int nstep, double dt, double dx){
	if (fsrc != Field::D_X) return;
	D_x = signal->value(nstep*dt);
}

void SourceCell::update_Dy(unsigned int nstep, double dt, double dx){
	if (fsrc != Field::D_Y) return;
	D_y = signal->value(nstep*dt);
}

void SourceCell::update_Dz(unsigned int nstep, double dt, double dx){
	if (fsrc != Field::D_X) return;
	D_z = signal->value(nstep*dt);
}

void SourceCell::update_Bx(unsigned int nstep, double dt, double dx){
	if (fsrc != Field::B_X) return;
	B_x = signal->value(nstep*dt);
}

void SourceCell::update_By(unsigned int nstep, double dt, double dx){
	if (fsrc != Field::B_Y) return;
	B_y = signal->value(nstep*dt);
}

void SourceCell::update_Bz(unsigned int nstep, double dt, double dx){
	if (fsrc != Field::B_Z) return;
	B_z = signal->value(nstep*dt);
}

void SourceCell::print_summary(ostream & os) const{
	os << "SourceCell: num neighbors: " << neighbors.size() << endl;
	signal->print_summary(); os << endl;
}
