#include "GhostCell.hpp"

GhostCell::GhostCell(double & buff){
	buffer_out = &buff;
	E_x = 0.0;
	E_y = 0.0;
	E_z = 0.0;
	H_x = 0.0;
	H_y = 0.0;
	H_z = 0.0;
	D_x = 0.0;
	D_y = 0.0;
	D_z = 0.0;
	B_x = 0.0;
	B_y = 0.0;
	B_z = 0.0;
	permittivity = 1.0;
	permeability = 1.0;
}

// updaters for boundary cell should do nothing
void GhostCell::update_Dx(unsigned int nstep, double dt, double dx){

}

void GhostCell::update_Dy(unsigned int nstep, double dt, double dx){

}

void GhostCell::update_Dz(unsigned int nstep, double dt, double dx){

}

void GhostCell::update_Bx(unsigned int nstep, double dt, double dx){

}

void GhostCell::update_By(unsigned int nstep, double dt, double dx){

}

void GhostCell::update_Bz(unsigned int nstep, double dt, double dx){

}

void GhostCell::fetch_data(Field f){
	if (f == Field::E_X) E_x = *buffer_out;
	else if (f == Field::E_Y) E_y = *buffer_out;
	else if (f == Field::E_Z) E_z = *buffer_out;
	else if (f == Field::H_X) H_x = *buffer_out;
	else if (f == Field::H_Y) H_y = *buffer_out;
	else if (f == Field::H_Z) H_z = *buffer_out;
}

