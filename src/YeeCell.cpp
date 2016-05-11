#include "YeeCell.hpp"

YeeCell::YeeCell(){
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

// update functions for the Yee algorithm
void YeeCell::update_Dx(unsigned int nstep, double dt, double dx){
	D_x += dt*( (H_z - neighbors[BoundaryLocation::Y_MIN]->H_z)/dx 
			  - (H_y - neighbors[BoundaryLocation::Z_MIN]->H_y)/dx );
}

void YeeCell::update_Dy(unsigned int nstep, double dt, double dx){
	D_y += dt*( (H_x - neighbors[BoundaryLocation::Z_MIN]->H_x)/dx 
			  - (H_z - neighbors[BoundaryLocation::X_MIN]->H_z)/dx );
}

void YeeCell::update_Dz(unsigned int nstep, double dt, double dx){
	D_z += dt*( (H_y - neighbors[BoundaryLocation::X_MIN]->H_y)/dx 
			  - (H_x - neighbors[BoundaryLocation::Y_MIN]->H_x)/dx );
}

void YeeCell::update_Bx(unsigned int nstep, double dt, double dx){
	B_x += -dt*( (neighbors[BoundaryLocation::Y_MAX]->E_z - E_z)/dx 
			  - (neighbors[BoundaryLocation::Z_MAX]->E_y - E_y)/dx );
}

void YeeCell::update_By(unsigned int nstep, double dt, double dx){
	B_y += -dt*( (neighbors[BoundaryLocation::Z_MAX]->E_x - E_x)/dx 
			  - (neighbors[BoundaryLocation::X_MAX]->E_z - E_z)/dx );
}

void YeeCell::update_Bz(unsigned int nstep, double dt, double dx){
	B_z += -dt*( (neighbors[BoundaryLocation::X_MAX]->E_y - E_y)/dx 
			  - (neighbors[BoundaryLocation::Y_MAX]->E_x - E_x)/dx );
}


// update using permittivity and permeability
void YeeCell::update_Ex(unsigned int nstep, double dt, double dx){
	E_x = D_x/(permittivity*m_eps0);
}

void YeeCell::update_Ey(unsigned int nstep, double dt, double dx){
	E_y = D_y/(permittivity*m_eps0);
}

void YeeCell::update_Ez(unsigned int nstep, double dt, double dx){
	E_z = D_z/(permittivity*m_eps0);
}

void YeeCell::update_Hx(unsigned int nstep, double dt, double dx){
	H_x = B_x/(permeability*m_mu0);
}

void YeeCell::update_Hy(unsigned int nstep, double dt, double dx){
	H_y = B_y/(permeability*m_mu0);
}

void YeeCell::update_Hz(unsigned int nstep, double dt, double dx){
	H_z = D_z/(permeability*m_mu0);
}


void YeeCell::set_neighbor(BoundaryLocation bl, YeeCell * neighbor){
	neighbors[bl] = neighbor;
}

void YeeCell::set_properties(double permittiv, double permeabil){
	permittivity = permittiv;
	permeability = permeabil;
}

void YeeCell::print_summary(ostream & os) const{
	os << "(eps, mu) = (" << permittivity << ", " << permeability << "), nneighbors: " << neighbors.size() ;
}

unsigned int YeeCell::num_neighbors() const{
	return neighbors.size(); 
}

map<BoundaryLocation, YeeCell *> YeeCell::get_neighbors(){
	return neighbors;
}
