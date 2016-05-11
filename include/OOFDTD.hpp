#ifndef _OOFDTD_H
#define _OOFDTD_H

#include "FDTDCommon.hpp"
#include "RegularMesh.hpp"
#include "RegularDomain.hpp"
#include "YeeCell.hpp"
#include "PECCell.hpp"
#include "SourceCell.hpp"
#include "GhostCell.hpp"

#include "GeometricObject.hpp"
#include "SimulationData.hpp"
//#include "SimulationDataHDF.hpp"
#include "Converter.hpp"
#include "mpitools.hpp"
#include "Timer.hpp"

using namespace std;

class OOFDTD{
public:

	OOFDTD(ParametricModel3D * model, RegularDomain domain,
		   double dx);
	
	// setters for the simulation
	void set_nsteps(unsigned int steps);
	void set_sources(RealPoint rp, Field f, SignalGenerator * sg);
	void capture_data(vector<unsigned int> timesteps, vector<Field> f);

	// running the simulation
	void run();


private:
	void update(Field f);
	void communicate(Field f);

	double delta_t, delta_x;
	double CFL=0.9;
	unsigned int nsteps;
	unsigned int tstep;

	RegularDomain dom;
	RegularMesh mesh;

	// data/output capture
	ofstream ofile;
	vector<Field> capture_fields;
	vector<unsigned int> capture_steps;
	SimulationData * simdata;

	// boundaries
	map<BoundaryLocation, BoundaryType> boundaries;
	map<BoundaryLocation, double *> boundary_buffer_out;
	map<BoundaryLocation, double *> boundary_buffer_in;

	// cells owned by this proc
	vector<YeeCell *> inner_cells;
	map<BoundaryLocation, vector<YeeCell *>> boundary_cells;
};


#endif