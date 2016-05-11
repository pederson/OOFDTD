#include "include/OOFDTD.hpp"
#include "include/mpitools.hpp"
#include "include/GeometricObject.hpp"
#include <iostream>


int main(int argc, char * argv[]){
	int ierr, num_procs, my_id;

	//int n_per_proc = atoi(argv[1]);
	mpi::init(&argc, &argv);
	my_id = mpi::rank();
	num_procs = mpi::size();


	// parameters
	double res = 50.0;			// number of cells per wavelength
	double srcfreq = 50.0e+9;
	double freqspread = 10.0e+9;
	double wavelen = 3.0e+8/srcfreq;
	double lx = 2.0*wavelen;
	double ly = 2.0*wavelen;			// domain length in y direction
	double lz = 2.0*wavelen;
	double dx = wavelen/res;
	unsigned int numcycles = 8;
	double srclocx = 0.0;
	double srclocy = 0.0;
	double srclocz = 0.0;
	//
	//
	//
	
	// calculate the resolution given the number of cells per proc
	//res = pow(pow(wavelen/lx, 3.0)*double(n_per_proc)*double(mpi::size()) , 1.0/3.0);
	//dx = wavelen/res;

	if (mpi::is_master()){
		cout << "freq(Hz): " << srcfreq << " wavelen(m): " << wavelen << endl;
		cout << "lx: " << lx << " ly: " << ly << endl;
		cout << "dx: " << dx << endl;
		cout << "nnodesx: " << lx/dx << endl;
	}

	// 3D
	
	ParametricModel3D paramodel;
	paramodel.set_model_name("Empty");
	paramodel.add_physical_property("m_index");
	paramodel.add_physical_property("metal_nodes");
	paramodel.add_material("Vacuum", {0.0, 0.0});
	paramodel.add_material("Metal", {1.0, 1.0});
	if (mpi::is_master()) paramodel.print_summary();
	//*/

	// define the simulation domain
	RealBounds rb = RealBounds(RealPoint(-lx/2.0, -ly/2.0, -lz/2.0), RealPoint(lx/2.0, ly/2.0, lz/2.0));
	RegularDomain dom = RegularDomain(rb, RealPoint(0,0,0), dx, dx, dx, mpi::size());

	// cout << "ABOUT TO CREATE SIMULATION OBJECT" << endl;
	// create and run the simulation
	OOFDTD sim(&paramodel, dom, dx);
	sim.set_nsteps(100);
	SignalSinusoid sig(srcfreq);
	//sim.set_sources(RealPoint(dx, dx, dx), Field::D_Z, &sig);
	// cout << "ABOUT TO RUN" << endl;
	sim.run();



	//cout << mpi::rank() << " FINISHED EVERYTHING... SHOULD BE EXITING NOW" << endl;
	mpi::finalize();

	//cout << "mpi finalized... about to return" << endl;
	return 0;
}
