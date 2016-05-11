#include "OOFDTD.hpp"

using namespace std;

OOFDTD::OOFDTD(ParametricModel3D * model, RegularDomain domain,
			   double dx){

	stringstream ss;
	ss << "proc" << mpi::rank() << ".log" ;
	ofile.open(ss.str());

	ofile << " AM HERE " << endl;
	// create a mesh from the given domain information
	// and from the model
	delta_x = dx;
	dom = domain;
	RealBounds rb = dom.get_proc_bounds(mpi::rank());
	ofile << "GOT PROC BOUNDS" << endl;
	mesh = build_simple_mesh_3d(*model, dx, rb, {0,1});
	ofile << "BUILT MESH" << endl;
	mesh.print_summary(ofile);

	if (mpi::is_master()){
		cout << "The proc grid looks like this: " << endl;
		dom.print_grid();
	}

	// set the time step using the CFL criterion
	delta_t = CFL*delta_x/m_c0;

	// allocate YeeCell objects for the interior of the domain
	inner_cells.resize(mesh.elementcount());
	for (auto i=0 ; i<mesh.elementcount(); i++){
		inner_cells[i] = new YeeCell();
	}
	ofile << "ALLOCATED INTERIOR CELLS (" << mesh.elementcount() << " cells)" << endl;

	// allocate boundary cells
	vector<BoundaryLocation> bloc = get_boundaries(EMMode::ThreeD);
	for (auto bl=bloc.begin(); bl!=bloc.end(); bl++){
		if (dom.is_proc_on_boundary(mpi::rank(), *bl)) boundaries[*bl] = BoundaryType::PEC;
		else boundaries[*bl] = BoundaryType::PARALLEL;
	}
	for (auto b=boundaries.begin(); b!=boundaries.end(); b++){
		// get number of cells on this boundary and resize vector
		IndexBounds ib = dom.get_proc_local_ibounds_cell(mpi::rank());
		ib.print_summary(ofile);
		IndexPlane ip = ib.get_boundary_plane(b->first);
		
		//ip.print_summary(ofile);
		ofile << "the boundary plane " << get_string(b->first) << " has " << ip.get_npts() << " cells" << endl;
		vector<YeeCell *> bcells (ip.get_npts());
		boundary_cells[b->first] = bcells;
		if (b->second == BoundaryType::PARALLEL){
			unsigned int neighbor = dom.get_proc_neighbor(mpi::rank(), b->first);
			boundary_buffer_out[b->first] = new double[bcells.size()];
			boundary_buffer_in[b->first] = new double[bcells.size()];
			unsigned int ct=0;
			for (auto c=boundary_cells[b->first].begin(); c!=boundary_cells[b->first].end(); c++){
				(*c) = new GhostCell(boundary_buffer_in[b->first][ct]);\
				ct++;
			}
		}
		else if (b->second == BoundaryType::PEC){
			for (auto c=boundary_cells[b->first].begin(); c!=boundary_cells[b->first].end(); c++){
				(*c) = new PECCell();
			}
		}
	}	
 
	for (auto b=boundaries.begin(); b!=boundaries.end(); b++){
		ofile << "boundary " << get_string(b->first) << " is of type " << get_string(b->second) << endl;
	}

	ofile << "ALLOCATED BOUNDARY CELLS" << endl;

	// link the interior YeeCells to their neighbors
	unsigned int ind;
	unsigned int indxm, indxp, indym, indyp, indzm, indzp;
	for (auto i=1; i<mesh.num_cells_x()-1; i++){
		for (auto j=1; j<mesh.num_cells_y()-1; j++){
			for (auto k=1; k<mesh.num_cells_z()-1; k++){
				// mapping to the global index
				ind = mesh.cell_map_to_global_ind(i,j,k);
				indxm = mesh.cell_map_to_global_ind(i-1,j,k);
				indxp = mesh.cell_map_to_global_ind(i+1,j,k);
				indym = mesh.cell_map_to_global_ind(i,j-1,k);
				indyp = mesh.cell_map_to_global_ind(i,j+1,k);
				indzm = mesh.cell_map_to_global_ind(i,j,k-1);
				indzp = mesh.cell_map_to_global_ind(i,j,k+1);


				inner_cells[ind]->set_neighbor(BoundaryLocation::X_MIN, inner_cells[indxm]);
				inner_cells[ind]->set_neighbor(BoundaryLocation::X_MAX, inner_cells[indxp]);
				inner_cells[ind]->set_neighbor(BoundaryLocation::Y_MIN, inner_cells[indym]);
				inner_cells[ind]->set_neighbor(BoundaryLocation::Y_MAX, inner_cells[indyp]);
				inner_cells[ind]->set_neighbor(BoundaryLocation::Z_MIN, inner_cells[indzm]);
				inner_cells[ind]->set_neighbor(BoundaryLocation::Z_MAX, inner_cells[indzp]);
			}
		}
	}

	ofile << "LINKED INTERIOR CELLS" << endl;

	// link the boundary YeeCells with their boundary neighbors
	IndexBounds ib = dom.get_proc_local_ibounds_cell(mpi::rank());
	ib.print_summary(ofile);
	for (auto b=boundaries.begin(); b!=boundaries.end(); b++){
		// get the plane of indices on this boundary
		IndexPlane ip = ib.get_boundary_plane(b->first);

		// get the boundary cells 
		vector<YeeCell *> bc = boundary_cells[b->first];
		ofile << "Boundary (" << get_string(b->first) << ") has " << bc.size() << " cells" << endl;
		ofile << "with indices: " ;
		ip.print_summary(ofile);

		unsigned int bci = 0;
		unsigned int inind;
		for (auto i=ip.bounds.min.i; i<=ip.bounds.max.i; i++){
			for (auto j=ip.bounds.min.j; j<=ip.bounds.max.j; j++){
				for (auto k=ip.bounds.min.k; k<=ip.bounds.max.k; k++){
					
					ind = mesh.cell_map_to_global_ind(i,j,k);
					indxm = mesh.cell_map_to_global_ind(i-1,j,k);
					indxp = mesh.cell_map_to_global_ind(i+1,j,k);
					indym = mesh.cell_map_to_global_ind(i,j-1,k);
					indyp = mesh.cell_map_to_global_ind(i,j+1,k);
					indzm = mesh.cell_map_to_global_ind(i,j,k-1);
					indzp = mesh.cell_map_to_global_ind(i,j,k+1);

					// ofile << "got mapped indices" << endl;
					// ofile << "ijk: (" << i << "," << j << "," << k << ") ... ind: " << ind << endl;

					if (i > 0){
						inner_cells[ind]->set_neighbor(BoundaryLocation::X_MIN, inner_cells[indxm]);
					}
					if (i < ib.max.i){
						inner_cells[ind]->set_neighbor(BoundaryLocation::X_MAX, inner_cells[indxp]);
					}
					if (j > 0){
						inner_cells[ind]->set_neighbor(BoundaryLocation::Y_MIN, inner_cells[indym]);
					}
					if (j < ib.max.j){
						inner_cells[ind]->set_neighbor(BoundaryLocation::Y_MAX, inner_cells[indyp]);
					}
					if (k > 0){
						inner_cells[ind]->set_neighbor(BoundaryLocation::Z_MIN, inner_cells[indzm]);
					}
					if (k < ib.max.k){
						inner_cells[ind]->set_neighbor(BoundaryLocation::Z_MAX, inner_cells[indzp]);
					}

					inner_cells[ind]->set_neighbor(b->first, bc[bci]);
					bci++;
				}
			}
		}

	}

	// check that linking succeeded
	for (auto c=inner_cells.begin(); c!=inner_cells.end(); c++){
		if ((*c)->num_neighbors() != 6){
			ofile << "WE HAVE A PROBLEM " << endl;
			(*c)->print_summary(ofile);
		}
	}


	ofile << "LINKED BOUNDARY CELLS" << endl;
}

// set how many time steps for the simulation
void OOFDTD::set_nsteps(unsigned int steps){
	nsteps = steps;
}

// set source points (hard source only for now)
void OOFDTD::set_sources(RealPoint rp, Field f, SignalGenerator * sg){
	// does the point exist in this subdomain?
	if (! dom.contains_point(rp, mpi::rank())) return;

	ofile << "index for source: " << mesh.nearest_element(rp) << endl;

	// replace the YeeCell at the point with a SourceCell
	SourceCell * sc = new SourceCell(sg, f);
	ofile << "created source cell" << endl;
	sc->print_summary(ofile);
	ofile << "need to replace neighbors " << inner_cells[mesh.nearest_element(rp)]->get_neighbors().size() << endl;
	map<BoundaryLocation, YeeCell *> neighb = inner_cells[mesh.nearest_element(rp)]->get_neighbors();
	for (auto n = neighb.begin(); n!=neighb.end(); n++){
		ofile << "Boundary Location: " << get_string(n->first) << endl;
		sc->set_neighbor(n->first, n->second);
	}
	for (auto n = sc->get_neighbors().begin(); n!=sc->get_neighbors().end(); n++){
		n->second->set_neighbor(opposite_boundary(n->first), sc);
	}
	sc->print_summary(ofile);
	YeeCell * yc = inner_cells[mesh.nearest_element(rp)];
	inner_cells[mesh.nearest_element(rp)] = sc;
	//delete[] yc;
	sc->print_summary(ofile);
}

void OOFDTD::capture_data(vector<unsigned int> timesteps, vector<Field> f){

}

// running the simulation
void OOFDTD::run(){
	// ofile << "YEP I AM RUNNING THE STUFF" << endl;
	Timer tot;
	Timer tE, tH;

	if (mpi::is_master()) tot.start();
	for (auto i=1; i<nsteps; i++){



		tstep = i;
		if (mpi::is_master()) cout << "on time step " << tstep << "/" << nsteps-1 << "\r" << flush;
		update(Field::D_X); update(Field::E_X); communicate(Field::E_X);
		//cout << "updated DX" << endl;
		update(Field::D_Y); update(Field::E_Y);	communicate(Field::E_Y);
		//cout << "updated DY" << endl;
		update(Field::D_Z); update(Field::E_Z); communicate(Field::E_Z);
		//cout << "updated DZ" << endl;

		update(Field::B_Z); update(Field::H_Z); communicate(Field::H_Z);
		//cout << "updated BZ" << endl;
		update(Field::B_Y); update(Field::H_Y); communicate(Field::H_Y);
		//cout << "update BY" << endl;
		update(Field::B_X); update(Field::H_X); communicate(Field::H_X);
		//cout << "update BX" << endl;

		tstep++;
	}
	if (mpi::is_master()){
		tot.stop();
		cout << "Simulation Time (s): " << tot.totaltime() << endl;
	}

	// do some quick and dirty cleanup
	ofile << "I finished iterating" << endl;
	ofile.close();

}

void OOFDTD::update(Field f){
	if (f == Field::D_X){
		//ofile << "about to update Dx" << endl;
		for (auto c=inner_cells.begin(); c!=inner_cells.end(); c++)	{
			// (*c)->print_summary(ofile);
			// map<BoundaryLocation, YeeCell *> neighb = (*c)->get_neighbors();
			// ofile << endl;
			// for (auto i=neighb.begin(); i!=neighb.end(); i++){
			// 	if (i->second == nullptr){
			// 		ofile << get_string(i->first) << " Neighbor pointer is null!" << endl;
			// 	}
			// }
			(*c)->update_Dx(tstep, delta_t, delta_x);
		}
	}
	else if (f == Field::D_Y){
		//ofile << "about to update Dy" << endl;
		for (auto c=inner_cells.begin(); c!=inner_cells.end(); c++) (*c)->update_Dy(tstep, delta_t, delta_x);
	}
	else if (f == Field::D_Z){
		//ofile << "about to update Dz" << endl;
		for (auto c=inner_cells.begin(); c!=inner_cells.end(); c++) (*c)->update_Dz(tstep, delta_t, delta_x);
	}
	else if (f == Field::B_X){
		for (auto c=inner_cells.begin(); c!=inner_cells.end(); c++) (*c)->update_Bx(tstep, delta_t, delta_x);
	}
	else if (f == Field::B_Y){
		for (auto c=inner_cells.begin(); c!=inner_cells.end(); c++) (*c)->update_By(tstep, delta_t, delta_x);
	}
	else if (f == Field::B_Z){
		for (auto c=inner_cells.begin(); c!=inner_cells.end(); c++) (*c)->update_Bz(tstep, delta_t, delta_x);
	}

	else if (f == Field::E_X){
		for (auto c=inner_cells.begin(); c!=inner_cells.end(); c++)	(*c)->update_Ex(tstep, delta_t, delta_x);
	}
	else if (f == Field::E_Y){
		for (auto c=inner_cells.begin(); c!=inner_cells.end(); c++) (*c)->update_Ey(tstep, delta_t, delta_x);
	}
	else if (f == Field::E_Z){
		for (auto c=inner_cells.begin(); c!=inner_cells.end(); c++) (*c)->update_Ez(tstep, delta_t, delta_x);
	}
	else if (f == Field::H_X){
		for (auto c=inner_cells.begin(); c!=inner_cells.end(); c++) (*c)->update_Hx(tstep, delta_t, delta_x);
	}
	else if (f == Field::H_Y){
		for (auto c=inner_cells.begin(); c!=inner_cells.end(); c++) (*c)->update_Hy(tstep, delta_t, delta_x);
	}
	else if (f == Field::H_Z){
		for (auto c=inner_cells.begin(); c!=inner_cells.end(); c++) (*c)->update_Hz(tstep, delta_t, delta_x);
	}

	
}

void OOFDTD::communicate(Field f){
	// loop over the boundaries
	for (auto b=boundaries.begin(); b!=boundaries.end(); b++){
		BoundaryLocation bl;
		
		if (dom.is_proc_red(mpi::rank())) bl = b->first;
		else bl = opposite_boundary(b->first);

		// if a boundary is a parallel boundary
		if (boundaries[bl] == BoundaryType::PARALLEL){
			// get the plane of indices on this boundary
			IndexBounds ib = dom.get_proc_local_ibounds_cell(mpi::rank());
			IndexPlane ip = ib.get_boundary_plane(bl);

			// get the boundary cells 
			vector<YeeCell *> bc = boundary_cells[bl];
			double * buffer = boundary_buffer_out[bl];

			// prepare buffer from the boundary cells
			unsigned int bci = 0;
			unsigned int ind;
			for (auto i=ip.bounds.min.i; i<=ip.bounds.max.i; i++){
				for (auto j=ip.bounds.min.j; j<=ip.bounds.max.j; j++){
					for (auto k=ip.bounds.min.k; k<=ip.bounds.max.k; k++){

						ind = mesh.cell_map_to_global_ind(i,j,k);
						buffer[bci] = inner_cells[ind]->Hx();
						bci++;
					}
				}
			}

			//ofile << "comm for boundary " << get_string(bl) << " with proc " << dom.get_proc_neighbor(mpi::rank(), bl) << endl;

			// communicate in a red-black manner
			MPI_Request req;
			MPI_Status status;
			mpi::sendrecv(boundary_buffer_out[bl], bc.size(),
					 MPI_DOUBLE, dom.get_proc_neighbor(mpi::rank(), bl), 0, 
					 boundary_buffer_in[bl], bc.size(),
					 MPI_DOUBLE, dom.get_proc_neighbor(mpi::rank(), bl), 0, MPI_COMM_WORLD,
					 &status);


			// grab data from buffer into the ghost cells
			GhostCell * gc;
			for (auto i=0; i<bc.size(); i++){
				gc = (GhostCell *) bc[i];
				gc->fetch_data(f);
			}
		}
	}
}
