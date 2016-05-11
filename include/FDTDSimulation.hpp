#ifndef _FDTDSIMULATION_H
#define _FDTDSIMULATION_H

#include "FDTDCommon.hpp"
#include "RegularMesh.hpp"
#include "RegularDomain.hpp"
#include "SignalGenerator.hpp"
#include "SignalProcessing.hpp"
#include "SimulationData.hpp"
#include "SimulationDataHDF.hpp"
#include "DispersiveMaterial.hpp"
#include "BoundaryUpdater.hpp"
#include "SourceUpdater.hpp"
#include "Timer.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <functional>
#include <algorithm>
#include <complex>
#include <memory>
#include <utility>

#include <stdlib.h>
#include <math.h>
#include <time.h>

// namespace FDTD{
class BoundaryUpdater;
class SourceUpdater;


struct PointProbe{
	RealPoint point;
	IndexPoint ip;

	unsigned int ind;
	std::vector<Field> fields;
	std::vector<std::vector<double>> values;

	void print_summary(std::ostream & os=std::cout) const{
		point.print_summary(os);
		ip.print_summary(os);
	}
};

struct PointDFT{
	RealPoint point;
	IndexPoint ip;

	unsigned int ind;
	std::vector<double> freqs;
	std::vector<Field> fields;
	std::vector<std::vector<StreamProcessing::dft>> values;

	void print_summary(std::ostream & os=std::cout) const{
		os << " (" << freqs.size() << " freqs)" ;
		point.print_summary(os);
		ip.print_summary(os);
	}
};

struct PortDFT{
	RealPlane rp;
	IndexPlane ip;

	std::vector<double> freqs;
	std::vector<complex<double>> transforms;

	void print_summary(std::ostream & os=std::cout) const{
		os << " (" << freqs.size() << " freqs)" ;
		// rp.print_summary(os);
		// ip.print_summary(os); // this causes a crash
	}
};

struct FluxPlane{
	IndexPlane ip;
	std::vector<double> freqs;
	std::vector<unsigned int> inds;
	std::vector<OutputField> fields;
	std::vector<std::vector<std::vector<StreamProcessing::dft>>> transforms;

	void print_summary(std::ostream & os=std::cout) const{
		//ip.print_summary(os); // this causes a crash in parallel??
		os << " (" << inds.size() << " points)" ;
		os << " (" << freqs.size() << " freqs)" ;
		os << " (Fields:" ;
		for (auto f=fields.begin(); f!= fields.end(); f++){
			os << ", " << get_string(*f) ;
		}
		os << ")" ;
	}
};

struct FluxBox{
	RealPoint center;
	double side;
	std::vector<FluxPlane> planes;
};


class FDTDSimulation{
	friend class BoundaryUpdater;
	friend class BoundaryParallel;
	friend class BoundaryPeriodic;
	friend class SourceUpdaterTFSFPlane;
public:

	FDTDSimulation();
	FDTDSimulation(EMMode mode, RegularDomain & rdom, const RegularMesh & mesh);
	~FDTDSimulation();

	// inspectors
	const double * field_ptr(Field f) const {return m_field_ptrs[(int)f];};
	const SimulationData & simdata() const {return *m_simdata;};
	double dt() const {return m_dt;};
	virtual void print_summary() const;

	// mutators
	void set_mode(EMMode mode);
	void bind_mesh(const RegularMesh & mesh);
	void set_boundary(BoundaryLocation loc, BoundaryCondition type);
	void set_CFL(double cfl);
	void set_PML(Direction d, bool onoff);
	void set_PML_thickness(unsigned int thick);
	void set_SF_thickness(unsigned int thick);
	void add_source(SourceType type, double central_freq, double power, double param1=0.0, double param2=0.0, double param3=0.0);
	void add_random_sources(SourceType type, double central_freq, unsigned int num);
	void set_waveform(const SignalGenerator & sig_gen);
	//void add_flux_region(Direction norm_dir, double width, double x0, double y0, double z0=0.0);
	void add_flux_plane(Direction norm, MinMax norm_side, RealPoint center, double width);
	//void add_flux_box(RealPoint center, double sidelen);
	void add_port_dft(Direction norm, MinMax norm_side, RealPoint center, double width);
	void add_point_probe(RealPoint rp);
	void add_point_dft(RealPoint rp);
	//void add_random_sample(unsigned int npts);
	virtual void add_capture_field(OutputField fieldname);
	void add_material(unsigned int mtl_key, std::shared_ptr<DispersiveMaterial<double>> electric_mtl, std::shared_ptr<DispersiveMaterial<double>> magnetic_mtl);
	void bind_material_keys(const unsigned int * mkeys);
	void bind_material_keys(std::function<unsigned int(unsigned int)> mkey_fn);
	void bind_current_density_x(const double * current_density_x);
	void bind_current_density_y(const double * current_density_y);
	void bind_current_density_z(const double * current_density_z);
	void bind_current_density_x(std::function<double(unsigned int)> current_density_x_fn);
	void bind_current_density_y(std::function<double(unsigned int)> current_density_y_fn);
	void bind_current_density_z(std::function<double(unsigned int)> current_density_z_fn);

	// new generation of mutators
	//void add_source(const SourceInfo & si);

	void set_num_iters(unsigned int num_iters) {m_num_iters = num_iters;};

	void run(int num_iters = -1);	// runs through the full simulation by default

protected:

	void calculate_dt();
	virtual void prerun_check();
	virtual void allocate_fields();
	virtual void allocate_coeffs();
	virtual void allocate_boundaries();
	virtual void allocate_materials();
	virtual void allocate_PML();
	virtual void allocate_simdata();
	virtual void allocate_probes();
	virtual void allocate_ports();
	virtual void allocate_sources();
	virtual void allocate_TFSF();
	virtual double calculate_phase_velocity(double angle);
	virtual double calculate_TFSF_aux_dx();
	void prepare_transforms();

	//virtual void run_internal(int num_iters = -1);
	virtual void run_internal(EMMode mode, int num_iters=-1);

	virtual void update_current_density_E();
	virtual void update_sources_E();
	virtual void update_sources_H();

	virtual void update_ports();
	virtual void update_probes();
	virtual void update_transforms();
	virtual void capture_simdata();

	// new generation of update equations
	virtual void update_curl(Field f, IndexBounds ib);
	virtual void update_sources(Field f);
	virtual void update_boundaries(Field f);
	virtual void update_PML(Field f);
	virtual void update_materials(Field f);
	virtual void update_real(Field f);
	virtual void update_derived();
	virtual void update_comm(Field f);

	virtual void print_local(std::ostream & os = std::cout) const;
	

	// user defined data
	const RegularMesh * m_mesh=nullptr;
	const double * m_current_density_x=nullptr;
	const double * m_current_density_y=nullptr;
	const double * m_current_density_z=nullptr;

	// Index bounds for overall mesh
	//IndexBounds m_ib_glob;

	// data output
	std::ofstream m_ofile;
	//std::ostream & m_ofile = std::cout;
	SimulationData * m_simdata;
	std::vector<OutputField> m_capture_fields;
	std::vector<RealPlane> m_flux_planes;
	//std::vector<FluxPlane> m_flux_boxes;
	std::vector<RealPlane> m_port_dfts;
	std::vector<RealPoint> m_point_probes;
	std::vector<RealPoint> m_point_dfts;


	// internal simulation variables
	double * m_alwaysnull=nullptr;

		//************** FDTD SIMULATION GLOBAL INFO ***************//
		// spatial and grid info
		RegularDomain * m_domain;
		IndexVolume m_local_volume = IndexVolume(0,0,0);		// the volume of this local grid 
		IndexBounds m_global_bounds = m_local_volume.get_bounds_offset(0);	// the bounds of this local volume in the overall grid
		double m_dx=0.0, m_dy=0.0, m_dz=0.0, m_dt=0.0;
		std::vector<std::reference_wrapper<double>> m_delta_x = {m_dx, m_dy, m_dz};
		
		// tfsf layers
		unsigned int m_nTFSF=14;

		// time stepping
		double m_CourantFactor=0.99;
		unsigned int m_num_iters=0, m_current_iter=0;
		double m_tcur=0.0;

		// simulation mode
		EMMode m_simmode;	// mode (TE, TM, or 3D)

		// source type
		SourceType m_sourcetype=SourceType::NONE; // source (TFSF, TFSF_Plane, Hard, Soft)
		//SourceWaveform m_waveform; // source time waveform
		double m_source_freq;		// central frequency in Hz
		double m_source_ampl;		// source amplitude

		// boundaries
		BoundaryCondition m_boundary_x_min=BoundaryCondition::PEC, m_boundary_x_max=BoundaryCondition::PEC;
		BoundaryCondition m_boundary_y_min=BoundaryCondition::PEC, m_boundary_y_max=BoundaryCondition::PEC;
		BoundaryCondition m_boundary_z_min=BoundaryCondition::PEC, m_boundary_z_max=BoundaryCondition::PEC;
		std::vector<std::reference_wrapper<BoundaryCondition>> m_boundaries = {m_boundary_x_min, m_boundary_x_max,
																			   m_boundary_y_min, m_boundary_y_max,
																			   m_boundary_z_min, m_boundary_z_max};

		// pml
		unsigned int m_nPML=10;
		bool m_pml_x=false, m_pml_y=false, m_pml_z=false;
		
		// sources
		std::vector<RealPoint> m_signal_loc;		// contains locations of soft/hard sources
	
		// waveform-specific data
		const SignalGenerator * m_source_signal;		// source signal

		// source type-specific data
		double m_TFSF_angle_phi;			// angle of the incident plane wave off the x axis
		double m_TFSF_angle_theta; 			// angle of the incident plane wave off the z axis (3D only);
		double m_TFSF_polarization_psi; 	// angle that the E field makes with k cross z
		TFSF_aux_sim * m_TFSF_aux;

		// complex materials
		std::map<unsigned int, std::shared_ptr<DispersiveMaterial<double>>> m_electric_mtl_list;	// summary list
		std::map<unsigned int, std::shared_ptr<DispersiveMaterial<double>>> m_magnetic_mtl_list;

		//*********************************************************//



		//************** FDTD SIMULATION LOCAL INFO ***************//
		// In other words, this is the information that exists in this
		// local simulation only (for parallel runs)

		// output flux planes
		std::vector<FluxPlane> m_flux_local;
		std::vector<PortDFT> m_port_local;
		std::map<unsigned int, PointProbe> m_probe_local;
		std::map<unsigned int, PointDFT> m_dft_local;
		std::vector<complex<double>> m_dft_kernel;

		// switches to be set
		bool m_is_allocated=false;
		bool m_fields_allocated=false;

		// boundaries
		std::map<BoundaryLocation, BoundaryUpdater *> m_boundary_updater;

		// pml
		std::vector<std::reference_wrapper<bool>> m_pml_on = {m_pml_x, m_pml_y, m_pml_z};
		std::vector<bool> m_pml_local = {false, false,
										 false, false,
										 false, false};

		// functions used for current density
		std::function<double(unsigned int)> m_current_density_x_fn=nullptr;
		std::function<double(unsigned int)> m_current_density_y_fn=nullptr;
		std::function<double(unsigned int)> m_current_density_z_fn=nullptr;
		
		// complex materials
		std::vector<unsigned int> m_mtl_keys;		// keys that associate the material in the material list to the material at each cell
		std::vector<std::shared_ptr<DispersiveMaterial<double>>> m_electric_mtls_x;	// used for actual computations
		std::vector<std::shared_ptr<DispersiveMaterial<double>>> m_magnetic_mtls_x;	
		std::vector<std::shared_ptr<DispersiveMaterial<double>>> m_electric_mtls_y;	// used for actual computations
		std::vector<std::shared_ptr<DispersiveMaterial<double>>> m_magnetic_mtls_y;	
		std::vector<std::shared_ptr<DispersiveMaterial<double>>> m_electric_mtls_z;	// used for actual computations
		std::vector<std::shared_ptr<DispersiveMaterial<double>>> m_magnetic_mtls_z;	
		std::vector<std::vector<std::shared_ptr<DispersiveMaterial<double>>>*> m_materials = {&m_electric_mtls_x, &m_electric_mtls_y, &m_electric_mtls_z,
																							 &m_electric_mtls_x, &m_electric_mtls_y, &m_electric_mtls_z,
																							 &m_magnetic_mtls_x, &m_magnetic_mtls_y, &m_magnetic_mtls_z,
																							 &m_magnetic_mtls_x, &m_magnetic_mtls_y, &m_magnetic_mtls_z};
		// TODO: delete these... they are deprecated
		std::vector<bool> m_pec_nodes;	// metal nodes receive special treatment
		std::vector<bool> m_pmc_nodes;

		// sources
		std::vector<SourceUpdater *> m_source_updater;


		// CPML parameters (electric)
		double * m_CPML_k_LR=nullptr; // left to right K
		double * m_CPML_k_TB=nullptr; // top to bottom K
		double * m_CPML_k_CS=nullptr; // charm to strange K
		double * m_CPML_sig_LR=nullptr; // left to right sigma
		double * m_CPML_sig_TB=nullptr; // left to right sigma
		double * m_CPML_sig_CS=nullptr; // left to right sigma
		double * m_CPML_a_LR=nullptr; // left to right a
		double * m_CPML_a_TB=nullptr; // left to right a
		double * m_CPML_a_CS=nullptr; // left to right a
		// CPML parameters (magnetic)
		double * m_CPML_k_m_LR=nullptr; // left to right K
		double * m_CPML_k_m_TB=nullptr; // top to bottom K
		double * m_CPML_k_m_CS=nullptr; // charm to strange K
		double * m_CPML_sig_m_LR=nullptr; // left to right sigma
		double * m_CPML_sig_m_TB=nullptr; // left to right sigma
		double * m_CPML_sig_m_CS=nullptr; // left to right sigma
		double * m_CPML_a_m_LR=nullptr; // left to right a
		double * m_CPML_a_m_TB=nullptr; // left to right a
		double * m_CPML_a_m_CS=nullptr; // left to right a
		// CPML integrated quantities
		double * m_CPML_I_Ex_y=nullptr;
		double * m_CPML_I_Ex_z=nullptr;
		double * m_CPML_I_Ey_x=nullptr;
		double * m_CPML_I_Ey_z=nullptr;
		double * m_CPML_I_Ez_y=nullptr;
		double * m_CPML_I_Ez_x=nullptr;
		double * m_CPML_I_Hx_y=nullptr;
		double * m_CPML_I_Hx_z=nullptr;
		double * m_CPML_I_Hy_x=nullptr;
		double * m_CPML_I_Hy_z=nullptr;
		double * m_CPML_I_Hz_y=nullptr;
		double * m_CPML_I_Hz_x=nullptr;
		// CPML derived quantities
		double * m_CPML_b_E_LR=nullptr; 
		double * m_CPML_b_E_TB=nullptr; 
		double * m_CPML_b_E_CS=nullptr; 
		double * m_CPML_c_E_LR=nullptr;
		double * m_CPML_c_E_TB=nullptr;
		double * m_CPML_c_E_CS=nullptr;
		double * m_CPML_b_H_LR=nullptr; 
		double * m_CPML_b_H_TB=nullptr; 
		double * m_CPML_b_H_CS=nullptr;
		double * m_CPML_c_H_LR=nullptr;
		double * m_CPML_c_H_TB=nullptr;
		double * m_CPML_c_H_CS=nullptr;

		

		
		std::vector<std::reference_wrapper<double *>> m_CPML_k = {m_CPML_k_LR, m_CPML_k_TB, m_CPML_k_CS,
										  m_CPML_k_m_LR, m_CPML_k_m_TB, m_CPML_k_m_CS};
		std::vector<std::reference_wrapper<double *>> m_CPML_sig = {m_CPML_sig_LR, m_CPML_sig_TB, m_CPML_sig_CS,
										  m_CPML_sig_m_LR, m_CPML_sig_m_TB, m_CPML_sig_m_CS};
		std::vector<std::reference_wrapper<double *>> m_CPML_a = {m_CPML_a_LR, m_CPML_a_TB, m_CPML_a_CS,
										  m_CPML_a_m_LR, m_CPML_a_m_TB, m_CPML_a_m_CS};
		std::vector<std::reference_wrapper<double *>> m_CPML_b = {m_CPML_b_E_LR, m_CPML_b_E_TB, m_CPML_b_E_CS,
										  m_CPML_b_H_LR, m_CPML_b_H_TB, m_CPML_b_H_CS};
		std::vector<std::reference_wrapper<double *>> m_CPML_c = {m_CPML_c_E_LR, m_CPML_c_E_TB, m_CPML_c_E_CS,
										  m_CPML_c_H_LR, m_CPML_c_H_TB, m_CPML_c_H_CS};
		std::vector<std::reference_wrapper<double *>> m_CPML_I_g1 = {m_alwaysnull, m_alwaysnull, m_alwaysnull,
																	 m_CPML_I_Hz_y, m_CPML_I_Hx_z, m_CPML_I_Hy_x,
																	 m_alwaysnull, m_alwaysnull, m_alwaysnull, 
																	 m_CPML_I_Ez_y, m_CPML_I_Ex_z, m_CPML_I_Ey_x};
		std::vector<std::reference_wrapper<double *>> m_CPML_I_g2 = {m_alwaysnull, m_alwaysnull, m_alwaysnull,
																	 m_CPML_I_Hy_z, m_CPML_I_Hz_x, m_CPML_I_Hx_y,
																	 m_alwaysnull, m_alwaysnull, m_alwaysnull, 
																	 m_CPML_I_Ey_z, m_CPML_I_Ez_x, m_CPML_I_Ex_y};

		double * cpml_g1(Field f){
			return m_CPML_I_g1[(int)f].get();
		}
		double * cpml_g2(Field f){
			return m_CPML_I_g2[(int)f].get();
		}
		double * cpml_k(Field f, Direction d){
			int off = (int)d;
			if (is_magnetic(f)) off+=3;
			return m_CPML_k[off].get();
		}
		double * cpml_sig(Field f, Direction d){
			int off = (int)d;
			if (is_magnetic(f)) off+=3;
			return m_CPML_sig[off].get();
		}
		double * cpml_a(Field f, Direction d){
			int off = (int)d;
			if (is_magnetic(f)) off+=3;
			return m_CPML_a[off].get();
		}
		double * cpml_b(Field f, Direction d){
			int off = (int)d;
			if (is_magnetic(f)) off+=3;
			return m_CPML_b[off].get();
		}
		double * cpml_c(Field f, Direction d){
			int off = (int)d;
			if (is_magnetic(f)) off+=3;
			return m_CPML_c[off].get();
		}

		
		// normalized fields used in the simulation
		// Electric polarization field
		double * m_Dn_x=nullptr;	// note that this is the normalized D field
		double * m_Dn_y=nullptr;
		double * m_Dn_z=nullptr;	

		// Electric fields
		double * m_En_x=nullptr;	// note that this is the normalized E field
		double * m_En_y=nullptr;
		double * m_En_z=nullptr;	
		
		// Magnetic flux density field
		double * m_Bn_x=nullptr;	// normalized
		double * m_Bn_y=nullptr;
		double * m_Bn_z=nullptr;

		// Magnetic fields
		double * m_Hn_x=nullptr;	// note that this is the normalized H field
		double * m_Hn_y=nullptr;
		double * m_Hn_z=nullptr;

		// derived/output fields
		// electric
		double * m_E_x=nullptr;	// real E fields
		double * m_E_y=nullptr;
		double * m_E_z=nullptr;
		double * m_E_magn=nullptr;
		// magnetic
		double * m_H_x=nullptr;	// real H fields
		double * m_H_y=nullptr;
		double * m_H_z=nullptr;
		double * m_H_magn=nullptr;
		// poynting
		double * m_S_x=nullptr;
		double * m_S_y=nullptr;
		double * m_S_z=nullptr;
		double * m_poynting=nullptr;

		std::vector<std::reference_wrapper<double *>> m_field_ptrs = {m_En_x, m_En_y, m_En_z,
									  m_Dn_x, m_Dn_y, m_Dn_z,
									  m_Hn_x, m_Hn_y, m_Hn_z,
									  m_Bn_x, m_Bn_y, m_Bn_z};
		std::vector<std::reference_wrapper<double *>> m_output_field_ptrs = {m_E_x, m_E_y, m_E_z, m_E_magn,
												   	  m_H_x, m_H_y, m_H_z, m_H_magn,
												   	  m_S_x, m_S_y, m_S_z,
												      m_poynting};
		//*********************************************************//


};

#endif