#ifndef _FDTDCOMMON_H
#define _FDTDCOMMON_H

#include "RegularMesh.hpp"
#include "SignalGenerator.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <functional>
#include <algorithm>
#include <complex>


	const double m_pi = 3.14159265358979323846264338327950288;
	const double m_eps0 = 8.854e-12;
	const double m_mu0 = m_pi*4.0e-7;
	const double m_c0 = 2.99792458e+8;
	const double m_imp0 = sqrt(m_mu0/m_eps0);

	// some enums and defines for FDTD code
	enum class EMMode : unsigned char {TE=0, TM, ThreeD};
	enum class Field : unsigned char {E_X=0, E_Y, E_Z, D_X, D_Y, D_Z, H_X, H_Y, H_Z, B_X, B_Y, B_Z, NO_FIELD};
	enum class OutputField : unsigned char{E_X, E_Y, E_Z, E_MAGN, H_X, H_Y, H_Z, H_MAGN, S_X, S_Y, S_Z, POYNTING};
	enum class BoundaryType : unsigned char {PARALLEL, PERIODIC, BLOCH_PERIODIC, PEC, PMC};
	enum class SourceType : unsigned char {NONE, HARD, SOFT, TFSF, TFSF_PLANE};
	enum class FluxRegionType : unsigned char {PLANE, RANDOM_POINTS};


	//////////// info associated with Mode
	const std::vector<std::string> mode_strings = {"TE","TM","3D"};
	const std::vector<Field> te_fields = {Field::D_X, Field::D_Y, Field:: B_Z};
	const std::vector<Field> tm_fields = {Field::D_Z, Field::B_X, Field:: B_Y};
	const std::vector<Field> threeD_fields = {Field::D_X, Field::D_Y, Field::D_Z,
											  Field::B_X, Field::B_Y, Field::B_Z};
	const std::vector<std::vector<Field>> mode_fields = {te_fields, tm_fields, threeD_fields};
	const std::vector<Field> te_allfields = {Field::E_X, Field::E_Y, 
											 Field::D_X, Field::D_Y,
											 Field::H_Z, Field:: B_Z};
	const std::vector<Field> tm_allfields = {Field::E_Z, Field::D_Z,
											 Field::H_X, Field::H_Y,
											 Field::B_X, Field:: B_Y};
	const std::vector<Field> threeD_allfields = {Field::E_X, Field::E_Y, Field::E_Z,
												 Field::D_X, Field::D_Y, Field::D_Z,
												 Field::H_X, Field::H_Y, Field::H_Z,
											  	 Field::B_X, Field::B_Y, Field::B_Z};
	const std::vector<std::vector<Field>> mode_allfields = {te_allfields, tm_allfields, threeD_allfields};

	const std::vector<Direction> te_directions = {Direction::X, Direction::Y};
	const std::vector<Direction> tm_directions = {Direction::X, Direction::Y};
	const std::vector<Direction> threeD_directions = {Direction::X, Direction::Y, Direction::Z};
	const std::vector<std::vector<Direction>> mode_directions = {te_directions, tm_directions, threeD_directions};

	const std::vector<OutputField> te_ofields = {OutputField::E_X, OutputField::E_Y, OutputField::E_MAGN, 
												 OutputField:: H_Z, OutputField::H_MAGN,
												 OutputField::S_X, OutputField::S_Y, OutputField::POYNTING};
	const std::vector<OutputField> tm_ofields = {OutputField::E_Z, OutputField::E_MAGN, 
												 OutputField::H_X, OutputField:: H_Y, OutputField::H_MAGN,
													   OutputField::S_X, OutputField::S_Y, OutputField::POYNTING};
	const std::vector<OutputField> threeD_ofields = {OutputField::E_X, OutputField::E_Y, OutputField::E_Z, OutputField::E_MAGN,
											  OutputField::H_X, OutputField::H_Y, OutputField::H_Z, OutputField::H_MAGN,
											  OutputField::S_X, OutputField::S_Y, OutputField::S_Z, OutputField::POYNTING};
	const std::vector<std::vector<OutputField>> mode_ofields = {te_ofields, tm_ofields, threeD_ofields};
	const std::vector<BoundaryLocation> te_bounds = {BoundaryLocation::X_MIN, BoundaryLocation::X_MAX,
													 BoundaryLocation::Y_MIN, BoundaryLocation::Y_MAX};
	const std::vector<BoundaryLocation> tm_bounds = {BoundaryLocation::X_MIN, BoundaryLocation::X_MAX,
													 BoundaryLocation::Y_MIN, BoundaryLocation::Y_MAX};
	const std::vector<BoundaryLocation> threeD_bounds = {BoundaryLocation::X_MIN, BoundaryLocation::X_MAX,
													 BoundaryLocation::Y_MIN, BoundaryLocation::Y_MAX,
													 BoundaryLocation::Z_MIN, BoundaryLocation::Z_MAX};
	const std::vector<std::vector<BoundaryLocation>> mode_boundaries = {te_bounds, tm_bounds, threeD_bounds};

	inline std::vector<Field> get_fields(EMMode mode){
		return mode_fields[(int) mode];
	}

	inline std::vector<Field> get_all_fields(EMMode mode){
		return mode_allfields[(int) mode];
	}

	inline std::vector<OutputField> get_ofields(EMMode mode){
		return mode_ofields[(int) mode];
	}

	inline std::vector<BoundaryLocation> get_boundaries(EMMode mode){
		return mode_boundaries[(int)mode];
	}

	inline std::vector<Direction> get_directions(EMMode mode){
		return mode_directions[(int)mode];
	}

	inline std::string get_string(EMMode mode){
		return mode_strings[(int)mode];
	}
	//////////////////////////////////////////////



	///////////// info associated with Field
	const std::vector<std::string> field_strings = {"E_x", "E_y", "E_z",
													"D_x", "D_y", "D_z",
													"H_x", "H_y", "H_z",
													"B_x", "B_y", "B_z"};
	// first boundary that the field rests on
	const std::vector<BoundaryLocation> field_boundary_1 = {BoundaryLocation::Z_MIN, BoundaryLocation::X_MIN, BoundaryLocation::Y_MIN,
													  BoundaryLocation::Z_MIN, BoundaryLocation::X_MIN, BoundaryLocation::Y_MIN,
													  BoundaryLocation::X_MIN, BoundaryLocation::Y_MIN, BoundaryLocation::Z_MIN,
													  BoundaryLocation::X_MIN, BoundaryLocation::Y_MIN, BoundaryLocation::Z_MIN};

	// second boundary that the field rests on
	// (might be the same as the first boundary)
	const std::vector<BoundaryLocation> field_boundary_2 = {BoundaryLocation::Y_MIN, BoundaryLocation::Z_MIN, BoundaryLocation::X_MIN,
												      BoundaryLocation::Y_MIN, BoundaryLocation::Z_MIN, BoundaryLocation::X_MIN,
												      BoundaryLocation::X_MIN, BoundaryLocation::Y_MIN, BoundaryLocation::Z_MIN,
												      BoundaryLocation::X_MIN, BoundaryLocation::Y_MIN, BoundaryLocation::Z_MIN};

	// subfield for a given field (e.g. E is the subfield of D, H is the subfield of B)
	const std::vector<Field> field_subfield = {Field::E_X, Field::E_Y, Field::E_Z,
												 Field::E_X, Field::E_Y, Field::E_Z,
												 Field::H_X, Field::H_Y, Field::H_Z,
											  	 Field::H_X, Field::H_Y, Field::H_Z};

	// output field (nonnormalized field) for a given field
	const std::vector<OutputField> field_ofield = {OutputField::E_X, OutputField::E_Y, OutputField::E_Z,
												 OutputField::E_X, OutputField::E_Y, OutputField::E_Z,
												 OutputField::H_X, OutputField::H_Y, OutputField::H_Z,
											  	 OutputField::H_X, OutputField::H_Y, OutputField::H_Z};


	// first field in curl term
	const std::vector<Field> mcurl_g1 = {Field::NO_FIELD, Field::NO_FIELD, Field::NO_FIELD,
								  Field::H_Z, Field::H_X, Field::H_Y, 
								  Field::NO_FIELD, Field::NO_FIELD, Field::NO_FIELD, 
								  Field::E_Z, Field::E_X, Field::E_Y};
	// second field in curl term
	const std::vector<Field> mcurl_g2 = {Field::NO_FIELD, Field::NO_FIELD, Field::NO_FIELD, 
								  Field::H_Y, Field::H_Z, Field::H_X, 
								  Field::NO_FIELD, Field::NO_FIELD, Field::NO_FIELD, 
								  Field::E_Y, Field::E_Z, Field::E_X};

	// first derivative direction in curl term
	const std::vector<Direction> mcurl_d1 = {Direction::NO_DIRECTION, Direction::NO_DIRECTION, Direction::NO_DIRECTION,
									  Direction::Y, Direction::Z, Direction::X,
									  Direction::NO_DIRECTION, Direction::NO_DIRECTION, Direction::NO_DIRECTION,
									  Direction::Y, Direction::Z, Direction::X};

	// second derivative direction in curl term
	const std::vector<Direction> mcurl_d2 = {Direction::NO_DIRECTION, Direction::NO_DIRECTION, Direction::NO_DIRECTION,
									  Direction::Z, Direction::X, Direction::Y,
									  Direction::NO_DIRECTION, Direction::NO_DIRECTION, Direction::NO_DIRECTION,
									  Direction::Z, Direction::X, Direction::Y};

	// describes the direction which each field points
	const std::vector<Direction> field_dir = {Direction::X, Direction::Y, Direction::Z,
											  Direction::X, Direction::Y, Direction::Z,
											  Direction::X, Direction::Y, Direction::Z,
											  Direction::X, Direction::Y, Direction::Z};

	
	// coefficient term in front of curl
	const std::vector<double> mcurl_a = {0.0, 0.0, 0.0,
								  1.0/m_mu0, 1.0/m_mu0, 1.0/m_mu0,
								  0.0, 0.0, 0.0,
								  -1.0/m_eps0, -1.0/m_eps0, -1.0/m_eps0};

	const std::vector<double> mnorm_factor = {m_imp0, m_imp0, m_imp0,
											  1.0/m_c0, 1.0/m_c0, 1.0/m_c0,
											  1.0/m_imp0, 1.0/m_imp0, 1.0/m_imp0,
											  1.0/m_c0, 1.0/m_c0, 1.0/m_c0};

	// staggered grid component offsets for derivatives
	const std::vector<int> dx_offset_l = {-1, -1, -1,
								  			-1, -1, -1,
								  			-1, 0, 0,
								  			-1, 0, 0};
	const std::vector<int> dx_offset_r = {1, 0, 0,
								  			1, 0, 0,
								  			1, 1, 1,
								  			1, 1, 1};
	const std::vector<int> dy_offset_l = {-1, -1, -1,
								  			-1, -1, -1,
								  			0, -1, 0,
								  			0, -1, 0};
	const std::vector<int> dy_offset_r = {0, 1, 0,
								  			0, 1, 0,
								  			1, 1, 1,
								  			1, 1, 1};
	const std::vector<int> dz_offset_l = {-1, -1, -1,
								  			-1, -1, -1,
								  			0, 0, -1,
								  			0, 0, -1};
	const std::vector<int> dz_offset_r = {0, 0, 1,
								  			0, 0, 1,
								  			1, 1, 1,
								  			1, 1, 1};
	const std::vector<int> zero_offset = {0, 0, 0,
								  			0, 0, 0,
								  			0, 0, 0,
								  			0, 0, 0};

	const std::vector<double> offset_x = {0, -0.5, -0.5,
							  			0, -0.5, -0.5,
							  			-0.5, 0, 0,
							  			-0.5, 0, 0};

	const std::vector<double> offset_y = {-0.5, 0, -0.5,
							  			-0.5, 0, -0.5,
							  			0, -0.5, 0,
							  			0, -0.5, 0};

	const std::vector<double> offset_z = {-0.5, -0.5, 0,
							  			-0.5, -0.5, 0,
							  			0, 0, -0.5,
							  			0, 0, -0.5};

	const std::vector<std::vector<int>> zoffset = {zero_offset,
													  zero_offset,
													  zero_offset,
													  };
	const std::vector<std::vector<int>> offset_l = {dx_offset_l,
													  dy_offset_l,
													  dz_offset_l,
													  };

	const std::vector<std::vector<int>> offset_r = {dx_offset_r,
													  dy_offset_r,
													  dz_offset_r,
													  };

	inline bool is_magnetic(Field f){
		int fi = (int) f;
		return (fi>5 && fi<=11);
	}
	inline bool is_electric(Field f){
		int fi = (int) f;
		return (fi>=0 && fi<6);
	}

	// }
	inline Field get_E(Direction d){
		return Field(int(Field::E_X) + (int)d);
	}

	inline Field get_H(Direction d){
		return Field(int(Field::H_X) + (int)d);
	}

	inline Field get_subfield(Field f){
		return field_subfield[(int)f];
	}

	inline OutputField get_output_field(Field f){
		return field_ofield[(int)f];
	}

	inline std::string get_string(Field f){
		return field_strings[(int)f];
	}

	inline double normalization_factor(Field f){
		return mnorm_factor[(int)f];
	}

	inline bool field_exists(EMMode mode, Field f){
		std::vector<Field> fields = mode_fields[(int)mode];
		for (auto i=0; i<fields.size(); i++) if (fields[i] == f) return true;
			return false;
	}

	inline Direction field_direction(Field f){
		return field_dir[(int)f];
	}

	inline bool on_boundary(Field f, BoundaryLocation bl){
		return (field_boundary_1[(int)f] == bl || field_boundary_2[(int)f] == bl);
	}

	inline bool normal_to_boundary(Field f, BoundaryLocation bl){
		return (on_boundary(f, bl) && field_dir[(int)f]==boundary_dir[(int)bl]);
	}
	
	inline bool tangent_to_boundary(Field f, BoundaryLocation bl){
		return (on_boundary(f, bl) && field_dir[(int)f]!=boundary_dir[(int)bl]);
	}

	inline RealPoint get_offset(Field f){
		return RealPoint(offset_x[(int)f], offset_y[(int)f], offset_z[(int)f]);
	}
	//////////////////////////////////


	///////////// info associated with OutputField
	const std::vector<std::string> ofield_strings = {"E_x", "E_y", "E_z", "E_magn",
													"H_x", "H_y", "H_z", "H_magn",
													"S_x", "S_y", "S_z",
													"Poynting"};
	
	inline std::string get_string(OutputField of){
		return ofield_strings[(int)of];
	}
	inline bool field_exists(EMMode mode, OutputField of){
		std::vector<OutputField> ofields = mode_ofields[(int)mode];
		for (auto i=0; i<ofields.size(); i++) if (ofields[i] == of) return true;
			return false;
	}
	//////////////////////////////////


	//////////// info associated with BoundaryCondition
	const std::vector<std::string> boundary_condition_strings = {"Parallel", "Periodic", "Bloch Periodic",
																 "PEC", "PMC"};
	inline std::string get_string(BoundaryType bc){
		return boundary_condition_strings[(int)bc];
	}
	/////////////////////////////////


	//////////// info associated with SourceType
	const std::vector<std::string> source_strings = {"None", "Hard", "Soft", "TFSF", "TFSF Plane"};
	inline std::string get_string(SourceType st){
		return source_strings[(int)st];
	}
	/////////////////////////////////


	class FDTDPoint{
	public:

		FDTDPoint(double x, double y=0.0, double z=0.0)
		: m_x(x)
		, m_y(y)
		, m_z(z)
		{
		}

		double x() const {return m_x;};
		double y() const {return m_y;};
		double z() const {return m_z;};

	private: 
		double m_x, m_y, m_z;
	};


	inline Field curl_g1(Field f){return mcurl_g1[(unsigned char)f];};
	inline Field curl_g2(Field f){return mcurl_g2[(unsigned char)f];};
	inline Direction curl_d1(Field f){return mcurl_d1[(unsigned char)f];};
	inline Direction curl_d2(Field f){return mcurl_d2[(unsigned char)f];};
	inline std::vector<Direction> curl_d(Field f){return {curl_d1(f), curl_d2(f)};};
	inline double curl_a(Field f){return mcurl_a[(unsigned char)f];};


	inline IndexPoint deriv_pos_offset(Field f, Direction d){
		std::vector<std::vector<int>> off(zoffset);
		// std::cout << "\n direction : " << get_string(d) << std::endl;
		// std::cout << " field: " << get_string(f) << std::endl;
		// std::cout << "pos: " << off[(int)d][(int)f] << " new: " << offset_r[(int)d][(int)f] << std::endl;
		off[(int)d][(int)f] = offset_r[(int)d][(int)f];
		int fi = (int) f;
		return IndexPoint(off[0][fi], off[1][fi], off[2][fi]);
	}

	inline IndexPoint deriv_neg_offset(Field f, Direction d){
		std::vector<std::vector<int>> off(zoffset);
		off[(int)d][(int)f] = offset_l[(int)d][(int)f];
		int fi = (int) f;
		return IndexPoint(off[0][fi], off[1][fi], off[2][fi]);
	}

	// inline IndexBounds get_updatable_bounds(IndexVolume iv, Field f){
	// 	IndexBounds ib = iv.get_bounds_offset(0);
	// 	Direction fdir = field_direction(f);
	// 	std::vector<Direction> dirs = {Direction::X, Direction::Y, Direction::Z};
	// 	IndexPoint dminus = IndexPoint(0,0,0);
	// 	IndexPoint dplus = IndexPoint(0,0,0);
	// 	// std::cout << "field dir: " << get_string(fdir) << std::endl;
	// 	// std::cout << "npts: " << iv.nx << ", " << iv.ny << ", " << iv.nz << std::endl;
		

	// 	for (auto d = 0; d<dirs.size(); d++){
	// 		if (dirs[d] == fdir) continue;
	// 		dminus = dminus + deriv_neg_offset(f, dirs[d]);
	// 		dplus = dplus + deriv_pos_offset(f, dirs[d]);
	// 	}
	// 	for (auto d = 0; d<dirs.size(); d++){
	// 		//std::cout << "npts in direction: " << get_string(dirs[d]) << " = " << iv.get_npts(dirs[d]) << std::endl;

	// 		if (iv.get_npts(dirs[d]) > 1) continue;
	// 		//std::cout << "direction: " << get_string(dirs[d]) << std::endl;
	// 		dminus[dirs[d]] = 0;
	// 		dplus[dirs[d]] = 0;
	// 	}

	// 	IndexPoint minp = ib.min - dminus;
	// 	IndexPoint maxp = ib.max - dplus;
	// 	return IndexBounds(minp, maxp, iv);

	// }

	// inline bool is_updatable(Field f, BoundaryLocation bl){
	// 	Direction d = boundary_direction(bl);
	// 	if (field_direction(f) == d) return true;
		
	// 	IndexPoint minus = deriv_neg_offset(f, d);
	// 	IndexPoint plus = deriv_pos_offset(f, d);
	// 	if (boundary_side(bl) == MinMax::MIN){
	// 		if (minus[d] < 0) return false;
	// 	}
	// 	else{
	// 		if (plus[d] > 0) return false;
	// 	}
	// 	return true;
	// }

	// inline IndexBounds get_updatable_bounds(IndexBounds ib, Field f, Direction d){

	// 	IndexPoint dminus = deriv_neg_offset(f, d);
	// 	IndexPoint dplus = deriv_pos_offset(f, d);
	// 	if (ib.min[d]==ib.max[d]){
	// 		dminus[d] = 0;
	// 		dplus[d] = 0;
	// 	}
	// 	return IndexBounds(ib.min-dminus, ib.max - dplus, ib.get_volume_context());
	// }



/**************************************************
*	Most of the following functions are operations done
*	over an IndexBounds region. These include but are
*	not limited to updates, setting values, and reduction
*	operations
/**************************************************/
	// // updates of the form f[i] = f[i] + a*g[i];
	// inline void update_value(double a,
	// 				  IndexBounds & ibf, double * f,
	// 				  IndexBounds & ibg, const double * g){

	// 	IndexVolume ctxt_f = ibf.get_volume_context();
	// 	IndexVolume ctxt_g = ibg.get_volume_context();
		
	// 	int find, gind;
	// 	// the k loop index is on the outside, because in 2D there is only 1 k layer
	// 	for (auto fk=ibf.min.k, gk=ibg.min.k; fk<=ibf.max.k, gk<=ibg.max.k; fk++, gk++){
	// 		for (auto fj=ibf.min.j, gj=ibg.min.j; fj<=ibf.max.j, gj<=ibg.max.j; fj++, gj++){
	// 			for (auto fi=ibf.min.i, gi=ibg.min.i; fi<=ibf.max.i, gi<=ibg.max.i; fi++, gi++){
							
	// 				find = ctxt_f.get_serial_index(fi,fj,fk);
	// 				gind = ctxt_g.get_serial_index(gi,gj,gk);
	// 				f[find] += a*g[gind];

	// 			}
	// 		}
	// 	}
	// }


	// // updates of the form f[i] = f[i] + a*b[d]*g[i];
	// // where the coefficient b varies in the direction d
	// inline void update_value(double a,
	// 				  IndexBounds & ibf, double * f,
	// 				  IndexBounds & ibg, const double * g,
	// 				  Direction d, const double * b){
	// 	int fi,fj,fk;
	// 	int gi,gj,gk;

	// 	std::vector<std::reference_wrapper<int>> idcs = {fi,fj,fk};
	// 	int &idxd = idcs[(int)d];

	// 	IndexVolume ctxt_f = ibf.get_volume_context();
	// 	IndexVolume ctxt_g = ibg.get_volume_context();
		
	// 	int find, gind;
	// 	// the k loop index is on the outside, because in 2D there is only 1 k layer
	// 	for (fk=ibf.min.k, gk=ibg.min.k; fk<=ibf.max.k, gk<=ibg.max.k; fk++, gk++){
	// 		for (fj=ibf.min.j, gj=ibg.min.j; fj<=ibf.max.j, gj<=ibg.max.j; fj++, gj++){
	// 			for (fi=ibf.min.i, gi=ibg.min.i; fi<=ibf.max.i, gi<=ibg.max.i; fi++, gi++){
							
	// 				find = ctxt_f.get_serial_index(fi,fj,fk);
	// 				gind = ctxt_g.get_serial_index(gi,gj,gk);
	// 				f[find] += a*b[idxd]*g[gind];

	// 			}
	// 		}
	// 	}
	// }

	// // updates of the form f[i] = f[i] + a*g1[i]*g2[i];
	// inline void update_multiply(double a,
	// 				  IndexBounds & ibf, double * f,
	// 				  IndexBounds & ibg1, const double * g1,
	// 				  IndexBounds & ibg2, const double * g2){

	// 	IndexVolume ctxt_f = ibf.get_volume_context();
	// 	IndexVolume ctxt_g1 = ibg1.get_volume_context();
	// 	IndexVolume ctxt_g2 = ibg2.get_volume_context();
		
	// 	int find, g1ind, g2ind;
	// 	// the k loop index is on the outside, because in 2D there is only 1 k layer
	// 	for (auto fk=ibf.min.k, g1k=ibg1.min.k, g2k=ibg2.min.k; fk<=ibf.max.k, g1k<=ibg1.max.k, g2k<=ibg2.max.k; fk++, g1k++, g2k++){
	// 		for (auto fj=ibf.min.j, g1j=ibg1.min.j, g2j=ibg2.min.j; fj<=ibf.max.j, g1j<=ibg1.max.j, g2j<=ibg2.max.j; fj++, g1j++, g2j++){
	// 			for (auto fi=ibf.min.i, g1i=ibg1.min.i, g2i=ibg2.min.i; fi<=ibf.max.i, g1i<=ibg1.max.i, g2i<=ibg2.max.i; fi++, g1i++, g2i++){
							
	// 				find = ctxt_f.get_serial_index(fi,fj,fk);
	// 				g1ind = ctxt_g1.get_serial_index(g1i,g1j,g1k);
	// 				g2ind = ctxt_g2.get_serial_index(g2i,g2j,g2k);
	// 				f[find] += a*g1[g1ind]*g2[g2ind];

	// 			}
	// 		}
	// 	}
	// }

	// // updates of the form f[i] = f[i] + a*(gr[i]-gl[i]);
	// inline void update_deriv(double a,
	// 				  IndexBounds & ibf, double * f,
	// 				  IndexBounds & ibgl, const double * gl,
	// 				  IndexBounds & ibgr, const double * gr){

	// 	int fi,fj,fk;
	// 	int gli,glj,glk;
	// 	int gri,grj,grk;

	// 	IndexVolume ctxt_f = ibf.get_volume_context();
	// 	IndexVolume ctxt_gl = ibgl.get_volume_context();
	// 	IndexVolume ctxt_gr = ibgr.get_volume_context();
		
	// 	int find, glind, grind;
	// 	// the k loop index is on the outside, because in 2D there is only 1 k layer
	// 	for (fk=ibf.min.k, glk=ibgl.min.k, grk=ibgr.min.k; fk<=ibf.max.k, glk<=ibgl.max.k, grk<=ibgr.max.k; fk++, glk++, grk++){
	// 		for (fj=ibf.min.j, glj=ibgl.min.j, grj=ibgr.min.j; fj<=ibf.max.j, glj<=ibgl.max.j, grj<=ibgr.max.j; fj++, glj++, grj++){
	// 			for (fi=ibf.min.i, gli=ibgl.min.i, gri=ibgr.min.i; fi<=ibf.max.i, gli<=ibgl.max.i, gri<=ibgr.max.i; fi++, gli++, gri++){
			
	// 				find = ctxt_f.get_serial_index(fi,fj,fk);
	// 				glind = ctxt_gl.get_serial_index(gli,glj,glk);
	// 				grind = ctxt_gr.get_serial_index(gri,grj,grk);
					
	// 				f[find] += a*(gr[grind]-gl[glind]);

	// 			}
	// 		}
	// 	}
	// }


	// // updates of the form f[i] = f[i] + a*cw[d]*(gr[i]-gl[i]);
	// // this is included to accomodate pml updates
	// inline void update_deriv(double a,
	// 				  Direction d, const double * cw,
	// 				  IndexBounds & ibf, double * f,
	// 				  IndexBounds & ibgl, const double * gl,
	// 				  IndexBounds & ibgr, const double * gr){

	// 	int fi=0,fj=0,fk=0;
	// 	int gli=0,glj=0,glk=0;
	// 	int gri=0,grj=0,grk=0;
	// 	std::vector<std::reference_wrapper<int>> idcs = {fi,fj,fk};
	// 	int &idx = idcs[(int)d];

	// 	IndexVolume ctxt_f = ibf.get_volume_context();
	// 	IndexVolume ctxt_gl = ibgl.get_volume_context();
	// 	IndexVolume ctxt_gr = ibgr.get_volume_context();

	// 	int ct=0;
	// 	int find, glind, grind;
	// 	// the k loop index is on the outside, because in 2D there is only 1 k layer
	// 	for (fk=ibf.min.k, glk=ibgl.min.k, grk=ibgr.min.k; fk<=ibf.max.k, glk<=ibgl.max.k, grk<=ibgr.max.k; fk++, glk++, grk++){
	// 		for (fj=ibf.min.j, glj=ibgl.min.j, grj=ibgr.min.j; fj<=ibf.max.j, glj<=ibgl.max.j, grj<=ibgr.max.j; fj++, glj++, grj++){
	// 			for (fi=ibf.min.i, gli=ibgl.min.i, gri=ibgr.min.i; fi<=ibf.max.i, gli<=ibgl.max.i, gri<=ibgr.max.i; fi++, gli++, gri++){
							
	// 				find = ctxt_f.get_serial_index(fi,fj,fk);
	// 				glind = ctxt_gl.get_serial_index(gli,glj,glk);
	// 				grind = ctxt_gr.get_serial_index(gri,grj,grk);
					
	// 				f[find] += a*cw[idx]*(gr[grind]-gl[glind]);

	// 			}
	// 		}
	// 	}
	// }

	// // updates of the form f[i] = f[i] + a/pmlk[i]*(gr[i]-gl[i]);
	// // this is included to accomodate inclusion of PML
	// inline void update_deriv(double a,
	// 				  IndexBounds & ibf, double * f,
	// 				  IndexBounds & ibgl, const double * gl,
	// 				  IndexBounds & ibgr, const double * gr,
	// 				  Direction d, const double * pmlk){

	// 	int fi=0,fj=0,fk=0;
	// 	int gli=0,glj=0,glk=0;
	// 	int gri=0,grj=0,grk=0;
	// 	std::vector<std::reference_wrapper<int>> idcs = {fi,fj,fk};
	// 	int &idxpml = idcs[(int)d];

	// 	IndexVolume ctxt_f = ibf.get_volume_context();
	// 	IndexVolume ctxt_gl = ibgl.get_volume_context();
	// 	IndexVolume ctxt_gr = ibgr.get_volume_context();

	// 	int ct=0;
	// 	int find, glind, grind;
	// 	// the k loop index is on the outside, because in 2D there is only 1 k layer
	// 	for (fk=ibf.min.k, glk=ibgl.min.k, grk=ibgr.min.k; fk<=ibf.max.k, glk<=ibgl.max.k, grk<=ibgr.max.k; fk++, glk++, grk++){
	// 		for (fj=ibf.min.j, glj=ibgl.min.j, grj=ibgr.min.j; fj<=ibf.max.j, glj<=ibgl.max.j, grj<=ibgr.max.j; fj++, glj++, grj++){
	// 			for (fi=ibf.min.i, gli=ibgl.min.i, gri=ibgr.min.i; fi<=ibf.max.i, gli<=ibgl.max.i, gri<=ibgr.max.i; fi++, gli++, gri++){
							
	// 				find = ctxt_f.get_serial_index(fi,fj,fk);
	// 				glind = ctxt_gl.get_serial_index(gli,glj,glk);
	// 				grind = ctxt_gr.get_serial_index(gri,grj,grk);
					
	// 				f[find] += a/pmlk[idxpml]*(gr[grind]-gl[glind]);

	// 			}
	// 		}
	// 	}

	// 	// ibf.print_summary();
	// 	// ibgl.print_summary();
	// 	// ibgr.print_summary();
	// 	// std::cout << "woo I looped " << ct << " times" << std::endl;
	// 	// throw -1;
	// }

	// // updates of the form f[i] = a*g[i];
	// inline void set_value(double a,
	// 				  IndexBounds & ibf, double * f,
	// 				  IndexBounds & ibg, const double * g){

	// 	// check that the bounds match


	// 	IndexVolume ctxt_f = ibf.get_volume_context();
	// 	IndexVolume ctxt_g = ibg.get_volume_context();
		
	// 	int find, gind;
	// 	// the k loop index is on the outside, because in 2D there is only 1 k layer
	// 	for (auto fk=ibf.min.k, gk=ibg.min.k; fk<=ibf.max.k, gk<=ibg.max.k; fk++, gk++){
	// 		for (auto fj=ibf.min.j, gj=ibg.min.j; fj<=ibf.max.j, gj<=ibg.max.j; fj++, gj++){
	// 			for (auto fi=ibf.min.i, gi=ibg.min.i; fi<=ibf.max.i, gi<=ibg.max.i; fi++, gi++){
							
	// 				find = ctxt_f.get_serial_index(fi,fj,fk);
	// 				gind = ctxt_g.get_serial_index(gi,gj,gk);
					
	// 				f[find] = a*g[gind];

	// 			}
	// 		}
	// 	}
	// }

	// // updates of the form f[i] = a*b[d]*g[i];
	// // where b only varies in 1 direction "d"
	// inline void set_value(double a,
	// 				  IndexBounds & ibf, double * f,
	// 				  IndexBounds & ibg, const double * g,
	// 				  Direction d, double * b){

	// 	// check that the bounds match
	// 	int fi=0,fj=0,fk=0;
	// 	int gi=0,gj=0,gk=0;
	// 	std::vector<std::reference_wrapper<int>> idcs = {fi,fj,fk};
	// 	int &idx = idcs[(int)d];


	// 	IndexVolume ctxt_f = ibf.get_volume_context();
	// 	IndexVolume ctxt_g = ibg.get_volume_context();
		
	// 	int find, gind;
	// 	// the k loop index is on the outside, because in 2D there is only 1 k layer
	// 	for (fk=ibf.min.k, gk=ibg.min.k; fk<=ibf.max.k, gk<=ibg.max.k; fk++, gk++){
	// 		for (fj=ibf.min.j, gj=ibg.min.j; fj<=ibf.max.j, gj<=ibg.max.j; fj++, gj++){
	// 			for (fi=ibf.min.i, gi=ibg.min.i; fi<=ibf.max.i, gi<=ibg.max.i; fi++, gi++){
							
	// 				find = ctxt_f.get_serial_index(fi,fj,fk);
	// 				gind = ctxt_g.get_serial_index(gi,gj,gk);
					
	// 				f[find] = a*b[idx]*g[gind];

	// 			}
	// 		}
	// 	}
	// }

	// // updates of the form f[i] = a;
	// inline void set_value(double a,
	// 				  const IndexBounds & ibf, double * f){

	// 	IndexVolume ctxt_f = ibf.get_volume_context();

	// 	int find;
	// 	// the k loop index is on the outside, because in 2D there is only 1 k layer
	// 	for (auto fk=ibf.min.k; fk<=ibf.max.k; fk++){
	// 		for (auto fj=ibf.min.j; fj<=ibf.max.j; fj++){
	// 			for (auto fi=ibf.min.i; fi<=ibf.max.i; fi++){
							
	// 				find = ctxt_f.get_serial_index(fi,fj,fk);
					
	// 				f[find] = a;

	// 			}
	// 		}
	// 	}
	// }

	// // updates of the form f[i] = a*sqrt(g[i]);
	// inline void set_sqrt(double a,
	// 				  IndexBounds & ibf, double * f,
	// 				  IndexBounds & ibg, const double * g){

	// 	// check that the bounds match


	// 	IndexVolume ctxt_f = ibf.get_volume_context();
	// 	IndexVolume ctxt_g = ibg.get_volume_context();
		
	// 	int find, gind;
	// 	// the k loop index is on the outside, because in 2D there is only 1 k layer
	// 	for (auto fk=ibf.min.k, gk=ibg.min.k; fk<=ibf.max.k, gk<=ibg.max.k; fk++, gk++){
	// 		for (auto fj=ibf.min.j, gj=ibg.min.j; fj<=ibf.max.j, gj<=ibg.max.j; fj++, gj++){
	// 			for (auto fi=ibf.min.i, gi=ibg.min.i; fi<=ibf.max.i, gi<=ibg.max.i; fi++, gi++){
							
	// 				find = ctxt_f.get_serial_index(fi,fj,fk);
	// 				gind = ctxt_g.get_serial_index(gi,gj,gk);
					
	// 				f[find] = a*sqrt(g[gind]);

	// 			}
	// 		}
	// 	}
	// }

	// // conditional update (only update if conditional is true)
	// // updates of the form f[i] = a;
	// inline void set_value_conditional(double a,
	// 				  IndexBounds & ibf, double * f,
	// 				  IndexBounds & ibc, const bool * c){

	// 	IndexVolume ctxt_f = ibf.get_volume_context();

	// 	int find;
	// 	// the k loop index is on the outside, because in 2D there is only 1 k layer
	// 	for (auto fk=ibf.min.k; fk<=ibf.max.k; fk++){
	// 		for (auto fj=ibf.min.j; fj<=ibf.max.j; fj++){
	// 			for (auto fi=ibf.min.i; fi<=ibf.max.i; fi++){
							
	// 				find = ctxt_f.get_serial_index(fi,fj,fk);
					
	// 				if (c[find]) f[find] = a;

	// 			}
	// 		}
	// 	}
	// }

	// // magnitude of a vector in cartesian coordinates
	// // of the form f = a*sqrt(g1^2 + g2^2 + g3^2)
	// // up to 2 null pointers are allowable
	// inline void set_magnitude(double a,
	// 				  IndexBounds & ibf, double * f,
	// 				  IndexBounds & ibg1, const double * g1,
	// 				  IndexBounds & ibg2, const double * g2,
	// 				  IndexBounds & ibg3, const double * g3){

	// 	// set to zero
	// 	set_value(0.0, ibf, f);

	// 	// add components
	// 	if (g1 != nullptr) update_multiply(1.0, ibf, f, ibg1, g1, ibg1, g1);
	// 	if (g2 != nullptr) update_multiply(1.0, ibf, f, ibg2, g2, ibg2, g2);
	// 	if (g3 != nullptr) update_multiply(1.0, ibf, f, ibg3, g3, ibg3, g3);

	// 	// take sqrt
	// 	set_sqrt(a, ibf, f, ibf, f);

	// }


	// // integration operation
	// // updates of the form a = sum(f[i]);
	// inline void integrate(double & a,
	// 				  IndexBounds & ibf, const double * f){

	// 	IndexVolume ctxt_f = ibf.get_volume_context();

	// 	int find;
	// 	a=0;
	// 	// the k loop index is on the outside, because in 2D there is only 1 k layer
	// 	for (auto fk=ibf.min.k; fk<=ibf.max.k; fk++){
	// 		for (auto fj=ibf.min.j; fj<=ibf.max.j; fj++){
	// 			for (auto fi=ibf.min.i; fi<=ibf.max.i; fi++){
							
	// 				find = ctxt_f.get_serial_index(fi,fj,fk);
					
	// 				a += f[find];

	// 			}
	// 		}
	// 	}
	// }


	// // integration operation on tangential fields
	// // updates of the form a = int(f x n dS);
	// inline void integrate_tangential(double & a,
	// 				  IndexBounds & ibf, const double * f){

	// 	IndexVolume ctxt_f = ibf.get_volume_context();

	// 	int find;
	// 	a=0;
	// 	// the k loop index is on the outside, because in 2D there is only 1 k layer
	// 	for (auto fk=ibf.min.k; fk<=ibf.max.k; fk++){
	// 		for (auto fj=ibf.min.j; fj<=ibf.max.j; fj++){
	// 			for (auto fi=ibf.min.i; fi<=ibf.max.i; fi++){
							
	// 				find = ctxt_f.get_serial_index(fi,fj,fk);
					
	// 				a += f[find];

	// 			}
	// 		}
	// 	}
	// }

	// // integration operation on normal fields
	// // updates of the form a = int(f . n dS);
	// inline void integrate_normal(double & a,
	// 				  IndexBounds & ibf, const double * f){

	// 	IndexVolume ctxt_f = ibf.get_volume_context();

	// 	int find;
	// 	a=0;
	// 	// the k loop index is on the outside, because in 2D there is only 1 k layer
	// 	for (auto fk=ibf.min.k; fk<=ibf.max.k; fk++){
	// 		for (auto fj=ibf.min.j; fj<=ibf.max.j; fj++){
	// 			for (auto fi=ibf.min.i; fi<=ibf.max.i; fi++){
							
	// 				find = ctxt_f.get_serial_index(fi,fj,fk);
					
	// 				a += f[find];

	// 			}
	// 		}
	// 	}
	// }



	// class TFSF_aux_sim{
	// public:
	// 	TFSF_aux_sim(){
	// 		m_D = nullptr;
	// 		m_B = nullptr;
	// 	}

	// 	TFSF_aux_sim(const SignalGenerator & sg, 
	// 				 double dx, unsigned int nnodes, 
	// 				 double dt, unsigned int ntime)
	// 	: m_np(nnodes)
	// 	, m_nt(ntime)
	// 	, m_dt(dt)
	// 	, m_dx(dx)
	// 	{
	// 		m_vphase = 1.0;
	// 		//print_summary();

	// 		// allocate
	// 		m_D = new double*[m_nt];
	// 		m_B = new double*[m_nt];
	// 		for (auto n=0; n<m_nt; n++){
	// 			m_D[n] = new double[m_np];
	// 			m_B[n] = new double[m_np];
	// 		}

	// 		// std::cout << "np: " << m_np << std::endl;

	// 		// perform 1D simulation
	// 		for (auto n=0; n<m_nt; n++){
	// 			for (auto m=0; m<m_np; m++){
	// 				m_D[n][m] = 0.0;
	// 				m_B[n][m] = 0.0;
	// 			}
	// 		}

	// 		// std::cout << "B[0][0]: " << m_B[0][0] << std::endl;

	// 		//vp_ratio = calculate_phase_velocity(0.0)/calculate_phase_velocity(m_TFSF_angle_phi);
	// 		//cout << "phase velocity ratio: " << vp_ratio << endl;
	// 		// cout << "phase velocity 0: " << calculate_phase_velocity(0.0) << endl;

	// 		m_D[0][0] = sg.value(0)/m_imp0;
	// 		for (auto n=1; n<m_nt; n++){
	// 			// time step D normalized
	// 			for (auto m=1; m<m_np-1; m++){
	// 				m_D[n][m] = m_D[n-1][m] - m_dt*m_c0/(m_vphase*m_imp0*m_dx)*(m_B[n-1][m] - m_B[n-1][m-1]);
	// 			}

	// 			// set E source at left end
	// 			m_D[n][0] = sg.value(n*m_dt)/m_imp0;	// divided by imp0 b/c of conversion from E to D

	// 			// time step B normalized
	// 			for (auto m=0; m<m_np-2; m++){
	// 				m_B[n][m] = m_B[n-1][m] - m_dt*m_c0*m_imp0/(m_vphase*m_dx)*(m_D[n][m+1] - m_D[n][m]);
	// 			}
	// 		}

	// 	}

	// 	// assignment operator
	// 	// TFSF_aux_sim operator=(const TFSF_aux_sim & tas){
	// 	// 	TFSF_aux_sim out;
	// 	// 	out.m_vphase = m_vphase;
	// 	// 	out.m_dt = m_dt;
	// 	// 	out.m_dx = m_dx;
	// 	// 	out.m_np = m_np;
	// 	// 	out.m_nt = m_nt;
	// 	// 	out.m_D = new double*[m_nt];
	// 	// 	out.m_B = new double*[m_nt];
	// 	// 	for (auto n=0; n<m_nt; n++){
	// 	// 		out.m_D[n] = new double[m_np];
	// 	// 		out.m_B[n] = new double[m_np];
	// 	// 	}
	// 	// 	for (auto n=0; n<m_nt; n++){
	// 	// 		for (auto m=0; m<m_np; m++){
	// 	// 			out.m_D[n][m] = m_D[n][m];
	// 	// 			out.m_B[n][m] = m_B[n][m];
	// 	// 		}
	// 	// 	}
	// 	// }

	// 	~TFSF_aux_sim(){
	// 		std::cout << "TFSF AUX GOT DESTROYED" << std::endl;
	// 		if (m_D == nullptr && m_B == nullptr) return;
	// 		for (auto n=0; n<m_nt; n++){
	// 			delete[] m_D[n];
	// 			delete[] m_B[n];
	// 		}
	// 		delete[] m_D;
	// 		delete[] m_B;
	// 	}

	// 	double get_value(Field f, double d, unsigned int tstep) const {
	// 		double dprime;
	// 		// std::cout << "dproj: " << d << " total len: " << m_dx*m_np << std::endl;
	// 		// std::cout << "timestep: " << tstep << std::endl;
	// 		if (is_electric(f)){
	// 			double dpp = d/m_dx;
	// 			dprime = dpp-int(dpp);
	// 			return ((1-dprime)*m_D[tstep][2+int(dpp)] + dprime*m_D[tstep][3+int(dpp)]);
	// 		}
	// 		else{
	// 			double dpp = d/m_dx + 0.5;
	// 			//if (m_B[0] == nullptr) std::cout << "Aux B is null!" << std::endl;
	// 			// std::cout << "int(dpp): " << int(dpp) << " dpp: " << dpp << std::endl;
	// 			// std::cout << "first term: " << m_B[0][0] << std::endl;
	// 			// std::cout << "second term: " << m_B[tstep][2+int(dpp)] << std::endl;
				
	// 			dprime = dpp - int(dpp);
	// 			return ((1-dprime)*m_B[tstep-1][1+int(dpp)] + dprime*m_B[tstep-1][2+int(dpp)]);
	// 		}
	// 	}

	// 	void print_summary() const{
	// 		std::cout << "TFSF Auxilliary Simulation: V_phase/c0=" << m_vphase;
	// 		std::cout << " dx: " << m_dx << " npoints: " << m_np << " ntime: " << m_nt << std::endl;
	// 	}

	// private:
	// 	double ** m_D;
	// 	double ** m_B;
	// 	double m_vphase;
	// 	double m_dt, m_dx;

	// 	unsigned int m_np, m_nt;
	// };

#endif
