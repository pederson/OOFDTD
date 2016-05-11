#ifndef _REGULARMESH_H
#define _REGULARMESH_H

#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <fstream>
#include <functional>

#include "Mesh.hpp"
#include "mpitools.hpp"

enum class Direction : unsigned char {X=0, Y, Z, NO_DIRECTION};
enum class BoundaryLocation : unsigned char {X_MIN=0, X_MAX, Y_MIN, Y_MAX, Z_MIN, Z_MAX};
enum class MinMax : signed char {MIN=-1, MAX=1};
const std::vector<std::string> direction_strings = {"x","y","z","none"};
const std::vector<std::string> boundary_loc_strings = {"xmin","xmax",
                               "ymin","ymax",
                               "zmin","zmax"};
const std::vector<std::string> minmax_strings = {"min","max"};

// specifies the boundary opposite to each boundary
const std::vector<BoundaryLocation> boundary_opp = {BoundaryLocation::X_MAX, BoundaryLocation::X_MIN, 
                        BoundaryLocation::Y_MAX, BoundaryLocation::Y_MIN,
                        BoundaryLocation::Z_MAX, BoundaryLocation::Z_MIN};

const std::vector<Direction> boundary_dir = {Direction::X, Direction::X, 
                      Direction::Y, Direction::Y,
                      Direction::Z, Direction::Z};

const std::vector<MinMax> boundary_minmax = {MinMax::MIN, MinMax::MAX,
                                             MinMax::MIN, MinMax::MAX,
                                             MinMax::MIN, MinMax::MAX};

const std::vector<int> boundary_offset_x = {-1, 1, 0, 0, 0, 0};
const std::vector<int> boundary_offset_y = {0, 0, -1, 1, 0, 0};
const std::vector<int> boundary_offset_z = {0, 0, 0, 0, -1, 1};


inline BoundaryLocation opposite_boundary(BoundaryLocation bl){
  return boundary_opp[(int)bl];
}

inline std::string get_string(Direction d){
  return direction_strings[(int)d];
}

inline std::string get_string(BoundaryLocation bl){
  return boundary_loc_strings[(int)bl];
}

inline std::string get_string(MinMax side){
  return minmax_strings[(side==MinMax::MIN? 0 : 1)];
}

inline MinMax boundary_side(BoundaryLocation bl){
  return boundary_minmax[(int)bl];
}

inline Direction boundary_direction(BoundaryLocation bl){
  return boundary_dir[(int)bl];
}


class IndexBounds;
class RealPlane;
class IndexPlane;

class RealPoint{
public:
  RealPoint(){};
  RealPoint(double xp, double yp, double zp)
  : x(xp)
  , y(yp)
  , z(zp)
  {
  }

  RealPoint operator+(RealPoint rp){
    return RealPoint(rp.x+x, rp.y+y, rp.z+z);
  }

  double & operator[](Direction d){
    std::vector<std::reference_wrapper<double>> pt= {std::reference_wrapper<double>(x),
						     std::reference_wrapper<double>(y),
						     std::reference_wrapper<double>(z)};
    return pt[(int)d].get();
  }

  const double & operator[](Direction d) const{
    std::vector<std::reference_wrapper<const double>> pt {std::reference_wrapper<const double>(x),
							  std::reference_wrapper<const double>(y),
							  std::reference_wrapper<const double>(z)};
    return pt[(int)d].get();
  }

  double dot(RealPoint rp) const{
    return rp.x*x + rp.y*y + rp.z*z;
  }

  void print_summary(std::ostream & os=std::cout) const {
    os << "(x,y,z): " << "(" << x << ", " << y << ", " << z << ")" ;// << std::endl;
  }

  double x,y,z;
};

class IndexPoint{
public:
  IndexPoint(){};
  IndexPoint(int ii, int ji, int ki)
  : i(ii)
  , j(ji)
  , k(ki)
  {
  }

  IndexPoint operator+(IndexPoint ip){
    return IndexPoint(ip.i+i, ip.j+j, ip.k+k);
  }

  IndexPoint operator-(IndexPoint ip){
    return IndexPoint(i-ip.i, j-ip.j, k-ip.k);
  }

  int & operator[](Direction d){
    std::vector<std::reference_wrapper<int>> pt = {std::reference_wrapper<int>(i),
						   std::reference_wrapper<int>(j),
						   std::reference_wrapper<int>(k)};
    return pt[(int)d].get();
  }

  const int & operator[](Direction d) const{
    std::vector<std::reference_wrapper<const int>> pt = {std::reference_wrapper<const int>(i),
							 std::reference_wrapper<const int>(j),
							 std::reference_wrapper<const int>(k)};
    return pt[(int)d].get();
  }

  void print_summary(std::ostream & os=std::cout) const {
    os << "(i,j,k): " << "(" << i << ", " << j << ", " << k << ")" ;// << std::endl;
  }

  int i,j,k;

};

inline IndexPoint boundary_offset(BoundaryLocation bl){
  int ibl = (int)bl;
  return IndexPoint(boundary_offset_x[ibl], boundary_offset_y[ibl], boundary_offset_z[ibl]);
}

class RealVolume{
public:
  RealVolume(double xvol, double yvol, double zvol)
  : xlen(xvol)
  , ylen(yvol)
  , zlen(zvol)
  , vtot(xlen*ylen*zlen)
  {
  }

  double xlen, ylen, zlen;
  double vtot;
};

class IndexVolume{
public:

  IndexVolume(int numx, int numy, int numz)
  : nx(numx)
  , ny(numy)
  , nz(numz)
  , ntot(numx*numy*numz)
  {
  }

  int & get_npts(Direction d){
    std::vector<std::reference_wrapper<int>> npt = {std::reference_wrapper<int>(nx), 
						    std::reference_wrapper<int>(ny), 
						    std::reference_wrapper<int>(nz)};
    int rf = npt[(int)d].get();
    return npt[(int)d].get();
  }

  int get_serial_index(IndexPoint ip) const{
    return get_serial_index(ip.i, ip.j, ip.k);
  }

  int get_serial_index(int i, int j, int k) const{
    return k*(nx*ny) + j*(nx) + i;
  }

  template<typename T>
  T * generate_volume_array(T init_val) const{
    T * v = new T[ntot];
    for (auto i=0; i<ntot; i++) v[i]=init_val;
    return v;
  }

  IndexBounds get_bounds_offset(int off) const;

  IndexBounds get_bounds_boundary(BoundaryLocation bl) const;

  IndexPlane get_boundary_plane(BoundaryLocation bl) const;

  void print_summary(std::ostream & os=std::cout) const{
    os << "(nx, ny, nz): " << "(" << nx << ", " << ny << ", " << nz << ")" ;
  }

  // IndexBounds get_bounds_boundary(BoundaryLocation bl, int off);

  int nx=0, ny=0, nz=0;
  int ntot=0;

};

class RealBounds{
public:
  RealBounds(RealPoint pmin, RealPoint pmax)
  : min(pmin)
  , max(pmax)
  {
  }

  RealBounds()
  : min(RealPoint(0,0,0))
  , max(RealPoint(0,0,0))
  {
  }

  bool contains_point(RealPoint rp) const{
    std::vector<Direction> alldirs = {Direction::X, Direction::Y, Direction::Z};
    Direction d;
    for (auto i=0; i<alldirs.size(); i++){
      d = alldirs[i];
      if (rp[d] < min[d]) return false;
      if (rp[d] > max[d]) return false;
    }
    return true;
  }

  RealPoint get_side(MinMax side) const {
    if (side==MinMax::MIN) return min;
    else return max;
  }

  RealBounds get_boundary(BoundaryLocation bl){
    Direction d = boundary_direction(bl);
    MinMax side = boundary_side(bl);
    RealBounds out = *this;
    out.min[d] = get_side(side)[d];
    out.max[d] = get_side(side)[d];
    return out;
  }

  RealPlane get_boundary_plane(BoundaryLocation bl) const;

  void print_summary(std::ostream & os=std::cout) const {
    os << "\tmin: " ;
    min.print_summary(os);
    os << " " << std::endl;
    os << "\tmax: " ;
    max.print_summary(os);
    os << " " << std::endl;
  }

  RealPoint min, max;

};

class IndexBounds{
public:
  
  IndexBounds()
  : min(IndexPoint(0,0,0))
  , max(IndexPoint(0,0,0))
  {
  }

  IndexBounds(IndexPoint imin, IndexPoint imax)
  : min(imin)
  , max(imax)
  {
    context = IndexVolume(imax.i-imin.i+1,
                          imax.j-imin.j+1,
                          imax.k-imin.k+1);
  }

  IndexBounds(IndexPoint imin, IndexPoint imax, IndexVolume ctxt)
  : min(imin)
  , max(imax)
  , context(ctxt)
  {
  }

  IndexPoint get_side(MinMax side) const {
    if (side==MinMax::MIN) return min;
    else return max;
  }

  IndexBounds get_boundary(BoundaryLocation bl, int offset=0) const{
    Direction d = boundary_direction(bl);
    MinMax mm = boundary_side(bl);
    IndexPoint minp = min;
    IndexPoint maxp = max;
    if (mm==MinMax::MIN) maxp[d] = minp[d] + offset;
    else minp[d] = maxp[d] - offset;
    return IndexBounds(minp, maxp, context);
  }

  IndexPlane get_boundary_plane(BoundaryLocation bl) const;

  IndexBounds offset(IndexPoint ip) const{
    return IndexBounds(IndexPoint(min.i+ip.i, min.j+ip.j, min.k+ip.k),
                       IndexPoint(max.i+ip.i, max.j+ip.j, max.k+ip.k),
                       context);
  }

  IndexVolume generate_new_volume() const{
    return IndexVolume(max.i-min.i+1,
                       max.j-min.j+1,
                       max.k-min.k+1);
  }

  IndexVolume get_volume_context() const {return context;};

  void print_summary(std::ostream & os=std::cout) const {
    os << "\tmin: " ;
    min.print_summary(os);
    os << " " << std::endl;
    os << "\tmax: " ;
    max.print_summary(os);
    os << " " << std::endl;
  }

  IndexPoint min, max;
private:
  IndexVolume context = IndexVolume(0,0,0);

};

class RealPlane{
public:
  RealPlane(){};
  RealPlane(RealBounds rb, double pt, Direction n, MinMax sd)
  : bounds(rb)
  , normal(n)
  , side(sd)
  {
    bounds.min[n] = pt;
    bounds.max[n] = pt;
  }

  // slides the plane in the normal direction
  // by the amount 'val'
  RealPlane offset(double val) const {
    return RealPlane(bounds, bounds.min[normal] + double(side)*val, normal, side);
  }

  // truncate the plane on the given boundary
  // by the given amount of index points
  RealPlane truncate(BoundaryLocation bl, double off){
    RealBounds bout = bounds;
    if (boundary_side(bl)== MinMax::MIN){
      bout.min[boundary_direction(bl)] += off;
    } 
    else{
      bout.max[boundary_direction(bl)] -= off;
    }
    return RealPlane(bout, bout.min[normal], normal, side);
  }

  bool intersects(RealBounds rb) const {
    if (rb.min[normal] > bounds.min[normal]) return false;
    if (rb.max[normal] < bounds.min[normal]) return false;
    return true;
  }

  void print_summary(std::ostream & os=std::cout) const{
    os << get_string(normal) << "-normal, side: " << get_string(side) << std::endl;
    bounds.print_summary(os);
    return;
  }

  RealBounds bounds;
  Direction normal;
  MinMax side;
};

class IndexPlane{
public:
  IndexPlane(){};
  IndexPlane(IndexBounds ib, int pt, Direction n, MinMax sd)
  : bounds(ib)
  , normal(n)
  , side(sd)
  {
    bounds.min[normal] = pt;
    bounds.max[normal] = pt;
  }

  // slides the plane in the normal direction
  // by the amount 'val'
  IndexPlane offset(int val) const {
    return IndexPlane(bounds, bounds.min[normal] + int(side)*val, normal, side);
  }

  // truncate the plane on the given boundary
  // by the given amount of index points
  IndexPlane truncate(BoundaryLocation bl, int ninds){
    IndexBounds bout = bounds;
    if (boundary_side(bl)== MinMax::MIN){
      bout.min[boundary_direction(bl)] += ninds;
    } 
    else{
      bout.max[boundary_direction(bl)] -= ninds;
    }
    return IndexPlane(bout, bout.min[normal], normal, side);
  }

  unsigned int get_npts() const{
    unsigned int npts=1;
    if (normal!=Direction::X) npts*=(bounds.max[Direction::X]-bounds.min[Direction::X]+1);
    if (normal!=Direction::Y) npts*=(bounds.max[Direction::Y]-bounds.min[Direction::Y]+1);
    if (normal!=Direction::Z) npts*=(bounds.max[Direction::Z]-bounds.min[Direction::Z]+1);
    return npts;
  }

  void print_summary(std::ostream & os=std::cout) const{
    os << get_string(normal) << "-normal, side: " << get_string(side) << std::endl;
    bounds.print_summary(os);
    return;
  }

  IndexBounds bounds;
  Direction normal;
  MinMax side;
};


class ProcGrid{
public:
  ProcGrid(){};

  ProcGrid(RealBounds rb, 
           int numx, int numy, int numz,
           double dx, double dy, double dz)
  : nx(numx)
  , ny(numy)
  , nz(numz)
  {
    // construct the grid
    grid.resize(nx);
    for (auto i=0; i<nx; i++) grid[i].resize(ny);
    for (auto i=0; i<nx; i++){
      for (auto j=0; j<ny; j++) grid[i][j].resize(nz);
    }
    for (auto i=0; i<nx; i++){
      for (auto j=0; j<ny; j++){
        for (auto k=0; k<nz; k++){
          grid[i][j][k] = k*nx*ny + j*nx + i;
        }
      }
    }

    // construct the realbounds
    RealPoint min = rb.min;
    RealPoint max = rb.max;
    double xcut, ycut, zcut;
    xcut = (max.x-min.x)/nx;
    ycut = (max.y-min.y)/ny;
    zcut = (max.z-min.z)/nz;
    double xmin, ymin, zmin;
    double xmax, ymax, zmax;

    // construct the grid
    int este;
    IndexPoint estep = IndexPoint(0,0,0);
    bounds.resize(nx);
    for (auto i=0; i<nx; i++) bounds[i].resize(ny);
    for (auto i=0; i<nx; i++){
      for (auto j=0; j<ny; j++) bounds[i][j].resize(nz);
    }
    for (auto i=0; i<nx; i++){
      for (auto j=0; j<ny; j++){
        for (auto k=0; k<nz; k++){

          este = grid[i][j][k];

          // determine which x, y, and z slice this is
          estep = get_index_proc(este);

          xmin = min.x + estep.i*xcut;
          xmax = xmin + xcut;
          ymin = min.y + estep.j*ycut;
          ymax = ymin + ycut;
          zmin = min.z + estep.k*zcut;
          zmax = zmin + zcut;

          // if (estep.i>0 && estep.i<nx) xmin -= dx;
          // if (estep.i>=0 && estep.i<nx-1)  xmax += dx;
          // if (estep.j>0 && estep.j<ny) ymin -= dy;
          // if (estep.j>=0 && estep.j<ny-1)  ymax += dy;
          // if (estep.k>0 && estep.k<nz) zmin -= dz;
          // if (estep.k>=0 && estep.k<nz-1)  zmax += dz;


          bounds[i][j][k] = RealBounds(RealPoint(xmin,ymin,zmin),RealPoint(xmax,ymax,zmax));

        }
      }
    }
  }

  RealBounds get_bounds_proc(int thisproc) const{
    IndexPoint ip = get_index_proc(thisproc);
    return bounds[ip.i][ip.j][ip.k];
  }

  bool contains_point(RealPoint rp, int proc) const{
    RealBounds rb = get_bounds_proc(proc);
    return rb.contains_point(rp);
  }

  IndexPoint get_index_proc(int thisproc) const{
    int zc, yc, xc;
    zc = thisproc/(nx*ny);
    yc = (thisproc-(zc*nx*ny))/nx;
    xc = (thisproc-(zc*nx*ny)-(yc*nx));
    return IndexPoint(xc,yc,zc);
  }

  int get_neighbor_proc(int thisproc, BoundaryLocation bl) const{
    IndexPoint tp = get_index_proc(thisproc);
    IndexPoint off = boundary_offset(bl);
    return get_proc_at(tp+off);
  }

  int get_proc_at(IndexPoint ip) const{
    return grid[(ip.i<0? nx-1 : ip.i%nx)][(ip.j<0? ny-1 : ip.j%ny)][(ip.k<0? nz-1 : ip.k%nz)];
  }

  void print_grid(std::ostream & os=std::cout) const{
    os << "Size (" << nx << " x " << ny << " x " << nz << ")" << std::endl;
    for (auto i=0; i<nx; i++){
      for (auto j=0; j<ny; j++){
        for (auto k=0; k<nz; k++){
          os << "(" << i << ", " << j << ", " << k << ") ===> Proc #" << get_proc_at(IndexPoint(i,j,k)) << std::endl;
        }
      }
    }
  }

  bool is_proc_on_boundary(int proc, BoundaryLocation bl) const{
    IndexPoint ploc = get_index_proc(proc);
    Direction d = boundary_direction(bl);
    int cval;
    std::vector<std::reference_wrapper<const int>> np = {std::reference_wrapper<const int>(nx),
							 std::reference_wrapper<const int>(ny),
							 std::reference_wrapper<const int>(nz)};
    if (boundary_side(bl) == MinMax::MIN) cval = 0;
    else if (boundary_side(bl) == MinMax::MAX) cval = np[(int)d]-1;
    return (ploc[d] == cval? true : false);
  }

  // split a volume into a number of processors
  static IndexPoint split_volume(RealBounds rb, int nprocs_tot,
                                 double dx, double dy, double dz){
    if (nprocs_tot==1) return IndexPoint(1,1,1);

    // if (nprocs_tot%2==1){
    //   std::cout << "ERROR: cannot split an odd number of processors" << std::endl;
    //   throw -1;
    // }
    std::cout << pow(nprocs_tot, 1.0/3.0) << std::endl;
    if (nprocs_tot == pow(int(pow(nprocs_tot, 1.0/3.0)),3.0)){
      int n3 = pow(nprocs_tot, 1.0/3.0);
      return IndexPoint(n3, n3, n3);
    }

    double nptsx = (rb.max.x - rb.min.x)/dx;
    double nptsy = (rb.max.y - rb.min.y)/dy;
    double nptsz = (rb.max.z - rb.min.z)/dz;

    int numz = nptsz/std::max(nptsx, nptsy) + 1;
    //int numx = sqrt(nptsx/nptsy*nprocs_tot/numz);
    int numy = sqrt(nptsy/nptsx*nprocs_tot/numz);
    int numx = nprocs_tot/(numz*numy);

    return IndexPoint(numx, numy, numz);
  }

  int nx, ny, nz;
  std::vector<std::vector<std::vector<int>>> grid; // maps a proc grid location to a proc
  std::vector<std::vector<std::vector<RealBounds>>> bounds;  // maps a proc grid location to a volume bound
  std::vector<std::reference_wrapper<int>> nprocs = {std::reference_wrapper<int>(nx),
						     std::reference_wrapper<int>(ny),
						     std::reference_wrapper<int>(nz)};
};



// regular mesh
class RegularMesh : public Mesh {
public:

  RegularMesh(); 
  RegularMesh(const Mesh & mesh);
  //RegularMesh(const RegularMesh & rmesh);
  //~RegularMesh();

  // inspectors
  MeshNode & regular_node(unsigned int i, unsigned int j=0, unsigned int k=0);
  unsigned int reg_num_nodes_x() const {return _num_nodes_x;};
  unsigned int reg_num_nodes_y() const {return _num_nodes_y;};
  unsigned int reg_num_nodes_z() const {return _num_nodes_z;};
  unsigned int num_cells_x() const {return _num_nodes_x-1;};
  unsigned int num_cells_y() const {return _num_nodes_y-1;};
  unsigned int num_cells_z() const {return _num_nodes_z-1;};
  double res() const {return _res;};
  //unsigned int & reg_nodes_boundary(BoundaryLocation loc);
  //CoreBoundaries core_boundaries(unsigned int core) const {return m_core_boundaries[core];};
  std::vector<unsigned int> left_inds() const;
  std::vector<unsigned int> right_inds() const;
  std::vector<unsigned int> top_inds() const;
  std::vector<unsigned int> bottom_inds() const;
  std::vector<unsigned int> front_inds() const;
  std::vector<unsigned int> back_inds() const;

  unsigned int cell_map_to_global_ind(unsigned int i, unsigned int j=0, unsigned int k=0) const;
  
  unsigned int reg_inds_to_glob_ind(unsigned int i, unsigned int j=0, unsigned int k=0) const;
  void glob_ind_to_reg_inds(unsigned int glob_ind, unsigned int & i, unsigned int & j, unsigned int & k) const;
  unsigned int nearest_node(double x_loc, double y_loc=0.0, double z_loc=0.0) const;
  unsigned int nearest_node(RealPoint rp) const {return nearest_node(rp.x, rp.y, rp.z);};
  unsigned int nearest_element(double x_loc, double y_loc=0.0, double z_loc=0.0) const;
  unsigned int nearest_element(RealPoint rp) const {return nearest_element(rp.x, rp.y, rp.z);};

  // grid generation and refinement
  static RegularMesh create_regular_grid(double res, RealBounds rb);
  static RegularMesh create_regular_grid_n(double res, unsigned int num_nodes_x, unsigned int num_nodes_y = 1, 
                      unsigned int num_nodes_z = 1); // create a regular grid of points and store it in the mesh
  static RegularMesh create_regular_grid_b(double res, double xmin, double xmax, double ymin=0.0, double ymax=0.0,
                      double zmin=0.0, double zmax=0.0);
  // void assign_core_groups(unsigned int ncores);

private:
  unsigned int _num_nodes_x, _num_nodes_y, _num_nodes_z;
  double _res;
  // std::vector<CoreBoundaries> m_core_boundaries;

  void create_regular_grid_internal(double res, unsigned int num_nodes_x, unsigned int num_nodes_y, 
                      unsigned int num_nodes_z,
                      double xcen=0.0, double ycen=0.0, double zcen=0.0);


};



#endif
