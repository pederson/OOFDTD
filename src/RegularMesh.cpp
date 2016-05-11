#include "RegularMesh.hpp"

using namespace std;

IndexBounds IndexVolume::get_bounds_offset(int off) const{
    return IndexBounds(IndexPoint(min(off, nx-1),min(off, ny-1),min(off, nz-1)),
                       IndexPoint(max(0, nx-off-1),max(0, ny-off-1),max(0,nz-off-1)),
                       *this);
}

IndexBounds IndexVolume::get_bounds_boundary(BoundaryLocation bl) const{
  // cout << "location: " << get_string(bl) << endl;
  // cout << "nx, ny, nz: " << nx << ", " << ny << ", " << nz << endl;
  int newval=0;
  vector<reference_wrapper<const int>> npts = {reference_wrapper<const int>(nx),
					       reference_wrapper<const int>(ny), 
					       reference_wrapper<const int>(nz)};
  if (boundary_side(bl)==MinMax::MIN) newval = 0;
  else newval = max(0,npts[(int)boundary_direction(bl)] -1);
  
  IndexPoint minpt = IndexPoint(0,0,0);
  IndexPoint maxpt = IndexPoint(nx-1,ny-1,nz-1);
  vector<reference_wrapper<int>> minp = {reference_wrapper<int>(minpt.i), 
					 reference_wrapper<int>(minpt.j), 
					 reference_wrapper<int>(minpt.k)};
  vector<reference_wrapper<int>> maxp = {reference_wrapper<int>(maxpt.i), 
					 reference_wrapper<int>(maxpt.j), 
					 reference_wrapper<int>(maxpt.k)};
  minp[(int)boundary_direction(bl)].get() = newval;
  maxp[(int)boundary_direction(bl)].get() = newval;


  return IndexBounds(minpt, maxpt, *this);
}

IndexPlane IndexVolume::get_boundary_plane(BoundaryLocation bl) const{
  Direction d = boundary_direction(bl);
  MinMax side = boundary_side(bl);
  IndexBounds ib = get_bounds_boundary(bl);
  return IndexPlane(ib, ib.get_side(side)[d], boundary_direction(bl), boundary_side(bl));
}

RealPlane RealBounds::get_boundary_plane(BoundaryLocation bl) const{
  Direction d = boundary_direction(bl);
  MinMax side = boundary_side(bl);
  return RealPlane(*this, get_side(side)[d], boundary_direction(bl), boundary_side(bl));
}

IndexPlane IndexBounds::get_boundary_plane(BoundaryLocation bl) const{
  Direction d = boundary_direction(bl);
  MinMax side = boundary_side(bl);
  return IndexPlane(*this, get_side(side)[d], boundary_direction(bl), boundary_side(bl));
} 

RegularMesh::RegularMesh(){
  
}

RegularMesh::RegularMesh(const Mesh & mesh){
  _mesh_type = REGULAR;

  // find res
  double x0 = mesh.node(0).x();
  double xc = mesh.node(0).x();
  double xdiff, xdiffmin;
  xdiffmin = x0*100000000 + 10000000000;
  for (auto i=1; i<mesh.nodecount(); i++){
    xc = mesh.node(i).x();
    if (xc != x0){
      xdiff = xc - x0;
      if (xdiff < xdiffmin) {
        _res = xdiff;
        xdiffmin = xdiff;
      }
    }
  }

  // find numx numy and numz
  _num_nodes_x = (mesh.xmax() - mesh.xmin())/_res + 1;
  _num_nodes_y = (mesh.ymax() - mesh.ymin())/_res + 1;
  _num_nodes_z = (mesh.zmax() - mesh.zmin())/_res + 1;

  // copy everything from the mesh data
  _xmax = mesh.xmax();
  _xmin = mesh.xmin();
  _ymax = mesh.ymax();
  _ymin = mesh.ymin();
  _zmax = mesh.zmax();
  _zmin = mesh.zmin();

  _num_dims = mesh.num_dims();

  set_nodecount(mesh.nodecount());
  set_elementcount(mesh.elementcount());

  for (auto i=0; i<mesh.nodecount(); i++){
    _nodes[i] = mesh.node(i);
  }
  for (auto i=0; i<mesh.elementcount(); i++){
    _elements[i] = mesh.element(i);
  }



}

RegularMesh RegularMesh::create_regular_grid(double res, RealBounds rb){
  unsigned int num_nodes_x, num_nodes_y, num_nodes_z;
  RegularMesh mesh_out;

  num_nodes_x = (unsigned int)((rb.max.x-rb.min.x)/res) + 1;
  num_nodes_y = (unsigned int)((rb.max.y-rb.min.y)/res) + 1;
  num_nodes_z = (unsigned int)((rb.max.z-rb.min.z)/res) + 1;

  //cout << "nx: " << num_nodes_x << "  ny: " << num_nodes_y << "  nz: " << num_nodes_z << endl;
  double xcen, ycen, zcen;
  xcen = (rb.max.x + rb.min.x)/2.0;
  ycen = (rb.max.y + rb.min.y)/2.0;
  zcen = (rb.max.z + rb.min.z)/2.0;
  mesh_out.create_regular_grid_internal(res, num_nodes_x, num_nodes_y, num_nodes_z, xcen, ycen, zcen);

  return mesh_out;
}

RegularMesh RegularMesh::create_regular_grid_n(double res, unsigned int num_nodes_x, unsigned int num_nodes_y, 
                    unsigned int num_nodes_z){
  RegularMesh mesh_out;
  mesh_out.create_regular_grid_internal(res, num_nodes_x, num_nodes_y, num_nodes_z, 0.0, 0.0, 0.0);

  return mesh_out;
}


RegularMesh RegularMesh::create_regular_grid_b(double res, double xmin, double xmax, double ymin, double ymax,
                    double zmin, double zmax){
  unsigned int num_nodes_x, num_nodes_y, num_nodes_z;
  RegularMesh mesh_out;

  num_nodes_x = (unsigned int)((xmax-xmin)/res) + 1;
  num_nodes_y = (unsigned int)((ymax-ymin)/res) + 1;
  num_nodes_z = (unsigned int)((zmax-zmin)/res) + 1;

  //cout << "nx: " << num_nodes_x << "  ny: " << num_nodes_y << "  nz: " << num_nodes_z << endl;
  double xcen, ycen, zcen;
  xcen = (xmin + xmax)/2.0;
  ycen = (ymin + ymax)/2.0;
  zcen = (zmin + zmax)/2.0;
  mesh_out.create_regular_grid_internal(res, num_nodes_x, num_nodes_y, num_nodes_z, xcen, ycen, zcen);

  return mesh_out;
}



void RegularMesh::create_regular_grid_internal(double res, unsigned int num_nodes_x, unsigned int num_nodes_y, 
                      unsigned int num_nodes_z,
                      double xcen, double ycen, double zcen){
  // declare vars

  // do some input checking

  // initialize things and fill in metadata
  _res = res;
  _num_nodes_x = num_nodes_x;
  _num_nodes_y = num_nodes_y;
  _num_nodes_z = num_nodes_z;
  _mesh_type = REGULAR;
  _nodes.resize(num_nodes_x*num_nodes_y*num_nodes_z);
  if (num_nodes_y==1 && num_nodes_z==1){
    _num_dims = 1;
  }
  else if (num_nodes_z==1){
    _num_dims = 2;
  }
  else{
    _num_dims = 3;
  }

  // create nodes
  unsigned int glidx, i, j, k;
  for (i=0; i<num_nodes_z; i++){
    for (j=0; j<num_nodes_y; j++){
      for (k=0; k<num_nodes_x; k++){
        //glidx = i*(num_nodes_x*num_nodes_y) + j*(num_nodes_x) + k;
        glidx = reg_inds_to_glob_ind(k,j,i);
        _nodes.at(glidx).set_z(double(i)*res);
        _nodes.at(glidx).set_y(double(j)*res);
        _nodes.at(glidx).set_x(double(k)*res);
      }
    }
  }
  // set boundary properties for nodes
  if (num_nodes_y==1 && num_nodes_z==1){
    i=0; j=0;
    for (k=0; k<num_nodes_x; k++){
      _nodes.at(k).set_num_connections(2);
    }
    _nodes.at(0).set_boundary(true);
    _nodes.at(0).set_num_connections(1);
    _nodes.at(num_nodes_x-1).set_boundary(true);
    _nodes.at(num_nodes_x-1).set_num_connections(1);
  }
  else if (num_nodes_z==1){
    i=0;
    for (j=0; j<num_nodes_y; j++){
      for (k=0; k<num_nodes_x; k++){
        glidx = i*(num_nodes_x*num_nodes_y) + j*(num_nodes_x) + k;
        if (k==0 || k==num_nodes_x-1 || j==0 || j==num_nodes_y-1){
          _nodes.at(glidx).set_boundary(true);
          _nodes.at(glidx).set_num_connections(3);
        }
      }
    }
    _nodes.at(0).set_num_connections(2);
    _nodes.at(num_nodes_x-1).set_num_connections(2);
    _nodes.at((num_nodes_y-1)*num_nodes_x).set_num_connections(2);
    _nodes.at((num_nodes_y-1)*num_nodes_x+num_nodes_x-1).set_num_connections(2);
  }
  else{
    for (i=0; i<num_nodes_z; i++){
      for (j=0; j<num_nodes_y; j++){
        for (k=0; k<num_nodes_x; k++){
          glidx = i*(num_nodes_x*num_nodes_y) + j*(num_nodes_x) + k;
          if (k==0 || k==num_nodes_x-1 || j==0 || j==num_nodes_y-1 || i==0 || i==num_nodes_z-1){
            _nodes.at(glidx).set_boundary(true);
            _nodes.at(glidx).set_num_connections(5);
          }
        }
      }
    }
    _nodes.at(0).set_num_connections(3);
    _nodes.at(num_nodes_x-1).set_num_connections(3);
    _nodes.at((num_nodes_y)*(num_nodes_x)-1).set_num_connections(3);
    _nodes.at((num_nodes_y-1)*(num_nodes_x)+1).set_num_connections(3);

    _nodes.at((num_nodes_z-1)*(num_nodes_x)*(num_nodes_y)-1).set_num_connections(3);
    _nodes.at((num_nodes_z-1)*(num_nodes_x)*(num_nodes_y)-1 + num_nodes_x-1).set_num_connections(3);
    _nodes.at((num_nodes_z-1)*(num_nodes_x)*(num_nodes_y)-1 + (num_nodes_y)*(num_nodes_x)-1).set_num_connections(3);
    _nodes.at((num_nodes_z)*(num_nodes_x)*(num_nodes_y)-1).set_num_connections(3);
  }
  

  // create elements
  unsigned int blf, tlf, brf, trf, blb, tlb, brb, trb, nex, ney, nez; 
  
  if (num_nodes_y==1 && num_nodes_z==1){
    nex = num_nodes_x-1;
    _elements.resize((num_nodes_x-1));
    i=0; j=0;
    for (k=0; k<_elements.size(); k++){
      glidx = i*(nex*ney) + j*(nex) + k;
      blf = k;
      brf = k+1;
      _elements.at(glidx).set_vertex_inds({blf, brf});
      _elements.at(glidx).set_element_type(ElementType::LINE_2);
    }
  }
  else if (num_nodes_z==1){
    nex = num_nodes_x-1;
    ney = num_nodes_y-1;
    _elements.resize((num_nodes_x-1)*(num_nodes_y-1));
    for (j=0; j<ney; j++){
      for (k=0; k<nex; k++){
        glidx = i*(nex*ney) + j*(nex) + k;
        blf = (nex+1)*(j) + k;
        brf = (nex+1)*(j) + k+1;
        trf = (nex+1)*(j+1) + k+1;
        tlf = (nex+1)*(j+1) + k;
        _elements.at(glidx).set_vertex_inds({blf, brf, trf, tlf});
        _elements.at(glidx).set_element_type(ElementType::QUAD_4);
      }
    }
  }
  else{
    _elements.resize((num_nodes_x-1)*(num_nodes_y-1)*(num_nodes_z-1));
    nex = num_nodes_x-1;
    ney = num_nodes_y-1;
    nez = num_nodes_z-1;
    for (i=0; i<nez; i++){
      for (j=0; j<ney; j++){
        for (k=0; k<nex; k++){
          glidx = i*(nex*ney) + j*(nex) + k;
          blf = (nex+1)*(ney+1)*(i) + (nex+1)*(j) + k;
          brf = (nex+1)*(ney+1)*(i) + (nex+1)*(j) + k+1;
          trf = (nex+1)*(ney+1)*(i) + (nex+1)*(j+1) + k+1;
          tlf = (nex+1)*(ney+1)*(i) + (nex+1)*(j+1) + k;
          blb = (nex+1)*(ney+1)*(i+1) + (nex+1)*(j) + k;
          brb = (nex+1)*(ney+1)*(i+1) + (nex+1)*(j) + k+1;
          trb = (nex+1)*(ney+1)*(i+1) + (nex+1)*(j+1) + k+1;
          tlb = (nex+1)*(ney+1)*(i+1) + (nex+1)*(j+1) + k;
          _elements.at(glidx).set_vertex_inds({blf, brf, trf, tlf, blb, brb, trb, tlb});
          _elements.at(glidx).set_element_type(ElementType::HEX_8);
        }
      }
    }
  }

  // translate the center to 0,0
  double xc, yc, zc;
  unsigned int centerind, ix, iy, iz; 
  if (num_nodes_x%2 == 1) ix = (num_nodes_x-1)/2;
  else ix = (num_nodes_x/2);
  if (num_nodes_y%2 == 1) iy = (num_nodes_y-1)/2;
  else iy = (num_nodes_y/2);
  if (num_nodes_z%2 == 1) iz = (num_nodes_z-1)/2;
  else iz = (num_nodes_z/2);
  centerind = reg_inds_to_glob_ind(ix, iy, iz);
  xc = _nodes.at(centerind).x();
  yc = _nodes.at(centerind).y();
  zc = _nodes.at(centerind).z();
  // xc = ((num_nodes_x-1)*res)/2.0;
  // yc = ((num_nodes_y-1)*res)/2.0;
  // zc = ((num_nodes_z-1)*res)/2.0;

  // translate to the correct centerpoint
  for (auto i=0; i<_nodes.size(); i++){
    //nd = _nodes.at(i);
    _nodes.at(i).set_x(_nodes.at(i).x() + xcen - xc);
    _nodes.at(i).set_y(_nodes.at(i).y() + ycen - yc);
    _nodes.at(i).set_z(_nodes.at(i).z() + zcen - zc);
  }
  calc_extents();

  return;
}

unsigned int RegularMesh::cell_map_to_global_ind(unsigned int i, unsigned int j, unsigned int k) const {
  return k*(_num_nodes_x-1)*(_num_nodes_y-1) + j*(_num_nodes_x-1) + i;
}

unsigned int RegularMesh::reg_inds_to_glob_ind(unsigned int i, unsigned int j, unsigned int k) const {
  return k*(_num_nodes_x*_num_nodes_y) + j*(_num_nodes_x) + i;
}

void RegularMesh::glob_ind_to_reg_inds(unsigned int glob_ind, unsigned int & i, unsigned int & j, unsigned int & k) const {
  k = glob_ind % (_num_nodes_z);
  j = glob_ind % (_num_nodes_y*_num_nodes_z);
  i = glob_ind % (_num_nodes_x*_num_nodes_y*_num_nodes_z);

  return;
}

unsigned int RegularMesh::nearest_node(double x_loc, double y_loc, double z_loc) const{

  unsigned int xind, yind, zind;
  xind = (unsigned int)(double(_num_nodes_x-1)*(x_loc-_xmin)/(_xmax-_xmin));
  if (_num_nodes_y==1) yind = 0;
  else yind = (unsigned int)(double(_num_nodes_y)*(y_loc-_ymin)/(_ymax-_ymin));
  if (_num_nodes_z==1) zind = 0;
  else zind = (unsigned int)(double(_num_nodes_z)*(z_loc-_zmin)/(_zmax-_zmin));
  return reg_inds_to_glob_ind(xind, yind, zind);
}

unsigned int RegularMesh::nearest_element(double x_loc, double y_loc, double z_loc) const{

  unsigned int xind, yind, zind;
  xind = (unsigned int)(double(_num_nodes_x-1)*(x_loc-_xmin)/(_xmax-_xmin));
  if (_num_nodes_y==1) yind = 0;
  else yind = (unsigned int)(double(_num_nodes_y)*(y_loc-_ymin)/(_ymax-_ymin));
  if (_num_nodes_z==1) zind = 0;
  else zind = (unsigned int)(double(_num_nodes_z)*(z_loc-_zmin)/(_zmax-_zmin));
  return cell_map_to_global_ind(xind, yind, zind);
}
