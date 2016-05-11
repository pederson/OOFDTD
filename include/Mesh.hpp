#ifndef _MESH_H
#define _MESH_H

#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <fstream>

#include "mpitools.hpp"

// these enums follow closely to the MSH format
// http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php
enum MeshType{REGULAR=0, UNSTRUCTURED_TRI=1, UNSTRUCTURED_QUAD=2, MESH_UNKNOWN};
enum class ElementType : unsigned int {EMPTY_0=0, POINT_1, LINE_2, TRI_3, QUAD_4, TET_4, HEX_8, PRISM_6, PYRAMID_5, UNKNOWN};

const std::vector<std::string> elnames = {"Empty", "Point_1", "Line_2", "Tri_3", "Quad_4", "Tet_4", "Hex_8", "Prism_6", "Pyramid_5", "Unknown"};
const std::vector<unsigned int> nvert = {0, 1, 2, 3, 4, 4, 8, 6, 5, 0};
inline std::string get_string(ElementType et){return elnames[(int)et];};
// inline ElementType get_element(std::string els){
//   for (auto e=ElementType::EMPTY_0; e<ElementType::UNKNOWN; e++){
//     if (get_string(e).compare(els)==0) return e;
//   }
//   return ElementType::UNKNOWN;
// }
inline unsigned int get_nvert(ElementType et){return nvert[(int)et];};


// forward declarations
class MeshNode;
class MeshElement;


class Mesh{
public:
  Mesh();   // constructor
  Mesh(const Mesh & mesh);  // copy constructor
  ~Mesh();    // destructor

  void print_summary(std::ostream & os=std::cout) const;
  void print_detailed(std::ostream & os=std::cout) const;
  

  // metadata access
  MeshType mesh_type() const {return _mesh_type;};
  unsigned int num_dims() const {return _num_dims;};
  unsigned int nodecount() const {return _nodes.size();};
  unsigned int elementcount() const {return _elements.size();};
  double xmin() const {return _xmin;};
  double ymin() const {return _ymin;};
  double zmin() const {return _zmin;};
  double xmax() const {return _xmax;};
  double ymax() const {return _ymax;};
  double zmax() const {return _zmax;};
  unsigned int nearest_node(double x_loc, double y_loc=0.0, double z_loc=0.0) const;

  // node and element access
  MeshNode & node(unsigned int i) {return _nodes.at(i);};
  MeshElement & element(unsigned int i) {return _elements.at(i);};
  const MeshNode & node(unsigned int i) const {return _nodes.at(i);};
  const MeshElement & element(unsigned int i) const {return _elements.at(i);};
  
  // property interaction and access
  const double & x();
  const double & y();
  const double & z();
  const bool & boundary();
  const unsigned int & core_group();
  const unsigned int & num_connections();
  const double & data(std::string fieldname) const;

  const double * nodedata(std::string fieldname) const;
  const double * elementdata(std::string fieldname) const;

  void set_nodecount(unsigned int count);
  void set_elementcount(unsigned int count);
  void set_num_dims(unsigned int ndims) {_num_dims = ndims;};

  void add_phys_property(std::string property_name, const double * property_vals);
  void add_phys_property(std::string proprety_name, double init_val);
  std::vector<std::string> get_phys_properties() const {return _phys_property_names;};

  void add_nodedata(std::string property_name, const double * values);
  void add_elementdata(std::string property_name, const double * values);

  void reset_property(std::string property_name, double reset_val=0.0);
  void set_phys_property(std::string property_name, unsigned int i, double val){_phys_properties.at(property_name).at(i) = val;};
  void increment_phys_property(std::string property_name, unsigned int i, double val){_phys_properties.at(property_name).at(i) += val;};

  // grid generation and refinement
  //static Mesh create_unstructured_tri_simple();

  // reading and writing files
  static Mesh read_MSH(std::string filename);
  static Mesh read_XML(std::string filename);
  static Mesh read_NEU(std::string filename, unsigned int byte_offset=0);
  static Mesh read_CAS(std::string filename, unsigned int byte_offset=0);
  void write_MSH(std::string filename) const;
  void write_NEU(std::string filename) const;
  void write_CAS(std::string filename) const;

  void calc_extents();

//private:
protected:
  // metadata
  MeshType _mesh_type;
  unsigned int _num_dims;
  double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;
  

  // nodes and elements
  std::vector<MeshNode> _nodes; // array of nodes
  std::vector<MeshElement> _elements; // array of elements

  // user-defined properties for the mesh
  // DYLAN_TODO: RENAME THIS NodeData and ElementData
  std::vector<std::string> _phys_property_names; // the name position in this vector corresponds with the position of the property value in the node
  std::map<std::string, std::vector<double>> _phys_properties;
  
  std::map<std::string, std::vector<double>> m_nodedata;
  std::map<std::string, std::vector<double>> m_elementdata;

  // other properties conveniently placed in arrays (on demand) in order to return data
  std::vector<bool> _boundary;
  std::vector<double> _x, _y, _z;
  std::vector<unsigned int> _core_group, _num_connections;

  void read_MSH_internal(std::string filename);

  

};


class MeshElement{
public:
  MeshElement();
  MeshElement(std::vector<unsigned int> vertex_inds);
  ~MeshElement();

  // utils
  void print_summary(std::ostream & os=std::cout) const;
  void print_detailed(std::ostream & os=std::cout) const;
  double area() const;
  double perimeter() const;

  // member data access
  ElementType type() const {return _element_type;};
  unsigned int num_vertices() const {return _vertex_inds.size();};
  unsigned int vertex_ind(unsigned int i) const {return _vertex_inds.at(i);};
  //std::vector<unsigned int> vertex_inds() const {return _vertex_inds;};
  const unsigned int & vertex_inds() const {return _vertex_inds.front();};

  // mutators
  void remove_vertex(unsigned int vert_ind);
  void add_vertex(unsigned int vert_ind, int position=-1);
  void set_vertex_inds(std::vector<unsigned int> vertex_inds) {_vertex_inds = vertex_inds;};
  void set_element_type(ElementType type) {_element_type = type;};

private:
  std::vector<unsigned int> _vertex_inds;
  ElementType _element_type; // 

  void recalc_type();

};





class MeshNode{
public:
  MeshNode();
  MeshNode(double x, double y, double z=0.0, bool boundary=false, unsigned int num_connections=0, unsigned int core_group=0);
  ~MeshNode();

  // utils
  void print_summary(std::ostream & os=std::cout) const;
  void print_detailed(std::ostream & os=std::cout) const;
  static double dist_sq(const MeshNode & node1, const MeshNode & node2);
  static double dist(const MeshNode & node1, const MeshNode & node2) {return sqrt(dist_sq(node1, node2));};
  static double area_tri(const MeshNode & node1, const MeshNode & node2, const MeshNode & node3);

  // inspectors
  double x() const {return _x;};
  double y() const {return _y;};
  double z() const {return _z;};
  unsigned int num_connections() const {return _num_connections;};
  unsigned int core_group() const {return _core_group;};
  bool boundary() const {return _boundary;};

  // mutators
  void add_connection() {_num_connections++;};
  void remove_connection() {_num_connections--;};
  void set_x(double xval) {_x = xval;};
  void set_y(double yval) {_y = yval;};
  void set_z(double zval) {_z = zval;};
  void set_boundary(bool boundary) {_boundary = boundary;};
  void set_core_group(unsigned int core_group) {_core_group = core_group;};
  void set_num_connections(unsigned int num_connections) {_num_connections = num_connections;};

private:
  double _x, _y, _z; 
  bool _boundary;
  unsigned int _core_group, _num_connections;

};


#endif
