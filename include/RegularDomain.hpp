#ifndef _REGULARDOMAIN_H
#define _REGULARDOMAIN_H

#include "Domain.hpp"
#include "RegularMesh.hpp"

class RegularDomain : public Domain{
public:
	RegularDomain(){};
	
	RegularDomain(RealBounds rb, RealPoint fixpt, double dx, double dy, double dz, int nprocs)
	: m_dx(dx)
	, m_dy(dy)
	, m_dz(dz)
	{

		// make adjustments to the bounds
		// because we are centering on a fixed point
		m_bounds = rb;
		m_bounds.min.x = fixpt.x-double(int((fixpt.x-m_bounds.min.x)/dx))*dx;
		m_bounds.min.y = fixpt.y-double(int((fixpt.y-m_bounds.min.y)/dy))*dy;
		m_bounds.min.z = fixpt.z-double(int((fixpt.z-m_bounds.min.z)/dz))*dz;

		m_bounds.max.x = fixpt.x+double(int((m_bounds.max.x-fixpt.x)/dx))*dx;
		m_bounds.max.y = fixpt.y+double(int((m_bounds.max.y-fixpt.y)/dy))*dy;
		m_bounds.max.z = fixpt.z+double(int((m_bounds.max.z-fixpt.z)/dz))*dz;
	
		// create the procgrid
		IndexPoint nproc = ProcGrid::split_volume(m_bounds, nprocs, dx, dy, dz);
		m_pg = ProcGrid(m_bounds, nproc.i, nproc.j, nproc.k, dx, dy, dz);

	}

	RealBounds get_bounds() const{
		return m_bounds;
	}

	RealBounds get_proc_bounds(int proc) const{
		return m_pg.get_bounds_proc(proc);
	}

	IndexBounds get_proc_ibounds(int proc) const{
		int imin, jmin, kmin;
		int imax, jmax, kmax;
		RealBounds rb = m_pg.get_bounds_proc(proc);
		imin = (rb.min.x - m_bounds.min.x)/m_dx;
		jmin = (rb.min.y - m_bounds.min.y)/m_dy;
		kmin = (rb.min.z - m_bounds.min.z)/m_dz;
		imax = imin + (rb.max.x - rb.min.x)/m_dx;
		jmax = jmin + (rb.max.y - rb.min.y)/m_dy;
		kmax = kmin + (rb.max.z - rb.min.z)/m_dz;
		return IndexBounds(IndexPoint(imin,jmin,kmin), IndexPoint(imax,jmax,kmax));
	}

	IndexBounds get_proc_local_ibounds(int proc){
		IndexBounds ib = get_proc_ibounds(proc);
		ib.max = ib.max - ib.min;
		ib.min = ib.min - ib.min;
		
		return ib;
	}

	IndexBounds get_proc_local_ibounds(int proc, RealBounds rb){
		RealBounds pb = get_proc_bounds(proc);
		RealBounds lb = rb;
		IndexPoint minp, maxp;
		vector<Direction> ds = {Direction::X, Direction::Y, Direction::Z};
		vector<double> dxi = {m_dx, m_dy, m_dz};
		for (auto i=0; i<ds.size(); i++){
			lb.min[ds[i]] = max(rb.min[ds[i]], pb.min[ds[i]]);
			lb.max[ds[i]] = min(rb.max[ds[i]], pb.max[ds[i]]);
		}
		for (auto i=0; i<ds.size(); i++){
			if (dxi[(int)ds[i]]==0) {
				minp[ds[i]] = 0;
				maxp[ds[i]] = 0;
				continue;
			}
			minp[ds[i]] = (lb.min[ds[i]]-pb.min[ds[i]])/dxi[(int)ds[i]];
			maxp[ds[i]] = (lb.max[ds[i]]-pb.min[ds[i]])/dxi[(int)ds[i]];
		}

		return IndexBounds(minp, maxp, get_proc_local_ibounds(proc).get_volume_context());
	}

	IndexBounds get_proc_ibounds_cell(int proc) const{
		int imin, jmin, kmin;
		int imax, jmax, kmax;
		RealBounds rb = m_pg.get_bounds_proc(proc);
		imin = (rb.min.x - m_bounds.min.x)/m_dx;
		jmin = (rb.min.y - m_bounds.min.y)/m_dy;
		kmin = (rb.min.z - m_bounds.min.z)/m_dz;
		imax = imin + (rb.max.x - rb.min.x)/m_dx -1;
		jmax = jmin + (rb.max.y - rb.min.y)/m_dy -1;
		kmax = kmin + (rb.max.z - rb.min.z)/m_dz -1;
		return IndexBounds(IndexPoint(imin,jmin,kmin), IndexPoint(imax,jmax,kmax));
	}

	IndexBounds get_proc_local_ibounds_cell(int proc){
		IndexBounds ib = get_proc_ibounds_cell(proc);
		ib.max = ib.max - ib.min;
		ib.min = ib.min - ib.min;
		//ib.max = ib.max - IndexPoint(1,1,1);
		return ib;
	}

	bool is_proc_red(int proc) const{
		IndexPoint ip = m_pg.get_index_proc(proc);
		if ((ip.i+ip.j+ip.k)%2 == 0) return true;
		return false;
	}


	int get_proc_neighbor(int proc, BoundaryLocation bl) const{
		return m_pg.get_neighbor_proc(proc, bl);
	}

	IndexPoint max_index() const{
		int imax, jmax, kmax;
		imax = (m_bounds.max.x-m_bounds.min.x)/m_dx;
		jmax = (m_bounds.max.y-m_bounds.min.y)/m_dy;
		kmax = (m_bounds.max.z-m_bounds.min.z)/m_dz;
		return IndexPoint(imax, jmax, kmax);
	}

	bool is_proc_on_boundary(int proc, BoundaryLocation bl) const{
		return m_pg.is_proc_on_boundary(proc, bl);
	}

	bool contains_point(RealPoint rp, int proc) const{
		return m_pg.contains_point(rp, proc);
	}

	unsigned int array_node_offset(int proc) const{
		unsigned int offset = 0;
		for (auto i=0; i<proc; i++){
			offset += nnodes(i);
		}
		return offset;
	}

	unsigned int array_element_offset(int proc) const{
		unsigned int offset = 0;
		for (auto i=0; i<proc; i++){
			offset += nelements(i);
		}
		return offset;
	}

	unsigned int nnodes(int proc) const{
		IndexBounds ib = get_proc_ibounds(proc);
		unsigned int nx, ny, nz;
		nz = (ib.max.k-ib.min.k+1);
		ny = (ib.max.j-ib.min.j+1);
		nx = (ib.max.i-ib.min.i+1);
		return nx*ny*nz;
	}

	unsigned int nelements(int proc) const{
		IndexBounds ib = get_proc_ibounds(proc);
		unsigned int nx, ny, nz;
		nz = (ib.max.k-ib.min.k);
		ny = (ib.max.j-ib.min.j);
		nx = (ib.max.i-ib.min.i);
		if (nz == 0) nz=1;
		if (ny == 0) ny=1;
		return nx*ny*nz;
	}

	unsigned int nnodes() const{
		int nproc = m_pg.nx*m_pg.ny*m_pg.nz;
		unsigned int ntot=0;
		for (auto i=0; i<nproc; i++){
			ntot += nnodes(i);
		}
		return ntot;
	}

	unsigned int nelements() const{
		int nproc = m_pg.nx*m_pg.ny*m_pg.nz;
		unsigned int ntot=0;
		for (auto i=0; i<nproc; i++){
			ntot += nelements(i);
		}
		return ntot;
	}

	void print_grid(ostream & os=std::cout) const{
		m_pg.print_grid(os);
		for (auto i=0; i<mpi::size(); i++){
			cout << "Proc #" << i << " ---> " << (this->is_proc_red(i)? "RED" : "BLACK") << endl;
		}
	}

private:
	RealBounds m_bounds;// = RealBounds(RealPoint(0,0,0), RealPoint(0,0,0));
	ProcGrid m_pg;
	double m_dx, m_dy, m_dz;
};


#endif