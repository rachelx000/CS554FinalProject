/*

Data structures for learnply

Eugene Zhang 2005

*/

#ifndef __POLYHEDRON_H__
#define __POLYHEDRON_H__

#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include <vector>
#include <Eigen/Dense>

const double EPS = 1.0e-6;
const double PI=3.1415926535898;

// forward declarations
class Vertex;
class Edge;
class Triangle;
class Corner;

// The vertex class
class Vertex 
{
public:
    double x,y,z;
    int index;
    
    int ntris;
    Triangle **tris;
    int max_tris;
    
    /* Project 1, Problem 2 */
    int ncorners;
    Corner **corners;
    
    /* Project 1, Problem 3 */
    int valence_deficit;
    double K; // gaussian curvature

    icVector3 normal;
    void *other_props;

    /* curvature */
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	double voronoi_area;
	double mean_curvature;
	double gauss_curvature;
    icVector3 local_x_axis;
	icVector3 local_y_axis;
    Eigen::Matrix2d local_tensor;
    Eigen::Matrix3d global_tensor;
    icVector3 pcurve_major;
	icVector3 pcurve_minor;

public:
    Vertex(double xx, double yy, double zz)
    {
        x = xx; y = yy; z = zz;
    
        /* Project 1, Problem 2 */
        ncorners = 0;
        corners = NULL;
    }
};

// The edge class
class Edge 
{
public:
    int index;
    Vertex *verts[2];
    int ntris;
    Triangle **tris;
    double length;

    /* Project 1, Problem 2 */
    Corner *corners[2];

	Vertex* other_vertex(Vertex* v)
	{
        if (verts[0] == v)
            return verts[1];
        else if (verts[1] == v)
            return verts[0];
        else
            return nullptr;
	}

    /* Final Project */
    icVector3 tangent_v1, tangent_v2;   // tangent vectors for two endpoints
    icVector3 silhouette_point;
    icVector3 silhouette_point_normal;
    bool is_new = false;
};

// The triangle class
class Triangle 
{
public:
    int index;
    int nverts;
    Vertex *verts[3];
    Edge *edges[3];

    /* Project 1, Problem 2 */
    Corner *corners[3];

    double angle[3];
    float area;

    icVector3 normal;
    void *other_props;

	/* curvature */
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    icVector3 local_x_axis;
    icVector3 local_y_axis;
    Eigen::Matrix2d local_tensor;
    Eigen::Matrix3d global_tensor;
    icVector3 pcurve_major;
    icVector3 pcurve_minor;

};

// The corner class
class Corner
{
public:
	int index;

	Edge *e;
	Vertex *v;
	Triangle *t;

	Corner *o; // opposite corner
	Corner *n; // next corner
	Corner *p; // previous corner

    double angle;
};

// The vertex list class
class VertexList
{
public:

    int num_verts;
    int max_verts;
    Vertex **verts;

    VertexList (int max) 
    {
        max_verts = max;
        verts = new Vertex *[max_verts];
        num_verts = 0;
    }
    void finalize() 
    {
        free(verts);
        free(this);
    }
    void append(Vertex *v)
    {
        int i;

        /* first make sure there is enough room for new vertex */
        if (num_verts >= max_verts) {
            max_verts += 10;
            Vertex **tlist = new Vertex *[max_verts];
            for (i = 0; i < num_verts; i++)
                tlist[i] = verts[i];
            delete (verts);
            verts = tlist;
        }

        /* add new vertex to list */
        verts[num_verts] = v;
        num_verts++;
    }
};

// The edge list class
class EdgeList 
{
public:

    int num_edges;
    int max_edges;
    Edge **edges;
    
    EdgeList (int max) 
    {
    max_edges = max;
    edges = new Edge *[max_edges];
    num_edges = 0;
    }
    void finalize() 
    {
        free(edges);
        free(this);
    }
    void append(Edge *e)
    {
        int i;

        /* first make sure there is enough room for new edge */
        if (num_edges >= max_edges) {
            max_edges += 10;
            Edge **tlist = new Edge *[max_edges];
            for (i = 0; i < num_edges; i++)
                tlist[i] = edges[i];
            delete (edges);
            edges = tlist;
        }

        /* add new edge to list */
        edges[num_edges] = e;
        num_edges++;
    }
};

// The corner list class
class CornerList 
{
public:
    int num_corners;
    int max_corners;
    Corner **corners;

    CornerList(int max)
    {
        max_corners = max;
        corners = new Corner *[max_corners];
        num_corners = 0;
    }
    void finalize() 
    {
        free(corners);
        free(this);
    }
    void append(Corner *c)
    {
        int i;

        if (num_corners >= max_corners)
        {
            max_corners += 16;
            Corner **clist = (Corner **)realloc(corners, sizeof(Corner *) * max_corners);
            if (!clist)
                abort();
            corners = clist;
        }

        corners[num_corners] = c;
        num_corners++;
    }
};

struct CornerTableEntry
{
	int c;
	int v_min;
	int v_max;
	int o;
};

struct LineSegment
{
    icVector3 start;
    icVector3 end;

	LineSegment(icVector3 s, icVector3 e)
	{
		start = s;
		end = e;
	}
    LineSegment(double sx, double sy, double sz, double ex, double ey, double ez)
    {
		start = icVector3(sx, sy, sz);
		end = icVector3(ex, ey, ez);
    }
};

/* Final Project */
struct RemeshData {
    std::vector<Vertex*> new_vertices;
    std::vector<Edge*> new_edges;
    std::vector<Triangle*> new_triangles;

    void clear() {
        new_vertices.clear();
        new_edges.clear();
        new_triangles.clear();
    }
};

// The polyhedron class
class Polyhedron 
{
public:

    int index;

    Triangle **tlist;   /* list of triangles */
    int ntris;
    int max_tris;

    Vertex **vlist;     /* list of vertices */
    int nverts;
    int max_verts;

    Edge **elist;       /* list of edges */
    int nedges;
    int max_edges;

    /* Final Project */
    RemeshData remesh_data;

    /* Project 1, Problem 2 */
    Corner **clist;     /* list of corners */
    int ncorners;
    int max_corners;
    CornerTableEntry *ctable;
    
    /* Project 1, Problem 3 */
    double K_min, K_max;
    int min_valence_deficit;

    icVector3 center;
    double radius;
    double area;
    int seed;

    PlyOtherProp *vert_other,*face_other;

    /* Project 2, Problem 1 */
    std::vector<LineSegment> silhouette;
    
    /* Project 2, Problem 2 */
    double sum_mean_curvature, min_mean_curvature, max_mean_curvature;
    double sum_gauss_curvature, min_gauss_curvature, max_gauss_curvature;
    
	/* Project 2, Problem 3 */
    std::vector<LineSegment> major_hatches;
	std::vector<LineSegment> minor_hatches;

    // constructors
    Polyhedron();
    Polyhedron(FILE *);
    
	// utility functions
    void initialize();
    void finalize();
    void write_file(FILE *);
    void create_pointers();
    void create_edge(Vertex *, Vertex *);
    void create_edges();
    int face_to_vertex_ref(Triangle *, Vertex *);
    void order_vertex_to_tri_ptrs(Vertex *);
    void vertex_to_tri_ptrs();
    Triangle *find_common_edge(Triangle *, Vertex *, Vertex *);
    Triangle *other_triangle(Edge *, Triangle *);
    void calc_bounding_sphere();
    int calc_face_normals_and_area();
    void calc_edge_length();
    
    // Project 1 functions
	/* 1c */
    void average_normals();
    /* 2 */
    void create_corners();
    /* 3a */
    void compute_euler_characteristic();
    /* 3b */
    void compute_gaussian_curvature_angle_deficit();
    /* 3c */
    void compute_gaussian_curvature_valence_deficit();
    /* 3d */
    void compute_handles();

    // Project 2 functions
    /* 1a */
    void compute_silhouette_edges(const icMatrix3x3& view, const icVector3& translate);
	/* 1b */
    void compute_silhouette_faces(const icMatrix3x3& view, const icVector3& translate);
    /* 2ab */
    void compute_vert_voronoi_areas();
    void compute_vert_mean_curvature();
	void compute_vert_gaussian_curvature();
    void compute_vert_curvature_tensor();
    void update_vert_global_tensors();
    void smooth_vert_curvature_tensors(int weight_scheme, double step_size, int iterations);
    void compute_vert_principal_curvatures();
	/* 3 */
    void compute_face_principal_curvatures();
    void build_curvature_hatch_lines(int principal_direction, int transition_scheme);
	void hatch_line_step(int principal_direction, int transition_scheme, int iterations, 
        Edge* current_edge, Triangle* current_face, icVector3& current_loc, icVector3& current_dir);

    /* Final Project */
    int find_odd_visibility(const bool visibility[3]);
    std::pair<Edge*, Edge*> compute_silhouette_bridges(Triangle* tri, const bool visibility[3], Vertex* odd_vertex);
    void compute_tangents_for_an_edge(Edge* edge);
    void compute_silhouette_point(Edge* edge, const icVector3& camera_pos);
    icVector3 hermite_interpolation(const icVector3& v1, const icVector3& v2, const icVector3& t1, const icVector3& t2, double u);
    icVector3 compute_silhouette_segment(Edge* e1, Edge* e2, float u);
    int compute_segment_num(Edge* v1, Edge* v2);
    Vertex* create_new_vertex(icVector3& pos, icVector3& normal);
    Edge* create_new_edge(Vertex* v1, Vertex* v2);
    void create_new_triangle(Vertex* v1, Vertex* v2, Vertex* v3, Edge* e1, Edge* e2, Edge* e3);
    void remesh_silhouette(const icMatrix3x3& view, const icVector3& translate);
};

#endif // __POLYHEDRON_H__