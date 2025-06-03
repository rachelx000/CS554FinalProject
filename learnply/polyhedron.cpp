/*

Data structures for learnply

Eugene Zhang, 2005

*/


#define _USE_MATH_DEFINES
#include <cassert>
#include "polyhedron.h"
#include "learnply_io.h"
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <map>

static PlyFile *in_ply;



;//////////////////////////////////////////////////////////////////////////////
;// Constructors
;//////////////////////////////////////////////////////////////////////////////

// read in a polyhedron from a file
Polyhedron::Polyhedron(FILE *file)
{
	int i,j;
	int elem_count;
	char *elem_name;

	/*** Read in the original PLY object ***/
	in_ply = read_ply (file);

	for (i = 0; i < in_ply->num_elem_types; i++) 
	{
		/* prepare to read the i'th list of elements */
		elem_name = setup_element_read_ply (in_ply, i, &elem_count);

		if (equal_strings ("vertex", elem_name)) {

			/* create a vertex list to hold all the vertices */
			nverts = max_verts = elem_count;
			vlist = new Vertex *[nverts];

			/* set up for getting vertex elements */
			setup_property_ply (in_ply, &vert_props[0]);
			setup_property_ply (in_ply, &vert_props[1]);
			setup_property_ply (in_ply, &vert_props[2]);
			vert_other = get_other_properties_ply (in_ply, 
								offsetof(Vertex_io,other_props));

			/* grab all the vertex elements */
			for (j = 0; j < nverts; j++) {
				Vertex_io vert;
				get_element_ply (in_ply, (void *) &vert);

				/* copy info from the "vert" structure */
				vlist[j] = new Vertex (vert.x, vert.y, vert.z);
				vlist[j]->other_props = vert.other_props;
			}
		}
		else if (equal_strings ("face", elem_name)) {

			/* create a list to hold all the face elements */
			ntris = max_tris = elem_count;
			tlist = new Triangle *[ntris];

			/* set up for getting face elements */
			setup_property_ply (in_ply, &face_props[0]);
			face_other = get_other_properties_ply (in_ply, offsetof(Face_io,other_props));

			/* grab all the face elements */
			for (j = 0; j < elem_count; j++) {
			Face_io face;
			get_element_ply (in_ply, (void *) &face);

			if (face.nverts != 3) {
				fprintf (stderr, "Face has %d vertices (should be three).\n",
						face.nverts);
				exit (-1);
			}

			/* copy info from the "face" structure */
			tlist[j] = new Triangle;
			tlist[j]->nverts = 3;
			tlist[j]->verts[0] = (Vertex *)(intptr_t)face.verts[0];
			tlist[j]->verts[1] = (Vertex *)(intptr_t)face.verts[1];
			tlist[j]->verts[2] = (Vertex *)(intptr_t)face.verts[2];
			tlist[j]->other_props = face.other_props;
			}
		}
		else
			get_other_element_ply (in_ply);
	}

	/* close the file */
	close_ply (in_ply);

	/* fix up vertex pointers in triangles */
	for (i = 0; i < ntris; i++) 
	{
		tlist[i]->verts[0] = vlist[(intptr_t) tlist[i]->verts[0]];
		tlist[i]->verts[1] = vlist[(intptr_t) tlist[i]->verts[1]];
		tlist[i]->verts[2] = vlist[(intptr_t) tlist[i]->verts[2]];
	}

	/* get rid of triangles that use the same vertex more than once */
	for (i = ntris-1; i >= 0; i--) 
	{
		Triangle *tri = tlist[i];
		Vertex *v0 = tri->verts[0];
		Vertex *v1 = tri->verts[1];
		Vertex *v2 = tri->verts[2];

		if (v0 == v1 || v1 == v2 || v2 == v0) 
		{
			free (tlist[i]);
			ntris--;
			tlist[i] = tlist[ntris];
		}
	}
}

// empty constructor
Polyhedron::Polyhedron()
{
	nverts = nedges = ntris = 0;
	max_verts = max_tris = 50;

	vlist = new Vertex * [max_verts];
	tlist = new Triangle * [max_tris];
}



;//////////////////////////////////////////////////////////////////////////////
;// Utility Functions
;//////////////////////////////////////////////////////////////////////////////

// write out a polyhedron to a file
void Polyhedron::write_file(FILE *file)
{
  int i;
  PlyFile *ply;
  char **elist;
  int num_elem_types;

  /*** Write out the transformed PLY object ***/

  elist = get_element_list_ply (in_ply, &num_elem_types);
  ply = write_ply (file, num_elem_types, elist, in_ply->file_type);

  /* describe what properties go into the vertex elements */
  /* TODO: casting literals to char * to surpress warnings */

  describe_element_ply (ply, (char *)"vertex", nverts);
  describe_property_ply (ply, &vert_props[0]);
  describe_property_ply (ply, &vert_props[1]);
  describe_property_ply (ply, &vert_props[2]);
//  describe_other_properties_ply (ply, vert_other, offsetof(Vertex_io,other_props));

  describe_element_ply (ply, (char *)"face", ntris);
  describe_property_ply (ply, &face_props[0]);

//  describe_other_properties_ply (ply, face_other,
//                                offsetof(Face_io,other_props));

//  describe_other_elements_ply (ply, in_ply->other_elems);

  copy_comments_ply (ply, in_ply);
	char mm[1024];
	snprintf(mm, sizeof(mm), "modified by learnply");
//  append_comment_ply (ply, "modified by simvizply %f");
	  append_comment_ply (ply, mm);
  copy_obj_info_ply (ply, in_ply);

  header_complete_ply (ply);

  /* set up and write the vertex elements */
  put_element_setup_ply (ply, (char *)"vertex");
  for (i = 0; i < nverts; i++) {
    Vertex_io vert;

    /* copy info to the "vert" structure */
    vert.x = vlist[i]->x;
    vert.y = vlist[i]->y;
    vert.z = vlist[i]->z;
    vert.other_props = vlist[i]->other_props;

    put_element_ply (ply, (void *) &vert);
  }

  /* index all the vertices */
  for (i = 0; i < nverts; i++)
    vlist[i]->index = i;

  /* set up and write the face elements */
  put_element_setup_ply (ply, (char *)"face");

  Face_io face;
  face.verts = new int[3];
  
  for (i = 0; i < ntris; i++) {

    /* copy info to the "face" structure */
    face.nverts = 3;
    face.verts[0] = tlist[i]->verts[0]->index;
    face.verts[1] = tlist[i]->verts[1]->index;
    face.verts[2] = tlist[i]->verts[2]->index;
    face.other_props = tlist[i]->other_props;

    put_element_ply (ply, (void *) &face);
  }
  put_other_elements_ply (ply);

  close_ply (ply);
  free_ply (ply);
}

void Polyhedron::initialize(){
	icVector3 v1, v2;

	create_pointers();
	calc_edge_length();
	seed = -1;
}

void Polyhedron::finalize(){
	int i;

	for (i=0; i<ntris; i++){
		free(tlist[i]->other_props);
		free(tlist[i]);
	}
	for (i=0; i<nedges; i++) {
		free(elist[i]->tris);
		free(elist[i]);
	}
	for (i=0; i<nverts; i++) {
		free(vlist[i]->tris);
		free(vlist[i]->other_props);
		free(vlist[i]);
	}

	/* 2 */
	for (i = 0; i < ncorners; ++i)
		delete clist[i];
	delete[] clist;
	delete[] ctable;

	free(tlist);
	free(elist);
	free(vlist);
	if (!vert_other)
		free(vert_other);
	if (!face_other)
		free(face_other);

	silhouette.clear();
	major_hatches.clear();
	minor_hatches.clear();
}

/******************************************************************************
Find out if there is another face that shares an edge with a given face.

Entry:
  f1    - face that we're looking to share with
  v1,v2 - two vertices of f1 that define edge

Exit:
  return the matching face, or NULL if there is no such face
******************************************************************************/
Triangle *Polyhedron::find_common_edge(Triangle *f1, Vertex *v1, Vertex *v2)
{
  int i,j;
  Triangle *f2;
  Triangle *adjacent = NULL;

  /* look through all faces of the first vertex */

  for (i = 0; i < v1->ntris; i++) {
    f2 = v1->tris[i];
    if (f2 == f1)
      continue;
    /* examine the vertices of the face for a match with the second vertex */
    for (j = 0; j < f2->nverts; j++) {

      /* look for a match */
      if (f2->verts[j] == v2) {

#if 0
	/* watch out for triple edges */

        if (adjacent != NULL) {

	  fprintf (stderr, "model has triple edges\n");

	  fprintf (stderr, "face 1: ");
	  for (k = 0; k < f1->nverts; k++)
	    fprintf (stderr, "%d ", f1->iverts[k]);
	  fprintf (stderr, "\nface 2: ");
	  for (k = 0; k < f2->nverts; k++)
	    fprintf (stderr, "%d ", f2->iverts[k]);
	  fprintf (stderr, "\nface 3: ");
	  for (k = 0; k < adjacent->nverts; k++)
	    fprintf (stderr, "%d ", adjacent->iverts[k]);
	  fprintf (stderr, "\n");

	}

	/* if we've got a match, remember this face */
        adjacent = f2;
#endif

#if 1
	/* if we've got a match, return this face */
        return (f2);
#endif

      }
    }
  }

  return (adjacent);
}

/******************************************************************************
Create an edge.

Entry:
  v1,v2 - two vertices of f1 that define edge
******************************************************************************/
void Polyhedron::create_edge(Vertex *v1, Vertex *v2)
{
  int i,j;
  Triangle *f;

  /* make sure there is enough room for a new edge */

  if (nedges >= max_edges) {

    max_edges += 100;
    Edge **list = new Edge *[max_edges];

    /* copy the old list to the new one */
    for (i = 0; i < nedges; i++)
      list[i] = elist[i];

    /* replace list */
    free (elist);
    elist = list;
  }

  /* create the edge */

  elist[nedges] = new Edge;
  Edge *e = elist[nedges];
  e->index = nedges;
  e->verts[0] = v1;
  e->verts[1] = v2;
  nedges++;

  /* 2 */
  e->corners[0] = e->corners[1] = NULL;

  /* count all triangles that will share the edge, and do this */
  /* by looking through all faces of the first vertex */

  e->ntris = 0;

  for (i = 0; i < v1->ntris; i++) {
    f = v1->tris[i];
    /* examine the vertices of the face for a match with the second vertex */
    for (j = 0; j < 3; j++) {
      /* look for a match */
      if (f->verts[j] == v2) {
        e->ntris++;
        break;
      }
    }
  }

  /* make room for the face pointers (at least two) */
  if (e->ntris < 2)
    e->tris = new Triangle *[2];
  else
    e->tris = new Triangle *[e->ntris];

  /* create pointers from edges to faces and vice-versa */

  e->ntris = 0; /* start this out at zero again for creating ptrs to tris */

  for (i = 0; i < v1->ntris; i++) {

    f = v1->tris[i];

    /* examine the vertices of the face for a match with the second vertex */
    for (j = 0; j < 3; j++)
      if (f->verts[j] == v2) {

        e->tris[e->ntris] = f;
        e->ntris++;

        if (f->verts[(j+1)%3] == v1)
          f->edges[j] = e;
        else if (f->verts[(j+2)%3] == v1)
          f->edges[(j+2)%3] = e;
        else {
          fprintf (stderr, "Non-recoverable inconsistancy in create_edge()\n");
          exit (-1);
        }

        break;  /* we'll only find one instance of v2 */
      }

  }
}

/******************************************************************************
Create edges.
******************************************************************************/
void Polyhedron::create_edges()
{
  int i,j;
  Triangle *f;
  Vertex *v1,*v2;
  double count = 0;

  /* count up how many edges we may require */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < f->nverts; j++) {
      v1 = f->verts[j];
      v2 = f->verts[(j+1) % f->nverts];
      Triangle *result = find_common_edge (f, v1, v2);
      if (result)
        count += 0.5;
      else
        count += 1;
    }
  }

  /*
  printf ("counted %f edges\n", count);
  */

  /* create space for edge list */

  max_edges = (int) (count + 10);  /* leave some room for expansion */
  elist = new Edge *[max_edges];
  nedges = 0;

  /* zero out all the pointers from faces to edges */

  for (i = 0; i < ntris; i++)
    for (j = 0; j < 3; j++)
      tlist[i]->edges[j] = NULL;

  /* create all the edges by examining all the triangles */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < 3; j++) {
      /* skip over edges that we've already created */
      if (f->edges[j])
        continue;
      v1 = f->verts[j];
      v2 = f->verts[(j+1) % f->nverts];
      create_edge (v1, v2);
    }
  }
}

/******************************************************************************
Create pointers from vertices to faces.
******************************************************************************/
void Polyhedron::vertex_to_tri_ptrs()
{
  int i,j;
  Triangle *f;
  Vertex *v;

  /* zero the count of number of pointers to faces */

  for (i = 0; i < nverts; i++)
    vlist[i]->max_tris = 0;

  /* first just count all the face pointers needed for each vertex */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < f->nverts; j++)
      f->verts[j]->max_tris++;
  }

  /* allocate memory for face pointers of vertices */

  for (i = 0; i < nverts; i++) {
    vlist[i]->tris = (Triangle **)
		      malloc (sizeof (Triangle *) * vlist[i]->max_tris);
    vlist[i]->ntris = 0;
  }

  /* now actually create the face pointers */

  for (i = 0; i < ntris; i++) {
    f = tlist[i];
    for (j = 0; j < f->nverts; j++) {
      v = f->verts[j];
      v->tris[v->ntris] = f;
      v->ntris++;
    }
  }
}

/******************************************************************************
Find the other triangle that is incident on an edge, or NULL if there is
no other.
******************************************************************************/
Triangle *Polyhedron::other_triangle(Edge *edge, Triangle *tri)
{
  /* search for any other triangle */

  for (int i = 0; i < edge->ntris; i++)
    if (edge->tris[i] != tri)
      return (edge->tris[i]);

  /* there is no such other triangle if we get here */
  return (NULL);
}

/******************************************************************************
Order the pointers to faces that are around a given vertex.

Entry:
  v - vertex whose face list is to be ordered
******************************************************************************/
void Polyhedron::order_vertex_to_tri_ptrs(Vertex *v)
{
  int i,j;
  Triangle *f;
  Triangle *fnext;
  int nf;
  int vindex;
  int boundary;
  int count;

  nf = v->ntris;
  f = v->tris[0];

  /* go backwards (clockwise) around faces that surround a vertex */
  /* to find out if we reach a boundary */

  boundary = 0;

  for (i = 1; i <= nf; i++) {

    /* find reference to v in f */
    vindex = -1;
    for (j = 0; j < f->nverts; j++)
      if (f->verts[j] == v) {
	vindex = j;
	break;
      }

    /* error check */
    if (vindex == -1) {
      fprintf (stderr, "can't find vertex #1\n");
      exit (-1);
    }

    /* corresponding face is the previous one around v */
    fnext = other_triangle (f->edges[vindex], f);

    /* see if we've reached a boundary, and if so then place the */
    /* current face in the first position of the vertice's face list */

    if (fnext == NULL) {
      /* find reference to f in v */
      for (j = 0; j < v->ntris; j++)
        if (v->tris[j] == f) {
	  v->tris[j] = v->tris[0];
	  v->tris[0] = f;
	  break;
	}
      boundary = 1;
      break;
    }

    f = fnext;
  }

  /* now walk around the faces in the forward direction and place */
  /* them in order */

  f = v->tris[0];
  count = 0;

  for (i = 1; i < nf; i++) {

    /* find reference to vertex in f */
    vindex = -1;
    for (j = 0; j < f->nverts; j++)
      if (f->verts[(j+1) % f->nverts] == v) {
	vindex = j;
	break;
      }

    /* error check */
    if (vindex == -1) {
      fprintf (stderr, "can't find vertex #2\n");
      exit (-1);
    }

    /* corresponding face is next one around v */
    fnext = other_triangle (f->edges[vindex], f);

    /* break out of loop if we've reached a boundary */
    count = i;
    if (fnext == NULL) {
      break;
    }

    /* swap the next face into its proper place in the face list */
    for (j = 0; j < v->ntris; j++)
      if (v->tris[j] == fnext) {
	v->tris[j] = v->tris[i];
	v->tris[i] = fnext;
	break;
      }

    f = fnext;
  }
}

/******************************************************************************
Find the index to a given vertex in the list of vertices of a given face.

Entry:
  f - face whose vertex list is to be searched
  v - vertex to return reference to

Exit:
  returns index in face's list, or -1 if vertex not found
******************************************************************************/
int Polyhedron::face_to_vertex_ref(Triangle *f, Vertex *v)
{
  int j;
  int vindex = -1;

  for (j = 0; j < f->nverts; j++)
    if (f->verts[j] == v) {
      vindex = j;
      break;
    }

  return (vindex);
}

/******************************************************************************
Create various face and vertex pointers.
******************************************************************************/
void Polyhedron::create_pointers()
{
  int i;

  /* index the vertices and triangles */

  for (i = 0; i < nverts; i++)
    vlist[i]->index = i;

  for (i = 0; i < ntris; i++) 
    tlist[i]->index = i;

  /* create pointers from vertices to triangles */
  vertex_to_tri_ptrs();

  /* make edges */
  create_edges();


  /* order the pointers from vertices to faces */
	for (i = 0; i < nverts; i++){
//		if (i %1000 == 0)
//			fprintf(stderr, "ordering %d of %d vertices\n", i, nverts);
    order_vertex_to_tri_ptrs(vlist[i]);
		
	}
  /* index the edges */

  for (i = 0; i < nedges; i++){
//		if (i %1000 == 0)
//			fprintf(stderr, "indexing %d of %d edges\n", i, nedges);
    elist[i]->index = i;
	}

}

void Polyhedron::calc_bounding_sphere()
{
	unsigned int i;
	icVector3 min, max;

	for (i=0; i<nverts; i++)
	{
		if (i==0)  {
			min.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
			max.set(vlist[i]->x, vlist[i]->y, vlist[i]->z);
		}
		else 
		{
			if (vlist[i]->x < min.entry[0])
				min.entry[0] = vlist[i]->x;
			if (vlist[i]->x > max.entry[0])
				max.entry[0] = vlist[i]->x;
			if (vlist[i]->y < min.entry[1])
				min.entry[1] = vlist[i]->y;
			if (vlist[i]->y > max.entry[1])
				max.entry[1] = vlist[i]->y;
			if (vlist[i]->z < min.entry[2])
				min.entry[2] = vlist[i]->z;
			if (vlist[i]->z > max.entry[2])
				max.entry[2] = vlist[i]->z;
		}
	}
	center = (min + max) * 0.5;
	radius = length(center - min);
}

void Polyhedron::calc_edge_length()
{
	int i;
	icVector3 v1, v2;

	for (i=0; i<nedges; i++) {
		v1.set(elist[i]->verts[0]->x, elist[i]->verts[0]->y, elist[i]->verts[0]->z);
		v2.set(elist[i]->verts[1]->x, elist[i]->verts[1]->y, elist[i]->verts[1]->z);
		elist[i]->length = length(v1-v2);
	}
}

int Polyhedron::calc_face_normals_and_area()
{
	unsigned int i, j;
	icVector3 v0, v1, v2;
	Triangle *temp_t;
	double length[3];

	area = 0.0;
	for (i=0; i<ntris; i++){
		for (j=0; j<3; j++)
			length[j] = tlist[i]->edges[j]->length;
		double temp_s = (length[0] + length[1] + length[2])/2.0;
		tlist[i]->area = sqrt(temp_s*(temp_s-length[0])*(temp_s-length[1])*(temp_s-length[2]));

		area += tlist[i]->area;
		temp_t = tlist[i];
		v1.set(vlist[tlist[i]->verts[0]->index]->x, vlist[tlist[i]->verts[0]->index]->y, vlist[tlist[i]->verts[0]->index]->z);
		v2.set(vlist[tlist[i]->verts[1]->index]->x, vlist[tlist[i]->verts[1]->index]->y, vlist[tlist[i]->verts[1]->index]->z);
		v0.set(vlist[tlist[i]->verts[2]->index]->x, vlist[tlist[i]->verts[2]->index]->y, vlist[tlist[i]->verts[2]->index]->z);
		tlist[i]->normal = cross(v0-v1, v2-v1);
		normalize(tlist[i]->normal);
	}

	double signedvolume = 0.0;
	icVector3 test = center;
	for (i=0; i<ntris; i++){
		icVector3 cent(vlist[tlist[i]->verts[0]->index]->x, vlist[tlist[i]->verts[0]->index]->y, vlist[tlist[i]->verts[0]->index]->z);
		signedvolume += dot(test-cent, tlist[i]->normal)*tlist[i]->area;
	}
	signedvolume /= area;
	if (signedvolume < 0)
		return 0; // return orientation
	else 
	{
		for (i=0; i<ntris; i++) // flip normals
			tlist[i]->normal *= -1.0;
		return 1; // return orientation
	}
}



;//////////////////////////////////////////////////////////////////////////////
;// Project 1 functions
;//////////////////////////////////////////////////////////////////////////////

/* 3a */
void Polyhedron::compute_euler_characteristic()
{
	int x = nverts - nedges + ntris;
	printf("Euler characteristic: %d\n", x);
}

/* 3b */
static double compute_angle(Corner *c)
{
	icVector3 v1(c->v->x, c->v->y, c->v->z);
	icVector3 v2(c->n->v->x, c->n->v->y, c->n->v->z);
	icVector3 v3(c->p->v->x, c->p->v->y, c->p->v->z);

	icVector3 v21 = v2 - v1;
	normalize(v21);

	icVector3 v31 = v3 - v1;
	normalize(v31);

	return acos(dot(v21, v31));
}

/* 3b */
void Polyhedron::compute_gaussian_curvature_angle_deficit()
{
	K_min = INFINITY, K_max = -INFINITY;

	double K_disc = 0.0;
	for (int i = 0; i < nverts; ++i)
	{
		double ang_tot = 0.0;

		Vertex *v = vlist[i];
		for (int j = 0; j < v->ncorners; ++j)
			ang_tot += v->corners[j]->angle;

		v->K = 2.0 * M_PI - ang_tot;
		K_disc += v->K;

		K_min = std::min(K_min, v->K);
		K_max = std::max(K_max, v->K);
	}

	printf("Discrete Gaussian Curvature (total angle deficit): %g\n", K_disc);
	printf("Minimum/Maximum Vertex Curvature: %g/%g\n", K_min, K_max);
}

/* 3c */
void Polyhedron::compute_gaussian_curvature_valence_deficit()
{
	min_valence_deficit = 6;

	int total_deficit = 0;
	for (int i = 0; i < nverts; ++i)
	{
		vlist[i]->valence_deficit = 6 - vlist[i]->ntris;
		//vlist[i]->K = M_PI / 3.0 * vlist[i]->valence_deficit;
		//K_disc += vlist[i]->K;
		total_deficit += vlist[i]->valence_deficit;

		min_valence_deficit = std::min(min_valence_deficit, vlist[i]->valence_deficit);
	}

	printf("Total Valence Deficit: %d\n", total_deficit);
	printf("Minimum Valence Deficit: %d\n", min_valence_deficit);
}

/* 3d */
void Polyhedron::compute_handles()
{
	int x = nverts - nedges + ntris;
	int nhandles = (2 - x) / 2;
	printf("Handle count: %d\n", nhandles);
}

/* 2 */
void Polyhedron::create_corners()
{
	int ci;

	/* always 3 corners per triangle */
	max_corners = ncorners = 3 * ntris;
	clist = new Corner *[max_corners];
	for (int i = 0; i < ncorners; ++i)
	{
		clist[i] = new Corner;
		clist[i]->index = i;
	}

	/* setup corners */
	ci = 0;
	for (int i = 0; i < ntris; ++i)
	{
		Triangle *t = tlist[i];
		for (int j = 0; j < 3; ++j)
		{
			Corner *c = clist[ci + j];

			/* assign edge/vertex/triangle */
			c->v = t->verts[j];
			c->e = t->edges[(j+1)%3]; // see edge creation for why this works
			c->t = t;
			t->corners[j] = c;

			/* compute next/previous */
			c->n = clist[ci + (j + 1) % 3];
			c->p = clist[ci + (j + 2) % 3];

			/* store corner in vertex */
			c->v->ncorners++;
			Corner **corners = (Corner **)realloc(c->v->corners, sizeof(Corner *) * c->v->ncorners);
			if (!corners)
				abort();
			corners[c->v->ncorners - 1] = c;
			c->v->corners = corners;
		}

		ci += 3;
	}

	/* construct corner table */
	ctable = new CornerTableEntry[ncorners];
	for (int i = 0; i < ncorners; ++i)
	{
		ctable[i].c = i;
		ctable[i].v_min = std::min(clist[i]->p->v->index, clist[i]->n->v->index);
		ctable[i].v_max = std::max(clist[i]->p->v->index, clist[i]->n->v->index);
		ctable[i].o = 0;
	}

	/* sort corner table */
	std::sort(ctable, ctable + ncorners, [](CornerTableEntry &a, CornerTableEntry &b) {
		/* returns an entry is "less than" another */
		if (a.v_min != b.v_min)
			return a.v_min < b.v_min;
		if (a.v_max != b.v_max)
			return a.v_max < b.v_max;
		return a.c < b.c;
	});

	/* link opposite corners */
	for (int i = 0; i < ncorners; i += 2)
	{
		CornerTableEntry &c = ctable[i];
		CornerTableEntry &o = ctable[i+1];

		/* link in corner table */
		c.o = o.c;
		o.o = c.c;

		/* create link in corner structure */
		clist[c.c]->o = clist[c.o];
		clist[o.c]->o = clist[o.o];

		/* track each corner in its edge */
		assert(clist[c.c]->e == clist[o.c]->e);
		clist[c.c]->e->corners[0] = clist[c.c];
		clist[o.c]->e->corners[1] = clist[c.o];
	}

	/* calculate corner angles */
	for (int i = 0; i < ncorners; ++i)
	{
		Corner *c = clist[i];
		c->angle = compute_angle(c);
	}

	/* print corner table */
	/*for (int i = 0; i < ncorners; ++i)
	{
		CornerTableEntry &c = ctable[i];
		printf("%5d %5d %5d %5d\n", c.c, c.v_min, c.v_max, c.o);
	}*/
}

/* 1c */
void Polyhedron::average_normals()
{
	int i, j;

	for (i = 0; i < nverts; i++) {
		vlist[i]->normal = icVector3(0.0);
		for (j = 0; j < vlist[i]->ntris; j++)
			vlist[i]->normal += vlist[i]->tris[j]->normal;
		normalize(vlist[i]->normal);
	}
}


;//////////////////////////////////////////////////////////////////////////////
;// Project 2 functions
;//////////////////////////////////////////////////////////////////////////////

/* Problem 1 */

void Polyhedron::compute_silhouette_edges(const icMatrix3x3& view, const icVector3& translate)
{
	// clear silhouette vector
	silhouette.clear();
	std::vector<bool> face_forward(ntris, false);

	// tag all faces as facing forward or backward
	Triangle* tri;
	icVector3 eye, norm;
	double dprod;
	for (int i = 0; i < ntris; i++) {
		tri = tlist[i];
		icVector3 loc(tri->verts[0]->x + tri->verts[1]->x + tri->verts[2]->x,
			tri->verts[0]->y + tri->verts[1]->y + tri->verts[2]->y,
			tri->verts[0]->z + tri->verts[1]->z + tri->verts[2]->z);
		loc /= 3.0;
		eye = view * loc + translate;
		norm = view * tri->normal;

		dprod = dot(norm, eye);
		face_forward[i] = (dprod < 0.0);
	}

	// find edges where the faces switch direction
	Edge* edge;
	for (int i = 0; i < nedges; i++) {
		edge = elist[i];
		if (face_forward[edge->tris[0]->index] != face_forward[edge->tris[1]->index]) {
			LineSegment segment = LineSegment(edge->verts[0]->x, edge->verts[0]->y, edge->verts[0]->z,
				edge->verts[1]->x, edge->verts[1]->y, edge->verts[1]->z);
			silhouette.push_back(segment);
		}
	}

}

void Polyhedron::compute_silhouette_faces(const icMatrix3x3& view, const icVector3& translate)
{
	// clear silhouette edge
	silhouette.clear();
	std::vector<double> face_forward(ntris, 0.0);

	// tag all the vertices as facing forward or backward
	Vertex* vert;
	icVector3 eye, norm;
	double dprod;
	for (int i = 0; i < nverts; i++) {
		vert = vlist[i];
		icVector3 loc(vert->x, vert->y, vert->z);

		eye = view * loc + translate;
		norm = view * vert->normal;

		dprod = dot(norm, eye);
		face_forward[i] = dprod;
	}

	// loop through faces to look for silhouette segments going through them
	Triangle* tri;
	Vertex* v0, * v1;
	std::vector<icVector3> crossings;
	for (int i = 0; i < ntris; i++) {
		tri = tlist[i];
		crossings.clear();

		// find edges where vertex normals switch direction
		for (int j = 0; j < 3; j++) {
			v0 = tri->edges[j]->verts[0];
			v1 = tri->edges[j]->verts[1];
			if ((face_forward[v0->index] > 0.0 && face_forward[v1->index] <= 0.0) ||
				(face_forward[v0->index] <= 0.0 && face_forward[v1->index] > 0.0) ) {
				// find where the silhouette crosses the edge using linear interpolation
				double t = face_forward[v0->index] / (face_forward[v0->index] - face_forward[v1->index]);
				icVector3 loc0(v0->x, v0->y, v0->z);
				icVector3 loc1(v1->x, v1->y, v1->z);
				icVector3 crossing = loc0 + t * (loc1 - loc0);
				crossings.push_back(crossing);
			}
		}

		if (crossings.size() == 2) {
			LineSegment segment(crossings[0], crossings[1]);
			silhouette.push_back(segment);
		}
	}

}

;//////////////////////////////////////////////////////////////////////////////
;// Final Project functions
;//////////////////////////////////////////////////////////////////////////////


int Polyhedron::find_odd_visibility(const bool visibility[3]) {
	if (visibility[1] == visibility[2] && visibility[0] != visibility[1]) {
		return 0;
	}
	if (visibility[0] == visibility[2] && visibility[1] != visibility[0]) {
		return 1;
	}
	if (visibility[0] == visibility[1] && visibility[2] != visibility[0]) {
		return 2;
	}

	return -1;	// should not happen if the remesh_silhouette logic is correct
}

void Polyhedron::compute_tangents_for_an_edge(Edge* edge) {
	icVector3 v1(edge->verts[0]->x, edge->verts[0]->y, edge->verts[0]->z);
	icVector3 v2(edge->verts[1]->x, edge->verts[1]->y, edge->verts[1]->z);

	icVector3 v1_normal = edge->verts[0]->normal;
	icVector3 v2_normal = edge->verts[1]->normal;

	icVector3 e_vector = v2 - v1;
	double edge_length = edge->length;

	// Tangent direction T1 = (V2−V1)−[(V2−V1)⋅N1]⋅N1
	icVector3 tangent_v1 = e_vector - dot(e_vector, v1_normal) * v1_normal;
	normalize(tangent_v1);

	// Tangent direction T2 = (V2−V1)−[(V2−V1)⋅N2]⋅N2
	icVector3 tangent_v2 = e_vector - dot(e_vector, v2_normal) * v2_normal;
	normalize(tangent_v2);

	// Compute the cos_theta between each tangent direction and the vector (V2-V1)
	normalize(e_vector);
	double cos_theta_v1 = dot(e_vector, tangent_v1);
	double cos_theta_v2 = dot(e_vector, tangent_v2);

	// Calculate the lengths of two tangent vectors for later Hermite interpolation
	double tangent_v1_len = 2 * edge_length / (1 + cos_theta_v1);
	double tangent_v2_len = 2 * edge_length / (1 + cos_theta_v2);

	// Update two tangent vectors for the edge
	edge->tangent_v1 = tangent_v1 * tangent_v1_len;
	edge->tangent_v2 = tangent_v2 * tangent_v2_len;
}

icVector3 Polyhedron::hermite_interpolation(const icVector3& v1, const icVector3& v2,
											const icVector3& t1, const icVector3& t2,
											double u) {
	// Perform Hermite interpolation on an edge to find the interpolated point
	// S(u) = (2V1−2V2+T1+T2)*u^3 - (3V1−3V2+2T1+T2)*u^2 + T1*u + V1
	icVector3 A = (2.f*v1 - 2.f*v2 + t1 + t2) * (u * u * u);
	icVector3 B = (3.f*v1 - 3.f*v2 + 2.f*t1 + t2) * (u * u);
	icVector3 C = t1 * u;
	icVector3 D = v1;

	return A - B + C + D;
}

void Polyhedron::compute_silhouette_point(Edge* edge, const icVector3& camera_pos) {
	//		N1(u)=(1−u)⋅N1+u⋅N2
	//		D(u)=(1−u)⋅D1+u⋅D2
	// Solve the quadratic equation given by D(u)*N1(u)=0 to locate the silhouette point:
	//		[(1−u)⋅D1+u⋅D2] * [(1−u)⋅N1+u⋅N2] = 0
	//		(1-u)^2*D1*N1 + (1-u)*u*(D2*N1 + D1*N2) + u^2*D2*N2 = 0
	icVector3 v1(edge->verts[0]->x, edge->verts[0]->y, edge->verts[0]->z);
	icVector3 v2(edge->verts[1]->x, edge->verts[1]->y, edge->verts[1]->z);

	icVector3 N1 = edge->verts[0]->normal;
	icVector3 N2 = edge->verts[1]->normal;

	icVector3 D1 = camera_pos - v1;
	icVector3 D2 = camera_pos - v2;

	double cA = dot(D1, N1);
	double cB = dot(D2, N1);
	double cC = dot(D1, N2);
	double cD = dot(D2, N2);

	//		u^2(A-B-C+D) + u(B+C-2A) + A = 0
	double a = cA - cB - cC + cD;
	double b = cB + cC - 2 * cA;
	double c = cA;

	double disc = sqrt( b*b - 4.f*a*c );
	if (disc < 0) {
		printf("Error: discriminant < 0\n");
		return;
	}

	double u = (-b + disc) / (2.0 * a);
	// check if interpolant u is valid (within the range [0, 1])
	// If not, use the other solution
	if ( u < 0.0 || u > 1.0 ) {
		u = (-b - disc) / (2.0 * a);
	}
	// raise error if no valid u
	if ( u < 0.0 || u > 1.0 ) {
		printf("invalid interpolant u for silhouette points");
		return;
	}

	icVector3 interpolated_normal = (1.0 - u) * N1 + u * N2;
	normalize(interpolated_normal);
	edge->silhouette_point_normal = interpolated_normal;
	edge->silhouette_point = hermite_interpolation(v1, v2, edge->tangent_v1, edge->tangent_v2, u);
}

std::pair<Edge*, Edge*> Polyhedron::compute_silhouette_bridges(Triangle* tri, const bool visibility[3], Vertex* odd_vertex)
{
	std::pair<Edge*, Edge*> silhouette_bridges;

	// Find the two edges with visibility (+, -)
	int j = 0;
	Edge* curr_edge;
	for (int i = 0; i < 3; i++) {
		curr_edge = tri->edges[i];
		// If edge[i] connects to odd_vertex, it is an edge for silhouette bridge
		if ( curr_edge->verts[0] == odd_vertex || curr_edge->verts[1] == odd_vertex ) {
			if (j == 0) {
				silhouette_bridges.first = curr_edge;
				compute_tangents_for_an_edge(curr_edge);
			}
			if (j == 1) {
				silhouette_bridges.second = curr_edge;
				compute_tangents_for_an_edge(curr_edge);
			}
			j++;
		}
	}
	return silhouette_bridges;
}

icVector3 Polyhedron::compute_silhouette_segment(Edge* e1, Edge* e2, float u) {
	// Interpolate along the silhouette segment connecting two silhouette bridges
	// using identified silhouette points along with their normals
	icVector3 v1 = e1->silhouette_point;
	icVector3 v2 = e2->silhouette_point;

	icVector3 n1 = e1->silhouette_point_normal;
	icVector3 n2 = e2->silhouette_point_normal;

	// Compute the tangent vectors at each silhouette point
	icVector3 e_vector = v2 - v1;
	double edge_length = length(e_vector);

	// Tangent direction T1 = (V2−V1)−[(V2−V1)⋅N1]⋅N1
	icVector3 tangent_v1 = e_vector - dot(e_vector, n1) * n1;
	normalize(tangent_v1);

	// Tangent direction T2 = (V2−V1)−[(V2−V1)⋅N2]⋅N2
	icVector3 tangent_v2 = e_vector - dot(e_vector, n2) * n2;
	normalize(tangent_v2);

	// Compute the cos_theta between each tangent direction and the vector (V2-V1)
	normalize(e_vector);
	double cos_theta_v1 = dot(e_vector, tangent_v1);
	double cos_theta_v2 = dot(e_vector, tangent_v2);

	// Calculate the lengths of two tangent vectors for later Hermite interpolation
	double tangent_v1_len = 2 * edge_length / (1 + cos_theta_v1);
	double tangent_v2_len = 2 * edge_length / (1 + cos_theta_v2);

	// Update two tangent vectors for the edge
	tangent_v1 = tangent_v1 * tangent_v1_len;
	tangent_v2 = tangent_v2 * tangent_v2_len;

	return hermite_interpolation(v1, v2, tangent_v1, tangent_v2, u);

}

int Polyhedron::compute_segment_num(Edge* e1, Edge* e2) {
	// Compute the midpoint of the chord formed by two silhouette points
	icVector3 v1 = e1->silhouette_point;
	icVector3 v2 = e2->silhouette_point;
	icVector3 mid = (v1 + v2) * 0.5f;

	// Compute the midpoint of the expected hermite curve
	icVector3 hermite_mid = compute_silhouette_segment(e1, e2, 0.5f);

	// Estimate arc height h and arc curvature r
	double h = length(hermite_mid - mid);
	icVector3 n1 = e1->silhouette_point_normal;
	icVector3 n2 = e2->silhouette_point_normal;
	double cos_theta = dot(n1, n2);
	double r = sqrt((1.0 - cos_theta) / 2.f);

	int ns = std::ceil(0.5 * sqrt(h * r));
	ns = std::max(2, ns);

	return ns;
}

Vertex* Polyhedron::create_new_vertex(icVector3& pos, icVector3& normal) {
	Vertex* v = new Vertex(pos.x, pos.y, pos.z);
	v->normal = normal;
	v->index = (int)remesh_data.new_vertices.size();
	remesh_data.new_vertices.push_back(v);
	return v;
}

Edge* Polyhedron::create_new_edge(Vertex* v1, Vertex* v2) {
	Edge* e = new Edge();
	e->verts[0] = v1;
	e->verts[1] = v2;
	e->index = (int)remesh_data.new_edges.size();
	e->is_new = true;
	remesh_data.new_edges.push_back(e);
	return e;
}

void Polyhedron::create_new_triangle(Vertex* v1, Vertex* v2, Vertex* v3, Edge* e1, Edge* e2, Edge* e3) {
	Triangle* tri = new Triangle();
	tri->nverts = 3;
	tri->verts[0] = v1;
	tri->verts[1] = v2;
	tri->verts[2] = v3;
	tri->edges[0] = e1;
	tri->edges[1] = e2;
	tri->edges[2] = e3;
	tri->index = (int)remesh_data.new_triangles.size();
	remesh_data.new_triangles.push_back(tri);
}

void Polyhedron::remesh_silhouette(const icMatrix3x3& view, const icVector3& translate)
{
	icVector3 camera_pos = inverse(view) * (-translate);

	// Clear previously stored remesh data
	remesh_data.clear();
	// Copy original vlist and elist to the remesh_data
	Vertex* vert;
	Edge* edge;
	for (int i = 0; i < nverts; i++) {
		vert = vlist[i];
		remesh_data.new_vertices.push_back(vert);
	}
	// DEBUG: printf( "vert_size: %d\n", (int)remesh_data.new_vertices.size());
	for (int i = 0; i < nedges; i++) {
		edge = elist[i];
		remesh_data.new_edges.push_back(edge);
	}
	// DEBUG: printf( "edge_size: %d\n", (int)remesh_data.new_edges.size());
	// DEBUG: printf( "tri_size: %d\n", (int)ntris);

	// Cache for shared silhouette vertices
	std::map<std::pair<int, int>, Vertex*> silhouette_vertex_cache;
	std::map<std::pair<int, int>, Vertex*> base_vertex_cache;

	// Iterate through all triangle faces in the original mesh
	Triangle* tri;
	for (int i = 0; i < ntris; i++) {
		/*
		 * Step #1: Perform visibility test to determine if it is a silhouette triangle
		 *			by checking if all its vertices has the same sign for visibility
		 */
		tri = tlist[i];
		bool visibility[3];
		Vertex* vert;
		icVector3 eye, norm;
		for (int j = 0; j < 3; j++) {
			vert = tri->verts[j];
			icVector3 loc(vert->x, vert->y, vert->z);
			eye = camera_pos - loc;
			normalize(eye);
			norm = vert->normal;
			normalize(norm);
			// From article: the vertex V is visible if and only if the viewpoint E is
			// above the plane  F(X;V)=0
			visibility[j] = dot(norm, eye) >= 0.f;
		}

		// Remeshing is needed only if it is silhouette triangle
		if (visibility[0] == visibility[1] && visibility[1] == visibility[2]) {
			tri->index = (int)remesh_data.new_triangles.size();
			remesh_data.new_triangles.push_back(tri);
		}
		else {
			/*
			 * Step #2: Find the vertex of the triangle with a different visibility compared
			 *			to other two vertices
			 */
			int odd_index = find_odd_visibility(visibility);
			if (odd_index == -1) {
				printf("Error: odd_index = -1, error logic happened for remesh_silhouette!");
				return;
			}
			Vertex* odd_vertex = tri->verts[odd_index];

			/*
			 * Step #3: Compute the silhouette bridges by Hermite interpolation for the two
			 *			edges with visibility (+, -) of the silhouette edge
			 */
			auto silhouette_bridges = compute_silhouette_bridges(tri, visibility, odd_vertex);
			Edge* e1 = silhouette_bridges.first;
			Edge* e2 = silhouette_bridges.second;

			/*
			 * Step #4: Compute the silhouette point on the silhouette bridge
			 */
			compute_silhouette_point(e1, camera_pos);
			compute_silhouette_point(e2, camera_pos);

			/*
			 * Step #5: Adaptively compute the silhouette segments based on silhouette points along with their normals
			 */

			 const int num_samples = 8;
			 // Compute the number of sample points based on local info
			// const int num_samples = compute_segment_num(e1, e2);

			const float step = 1.f / num_samples;
			// Vertex* center_v = tri->verts[odd_index];
			Vertex* vA = (e1->verts[0] == odd_vertex) ? e1->verts[1] : e1->verts[0];
			icVector3 vA_pos = (vA->x, vA->y, vA->z);
			icVector3 vA_normal = vA->normal;
			Vertex* vB = (e2->verts[0] == odd_vertex) ? e2->verts[1] : e2->verts[0];
			icVector3 vB_pos = (vB->x, vB->y, vB->z);
			icVector3 vB_normal = vA->normal;
			// Vertex* vA = tri->verts[(odd_index + 1) % 3];
			// Vertex* vB = tri->verts[(odd_index + 2) % 3];

			// Determine remeshing strategy: compare depth of silhouette segment vs edge (vA-vB)
			icVector3 p_cam = inverse(view) * (e1->silhouette_point + e2->silhouette_point) * 0.5f;
			icVector3 e_cam = inverse(view) * ((icVector3(vA->x, vA->y, vA->z) + icVector3(vB->x, vB->y, vB->z)) * 0.5f);

			bool segment_behind_edge = p_cam.z > e_cam.z;

			for (int k = 0; k < num_samples; k++) {
				float t0 = step * k;
				float t1 = step * (k + 1);

				icVector3 p0 = compute_silhouette_segment(e1, e2, t0);
				icVector3 p1 = compute_silhouette_segment(e1, e2, t1);

				icVector3 n1 = e1->silhouette_point_normal;
				icVector3 n2 = e2->silhouette_point_normal;
				icVector3 n0_interp = (1.0 - t0) * n1 + t0 * n2;
				icVector3 n1_interp = (1.0 - t1) * n1 + t1 * n2;
				normalize(n0_interp);
				normalize(n1_interp);

				std::pair<int, int> key0(i, k);
				std::pair<int, int> key1(i, k + 1);
				Vertex* v0;
				Vertex* v1;

				if (silhouette_vertex_cache.count(key0)) v0 = silhouette_vertex_cache[key0];
				else {
					v0 = new Vertex(p0.x, p0.y, p0.z);
					v0->normal = n0_interp;
					v0->index = (int)remesh_data.new_vertices.size();
					remesh_data.new_vertices.push_back(v0);
					silhouette_vertex_cache[key0] = v0;
				}

				if (silhouette_vertex_cache.count(key1)) v1 = silhouette_vertex_cache[key1];
				else {
					v1 = new Vertex(p1.x, p1.y, p1.z);
					v1->normal = n1_interp;
					v1->index = (int)remesh_data.new_vertices.size();
					remesh_data.new_vertices.push_back(v1);
					silhouette_vertex_cache[key1] = v1;
				}

				auto make_edge = [&](Vertex* a, Vertex* b) {
					Edge* e = new Edge();
					e->verts[0] = a;
					e->verts[1] = b;
					e->index = (int)remesh_data.new_edges.size();
					e->is_new = true;
					remesh_data.new_edges.push_back(e);
					return e;
					};
				
				tri->index = (int)remesh_data.new_triangles.size();
				remesh_data.new_triangles.push_back(tri);

				if (segment_behind_edge) {

					icVector3 b_p0 = (1.0 - t0) * vA_pos + t0 * vB_pos;
					icVector3 b_p1 = (1.0 - t1) * vA_pos + t1 * vB_pos;

					icVector3 b_n0_interp = (1.0 - t0) * vA_normal + t0 * vB_normal;
					icVector3 b_n1_interp = (1.0 - t1) * vA_normal + t1 * vB_normal;
					normalize(b_n0_interp);
					normalize(b_n1_interp);

					Vertex* b_v0;
					if (k == 0) {
						b_v0 = vA;
					}
					else {
						if (base_vertex_cache.count(key0)) b_v0 = base_vertex_cache[key0];
						else {
							b_v0 = new Vertex(p0.x, p0.y, p0.z);
							b_v0->normal = b_n0_interp;
							b_v0->index = (int)remesh_data.new_vertices.size();
							remesh_data.new_vertices.push_back(b_v0);
							base_vertex_cache[key0] = b_v0;
						}
					}
					Vertex* b_v1;
					if (k == num_samples - 1) {
						b_v1 = vB;
					}
					else {
						if (base_vertex_cache.count(key1)) b_v1 = base_vertex_cache[key1];
						else {
							b_v1 = new Vertex(p1.x, p1.y, p1.z);
							b_v1->normal = b_n1_interp;
							b_v1->index = (int)remesh_data.new_vertices.size();
							remesh_data.new_vertices.push_back(b_v1);
							base_vertex_cache[key1] = b_v1;
						}
					}
					
					// form the triangulation between each silhouette segment and the base segment
					Edge* sil_base_e0 = make_edge(b_v0, v0);
					Edge* sil_base_e1 = make_edge(v0, b_v1);
					Edge* base_e = make_edge(b_v1, b_v0);
					create_new_triangle(b_v0, v0, b_v1, sil_base_e0, sil_base_e1, base_e);

					Edge* sil_base_e2 = make_edge(b_v1, v1);
					Edge* sil_e = make_edge(v1, v0);
					create_new_triangle(b_v1, v1, v0, sil_base_e2, sil_e, sil_base_e1);

					// ensure local convexity on the current triangle
					if (k == 0) {
						Edge* tri_e0 = make_edge(odd_vertex, b_v0);
						Edge* tri_e1 = make_edge(v0, odd_vertex);
						create_new_triangle(odd_vertex, b_v0, v0, tri_e0, sil_base_e0, tri_e1);
					}
					Edge* tri_e0 = make_edge(v0, odd_vertex);
					Edge* tri_e1 = make_edge(odd_vertex, v1);
					create_new_triangle(v0, odd_vertex, v1, tri_e0, tri_e1, sil_e);
					if (k == num_samples - 1) {
						Edge* tri_e0 = make_edge(v1, odd_vertex);
						Edge* tri_e1 = make_edge(odd_vertex, b_v1);
						// Edge* sil_base_e2 = make_edge(b_v1, v1);
						create_new_triangle(v1, odd_vertex, b_v1, tri_e0, tri_e1, sil_base_e2);
					}


				}
				else {
					if (k == 0) {
						// form the triangle between vA and the silhouette point
						Edge* init_tri_e0 = make_edge(odd_vertex, vA);
						Edge* init_tri_e1 = make_edge(vA, v0);
						Edge* init_tri_e2 = make_edge(v0, odd_vertex);

						create_new_triangle(odd_vertex, vA, v0, init_tri_e0, init_tri_e1, init_tri_e2);
					}
					// The sample points on the segment S(t) are connected to the odd vertex
					Edge* sil_tri_e0 = make_edge(odd_vertex, v0);
					Edge* sil_tri_e1 = make_edge(v0, v1);
					Edge* sil_tri_e2 = make_edge(v1, odd_vertex);

					create_new_triangle(odd_vertex, v0, v1, sil_tri_e0, sil_tri_e1, sil_tri_e2);

					if (k == num_samples - 1) {
						// form the triangle between vB and the silhouette point
						Edge* end_tri_e0 = make_edge(odd_vertex, vB);
						Edge* end_tri_e1 = make_edge(vB, v1);
						Edge* end_tri_e2 = make_edge(v1, odd_vertex);

						create_new_triangle(odd_vertex, vB, v1, end_tri_e0, end_tri_e1, end_tri_e2);
					} 
				}
			}
		}
		// DEBUG: printf( "remesh_vert_size: %d\n", (int)remesh_data.new_vertices.size());
		// DEBUG: printf( "remesh_tri_size: %d\n", (int)remesh_data.new_triangles.size());
	}
}