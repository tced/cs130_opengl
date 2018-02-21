/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

//*****GLOBAL VARIABLES AND FUNCTIONS****//////// 
vec4 curr_color;
MGLpoly_mode curr_type;
MGLmatrix_mode curr_matrix;  
mat4 projection_matrix;
mat4 modelview_matrix;   
//creating a vertex list 
struct Vertex {
   vec4 position; 
   vec4 color;
};

struct Triangle {
   Vertex A;
   Vertex B; 
   Vertex C; 
};  
vector<Vertex> list_vertices; 
vector<Triangle>list_triangles; 
/**
 * Standard macro to report errors
 */
inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

//*******HELPER FUNCTIONS******// 
mat4& current_matrix() { 
   if(curr_matrix == MGL_PROJECTION) {
       return projection_matrix; 
   }
   else {
     return modelview_matrix; 
   }
}
 
//CREATED: helper function for rasterize_triangle functn
MGLfloat GetArea(vec2 a, vec2 b, vec2 c) { 
   return (a[0] * (b[1]-c[1]) + a[1] * (c[0]-b[0]) + (b[0]*c[1]-b[1]*c[0])); 
}

 
void Rasterize_Triangle(const Triangle& tri, int width, int height, MGLpixel* data) {
    vec2 a, b, c; 
    float area, alpha, beta, gamma; 
    a[0] = (tri.A.position[0] + 1)*width/2-0.5; 
    a[1] = (tri.A.position[1] + 1)*height/2-0.5; 

    b[0] = (tri.B.position[0] + 1)*width/2-0.5; 
    b[1] = (tri.B.position[1] + 1)*height/2-0.5; 

    c[0] = (tri.C.position[0] + 1)*width/2-0.5; 
    c[1] = (tri.C.position[1] + 1)*height/2-0.5; 

    area = GetArea(a,b,c); 
 
    for(int i = 0; i < width; ++i) { 
       for (int j = 0; j < height; ++j) {
           vec2 p; 
           p[0] = i; 
           p[1] = j; 
   
           alpha = GetArea(p,b,c)/area; 
           beta = GetArea(a,p,c)/area; 
           gamma = GetArea(a,b,p)/area; 

           if (alpha >= 0 && beta >=0 && gamma >= 0) {
              data[i+j*width] = Make_Pixel(tri.A.color[0]*255, tri.A.color[1]*255,tri.A.color[2]*255); 
           }
       }
    }
 
}
//******END OF HELPER FUNCTIONS****//// 

/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
   //fill whole pixel data black  
   Make_Pixel(0,0,0); 

   //call rasterize_Triangle to raster each triangle on data then clear the triangle list 
   for (unsigned int i = 0; i < list_triangles.size(); ++i) {
       Rasterize_Triangle(list_triangles[i], width, height, data); 
   }
   list_triangles.clear(); 
}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
    curr_type = mode; 
}


/**
 * /Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
  //create a triangle using each triplet of vertices 
  if (curr_type == MGL_TRIANGLES) {
     for (unsigned int i = 0; i < list_vertices.size(); i=+3) {
          Triangle t;
          t.A = list_vertices.at(i);
          t.B = list_vertices.at(i+1); 
          t.C = list_vertices.at(i+2);  
          list_triangles.push_back(t);

     } 
  }
  //create a quad, using 3 triangles 
  if (curr_type == MGL_QUADS) {
     for (unsigned int j = 0; j < list_vertices.size(); j=+4) {
         unsigned int TopRight, TopLeft, BottomRight, BottomLeft;
         Triangle first, second; 
 
         TopRight = j + 2; 
         TopLeft = j + 3; 
         BottomLeft = j; 
         BottomRight = j + 1; 

         first.A = list_vertices.at(BottomLeft);
         first.B = list_vertices.at(BottomRight); 
         first.C = list_vertices.at(TopRight);

	 second.A = list_vertices.at(BottomLeft); 
         second.B = list_vertices.at(TopLeft); 
         second.C = list_vertices.at(TopRight);  
         
         list_triangles.push_back(first); 
         list_triangles.push_back(second);  

      } 
  }  
  list_vertices.clear(); 
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
   
   mglVertex3(x, y, 0);   
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{
   Vertex position; 
   vec4 position_4d; 
   position_4d[0] = x; 
   position_4d[1] = y; 
   position_4d[2] = z; 
   position_4d[3] = 1;
   position_4d = projection_matrix * modelview_matrix * position_4d; 
   position.position = position_4d; 
   position.color = curr_color;  
   
   list_vertices.push_back(position);   
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
    curr_matrix = mode; 
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
    
   mat4 identity_matrix = {{1,0,0,0,
                     0,1,0,0,
                     0,0,1,0,
                     0,0,0,1}};

   mat4& curr_mat = current_matrix(); 
   curr_mat = identity_matrix;
     
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
   float t_x, t_y, t_z; 
   t_x = -(right+left)/(right-left);
   t_y = -(top+bottom)/(top-bottom); 
   t_z = -(far+near)/(far-near); 

   mat4& curr_matrix = current_matrix(); 

   mat4 temp_matrix = {{2/(right-left),0,0,0, 
                        0,2/(top-bottom),0,0,
                        0,0,-2/(far-near),0, 
                        t_x, t_y, t_z, 1}}; 

   curr_matrix = temp_matrix * curr_matrix; 

}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
   float t_x, t_y, t_z; 
   t_x = -(right+left)/(right-left);
   t_y = -(top+bottom)/(top-bottom); 
   t_z = -(far+near)/(far-near); 

   mat4& curr_matrix = current_matrix(); 

   mat4 temp_matrix = {{2/(right-left),0,0,0,
                 0,2/(top-bottom),0,0,
                 0,0,-2/(far-near),0,
                 t_x, t_y, t_z, 1}}; 
  
   curr_matrix = temp_matrix * curr_matrix; 
}


/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
  curr_color[0] = red; 
  curr_color[1] = green; 
  curr_color[2] = blue;
}
