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
vec3 curr_color;
MGLpoly_mode curr_type;
MGLmatrix_mode mode_matrix;  
mat4 projection_matrix;
mat4 modelview_matrix;   
std::vector <mat4> projection_stack(1);
std::vector <mat4> modelview_stack(1); 
//creating a vertex list 
struct Vertex {
   vec4 position; 
   vec3 color;
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
   if(mode_matrix == MGL_PROJECTION) {
       return projection_matrix; 
   }
   else {
     return modelview_matrix; 
   } 
} 

mat4& top_of_active_matrix_stack(){
   if(mode_matrix == MGL_PROJECTION) {
     return projection_stack.back(); 
   }
   else if(mode_matrix == MGL_MODELVIEW) {
     return modelview_stack.back();
   }
   else {} 

}
 
 
//CREATED: helper function for rasterize_triangle functn
MGLfloat GetArea(vec2 a, vec2 b, vec2 c) { 
   return (a[0] * (b[1]-c[1]) + a[1] * (c[0]-b[0]) + (b[0]*c[1]-b[1]*c[0])); 
}

 
void Rasterize_Triangle(const Triangle& tri, int width, int height, MGLpixel* data) {
    //vec2 a, b, c; 
    MGLfloat i1, j1;  
    MGLfloat i2, j2; 
    MGLfloat i3, j3; 
    float area, alpha, beta, gamma; 
    
    i1 = ((tri.A.position[0]/tri.A.position[3] + 1)*width)/2-0.5; 
    j1 = ((tri.A.position[1]/tri.A.position[3] + 1)*height)/2-0.5; 
    vec2 a(i1, j1); 

    i2 = ((tri.B.position[0]/tri.B.position[3] + 1)*width)/2-0.5; 
    j2 = ((tri.B.position[1]/tri.B.position[3] + 1)*height)/2-0.5; 
    vec2 b(i2, j2); 
 
    i3 = ((tri.C.position[0]/tri.C.position[3] + 1)*width)/2-0.5; 
    j3 = ((tri.C.position[1]/tri.C.position[3] + 1)*height)/2-0.5; 
    vec2 c(i3, j3);     
    
    area = GetArea(a,b,c); 
 
    for(int i = 0; i < width; ++i) { 
       for (int j = 0; j < height; ++j) {
           vec2 p; 
           p[0] = i; 
           p[1] = j; 
   
           alpha = GetArea(p,b,c)/area; 
           beta = GetArea(a,p,c)/area; 
           gamma = GetArea(a,b,p)/area; 

           if ((alpha >= 0 && beta >=0 && gamma >= 0) 
              && (gamma >= 0 && gamma <= 1)){
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
   for(unsigned int i = 0; i< width; ++i) {
     for(unsigned int j = 0; j < height; ++j) {
        data[i+j*width] = Make_Pixel(0,0,0); 
     }
   }
   //Make_Pixel(0,0,0); 

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
  else if (curr_type == MGL_QUADS) {
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
   Vertex vertex;
   vertex.color = curr_color; 
   vertex.position[0] = x; 
   vertex.position[1] = y; 
   vertex.position[2] = z; 
   vertex.position[3] = 1.0f;
   vertex.position = projection_stack.back() * modelview_stack.back() * vertex.position; 
   list_vertices.push_back(vertex);

}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
    mode_matrix = mode;
    cout << "this is the matrix mode: " << mode_matrix;  
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
   mat4& pos = top_of_active_matrix_stack();   
   if (mode_matrix == MGL_PROJECTION) {
      //if(!projection_stack.empty()) {
        projection_stack.push_back(pos); 
      //}
   }  

   else {
      //if(!modelview_stack.empty()) { 
      	modelview_stack.push_back(pos);
      //} 
   }

 
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
    
   if (mode_matrix == MGL_PROJECTION) {
      if(!projection_stack.empty()) { 
         projection_stack.pop_back(); 
      }
   }
   else {
      if(!modelview_stack.empty()){  
       modelview_stack.pop_back(); 
      }
    }
    
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

  mat4& curr_stack = top_of_active_matrix_stack(); 
  curr_stack = identity_matrix; 
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
   if (mode_matrix == MGL_PROJECTION) {
      for (uint i = 0; i < 16; ++i) {
          projection_matrix.values[i] = matrix[i]; 
      }
   }
  
   if (mode_matrix == MGL_MODELVIEW) {
      for (uint i = 0; i < 16; ++i) {
          modelview_matrix.values[i] = matrix[i]; 
      }
   }
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
    mat4 mult_matrix; 
    mult_matrix.make_zero(); 
    mat4& curr_matrix = top_of_active_matrix_stack(); 
    for(uint i = 0; i < 16; ++i) {
       mult_matrix.values[i] = matrix[i]; 
    }
    curr_matrix = mult_matrix * curr_matrix; 
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
   mat4 translate_matrix = {{1,0,0,0,
  			     0,1,0,0,
                             0,0,1,0,
  			     x,y,z,1}}; 
   mat4& curr_stack = top_of_active_matrix_stack();
   curr_stack = curr_stack * translate_matrix; 
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
    mat4& curr_stack = top_of_active_matrix_stack(); 
    MGLfloat c = cos(angle * M_PI/180); 
    MGLfloat s = sin(angle * M_PI/180);
    MGLfloat normal = sqrt(x*x + y*y + z*z); 
    x = x/normal; 
    y = y/normal;
    z = z/normal;  

    mat4 rotate_matrix = {{x*x*(1-c)+c,y*x*(1-c)+z*s, x*z*(1-c)-y*s, 0,
			   x*y*(1-c)-z*s,y*y*(1-c)+c, y*z*(1-c)+x*s, 0,
   			   x*z*(1-c)+y*s, y*z*(1-c)-x*s, z*z*(1-c)+c, 0,
			   0,0,0,1}};
    curr_stack = curr_stack * rotate_matrix;  
}

/**
  //matrix_stack = identity_matrix; 
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{  
   mat4 scale_matrix = {{ x, 0, 0, 0, 
			  0, y, 0, 0, 
			  0, 0, z, 0,
		          0, 0, 0, 1}}; 
  mat4& curr_matrix = top_of_active_matrix_stack(); 
  curr_matrix = curr_matrix * scale_matrix;  
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
   MGLfloat A, B, C, D; 
   A = (right+left)/(right-left);
   B = (top+bottom)/(top-bottom); 
   C = -(far+near)/(far-near); 
   D = -2*(far*near)/(far-near); 

   mat4 frustum_matrix = {{2*near/(right-left),0,0,0, 
                        0,2*near/(top-bottom),0,0,
                        A,B,C,-1, 
                        0, 0, D, 0}}; 
   //mat4& matrix = current_matrix(); 
   //matrix = frustum_matrix * matrix; 
   mat4& curr_matrix = top_of_active_matrix_stack(); 
   curr_matrix = curr_matrix * frustum_matrix;
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
   MGLfloat t_x, t_y, t_z; 
   t_x = -(right+left)/(right-left);
   t_y = -(top+bottom)/(top-bottom); 
   t_z = -(far+near)/(far-near); 


   mat4 ortho_matrix = {{2/(right-left),0,0,0,
                 0,2/(top-bottom),0,0,
                 0,0,-2/(far-near),0,
                 t_x, t_y, t_z, 1}}; 
  
   //mat4& matrix = current_matrix(); 
   //matrix = ortho_matrix * matrix; 
   mat4& curr_matrix = top_of_active_matrix_stack(); 
   curr_matrix = curr_matrix * ortho_matrix; 
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
