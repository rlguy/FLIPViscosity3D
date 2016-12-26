#ifndef ARRAY3_UTILS_H
#define ARRAY3_UTILS_H

#include "array3.h"
#include "util.h"

template<class T>
T interpolate_value(vmath::vec3 point, const Array3<T, Array1<T> >& grid) {
   int i,j,k;
   float fi,fj,fk;

   get_barycentric(point[0], i, fi, 0, grid.ni);
   get_barycentric(point[1], j, fj, 0, grid.nj);
   get_barycentric(point[2], k, fk, 0, grid.nk);

   return trilerp(
         grid(i,j,k), grid(i+1,j,k), grid(i,j+1,k), grid(i+1,j+1,k), 
         grid(i,j,k+1), grid(i+1,j,k+1), grid(i,j+1,k+1), grid(i+1,j+1,k+1), 
         fi,fj,fk);
}

float interpolate_value(vmath::vec3 point, Array3d<float> &grid) {
   int i,j,k;
   float fi,fj,fk;

   get_barycentric(point[0], i, fi, 0, grid.width);
   get_barycentric(point[1], j, fj, 0, grid.height);
   get_barycentric(point[2], k, fk, 0, grid.depth);

   return trilerp(
         grid(i,j,k), grid(i+1,j,k), grid(i,j+1,k), grid(i+1,j+1,k), 
         grid(i,j,k+1), grid(i+1,j,k+1), grid(i,j+1,k+1), grid(i+1,j+1,k+1), 
         fi,fj,fk);
}


template<class T>
T interpolate_gradient(vmath::vec3 &gradient, vmath::vec3 point, const Array3<T, Array1<T> >& grid) {
   int i,j,k;
   float fx,fy,fz;
   
   get_barycentric(point[0], i, fx, 0, grid.ni);
   get_barycentric(point[1], j, fy, 0, grid.nj);
   get_barycentric(point[2], k, fz, 0, grid.nk);
   
   T v000 = grid(i,j,k);
   T v001 = grid(i,j,k+1);
   T v010 = grid(i,j+1,k);
   T v011 = grid(i,j+1,k+1);
   T v100 = grid(i+1,j,k);
   T v101 = grid(i+1,j,k+1);
   T v110 = grid(i+1,j+1,k);
   T v111 = grid(i+1,j+1,k+1);

   T ddx00 = (v100 - v000);
   T ddx10 = (v110 - v010);
   T ddx01 = (v101 - v001);
   T ddx11 = (v111 - v011);
   T dv_dx = bilerp(ddx00,ddx10,ddx01,ddx11, fy,fz);

   T ddy00 = (v010 - v000);
   T ddy10 = (v110 - v100);
   T ddy01 = (v011 - v001);
   T ddy11 = (v111 - v101);
   T dv_dy = bilerp(ddy00,ddy10,ddy01,ddy11, fx,fz);

   T ddz00 = (v001 - v000);
   T ddz10 = (v101 - v100);
   T ddz01 = (v011 - v010);
   T ddz11 = (v111 - v110);
   T dv_dz = bilerp(ddz00,ddz10,ddz01,ddz11, fx,fy);

   gradient.x = dv_dx;
   gradient.y = dv_dy;
   gradient.z = dv_dz;
   
   //return value for good measure.
   return trilerp(
      v000, v100,
      v010, v110, 
      v001, v101,
      v011, v111,
      fx, fy, fz);
}

float interpolate_gradient(vmath::vec3 &gradient, vmath::vec3 point, Array3d<float> &grid) {
   int i,j,k;
   float fx,fy,fz;
   
   get_barycentric(point[0], i, fx, 0, grid.width);
   get_barycentric(point[1], j, fy, 0, grid.height);
   get_barycentric(point[2], k, fz, 0, grid.depth);
   
   float v000 = grid(i,j,k);
   float v001 = grid(i,j,k+1);
   float v010 = grid(i,j+1,k);
   float v011 = grid(i,j+1,k+1);
   float v100 = grid(i+1,j,k);
   float v101 = grid(i+1,j,k+1);
   float v110 = grid(i+1,j+1,k);
   float v111 = grid(i+1,j+1,k+1);

   float ddx00 = (v100 - v000);
   float ddx10 = (v110 - v010);
   float ddx01 = (v101 - v001);
   float ddx11 = (v111 - v011);
   float dv_dx = bilerp(ddx00,ddx10,ddx01,ddx11, fy,fz);

   float ddy00 = (v010 - v000);
   float ddy10 = (v110 - v100);
   float ddy01 = (v011 - v001);
   float ddy11 = (v111 - v101);
   float dv_dy = bilerp(ddy00,ddy10,ddy01,ddy11, fx,fz);

   float ddz00 = (v001 - v000);
   float ddz10 = (v101 - v100);
   float ddz01 = (v011 - v010);
   float ddz11 = (v111 - v110);
   float dv_dz = bilerp(ddz00,ddz10,ddz01,ddz11, fx,fy);

   gradient.x = dv_dx;
   gradient.y = dv_dy;
   gradient.z = dv_dz;
   
   //return value for good measure.
   return trilerp(
      v000, v100,
      v010, v110, 
      v001, v101,
      v011, v111,
      fx, fy, fz);
}

#endif