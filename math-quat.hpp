//LQ
#include <stdio.h>
#include <cstdlib>
#include <memory>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <vector>

#define THIRD 0.333333333333333
#define ROOTTHREE 1.73205080756888       
#define PI 3.141592653589793238462643383279502884197

double* quatrot (double q[], double p[])
{
  //Rotate a point (double[3]) w.r.t. a quaternion (double[4])

  double* rotate = new double[3];
  for (int i=0; i<3; i++)
    {
      rotate[i] = 0.0;
    }

  double matrix[3][3];

  matrix[0][0] = 1.0 - 2*q[2]*q[2]- 2*q[3]*q[3];
  matrix[0][1] = 2*q[1]*q[2] - 2*q[0]*q[3];
  matrix[0][2] = 2*q[1]*q[3] + 2*q[0]*q[2];

  matrix[1][0] = 2*q[1]*q[2] + 2*q[0]*q[3];
  matrix[1][1] = 1.0 - 2*q[1]*q[1]- 2*q[3]*q[3];
  matrix[1][2] = 2*q[2]*q[3] - 2*q[0]*q[1];


  matrix[2][0] = 2*q[1]*q[3] - 2*q[0]*q[2];
  matrix[2][1] = 2*q[2]*q[3] + 2*q[0]*q[1];
  matrix[2][2] = 1.0 - 2*q[1]*q[1]- 2*q[2]*q[2];

 
  rotate[0] = matrix[0][0]*p[0] + matrix[0][1]*p[1] + matrix[0][2]*p[2];
  rotate[1] = matrix[1][0]*p[0] + matrix[1][1]*p[1] + matrix[1][2]*p[2];
  rotate[2] = matrix[2][0]*p[0] + matrix[2][1]*p[1] + matrix[2][2]*p[2];

  return rotate;
}




double* invquatrot (double q[], double p[])
{
  //Inverse rotate a point (double[3]) w.r.t. a quaternion (double[4])

  double* rotate = new double[3];
  for (int i=0; i<3; i++)
    {
      rotate[i] = 0.0;
    }

  double matrix[3][3];

  matrix[0][0] = 1.0 - 2*q[2]*q[2]- 2*q[3]*q[3];
  matrix[0][1] = 2*q[1]*q[2] + 2*q[0]*q[3];
  matrix[0][2] = 2*q[1]*q[3] - 2*q[0]*q[2];

  matrix[1][0] = 2*q[1]*q[2] - 2*q[0]*q[3];
  matrix[1][1] = 1.0 - 2*q[1]*q[1]- 2*q[3]*q[3];
  matrix[1][2] = 2*q[2]*q[3] + 2*q[0]*q[1];


  matrix[2][0] = 2*q[1]*q[3] + 2*q[0]*q[2];
  matrix[2][1] = 2*q[2]*q[3] - 2*q[0]*q[1];
  matrix[2][2] = 1.0 - 2*q[1]*q[1]- 2*q[2]*q[2];

 
  rotate[0] = matrix[0][0]*p[0] + matrix[0][1]*p[1] + matrix[0][2]*p[2];
  rotate[1] = matrix[1][0]*p[0] + matrix[1][1]*p[1] + matrix[1][2]*p[2];
  rotate[2] = matrix[2][0]*p[0] + matrix[2][1]*p[1] + matrix[2][2]*p[2];

  //  std::cout << "[invquatrot] " << rotate[0] << " " << rotate[1] << " " << rotate[2] << std::endl;

  return rotate;
}



double* quatmul (double q1[], double q2[])
{
  //Multiply two quaternions

  double* qprod = new double[4];
  for (int i=0; i<4; i++)
    {
      qprod[i] = 0.0;
    }

  qprod[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
  qprod[1] = q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2];
  qprod[2] = q1[0]*q2[2] - q1[1]*q2[3] + q1[2]*q2[0] + q1[3]*q2[1];
  qprod[3] = q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1] + q1[3]*q2[0];

  return qprod;
}



double* invquat (double q[])
{
  //Compute inverse of a quaternion

  double* qinv = new double[4];
  for (int i=0; i<4; i++)
    {
      qinv[i] = 0.0;
    }

  double qnorm = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];

  qinv[0] =  q[0]/qnorm;
  qinv[1] = -q[1]/qnorm;
  qinv[2] = -q[2]/qnorm;
  qinv[3] = -q[3]/qnorm;

  return qinv;
}


double* quatsum (double q1[], double q2[])
{
  //Compute sum between two quaternions

  double* qsum = new double[4];
  for (int i=0; i<4; i++)
    qsum[i] = 0.0;

  qsum[0] = q1[0] + q2[0];
  qsum[1] = q1[1] + q2[1];
  qsum[2] = q1[2] + q2[2];
  qsum[3] = q1[3] + q2[3];

  return qsum;
}


double* quatdiff (double q1[], double q2[])
{
  //Compute difference between two quaternions (q1 - q2)

  double* qdiff = new double[4];
  for (int i=0; i<4; i++)
    qdiff[i] = 0.0;

  qdiff[0] = q1[0] - q2[0];
  qdiff[1] = q1[1] - q2[1];
  qdiff[2] = q1[2] - q2[2];
  qdiff[3] = q1[3] - q2[3];

  return qdiff;
}


double* quatprod (double q1[], double q2[])
{
  //Compute product between two quaternions

  double* qprod = new double[4];
  for (int i=0; i<4; i++)
    qprod[i] = 0.0;

  qprod[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
  qprod[1] = q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2];
  qprod[2] = q1[0]*q2[2] + q1[2]*q2[0] + q1[3]*q2[1] - q1[1]*q2[3];
  qprod[3] = q1[0]*q2[3] + q1[3]*q2[0] + q1[1]*q2[2] - q1[2]*q2[1];

  return qprod;
}


double* quatprodscalar (double a, double q[])
{
  //Compute product between a scalar and a quaternion

  double* qprodscalar = new double[4];
  for (int i=0; i<4; i++)
    qprodscalar[i] = 0.0;

  qprodscalar[0] = a*q[0];
  qprodscalar[1] = a*q[1];
  qprodscalar[2] = a*q[2];
  qprodscalar[3] = a*q[3];

  return qprodscalar;
}


double quatnorm (double q[])
{
  //Compute norm of a quaternion

  double qnorm;

  qnorm = sqrt(q[0]*q[0] + q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);

  return qnorm;
}



double quatdot (double q1[], double q2[])
{
  //Dot product of two quaternions

  double qprod;

  qprod = q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2] + q1[3]*q2[3];

  return qprod;
}



double* quatMATprod (double matrix[4][4], double q[])
{
  //Multiply a 4x4 matrix by a quaternion

  double* r = new double[4];
  for (int i=0; i<4; i++)
    r[i] = 0.0;

  r[0] = matrix[0][0]*q[0] + matrix[0][1]*q[1]  + matrix[0][2]*q[2]  + matrix[0][3]*q[3];
  r[1] = matrix[1][0]*q[0] + matrix[1][1]*q[1]  + matrix[1][2]*q[2]  + matrix[1][3]*q[3];
  r[2] = matrix[2][0]*q[0] + matrix[2][1]*q[1]  + matrix[2][2]*q[2]  + matrix[2][3]*q[3];
  r[3] = matrix[3][0]*q[0] + matrix[3][1]*q[1]  + matrix[3][2]*q[2]  + matrix[3][3]*q[3];

  return r;
}


double* qRotX (double q[], double theta)
{
  //Rotation of a quaternion of an theta around laboratory X axis
  double* r = new double[4];
  for (int i=0; i<4; i++)
    r[i] = 0.0;

  double ct = cos(0.5*theta);
  double st = sin(0.5*theta);

  double RotMat[4][4];
  
  RotMat[0][0] =  ct;
  RotMat[0][1] = -st;
  RotMat[0][2] = 0.0;
  RotMat[0][3] = 0.0;

  RotMat[1][0] =  st;
  RotMat[1][1] =  ct;
  RotMat[1][2] = 0.0;
  RotMat[1][3] = 0.0;

  RotMat[2][0] = 0.0;
  RotMat[2][1] = 0.0;
  RotMat[2][2] =  ct;
  RotMat[2][3] = -st;

  RotMat[3][0] = 0.0;
  RotMat[3][1] = 0.0;
  RotMat[3][2] =  st;
  RotMat[3][3] =  ct;

  r = quatMATprod(RotMat, q);

  return r;
}


double* qRotY (double q[], double theta)
{
  //Rotation of a quaternion of an theta around laboratory Y axis
  double* r = new double[4];
  for (int i=0; i<4; i++)
    r[i] = 0.0;

  double ct = cos(0.5*theta);
  double st = sin(0.5*theta);

  double RotMat[4][4];
  
  RotMat[0][0] =  ct;
  RotMat[0][1] = 0.0;
  RotMat[0][2] = -st;
  RotMat[0][3] = 0.0;

  RotMat[1][0] = 0.0;
  RotMat[1][1] =  ct;
  RotMat[1][2] = 0.0;
  RotMat[1][3] =  st;

  RotMat[2][0] =  st;
  RotMat[2][1] = 0.0;
  RotMat[2][2] =  ct;
  RotMat[2][3] = 0.0;

  RotMat[3][0] = 0.0;
  RotMat[3][1] = -st;
  RotMat[3][2] = 0.0;
  RotMat[3][3] =  ct;

  r = quatMATprod(RotMat, q);

  return r;
}



double* qRotZ (double q[], double theta)
{
  //Rotation of a quaternion of an theta around laboratory Z axis
  double* r = new double[4];
  for (int i=0; i<4; i++)
    r[i] = 0.0;

  double ct = cos(0.5*theta);
  double st = sin(0.5*theta);

  double RotMat[4][4];
  
  RotMat[0][0] =  ct;
  RotMat[0][1] = 0.0;
  RotMat[0][2] = 0.0;
  RotMat[0][3] = -st;

  RotMat[1][0] = 0.0;
  RotMat[1][1] =  ct;
  RotMat[1][2] = -st;
  RotMat[1][3] = 0.0;

  RotMat[2][0] = 0.0;
  RotMat[2][1] =  st;
  RotMat[2][2] =  ct;
  RotMat[2][3] = 0.0;

  RotMat[3][0] =  st;
  RotMat[3][1] = 0.0;
  RotMat[3][2] = 0.0;
  RotMat[3][3] =  ct;

  r = quatMATprod(RotMat, q);

  return r;
}

double* qrot (double q[], double p[])
{
  //Rotate a point in a quaternion frame (quatrot/invquatrot)
 
  double* rotate = new double[3];
  for (int i=0; i<3; i++)
    {
      rotate[i] = 0.0;
    }

  double v1 =  p[0]*q[1] + p[1]*q[2] + p[2]*q[3];
  double v2 =  p[0]*q[0] - p[1]*q[3] + p[2]*q[2];
  double v3 =  p[0]*q[3] + p[1]*q[0] - p[2]*q[1];
  double v4 = -p[0]*q[2] + p[1]*q[1] + p[2]*q[0];

  rotate[0] = q[1]*v1 + q[0]*v2 - q[3]*v3 + q[2]*v4;
  rotate[1] = q[2]*v1 + q[3]*v2 + q[0]*v3 - q[1]*v4;
  rotate[2] = q[3]*v1 - q[2]*v2 + q[1]*v3 + q[0]*v4;

  return rotate;
}





double* qDirCosX (double q[])
{
  //Director cosine corresponding to scalar product between x and (X, Y, Z) --> 1st row of RotMat(q)
  double* r = new double[3];
  for (int i=0; i<3; i++)
    r[i] = 0.0;

  double qxx = q[1]*q[1];
  double qyy = q[2]*q[2];
  double qzz = q[3]*q[3];

  double qxy = q[1]*q[2];
  double qxz = q[1]*q[3];
  double qyz = q[2]*q[3];

  double qwx = q[0]*q[1];
  double qwy = q[0]*q[2];
  double qwz = q[0]*q[3];


  r[0] = 1.0 - 2.0*(qyy+qzz);
  r[1] = 2.0*(qxy+qwz);
  r[2] = 2.0*(qxz-qwy);

  return r;
}



double* qDirCosY (double q[])
{
  //Director cosine corresponding to scalar product between y and (X, Y, Z) --> 2nd row of RotMat(q)
  double* r = new double[3];
  for (int i=0; i<3; i++)
    r[i] = 0.0;

  double qxx = q[1]*q[1];
  double qyy = q[2]*q[2];
  double qzz = q[3]*q[3];

  double qxy = q[1]*q[2];
  double qxz = q[1]*q[3];
  double qyz = q[2]*q[3];

  double qwx = q[0]*q[1];
  double qwy = q[0]*q[2];
  double qwz = q[0]*q[3];

  r[0] = 2.0*(qxy-qwz);
  r[1] = 1.0 - 2.0*(qxx+qzz);
  r[2] = 2.0*(qyz+qwx);

  return r;
}



double* qDirCosZ (double q[])
{
  //Director cosine corresponding to scalar product between z and (X, Y, Z) --> 3rd row of RotMat(q)
  double* r = new double[3];
  for (int i=0; i<3; i++)
    r[i] = 0.0;

  double qxx = q[1]*q[1];
  double qyy = q[2]*q[2];
  double qzz = q[3]*q[3];

  double qxy = q[1]*q[2];
  double qxz = q[1]*q[3];
  double qyz = q[2]*q[3];

  double qwx = q[0]*q[1];
  double qwy = q[0]*q[2];
  double qwz = q[0]*q[3];

  r[0] = 2.0*(qxz+qwy);
  r[1] = 2.0*(qyz-qwx);
  r[2] = 1.0 - 2.0*(qxx+qyy);

  return r;
}


double qZx (double q[])
{
  //Director cosine corresponding to scalar product between x and (0,0,1)
  double r;

  double qxx = q[1]*q[1];
  double qyy = q[2]*q[2];
  double qzz = q[3]*q[3];

  double qxy = q[1]*q[2];
  double qxz = q[1]*q[3];
  double qyz = q[2]*q[3];

  double qwx = q[0]*q[1];
  double qwy = q[0]*q[2];
  double qwz = q[0]*q[3];

  r = 2.0*(qxz+qwy);

  return r;
}



double qZy (double q[])
{
  //Director cosine corresponding to scalar product between y and (0,0,1)
  double r;

  double qxx = q[1]*q[1];
  double qyy = q[2]*q[2];
  double qzz = q[3]*q[3];

  double qxy = q[1]*q[2];
  double qxz = q[1]*q[3];
  double qyz = q[2]*q[3];

  double qwx = q[0]*q[1];
  double qwy = q[0]*q[2];
  double qwz = q[0]*q[3];

  r = 2.0*(qyz-qwx);

  return r;
}



double qZz (double q[])
{
  //Director cosine corresponding to scalar product between z and (0,0,1)
  double r;

  double qxx = q[1]*q[1];
  double qyy = q[2]*q[2];
  double qzz = q[3]*q[3];

  double qxy = q[1]*q[2];
  double qxz = q[1]*q[3];
  double qyz = q[2]*q[3];

  double qwx = q[0]*q[1];
  double qwy = q[0]*q[2];
  double qwz = q[0]*q[3];

  r = 1.0-2.0*(qxx+qyy);

  return r;
}





double* MatrixToQuaternion (double matrix[3][3], bool active)
{
  std::vector<double> quat;

  //Convert a rotation matrix to the corresponding quaternion

//  std::cout << "[MatrixToQuaternion] beginning " << std::endl; 
//  std::cout << matrix[0][0] << " "<< matrix[0][1] << " " << matrix[0][2] << " " << std::endl; 
//  std::cout << matrix[1][0] << " "<< matrix[1][1] << " " << matrix[1][2] << " " << std::endl; 
//  std::cout << matrix[2][0] << " "<< matrix[2][1] << " " << matrix[2][2] << " " << std::endl; 


  //Check matrix to be normalized
  double n0 = sqrt(matrix[0][0]*matrix[0][0] + matrix[0][1]*matrix[0][1] + matrix[0][2]*matrix[0][2]);
  double n1 = sqrt(matrix[1][0]*matrix[1][0] + matrix[1][1]*matrix[1][1] + matrix[1][2]*matrix[1][2]);
  double n2 = sqrt(matrix[2][0]*matrix[2][0] + matrix[2][1]*matrix[2][1] + matrix[2][2]*matrix[2][2]);
  matrix[0][0] = matrix[0][0]/n0;
  matrix[0][1] = matrix[0][1]/n0;
  matrix[0][2] = matrix[0][2]/n0;
  matrix[1][0] = matrix[1][0]/n1;
  matrix[1][1] = matrix[1][1]/n1;
  matrix[1][2] = matrix[1][2]/n1;
  matrix[2][0] = matrix[2][0]/n2;
  matrix[2][1] = matrix[2][1]/n2;
  matrix[2][2] = matrix[2][2]/n2;


  //Symmetrize Matrix Elements
//  matrix[1][0] = matrix[1][0] - matrix[0][0]*(matrix[0][0]*matrix[1][0]);
//  matrix[1][1] = matrix[1][1] - matrix[0][1]*(matrix[0][1]*matrix[1][1]);
//  matrix[1][2] = matrix[1][2] - matrix[0][0]*(matrix[2][2]*matrix[1][2]);
//  matrix[2][0] = pow(matrix[0][0], matrix[1][0]);
//  matrix[2][1] = pow(matrix[0][1], matrix[1][1]);
//  matrix[2][2] = pow(matrix[0][2], matrix[1][2]);


//  std::cout << "[MatrixToQuaternion] symm " << std::endl; 
//  std::cout << matrix[0][0] << " "<< matrix[0][1] << " " << matrix[0][2] << " " << std::endl; 
//  std::cout << matrix[1][0] << " "<< matrix[1][1] << " " << matrix[1][2] << " " << std::endl; 
//  std::cout << matrix[2][0] << " "<< matrix[2][1] << " " << matrix[2][2] << " " << std::endl; 



//  std::cout << "[MatrixToQuaternion] " << std::endl; 
//  std::cout << matrix[0][0] << " "<< matrix[0][1] << " " << matrix[0][2] << " " << std::endl; 
//  std::cout << matrix[1][0] << " "<< matrix[1][1] << " " << matrix[1][2] << " " << std::endl; 
//  std::cout << matrix[2][0] << " "<< matrix[2][1] << " " << matrix[2][2] << " " << std::endl; 

  double trace = matrix[0][0]+matrix[1][1]+matrix[2][2];
  //  std::cout << "[MatrixToQuaternion] trace " << trace << std::endl; 
  double qw, qx, qy, qz;

  if (trace > 0.0)
    {
      qw = sqrt(1+trace)/2.0;
      qx = (matrix[2][1] - matrix[1][2])/(4*qw);
      qy = (matrix[0][2] - matrix[2][0])/(4*qw);
      qz = (matrix[1][0] - matrix[0][1])/(4*qw);
    }
  else if (matrix[0][0] > matrix[1][1] && matrix[0][0] > matrix[2][2])
    {
      //      std::cout << "a" << std::endl;
      double tmp = matrix[0][0] - matrix[1][1] - matrix[2][2];
      qx = sqrt(1.0+tmp)/2.0;
      qw = (matrix[2][1] - matrix[1][2])/(4*qx);
      qy = (matrix[0][1] + matrix[1][0])/(4*qx);
      qz = (matrix[0][2] + matrix[2][0])/(4*qx);
    }
  else if (matrix[1][1] > matrix[2][2])
    {
      //      std::cout << "b" << std::endl;
      double tmp = matrix[1][1] - matrix[0][0] - matrix[2][2];
      qy = sqrt(1.0+tmp)/2.0;
      qw = (matrix[0][2] - matrix[2][0])/(4*qy);
      qx = (matrix[0][1] + matrix[1][0])/(4*qy);
      qz = (matrix[1][2] + matrix[2][1])/(4*qy);
    }
  else 
    {
      //      std::cout << "c" << std::endl;
      double tmp = matrix[2][2] - matrix[0][0] - matrix[1][1];
      qz = sqrt(1.0+tmp)/2.0;
      qw = (matrix[1][0] - matrix[0][1])/(4*qz);
      qx = (matrix[0][2] + matrix[2][0])/(4*qz);
      qy = (matrix[1][2] + matrix[2][1])/(4*qz);
    }

  if (active == 0)
    {
      qw = -qw;
      qx = qx;
      qy = qy;
      qz = qz;
    }

  //normalize
//  double magnitude = sqrt(qw*qw + qx*qx + qy*qy + qz*qz);
//  qw = qw/magnitude;
//  qx = qx/magnitude;
//  qy = qy/magnitude;
//  qz = qz/magnitude;

  quat.push_back(qw);
  quat.push_back(qx);
  quat.push_back(qy);
  quat.push_back(qz);

  double* q = new double[4];
  for (int i=0; i<4; i++)
    {
      q[i] = 0.0;
    }

  q[0] = quat[0];
  q[1] = quat[1];
  q[2] = quat[2];
  q[3] = quat[3];

  return q;
}



double** QuaternionToMatrix (double q[])
{
  //Convert a quaternion into a rotation matrix

  double** mat = new double*[3];
  for (int i=0; i<3; i++)
    mat[i] = new double[3];

  double qxx = q[1]*q[1];
  double qyy = q[2]*q[2];
  double qzz = q[3]*q[3];

  double qxy = q[1]*q[2];
  double qxz = q[1]*q[3];
  double qyz = q[2]*q[3];

  double qwx = q[0]*q[1];
  double qwy = q[0]*q[2];
  double qwz = q[0]*q[3];

  mat[0][0] = 1.0 - 2.0*(qyy+qzz);
  mat[0][1] = 2.0*(qxy+qwz);
  mat[0][2] = 2.0*(qxz-qwy);

  mat[1][0] = 2.0*(qxy-qwz);
  mat[1][1] = 1.0 - 2.0*(qxx+qzz);
  mat[1][2] = 2.0*(qyz+qwx);

  mat[2][0] = 2.0*(qxz+qwy);
  mat[2][1] = 2.0*(qyz-qwx);
  mat[2][2] = 1.0 - 2.0*(qxx+qyy);

  return mat;
}



double* eul2quat(double alpha, double beta, double gamma)
{
  //Convert 3 Euler angles into a quaternion 

  double* q = new double[4];
  for (int i=0; i<4; i++)
      q[i] = 0.0;

  double hbeta = 0.5*beta;
  double halpga = 0.5*(alpha+gamma);
  double halmga = 0.5*(alpha-gamma);
  double cohbe = cos(hbeta);
  double sihbe = sin(hbeta);

  q[0] =  cohbe*cos(halpga);
  q[1] = -sihbe*sin(halmga);
  q[2] =  sihbe*cos(halmga);
  q[0] =  cohbe*sin(halpga);

  return q;
}


double* vec2quat (double v[3])
{
  //Convert a vector into a quaternion
 
  double* q = new double[4];
  for (int i=0; i<4; i++)
      q[i] = 0.0;

  double* vrs = versor(v);
  double alpha, beta, gamma;

  if ( (1.0-abs(vrs[2])) > 0.001 )
    {
      beta = acos(vrs[2]);
      double s = sin(beta);
      alpha = acos(vrs[0]/s);
      
      if ( vrs[1] < 0.0 )
	alpha = 2.0*PI - alpha;

    }
  else
    {
	
      if ( vrs[2] > 0.0 )
	beta = 0.0;
      else
	beta = PI;

      alpha = 0.0;
      gamma = 0.0;
    }
    
  q = eul2quat(alpha, beta, gamma);
    
  return q;
}

