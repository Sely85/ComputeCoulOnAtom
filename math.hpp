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



double* distance(double crd1[], double crd2[])
{
  //Compute distance vector between two vectors (double[3])
 
  double* distance = new double[3];
  for (int i=0; i<3; i++)
    {
      distance[i]=0.0;
    }

  distance[0] = crd2[0] - crd1[0];
  distance[1] = crd2[1] - crd1[1];
  distance[2] = crd2[2] - crd1[2];
  //std::cout << " distance " << distance[0] << " " << distance[1] << " " << distance[2] << std::endl;

  return distance;
}



double norm(double crd[])
{
  //Compute norm value of a vector (double[3])

  double norm=0.0;

  norm = sqrt(crd[0]*crd[0] + crd[1]*crd[1] + crd[2]*crd[2]);

  return norm;
}




double* crossprod(double crd1[], double crd2[])
{
  //Compute cross product between two vectors (double[3])

  double* cross = new double[3];
  for (int i=0; i<3; i++)
    {
      cross[i]=0.0;
    }

  cross[0] = crd1[1]*crd2[2]-crd1[2]*crd2[1];
  cross[1] = crd1[2]*crd2[0]-crd1[0]*crd2[2];
  cross[2] = crd1[0]*crd2[1]-crd1[1]*crd2[0];

  return cross;
}


double dotprod(double crd1[], double crd2[])
{
  //Compute dot product between two vectors (double[3])

  double dot;

  dot = crd1[0]*crd2[0] +  crd1[1]*crd2[1] + crd1[2]*crd2[2];

  return dot;
}



double* matvecprod(double matrix[3][3], double vector[])
{
  //Compute product between a 3x3 matrix and a vector 

  double* prod = new double[3];
  for (int i=0; i<3; i++)
    {
      prod[i] = 0.0;
    }

  prod[0] = matrix[0][0]*vector[0] + matrix[0][1]*vector[1] + matrix[0][2]*vector[2];
  prod[1] = matrix[1][0]*vector[0] + matrix[1][1]*vector[1] + matrix[1][2]*vector[2];
  prod[2] = matrix[2][0]*vector[0] + matrix[2][1]*vector[1] + matrix[2][2]*vector[2];

  return prod;
}



double* versor(double vector[])
{
  //Compute distance vector between two vectors (double[3])

  double* versor = new double[3];
  for (int i=0; i<3; i++)
      versor[i]=0.0;

  double nrm = norm(vector);
  
  versor[0] = vector[0]/nrm;
  versor[1] = vector[1]/nrm;
  versor[2] = vector[2]/nrm;

  return versor;
}


double compute_angle(double x1[], double x2[])
{
  double angle = 0.0;

  double num = dotprod(x1, x2);
  double norm1 = norm(x1);
  double norm2 = norm(x2);

  angle = num/(norm1*norm2);

  return angle;
}

