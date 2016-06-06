#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <map>
#include <utility>
#include "Source/Colsamm.h"
using namespace ::_COLSAMM_;


typedef struct{
  double x;
  double y;
} point; 

typedef struct{
  uint16_t v[3];
} face; 

//typedef double (*func_t)(double, double);

void readInput( point* points, face* faces ){
  double i, x, y;
  double v0, v1, v2;
  
  // read points
  std::string line;
  std::ifstream infile ("unit_circle.txt");
  if (infile.is_open())
  {
    std::getline(infile, line);
    std::getline(infile, line);
    for(int j=0; j<= 1038; j++) 
    {
      infile >> i >> x >> y;
      points[j].x = x;
      points[j].y = y;
    }
    
    // read faces
    std::getline(infile, line);
    std::getline(infile, line);
        std::getline(infile, line);


    
    for(int j=0; j< 1976; j++) 
    {
      infile >> v0 >> v1 >> v2;
      faces[j].v[0] = v0;
      faces[j].v[1] = v1;
      faces[j].v[2] = v2;
    }
    
    infile.close();
  }

  else std::cout << "Unable to open file"; 
}

void setK( point* points, double* k, double delta ){
  
  for( int i =0; i <= 1039; i++ ){
      k[i] = (100.0 + delta )* std::exp( -50.0 * ((points[i].x*points[i].x) + (points[i].y*points[i].y)) ) - 100.0;
  }
}

void testK( point* points, double* k ){
    double x, y, k2;
     std::ifstream infile ("ksq-ref.txt");
    for(int j=0; j<= 1038; j++) 
    {
      infile >> x >> y >> k2;
      if( points[j].x != x || points[j].y != y || k[j] != k2 ){
        std::cout << "Error: x: " << x << " y: " << y << " k: " << k2 << std::endl;
        std::cout << "Expected: x: " << points[j].x << " y: " << points[j].y << " k: " << k[j] << std::endl;
      }
    }
}

double my_func(double x, double y){
    double k = 0;
    k = (100.0 + 0.01 )* std::exp( -50.0 * ((x*x) + (y*y)) ) - 100.0;
    return k;
}


void computeA_M( point* points, face* faces, std::vector<std::map<uint16_t, double>>& matrixA, std::vector<std::map<uint16_t, double>>& matrixM){

    ELEMENTS::Triangle my_element;
    std::vector< std::vector< double > > my_local_matrixA;
    std::vector< std::vector< double > > my_local_matrixM;

    std::vector<double> corners(6, 0.0);

    for(int i=0; i<1976; i++){
        corners[0] = points[faces[i].v[0]].x; corners[1] = points[faces[i].v[0]].y;
        corners[2] = points[faces[i].v[1]].x; corners[3] = points[faces[i].v[1]].y;
        corners[4] = points[faces[i].v[2]].x; corners[5] = points[faces[i].v[2]].y;
        my_element(corners);


        my_local_matrixA = my_element.integrate(grad(v_()) * grad(w_()) - func<double>(my_func)*w_()*v_());
        my_local_matrixM = my_element.integrate( v_() * w_());

        for(int k=0; k<3; k++){
	   for(int l=0; l<3; l++){
	     if(faces[i].v[l] == 0 && faces[i].v[k] == 0 ){
		std::cout << "i: " << i << " localmatrix : " << my_local_matrixA[0][0] << std::endl;
	     }
		matrixA[faces[i].v[l]][faces[i].v[k]] += my_local_matrixA[l][k];

		matrixM[faces[i].v[l]][faces[i].v[k]] += my_local_matrixM[l][k];

	   }
        }

    }



}

void Jaccobi(  std::vector<std::map<uint16_t, double>>& &u,
	       std::vector<std::map<uint16_t, double>>& matrixA,
	       std::vector<std::map<uint16_t, double>>& u_neu, 
	       std::vector<std::map<uint16_t, double>>& f, 
	       double h, double epsilon){

  double tmp = 0.0;
  int count = 0;
	for (int iterations = 0; iterations<numIterations; iterations++){
		for (int m = 1; m<1039; m++){
			for (std::map<uint16_t, double>::iterator it=matrixA[i].begin(); it!=matrixA[i].end(); ++it){
			  tmp += it->second * u[it->first];
			  count++;
				//u_neu(q, m) = (1.0 / 4.0) * (h*h*f(q, m) + (u(q - 1, m) + u(q + 1, m) + u(q, m - 1) + u(q, m + 1)));
			}
			tmp += f[m];
			u_neu[m] = tmp / count;
			tmp =0.0;
			count= 0;
		}
	}


    

}


void printMat( std::string name, int dim, std::vector<std::map<uint16_t, double>>& matrix ){
  std::ofstream out( name );

  for(int i=0; i<dim; i++ ){
    for (std::map<uint16_t, double>::iterator it=matrix[i].begin(); it!=matrix[i].end(); ++it){
      out << i << " " << it->first << " " << it->second << std::endl;
    }
  }
}


int main( int args, char** argv ){
  
  if( args != 3 ){
    std::cout << "Usage: ./waveguide delta epsilon" << std::endl;
    exit( EXIT_SUCCESS );
  }
  
  double delta = 0.0;
  //double epsilon = 0.0;
  point points[1039];
  face faces[1977];
  double k[1039];

  // matrices
  std::vector< std::map<uint16_t, double> > matrixA(1039);
  std::vector< std::map<uint16_t, double> > matrixM(1039);

  
  delta = atof( argv[1] );
  //epsilon = atof( argv[2] );
  
  readInput( points, faces );

   setK( points, k, delta );

  computeA_M( points, faces, matrixA, matrixM);
  printMat( std::string("A.txt"), 1039, matrixA );
  printMat( std::string("M.txt"), 1039, matrixM );
  

  //testK( points, k );
  exit( EXIT_SUCCESS );
  
}
