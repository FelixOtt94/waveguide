 
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "Source/Colsamm.h"
using namespace ::_COLSAMM_;


typedef struct{
  double x;
  double y;
} point; 

typedef struct{
  uint16_t v0;
  uint16_t v1;
  uint16_t v2;
} face; 

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
    for(int j=0; j<= 1976; j++) 
    {
      infile >> v0 >> v1 >> v2;
      faces[j].v0 = v0;
      faces[j].v1 = v1;
      faces[j].v2 = v2;
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

void computeA_M( point* points, face* faces ){

    ELEMENTS::Triangle my_element;
    /*std::vector< std::vector< double > > my_local_matrix;

    std::vector<double> corners(6, 0.0);

    for(int i=0; i<1977; i++){
        corners[0] = points[faces[i].v0].x; corners[1] = points[faces[i].v0].y;
        corners[2] = points[faces[i].v1].x; corners[3] = points[faces[i].v0].y;
        corners[4] = points[faces[i].v2].x; corners[5] = points[faces[i].v0].y;
        my_element(corners);


        //my_local_matrix = my_element.integrate(grad(v_()) * grad(w_()));


    }

*/

}

//double my_func(double x, double y);

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
  
  delta = atof( argv[1] );
  //epsilon = atof( argv[2] );
  
  readInput( points, faces );
  
  setK( points, k, delta );
  
  testK( points, k );
  exit( EXIT_SUCCESS );
  
}
