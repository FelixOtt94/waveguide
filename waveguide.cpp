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



static double delta = 0.0;

typedef struct{
  double x;
  double y;
} point; 

typedef struct{
  uint16_t v[3];
} face; 



// computes the skalarprodukt of a and b.
double skalarprodukt(std::vector<double>& a, std::vector<double>&  b){
    double skalarprodukt = 0.0;

    for(int i=0; i<1039; i++){
            skalarprodukt += a[i] * b[i];
        
    }
    return skalarprodukt;
}

// calculates the l2-norm of the given vector r
double residuum( std::vector<double> &r) {
	double residuum = 0.0;
	double sum = 0.0;
	for (int i = 0; i < 1039; i++) {
			sum += r[i] * r[i];
		
	}
	residuum = sqrt( sum) / 1039;
	return residuum;
}

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



void setK( point* points, std::vector<double>& k){
  
  for( int i =0; i < 1039; i++ ){
      k[i] = (100.0 + delta )* std::exp( -50.0 * ((points[i].x*points[i].x) + (points[i].y*points[i].y)) ) - 100.0;
  }
}




//Calculates k^2. Needed for function pointer in function computeA_M()
double my_func(double x, double y){
    double k = 0;
    k = (100.0 + delta )* std::exp( -50.0 * ((x*x) + (y*y)) ) - 100.0;
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
		matrixA[faces[i].v[l]][faces[i].v[k]] += my_local_matrixA[l][k];

		matrixM[faces[i].v[l]][faces[i].v[k]] += my_local_matrixM[l][k];

	   }
        }

    }



}


void cg(std::vector<double>&  f_x_y,
	std::vector<double>& grid, 
	std::vector<std::map<uint16_t, double>>& matrixA,
	double epsilon){
  
    double alpha = 0.0;
    double delta_1 = 0.0;
    double betta = 0.0;
    double delta_0 = 0.0;
    double tmp = 0.0;
    std::vector<double> r(1039);
    std::vector<double> d(1039);
    std::vector<double> z(1039);

    // 1)
    for(int i=0; i<1039; i++){
      for (std::map<uint16_t, double>::iterator it=matrixA[i].begin(); it!=matrixA[i].end(); ++it){
			  tmp += it->second * grid[it->first];
      }
      r[i] = f_x_y[i] - tmp;
      tmp = 0.0;
    }

    // 2)
    delta_0 = skalarprodukt(r, r);

    // 3)
    if(residuum(r) < epsilon){
        //std::cout << "number of needed iterations: 0" << std::endl;
        return;
    }

    // 4)
    for(int i=0; i<1039; i++){
            d[i] = r[i];

    }

    // 5-15)
    for(int k=1; k<=1000000; k++){
        // 6)
        for(int i=0; i<1039; i++){
	            for (std::map<uint16_t, double>::iterator it=matrixA[i].begin(); it!=matrixA[i].end(); ++it){
			  tmp += it->second * d[it->first];
	    }
                z[i] = tmp;
                tmp = 0.0;
           
        }

        // 7)
        alpha = delta_0 / skalarprodukt(d, z);

        // 8)
        for(int i=0; i<1039; i++){
                grid[i] = grid[i] + alpha * d[i];
            
        }

        // 9)
        for(int i=0; i<1039; i++){
                r[i] = r[i] - alpha * z[i];
            
        }

        // 10)
        delta_1 = skalarprodukt(r, r);

        // 11)
        if(residuum(r) < epsilon){
            //std::cout << "number of needed iterations: " << k << std::endl;
            return;
        }

        // 12)
        betta = delta_1/delta_0;

        // 13)
        for(int i=0; i<1039; i++){
                d[i] = r[i] + betta * d[i];
            
        }

        // 14)
        delta_0 = delta_1;
    }
    //std::cout << "number of needed iterations: " << 1000000 << std::endl;
}

void inverse_power(std::vector<std::map<uint16_t, double>>& matrixA,std::vector<std::map<uint16_t, double>>& matrixM, 	std::vector<double>& u, double epsilon){
    double lambda = 1.1;
    double lambda_old = 1;
    std::vector<double> f (1039);
    double abort = 1e-10;
    double tmp = 0.0;
    double tmp_lambda = 0.0;
    double length = 0.0;
    int counter = 0;


    while(fabs((lambda-lambda_old)/lambda_old) > abort ){
        //2
        lambda_old = lambda;

        //3
        for(int i = 0; i<1039; ++i){
            for (std::map<uint16_t, double>::iterator it=matrixM[i].begin(); it!=matrixM[i].end(); ++it){
                tmp += (it-> second)*u[it->first];
            }
            f[i] = tmp;
            tmp = 0.0;
        }

        //4
        cg(f, u, matrixA, epsilon);

        //5
        for(int i = 0; i < 1039; ++i){
            length += u[i]*u[i];
        }
        length = sqrt(length);
        for(int i = 0; i < 1039; ++i){
            u[i] = u[i]/length;
        }

        //6
        for(int i = 0; i<1039; ++i){
            for (std::map<uint16_t, double>::iterator it=matrixA[i].begin(); it!=matrixA[i].end(); ++it){
                tmp += (it-> second)*u[it->first];
            }
            f[i] = tmp;
            tmp = 0.0;
        }
        for(int i = 0; i<1039; ++i){
            tmp_lambda += u[i]*f[i];
        }
        for(int i = 0; i<1039; ++i){
            for (std::map<uint16_t, double>::iterator it=matrixM[i].begin(); it!=matrixM[i].end(); ++it){
                tmp += (it-> second)*u[it->first];
            }
            f[i] = tmp;
            tmp = 0.0;
        }
        for(int i = 0; i<1039; ++i){
            tmp += u[i]*f[i];
        }
        lambda = tmp_lambda/tmp;


        tmp = 0.0;
        tmp_lambda = 0.0;
        length = 0.0;
        std::cout << "Iteration: " << counter++ << ", lambda: " << lambda << std::endl;

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

void printSolution(std::vector<double>& u, point * points, std::string name){
    std::ofstream out (name);
    for(int i = 0; i < 1039; i++){
        out << points[i].x << " " << points[i].y << " " << u[i] << std::endl;
    }
}

int main( int args, char** argv ){
  
  if( args != 3 ){
    std::cout << "Usage: ./waveguide delta epsilon" << std::endl;
    exit( EXIT_SUCCESS );
  }
  

  double epsilon = 0.0;
  point points[1039];
  face faces[1977];
  std::vector<double> k(1039);
  std::vector<double> u(1039, 1.0);

  // matrices
  std::vector< std::map<uint16_t, double> > matrixA(1039);
  std::vector< std::map<uint16_t, double> > matrixM(1039);

  
  delta = atof( argv[1] );
  epsilon = atof( argv[2] );
  
  readInput( points, faces );

  setK( points, k);

  computeA_M( points, faces, matrixA, matrixM);
  printMat( std::string("A.txt"), 1039, matrixA );
  printMat( std::string("M.txt"), 1039, matrixM );

  inverse_power(matrixA, matrixM, u, epsilon);
  
  printSolution(u, points, "eigenmode.txt");
  printSolution(k, points, "ksq.txt");

  exit( EXIT_SUCCESS );
  
}
