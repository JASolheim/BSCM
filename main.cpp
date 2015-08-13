/*
This is a simple test program whose sole purpose is to exercise the class BSCM::Spline.
Author: Jeff Solheim   Date: 03 JAN 2013

  In this example :
  M = order = 3
  N = # of collocation points = 3
  *Physical* boundaries extend from Xmin = 2 to Xmax = 5 :
        __________      // line over physical region
    0  1  2  3  4  5  6  7
  // Knots are at x = 0, 1, 2, 3, 4, 5, 6, & 7.
  The knots at x = 0, 1, 6, & 7 fall outside the *physical* boundaries.
  Collocation points will be located at 2.5, 3.5, 4.5
  (halfway between each consecutive pair of physical boundary knots).
*/

#include <iostream>
#include <iomanip>
#include "Spline.h"
#include <unsupported/Eigen/MatrixFunctions>

using namespace std ;

int main ( int argc, char *argv[] )
  {
  const unsigned int         splineOrder  =  3;//5;//3 ;
  const double               KNOT_ARRAY []  =  { 0, 1, 2, 3, 4, 5, 6, 7 } ;
  //const double               KNOT_ARRAY []  =  { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 } ;
  const std::vector<double>  knotVector ( KNOT_ARRAY, KNOT_ARRAY + sizeof(KNOT_ARRAY) / sizeof(double) ) ;

   Eigen::MatrixXi  boundaryConditionsMatrix ( 2, 3 ) ;
   boundaryConditionsMatrix  <<  1, 0, 0,   //  0th derivative at left  boundary to be zero
                         0, 1, 0 ;  //  1st derivative at right boundary to be zero
  //Eigen::MatrixXi  boundaryConditionsMatrix ( 4, 5 ) ;
  //boundaryConditionsMatrix  <<  1, 0, 0, 0, 0,
  //                      0, 1, 0, 0, 0,
  //                      1, 0, 0, 0, 0,
  //                      0, 1, 0, 0, 0 ;
  // instantiate the Spline class ...
  BSCM::Spline     testSpline ( splineOrder, knotVector, boundaryConditionsMatrix ) ;
  // obtain the matrix representation of the second derivative ...
  Eigen::MatrixXd  D  =  testSpline.operatorMatrix ( 2 ) ;
  double thermalDiffusivity  =  0.5 ;
  Eigen::MatrixXd  A  =  (thermalDiffusivity * testSpline.operatorMatrix(2)).exp() ;
  Eigen::VectorXd  u ;
  u.resize ( testSpline.numKnots - 2 * testSpline.order + 1 ) ; // 8 - 2*3 + 1 = 3
  u << 1, 0, 0.5 ; // initial temperatures at collocation points
  // iterate through several time periods to simulate diffusion of heat within rod ...
  for ( int t = 1 ; t <= 5 ; t++ )
    {
    u = A * u ;
	cout << u << endl << endl ;
    }

exit(0);

  cout << "=====================================================================\n" ;
  cout << "Collocation points collocationX.at(i) ..." ;
  for ( unsigned int i = 0 ; i < testSpline.collocationX.size() ; i ++ )
   cout << testSpline.collocationX.at(i) << endl ;

  cout << "=====================================================================\n" ;
  cout << "B(k,i,x) = ..." ;
  for ( size_t k = 1 ; k <= testSpline.order ; k ++ )
    for ( size_t i = 0 ; i < (testSpline.knotX.size() - k - 1) ; i ++ )
      for ( double x = 1.41 ; x <= 6.99 ; x += 0.40 )
        cout << "B(" << k << "," << i << "," << x << ") \t" << testSpline.B(k,i,x) << endl ;

  cout << "=====================================================================\n" ;
  cout << "B(k,i,alpha) = ..." ;
  for ( size_t k = 1 ; k <= testSpline.order ; k ++ )
    for ( size_t i = 0 ; i < (testSpline.knotX.size() - k - 1) ; i ++ )
      {
      for ( size_t alpha = 0 ; alpha < testSpline.collocationX.size() ; alpha ++ )
        {
        cout << "B(" << k << "," << i << "," << alpha << ") [" << testSpline.collocationX.at(alpha) << "] \t"
             << testSpline.B(k,i,alpha) << endl ;
        }
      cout << "--------------------------" << endl ;
      }

  cout << "=====================================================================\n" ;
  cout << "D_B(p,k,i,x) = ..." ;
  for ( size_t p = 0 ; p < testSpline.order ; p ++ )
    {
    for ( size_t k = (p+1) ; k <= testSpline.order ; k ++ )
      for ( size_t i = 0 ; i < (testSpline.knotX.size() - k - 1) ; i ++ )
        for ( double x = 1.41 ; x <= 6.99 ; x += 0.40 )
          cout << "D_B(" << p << "," << k << "," << i << "," << x << ") \t"
               << testSpline.D_B(p,k,i,x) << endl ;
    cout << "=============================================================================" << endl ;
    }

  std::cout << std::setprecision(2) << setw(6) ;
  cout << "=====================================================================\n" ;
  cout << "B_matrix is\n" << testSpline.B_matrix << endl ;

  cout << "=====================================================================\n" ;
  cout << "K_matrix is\n" << testSpline.K_matrix << endl ;

  cout << "=====================================================================\n" ;
  cout << "beta_matrix is\n" << testSpline.beta_matrix << endl ;

  cout << "=====================================================================\n" ;
  cout << "test value D_B == "
       << testSpline.D_B ( testSpline.order-1, testSpline.order, 8, 10.0 ) << endl ;
// failed:  assert ( i < (knotX.size() - k) ) ;  //  i can range 0 to # of knots - k


  cout << "=====================================================================\n" ;
  for ( double x = 1.0 ; x <= 8.0 ; x +=0.5 )
    cout << x << "\t" << testSpline.D_B ( 2, testSpline.order, 1, x ) << endl ;
  // Maple equivalent is `&Delta;B`(2, M, 2, x)

  cout << "=====================================================================\n" ;
  cout << "operatorMatrix ( 2 ) is\n" << testSpline.operatorMatrix(2) << endl ;

  cout << "=====================================================================\n" ;

for ( size_t k = 1 ; k <= testSpline.order ; k ++ )
  {
  for ( size_t i = 0 ; i <= (testSpline.N + 2*testSpline.order - k - 2) ; i ++ )
    {
    for ( size_t alpha = 0 ; alpha <= (testSpline.N - 1) ; alpha ++ )
      {
      double  a  =  testSpline.B ( k, i, testSpline.collocationX[alpha] ) ;
      double  b  =  testSpline.B ( k, i, alpha ) ;
      //double  b  =  testSpline.B_k_i_alpha.at(k).at(i).at(alpha) ;
      double  c  =  a - b ;
      cout << k << i << alpha << " (" ;
      cout << c ;
      cout << ")" << endl ;
      }
    cout << endl ;
    }
  cout << "----------------------------------------------" ;
  cout << endl ;
  }
  cout << "=====================================================================\n" ;
  } // end main
