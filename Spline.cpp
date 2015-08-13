/**
 * @file    Spline.cpp
 * @author  Jeff Solheim <JASolheim@FHSU.edu>
 * @version  1.0
 *
 * @section LICENSE
 * This program is distributed WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * @section DESCRIPTION
 * File Spline.cpp contains the declaration of the Spline class
 * of the Basis Spline Collocation Method (BSCM).
 */

#include <cassert>
#include <cmath>
#include "Spline.h"

using namespace BSCM ;

// ================================================================================================

// constructor
Spline::Spline ( size_t order, std::vector<double> knotX, Eigen::MatrixXi K_matrix )
  {
  assert ( (order % 2) == 1 ) ;
  assert ( (MIN_ORDER <= order) && (order <= MAX_ORDER) ) ;
  assert ( ((2 * order) <= knotX.size()) && (knotX.size() <= MAX_NUMBER_KNOTS) ) ;

  // Save parameters' values as instance variables.
  this->order     =  order ;
  this->knotX     =  knotX ;
  this->numKnots  =  knotX.size() ;
  this->K_matrix  =  K_matrix ;
  this->xMin      =  knotX [ order - 1 ] ;
  this->xMax      =  knotX [ numKnots - order ] ;

  // Determine collocation points within physical boundaries.
  // Also, confirm that, *within physical boundaries*, each knot is strictly less than its successor.
  for ( int i = (order - 1) ; i < (numKnots - order) ; i ++ )
    {
    assert ( knotX[i] < knotX[i+1] ) ;
    collocationX.push_back ( (knotX[i] + knotX[i+1]) / 2 ) ;
    } // end for i loop
  assert ( collocationX.size() == (numKnots - (2 * order) + 1) ) ; // # of knots = N + 2M - 1
  this->N  =  collocationX.size() ;

  // Prepare vector that is to store values of B(k,i,alpha).
  B_k_i_alpha  =  std::vector< std::vector< std::vector<double> > >
    ( order + 1, // index k = 0 will not be used; values of k = 1 .. M will be used
      std::vector< std::vector<double> >(0) // initially empty vectors for all indices k
     ) ;
  // Values in B_k_i_alpha all initially Not-a-Number (NaN), but those values
  // are replaced as they are calculated.
  for ( size_t k = 1 ; k <= order ; k ++ )
    B_k_i_alpha.at(k)
      =  std::vector< std::vector<double> >
         ( N + 2*order - k - 1, std::vector<double>(N,std::numeric_limits<double>::quiet_NaN()) ) ;

  // Assign values of B(M,i,alpha) to B_matrix.
  B_matrix  =  Eigen::MatrixXd::Zero ( N, (N + order - 1) ) ;
  for ( size_t alpha = 0 ; alpha < static_cast<size_t>(B_matrix.rows()) ; alpha ++ )
    for ( size_t i = 0 ; i < static_cast<size_t>(B_matrix.cols()) ; i ++ )
      B_matrix(alpha,i)  =  B ( order, i, alpha ) ;

  // Assign values within beta_matrix according to Umar's Equation (18), p. 432.
  beta_matrix  =  Eigen::MatrixXd::Zero ( (order - 1), (order + N - 1) ) ;
  for ( int r = 0 ; r < beta_matrix.rows() ; r ++ )
    {
    // First half of rows are evaluated at left boundary; second half at right boundary.
    // xMin would be left boundary of physical region; use xMax for right boundary.
    double  x  =  ( (r < (order / 2)) ? (xMin) : (xMax) ) ;
    for ( int i = 0 ; i < beta_matrix.cols() ; i ++ )
      {
      double  sum  =  0.0 ;
      for ( size_t p = 0 ; p < order ; p ++ )
        sum  +=  ( K_matrix(r,p) * D_B ( p, order, i, x ) ) ;
      beta_matrix ( r, i )  =  sum ;
      } // end for i loop
    } // end for r loop

  // Create & assign B_tilde_matrix of Umar's Equation (20), p. 433.
  B_tilde_matrix  =  Eigen::MatrixXd::Zero ( (N + order - 1), (N + order - 1) ) ;
  B_tilde_matrix << B_matrix, beta_matrix ;

  // Create & assign C_tilde_matrix of Umar's Equation (22), p. 433.
  // Note: one could also consider using Eigen::FullPivLU
  //       if one wants to confirm that B_tilde_matrix is invertible.
  Eigen::PartialPivLU<Eigen::MatrixXd>  thePartialPivLU ( B_tilde_matrix ) ;
  C_tilde_matrix  =  thePartialPivLU.inverse() ;
  } // end constructor

// ================================================================================================

double Spline::B ( size_t k, size_t i, size_t alpha )
  // Evaluate basis function B(k,i) at x = the alpha'th collocation point.
  {
  assert ( (1 <= k) && (k <= order) ) ;  //  k can range 1 to order M
  assert ( i < (numKnots - k) ) ;        //  i can range 0 to # of knots - k - 1
  assert ( alpha < N ) ;                 //  alpha can range 0 to (N-1)

  // If a # has already been calculated for B(k,i,alpha), then just return that number.
  if ( ! std::isnan( B_k_i_alpha[k][i][alpha] ) )
    return  ( B_k_i_alpha[k][i][alpha] ) ;

  // If execution reaches this point, then no # has yet been calculated for B(k,i,alpha).
  assert ( std::isnan(B_k_i_alpha[k][i][alpha]) ) ;

  // B(k,i) will be evaluated at x = the alpha'th collocation point.
  double x  =  collocationX[alpha] ;

  if ( k == 1 )
    return (  B_k_i_alpha[1][i][alpha]  =  B ( 1, i, x )  ) ;

  assert ( k >= 2 ) ;

  double  firstTermB   =  B ( k-1, i  , alpha ) ;
  double  secondTermB  =  B ( k-1, i+1, alpha ) ;

  double firstTerm   =  firstTermB   *  ( (x - knotX[i]) / (knotX[k+i-1] - knotX[i]) ) ;
  double secondTerm  =  secondTermB  *  ( (knotX[k+i] - x) / (knotX[k+i] - knotX[i+1]) ) ;

  return  (  B_k_i_alpha[k][i][alpha]  =  ( firstTerm + secondTerm )  ) ;
  } // end B(k,i,alpha)

// ================================================================================================

double Spline::B ( size_t k, size_t i, double x )
  // Evaluate basis function B(k,i) at arbitrary real x.
  {
  assert ( (1 <= k) && (k <= order) ) ;   //  k can range 1 to order M
  assert ( i < (numKnots - k) ) ;         //  i can range 0 to # of knots - k
  assert ( (knotX.front() <= x) && (x <= knotX.back() ) ) ;
  //  B is defined only for those x between the first & last knots.

  if ( (x < knotX[i]) || (x > knotX[i+k]) )  //  B(k,i,x) falls off to zero to its left & right
    return  0.0 ;

  if ( k >= 2 )
    {
    // see Umar Equation (1), p. 428
    double firstTerm   =  B(k-1,i  ,x) * ( (x - knotX[i]) / (knotX[k+i-1] - knotX[i]) ) ;
    double secondTerm  =  B(k-1,i+1,x) * ( (knotX[k+i] - x) / (knotX[k+i] - knotX[i+1]) ) ;
    return  ( firstTerm + secondTerm ) ;
    } // end if

  // When k == 1, B(k,i) is a step function; see Umar Equation (2), p. 428.
  return  ( ( (knotX[i] <= x) && (x < knotX[i+1]) ) ? (1.0) : (0.0) ) ;
  } // end function B

// ================================================================================================

double Spline::C ( size_t k, size_t i, double x )
  // Implements Umar's Equation (5), p. 429.
  {
  assert ( k >= 1 ) ;

  if ( k == 1 ) // 0th derivative, which is the basis function B, itself
    return B ( 1, i, x ) ;

  double firstTerm   =  C(k-1,i  ,x) / ( knotX[k+i-1] - knotX[i  ] ) ;
  double secondTerm  =  C(k-1,i+1,x) / ( knotX[k+i  ] - knotX[i+1] ) ;
  return  (k-1) * ( firstTerm - secondTerm ) ;
  } // end function C

// ================================================================================================

double Spline::D_B ( size_t p, size_t k, size_t i, double x )
  // Evaluate p'th derivative of B(k,i) at arbitrary real x,
  // for k in the range (p+1) .. M.
  {
  assert ( p < order ) ;                     // p can range 0 .. (M-1)
  assert ( ((p+1) <= k) && (k <= order) ) ;  // k can range (p+1) .. M
  assert ( i < (numKnots - k) ) ;            // i can range 0 .. (numKnots - k - 1)
  assert ( (knotX.front() <= x) && (x <= knotX.back() ) ) ;
  //  D_B is defined only for those x between the first & last knots.

  if ( k >= (p+2) )  //  use Umar's Equation (4), p. 428
    {
    double  firstTermRatio  =  ( x - knotX[i] ) / ( knotX[k+i-1] - knotX[i] ) ;
    double  secondTermRatio  =  ( knotX[k+i] - x ) / ( knotX[k+i] - knotX[i+1] ) ;
    double  firstTermDeriv   =  D_B ( p, k-1, i  , x ) ;
    double  secondTermDeriv  =  D_B ( p, k-1, i+1, x ) ;
    return  ( static_cast<double>(k-1) / static_cast<double>(k-p-1) )
          * ( ( firstTermRatio * firstTermDeriv ) + ( secondTermRatio * secondTermDeriv ) ) ;
    } // end if

  assert ( k == (p+1) ) ;  //  in other words, p = (k-1)

  // Use Umar's Equations (5) & (6) to calculate D_B ( k-1, k, i, x ) = C ( k, i, x ).
  return  C ( k, i, x ) ;
  } // end function D_B

// ================================================================================================

Eigen::MatrixXd  Spline::operatorMatrix ( size_t derivativeOrder )
  // Determine the matrix representation of differentiation operator.
  {
  Eigen::MatrixXd  returnMatrix  =  Eigen::MatrixXd::Zero ( N, N ) ;
  double   sum ;
  for ( size_t alpha = 0 ; alpha < N ; alpha ++ )
    for ( size_t beta = 0 ; beta < N ; beta ++ )
      {
      sum  =  0.0 ;
      for ( int i = 0 ; i < (order + N - 1) ; i ++ )
        sum  +=  ( C_tilde_matrix(i,beta) * D_B(derivativeOrder,order,i,collocationX[alpha]) ) ;
      returnMatrix ( alpha, beta )  =  sum ;
      } // end for beta loop
  return  returnMatrix ;
  } // end operatorMatrix function

// ================================================================================================
