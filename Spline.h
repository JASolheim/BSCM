/**
 * @file    Spline.h
 * @author  Jeff Solheim <JASolheim@FHSU.edu>
 * @version  1.0
 *
 * @section LICENSE
 * This program is distributed WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * @section DESCRIPTION
 * File Spline.h contains the declaration of the Spline class
 * of the Basis Spline Collocation Method (BSCM).
 */

#ifndef  SPLINE_H
#define  SPLINE_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/LU>

namespace BSCM
  {

  /**
   * @brief
   * Class %Spline implements the <em>Basis %Spline Collocation Method</em> (BSCM).
   *
   * See the paper entitled <b>Basis %Spline Collocation Method for the Lattice Solution
   * of Boundary Value Problems</b>, by A.S. Umar, J. Wu, M.R. Strayer, & C. Bottcher.
   * This paper appeared in the <em>Journal of Computational Physics</em>, 93, 426-448, 1991,
   * and, in the following documentation, is referred to as "Umar".
   */

  class  Spline
    {

    public :  //  ----------------------------------  Data Members  ----------------------------------------------

      /**
        * @brief
        * <b><em>Order&nbsp; M</em></b> &nbsp;of interpolating spline function
        *
        * Each basis spline function &nbsp;<b><em>B<sub>&nbsp;i</sub><sup>M</sup></em></b>&nbsp;
        * is defined piecewise via polynomials of degree <b><em>M</em></b> &minus; 1.\n\n
        * <b><em>M</em></b> &nbsp;must be odd so that <b><em>K_matrix</em></b> will have an even
        * number of rows.\n
        * See Umar's p. 427.
        */
      size_t  order ;

      /**
        * @brief
        * <b><em>Sequence of knot points</em></b>&nbsp;&nbsp;
        * {&nbsp;<b><em>x<sub>&nbsp;i</sub></em></b>&nbsp;}&nbsp;&nbsp; along the x-axis
        * 
        *
        * Uses zero-based indexing, <b><em>in contrast to Umar's one-based indexing</em></b>.\n
        * Therefore, the knots themselves are&nbsp; <em>x</em><sub>&nbsp;0</sub> .. 
        * <em>x</em><sub>&nbsp;<b><em>N</em></b>+2<b><em>M</em></b>-2</sub> .
        * See p. 427 of Umar et al.
        */
      std::vector<double>  knotX ;

      /**
        * @brief
        * <b><em>Number of knot points</em></b> along the <em>x</em>-axis
        *
        * Denoted <b><em>N</em></b> + 2<b><em>M</em></b> &minus; 1 by Umar et al (p. 430).
        */
      size_t  numKnots ;

      /**
        * @brief
        * <b><em>Sequence of collocation points</em></b>
        * &nbsp;&nbsp;{&nbsp;<b><em>x<sub>&nbsp;&alpha;</sub></em></b>&nbsp;}&nbsp;&nbsp;
        * within the physical boundaries along the x-axis
        *
        * Following Umar's suggestion (Equation (13), p. 431), each collocation point is located
        * midway between two knots.\n
        * Uses zero-based indexing, <b><em>in contrast to Umar's one-based indexing</em></b>.\n
        */
      std::vector<double>  collocationX ;

      /**
        * @brief
        * <b><em>Number of collocation points</em></b> along the <em>x</em>-axis
        *
        * Denoted <b><em>N</em></b> Umar et al (p. 431, Equation (13)).
        */
      size_t  N ;

      /**
        * @brief
        * The matrix&nbsp; <b><em>B<sub>&nbsp;&alpha;&nbsp;i</sub></em></b>
        * &nbsp;of Umar's Equation (14), p. 431
        *
        * <em><b>B<sub>&nbsp;&alpha;&nbsp;i</sub></b></em> &nbsp;is defined as&nbsp;
        * <em><b>B<sub>&nbsp;&alpha;&nbsp;i</sub></b></em> &nbsp;=&nbsp; 
        * <em><b>B<sub>&nbsp;i</sub><sup>M</sup>&nbsp;(&nbsp;x<sub>&alpha;&nbsp;</sub>)</b></em>.\n
        * @par  Subscript <b><em>&nbsp;&alpha;&nbsp;</em></b>
        *   Subscript <b><em>&nbsp;&alpha;&nbsp;</em></b> varies over collocation points,
        *   0 .. (<b><em>N</em></b>&minus;1).\n
        *   <b><em>This differs from Umar's paper</em></b>,
        *   in which <b><em>&nbsp;&alpha;&nbsp;</em></b> varies 1 .. <b><em>N</em></b>.\n
        * @par  Subscript <b><em>&nbsp;i&nbsp;</em></b>
        *   Subscript <b><em>&nbsp;i&nbsp;</em></b> is an index selecting basis function
        *   &nbsp;<em><b>B<sub>&nbsp;i</sub><sup>M</sup>&nbsp;</b></em>.\n
        */
      Eigen::MatrixXd  B_matrix ;

      /**
        * @brief
        * The matrix&nbsp; <b><em>K<sub>&nbsp;r&nbsp;p</sub></em></b>
        * &nbsp;of Umar's Equation (16), p. 432
        *
        *  Determines which linear combinations of derivatives of
        * <b><em>&nbsp;f&nbsp;</em></b> shall be forced to be zero.\n
        * @par  Subscript <b><em>&nbsp;r&nbsp;</em></b>
        *   Subscript <b><em>&nbsp;r&nbsp;</em></b> enumerates the fixed boundary conditions
        *   and varies 0 .. (<b><em>M</em></b>&minus;2).\n
        *   <b><em>This differs from Umar's paper</em></b>, in which <b><em>&nbsp;r&nbsp;</em></b> varies
        *   1 .. (<b><em>M</em></b>&minus;1).\n\n
        *   In the first (<em><b>M</b></em>&minus;1)/2 rows of&nbsp;
        *   <b><em>K<sub>&nbsp;r&nbsp;p</sub></em></b>,
        *   derivatives are evaluated at the left boundary,\n
        *   while in the last (<em><b>M</b></em>&minus;1)/2 rows of&nbsp;
        *   <b><em>K<sub>&nbsp;r&nbsp;p</sub></em></b>,
        *   derivatives are evaluated at the right boundary.\n
        * @par  Subscript <b><em>&nbsp;p&nbsp;</em></b>
        *   Subscript <b><em>&nbsp;p&nbsp;</em></b> enumerates derivatives; <b><em>p</em></b> varies 0 .. (<b><em>M</em></b>&minus;1).\n
        */
      Eigen::MatrixXi K_matrix ;

      /**
        * @brief
        * The matrix&nbsp; <b><em>&beta;<sub>&nbsp;r&nbsp;i</sub></em></b>
        * &nbsp;of Umar's Equation (18), p. 432
        * @par  Subscript <b><em>&nbsp;r&nbsp;</em></b>
        *   Subscript <b><em>&nbsp;r&nbsp;</em></b> varies 0 .. (<b><em>M</em></b>&minus;2).\n
        *   <b><em>This differs from Umar's paper</em></b>, in which <b><em>&nbsp;r&nbsp;</em></b> varies
        *   1 .. (<b><em>M</em></b>&minus;1).\n
        * @par  Subscript <b><em>&nbsp;i&nbsp;</em></b>
        *   Subscript <b><em>&nbsp;i&nbsp;</em></b> varies 0 ..
        *   (&nbsp;<b><em><b><em>N</em></b></em></b> + <em><b>M</b></em> &minus; 2&nbsp;).\n
        *   <b><em>This differs from Umar's paper</em></b>, in which <b><em>&nbsp;i&nbsp;</em></b> varies
        *   1 .. (&nbsp;<em><b><b><em>N</em></b></b></em> + <b><em>M</em></b> &minus; 1&nbsp;).\n
        */
      Eigen::MatrixXd  beta_matrix ;

    public :  //  ----------------------------------  Member Functions  ------------------------------------------

      /**
        * @brief
        * Construct a BSCM %Spline object
        *
        * @param order    %Spline order (denoted <em><b>M</b></em> by Umar et al, p. 427)
        * @param knotX    Collection of knots &nbsp;<em><b>x<sub>0&nbsp;</sub></b></em>,
        *                 <em><b>x<sub>1&nbsp;</sub></b></em>, <em><b>x<sub>2&nbsp;</sub></b></em>, ...
        * @param K_matrix Specifies fixed boundary conditions
        *                 (denoted&nbsp; <em><b>K<sub>&nbsp;r&nbsp;p</sub></b></em>
        *                 &nbsp;in Umar's Equation (16), p. 432)
        */
      Spline ( size_t order, std::vector<double> knotX, Eigen::MatrixXi K_matrix ) ;
    
      /**
        * @brief
        * <em><b>B<sub>&nbsp;i</sub><sup>k</sup>&nbsp;(x)</b></em>,
        * &nbsp;basis function &nbsp;<em><b>B<sub>&nbsp;i</sub><sup>k</sup></b></em>,
        * &nbsp;evaluated at <b><em>&nbsp;x&nbsp;</em></b>
        *
        * <em>B<sub>&nbsp;i</sub><sup>k</sup>(x)</em> is defined only when <em>x</em>
        * is between the first &amp; last knots.\n
        * Also, <em>B<sub>i</sub><sup>k</sup>(x)</em> &ne; 0 only when <em>x</em> is
        * between knots <em>x<sub>i</sub></em> &amp; <em>x<sub>i+k</sub></em>
        * (i.e., when <em>x<sub>i</sub></em> &le; <em>x</em> &lt; <em>x<sub>i+k</sub></em>).\n
        * See Umar's Equations (1), (2), and (3), p. 428.
        * @param  k  %Spline order, ranging from 1, ..., <b><em>M</em></b>
        * @param  i  %Spline index, ranging 0, ..., (<em>numKnots</em> &minus; <em>k</em> &minus; 1).
        * @param  x  Location along horizontal axis
        * @return  double
        * @par     Note
        *          This version of&nbsp; <b><em>B</em></b> &nbsp;has third parameter of type double,
        *          whereas another version of overloaded&nbsp; <b><em>B</em></b> &nbsp;has
        *          third parameter of type size_t.
        */
      double B ( size_t k, size_t i, double x ) ;

      /**
        * @brief
        * <em><b>B<sub>&nbsp;i</sub><sup>k</sup>&nbsp;
        * (&nbsp;x<sub>&alpha;&nbsp;</sub>)</b></em>,
        * &nbsp;basis function <em><b>B<sub>&nbsp;i</sub><sup>k</sup>&nbsp;</b></em>
        * &nbsp;evaluated at a collocation point <b><em>&nbsp;x<sub>&nbsp;&alpha;&nbsp;</sub></em></b>
        * 
        *
        * @param  k     %Spline order, ranging from 1, ..., <b><em>M</em></b>
        * @param  i     %Spline index, ranging 0, ..., (<em>numKnots</em> &minus; <em>k</em> &minus; 1).
        * @param  alpha Index specifying collocation point
        *         <b><em>&nbsp;x<sub>&nbsp;&alpha;</sub></em></b>\n
        *         <b><em>&alpha;</em></b> &nbsp;is an index into
        *         <b><em>collocationX</em></b>, the vector of collocation points\n
        *         <b><em>&alpha;</em></b> &nbsp; ranges 0 .. (<b><em>N</em></b> &minus; 1)
        * @return  double
        * @par     Note
        *          This version of&nbsp; <b><em>B</em></b> &nbsp;has third parameter of type size_t,
        *          whereas another version of overloaded&nbsp; <b><em>B</em></b> &nbsp;has
        *          third parameter of type double.
        */
      double B ( size_t k, size_t i, size_t alpha ) ;

      /**
        * @brief
        * <b><em>&part;<sup>&nbsp;p</sup></em></b>
        * &nbsp;<em><b>B<sub>&nbsp;i</sub><sup>k</sup>&nbsp;(&nbsp;x&nbsp;)</b></em>,
        * &nbsp;the <em>&nbsp;<b>p<sup>&nbsp;th</sup></b>&nbsp;</em> derivative of
        * &nbsp;<em><b>B<sub>&nbsp;i</sub><sup>k</sup></b></em>, 
        * &nbsp;evaluated at&nbsp; <b><em>x</em></b>
        * 
        * See Umar's Equations (4), (5), (6), &amp; (7), pp. 428-429.
        * @param  p  Derivative order, ranging 0, ..., (<b><em>M</em></b>&minus;1)\n
        *            <em>p</em> = 0 is 0<sup>th</sup> derivative,&nbsp; 
        *            <em>p</em> = 1 is 1<sup>st</sup> derivative, etc.
        * @param  k  %Spline order, ranging (<em>p</em>+1), ..., <b><em>M</em></b>\n
        * @param  i  %Spline index, ranging 0, ...,
        *            (<b><em>numKnots</em></b> &minus; <b><em>k</em></b> &minus; 1).
        * @param  x  Location along horizontal axis
        * @return  double
        * @par     Note
        *          If <em>p</em> &ge; <em>k</em>, &nbsp;
        *          <em>&part;<sup>&nbsp;p</sup></em> &nbsp;<em>B<sub>&nbsp;i</sub><sup>k</sup>&nbsp;(&nbsp;x&nbsp;)</em> would be 0
        *          (where it is defined) since
        *          &nbsp;<em>B<sub>&nbsp;i</sub><sup>k</sup>&nbsp;(&nbsp;x&nbsp;)</em>
        *          &nbsp; is comprised (piecewise) of polynomials of degree <em>k</em> &minus; 1.
        */
      double D_B ( size_t p, size_t k, size_t i, double x ) ;

      /**
        * @brief Matrix representation of differentiation operator
        *
        * Denoted <b><em>O<sub>&nbsp;&alpha;</sub><sup>&nbsp;&beta;</sup></em></b> &nbsp;by Umar,
        * this is the collocation space matrix representation of a differentiation operator
        * as defined by Umar's Equation (28), p. 434.\n
        * <em>O<sub>&nbsp;&alpha;</sub><sup>&nbsp;&beta;</sup></em> &nbsp;is of dimensions
        * <b><em>N</em></b> by <b><em>N</em></b>, and row &amp; column indices are
        * zero-based.
        * @param   derivativeOrder 1 indicates &part;/&part;<em>x</em>,&nbsp; 
        *          2 indicates &part;<sup>2</sup>/&part;<em>x</em><sup>2</sup>,&nbsp; etc.
        * @return  Eigen::MatrixXd
        * @par     Note
        *          Assumes that&nbsp; <em>f</em>(<em>N</em>), ...,
        *          <em>f</em>(<em>M</em>+<em>N</em>&minus;2) &nbsp;have all been set equal zero,
        *          as shown in Umar's Equation (21), p. 433.
        */
      Eigen::MatrixXd  operatorMatrix ( size_t derivativeOrder ) ;

    private :  //  -----------------------------------------------------------------------------------------------

      static const size_t MIN_ORDER      =   3 ;
      static const size_t MAX_ORDER      =  15 ;

      /**
        * This is the matrix "B tilde" of Umar, Equation (20), p. 433.
        */
      Eigen::MatrixXd  B_tilde_matrix ;

      /**
        * This is the matrix "C tilde" of Umar, Equation (22), p. 433.
        */
      Eigen::MatrixXd  C_tilde_matrix ;

      /**
        * This is the leftmost physical boundary; see Umar p 430.
        */
     double           xMin ;

      /**
        * This is the rightmost physical boundary; see Umar p 430.
        */
     double           xMax ;

      static const size_t MAX_NUMBER_KNOTS  =  100 ;

      /**
        * This function implements the recursion relation to find lower order derivatives
        * as shown in Umar p. 429, Equations (5, 6, and 7).
        * @return  double
        */
      double C ( size_t k, size_t i, double x ) ;

      /**
        * Vector storing values of B( size_t k, size_t i, size_t alpha )
        * B_k_i_alpha[k][i][alpha] = B ( k, i, alpha )
        */
      std::vector< std::vector< std::vector<double> > >  B_k_i_alpha ;

    } ; // end Spline class

  } // end namespace BSCM

#endif  //  SPLINE_H
