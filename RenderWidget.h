/*
  RenderWidget.h    Jeffery Solheim
  Portion of a Qt application illustrating use of the BSCM C++ library.
*/

#ifndef RENDERWIDGET_H
#define RENDERWIDGET_H

#include <QWidget>
#include <QPainter>
#include <QPen>
#include <vector>
#include <cmath>

//#include "MainWindow.h"

#include "Spline.h"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <unsupported/Eigen/MatrixFunctions>

class RenderWidget : public QWidget
  {
  Q_OBJECT // must include using Qt signals/slots

  private:
    int argc ;
    char **argv ;

  public:  //  ====================  constructor  ====================
    double f ( double x )
      {
        return  ( 100.0 * sin( 4.0 * atan(1.0) * x / 800.0 ) ) ;
      } // end f

    RenderWidget ( int argc, char **argv, QWidget *parent = 0 )
        : QWidget (parent)
      {
      order_M       =   3 ;
      num_Coll_Pts  =   3 ;
      left_Knot     = 100.0 ;
      right_Knot    = 800.0 ;
      this->argc  =  argc ;
      this->argv  =  argv ;
      } // end RenderWidget constructor

  public slots :
    void order_M_Changed ( int new_M )
      {
      order_M  =  new_M ;
      this->repaint();
      } // end order_M_Changed
    void num_Coll_Pts_Changed ( int new_N )
      {
      num_Coll_Pts  =  new_N ;
      this->repaint();
      } // end num_Coll_Pts_Changed
    void left_Knot_Changed ( const QString & new_Left_Knot )
      {
      left_Knot  =  new_Left_Knot.toDouble() ;
      this->repaint();
      } // end left_Knot_Changed
    void right_Knot_Changed ( const QString & new_Right_Knot )
      {
      right_Knot  =  new_Right_Knot.toDouble() ;
      this->repaint();
      } // end left_Knot_Changed

  protected :
    void paintEvent ( QPaintEvent * pep )
      {
      //{ // ---------------------------------------------------------------------------------
      const unsigned int         SPLINE_ORDER  =  3 ;
      const double               KNOT_ARRAY []  =  { 100, 200, 300, 400, 500, 600, 700, 800 } ;
      const std::vector<double>  KNOT_VECTOR ( KNOT_ARRAY, KNOT_ARRAY + sizeof(KNOT_ARRAY) / sizeof(double) ) ;

      Eigen::MatrixXi  constraintMatrix ( 2, 3 ) ;
      constraintMatrix  <<  1, 0, 0,   //  0th derivative at left  boundary to be zero
                            1, 0, 0 ;  //  0th derivative at right boundary to be zero
      //} // ---------------------------------------------------------------------------------

      // Spline ( size_t order, std::vector<double> knotX, Eigen::MatrixXi K_matrix ) ;
      unsigned int num_Knots  =  (num_Coll_Pts + 2*order_M - 1) ;
      std::vector<double> knotX ;
      double  length  =  right_Knot - left_Knot ;
      knotX.push_back ( left_Knot ) ;
      double  delta_x  =  length / (num_Knots - 1) ;
      double  x        =  left_Knot + delta_x ;
      for ( unsigned int i = 1 ; i <= (num_Coll_Pts + 2*(order_M) - 3) ; i ++ )
        {
        knotX.push_back ( x );
        x  +=  delta_x ;
        }
      knotX.push_back ( right_Knot ) ;
      // instantiate the Spline class ...
      BSCM::Spline  test_Spline ( order_M, knotX, constraintMatrix ) ;
      Eigen::MatrixXd  opMtrx  =  test_Spline.operatorMatrix(2) ;
      Eigen::MatrixXd  E_to_the_D  =  opMtrx.exp() ;
      Eigen::Vector3d  f_alpha ;
      f_alpha <<  f(test_Spline.collocationX[0]),
                  f(test_Spline.collocationX[1]),
                  f(test_Spline.collocationX[2]) ;


/*
      A ←  e^[D∙α^2 ]
      f ←  f_α  (initial conditions u(x,0))
      for t = 1, 2, 3, 4, 5, …
          f←  A ∙f   (new f obtained by multiplying preceding f by A)
*/

      QPainter  painter ( this ) ;
      QPen  blackPen ( QColor(0,0,0)) ;
      QPen  redPen ( QColor(255,0,0)) ;
      QPen  bluePen ( QColor(0,0,255)) ;
      redPen.setWidth(10);
      bluePen.setWidth(5);
      painter.setPen(blackPen);
      painter.drawLine(0,this->height()/2,this->width(),this->height()/2);
      painter.setPen(redPen);
      for
        (
          std::vector<double>::iterator it = knotX.begin() ;
          it != knotX.end() ;
          ++it
        )
          {
          double  x  =  *it ;
          painter.drawPoint ( x, this->height()/2 ) ;
          }// end for loop

      for ( int x = 0 ; x <= 800 ; x ++ )
        painter.drawPoint( x, (this->height()/2) + floor(f(x)) );
      painter.setPen(bluePen);
      for ( int alpha = 0 ; alpha < test_Spline.collocationX.size() ; alpha ++ )
        {
        double x = test_Spline.collocationX[alpha] ;
        double y = f(x) ;
        painter.drawPoint( x, (this->height()/2) + floor(y) );
        }

/*
      for ( int t = 0 ; t < 3 ; t ++ )
        {
        f_alpha  =  E_to_the_D * f_alpha ;
        for ( int alpha = 0 ; alpha < test_Spline.collocationX.size() ; alpha ++ )
          {
          double x = test_Spline.collocationX[alpha] ;
          double y = f(x) ;
          painter.drawPoint( floor(x), floor(y) );
          }
        } // end for loop
*/


      /*
      QPainter  painter ( this ) ;
      QPen  myPen ( QColor(150,150,0)) ;
      painter.setPen(myPen);
      //painter.setWindow(0,-3,9,6); // x varies 0 to 9, y varies -3 to +3
      //((MainWindow *)(((RenderBox*)parent)->parent))->controlBox->leftKnotEdit->value();

      for( int i = 0 ; i <= 999 ; i+=1 )
          for ( int j = 0 ; j <= 679 ; j+=1 )
              {
              myPen.setColor( QColor( (150+10*i)%256, (150+(j+2)*10)%256, 0 ) );
              painter.setPen(myPen);
              painter.drawPoint(i,j);
              }
*/
/*
      QPainter  painter ( this ) ;
      painter.translate( QString(argv[1]).toInt(), QString(argv[2]).toInt() );
      painter.rotate ( QString(argv[3]).toInt() ) ;
      QPen  myPen ( QColor(150,150,0)) ;
      //myPen.setWidthF(0.25);
      painter.setPen(myPen);
      QPen  redPen ( QColor(255,0,0) ) ;
      //redPen.setWidthF( 0.01 * length );
      //redPen.setWidth(1); // in "pixels"
      QPen  greenPen ( QColor(0,255,0) ) ;
      //greenPen.setWidthF( 0.01 * length );
      //greenPen.setWidth(1); // in "pixels"
      //painter.setViewport(0,0,1000,680);
      painter.setWindow(0,-3,9,6); // x varies 0 to 9, y varies -3 to +3
//      painter.drawPoint(0,0);
      for( int i = 0 ; i <= 9 ; i+=1 )
          for ( int j = -3 ; j <= +3 ; j+=1 )
              {
//              if ( j < 0 )
//                painter.setPen(redPen);
//              if ( j > 0 )
//                painter.setPen ( greenPen ) ;
              myPen.setColor( QColor( 150+10*i, 150+(j+2)*10, 0 ) );
              painter.setPen(myPen);
              painter.drawPoint(i,j);
              }
//      myPen.setColor ( QColor(200,150,100)) ;
//      painter.drawPoint(8,2);
*/
/*
      painter.setWindow ( floor(left_Knot),   0,
                          ceil(length),     100 ) ;
      // if window measures (x,y,width,height) = (1,0,7,100)
      // x varies from 1 to 8
      // y varies from 0 to 100
      //painter.drawPoint(floor(left_Knot), 50) ;
      for ( int i = 2 ; i <= 4 ; i ++ )
          for ( int j = 50 ; j <= 52 ; j ++ )
              painter.drawPoint(i,j);
//      painter.drawLine ( floor(left_Knot), 50,
//                         ceil(right_Knot), 50 ) ;
 **/
/*
unsigned int  order_M ;
unsigned int  num_Coll_Pts ;
double        left_Knot ;
double        right_Knot ;
      // y = 200 + 100 * sin ( x / 100 )
      painter.setPen ( QColor(0,0,0) ) ;

      // draw x-axis ...
      painter.drawLine ( 0, 200, 628, 200 ) ;
      for ( int x = 0 ; x <= 628 ; x += (628 / 4) )
        painter.drawLine ( x, 200-20, x, 200+20 ) ;

      // draw several periods of sine wave
      painter.setPen ( QColor(255,0,0) ) ;
      for ( int x1 = 0 ; x1 <= (628 - 4) ; x1 += 4 )
        {
        int x2 = x1 + 4 ;
        int y1 = 200 + sineAmplitude * sin ( (x1 * numberPeriods) / 100.0 ) ;
        int y2 = 200 + sineAmplitude * sin ( (x2 * numberPeriods) / 100.0 ) ;
        painter.drawLine ( x1, y1, x2, y2 ) ;
        } // end for x loop
*/
      } // end paintEvent

  private:  //  ================  instance variables  ================
    unsigned int  order_M ;
    unsigned int  num_Coll_Pts ;
    double        left_Knot ;
    double        right_Knot ;

  }; // end RenderWidget class

#endif // RENDERWIDGET_H
