/*
  MainWindow.h    Jeffery Solheim
  Portion of a Qt application illustrating use of the BSCM C++ library.
*/

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "RenderBox.h"
#include "ControlBox.h"

class MainWindow : public QMainWindow
  {

  public:  //  ====================  constructor  ====================
    MainWindow (  int argc, char **argv  )
      {
      renderBox   =  new RenderBox ( argc, argv, "  Heat Profile ", this ) ;
      controlBox  =  new ControlBox ( "  BSCM Parameters ", this ) ;

      // connect signals to slots ...
      this->connect
          (    controlBox->orderSpinBox ,
               SIGNAL ( valueChanged ( int ) ) ,
               renderBox->renderWidget ,
               SLOT   ( order_M_Changed ( int ) )    ) ;
      this->connect
          (    controlBox->nSpinBox ,
               SIGNAL ( valueChanged ( int ) ) ,
               renderBox->renderWidget ,
               SLOT   ( num_Coll_Pts_Changed ( int ) )    ) ; // change this line
      this->connect
          (    controlBox->leftKnotEdit ,
               SIGNAL ( textEdited ( const QString & ) ) ,
               renderBox->renderWidget ,
               SLOT   ( left_Knot_Changed ( const QString & ) )    ) ;
      this->connect
          (    controlBox->rightKnotEdit ,
               SIGNAL ( textEdited ( const QString & ) ) ,
               renderBox->renderWidget ,
               SLOT   ( right_Knot_Changed ( const QString & ) )    ) ;
      } ; // end MainWindow constructor

  public:  //  ================  instance variables  ================
    RenderBox *   renderBox ;
    ControlBox *  controlBox ;

  } ; // end MainWindow class

#endif // MAINWINDOW_H
