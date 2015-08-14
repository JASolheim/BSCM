/*
  RenderBox.h    Jeffery Solheim
  Portion of a Qt application illustrating use of the BSCM C++ library.
*/

#ifndef RENDERBOX_H
#define RENDERBOX_H

#include <QGroupBox>
#include <QVBoxLayout>
#include "RenderWidget.h"

class RenderBox : public QGroupBox // inherits from QWidget
  {

  public:  //  ================  instance variables  =================
    RenderWidget *  renderWidget ;

  public:  //  ====================  constructor  ====================
    RenderBox ( int argc, char **argv, const QString & title, QWidget *parent = 0 )
        : QGroupBox(title,parent)
      {
      this->setStyleSheet
       ("QGroupBox { margin:5px; border:5px solid gray; border-radius:5px; padding:5px; }");
      QVBoxLayout * layout  =  new QVBoxLayout ;
      this->renderWidget    =  new RenderWidget ( argc, argv, this ) ;
      layout->addWidget ( this->renderWidget ) ;
      this->setLayout ( layout ) ;
      ((QMainWindow *)parent)->setCentralWidget ( this ) ;
      } // end RenderBox constructor

  }; // end RenderBox class

#endif // RENDERBOX_H
