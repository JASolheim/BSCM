/*
  ControlBox.h    Jeffery Solheim
  Portion of a Qt application illustrating use of the BSCM C++ library.
*/

#ifndef CONTROLBOX_H
#define CONTROLBOX_H

#include <QGroupBox>
#include <QDockWidget>
#include <QGridLayout>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QSpinBox>
#include <QLineEdit>
#include <QCheckBox>
#include <QLabel>
#include <QSpacerItem>
#include <QPushButton>

class ControlBox : public QGroupBox // inherits from QWidget
  {

  public:  //  ================  instance variables  =================
    QSpinBox *     orderSpinBox ;
    QSpinBox *     nSpinBox ;
    QLineEdit *    leftKnotEdit ;
    QLineEdit *    rightKnotEdit ;
    QLineEdit *    thermDiffEdit ;
    QSpinBox *     initProfileSpinBox ;
    QPushButton *  goButton ;

  public:  //  ====================  constructor  ====================
    ControlBox ( const QString & title, QWidget *parent = 0 )
        : QGroupBox(title,parent)
      {
      orderSpinBox  =  new QSpinBox ( this ) ;
      orderSpinBox->setMinimum ( 3 ) ;
      orderSpinBox->setMaximum ( 7 ) ;
      orderSpinBox->setSingleStep ( 2 ) ;

      nSpinBox  =  new QSpinBox ( this ) ;
      nSpinBox->setMinimum (  3 ) ;
      nSpinBox->setMaximum ( 17 ) ;

      leftKnotEdit   =  new QLineEdit ( "100.0", this ) ;
      rightKnotEdit  =  new QLineEdit ( "800.0", this ) ;
      thermDiffEdit  =  new QLineEdit ( "1.0", this ) ;

      initProfileSpinBox  =  new QSpinBox ( this ) ;
      initProfileSpinBox->setMinimum ( 1 ) ;
      initProfileSpinBox->setMaximum ( 3 ) ;

      goButton  =  new QPushButton ( "Start Animation", this ) ;

      // boundary conditions specifier ...
      QGroupBox *  leftGB  =  new QGroupBox (" Left Boundary Conditions -- Zero Derivatives" );
      leftGB->setStyleSheet ( "QGroupBox { border:2px solid gray; border-radius:2px; }" ) ;
      QHBoxLayout *  leftGBLayout = new QHBoxLayout;
      leftGB->setLayout(leftGBLayout);

      for ( int p = 0 ; p < 7 ; p ++ )
        {
        QCheckBox * lcb = new QCheckBox(QString(p+'0'));
        leftGBLayout->addWidget(lcb);
        }

      QGroupBox *  rightGB  =  new QGroupBox (" Right Boundary Conditions -- Zero Derivatives");
      rightGB->setStyleSheet ( "QGroupBox { border:2px solid gray; border-radius:2px; }" ) ;
      QHBoxLayout *  rightGBLayout = new QHBoxLayout;
      rightGB->setLayout(rightGBLayout);

      for ( int p = 0 ; p < 7 ; p ++ )
        {
        QCheckBox * rcb = new QCheckBox(QString(p+'0'));
        rightGBLayout->addWidget(rcb);
        }

      QGridLayout *  layout = new QGridLayout ;
      layout->setColumnMinimumWidth(0,1);
      layout->setColumnMinimumWidth(1,1);
      layout->setColumnMinimumWidth(2,20);
      layout->setColumnMinimumWidth(3,1);
      layout->setColumnMinimumWidth(4,1);
      layout->setColumnMinimumWidth(5,20);
      layout->setColumnMinimumWidth(6,1);
      layout->setColumnMinimumWidth(7,20);
      layout->setColumnMinimumWidth(8,1);
      layout->setColumnMinimumWidth(9,1);
      layout->setColumnMinimumWidth(10,20);
      layout->setColumnMinimumWidth(11,1);
      layout->addWidget( new QLabel("Spline Order M: "),0,0);
      layout->addWidget( orderSpinBox,0,1);
      layout->addWidget( new QLabel("# Coll Pts N: "),1,0);
      layout->addWidget( nSpinBox,1,1);
      layout->addWidget( new QLabel("Leftmost Knot: "),0,3);
      layout->addWidget( leftKnotEdit,0,4);
      layout->addWidget( new QLabel("Rightmost Knot: "),1,3);
      layout->addWidget( rightKnotEdit,1,4);
      layout->addWidget(leftGB,0,6);
      layout->addWidget(rightGB,1,6);
      layout->addWidget( new QLabel("Thermal Diffusivity: "),0,8);
      layout->addWidget( thermDiffEdit,0,9);
      layout->addWidget( new QLabel("Initial Profile: "),1,8);
      layout->addWidget( initProfileSpinBox,1,9);
      layout->addWidget( goButton, 0, 11, 2, 1 );
      this->setLayout(layout);

      this->setStyleSheet
       ("QGroupBox { margin:5px; border:5px solid gray; border-radius:5px; padding:5px; }");
      QDockWidget *  dockWidget  =  new QDockWidget ( parent ) ;
      dockWidget->setAllowedAreas ( Qt::BottomDockWidgetArea ) ;
      dockWidget->setWidget ( this ) ;
      ((QMainWindow *)parent)->addDockWidget ( Qt::BottomDockWidgetArea, dockWidget ) ;
      } // end ControlBox constructor

  }; // end ControlBox class

#endif // CONTROLBOX_H
