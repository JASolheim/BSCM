/*
  mainQt.cpp    Jeffery Solheim
  Portion of a Qt application illustrating use of the BSCM C++ library.
*/

#include <QApplication>
#include "MainWindow.h"

#include <iostream>

int main ( int argc, char *argv[] )
    {
    QApplication  app ( argc, argv ) ;
    for ( int i = 0 ; i < argc ; i ++ )
        std::cout << argv[i] << std::endl ;
    app.setApplicationName ( "Heat Equation via BSCM" ) ;
    MainWindow  mainWin ( argc, argv )  ;
    mainWin.resize ( 1000, 680 ) ;
    mainWin.show ( ) ;
    return  app.exec ( ) ;
    } // end main
