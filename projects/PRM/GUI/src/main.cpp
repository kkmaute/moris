#include "main_gui.hpp"


moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

int main( int argc, char *argv[] )
{
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    QApplication app( argc, argv );

    Moris_Gui widget;
    widget.show();

    return app.exec();
}