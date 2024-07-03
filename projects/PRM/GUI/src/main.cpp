#include "main_gui.hpp"

int main( int argc, char *argv[] )
{
    QApplication app( argc, argv );
    
    Moris_Gui widget;
    widget.show();

    return app.exec();
}