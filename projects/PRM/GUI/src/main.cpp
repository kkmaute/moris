#include "main_gui.hpp"

int main( int argc, char *argv[] )
{
    QApplication app( argc, argv );

    // to test access to parameter lists
    // moris::Parameter_List tList = moris::prm::create_property_parameter_list();

    // for ( auto it = tList.begin(); it != tList.end(); ++it )
    // {
    //     std::cout << it->first << "\n";
    // }

    // simple window
    Moris_Gui widget;
    widget.show();

    return app.exec();
}