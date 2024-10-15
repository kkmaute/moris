#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include <QPushButton>

#include "cl_line_edit.hpp"
#include "cl_combo_box.hpp"
#include "main_gui.hpp"
#include "TestWindow.hpp"
#include "cl_Parameter_List.hpp"

#include "main.moc"

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

int main( int argc, char *argv[] )
{
    gMorisComm = moris::Comm_Manager( &argc, &argv );
    QApplication app( argc, argv );

    // Create a Parameter List and populate it with initial test values
    moris::Parameter_List parameterList;
    parameterList.insert( "lineEdit", "Initial text", {} );
    parameterList.insert_enum( "comboBox", { "Option 1", "Option 2", "Option 3" } );
    parameterList.insert( "doubleSpinBox", 0.5, {} );
    parameterList.insert( "intSpinBox", 10, {} );
    parameterList.insert_enum( "pairBox", { "Example 1", "Example 2", "Example 3" } );
    parameterList.insert_enum( "boolBox", { "True", "False" } );

    moris::Moris_Gui widget;
    widget.show();

    // TestWindow mainWindow( nullptr, parameterList );
    // mainWindow.show();

    return app.exec();
}
