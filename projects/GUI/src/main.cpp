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

    moris::Moris_Gui widget;
    widget.show();

    // TestWindow mainWindow( nullptr, parameterList );
    // mainWindow.show();

    return app.exec();
}
