#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include <QPushButton>

#include "moris_line_edit.hpp"
#include "moris_combo_box.hpp"
#include "main_gui.hpp"
#include "TestWindow.hpp"


#include "main.moc"

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

int main(int argc, char *argv[])
{
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    QApplication app(argc, argv);

    TestWindow mainWindow;
    mainWindow.show();

    return app.exec();
}


//Moris_Gui moris_gui;
//moris_gui.show();