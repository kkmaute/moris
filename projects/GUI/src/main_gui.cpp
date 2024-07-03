#include <QApplication>
#include <QLineEdit>
#include <QPushButton>
#include <QVBoxLayout>
#include <QWidget>
#include <iostream>

#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Logger.hpp"                   // MRS/IOS/src

#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

class MyWidget : public QWidget
{
    Q_OBJECT

  public:
    MyWidget( QWidget *parent = nullptr )
            : QWidget( parent )
    {
        // Create the QLineEdit and QPushButton
        lineEdit            = new QLineEdit( this );
        QPushButton *button = new QPushButton( "Store Name", this );

        // Set up the layout
        QVBoxLayout *layout = new QVBoxLayout( this );
        layout->addWidget( lineEdit );
        layout->addWidget( button );
        setLayout( layout );

        // Connect the button's clicked signal to the appropriate slot
        connect( button, &QPushButton::clicked, this, &MyWidget::storeName );
    }

  public slots:
    void storeName()
    {
        // Store the name and print it
        QString name = lineEdit->text();
        std::cout << "Name stored: " << name.toStdString() << std::endl;
    }

  private:
    QLineEdit *lineEdit;
};

int main( int argc, char *argv[] )
{
    QApplication app( argc, argv );

    // to test access to parameter lists
    moris::Parameter_List tList = moris::prm::create_property_parameter_list();

    for ( auto it = tList.begin(); it != tList.end(); ++it )
    {
        std::cout << it->first << "\n";
    }

    // simple window
    MyWidget widget;
    widget.show();

    return app.exec();
}

#include "main_gui.moc"
