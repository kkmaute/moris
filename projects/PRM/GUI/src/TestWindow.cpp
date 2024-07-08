#include <QVBoxLayout>
#include <QPushButton>
#include <QDebug>

#include "TestWindow.hpp"
#include "moris_line_edit.hpp"
#include "moris_combo_box.hpp"
#include "moris_double_spin_box.hpp"
#include "moris_int_spin_box.hpp"

TestWindow::TestWindow( QWidget *parent )
        : QMainWindow( parent )
{
    auto *centralWidget = new QWidget( this );
    auto *layout        = new QVBoxLayout( centralWidget );

    // Create line edits and add them to the layout
    for ( int i = 1; i < 6; ++i )
    {
        auto *lineEdit = new Moris_Line_Edit( this );
        lineEdit->setObjectName( QString( "lineEdit%1" ).arg( i ) );
        layout->addWidget( lineEdit );

        // Connect lineEdit's textChanged signal to updateInput slot
        connect( lineEdit, &Moris_Line_Edit::textChanged, this, &TestWindow::updateInput );

        // Store lineEdit in lineEdits map
        lineEdits.insert( lineEdit->objectName(), lineEdit );

        // Insert initial empty parameters into parameterList
        parameterList.insert( lineEdit->objectName().toStdString(), "" );
    }

    // Create combo box and add it to the layout
    comboBox = new Moris_Combo_Box( this );
    comboBox->setObjectName( "comboBox1" );
    comboBox->addItem( "Option 1" );
    comboBox->addItem( "Option 2" );
    comboBox->addItem( "Option 3" );
    layout->addWidget( comboBox );

    // Connect comboBox's currentIndexChanged signal to updateComboBox slot
    connect( comboBox, &Moris_Combo_Box::currentIndexChanged, this, &TestWindow::updateComboBox );

    // Insert combo box into parameterList with an empty initial value
    parameterList.insert( comboBox->objectName().toStdString(), "" );

    // Create double spin box and add it to the layout
    doubleSpinBox = new Moris_Double_Spin_Box( this );
    doubleSpinBox->setObjectName( "doubleSpinBox1" );
    layout->addWidget( doubleSpinBox );

    // Connect doubleSpinBox's valueChanged signal to updateDoubleSpinBox slot
    connect( doubleSpinBox, QOverload< const QString &, const QVariant & >::of( &Moris_Double_Spin_Box::valueChanged ), this, &TestWindow::updateDoubleSpinBox );

    // Insert double spin box into parameterList with an empty initial value
    parameterList.insert( doubleSpinBox->objectName().toStdString(), "" );

    // Create integer spin box and add it to the layout
    intSpinBox = new Moris_Int_Spin_Box( this );
    intSpinBox->setObjectName( "intSpinBox1" );
    layout->addWidget( intSpinBox );

    // Connect intSpinBox's valueChanged signal to updateIntSpinBox slot
    connect( intSpinBox, QOverload< const QString &, const QVariant & >::of( &Moris_Int_Spin_Box::valueChanged ), this, &TestWindow::updateIntSpinBox );

    // Insert integer spin box into parameterList with an empty initial value
    parameterList.insert( intSpinBox->objectName().toStdString(), "" );

    // Create button to save and print inputs
    auto *printButton = new QPushButton( "Save and Print Inputs", this );
    layout->addWidget( printButton );

    // Connect button's clicked signal to saveAndPrintInputs slot
    connect( printButton, &QPushButton::clicked, this, &TestWindow::saveAndPrintInputs );

    // Set central widget for the main window
    setCentralWidget( centralWidget );
}

// Slot to update currentInputs with line edit text changes
void TestWindow::updateInput( const QString &name, const QString &text )
{
    currentInputs.setParameter( name, text );    // Set parameter in currentInputs
}

// Slot to update currentInputs with combo box selection changes
void TestWindow::updateComboBox( const QString &name, int index )
{
    QString selectedItem = comboBox->itemText( index );    // Get selected item text
    currentInputs.setParameter( name, selectedItem );      // Set parameter in currentInputs
}

// Slot to update currentInputs with double spin box value changes
void TestWindow::updateDoubleSpinBox( const QString &name, const QVariant &value )
{
    currentInputs.setParameter( name, value );    // Set parameter in currentInputs
}

// Slot to update currentInputs with integer spin box value changes
void TestWindow::updateIntSpinBox( const QString &name, const QVariant &value )
{
    currentInputs.setParameter( name, value );    // Set parameter in currentInputs
}

// Slot to save inputs from currentInputs to parameterList and print them
void TestWindow::saveAndPrintInputs()
{
    bool allInputsAdded = true;

    // Transfer inputs from currentInputs to parameterList
    for ( const auto &key : currentInputs.allParameters().keys() )
    {
        QVariant variantValue = currentInputs.parameter( key );

        // Convert QVariant to the appropriate type
        QString stringValue;
        if ( variantValue.canConvert< QString >() )
        {
            stringValue = variantValue.toString();
        }
        else if ( variantValue.canConvert< int >() )
        {
            int intValue = variantValue.toInt();
            stringValue  = QString::number( intValue );
        }
        else if ( variantValue.canConvert< double >() )
        {
            double doubleValue = variantValue.toDouble();
            stringValue        = QString::number( doubleValue, 'g', 15 );    // 'g' format for double
        }
        else
        {
            qDebug() << "Unsupported data type for parameter" << key;
            continue;    // Skip unsupported types
        }

        // Set the parameter in parameterList, unlocking after setting
        try
        {
            parameterList.set( key.toStdString(), stringValue.toStdString(), false );
        } catch ( const std::runtime_error &e )
        {
            qDebug() << "Error setting parameter" << key << ": " << e.what();
            allInputsAdded = false;
        }
    }

    // Print all parameters in parameterList
    qDebug() << "Saved Input Set:";
    for ( const auto &param : parameterList )
    {
        std::string paramName  = param.first;
        std::string paramValue = param.second.get_string();

        // Remove surrounding double quotes if they exist
        if ( !paramValue.empty() && paramValue.front() == '"' && paramValue.back() == '"' )
        {
            paramValue = paramValue.substr( 1, paramValue.size() - 2 );
        }

        qDebug() << QString::fromStdString( paramName ) << ": " << QString::fromStdString( paramValue );
    }

    // Check if all inputs from currentInputs have been successfully allocated to parameterList
    for ( const auto &key : currentInputs.allParameters().keys() )
    {
        if ( !parameterList.exists( key.toStdString() ) )
        {
            qDebug() << "Input " << key << " is missing in parameterList.";
            allInputsAdded = false;
        }
    }

    if ( allInputsAdded )
    {
        qDebug() << "All inputs have been successfully allocated to parameterList.";
    }
    else
    {
        qDebug() << "Inputs are missing in parameterList.";
    }
}
