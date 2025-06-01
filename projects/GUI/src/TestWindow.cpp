#include "TestWindow.hpp"

namespace moris
{

    // Constructor for the TestWindow class
    TestWindow::TestWindow( QWidget *parent, moris::Parameter_List &parameterList )
            : QMainWindow( parent )
            , parameterList( parameterList )
    {
        // Create central widget and layout
        auto *centralWidget = new QWidget( this );
        auto *layout        = new QVBoxLayout( centralWidget );

        // Iterate over each parameter in the parameter list to initialize GUI elements
        for ( auto paramEntry : parameterList )
        {
            const auto &paramName = paramEntry.get_name();

            // Initialize Moris_Line_Edit if the parameter name matches
            if ( paramName == "lineEdit" )
            {
                lineEdit = new Moris_Line_Edit( this, paramEntry.get_parameter() );
                lineEdit->setObjectName( QString::fromStdString( paramName ) );
                lineEdit->setText( QString::fromStdString( paramEntry.get_parameter().get_value< std::string >() ) );
                layout->addWidget( lineEdit );
                connect( lineEdit, &Moris_Line_Edit::textChanged, this, &TestWindow::updateInput );
            }
            // Initialize Moris_Combo_Box if the parameter name matches
            else if ( paramName == "comboBox" )
            {
                comboBox = new Moris_Combo_Box( this, paramEntry.get_parameter() );
                comboBox->setObjectName( QString::fromStdString( paramName ) );
                comboBox->setCurrentIndex( paramEntry.get_parameter().get_value< moris::uint >() );
                layout->addWidget( comboBox );
                connect( comboBox, &Moris_Combo_Box::index_changed, this, &TestWindow::updateComboBox );
            }
            // Initialize Moris_Double_Spin_Box if the parameter name matches
            else if ( paramName == "doubleSpinBox" )
            {
                doubleSpinBox = new Moris_Double_Spin_Box( this, paramEntry.get_parameter() );
                doubleSpinBox->setObjectName( QString::fromStdString( paramName ) );
                doubleSpinBox->setValue( paramEntry.get_parameter().get_value< moris::real >() );
                layout->addWidget( doubleSpinBox );
                connect( doubleSpinBox, &Moris_Double_Spin_Box::value_changed, this, &TestWindow::updateDoubleSpinBox );
            }
            // Initialize Moris_Int_Spin_Box if the parameter name matches
            else if ( paramName == "intSpinBox" )
            {
                intSpinBox = new Moris_Int_Spin_Box( this, paramEntry.get_parameter() );
                intSpinBox->setObjectName( QString::fromStdString( paramName ) );
                intSpinBox->setValue( paramEntry.get_parameter().get_value< int >() );
                layout->addWidget( intSpinBox );
                connect( intSpinBox, &Moris_Int_Spin_Box::value_changed, this, &TestWindow::updateIntSpinBox );
            }
            // Initialize Moris_Pair_Box if the parameter name matches
            // else if ( paramName == "pairBox" )
            // {
            //     pairBox = new Moris_Pair_Box( this, paramEntry.get_parameter(), QStringList() << "Example 1" << "Example 2" << "Example 3" );
            //     pairBox->setObjectName( QString::fromStdString( paramName ) );
            //     layout->addWidget( pairBox );
            //     connect( pairBox, &Moris_Pair_Box::combo_box_text_changed, this, &TestWindow::updatePairBox );
            //     connect( pairBox, &Moris_Pair_Box::line_edit_text_changed, this, &TestWindow::updatePairBox );
            // }
        }

        // Create and add Save and Print Inputs button
        auto *printButton = new QPushButton( "Save and Print Inputs", this );
        layout->addWidget( printButton );
        connect( printButton, &QPushButton::clicked, this, &TestWindow::saveAndPrintInputs );

        // Set the central widget and layout
        centralWidget->setLayout( layout );
        setCentralWidget( centralWidget );
    }

    // Slot to update the text of Moris_Line_Edit
    void TestWindow::updateInput( const QString &name, const QString &text )
    {
        if ( lineEdit && lineEdit->objectName() == name )
        {
            lineEdit->setText( text );
        }
    }

    // Slot to update the index of Moris_Combo_Box
    void TestWindow::updateComboBox( const QString &name, int index )
    {
        if ( comboBox && comboBox->objectName() == name )
        {
            comboBox->setCurrentIndex( index );
        }
    }

    // Slot to update the value of Moris_Double_Spin_Box
    void TestWindow::updateDoubleSpinBox( const QString &name, double value )
    {
        if ( doubleSpinBox && doubleSpinBox->objectName() == name )
        {
            doubleSpinBox->setValue( value );
        }
    }

    // Slot to update the value of Moris_Int_Spin_Box
    void TestWindow::updateIntSpinBox( const QString &name, int value )
    {
        if ( intSpinBox && intSpinBox->objectName() == name )
        {
            intSpinBox->setValue( value );
        }
    }

    // Slot to update the Moris_Pair_Box values
    // void TestWindow::updatePairBox( const QString &name, const QString &text )
    // {
    //     if ( pairBox )
    //     {
    //         if ( pairBox->moris_pair_combo_box->objectName() == name )
    //         {
    //             pairBox->moris_pair_combo_box->setCurrentText( text );
    //         }
    //         else if ( pairBox->moris_pair_line_edit->objectName() == name )
    //         {
    //             pairBox->moris_pair_line_edit->setText( text );
    //         }
    //     }
    // }

    // Save all inputs to XML file
    void TestWindow::saveInputsToXML( const std::string &filePath )
    {
        // Initialize the XML parser with write mode
        moris::XML_Parser xmlParser( filePath, moris::XML_Mode::WRITE );

        // Set the root element for XML
        xmlParser.set( "UserInputs", "" );

        // Iterate over each parameter in the parameter list and save the values to XML
        for ( const auto &paramEntry : parameterList )
        {
            const auto &paramName = paramEntry.get_name();

            // Save value of Moris_Line_Edit
            if ( paramName == "lineEdit" )
            {
                if ( lineEdit )
                {
                    xmlParser.set( "UserInputs.LineEdit", lineEdit->text().toStdString() );
                }
            }
            // Save value of Moris_Combo_Box
            else if ( paramName == "comboBox" )
            {
                if ( comboBox )
                {
                    xmlParser.set( "UserInputs.ComboBox", comboBox->currentText().toStdString() );
                }
            }
            // Save value of Moris_Double_Spin_Box
            else if ( paramName == "doubleSpinBox" )
            {
                if ( doubleSpinBox )
                {
                    xmlParser.set( "UserInputs.DoubleSpinBox", std::to_string( doubleSpinBox->value() ) );
                }
            }
            // Save value of Moris_Int_Spin_Box
            else if ( paramName == "intSpinBox" )
            {
                if ( intSpinBox )
                {
                    xmlParser.set( "UserInputs.IntSpinBox", std::to_string( intSpinBox->value() ) );
                }
            }
            // Save combined values of Moris_Pair_Box
            else if ( paramName == "pairBox" )
            {
                // if ( pairBox )
                // {
                //     // Retrieve combined values from Moris_Pair_Box
                //     std::pair< std::string, std::string > pairValue;
                //     pairValue.first  = pairBox->moris_pair_combo_box->currentText().toStdString();
                //     pairValue.second = pairBox->moris_pair_line_edit->text().toStdString();
                //     xmlParser.set( "UserInputs.PairBox.LineEdit", pairValue.second );
                //     xmlParser.set( "UserInputs.PairBox.ComboBox", pairValue.first );
                // }
            }
        }

        // Save the XML file
        xmlParser.save();
        qDebug() << "Successfully saved inputs to XML.";
    }

    // Save and print all inputs
    void TestWindow::saveAndPrintInputs()
    {
        std::string xmlFilePath = "user_inputs.xml";

        // Save inputs to XML
        saveInputsToXML( xmlFilePath );

        // Print saved inputs to the terminal
        qDebug() << "Saved Input Set:";

        for ( const auto &paramEntry : parameterList )
        {
            const auto &paramName = paramEntry.get_name();

            // Print value of Moris_Line_Edit
            if ( paramName == "lineEdit" )
            {
                if ( lineEdit )
                {
                    qDebug() << QString::fromStdString( paramName ) << ": " << lineEdit->text();
                }
                else
                {
                    qDebug() << "LineEdit is null";
                }
            }
            // Print index of Moris_Combo_Box
            else if ( paramName == "comboBox" )
            {
                if ( comboBox )
                {
                    qDebug() << QString::fromStdString( paramName ) << ": " << comboBox->currentText();
                }
                else
                {
                    qDebug() << "ComboBox is null";
                }
            }
            // Print value of Moris_Double_Spin_Box
            else if ( paramName == "doubleSpinBox" )
            {
                if ( doubleSpinBox )
                {
                    qDebug() << QString::fromStdString( paramName ) << ": " << doubleSpinBox->value();
                }
                else
                {
                    qDebug() << "DoubleSpinBox is null";
                }
            }
            // Print value of Moris_Int_Spin_Box
            else if ( paramName == "intSpinBox" )
            {
                if ( intSpinBox )
                {
                    qDebug() << QString::fromStdString( paramName ) << ": " << intSpinBox->value();
                }
                else
                {
                    qDebug() << "IntSpinBox is null";
                }
            }
            // Print values of Moris_Pair_Box
            else if ( paramName == "pairBox" )
            {
                // if ( pairBox )
                // {
                //     qDebug() << QString::fromStdString( paramName ) << ": LineEdit:" << pairBox->moris_pair_line_edit->text() << ", ComboBox:" << pairBox->moris_pair_combo_box->currentText();
                // }
                // else
                // {
                //     qDebug() << "PairBox is null";
                // }
            }
        }
    }
}    // namespace moris