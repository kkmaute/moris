#include "TestWindow.hpp"


TestWindow::TestWindow( QWidget *parent, moris::Parameter_List &parameterList )
        : QMainWindow( parent )
{
    // Create central widget and layout for the main window
    auto *centralWidget = new QWidget( this );
    auto *layout        = new QVBoxLayout( centralWidget );
    // std::cout << "here" << std::endl;

    // Loop through parameterList and initialize widgets
    for ( std::pair< std::string, moris::Parameter >  &paramEntry : parameterList )
    {
        const auto &paramName         = paramEntry.first;

        if ( paramName == "lineEdit" )
        {
            std::cout << "line" << std::endl;
            lineEdit = new Moris_Line_Edit( this, paramEntry.second );
            lineEdit->setObjectName( QString::fromStdString( paramName ) );
            lineEdit->setText( QString::fromStdString( paramEntry.second.get_value< std::string >(  ) ) );
            layout->addWidget( lineEdit );
            lineEdits.insert( lineEdit->objectName(), lineEdit );
            connect( lineEdit, &Moris_Line_Edit::textChanged, this, &TestWindow::updateInput );

            
        }
        else if ( paramName == "comboBox" )
        {
            std::cout << "comboBox" << std::endl;
            comboBox = new Moris_Combo_Box( this, paramEntry.second );
            comboBox->setObjectName( QString::fromStdString( paramName ) );
            comboBox->addItem( "Option 1" );
            comboBox->addItem( "Option 2" );
            comboBox->addItem( "Option 3" );
            comboBox->setCurrentIndex( paramEntry.second.get_value< moris::uint >(  ) );
            layout->addWidget( comboBox );
            connect( comboBox, &Moris_Combo_Box::indexChanged, this, &TestWindow::updateComboBox );
        }
        else if ( paramName == "doubleSpinBox" )
        {
            std::cout << "dspinbox" << std::endl;
            doubleSpinBox = new Moris_Double_Spin_Box( this,  paramEntry.second);
            doubleSpinBox->setObjectName( QString::fromStdString( paramName ) );
            doubleSpinBox->setValue( paramEntry.second.get_value< moris::real >( ) );
            layout->addWidget( doubleSpinBox );
            doubleSpinBoxes.insert( doubleSpinBox->objectName(), doubleSpinBox );
            connect( doubleSpinBox, &Moris_Double_Spin_Box::valueChanged, this, &TestWindow::updateDoubleSpinBox );
        }
        else if ( paramName == "intSpinBox" )
        {
            std::cout << "intSpinbox" << std::endl;
            intSpinBox = new Moris_Int_Spin_Box( this, paramEntry.second );
            intSpinBox->setObjectName( QString::fromStdString( paramName ) );
            intSpinBox->setValue( paramEntry.second.get_value< int >(  ) );
            layout->addWidget( intSpinBox );
            intSpinBoxes.insert( intSpinBox->objectName(), intSpinBox );
            connect( intSpinBox, &Moris_Int_Spin_Box::valueChanged, this, &TestWindow::updateIntSpinBox );
        }
        else if ( paramName == "pairBox" )
        {
            pairBox = new Moris_Pair_Box( this, paramEntry.second, QStringList() << "Example 1" << "Example 2" << "Example 3" );
            pairBox->setObjectName( QString::fromStdString( paramName ) );
            layout->addWidget( pairBox );
            pairBoxes.insert( pairBox->objectName(), pairBox );
            connect( pairBox, &Moris_Pair_Box::comboBoxTextChanged, this, &TestWindow::updatePairBox );
            connect( pairBox, &Moris_Pair_Box::lineEditTextChanged, this, &TestWindow::updatePairBox );
        }
    }

    // Create a button to save and print inputs
    auto *printButton = new QPushButton( "Save and Print Inputs", this );
    layout->addWidget( printButton );
    connect( printButton, &QPushButton::clicked, this, &TestWindow::saveAndPrintInputs );

    // Set the layout and central widget for the main window
    centralWidget->setLayout( layout );
    setCentralWidget( centralWidget );
}

void TestWindow::updateInput( const QString &name, const QString &text )
{
    auto lineEdit = lineEdits.value( name );
    if ( lineEdit )
    {
        lineEdit->setText( text );
    }
}

void TestWindow::updateComboBox( const QString &name, int index )
{
    if ( comboBox && comboBox->objectName() == name )
    {
        comboBox->setCurrentIndex( index );
    }
}

void TestWindow::updateDoubleSpinBox( const QString &name, double value )
{
    auto doubleSpinBox = doubleSpinBoxes.value( name );
    if ( doubleSpinBox )
    {
        doubleSpinBox->setValue( value );
    }
}

void TestWindow::updateIntSpinBox( const QString &name, int value )
{
    auto intSpinBox = intSpinBoxes.value( name );
    if ( intSpinBox )
    {
        intSpinBox->setValue( value );
    }
}

void TestWindow::updatePairBox( const QString &name, const QString &text )
{
    auto pairBox = pairBoxes.value( name );
    if ( pairBox )
    {
        if ( pairBox->morisPairComboBox->objectName() == name )
        {
            pairBox->morisPairComboBox->setCurrentText( text );
        }
        else if ( pairBox->morisPairLineEdit->objectName() == name )
        {
            pairBox->morisPairLineEdit->setText( text );
        }
        qDebug() << "Pair Box updated:" << name << "=" << text;
    }
}

void TestWindow::saveAndPrintInputs()
{
    qDebug() << "Saved Input Set:";

    // Print inputs for LineEdits
    for ( auto it = lineEdits.begin(); it != lineEdits.end(); ++it )
    {
        const QString   &name     = it.key();
        Moris_Line_Edit *lineEdit = it.value();
        if ( lineEdit )
        {
            const QString &text = lineEdit->text();
            qDebug() << name << ": " << text;
        }
        else
        {
            qDebug() << name << ": LineEdit is null";
        }
    }

    // Print input for ComboBox
    if ( comboBox )
    {
        const QString &name  = comboBox->objectName();
        const int      index = comboBox->currentIndex();
        qDebug() << name << ": " << index;
    }
    else
    {
        qDebug() << "ComboBox is null";
    }

    // Print inputs for DoubleSpinBoxes
    for ( auto it = doubleSpinBoxes.begin(); it != doubleSpinBoxes.end(); ++it )
    {
        const QString         &name          = it.key();
        Moris_Double_Spin_Box *doubleSpinBox = it.value();
        if ( doubleSpinBox )
        {
            const double value = doubleSpinBox->value();
            qDebug() << name << ": " << value;
        }
        else
        {
            qDebug() << name << ": DoubleSpinBox is null";
        }
    }

    // Print inputs for IntSpinBoxes
    for ( auto it = intSpinBoxes.begin(); it != intSpinBoxes.end(); ++it )
    {
        const QString      &name       = it.key();
        Moris_Int_Spin_Box *intSpinBox = it.value();
        if ( intSpinBox )
        {
            const int value = intSpinBox->value();
            qDebug() << name << ": " << value;
        }
        else
        {
            qDebug() << name << ": IntSpinBox is null";
        }
    }

    // Print inputs for PairBoxes (Moris_Pair_Box)
    for ( auto it = pairBoxes.begin(); it != pairBoxes.end(); ++it )
    {
        const QString  &name    = it.key();
        Moris_Pair_Box *pairBox = it.value();
        if ( pairBox )
        {
            const QString &text = pairBox->morisPairLineEdit->text();
            qDebug() << name << ": " << text;
        }
        else
        {
            qDebug() << name << ": PairBox is null";
        }
    }
}
