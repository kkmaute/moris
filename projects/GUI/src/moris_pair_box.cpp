#include "moris_pair_box.hpp"

// Constructor for Moris_Pair_Box
// Inputs:
// - parent: Pointer to the parent widget (default is nullptr).
// - param: Reference to a moris::Parameter object to be linked with this widget.
// - options: QStringList containing options for the combo box.
Moris_Pair_Box::Moris_Pair_Box( QWidget *parent, moris::Parameter &param, const QStringList &options )
        : QWidget( parent )
        , morisPairComboBox( new QComboBox( this ) )
        , morisPairLineEdit( new QLineEdit( this ) )
        , parameter( param )
{
    // Set up the combo box with provided options
    morisPairComboBox->addItems( options );
    connect( morisPairComboBox, &QComboBox::currentTextChanged, this, &Moris_Pair_Box::onComboBoxTextChanged );

    // Set up the line edit
    connect( morisPairLineEdit, &QLineEdit::textChanged, this, &Moris_Pair_Box::onLineEditTextChanged );

    // Layout the widgets in a horizontal layout
    auto *layout = new QHBoxLayout( this );
    layout->addWidget( morisPairComboBox );
    layout->addWidget( morisPairLineEdit );
    setLayout( layout );
}

// Slot to handle changes in the combo box text
// Inputs:
// - text: New text selected in the combo box.
void Moris_Pair_Box::onComboBoxTextChanged( const QString &text )
{
    emit comboBoxTextChanged( objectName(), text );
}

// Slot to handle changes in the line edit text
// Inputs:
// - text: New text entered in the line edit.
void Moris_Pair_Box::onLineEditTextChanged( const QString &text )
{
    emit lineEditTextChanged( objectName(), text );
}

// Set the associated moris::Parameter with the combined value of the combo box and line edit
void Moris_Pair_Box::setParameter()
{
    std::pair<std::string, std::string> pairValue;
    pairValue.first = morisPairComboBox->currentText().toStdString();
    pairValue.second = morisPairLineEdit->text().toStdString();
    parameter.set_value(objectName().toStdString(), pairValue, false);
}
