#include <QDebug>

#include "moris_line_edit.hpp"

// Constructor for Moris_Line_Edit
// Initializes the line edit widget and sets up the signal-slot connection for text changes
// Inputs:
// - parent: Pointer to the parent widget, default is nullptr
// Outputs: None
Moris_Line_Edit::Moris_Line_Edit( QWidget *parent )
        : QLineEdit( parent )
{
    // Connect the QLineEdit's textChanged signal to the onTextChanged slot
    connect( this, &QLineEdit::textChanged, this, &Moris_Line_Edit::onTextChanged );

    // Setting an example parameter for demonstration
    setParameter( "exampleKey", "exampleValue" );
}

// Destructor for Moris_Line_Edit
// Currently, it uses the default destructor
// Inputs: None
// Outputs: None
Moris_Line_Edit::~Moris_Line_Edit() = default;

// Sets a parameter key-value pair
// Inputs:
// - key: The parameter key (QString)
// - value: The parameter value (QVariant)
// Outputs: None
void Moris_Line_Edit::setParameter( const QString &key, const QVariant &value )
{
    parameters.insert( key, value );
}

// Retrieves the value for a given parameter key
// Inputs:
// - key: The parameter key (QString)
// Outputs:
// - Returns the value associated with the key (QVariant)
QVariant Moris_Line_Edit::parameter( const QString &key ) const
{
    return parameters.value( key );
}

// Slot that handles the line edit's text change
// Inputs:
// - text: The new text (QString)
// Outputs: None
void Moris_Line_Edit::onTextChanged( const QString &text )
{
    // Emit the custom textChanged signal with the object name and the new text
    emit textChanged( objectName(), text );
}
