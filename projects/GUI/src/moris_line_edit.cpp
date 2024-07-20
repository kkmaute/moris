#include "moris_line_edit.hpp"

// Constructor for Moris_Line_Edit
// Inputs:
// - parent: Pointer to the parent widget (default is nullptr).
// - parameter: Pointer to a moris::Parameter object to be linked with this widget (default is nullptr).
Moris_Line_Edit::Moris_Line_Edit( QWidget *parent, moris::Parameter &parameter )
        : QLineEdit( parent )
        , mParameter( parameter )
{
    // Connect the textChanged(const QString &) signal to the onTextChanged slot
    connect( this, &QLineEdit::textChanged, this, &Moris_Line_Edit::onTextChanged );

    // If parameter is not null, set the initial text from the parameter value
        setText( QString::fromStdString( mParameter.get_string() ) );
}

// Destructor for Moris_Line_Edit
Moris_Line_Edit::~Moris_Line_Edit() = default;

// Getter for the associated moris::Parameter object
moris::Parameter &Moris_Line_Edit::getParameter() 
{
    return mParameter;
}

// Slot to handle text changes
// This slot is connected to the textChanged signal of QLineEdit and updates the linked moris::Parameter object.
// Inputs:
// - new_text: New text input in the widget.
void Moris_Line_Edit::onTextChanged( const QString &new_text )
{
    // If parameter is not null, update the parameter value with the new text

        mParameter.set_value( objectName().toStdString(), new_text.toStdString(), false );


    // Emit the custom textChanged signal with the widget's name and new text
    emit textChanged( objectName(), new_text );
}
