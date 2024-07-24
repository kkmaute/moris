#include "moris_line_edit.hpp"

// Constructor for Moris_Line_Edit
// Inputs:
// - parent: Pointer to the parent widget (default is nullptr).
// - parameter: Reference to a moris::Parameter object to be linked with this widget.
Moris_Line_Edit::Moris_Line_Edit( QWidget *parent, moris::Parameter &parameter )
        : QLineEdit( parent )
        , mParameter( parameter )
{
    // If parameter is not null, set the initial text from the parameter value
        setPlaceholderText( QString::fromStdString( mParameter.get_string() ) );
        connect( this, &QLineEdit::textChanged, this, &Moris_Line_Edit::onTextChanged );

    // Connect the textChanged signal of QLineEdit to the onTextChanged slot
    connect( this, &QLineEdit::textChanged, this, &Moris_Line_Edit::onTextChanged );
}

// Destructor for Moris_Line_Edit
// No specific cleanup is required, so the destructor is defaulted.
Moris_Line_Edit::~Moris_Line_Edit() = default;

// Getter for the associated moris::Parameter object
// Returns:
// - Reference to the moris::Parameter object linked with this widget.
moris::Parameter &Moris_Line_Edit::getParameter()
{
    return mParameter;
}

// Setter for the associated moris::Parameter object
// Inputs:
// - parameter: Reference to a moris::Parameter object to be linked with this widget.
void Moris_Line_Edit::setParameter( moris::Parameter &parameter )
{
    mParameter = parameter;
}

// Slot to handle text changes
// This slot is connected to the textChanged signal of QLineEdit and updates the linked moris::Parameter object.
// Inputs:
// - new_text: New text input in the widget.
void Moris_Line_Edit::onTextChanged( const QString &new_text )
{

    // Switch between different types based on the parameter index and update the parameter value
    if ( mParameter.index() == moris::variant_index< std::string >() )
    {
        mParameter.set_value( objectName().toStdString(), new_text.toStdString(), false );
    }
    else if ( mParameter.index() == moris::variant_index< std::pair< std::string, std::string > >() )
    {
        moris::Vector< std::string >          tVec  = moris::string_to_cell< std::string >( new_text.toStdString() );
        std::pair< std::string, std::string > tPair = std::make_pair( tVec( 0 ), tVec( 1 ) );
        mParameter.set_value( objectName().toStdString(), tPair, false );
    }
    else if ( mParameter.index() == moris::variant_index< moris::Vector< uint > >() )
    {
        moris::Vector< uint > tVec = moris::string_to_cell< uint >( new_text.toStdString() );
        mParameter.set_value( objectName().toStdString(), tVec, false );
    }
    else if ( mParameter.index() == moris::variant_index< moris::Vector< moris::sint > >() )
    {
        moris::Vector< moris::sint > tVec = moris::string_to_cell< moris::sint >( new_text.toStdString() );
        mParameter.set_value( objectName().toStdString(), tVec, false );
    }
    else if ( mParameter.index() == moris::variant_index< moris::Vector< moris::real > >() )
    {
        moris::Vector< moris::real > tVec = moris::string_to_cell< moris::real >( new_text.toStdString() );
        mParameter.set_value( objectName().toStdString(), tVec, false );
    }
    else if ( mParameter.index() == moris::variant_index< moris::Vector< std::string > >() )
    {
        moris::Vector< std::string > tVec = moris::string_to_cell< std::string >( new_text.toStdString() );
        mParameter.set_value( objectName().toStdString(), tVec, false );
    }
    else
    {
        moris::Vector< moris::real > tVec             = moris::string_to_cell< moris::real >( new_text.toStdString() );
        moris::Design_Variable       tDesign_Variable = moris::Design_Variable( tVec( 0 ), tVec( 1 ), tVec( 2 ) );
        mParameter.set_value( objectName().toStdString(), tDesign_Variable, false );
    }

    // Emit the custom textChanged signal with the widget's name and new text
    emit textChanged( objectName(), new_text );
}
