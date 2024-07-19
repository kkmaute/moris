#include "moris_line_edit.hpp"

// Constructor for Moris_Line_Edit
// Inputs:
// - parent: Pointer to the parent widget (default is nullptr).
// - parameter: Pointer to a moris::Parameter object to be linked with this widget (default is nullptr).
Moris_Line_Edit::Moris_Line_Edit( QWidget *parent, moris::Parameter *parameter )
        : QLineEdit( parent )
        , mParameter( parameter )
{
    // If parameter is not null, set the initial text from the parameter value
    if ( mParameter )
    {
        if ( mParameter->get_string() != "\"\"" )
        {
            setPlaceholderText( QString::fromStdString( mParameter->get_string() ) );
        }
    }

    // Connect the textChanged(const QString &) signal to the onTextChanged slot
    connect( this, &QLineEdit::textChanged, this, &Moris_Line_Edit::onTextChanged );

}

// Destructor for Moris_Line_Edit
Moris_Line_Edit::~Moris_Line_Edit() = default;

// Getter for the associated moris::Parameter object
moris::Parameter *Moris_Line_Edit::getParameter() const
{
    return mParameter;
}

void Moris_Line_Edit::setParameter( moris::Parameter *parameter )
{
    mParameter = parameter;
}

// Slot to handle text changes
// This slot is connected to the textChanged signal of QLineEdit and updates the linked moris::Parameter object.
// Inputs:
// - new_text: New text input in the widget.
void Moris_Line_Edit::onTextChanged(const QString &new_text )
{
    // If parameter is not null, update the parameter value with the new text

    // Parameter needs to check with index and save the type

    // Switch between the different indices of pararmeters except bool, uint, sint, real
    // Based on the index, call string_to_cell to create the Vector of the specific type
    // Pass the vector into the set_value function
    // iElements.second.index() == variant_index <sint> ()

    if ( mParameter )
    {
        std::cout << mParameter->index() << std::endl;
        if ( mParameter->index() == moris::variant_index< std::string >() )
        {
            std::cout << "string" << std::endl;
            mParameter->set_value( objectName().toStdString(), new_text.toStdString(), false );
        }
        else if ( mParameter->index() == moris::variant_index< std::pair< std::string, std::string > >() )
        {
            std::cout << "pair" << std::endl;
            moris::Vector< std::string >          tVec  = moris::string_to_cell< std::string >( new_text.toStdString() );
            std::pair< std::string, std::string > tPair = std::make_pair( tVec( 0 ), tVec( 1 ) );
            mParameter->set_value( objectName().toStdString(), tPair, false );
        }
        else if ( mParameter->index() == moris::variant_index< moris::Vector< uint > >() )
        {
            std::cout << "Vector: uint" << std::endl;
            moris::Vector< uint > tVec = moris::string_to_cell< uint >( new_text.toStdString() );
            mParameter->set_value( objectName().toStdString(), tVec, false );
        }
        else if ( mParameter->index() == moris::variant_index< moris::Vector< moris::sint > >() )
        {
            std::cout << "Vector: sint" << std::endl;
            moris::Vector< moris::sint > tVec = moris::string_to_cell< moris::sint >( new_text.toStdString() );
            mParameter->set_value( objectName().toStdString(), tVec, false );
        }
        else if ( mParameter->index() == moris::variant_index< moris::Vector< moris::real > >() )
        {
            std::cout << "Vector: real" << std::endl;
            moris::Vector< moris::real > tVec = moris::string_to_cell< moris::real >( new_text.toStdString() );
            mParameter->set_value( objectName().toStdString(), tVec, false );
        }
        else if ( mParameter->index() == moris::variant_index< moris::Vector< std::string > >() )
        {
            std::cout << "Vector: string" << std::endl;
            moris::Vector< std::string > tVec = moris::string_to_cell< std::string >( new_text.toStdString() );
            mParameter->set_value( objectName().toStdString(), tVec, false );
        }
        else
        {
            std::cout << "Design_Variable" << std::endl;
            moris::Vector< moris::real > tVec             = moris::string_to_cell< moris::real >( new_text.toStdString() );
            moris::Design_Variable       tDesign_Variable = moris::Design_Variable( tVec( 0 ), tVec( 1 ), tVec( 2 ) );
            mParameter->set_value( objectName().toStdString(), tDesign_Variable, false );
        }
        std::cout << mParameter->get_string() << std::endl;
    }

    // Emit the custom textChanged signal with the widget's name and new text
    emit textChanged( objectName(), new_text );
}
