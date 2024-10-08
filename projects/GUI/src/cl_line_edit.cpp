#include "cl_line_edit.hpp"

namespace moris
{
    // Constructor for Moris_Line_Edit
    // Inputs:
    // - parent: Pointer to the parent widget (default is nullptr).
    // - parameter: Reference to a Parameter object to be linked with this widget.
    Moris_Line_Edit::Moris_Line_Edit( QWidget *parent, Parameter &parameter )
            : QLineEdit( parent )
            , mParameter( parameter )
    {

        if ( mParameter.index() == variant_index< std::string >() )
        {
            setText( QString::fromStdString( mParameter.get_value< std::string >() ) );
        }
        else
        {
            setText( QString::fromStdString( mParameter.get_string() ) );
        }

        // If parameter is not null, set the initial text from the parameter value

        // Connect the textChanged signal of QLineEdit to the onTextChanged slot
        if ( mParameter.is_locked() )
        {
            setReadOnly( true );
        }
        else
        {
            connect( this, &QLineEdit::textChanged, this, &Moris_Line_Edit::onTextChanged );
        }
    }

    // Destructor for Moris_Line_Edit
    // No specific cleanup is required, so the destructor is defaulted.
    Moris_Line_Edit::~Moris_Line_Edit() = default;

    // Getter for the associated Parameter object
    // Returns:
    // - Reference to the Parameter object linked with this widget.
    Parameter &Moris_Line_Edit::getParameter()
    {
        return mParameter;
    }

    // Setter for the associated Parameter object
    // Inputs:
    // - parameter: Reference to a Parameter object to be linked with this widget.
    void Moris_Line_Edit::setParameter( Parameter &parameter )
    {
        mParameter = parameter;
    }

    // Slot to handle text changes
    // This slot is connected to the textChanged signal of QLineEdit and updates the linked Parameter object.
    // Inputs:
    // - new_text: New text input in the widget.
    void Moris_Line_Edit::onTextChanged( const QString &new_text )
    {

        // Switch between different types based on the parameter index and update the parameter value
        if ( mParameter.index() == variant_index< std::string >() )
        {
            mParameter.set_value( objectName().toStdString(), new_text.toStdString(), false );
        }
        else if ( mParameter.index() == variant_index< std::pair< std::string, std::string > >() )
        {
            Vector< std::string >                 tVec  = split_string( new_text.toStdString(), "," );
            std::pair< std::string, std::string > tPair = std::make_pair( tVec( 0 ), tVec( 1 ) );
            mParameter.set_value( objectName().toStdString(), tPair, false );
        }
        else if ( mParameter.index() == variant_index< Vector< uint > >() )
        {
            Vector< uint > tVec = string_to_vector< uint >( new_text.toStdString() );
            mParameter.set_value( objectName().toStdString(), tVec, false );
        }
        else if ( mParameter.index() == variant_index< Vector< sint > >() )
        {
            Vector< sint > tVec = string_to_vector< sint >( new_text.toStdString() );
            mParameter.set_value( objectName().toStdString(), tVec, false );
        }
        else if ( mParameter.index() == variant_index< Vector< real > >() )
        {
            Vector< real > tVec = string_to_vector< real >( new_text.toStdString() );
            mParameter.set_value( objectName().toStdString(), tVec, false );
        }
        else if ( mParameter.index() == variant_index< Vector< std::string > >() )
        {
            Vector< std::string > tVec = string_to_vector< std::string >( new_text.toStdString() );
            mParameter.set_value( objectName().toStdString(), tVec, false );
        }
        else
        {
            // Geometry center_x variable spazzes out when empty
            Vector< real >  tVec             = string_to_vector< real >( new_text.toStdString() );
            Design_Variable tDesign_Variable = Design_Variable( tVec( 0 ) );
            mParameter.set_value( objectName().toStdString(), tDesign_Variable, false );
        }

        // Emit the custom textChanged signal with the widget's name and new text
        emit textChanged( objectName(), new_text );
    }
}    // namespace moris