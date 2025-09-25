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

    // Template function to set the parameter value based on its type
    // Inputs:
    // - aParameterName: Name of the parameter to be set.
    // - aValue: Value to be set for the parameter.
    // - aNewText: New text input in the widget.
    template<typename T>
    void Moris_Line_Edit::mTrySetParameter(const std::string& aParameterName, const T& aValue, const QString& aNewText)
    {
        mParameter.set_value(aParameterName, aValue, false);
        emit textChanged(objectName(), aNewText);
    }

    // Slot to handle text changes
    // This slot is connected to the textChanged signal of QLineEdit and updates the linked Parameter object.
    // Inputs:
    // - new_text: New text input in the widget.
    void Moris_Line_Edit::onTextChanged( const QString &new_text )
    {   
        // early escape if the parameter is locked
        if(mParameter.is_locked()) return;

        // convert to string so to avoid unnecessary conversions later
        const std::string tParameterName = objectName().toStdString();
        const std::string tNewTextString = new_text.toStdString();

        if(mParameter.index() == variant_index<std::string>())
        {
           mTrySetParameter(tParameterName, tNewTextString, new_text);
           return;
        }
        
        if(mParameter.index() == variant_index<std::pair< std::string, std::string > >())
        {
            Vector< std::string > tVec = split_string( tNewTextString, "," );
            if ( tVec.size() < 2 ) return;
            auto tPair = std::make_pair( tVec( 0 ), tVec( 1 ) );
            mTrySetParameter(tParameterName, tPair, new_text);
            return;
        }

        if(mParameter.index() == variant_index< Vector< uint > >() )
        {
            auto tVec = string_to_vector< uint >( tNewTextString );
            if(tVec.empty()) return;
            mTrySetParameter(tParameterName, tVec, new_text);
            return;
        }

        if(mParameter.index() == variant_index< Vector< sint > >() )
        {
            auto tVec = string_to_vector< sint >( tNewTextString );
            if(tVec.empty()) return;
            mTrySetParameter(tParameterName, tVec, new_text);
            return;
        }

        if(mParameter.index() == variant_index< Vector< real > >() )
        {
            auto tVec = string_to_vector< real >( tNewTextString );
            if(tVec.empty()) return;
            mTrySetParameter(tParameterName, tVec, new_text);
            return;
        }

        if(mParameter.index() == variant_index< Vector< std::string > >() )
        {
            auto tVec = string_to_vector< std::string >( tNewTextString );
            if(tVec.empty()) return;
            mTrySetParameter(tParameterName, tVec, new_text);
            return;
        }

        
        {
            // Geometry center_x variable spazzes out when empty
            Vector< real >  tVec      = string_to_vector< real >( tNewTextString );
            if ( tVec.empty() ) return;
            Design_Variable tDesign_Variable = Design_Variable( tVec( 0 ) );
            mTrySetParameter(tParameterName, tDesign_Variable, new_text);
        }

    }
}    // namespace moris