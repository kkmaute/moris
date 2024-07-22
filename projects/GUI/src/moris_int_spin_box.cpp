#include "moris_int_spin_box.hpp"

// Constructor for Moris_Int_Spin_Box
// Inputs:
// - parent: Pointer to the parent widget (default is nullptr).
// - parameter: Pointer to a moris::Parameter object to be linked with this widget (default is nullptr).
Moris_Int_Spin_Box::Moris_Int_Spin_Box( QWidget *parent, moris::Parameter &parameter )
        : QSpinBox( parent )
        , mParameter( parameter )
{

    // If parameter is not null, set the initial value from the parameter value
    if (mParameter.index() == moris::variant_index<uint>()) {
        setValue( mParameter.get_value< uint >() );
    }
    else {
        setValue( mParameter.get_value< int >() );
        
    }
    // Connect the valueChanged(int) signal of QSpinBox to the onValueChanged slot
    connect( this, QOverload< int >::of( &QSpinBox::valueChanged ), this, &Moris_Int_Spin_Box::onValueChanged );
}

// Destructor for Moris_Int_Spin_Box
Moris_Int_Spin_Box::~Moris_Int_Spin_Box() = default;

// Getter for the associated moris::Parameter object
moris::Parameter &Moris_Int_Spin_Box::getParameter()
{
    return mParameter;
}

// Slot to handle value changes
// This slot is connected to the valueChanged(int) signal of QSpinBox and updates the linked moris::Parameter object.
// Inputs:
// - value: New integer value input in the widget.
void Moris_Int_Spin_Box::onValueChanged( int value )
{

        mParameter.set_value( objectName().toStdString(), value, false );


    // Emit the custom valueChanged signal with the widget's name and the new integer value
    emit valueChanged( objectName(), value );
}
