#include "moris_int_spin_box.hpp"

// Constructor for Moris_Int_Spin_Box.
// Initializes the spin box widget and sets its initial value based on the provided parameter.
// Connects the valueChanged signal to the appropriate slot for handling changes.
// Inputs:
// - parent: Pointer to the parent widget (default is nullptr).
// - parameter: Reference to a moris::Parameter object to be linked with this widget.
Moris_Int_Spin_Box::Moris_Int_Spin_Box( QWidget *parent, moris::Parameter &parameter )
        : QSpinBox( parent )
        , mParameter( parameter )
{
    // Set the initial value from the parameter value based on its type
    if ( mParameter.index() == moris::variant_index< uint >() )
    {
        setValue( mParameter.get_value< uint >() );
        setRange( 0, INT_MAX );
    }
    else
    {
        setRange( -INT_MAX, INT_MAX );
        setValue( mParameter.get_value< moris::sint >() );
    }
    // Connect the valueChanged(int) signal of QSpinBox to the onValueChanged slot
    connect( this, QOverload< int >::of( &QSpinBox::valueChanged ), this, &Moris_Int_Spin_Box::onValueChanged );
}

// Destructor for Moris_Int_Spin_Box.
// The destructor is defaulted as there are no specific cleanup requirements.
Moris_Int_Spin_Box::~Moris_Int_Spin_Box() = default;

// Getter for the associated moris::Parameter object.
// Returns the reference to the parameter linked with this widget.
// Output:
// - Reference to the moris::Parameter object.
moris::Parameter &Moris_Int_Spin_Box::getParameter()
{
    return mParameter;
}

// Slot to handle value changes in the spin box.
// Updates the linked moris::Parameter object with the new value whenever the spin box value changes.
// Inputs:
// - value: The new integer value input in the widget.
void Moris_Int_Spin_Box::onValueChanged( int value )
{
    // Update the parameter with the new value
    mParameter.set_value( objectName().toStdString(), value, false );

    // Emit the custom valueChanged signal with the widget's name and the new integer value
    emit valueChanged( objectName(), value );
}
