#include "moris_double_spin_box.hpp"

// Constructor for Moris_Double_Spin_Box
// Initializes the double spin box widget and sets up its signal-slot connections.
// Inputs:
// - parent: Pointer to the parent widget (default is nullptr).
// - parameter: Reference to a moris::Parameter object to be linked with this widget.
Moris_Double_Spin_Box::Moris_Double_Spin_Box( QWidget *parent, moris::Parameter &parameter )
        : QDoubleSpinBox( parent )
        , mParameter( parameter )
{
    
    setRange( -MORIS_REAL_MAX, MORIS_REAL_MAX );
    setValue( mParameter.get_value< moris::real  >() );
    // Connect the valueChanged(double) signal of QDoubleSpinBox to the onValueChanged slot
    connect( this, QOverload< double >::of( &QDoubleSpinBox::valueChanged ), this, &Moris_Double_Spin_Box::onValueChanged );
}

// Destructor for Moris_Double_Spin_Box
// The destructor is defaulted as there are no specific cleanup requirements.
Moris_Double_Spin_Box::~Moris_Double_Spin_Box() = default;

// Getter for the associated moris::Parameter object
// Returns the reference to the parameter linked with this widget.
// Outputs:
// - Reference to the moris::Parameter object.
moris::Parameter &Moris_Double_Spin_Box::getParameter()
{
    return mParameter;
}

// Slot to handle value changes
// This slot is connected to the valueChanged(double) signal of QDoubleSpinBox and updates the linked moris::Parameter object.
// Inputs:
// - value: New double value input in the widget.
void Moris_Double_Spin_Box::onValueChanged( double value )
{
    // Update the parameter value with the new double value
    mParameter.set_value( objectName().toStdString(), value, false );

    // Emit the custom valueChanged signal with the widget's name and the new double value
    emit valueChanged( objectName(), value );
}
