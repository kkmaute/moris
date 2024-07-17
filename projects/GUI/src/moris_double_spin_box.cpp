#include "moris_double_spin_box.hpp"

// Constructor for Moris_Double_Spin_Box
// Inputs:
// - parent: Pointer to the parent widget (default is nullptr).
// - parameter: Pointer to a moris::Parameter object to be linked with this widget (default is nullptr).
Moris_Double_Spin_Box::Moris_Double_Spin_Box( QWidget *parent, moris::Parameter *parameter )
        : QDoubleSpinBox( parent )
        , mParameter( parameter )
{
    // Connect the valueChanged(double) signal of QDoubleSpinBox to the onValueChanged slot
    connect( this, QOverload< double >::of( &QDoubleSpinBox::valueChanged ), this, &Moris_Double_Spin_Box::onValueChanged );

    // If parameter is not null, set the initial value from the parameter value
    if ( mParameter )
    {
        setValue( mParameter->get_value< double >() );
    }
}

// Destructor for Moris_Double_Spin_Box
Moris_Double_Spin_Box::~Moris_Double_Spin_Box() = default;

// Getter for the associated moris::Parameter object
moris::Parameter* Moris_Double_Spin_Box::getParameter() const
{
    return mParameter;
}

// Slot to handle value changes
// This slot is connected to the valueChanged(double) signal of QDoubleSpinBox and updates the linked moris::Parameter object.
// Inputs:
// - value: New double value input in the widget.
void Moris_Double_Spin_Box::onValueChanged( double value )
{
    // If parameter is not null, update the parameter value with the new double value
    if ( mParameter )
    {
        mParameter->set_value( objectName().toStdString(), value, false );
    }

    // Emit the custom valueChanged signal with the widget's name and the new double value
    emit valueChanged( objectName(), value );
}
