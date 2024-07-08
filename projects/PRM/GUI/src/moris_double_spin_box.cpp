#include "moris_double_spin_box.hpp"

// Moris_Double_Spin_Box constructor
// Inputs:
// - parent: The parent widget (QWidget *)
Moris_Double_Spin_Box::Moris_Double_Spin_Box( QWidget *parent )
        : QDoubleSpinBox( parent )
{
    // Set the range for the spin box using setMinimum and setMaximum
    setMinimum( 0 );
    setMaximum( static_cast< double >( std::numeric_limits< long long >::max() ) );

    // Set the number of decimal places to display
    setDecimals( 4 );    // Display 4 digits after the decimal point

    // Connect the valueChanged signal to the onValueChanged slot
    connect( this, QOverload< double >::of( &QDoubleSpinBox::valueChanged ), this, &Moris_Double_Spin_Box::onValueChanged );
}

// Moris_Double_Spin_Box destructor
Moris_Double_Spin_Box::~Moris_Double_Spin_Box()
{
    // Destructor
}

// Sets a parameter key-value pair in the parameters map
// Inputs:
// - key: The parameter key (QString)
// - value: The parameter value (QVariant)
// Outputs: None
void Moris_Double_Spin_Box::setParameter( const QString &key, const QVariant &value )
{
    // Insert the key-value pair into the parameters map
    parameters.insert( key, value );
}

// Retrieves the value for a given parameter key from the parameters map
// Inputs:
// - key: The parameter key (QString)
// Outputs:
// - Returns the value associated with the key (QVariant)
QVariant Moris_Double_Spin_Box::parameter( const QString &key ) const
{
    // Return the value associated with the provided key from the parameters map
    return parameters.value( key );
}

// Slot to handle value changes in the spin box
// Inputs:
// - new_value: The new value of the spin box (double)
void Moris_Double_Spin_Box::onValueChanged( double new_value )
{
    QVariant value;
    if ( new_value == static_cast< int >( new_value ) )
    {
        value = static_cast< int >( new_value );    // Save as integer if it's a whole number
    }
    else
    {
        value = new_value;    // Save as double otherwise
    }

    // Emit the valueChanged signal with the object name and the new value
    emit valueChanged( objectName(), value );
}
