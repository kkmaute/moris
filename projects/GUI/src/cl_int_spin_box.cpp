#include "cl_int_spin_box.hpp"

namespace moris
{
    // Constructor for Moris_Int_Spin_Box.
    // Initializes the spin box widget and sets its initial value based on the provided parameter.
    // Connects the valueChanged signal to the appropriate slot for handling changes.
    // Inputs:
    // - a_parent: Pointer to the parent widget (default is nullptr).
    // - a_parameter: Reference to a Parameter object to be linked with this widget.
    Moris_Int_Spin_Box::Moris_Int_Spin_Box( QWidget *a_parent, Parameter &a_parameter )
            : QSpinBox( a_parent )
            , m_parameter( a_parameter )
    {
        // Set the initial value from the parameter value based on its type
        if ( m_parameter.index() == variant_index< uint >() )
        {
            setValue( m_parameter.get_value< uint >() );
            setRange( 0, INT_MAX );
        }
        else
        {
            setRange( -INT_MAX, INT_MAX );
            setValue( m_parameter.get_value< sint >() );
        }
        // Connect the valueChanged(int) signal of QSpinBox to the on_value_changed slot
        connect( this, QOverload< int >::of( &QSpinBox::valueChanged ), this, &Moris_Int_Spin_Box::on_value_changed );
    }

    // Destructor for Moris_Int_Spin_Box.
    // The destructor is defaulted as there are no specific cleanup requirements.
    Moris_Int_Spin_Box::~Moris_Int_Spin_Box() = default;

    // Getter for the associated Parameter object.
    // Returns the reference to the parameter linked with this widget.
    // Outputs:
    // - Reference to the Parameter object.
    Parameter &Moris_Int_Spin_Box::get_parameter()
    {
        return m_parameter;
    }

    // Slot to handle value changes in the spin box.
    // Updates the linked Parameter object with the new value whenever the spin box value changes.
    // Inputs:
    // - a_value: The new integer value input in the widget.
    void Moris_Int_Spin_Box::on_value_changed( int a_value )
    {
        // Update the parameter with the new value
        m_parameter.set_value( objectName().toStdString(), a_value, false );

        // Emit the custom value_changed signal with the widget's name and the new integer value
        emit value_changed( objectName(), a_value );
    }
}    // namespace moris
