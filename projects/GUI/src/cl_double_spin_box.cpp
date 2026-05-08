#include "cl_double_spin_box.hpp"

namespace moris
{
    // Constructor for Moris_Double_Spin_Box
    // Initializes the double spin box widget and sets up its signal-slot connections.
    // Inputs:
    // - a_parent: Pointer to the parent widget (default is nullptr).
    // - a_parameter: Reference to a Parameter object to be linked with this widget.
    Moris_Double_Spin_Box::Moris_Double_Spin_Box( QWidget *a_parent, Parameter &a_parameter )
            : QDoubleSpinBox( a_parent )
            , mParameter( &a_parameter )
    {
        // Connect the valueChanged(double) signal of QDoubleSpinBox to the on_value_changed slot
        connect( this,
                 QOverload< double >::of( &QDoubleSpinBox::valueChanged ),
                 this, 
                 &Moris_Double_Spin_Box::on_value_changed );
        
        refreshDataParameter();

        if ( mParameter->is_locked() )
        {
            setReadOnly( true );
        }
        else
        {
            setReadOnly( false );
        }
    }

    // Destructor for Moris_Double_Spin_Box
    // The destructor is defaulted as there are no specific cleanup requirements.
    Moris_Double_Spin_Box::~Moris_Double_Spin_Box() = default;

    // Getter for the associated Parameter object
    // Returns the reference to the parameter linked with this widget.
    // Outputs:
    // - Reference to the Parameter object.
    Parameter &Moris_Double_Spin_Box::get_parameter()
    {
        MORIS_ERROR(mParameter, "Moris_Double_Spin_Box::get_parameter() called with null mParameter.");
        return *mParameter;
    }

    void Moris_Double_Spin_Box::setParameter( Parameter &parameter )
    {
        // qDebug() << "[DoubleSpinBox::setParameter]"
        //         << "widget =" << this
        //         << "name =" << objectName()
        //         << "old mParameter =" << mParameter
        //         << "new mParameter =" << &parameter;
        mParameter = &parameter;

        if( mParameter && mParameter->is_locked() )
        {
            setReadOnly( true );
        }
        else
        {
            setReadOnly( false );
        }
        refreshDataParameter();
    }

    void Moris_Double_Spin_Box::refreshDataParameter()
    {
        QSignalBlocker tBlocker( this ); // block signals to prevent recursive updates

        if(!mParameter)
        {
            setRange(-MORIS_REAL_MAX, MORIS_REAL_MAX);
            setValue( 0.0 );
            return;
        }

        if( mParameter->index() == variant_index< real >() )
        {
            setRange(-MORIS_REAL_MAX, MORIS_REAL_MAX);
            setValue( mParameter->get_value< real >() );
        }
        else
        {
            setRange(-MORIS_REAL_MAX, MORIS_REAL_MAX);
            bool tOk = false;
            double tValue = QString::fromStdString( mParameter->get_string() ).toDouble(&tOk);

            setValue( tOk ? tValue : 0.0 );
        }
    }

    // Slot to handle value changes
    // This slot is connected to the valueChanged(double) signal of QDoubleSpinBox and updates the linked Parameter object.
    // Inputs:
    // - a_value: New double value input in the widget.
    void Moris_Double_Spin_Box::on_value_changed( double a_value )
    {

        if(!mParameter || mParameter->is_locked()) return;
        // Update the parameter value with the new double value
        mParameter->set_value( objectName().toStdString(), a_value, false );

        // Emit the custom value_changed signal with the widget's name and the new double value
        emit value_changed( objectName(), a_value );
    }
}    // namespace moris
