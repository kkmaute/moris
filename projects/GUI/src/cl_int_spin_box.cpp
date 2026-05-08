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
            , mParameter( &a_parameter )
    {
        connect(
            this, 
            QOverload< int >::of( &QSpinBox::valueChanged ), 
            this, 
            &Moris_Int_Spin_Box::on_value_changed); // always connect once

        refreshDataParameter();

        // Connect the valueChanged(int) signal of QSpinBox to the on_value_changed slot
        if ( mParameter->is_locked() )
        {
            setReadOnly( true );
        }
        else
        {
            setReadOnly( false );
        }
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
        MORIS_ERROR(mParameter, "Moris_Int_Spin_Box::get_parameter() called with null mParameter.");
        return *mParameter;
    }

    void Moris_Int_Spin_Box::setParameter( Parameter &parameter )
    {
        // qDebug() << "[IntSpinBox::setParameter]"
        //         << "widget =" << this
        //         << "name =" << objectName()
        //         << "old mParameter =" << mParameter
        //         << "new mParameter =" << &parameter;

        mParameter = &parameter;


        refreshDataParameter();

        if ( mParameter && mParameter->is_locked() )
        {
            setReadOnly( true );
        }
        else
        {
            setReadOnly( false );
        }
    }

    void Moris_Int_Spin_Box::refreshDataParameter()
    {
        QSignalBlocker tBlocker( this ); // block signals to prevent recursive updates

        if(!mParameter)
        {
            setRange(INT_MIN, INT_MAX);
            setValue( 0 );
            return;
        }

        if(mParameter->index() == variant_index< uint >() )
        {
            if ( objectName() == "NLA_Solver_Implementation" )
                {
                    setRange( 0, 2);
                    uint tValue = mParameter->get_value<uint>();
                    if(tValue > static_cast<uint>(2))
                    {
                        tValue = static_cast<uint>(2);
                    }
                    setValue( static_cast<int>(tValue) );
                    return;
                }
            setRange( 0, INT_MAX ); 
            uint tValue = mParameter->get_value<uint>();
            if(tValue > static_cast<uint>(INT_MAX))
            {
                tValue = static_cast<uint>(INT_MAX);
            }
            setValue( static_cast<int>(tValue) );

        }
        else if(mParameter->index() == variant_index< sint >() )
        {

            setRange( INT_MIN, INT_MAX );
            setValue( static_cast<int>(mParameter->get_value<sint>()) );
        }
        else
        {
            setRange( INT_MIN, INT_MAX );
            bool tOk = false;
            int tValue = QString::fromStdString( mParameter->get_string() ).toInt(&tOk);

            setValue( tOk ? tValue : 0 );
        }

    }

    // Slot to handle value changes in the spin box.
    // Updates the linked Parameter object with the new value whenever the spin box value changes.
    // Inputs:
    // - a_value: The new integer value input in the widget.
    void Moris_Int_Spin_Box::on_value_changed( int a_value )
    {

        if(!mParameter || mParameter->is_locked()) return;

        const std::string tName = objectName().toStdString();
        if(mParameter->index() == variant_index<uint>() )
        {
            if(a_value < 0)
            {
                return;
            }
            mParameter->set_value(  tName, static_cast<uint>(a_value), false );
        }
        else
        {
            mParameter->set_value( tName, static_cast<sint>(a_value), false);
        }


        // Update the parameter with the new value


        // Emit the custom value_changed signal with the widget's name and the new integer value
        emit value_changed( objectName(), a_value );
    }
}    // namespace moris
