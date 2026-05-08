#include "cl_bool_combo_box.hpp"

namespace moris
{

    // Constructor for Moris_Bool_Combo_Box
    // Initializes the combo box widget and sets up its items and signal-slot connections.
    // Inputs:
    // - a_parent: Pointer to the parent widget (default is nullptr).
    // - a_parameter: Reference to a Parameter object to be linked with this widget.
    Moris_Bool_Combo_Box::Moris_Bool_Combo_Box( QWidget *a_parent, Parameter &a_parameter )
            : QComboBox( a_parent )
            , mParameter( &a_parameter )
    {
        connect(
                this, 
                QOverload< int >::of( &QComboBox::currentIndexChanged ), 
                this, 
                &Moris_Bool_Combo_Box::on_index_changed
        );

        refreshDataParameter();

        // Connect the currentIndexChanged(int) signal of QComboBox to the on_index_changed slot
        if ( mParameter && mParameter->is_locked() )
        {
            setDisabled( true );
        }
        else
        {
            setDisabled( false ); 
        }
    }

    // Destructor for Moris_Bool_Combo_Box
    // The destructor is defaulted as there are no specific cleanup requirements.
    Moris_Bool_Combo_Box::~Moris_Bool_Combo_Box() = default;

    // Getter for the linked Parameter object.
    // Returns the Parameter currently associated with this widget.
    // Checks that the stored Parameter pointer is valid before dereferencing.
    // Inputs:
    // - None.
    // Outputs:
    // - Reference to the linked Parameter object.
    Parameter &Moris_Bool_Combo_Box::get_parameter()
    {
        MORIS_ERROR(mParameter, "Moris_Line_Edit::getParameter() called with null mParameter.");
        return *mParameter;
    }
    // Setter for the linked Parameter object.
    // Reassigns this widget to a new Parameter, refreshes the displayed value,
    // and reapplies the locked/read-only state based on the new Parameter.
    // Useful when widgets are dynamically created, reused, or rebound.
    // Inputs:
    // - a_parameter: Reference to the new Parameter object to link.
    // Outputs:
    // - None.
    void Moris_Bool_Combo_Box::setParameter(Parameter &parameter )
    {

        // qDebug() << "[IntSpinBox::setParameter]"
        //         << "widget =" << this
        //         << "name =" << objectName()
        //         << "old mParameter =" << mParameter
        //         << "new mParameter =" << &parameter;
        mParameter = &parameter;
        refreshDataParameter();

        if( mParameter && mParameter->is_locked() )
        {
            setDisabled(true);
        }
        else
        {
            setDisabled(false);
        }
    }

    // Refreshes the widget display from the linked Parameter object.
    // Reads the current value stored in the Parameter and updates the visible UI.
    // Signals should be blocked during this update to avoid triggering write-back
    // slots during initialization or rebinding.
    // Inputs:
    // - None.
    // Outputs:
    // - None.
    void Moris_Bool_Combo_Box::refreshDataParameter()
    {
        QSignalBlocker tBlocker( this );

        clear();
        addItem( "True" );
        addItem( "False");

        if (!mParameter)
        {
            setCurrentIndex( 1 );
            return;
        }
        bool tValue = false;

        if( mParameter->index() == variant_index< bool >() )
        {
            tValue = mParameter->get_value< bool >();
        }
        else // fallback for string like input, need to refine later for expected string inputs
        {
            std::string tValueString = mParameter->get_string();
            
            if ( tValueString == "true" || tValueString == "1" || tValueString == "True")
            {
                tValue = true;
            }
            else
            {
                tValue = false;
            }
        }
        // if tValue, set index 0, else 1
        setCurrentIndex( tValue ? 0 : 1 );

    }

    // Slot for handling user edits or selection changes.
    // Reads the current widget value, converts it to the correct Parameter type,
    // and writes the updated value back to the linked Parameter object.
    // Returns early if the Parameter is null, locked, or has an unsupported type.
    // Inputs:
    // - a_index: New index emitted by the widget signal.
    // Outputs:
    // - None.
    void Moris_Bool_Combo_Box::on_index_changed( int a_index )
    {
        // Update the parameter with the new value based on the selected index
        if(!mParameter || mParameter->is_locked() )
        {
            return;
        }
        if ( a_index == 0 )
        {
            mParameter->set_value( objectName().toStdString(), true, false );
        }
        else
        {
            mParameter->set_value( objectName().toStdString(), false, false );
        }

        // Emit the custom index_changed signal with the widget's name and the new index
        emit index_changed( objectName(), a_index );
    }

}    // namespace moris
