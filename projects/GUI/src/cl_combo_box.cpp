#include "cl_combo_box.hpp"

namespace moris
{
    // Constructor for Moris_Combo_Box
    // Initializes the combo box widget and sets up its items and signal-slot connections.
    // Inputs:
    // - a_parent: Pointer to the parent widget (default is nullptr).
    // - a_parameter: Reference to a Parameter object to be linked with this widget.
    Moris_Combo_Box::Moris_Combo_Box( QWidget *a_parent, Parameter &a_parameter )
            : QComboBox( a_parent )
            , mParameter( &a_parameter )
    {
        // Add items to the combo box from the parameter's selection names
        if(!mParameter->is_locked()) 
        {
            if(mParameter->get_entry_type() == Entry_Type::SELECTION) 
            {
                for ( const std::string &selection_option : mParameter->get_selection_names() )
                {
                    addItem( QString::fromStdString( selection_option ) );
                }
            }
        }   
        else 
        {
            addItem(QString::fromStdString(mParameter->get_string()));
        }
        this->blockSignals(true);

        if ( mParameter->index() == variant_index< uint >() )
        {
            setCurrentIndex( mParameter->get_value< uint >() );
        }
        this->blockSignals(false);
        
        // Connect the currentIndexChanged(int) signal of QComboBox to the on_index_changed slot
        if ( mParameter->is_locked() )
        {
            setDisabled( true );
        }
        else
        {
            connect( this, QOverload< int >::of( &QComboBox::currentIndexChanged ), this, &Moris_Combo_Box::on_index_changed );
        }
    }

    // Overload constructor that gives the option to set the combo box items
    // Inputs:
    // - a_parent: Pointer to the parent widget (default is nullptr).
    // - a_parameter: Reference to a Parameter object to be linked with this widget.
    // - a_options: QStringList containing the options to be set in the combo box.
    Moris_Combo_Box::Moris_Combo_Box( QWidget *a_parent, Parameter &a_parameter, QStringList &a_options )
            : QComboBox( a_parent )
            , mParameter( &a_parameter )
            , m_options( a_options )
    {
        // Set up the combo box with the provided options
        set_options_list( a_options );
        if ( mParameter->index() == variant_index< uint >() )
        {
            this->blockSignals(true);
            setCurrentIndex( mParameter->get_value< uint >() );
            this->blockSignals(false);
        }
        // Connect the currentIndexChanged(int) signal of QComboBox to the on_index_changed slot
        if ( mParameter->is_locked() )
        {
            setDisabled( true );
        }
        else
        {
            connect( this, QOverload< int >::of( &QComboBox::currentIndexChanged ), this, &Moris_Combo_Box::on_index_changed );
        }
    }

    // Destructor for Moris_Combo_Box
    // The destructor is defaulted as there are no specific cleanup requirements.
    Moris_Combo_Box::~Moris_Combo_Box() = default;

    // Getter for the associated Parameter object
    // Returns the reference to the parameter linked with this widget.
    // Outputs:
    // - Reference to the Parameter object.
    Parameter &Moris_Combo_Box::get_parameter()
    {
        MORIS_ERROR(mParameter, "Moris_Line_Edit::getParameter() called with null mParameter.");
        return *mParameter;
    }


    void Moris_Combo_Box::setParameter(Parameter &parameter )
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



    // Slot to handle index changes in the combo box
    // Updates the linked Parameter object with the new value based on the selected index.
    // Inputs:
    // - a_index: The new index selected in the widget.
    // Outputs:
    // - None.
    void Moris_Combo_Box::on_index_changed( int a_index )
    {
        if(mParameter->index() == variant_index<uint>()){
            if(static_cast< uint >(a_index) == mParameter->get_value< uint >()) return;
        } else if(mParameter->index() == variant_index<std::string>()){
            if(currentText().toStdString() == mParameter->get_value<std::string>()) return;
        }

        
        // Update the parameter with the new value based on the selected index

        mParameter->set_value( objectName().toStdString(), currentText().toStdString(), false );

        // Emit the custom index_changed signal with the widget's name and the new index
        emit index_changed( objectName(), a_index );
    }


    // Refreshes the widget display from the linked Parameter object.
    // Reads the current value stored in the Parameter and updates the visible UI.
    // Signals should be blocked during this update to avoid triggering write-back
    // slots during initialization or rebinding.
    // Inputs:
    // - None.
    // Outputs:
    // - None.
    void Moris_Combo_Box::refreshDataParameter()
    {
        QSignalBlocker tBlocker( this ); // block signals to prevent unwanted signal emissions during refresh

        if ( !mParameter )
        {
            clear();
            return;
        }
        // if locked, show only the current stored value and disable further interaction
        if(mParameter->is_locked())
        {
            clear();
            addItem( QString::fromStdString(mParameter->get_string()) );
            setCurrentIndex( 0 );
            return;
        }
        // rebuild selection options if this parameter is a selection entry
        if ( mParameter->get_entry_type() == Entry_Type::SELECTION )
        {
            clear();
            for ( const std::string &selection_option : mParameter->get_selection_names() )
            {
                addItem(QString::fromStdString(selection_option));
            }
        }

        // sync the current selected value
        if (mParameter->index() == variant_index< uint >())
        {
            uint tIndex = mParameter->get_value< uint >();
            if ( tIndex < (uint)this->count() )
            {
                setCurrentIndex( static_cast<int>(tIndex) );
            }
            else
            {
                setCurrentIndex( 0 );
            }
        }
        else if ( mParameter->index() == variant_index< std::string >() )
        {
            QString tValue = QString::fromStdString( mParameter->get_value< std::string >() );
            int tIndex = findText( tValue );
            if ( tIndex >= 0 )
            {
                setCurrentIndex( tIndex );
            }
            else
            {
                // if current value is not in box, add it
                addItem ( tValue );
                setCurrentIndex( count() - 1 );
            }
        }
        else 
        {
             QString tValue = QString::fromStdString( mParameter->get_string() );
             int tIndex = findText( tValue );

            if ( tIndex >= 0 )
            {
                setCurrentIndex( tIndex );
            }
            else
            {
                // if current value is not in box, add it
                addItem ( tValue );
                setCurrentIndex( count() - 1 );
            }
        }
    }
}    // namespace moris