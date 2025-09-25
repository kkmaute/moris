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
            , m_parameter( a_parameter )
    {
        // Add items to the combo box from the parameter's selection names
        if(!m_parameter.is_locked()) {
            if(m_parameter.get_entry_type() == Entry_Type::SELECTION) {
                for ( const std::string &selection_option : m_parameter.get_selection_names() )
                {
                    addItem( QString::fromStdString( selection_option ) );
                }
            }
        }   else {
            addItem(QString::fromStdString(m_parameter.get_string()));
        }
        this->blockSignals(true);
        if ( m_parameter.index() == variant_index< uint >() )
        {
            setCurrentIndex( m_parameter.get_value< uint >() );
        }
        this->blockSignals(false);
        // Connect the currentIndexChanged(int) signal of QComboBox to the on_index_changed slot
        if ( m_parameter.is_locked() )
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
            , m_parameter( a_parameter )
            , m_options( a_options )
    {
        // Set up the combo box with the provided options
        set_options_list( a_options );
        if ( m_parameter.index() == variant_index< uint >() )
        {
            this->blockSignals(true);
            setCurrentIndex( m_parameter.get_value< uint >() );
            this->blockSignals(false);
        }
        // Connect the currentIndexChanged(int) signal of QComboBox to the on_index_changed slot
        if ( m_parameter.is_locked() )
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
        return m_parameter;
    }

    // Slot to handle index changes in the combo box
    // Updates the linked Parameter object with the new value based on the selected index.
    // Inputs:
    // - a_index: The new index selected in the widget.
    // Outputs:
    // - None.
    void Moris_Combo_Box::on_index_changed( int a_index )
    {
        if(m_parameter.index() == variant_index<uint>()){
            if(static_cast< uint >(a_index) == m_parameter.get_value< uint >()) return;
        } else if(m_parameter.index() == variant_index<std::string>()){
            if(currentText().toStdString() == m_parameter.get_value<std::string>()) return;
        }

        
        // Update the parameter with the new value based on the selected index

        std::cout << "objectName" << objectName().toStdString() << "\n";
        std::cout << "currentText" << currentText().toStdString();
        m_parameter.set_value( objectName().toStdString(), currentText().toStdString(), false );

        // Emit the custom index_changed signal with the widget's name and the new index
        emit index_changed( objectName(), a_index );
    }
}    // namespace moris
