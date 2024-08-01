#include "cl_combo_box.hpp"

namespace moris
{
    // Constructor for Moris_Combo_Box
    // Initializes the combo box widget and sets up its items and signal-slot connections.
    // Inputs:
    // - a_parent: Pointer to the parent widget (default is nullptr).
    // - a_parameter: Reference to a moris::Parameter object to be linked with this widget.
    Moris_Combo_Box::Moris_Combo_Box( QWidget *a_parent, moris::Parameter &a_parameter )
            : QComboBox( a_parent )
            , m_parameter( a_parameter )
    {
        // Add items to the combo box from the parameter's selection names
        for ( const std::string &selection_option : m_parameter.get_selection_names() )
        {
            addItem( QString::fromStdString( selection_option ) );
        }
        // Connect the currentIndexChanged(int) signal of QComboBox to the on_index_changed slot
        connect( this, QOverload< int >::of( &QComboBox::currentIndexChanged ), this, &Moris_Combo_Box::on_index_changed );
    }

    // Destructor for Moris_Combo_Box
    // The destructor is defaulted as there are no specific cleanup requirements.
    Moris_Combo_Box::~Moris_Combo_Box() = default;

    // Getter for the associated moris::Parameter object
    // Returns the reference to the parameter linked with this widget.
    // Outputs:
    // - Reference to the moris::Parameter object.
    moris::Parameter &Moris_Combo_Box::get_parameter()
    {
        return m_parameter;
    }

    // Slot to handle index changes in the combo box
    // Updates the linked moris::Parameter object with the new value based on the selected index.
    // Inputs:
    // - a_index: The new index selected in the widget.
    // Outputs:
    // - None.
    void Moris_Combo_Box::on_index_changed( int a_index )
    {
        // Update the parameter with the new value based on the selected index
        m_parameter.set_value( objectName().toStdString(), currentText().toStdString(), false );

        // Emit the custom index_changed signal with the widget's name and the new index
        emit index_changed( objectName(), a_index );
    }
}    // namespace moris
