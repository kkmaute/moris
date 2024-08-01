#include "cl_bool_combo_box.hpp"

namespace moris
{

    // Constructor for Moris_Bool_Combo_Box
    // Initializes the combo box widget and sets up its items and signal-slot connections.
    // Inputs:
    // - a_parent: Pointer to the parent widget (default is nullptr).
    // - a_parameter: Reference to a moris::Parameter object to be linked with this widget.
    Moris_Bool_Combo_Box::Moris_Bool_Combo_Box( QWidget *a_parent, moris::Parameter &a_parameter )
            : QComboBox( a_parent )
            , m_Parameter( a_parameter )
    {
        // Add boolean options to the combo box
        addItem( "true" );
        addItem( "false" );

        // Connect the currentIndexChanged(int) signal of QComboBox to the on_index_changed slot
        connect( this, QOverload< int >::of( &QComboBox::currentIndexChanged ), this, &Moris_Bool_Combo_Box::on_index_changed );
    }

    // Destructor for Moris_Bool_Combo_Box
    // The destructor is defaulted as there are no specific cleanup requirements.
    Moris_Bool_Combo_Box::~Moris_Bool_Combo_Box() = default;

    // Getter for the associated moris::Parameter object
    // Returns the reference to the parameter linked with this widget.
    // Output:
    // - Reference to the moris::Parameter object.
    moris::Parameter &Moris_Bool_Combo_Box::get_parameter()
    {
        return m_Parameter;
    }

    // Slot to handle index changes in the combo box
    // Updates the linked moris::Parameter object with the new value based on the selected index.
    // Inputs:
    // - a_index: The new index selected in the widget.
    void Moris_Bool_Combo_Box::on_index_changed( int a_index )
    {
        // Debugging output to check the parameter index
        std::cout << m_Parameter.index() << std::endl;

        // Update the parameter with the new value based on the selected index
        if ( a_index == 0 )
        {
            m_Parameter.set_value( objectName().toStdString(), true, false );
        }
        else
        {
            m_Parameter.set_value( objectName().toStdString(), false, false );
        }

        // Emit the custom index_changed signal with the widget's name and the new index
        emit index_changed( objectName(), a_index );
    }

}    // namespace moris
