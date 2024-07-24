#include "moris_combo_box.hpp"

// Constructor for Moris_Combo_Box
// Initializes the combo box widget and sets up its items and signal-slot connections.
// Inputs:
// - parent: Pointer to the parent widget (default is nullptr).
// - parameter: Reference to a moris::Parameter object to be linked with this widget.
Moris_Combo_Box::Moris_Combo_Box( QWidget *parent, moris::Parameter &parameter )
        : QComboBox( parent )
        , mParameter( parameter )
{
    // Add items to the combo box from the parameter's selection names
    for ( const std::string &iSelectionOption : mParameter.get_selection_names() )
    {
        addItem( QString::fromStdString( iSelectionOption ) );
    }
    // Connect the currentIndexChanged(int) signal of QComboBox to the onIndexChanged slot
    connect( this, QOverload< int >::of( &QComboBox::currentIndexChanged ), this, &Moris_Combo_Box::onIndexChanged );
}

// Destructor for Moris_Combo_Box
// The destructor is defaulted as there are no specific cleanup requirements.
Moris_Combo_Box::~Moris_Combo_Box() = default;

// Getter for the associated moris::Parameter object
// Returns the reference to the parameter linked with this widget.
// Outputs:
// - Reference to the moris::Parameter object.
moris::Parameter &Moris_Combo_Box::getParameter()
{
    return mParameter;
}

// Slot to handle index changes in the combo box
// Updates the linked moris::Parameter object with the new value based on the selected index.
// Inputs:
// - index: The new index selected in the widget.
// Outputs:
// - None.
void Moris_Combo_Box::onIndexChanged( int index )
{
    // Update the parameter with the new value based on the selected index
    mParameter.set_value( objectName().toStdString(), currentText().toStdString(), false );

    // Emit the custom indexChanged signal with the widget's name and the new index
    emit indexChanged( objectName(), index );
}
