#include "moris_bool_combo_box.hpp"

// Constructor for Moris_Bool_Combo_Box
// Initializes the combo box widget and sets up its items and signal-slot connections.
// Inputs:
// - parent: Pointer to the parent widget (default is nullptr).
// - parameter: Reference to a moris::Parameter object to be linked with this widget.
Moris_Bool_Combo_Box::Moris_Bool_Combo_Box( QWidget *parent, moris::Parameter &parameter )
        : QComboBox( parent )
        , mParameter( parameter )
{
    
    // Add boolean options to the combo box
    addItem( "true" );
    addItem( "false" );

    // Connect the currentIndexChanged(int) signal of QComboBox to the onIndexChanged slot
    connect( this, QOverload< int >::of( &QComboBox::currentIndexChanged ), this, &Moris_Bool_Combo_Box::onIndexChanged );
}

// Destructor for Moris_Bool_Combo_Box
// The destructor is defaulted as there are no specific cleanup requirements.
Moris_Bool_Combo_Box::~Moris_Bool_Combo_Box() = default;

// Getter for the associated moris::Parameter object
// Returns the reference to the parameter linked with this widget.
// Output:
// - Reference to the moris::Parameter object.
moris::Parameter &Moris_Bool_Combo_Box::getParameter()
{
    return mParameter;
}

// Slot to handle index changes in the combo box
// Updates the linked moris::Parameter object with the new value based on the selected index.
// Inputs:
// - index: The new index selected in the widget.
void Moris_Bool_Combo_Box::onIndexChanged( int index )
{
    // Debugging output to check the parameter index
    std::cout << mParameter.index() << std::endl;

    // Update the parameter with the new value based on the selected index
    if ( index == 0 )
    {
        mParameter.set_value( objectName().toStdString(), true, false );
    }
    else
    {
        mParameter.set_value( objectName().toStdString(), false, false );
    }

    // Emit the custom indexChanged signal with the widget's name and the new index
    emit indexChanged( objectName(), index );
}
