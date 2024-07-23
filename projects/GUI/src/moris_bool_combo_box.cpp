#include "moris_bool_combo_box.hpp"

// Constructor for Moris_Bool_Combo_Box
// Inputs:
// - parent: Pointer to the parent widget (default is nullptr).
// - parameter: Pointer to a moris::Parameter object to be linked with this widget (default is nullptr).
Moris_Bool_Combo_Box::Moris_Bool_Combo_Box( QWidget* parent, moris::Parameter& parameter )
        : QComboBox( parent )
        , mParameter( parameter )
{
    std::cout << "Construction index:" << mParameter.index() << std::endl;

    addItem( "true" );
    addItem( "false" );
    // Connect the currentIndexChanged(int) signal of QComboBox to the onIndexChanged slot
    connect( this, &QComboBox::currentIndexChanged, this, &Moris_Bool_Combo_Box::onIndexChanged );
}

// Destructor for Moris_Bool_Combo_Box
Moris_Bool_Combo_Box::~Moris_Bool_Combo_Box() = default;

// Slot to handle index changes
// This slot is connected to the currentIndexChanged(int) signal of QComboBox and updates the linked moris::Parameter object.
// Inputs:
// - index: New index selected in the widget.
void Moris_Bool_Combo_Box::onIndexChanged( int index )
{

    std::cout << "Index changed:" << mParameter.index() << std::endl;

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
