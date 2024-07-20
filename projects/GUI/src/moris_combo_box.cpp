#include "moris_combo_box.hpp"

// Constructor for Moris_Combo_Box
// Inputs:
// - parent: Pointer to the parent widget (default is nullptr).
// - parameter: Pointer to a moris::Parameter object to be linked with this widget (default is nullptr).
Moris_Combo_Box::Moris_Combo_Box( QWidget *parent, moris::Parameter &parameter )
        : QComboBox( parent )
        , mParameter( parameter )
{
    // Connect the currentIndexChanged(int) signal of QComboBox to the onIndexChanged slot
    connect( this, QOverload< int >::of( &QComboBox::currentIndexChanged ), this, &Moris_Combo_Box::onIndexChanged );

    // If parameter is not null, set the initial index from the parameter value

        setCurrentIndex( mParameter.get_value< uint >() );

}

// Destructor for Moris_Combo_Box
Moris_Combo_Box::~Moris_Combo_Box() = default;

// Getter for the associated moris::Parameter object
moris::Parameter &Moris_Combo_Box::getParameter()
{
    return mParameter;
}

// Slot to handle index changes
// This slot is connected to the currentIndexChanged(int) signal of QComboBox and updates the linked moris::Parameter object.
// Inputs:
// - index: New index selected in the widget.
void Moris_Combo_Box::onIndexChanged( int index )
{
    // If parameter is not null, update the parameter value with the new index
        mParameter.set_value( objectName().toStdString(), index, false );

    // Emit the custom indexChanged signal with the widget's name and the new index
    emit indexChanged( objectName(), index );
}
