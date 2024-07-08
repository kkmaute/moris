#include "moris_combo_box.hpp"

// Constructor for Moris_Combo_Box
// Initializes the combo box and sets up the signal-slot connection for index changes
Moris_Combo_Box::Moris_Combo_Box(QWidget *parent)
    : QComboBox(parent)
{
    // Connect the QComboBox's currentIndexChanged signal to the onCurrentIndexChanged slot
    connect(this, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &Moris_Combo_Box::onCurrentIndexChanged);

    // Setting a parameter
    setParameter("Key", "Value");
}

// Destructor for Moris_Combo_Box
// Uses the default destructor
Moris_Combo_Box::~Moris_Combo_Box() = default;

// Sets a parameter key-value pair
// Inputs:
// - key: The parameter key (QString)
// - value: The parameter value (QVariant)
// Outputs: None
void Moris_Combo_Box::setParameter(const QString &key, const QVariant &value)
{
    parameters.insert(key, value);
}

// Retrieves the value for a given parameter key
// Inputs:
// - key: The parameter key (QString)
// Outputs:
// - Returns the value associated with the key (QVariant)
QVariant Moris_Combo_Box::parameter(const QString &key) const
{
    return parameters.value(key);
}

// Slot that handles the combo box's index change
// Inputs:
// - index: The new selected index (int)
// Outputs: None
void Moris_Combo_Box::onCurrentIndexChanged(int index)
{
    // Emit the custom currentIndexChanged signal with the object name and the new index
    emit currentIndexChanged(objectName(), index);
}
