#include "moris_int_spin_box.hpp"

// Constructor
Moris_Int_Spin_Box::Moris_Int_Spin_Box(QWidget *parent)
    : QSpinBox(parent)
{
    // Set the range for the spin box
    setRange(0, 100000); // Set the minimum to 0 and the maximum to 1000

    // Connect the valueChanged(int) signal to the onValueChanged(int) slot
    connect(this, QOverload<int>::of(&QSpinBox::valueChanged), this, &Moris_Int_Spin_Box::onValueChanged);
}

// Destructor
Moris_Int_Spin_Box::~Moris_Int_Spin_Box() = default;

// Sets a parameter key-value pair in the parameters map
void Moris_Int_Spin_Box::setParameter(const QString &key, const QVariant &value)
{
    parameters.insert(key, value);
}

// Retrieves the value for a given parameter key from the parameters map
QVariant Moris_Int_Spin_Box::parameter(const QString &key) const
{
    return parameters.value(key);
}

// Slot to handle value changes in the spin box
void Moris_Int_Spin_Box::onValueChanged(int new_value)
{
    emit valueChanged(objectName(), QVariant(new_value)); // Emit signal with new QVariant value
}
