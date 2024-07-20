#include "moris_pair_box.hpp"

Moris_Pair_Box::Moris_Pair_Box(QWidget *parent, moris::Parameter &param, const QStringList &options)
    : QWidget(parent), morisPairComboBox(new QComboBox(this)), morisPairLineEdit(new QLineEdit(this)), parameter(param)
{
    // Set up the combo box with provided options
    morisPairComboBox->addItems(options);
    connect(morisPairComboBox, &QComboBox::currentTextChanged, this, &Moris_Pair_Box::onComboBoxTextChanged);

    // Set up the line edit
    connect(morisPairLineEdit, &QLineEdit::textChanged, this, &Moris_Pair_Box::onLineEditTextChanged);

    // Layout the widgets
    auto *layout = new QHBoxLayout(this);
    layout->addWidget(morisPairComboBox);
    layout->addWidget(morisPairLineEdit);
    setLayout(layout);
}

// // Public getter functions
// QComboBox* Moris_Pair_Box::getComboBox() const {
//     return morisPairComboBox;
// }

// QLineEdit* Moris_Pair_Box::getLineEdit() const {
//     return morisPairLineEdit;
// }

// Private slots for handling changes
void Moris_Pair_Box::onComboBoxTextChanged(const QString &text) {
    emit comboBoxTextChanged(objectName(), text);
}

void Moris_Pair_Box::onLineEditTextChanged(const QString &text) {
    emit lineEditTextChanged(objectName(), text);
}
