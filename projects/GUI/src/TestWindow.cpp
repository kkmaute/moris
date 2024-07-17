#include "TestWindow.hpp"

// Constructor for TestWindow
// Inputs:
// - parent: Pointer to the parent widget (default is nullptr).
// - parameterList: Reference to a moris::Parameter_List object containing parameters to be used in the UI.
TestWindow::TestWindow(QWidget *parent, const moris::Parameter_List &parameterList)
    : QMainWindow(parent), comboBox(nullptr)
{
    // Create central widget and layout for the main window
    auto *centralWidget = new QWidget(this);
    auto *layout = new QVBoxLayout(centralWidget);

    // Iterate through the parameter list to create UI elements
    for (auto it = parameterList.begin(); it != parameterList.end(); ++it)
    {
        const std::string &paramName = it->first;        // Parameter name
        const moris::Parameter &param = it->second;      // Corresponding moris::Parameter object

        if (paramName == "comboBox1")
        {
            // Create a combo box if the parameter name matches "comboBox1"
            comboBox = new Moris_Combo_Box(this, const_cast<moris::Parameter *>(&param));
            comboBox->setObjectName(QString::fromStdString(paramName));
            comboBox->addItem("Option 1");
            comboBox->addItem("Option 2");
            comboBox->addItem("Option 3");
            layout->addWidget(comboBox);

            // Connect combo box signal to updateComboBox slot
            connect(comboBox, &Moris_Combo_Box::indexChanged, this, &TestWindow::updateComboBox);
        }
        else if (paramName == "doubleSpinBox1")
        {
            // Create a double spin box if the parameter name matches "doubleSpinBox1"
            auto *doubleSpinBox = new Moris_Double_Spin_Box(this, const_cast<moris::Parameter *>(&param));
            doubleSpinBox->setObjectName(QString::fromStdString(paramName));
            layout->addWidget(doubleSpinBox);

            // Store the double spin box in a map for later access
            doubleSpinBoxes.insert(doubleSpinBox->objectName(), doubleSpinBox);
            // Connect double spin box signal to updateDoubleSpinBox slot
            connect(doubleSpinBox, &Moris_Double_Spin_Box::valueChanged, this, &TestWindow::updateDoubleSpinBox);
        }
        else if (paramName == "intSpinBox1")
        {
            // Create an integer spin box if the parameter name matches "intSpinBox1"
            auto *intSpinBox = new Moris_Int_Spin_Box(this, const_cast<moris::Parameter *>(&param));
            intSpinBox->setObjectName(QString::fromStdString(paramName));
            layout->addWidget(intSpinBox);

            // Store the integer spin box in a map for later access
            intSpinBoxes.insert(intSpinBox->objectName(), intSpinBox);
            // Connect integer spin box signal to updateIntSpinBox slot
            connect(intSpinBox, &Moris_Int_Spin_Box::valueChanged, this, &TestWindow::updateIntSpinBox);
        }
        else
        {
            // Create a line edit for other parameters
            auto *lineEdit = new Moris_Line_Edit(this, const_cast<moris::Parameter *>(&param));
            lineEdit->setObjectName(QString::fromStdString(paramName));
            layout->addWidget(lineEdit);

            // Store the line edit in a map for later access
            lineEdits.insert(lineEdit->objectName(), lineEdit);
            // Connect line edit signal to updateInput slot
            connect(lineEdit, &Moris_Line_Edit::textChanged, this, &TestWindow::updateInput);
        }
    }

    // Create a button to save and print inputs
    auto *printButton = new QPushButton("Save and Print Inputs", this);
    layout->addWidget(printButton);

    // Connect button click signal to saveAndPrintInputs slot
    connect(printButton, &QPushButton::clicked, this, &TestWindow::saveAndPrintInputs);

    // Set the central widget for the main window
    setCentralWidget(centralWidget);
}

// Slot to update text of a LineEdit widget based on its name
// Inputs:
// - name: Name of the LineEdit widget.
// - text: New text input in the LineEdit widget.
void TestWindow::updateInput(const QString &name, const QString &text)
{
    auto lineEdit = lineEdits.value(name);
    if (lineEdit)
    {
        lineEdit->setText(text);
    }
}

// Slot to update the current index of the combo box widget based on its name
// Inputs:
// - name: Name of the combo box widget.
// - index: New index selected in the combo box.
void TestWindow::updateComboBox(const QString &name, int index)
{
    if (comboBox && comboBox->objectName() == name)
    {
        comboBox->setCurrentIndex(index);
    }
}

// Slot to update the value of a double spin box widget based on its name
// Inputs:
// - name: Name of the double spin box widget.
// - value: New value set in the double spin box.
void TestWindow::updateDoubleSpinBox(const QString &name, double value)
{
    auto doubleSpinBox = doubleSpinBoxes.value(name);
    if (doubleSpinBox)
    {
        doubleSpinBox->setValue(value);
    }
}

// Slot to update the value of an integer spin box widget based on its name
// Inputs:
// - name: Name of the integer spin box widget.
// - value: New value set in the integer spin box.
void TestWindow::updateIntSpinBox(const QString &name, int value)
{
    auto intSpinBox = intSpinBoxes.value(name);
    if (intSpinBox)
    {
        intSpinBox->setValue(value);
    }
}

// Slot to save and print input values of all widgets
void TestWindow::saveAndPrintInputs()
{
    qDebug() << "Saved Input Set:";

    // Print inputs for LineEdits
    for (auto it = lineEdits.begin(); it != lineEdits.end(); ++it)
    {
        const QString &name = it.key();
        const QString &text = it.value()->text();
        qDebug() << name << ": " << text;
    }

    // Print input for ComboBox
    if (comboBox)
    {
        const QString &name = comboBox->objectName();
        const int index = comboBox->currentIndex();
        qDebug() << name << ": " << index;
    }

    // Print inputs for DoubleSpinBoxes
    for (auto it = doubleSpinBoxes.begin(); it != doubleSpinBoxes.end(); ++it)
    {
        const QString &name = it.key();
        const double value = it.value()->value();
        qDebug() << name << ": " << value;
    }

    // Print inputs for IntSpinBoxes
    for (auto it = intSpinBoxes.begin(); it != intSpinBoxes.end(); ++it)
    {
        const QString &name = it.key();
        const int value = it.value()->value();
        qDebug() << name << ": " << value;
    }
}
