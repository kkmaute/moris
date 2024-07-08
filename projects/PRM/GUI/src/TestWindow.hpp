#ifndef TESTWINDOW_HPP
#define TESTWINDOW_HPP

#include "moris_line_edit.hpp"
#include "moris_combo_box.hpp"
#include "moris_double_spin_box.hpp"
#include "moris_int_spin_box.hpp" 
#include "cl_Parameter_List.hpp"
#include "cl_Library_IO.hpp"
#include "input_parameters.hpp"

#include <QMainWindow>
#include <QVBoxLayout>
#include <QPushButton>
#include <QString>

class TestWindow : public QMainWindow
{
    Q_OBJECT

public:
    // Constructor
    TestWindow(QWidget *parent = nullptr);

private slots:
    // Slots to handle user input changes
    void updateInput(const QString &name, const QString &text);
    void updateComboBox(const QString &name, int index);
    void updateDoubleSpinBox(const QString &name, const QVariant &value); // Updated to QVariant
    void updateIntSpinBox(const QString &name, const QVariant &value); // Slot for integer spin box
    void saveAndPrintInputs();

private:
    Input_Parameters currentInputs;         // Stores current user inputs
    moris::Parameter_List parameterList;    // Stores final parameters to save and print

    Moris_Combo_Box *comboBox;              // Pointer to the combo box instance
    QMap<QString, Moris_Line_Edit*> lineEdits; // Map to store line edits dynamically
    Moris_Double_Spin_Box *doubleSpinBox;   // Pointer to the double spin box instance
    Moris_Int_Spin_Box *intSpinBox;         // Pointer to the integer spin box instance
};

#endif // TESTWINDOW_HPP
