#ifndef TESTWINDOW_HPP
#define TESTWINDOW_HPP

#include <QMainWindow>
#include <QVBoxLayout>
#include <QPushButton>
#include <QDebug>

#include "cl_Parameter_List.hpp"
#include "moris_line_edit.hpp"
#include "moris_combo_box.hpp"
#include "moris_double_spin_box.hpp"
#include "moris_int_spin_box.hpp"

class Moris_Line_Edit;
class Moris_Combo_Box;
class Moris_Double_Spin_Box;
class Moris_Int_Spin_Box;

class TestWindow : public QMainWindow
{
    Q_OBJECT

  public:
    // Constructor
    TestWindow( QWidget *parent = nullptr, const moris::Parameter_List &parameterList = {} );

  private slots:
    // Slots for updating input values
    void updateInput( const QString &name, const QString &text );
    void updateComboBox( const QString &name, int index );
    void updateDoubleSpinBox( const QString &name, double value );
    void updateIntSpinBox( const QString &name, int value );

    // Slot for saving and printing inputs
    void saveAndPrintInputs();

  private:
    Moris_Combo_Box                         *comboBox = nullptr;    // Pointer to ComboBox widget
    QMap< QString, Moris_Line_Edit * >       lineEdits;             // Map to store LineEdit widgets
    QMap< QString, Moris_Double_Spin_Box * > doubleSpinBoxes;       // Map to store DoubleSpinBox widgets
    QMap< QString, Moris_Int_Spin_Box * >    intSpinBoxes;          // Map to store IntSpinBox widgets
};

#endif    // TESTWINDOW_HPP
