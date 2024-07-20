#ifndef TESTWINDOW_HPP
#define TESTWINDOW_HPP

#include <QMainWindow>
#include <QVBoxLayout>
#include <QPushButton>
#include <QMap>
#include <QDebug>
#include <variant>

#include "moris_line_edit.hpp"
#include "moris_int_spin_box.hpp"
#include "moris_double_spin_box.hpp"
#include "moris_combo_box.hpp"
#include "moris_pair_box.hpp"
#include "cl_Parameter_List.hpp"

class TestWindow : public QMainWindow
{
    Q_OBJECT

  public:
    explicit TestWindow( QWidget *parent, moris::Parameter_List &parameterList);

  private slots:
    void updateInput( const QString &name, const QString &text );
    void updateComboBox( const QString &name, int index );
    void updateDoubleSpinBox( const QString &name, double value );
    void updateIntSpinBox( const QString &name, int value );
    void updatePairBox( const QString &name, const QString &text );
    void saveAndPrintInputs();

  private:
    Moris_Line_Edit       *lineEdit;
    Moris_Int_Spin_Box    *intSpinBox;
    Moris_Double_Spin_Box *doubleSpinBox;
    Moris_Combo_Box       *comboBox;
    Moris_Pair_Box        *pairBox;

    QMap< QString, Moris_Line_Edit * >       lineEdits;
    QMap< QString, Moris_Int_Spin_Box * >    intSpinBoxes;
    QMap< QString, Moris_Double_Spin_Box * > doubleSpinBoxes;
    QMap< QString, Moris_Pair_Box * >        pairBoxes;
};

#endif    // TESTWINDOW_HPP
