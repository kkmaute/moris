#ifndef TESTWINDOW_HPP
#define TESTWINDOW_HPP

#include <QMainWindow>
#include <QVBoxLayout>
#include <QPushButton>
#include <QMap>
#include <QDebug>
#include <variant>

#include "cl_line_edit.hpp"
#include "cl_int_spin_box.hpp"
#include "cl_double_spin_box.hpp"
#include "cl_combo_box.hpp"
#include "cl_pair_box.hpp"
#include "cl_Parameter_List.hpp"
#include "cl_XML_Parser.hpp"


namespace moris
{
    class TestWindow : public QMainWindow
    {
        Q_OBJECT

      public:
        explicit TestWindow( QWidget *parent, moris::Parameter_List &parameterList );

      private slots:
        void updateInput( const QString &name, const QString &text );
        void updateComboBox( const QString &name, int index );
        void updateDoubleSpinBox( const QString &name, double value );
        void updateIntSpinBox( const QString &name, int value );
        //void updatePairBox( const QString &name, const QString &text );
        void saveAndPrintInputs();
        void saveInputsToXML( const std::string &filePath );

      private:
        Moris_Line_Edit       *lineEdit;
        Moris_Int_Spin_Box    *intSpinBox;
        Moris_Double_Spin_Box *doubleSpinBox;
        Moris_Combo_Box       *comboBox;
        Moris_Pair_Box        *pairBox;

        moris::Parameter_List &parameterList;
    };
}    // namespace moris
#endif    // TESTWINDOW_HPP
