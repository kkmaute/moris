#ifndef CL_PAIR_BOX_HPP
#define CL_PAIR_BOX_HPP

#include "cl_Parameter.hpp"
#include "cl_Variant.hpp"

#include <QWidget>
#include <QComboBox>
#include <QLineEdit>
#include <QStringList>
#include <QHBoxLayout>

namespace moris
{
    class Parameter;
}

namespace moris
{
    // Moris_Pair_Box
    // Custom widget that contains a QComboBox and a QLineEdit.
    // It links to a Parameter object and emits signals when the text or selection changes.
    class Moris_Pair_Box : public QWidget
    {
        Q_OBJECT

      public:
        // Constructor
        // Inputs:
        // - a_parent: Pointer to the parent widget (default is nullptr).
        // - a_param: Reference to a Parameter object to be linked with this widget.
        // - a_options: QStringList containing options for the combo box (default is an empty list).

        // explicit Moris_Pair_Box( QWidget *a_parent, Parameter &a_param, const QStringList &a_options = {} );
        explicit Moris_Pair_Box( QWidget *a_parent, Parameter &a_param );

        // Public member variables
        // QComboBox *moris_pair_combo_box;
        QLineEdit *moris_pair_line_edit;
        QLineEdit *moris_pair_line_edit_2;

      public slots:

        // Slot to handle text changes in the line edit
        // Inputs:
        // - a_text: New text entered in the line edit.
        void on_line_edit_text_changed( const QString &a_text );

        // Slot to set the associated Parameter with the combined value of the combo box and line edit
        // void set_parameter();

      private:
        Parameter                            &mParameter;    // Reference to the Parameter object
        std::pair< std::string, std::string > mPairValue;
    };

}    // namespace moris

#endif    // CL_PAIR_BOX_HPP
