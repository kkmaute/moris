#ifndef CL_COMBO_BOX_HPP
#define CL_COMBO_BOX_HPP

#include <QComboBox>
#include "cl_Parameter.hpp"

namespace moris
{
    // Moris_Combo_Box
    // Custom QComboBox widget for handling index changes and linking to moris::Parameter objects.
    // This class extends QComboBox to provide additional functionality for managing combo box index parameters.
    // It emits a custom signal when the index changes, which includes the name and the new index value.
    class Moris_Combo_Box : public QComboBox
    {
        Q_OBJECT

      public:
        // Constructor for Moris_Combo_Box.
        // Initializes the combo box widget and sets up its items and signal-slot connections.
        // Inputs:
        // - a_parent: Pointer to the parent widget (default is nullptr).
        // - a_parameter: Reference to a moris::Parameter object to be linked with this widget.
        explicit Moris_Combo_Box( QWidget *a_parent, moris::Parameter &a_parameter );

        // Destructor for Moris_Combo_Box.
        // The destructor is defaulted as there are no specific cleanup requirements.
        // Inputs:
        // - None.
        // Outputs:
        // - None.
        ~Moris_Combo_Box() override;

        // Getter for the associated moris::Parameter object.
        // Returns the reference to the parameter linked with this widget.
        // Outputs:
        // - Reference to the moris::Parameter object.
        moris::Parameter &get_parameter();

      signals:
        // Signal emitted when the index changes.
        // Inputs:
        // - a_name: Name associated with the widget.
        // - a_index: New index selected in the widget.
        void index_changed( const QString &a_name, int a_index );

      private slots:
        // Slot to handle index changes.
        // This slot is connected to the currentIndexChanged(int) signal of QComboBox and updates the linked moris::Parameter object.
        // Inputs:
        // - a_index: New index selected in the widget.
        // Outputs:
        // - None.
        void on_index_changed( int a_index );

      private:
        moris::Parameter &m_parameter;    // Reference to the associated moris::Parameter object
    };

}    // namespace moris
#endif    // CL_COMBO_BOX_HPP