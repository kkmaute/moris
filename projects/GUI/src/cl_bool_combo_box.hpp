#pragma once

#include <QComboBox>
#include "cl_Parameter.hpp"

namespace moris
{

    // Moris_Bool_Combo_Box
    // Custom QComboBox widget for handling index changes and linking to Parameter objects.
    // This class extends QComboBox to provide additional functionality for managing combo box index parameters.
    // It emits a custom signal when the index changes, which includes the name and the new index value.
    class Moris_Bool_Combo_Box : public QComboBox
    {
        Q_OBJECT

      public:
        // Constructor for Moris_Bool_Combo_Box
        // Initializes the combo box widget and sets up its items and signal-slot connections.
        // Inputs:
        // - a_parent: Pointer to the parent widget (default is nullptr).
        // - a_parameter: Reference to a Parameter object to be linked with this widget.
        explicit Moris_Bool_Combo_Box( QWidget *a_parent, Parameter &a_parameter );

        // Destructor for Moris_Bool_Combo_Box
        // The destructor is defaulted as there are no specific cleanup requirements.
        // Inputs:
        // - None.
        // Outputs:
        // - None.
        ~Moris_Bool_Combo_Box() override;

        // Getter for the associated Parameter object
        // Returns the reference to the parameter linked with this widget.
        // Outputs:
        // - Reference to the Parameter object.
        Parameter &get_parameter();

      signals:
        // Signal emitted when the index changes
        // Inputs:
        // - a_name: Name associated with the widget.
        // - a_index: New index selected in the widget.
        void index_changed( const QString &a_name, int a_index );

      private slots:
        // Slot to handle index changes
        // This slot is connected to the currentIndexChanged(int) signal of QComboBox and updates the linked Parameter object.
        // Inputs:
        // - a_index: New index selected in the widget.
        // Outputs:
        // - None.
        void on_index_changed( int a_index );

      private:
        Parameter &m_Parameter;    // Reference to the associated Parameter object
    };

}    // namespace moris
