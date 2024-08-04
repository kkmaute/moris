#pragma once

#include <QSpinBox>
#include "cl_Parameter.hpp"

namespace moris
{
    // Moris_Int_Spin_Box
    // Custom QSpinBox widget for handling integer value input and linking to Parameter objects.
    // This class extends QSpinBox to provide additional functionality for managing integer input parameters.
    // It emits a custom signal when the value changes, including the widget's name and the new integer value.
    class Moris_Int_Spin_Box : public QSpinBox
    {
        Q_OBJECT

      public:
        // Constructor for Moris_Int_Spin_Box.
        // Inputs:
        // - a_parent: Pointer to the parent widget (default is nullptr).
        // - a_parameter: Reference to a Parameter object to be linked with this widget.
        explicit Moris_Int_Spin_Box( QWidget *a_parent, Parameter &a_parameter );

        // Destructor for Moris_Int_Spin_Box.
        // No special cleanup is required, so the destructor is defaulted.
        ~Moris_Int_Spin_Box() override;

        // Getter for the associated Parameter object.
        // Returns:
        // - Reference to the Parameter object linked with this widget.
        Parameter &get_parameter();

      signals:
        // Signal emitted when the value changes.
        // Inputs:
        // - a_name: Name associated with the widget.
        // - a_value: New integer value input in the widget.
        void value_changed( const QString &a_name, int a_value );

      private slots:
        // Slot to handle value changes.
        // This slot is connected to the valueChanged(int) signal of QSpinBox and updates the linked Parameter object.
        // Inputs:
        // - a_value: New integer value input in the widget.
        void on_value_changed( int a_value );

      private:
        Parameter &m_parameter;    // Reference to the associated Parameter object
    };

}    // namespace moris