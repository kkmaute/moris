#pragma once

#include <QLineEdit>
#include "cl_Parameter.hpp"

namespace moris
{
    // Moris_Line_Edit
    // Custom QLineEdit widget for handling text input and linking to Parameter objects.
    // This class extends QLineEdit to provide additional functionality for managing text input parameters.
    // It emits custom signals when the text changes or when the selection changes.
    class Moris_Line_Edit : public QLineEdit
    {
        Q_OBJECT

      public:
        // Constructor for Moris_Line_Edit.
        // Initializes the line edit widget and links it with a Parameter object.
        // Inputs:
        // - parent: Pointer to the parent widget (default is nullptr).
        // - parameter: Reference to a Parameter object to be linked with this widget.
        explicit Moris_Line_Edit( QWidget *parent, Parameter &parameter );

        // Destructor for Moris_Line_Edit.
        // Defaulted as there are no specific cleanup requirements.
        ~Moris_Line_Edit() override;

        // Getter for the associated Parameter object.
        // Returns:
        // - Reference to the Parameter object linked with this widget.
        Parameter &getParameter();

        // Setter for the associated Parameter object.
        // Inputs:
        // - parameter: Reference to a Parameter object to be linked with this widget.
        void setParameter( Parameter &parameter );

      signals:
        // Signal emitted when the text changes.
        // Inputs:
        // - name: Name associated with the widget.
        // - new_text: New text input in the widget.
        void textChanged( const QString &name, const QString &new_text );

        // Signal emitted when the selection changes.
        // No additional inputs are required.
        void selectionChanged();

      private slots:
        // Slot to handle text changes.
        // This slot is connected to the textChanged signal of QLineEdit and updates the linked Parameter object.
        // Inputs:
        // - new_text: New text input in the widget.
        void onTextChanged( const QString &new_text );

        // Placeholder for a slot to handle item changes, if needed in the future.
        // void onItemChanged();

      private:
        Parameter &mParameter;    // Reference to the associated Parameter object
        template<typename T>
        void mTrySetParameter(const std::string& aParameterName, const T& aValue, const QString& aNewText);
    };

}    // namespace moris

