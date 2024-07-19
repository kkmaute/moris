#ifndef MORIS_LINE_EDIT_HPP
#define MORIS_LINE_EDIT_HPP

#include <QLineEdit>
#include "cl_Parameter.hpp"

// Moris_Line_Edit
// Custom QLineEdit widget for handling text input and linking to moris::Parameter objects.
// This class extends QLineEdit to provide additional functionality for managing input parameters.
// It emits a custom signal when the text changes, which includes the name and the new text.
class Moris_Line_Edit : public QLineEdit
{
    Q_OBJECT

  public:
    // Constructor for Moris_Line_Edit.
    // Inputs:
    // - parent: Pointer to the parent widget (default is nullptr).
    // - parameter: Pointer to a moris::Parameter object to be linked with this widget (default is nullptr).
    explicit Moris_Line_Edit( QWidget *parent = nullptr, moris::Parameter *parameter = nullptr );

    // Destructor for Moris_Line_Edit.
    ~Moris_Line_Edit() override;

    // Getter for the associated moris::Parameter object
    moris::Parameter *getParameter() const;
    
    void setParameter( moris::Parameter *parameter );

  signals:
    // Signal emitted when the text changes.
    // Inputs:
    // - name: Name associated with the widget.
    // - new_text: New text input in the widget.
    void textChanged( const QString &name, const QString &new_text );
    void selectionChanged();

  private slots:
    // Slot to handle text changes.
    // This slot is connected to the textChanged signal of QLineEdit and updates the linked moris::Parameter object.
    // Inputs:
    // - new_text: New text input in the widget.
    void onTextChanged( const QString &new_text );
    //void onItemChanged();

  private:
    moris::Parameter *mParameter;    // Pointer to the associated moris::Parameter object
};

#endif    // MORIS_LINE_EDIT_HPP
