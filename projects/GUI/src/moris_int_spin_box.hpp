#ifndef MORIS_INT_SPIN_BOX_HPP
#define MORIS_INT_SPIN_BOX_HPP

#include <QSpinBox>
#include "cl_Parameter.hpp"

// Moris_Int_Spin_Box
// Custom QSpinBox widget for handling integer value input and linking to moris::Parameter objects.
// This class extends QSpinBox to provide additional functionality for managing integer input parameters.
// It emits a custom signal when the value changes, which includes the name and the new integer value.
class Moris_Int_Spin_Box : public QSpinBox
{
    Q_OBJECT

  public:
    // Constructor for Moris_Int_Spin_Box.
    // Inputs:
    // - parent: Pointer to the parent widget (default is nullptr).
    // - parameter: Pointer to a moris::Parameter object to be linked with this widget (default is nullptr).
    explicit Moris_Int_Spin_Box( QWidget *parent, moris::Parameter &parameter );

    // Destructor for Moris_Int_Spin_Box.
    // Inputs:
    // - None.
    // Outputs:
    // - None.
    ~Moris_Int_Spin_Box() override;

    // Getter for the associated moris::Parameter object
    moris::Parameter &getParameter() ;

  signals:
    // Signal emitted when the value changes.
    // Inputs:
    // - name: Name associated with the widget.
    // - value: New integer value input in the widget.
    void valueChanged( const QString &name, int value );

  private slots:
    // Slot to handle value changes.
    // This slot is connected to the valueChanged(int) signal of QSpinBox and updates the linked moris::Parameter object.
    // Inputs:
    // - value: New integer value input in the widget.
    // Outputs:
    // - None.
    void onValueChanged( int value );

  private:
    moris::Parameter &mParameter;    // Pointer to the associated moris::Parameter object
};

#endif    // MORIS_INT_SPIN_BOX_HPP
