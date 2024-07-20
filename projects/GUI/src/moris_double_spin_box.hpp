#ifndef MORIS_DOUBLE_SPIN_BOX_HPP
#define MORIS_DOUBLE_SPIN_BOX_HPP

#include <QDoubleSpinBox>
#include "cl_Parameter.hpp"    

// Moris_Double_Spin_Box
// Custom QDoubleSpinBox widget for handling double value input and linking to moris::Parameter objects.
// This class extends QDoubleSpinBox to provide additional functionality for managing double input parameters.
// It emits a custom signal when the value changes, which includes the name and the new double value.
class Moris_Double_Spin_Box : public QDoubleSpinBox
{
    Q_OBJECT

  public:
    // Constructor for Moris_Double_Spin_Box.
    // Inputs:
    // - parent: Pointer to the parent widget (default is nullptr).
    // - parameter: Pointer to a moris::Parameter object to be linked with this widget (default is nullptr).
    explicit Moris_Double_Spin_Box( QWidget *parent, moris::Parameter &parameter );

    // Destructor for Moris_Double_Spin_Box.
    // Inputs:
    // - None.
    // Outputs:
    // - None.
    ~Moris_Double_Spin_Box() override;

    // Getter for the associated moris::Parameter object
    moris::Parameter &getParameter();

  signals:
    // Signal emitted when the value changes.
    // Inputs:
    // - name: Name associated with the widget.
    // - value: New double value input in the widget.
    void valueChanged( const QString &name, double value );

  private slots:
    // Slot to handle value changes.
    // This slot is connected to the valueChanged(double) signal of QDoubleSpinBox and updates the linked moris::Parameter object.
    // Inputs:
    // - value: New double value input in the widget.
    // Outputs:
    // - None.
    void onValueChanged( double value );

  private:
    moris::Parameter &mParameter;    // Pointer to the associated moris::Parameter object
};

#endif    // MORIS_DOUBLE_SPIN_BOX_HPP
