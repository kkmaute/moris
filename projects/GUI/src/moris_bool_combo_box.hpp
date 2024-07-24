#ifndef MORIS_BOOL_COMBO_BOX_HPP
#define MORIS_BOOL_COMBO_BOX_HPP

#include <QComboBox>
#include "cl_Parameter.hpp"

// Moris_Bool_Combo_Box
// Custom QComboBox widget for handling index changes and linking to moris::Parameter objects.
// This class extends QComboBox to provide additional functionality for managing combo box index parameters.
// It emits a custom signal when the index changes, which includes the name and the new index value.
class Moris_Bool_Combo_Box : public QComboBox
{
    Q_OBJECT

  public:
    // Constructor for Moris_Bool_Combo_Box
    // Initializes the combo box widget and sets up its items and signal-slot connections.
    // Inputs:
    // - parent: Pointer to the parent widget (default is nullptr).
    // - parameter: Reference to a moris::Parameter object to be linked with this widget.
    explicit Moris_Bool_Combo_Box( QWidget *parent, moris::Parameter &parameter );

    // Destructor for Moris_Bool_Combo_Box
    // The destructor is defaulted as there are no specific cleanup requirements.
    // Inputs:
    // - None.
    // Outputs:
    // - None.
    ~Moris_Bool_Combo_Box() override;

    // Getter for the associated moris::Parameter object
    // Returns the reference to the parameter linked with this widget.
    // Outputs:
    // - Reference to the moris::Parameter object.
    moris::Parameter &getParameter();

  signals:
    // Signal emitted when the index changes
    // Inputs:
    // - name: Name associated with the widget.
    // - index: New index selected in the widget.
    void indexChanged( const QString &name, int index );

  private slots:
    // Slot to handle index changes
    // This slot is connected to the currentIndexChanged(int) signal of QComboBox and updates the linked moris::Parameter object.
    // Inputs:
    // - index: New index selected in the widget.
    // Outputs:
    // - None.
    void onIndexChanged( int index );

  private:
    moris::Parameter &mParameter;    // Reference to the associated moris::Parameter object
};

#endif    // MORIS_BOOL_COMBO_BOX_HPP
