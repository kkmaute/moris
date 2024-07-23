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
    // Constructor for Moris_Bool_Combo_Box.
    // Inputs:
    // - parent: Pointer to the parent widget (default is nullptr).
    // - parameter: Pointer to a moris::Parameter object to be linked with this widget (default is nullptr).
    explicit Moris_Bool_Combo_Box( QWidget* parent, moris::Parameter& parameter );

    // Destructor for Moris_Bool_Combo_Box.
    // Inputs:
    // - None.
    // Outputs:
    // - None.
    ~Moris_Bool_Combo_Box() override;

  signals:
    // Signal emitted when the index changes.
    // Inputs:
    // - name: Name associated with the widget.
    // - index: New index selected in the widget.
    void indexChanged( const QString &name, int index );

  private slots:
    // Slot to handle index changes.
    // This slot is connected to the currentIndexChanged(int) signal of QComboBox and updates the linked moris::Parameter object.
    // Inputs:
    // - index: New index selected in the widget.
    // Outputs:
    // - None.
    void onIndexChanged( int index );

  private:
    moris::Parameter &mParameter;    // Pointer to the associated moris::Parameter object
};

#endif    // Moris_Bool_Combo_Box_HPP
