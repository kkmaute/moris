#ifndef MORIS_PAIR_BOX_HPP
#define MORIS_PAIR_BOX_HPP

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

// Moris_Pair_Box
// Custom widget that contains a QComboBox and a QLineEdit.
// It links to a moris::Parameter object and emits signals when the text or selection changes.
class Moris_Pair_Box : public QWidget
{
    Q_OBJECT

  public:
    // Constructor
    // Inputs:
    // - parent: Pointer to the parent widget (default is nullptr).
    // - param: Reference to a moris::Parameter object to be linked with this widget.
    // - options: QStringList containing options for the combo box (default is an empty list).
    explicit Moris_Pair_Box( QWidget *parent, moris::Parameter &param, const QStringList &options = {} );

    // Public member variables
    QComboBox *morisPairComboBox;
    QLineEdit *morisPairLineEdit;

  signals:
    // Signal emitted when the text in the combo box changes
    // Inputs:
    // - name: Name associated with the widget.
    // - text: New text selected in the combo box.
    void comboBoxTextChanged( const QString &name, const QString &text );

    // Signal emitted when the text in the line edit changes
    // Inputs:
    // - name: Name associated with the widget.
    // - text: New text entered in the line edit.
    void lineEditTextChanged( const QString &name, const QString &text );

  public slots:
    // Slot to handle text changes in the combo box
    // Inputs:
    // - text: New text selected in the combo box.
    void onComboBoxTextChanged( const QString &text );

    // Slot to handle text changes in the line edit
    // Inputs:
    // - text: New text entered in the line edit.
    void onLineEditTextChanged( const QString &text );

    // Slot to set the associated moris::Parameter with the combined value of the combo box and line edit
    void setParameter();

  private:
    moris::Parameter &parameter;    // Reference to the moris::Parameter object
};

#endif    // MORIS_PAIR_BOX_HPP
