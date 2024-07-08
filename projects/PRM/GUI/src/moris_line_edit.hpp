#ifndef MORIS_LINE_EDIT_HPP
#define MORIS_LINE_EDIT_HPP

#include <QLineEdit>    
#include <QVariant>     
#include <QMap>         
#include <QDebug>     

// Moris_Line_Edit class declaration
// Inherits from QLineEdit
class Moris_Line_Edit : public QLineEdit
{
    Q_OBJECT

public:
    // Constructor
    // Initializes the line edit widget with the given parent widget
    // Inputs:
    // - parent: Pointer to the parent widget, default is nullptr
    // Outputs: None
    explicit Moris_Line_Edit(QWidget *parent = nullptr);

    // Destructor
    // Outputs: None
    ~Moris_Line_Edit();

    // Sets a parameter key-value pair
    // Inputs:
    // - key: The parameter key (QString)
    // - value: The parameter value (QVariant)
    // Outputs: None
    void setParameter(const QString &key, const QVariant &value);

    // Retrieves the value for a given parameter key
    // Inputs:
    // - key: The parameter key (QString)
    // Outputs:
    // - Returns the value associated with the key (QVariant)
    QVariant parameter(const QString &key) const;

signals:
    // Signal emitted on text change
    // Inputs:
    // - name: The object name of the line edit widget (QString)
    // - text: The new text entered in the line edit (QString)
    // Outputs: None
    void textChanged(const QString &name, const QString &text);

private slots:
    // Slot to handle text change
    // Inputs:
    // - text: The new text entered in the line edit (QString)
    // Outputs: None
    void onTextChanged(const QString &text);

private:
    QMap<QString, QVariant> parameters;   // Map to store parameters
};

#endif // MORIS_LINE_EDIT_HPP
