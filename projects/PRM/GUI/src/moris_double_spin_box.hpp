#ifndef MORIS_DOUBLE_SPIN_BOX_HPP
#define MORIS_DOUBLE_SPIN_BOX_HPP

#include <QDoubleSpinBox>
#include <QVariant>
#include <QMap>

class Moris_Double_Spin_Box : public QDoubleSpinBox
{
    Q_OBJECT

public:
    // Constructor
    // Inputs:
    // - parent: The parent widget (QWidget *)
    explicit Moris_Double_Spin_Box(QWidget *parent = nullptr);

    // Destructor
    ~Moris_Double_Spin_Box();

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
    // Signal emitted when the value changes
    // Outputs:
    // - name: The name of the spin box object (QString)
    // - new_value: The new value of the spin box (QVariant)
    void valueChanged(const QString &name, const QVariant &new_value);

private slots:
    // Slot to handle value changes in the spin box
    // Inputs:
    // - new_value: The new value of the spin box (double)
    void onValueChanged(double new_value);

private:
    QMap<QString, QVariant> parameters; // Store parameters as key-value pairs
};

#endif // MORIS_DOUBLE_SPIN_BOX_HPP
