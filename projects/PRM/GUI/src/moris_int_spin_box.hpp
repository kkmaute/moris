#ifndef MORIS_INT_SPIN_BOX_HPP
#define MORIS_INT_SPIN_BOX_HPP

#include <QSpinBox>
#include <QVariant>
#include <QMap>

class Moris_Int_Spin_Box : public QSpinBox
{
    Q_OBJECT

public:
    // Constructor
    explicit Moris_Int_Spin_Box(QWidget *parent = nullptr);

    // Destructor
    ~Moris_Int_Spin_Box();

    // Sets a parameter key-value pair in the parameters map
    void setParameter(const QString &key, const QVariant &value);

    // Retrieves the value for a given parameter key from the parameters map
    QVariant parameter(const QString &key) const;

signals:
    // Signal emitted when the value changes
    void valueChanged(const QString &name, const QVariant &new_value);

private slots:
    // Slot to handle value changes in the spin box
    void onValueChanged(int new_value);

private:
    // Store parameters as key-value pairs
    QMap<QString, QVariant> parameters;
};

#endif // MORIS_INT_SPIN_BOX_HPP
