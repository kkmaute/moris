#ifndef MORIS_COMBO_BOX_HPP
#define MORIS_COMBO_BOX_HPP

#include <QComboBox>
#include <QVariant>
#include <QMap>
#include <QDebug>

class Moris_Combo_Box : public QComboBox
{
    Q_OBJECT

  public:
    explicit Moris_Combo_Box( QWidget *parent = nullptr );    // Constructor with optional parent widget
    ~Moris_Combo_Box();                                       // Destructor

    // Sets a parameter key-value pair
    // Inputs:
    // - key: The parameter key (QString)
    // - value: The parameter value (QVariant)
    // Outputs: None
    void setParameter( const QString &key, const QVariant &value );

    // Retrieves the value for a given parameter key
    // Inputs:
    // - key: The parameter key (QString)
    // Outputs:
    // - Returns the value associated with the key (QVariant)
    QVariant parameter( const QString &key ) const;

  signals:
    // Signal emitted when the combo box index changes
    // Inputs:
    // - name: The name associated with the combo box
    // - index: The new index of the combo box
    void currentIndexChanged( const QString &name, int index );

  private slots:
    // Slot to handle index changes in the combo box
    // Inputs:
    // - index: The new selected index (int)
    // Outputs: None
    void onCurrentIndexChanged( int index );

  private:
    QMap< QString, QVariant > parameters;    // Map to store parameters
};

#endif    // MORIS_COMBO_BOX_HPP
