#ifndef INPUT_PARAMETERS_HPP
#define INPUT_PARAMETERS_HPP

#include <QString>
#include <QVariant>
#include <QMap>

class Input_Parameters
{
  public:
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

    // Retrieves all parameters as a map
    // Outputs:
    // - Returns a map of all parameters (QMap<QString, QVariant>)
    QMap< QString, QVariant > allParameters() const;

    // Comparison operator for equality
    // Inputs:
    // - other: The other Input_Parameters object to compare with
    // Outputs:
    // - Returns true if both objects have the same parameters, false otherwise (bool)
    bool operator==( const Input_Parameters &other ) const;

    // Comparison operator for inequality
    // Inputs:
    // - other: The other Input_Parameters object to compare with
    // Outputs:
    // - Returns true if objects have different parameters, false if they are equal (bool)
    bool operator!=( const Input_Parameters &other ) const;

  private:
    QMap< QString, QVariant > parameters;    // Map to store parameters
};

#endif    // INPUT_PARAMETERS_HPP
