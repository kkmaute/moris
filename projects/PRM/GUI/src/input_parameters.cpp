#include "input_parameters.hpp"

// Sets a parameter key-value pair in the parameters map
// Inputs:
// - key: The parameter key (QString)
// - value: The parameter value (QVariant)
// Outputs: None
void Input_Parameters::setParameter(const QString &key, const QVariant &value)
{
    // Insert the key-value pair into the parameters map
    parameters.insert(key, value);
}

// Retrieves the value for a given parameter key
// Inputs:
// - key: The parameter key (QString)
// Outputs:
// - Returns the value associated with the key (QVariant)
QVariant Input_Parameters::parameter(const QString &key) const
{
    // Return the value associated with the provided key from the parameters map
    return parameters.value(key);
}

// Retrieves all parameters as a map
// Outputs:
// - Returns a map of all parameters (QMap<QString, QVariant>)
QMap<QString, QVariant> Input_Parameters::allParameters() const
{
    // Return the entire parameters map
    return parameters;
}

// Comparison operator for equality
// Inputs:
// - other: The other Input_Parameters object to compare with
// Outputs:
// - Returns true if both objects have the same parameters, false otherwise (bool)
bool Input_Parameters::operator==(const Input_Parameters &other) const
{
    // Compare the parameters map of the current object with the other object
    return parameters == other.parameters;
}

// Comparison operator for inequality
// Inputs:
// - other: The other Input_Parameters object to compare with
// Outputs:
// - Returns true if objects have different parameters, false if they are equal (bool)
bool Input_Parameters::operator!=(const Input_Parameters &other) const
{
    // Return the negation of the equality comparison
    return !(*this == other);
}
