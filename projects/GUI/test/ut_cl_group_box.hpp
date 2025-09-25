#pragma once

#include <QObject>
#include <QtTest>

namespace moris
{
    class Parameter;
    class Moris_Group_Box;
}

class ut_cl_group_box : public QObject
{
    Q_OBJECT

private slots:
    void initTestCase();
    void testInitialState();
    void testParameterChange();
    void testReadOnlyParameter();

};