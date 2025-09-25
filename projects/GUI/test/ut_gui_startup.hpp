#pragma once

#include <QObject>
#include <QtTest>

class ut_gui_startup : public QObject
{
    Q_OBJECT

private slots:

    void testGuiLaunch();
};