
#include "ut_gui_startup.hpp"

#include <QtTest> // may need to update this to QtTest in future, if extra libraries are necessary
#include <QProcess> // to launch and kill gui process (idk if theres some alternative)
#include <QCoreApplication>
#include <QDir>


void ut_gui_startup::testGuiLaunch()
{
    QProcess tGUILaunchProcess; 
    QString tBaseDirectory = QCoreApplication::applicationDirPath();

    QDir tDirectory(tBaseDirectory); 
    tDirectory.cdUp();
    tDirectory.cdUp();
    tDirectory.cd("src");
    QString tPathGUIexe = tDirectory.absoluteFilePath("moris_gui");
    QFileInfo info(tPathGUIexe); // tInfoGUIexe
    QVERIFY2(info.exists() && info.isExecutable(),
            qPrintable("Cannot find or execute: " + tPathGUIexe));

    tGUILaunchProcess.start(tPathGUIexe);

        // Ensure the process started successfully
    QVERIFY(tGUILaunchProcess.waitForStarted(3000));
    QVERIFY(tGUILaunchProcess.state() == QProcess::Running);

    tGUILaunchProcess.kill();
    tGUILaunchProcess.waitForFinished();
    QTest::qWait(2000); // wait for the GUI to close
    
      // Check to make sure process is no longer running
    QVERIFY(tGUILaunchProcess.state() != QProcess::Running);
}
