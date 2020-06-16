#-------------------------------------------------
#
# Project created by QtCreator 2020-01-23T10:43:01
#
#-------------------------------------------------

QT += core gui widgets
QT += printsupport          #Needed for making qcustomplot work

TARGET = PLIS
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11

win32:RC_ICONS = logo.ico
macx:ICON = logo.icns


SOURCES += \
    help.cpp \
        main.cpp \
    modify_text.cpp \
    optimize.cpp \
    plis.cpp \
    data.cpp \
    plot.cpp \
    qcustomplot.cpp \
    calculations.cpp \
    result.cpp \
    simulatedialog.cpp \
    save_open.cpp \
    jackknife.cpp \
    legend.cpp \
    modify_data.cpp \
    show_result.cpp \
    guess_parameters.cpp \
    table.cpp

HEADERS += \
    help.h \
    modify_text.h \
    optimize.h \
    plis.h \
    data.h \
    plot.h \
    qcustomplot.h \
    calculations.h \
    fitmrq.h \
    result.h \
    simulatedialog.h \
    nr3.h \
    save_open.h \
    jackknife.h \
    legend.h \
    modify_data.h \
    show_result.h \
    guess_parameters.h \
    table.h

FORMS += \
    help.ui \
    modify_text.ui \
    optimize.ui \
    plis.ui \
    simulatedialog.ui \
    modify_data.ui \
    show_result.ui \
    guess_parameters.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

RESOURCES += \
    Resources.qrc \

DISTFILES += \
    img/background.png \
    img/bgui.png \
    img/color.png \
    img/deleteGraph.png \
    img/exit_icn.png \
    img/expand.png \
    img/export.png \
    img/file.png \
    img/fit.png \
    img/ftest.png \
    img/help.png \
    img/logo.png \
    img/logo2(100).png \
    img/logo2(50).png \
    img/logo2(75).png \
    img/merge.png \
    img/modelbutton.png \
    img/open.png \
    img/pen.png \
    img/point.png \
    img/preset.png \
    img/printer.png \
    img/rc_icon.png \
    img/save.png \
    img/trbutton.png \
    img/rc_icon.icns \
    img/shapes.png \
    img/shapes.png
