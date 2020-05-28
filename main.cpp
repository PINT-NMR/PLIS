#include <QApplication>
#include "plis.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Plis* p=new Plis{};
    QPixmap pixmap(":/img/img/logo2(100).png");
    QSplashScreen *splash=new QSplashScreen;
    splash->setPixmap(pixmap);
    splash->show();
    QTimer::singleShot(3000,splash,SLOT(close()));
    QTimer::singleShot(3000,p,SLOT(showMaximized()));
    //p->showMaximized();

    return a.exec();
}
