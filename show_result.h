#ifndef SHOW_RESULT_H
#define SHOW_RESULT_H

#include <QDialog>
#include "data.h"

namespace Ui {
class Show_result;
}

class Show_result : public QDialog
{
    Q_OBJECT

public:
    explicit Show_result(QWidget *parent, Data* &d);
    ~Show_result();

private:
    Ui::Show_result *ui;
    Data* data{};

    void showResultOneSite();
    void showResultTwoSite();
    void showResultFourSite();
    void showResultComp();
    void showResultCPMG();
    void writeLine(QLineEdit *line, double res, double err, QString unit);
    void writeLine(QLineEdit *valbox, QLineEdit *errbox, double res, double err, QString unit);

};

#endif // SHOW_RESULT_H
