#ifndef HELP_H
#define HELP_H

#include <QDialog>

namespace Ui {
class Help;
}

class Help : public QDialog
{
    Q_OBJECT

public:
    explicit Help(QWidget *parent = nullptr);
    ~Help();

    void setup(int i=0);

private slots:
    void on_listWidget_currentTextChanged(const QString &currentText);

    void on_pushButton_2_clicked();

private:
    Ui::Help *ui;
};

#endif // HELP_H
