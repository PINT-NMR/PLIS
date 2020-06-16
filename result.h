#ifndef RESULT_H
#define RESULT_H

#include "plis.h"


class Result
{
public:
    Result();
    Result(QLineEdit *kd1_box, QLineEdit *kd1_error, QLabel *kd1_label, QLabel *plus_minus1,
           QLineEdit *kd2_box, QLineEdit *kd2_error, QLabel *Kd2_label, QLabel *plus_minus2,
           QLineEdit *kd3_box, QLineEdit *kd3_error, QLabel *Kd3_label, QLabel *plus_minus3,
           QLineEdit *kd4_box, QLineEdit *kd4_error, QLabel *Kd4_label, QLabel *plus_minus4,
           QLineEdit *name, QLineEdit *model, QLineEdit *num_data, QLineEdit *chi2_box);
    void hideKd2();
    void hideKd3();
    void hideKd4();
    void showKd2();
    void showKd3();
    void showKd4();
    void clearAll();
    void clearNameModelNdataChi2();
    void clearKd1();
    void clearKd2();
    void writeToKdBox(QLineEdit *kdbox, double d);
    void writeToKdErrorBox(QLineEdit *kdbox, double d, QString unit);
    void writeNameModelNdataChi2(QString mode, Data *dataSet);
    void writeKd1(Data *dataSet);
    void writeKd1_CPMG(Data *dataSet);
    void writeKd2(Data *dataSet);
    void writeKd3(Data *dataSet);
    void writeKd4(Data *dataSet);

private:
    QLineEdit *_kd1_box, *_kd1_error;
    QLabel *_kd1_label, *_plus_minus1;
    QLineEdit *_kd2_box, *_kd2_error;
    QLabel *_kd2_label, *_plus_minus2;
    QLineEdit *_kd3_box, *_kd3_error;
    QLabel *_Kd3_label, *_plus_minus3;
    QLineEdit *_kd4_box, *_kd4_error;
    QLabel *_Kd4_label, *_plus_minus4;
    QLineEdit *_name, *_model,  *_num_data, *_chi2_box;
};

#endif // RESULT_H
