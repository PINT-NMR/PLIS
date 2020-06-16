#include "result.h"

Result::Result()
{
}

Result::Result(QLineEdit *kd1_box, QLineEdit *kd1_error, QLabel *kd1_label, QLabel *plus_minus1,
               QLineEdit *kd2_box, QLineEdit *kd2_error, QLabel *kd2_label, QLabel *plus_minus2,
               QLineEdit *kd3_box, QLineEdit *kd3_error, QLabel *Kd3_label, QLabel *plus_minus3,
               QLineEdit *kd4_box, QLineEdit *kd4_error, QLabel *Kd4_label, QLabel *plus_minus4,
               QLineEdit *name, QLineEdit *model, QLineEdit *num_data, QLineEdit *chi2_box) :
    _kd1_box(kd1_box), _kd1_error(kd1_error), _kd1_label(kd1_label), _plus_minus1(plus_minus1),
    _kd2_box(kd2_box), _kd2_error(kd2_error), _kd2_label(kd2_label), _plus_minus2(plus_minus2),
    _kd3_box(kd3_box), _kd3_error(kd3_error), _Kd3_label(Kd3_label), _plus_minus3(plus_minus3),
    _kd4_box(kd4_box), _kd4_error(kd4_error), _Kd4_label(Kd4_label), _plus_minus4(plus_minus4),
    _name(name),  _model(model), _num_data(num_data), _chi2_box(chi2_box)
{
    _kd1_label->setText("Kd");
}

void Result::hideKd2()
{
    _kd2_box->setVisible(false);
    _kd2_error->setVisible(false);
    _kd2_label->setVisible(false);
    _plus_minus2->setVisible(false);
}

void Result::hideKd3()
{
    _kd3_box->setVisible(false);
    _kd3_error->setVisible(false);
    _Kd3_label->setVisible(false);
    _plus_minus3->setVisible(false);
}

void Result::hideKd4()
{
    _kd4_box->setVisible(false);
    _kd4_error->setVisible(false);
    _Kd4_label->setVisible(false);
    _plus_minus4->setVisible(false);
}

void Result::showKd2()
{
    _kd2_box->setVisible(true);
    _kd2_error->setVisible(true);
    _kd2_label->setVisible(true);
    _plus_minus2->setVisible(true);
}

void Result::showKd3()
{
    _kd3_box->setVisible(true);
    _kd3_error->setVisible(true);
    _Kd3_label->setVisible(true);
    _plus_minus3->setVisible(true);
}

void Result::showKd4()
{
    _kd4_box->setVisible(true);
    _kd4_error->setVisible(true);
    _Kd4_label->setVisible(true);
    _plus_minus4->setVisible(true);
}

void Result::clearNameModelNdataChi2()
{
    _name->setText("");
    _model->setText("");
    _num_data->setText("");
    _chi2_box->setText("");
}

void Result::clearKd1()
{
    _kd1_error->setText("");
    _kd1_box->setText("");
}

void Result::clearKd2()
{
    _kd2_error->setText("");
    _kd2_box->setText("");
}

void Result::writeToKdBox(QLineEdit *kdbox, double d)
{
    if (d<0.001)
        kdbox->setText(QString::number(d,'e', 3));
    if (d<0.01)
        kdbox->setText(QString::number(d,'e', 6));
    else if (d<0.1)
        kdbox->setText(QString::number(d,'f', 5));
    else if (d<1)
        kdbox->setText(QString::number(d,'f', 4));
    else if (d<10)
        kdbox->setText(QString::number(d,'f', 3));
    else if (d<100)
        kdbox->setText(QString::number(d,'f', 2));
    else if (d<1000)
        kdbox->setText(QString::number(d,'f', 1));
    else if (d<10000)
        kdbox->setText(QString::number(d,'f', 0));
    else
        kdbox->setText(QString::number(d,'e', 3));
}

void Result::writeToKdErrorBox(QLineEdit *kdbox, double d, QString unit)
{
    if (d<0.01)
        kdbox->setText(QString::number(d,'e', 2) + " " + unit);
    else if (d<0.1)
        kdbox->setText(QString::number(d,'f', 4) + " " + unit);
    else if (d<1)
        kdbox->setText(QString::number(d,'f', 3) + " " + unit);
    else if (d<10)
        kdbox->setText(QString::number(d,'f', 2) + " " + unit);
    else if (d<100)
        kdbox->setText(QString::number(d,'f', 1) + " " + unit);
    else if (d<1000)
        kdbox->setText(QString::number(d,'f', 0) + " " + unit);
    else
        kdbox->setText(QString::number(d,'e', 2) + " " + unit);
}

void Result::writeNameModelNdataChi2(QString mode, Data *dataSet)
{
    hideKd2();
    hideKd3();
    hideKd4();
    _name->setText(dataSet->name);
    _model->setText(dataSet->model);
    if(mode=="standard" || mode=="chemshift")
        _num_data->setText(QString::number(dataSet->responceVector.size()));
    else
        _num_data->setText(QString::number(dataSet->n_cpmgVector.size()));
    _chi2_box->setText(QString::number(dataSet->chi2));
}

void Result::writeKd1(Data *dataSet)
{
    _kd1_label->setText("Kd");
    writeToKdBox(_kd1_box, dataSet->a[2]);
    writeToKdErrorBox(_kd1_error, dataSet->kd_error, dataSet->ligand_unit);
}

void Result::writeKd1_CPMG(Data *dataSet)
{
    _kd1_label->setText("Kd");
    writeToKdBox(_kd1_box, dataSet->kd_cpmg);
    writeToKdErrorBox(_kd1_error, dataSet->kd_error, dataSet->ligand_unit);
}

void Result::writeKd2(Data *dataSet)
{
    _kd1_label->setText("Kd1");
    showKd2();
    writeToKdBox(_kd1_box, dataSet->a[3]);
    writeToKdBox(_kd2_box, dataSet->a[4]);
    writeToKdErrorBox(_kd1_error, dataSet->kd1_error, dataSet->ligand_unit);
    writeToKdErrorBox(_kd2_error, dataSet->kd2_error, dataSet->ligand_unit);
}

void Result::writeKd3(Data *dataSet)
{
    _kd1_label->setText("Kd1");
    _Kd3_label->setText("Kd(C)");
    showKd2();
    showKd3();
    writeToKdBox(_kd1_box, dataSet->a[2]);
    writeToKdBox(_kd2_box, dataSet->a[3]);
    writeToKdBox(_kd3_box, dataSet->a[4]);
    writeToKdErrorBox(_kd1_error, dataSet->kd1_error, dataSet->ligand_unit);
    writeToKdErrorBox(_kd2_error, dataSet->kd2_error, dataSet->ligand_unit);
    writeToKdErrorBox(_kd3_error, dataSet->kdc_error, dataSet->ligand_unit);
}

void Result::writeKd4(Data *dataSet)
{
    _kd1_label->setText("Kd1");
    _Kd3_label->setText("Kd3");
    showKd2();
    showKd3();
    showKd4();
    writeToKdBox(_kd1_box, dataSet->a[5]);
    writeToKdBox(_kd2_box, dataSet->a[6]);
    writeToKdBox(_kd3_box, dataSet->a[7]);
    writeToKdBox(_kd4_box, dataSet->a[8]);
    writeToKdErrorBox(_kd1_error, dataSet->kd1_error, dataSet->ligand_unit);
    writeToKdErrorBox(_kd2_error, dataSet->kd2_error, dataSet->ligand_unit);
    writeToKdErrorBox(_kd3_error, dataSet->kd3_error, dataSet->ligand_unit);
    writeToKdErrorBox(_kd4_error, dataSet->kd4_error, dataSet->ligand_unit);
}

void Result::clearAll() {
    clearNameModelNdataChi2();
    clearKd1();
    clearKd2();
}
