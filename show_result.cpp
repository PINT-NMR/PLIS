#include "show_result.h"
#include "ui_show_result.h"

Show_result::Show_result(QWidget *parent, Data* &d) :
    QDialog(parent), ui(new Ui::Show_result), data(d)
{
    ui->setupUi(this);
    if (data->num_bind_site==1)
        showResultOneSite();
    else if (data->num_bind_site==2)
        showResultTwoSite();
    else if (data->num_bind_site==5)
        showResultFourSite();
    else if (data->num_bind_site==3)
        showResultComp();
    else if (data->num_bind_site==4)
        showResultCPMG();
}

void Show_result::writeLine(QLineEdit *line, double res, double err, QString unit="")
{
    QString result{}, error{};
    if (fabs(res)<0.001)
        result = QString::number(res,'e',3);
    else if (fabs(res)<0.01)
        result = QString::number(res,'f',6);
    else if (fabs(res)<0.1)
        result = QString::number(res,'f',5);
    else if (fabs(res)<1)
        result = QString::number(res,'f',4);
    else if (fabs(res)<10)
        result = QString::number(res,'f',3);
    else if (fabs(res)<100)
        result = QString::number(res,'f',2);
    else if (fabs(res)<1000)
        result = QString::number(res,'f',1);
    else if (fabs(res)<10000)
        result = QString::number(res,'f',0);
    else
        result = QString::number(res,'e',3);

    if (err<0.01)
        error = QString::number(err,'e',2);
    else if (err<0.1)
        error = QString::number(err,'f',4);
    else if (err<1)
        error = QString::number(err,'f',3);
    else if (err<10)
        error = QString::number(err,'f',2);
    else if (err<100)
        error = QString::number(err,'f',1);
    else if (err<1000)
        error = QString::number(err,'f',0);
    else
        error = QString::number(err,'e',2);

    if (unit != "")
        line->setText(result  + " ± " + error + " " + unit);
    else
        line->setText(result  + " ± " + error);
}

void Show_result::writeLine(QLineEdit *valbox, QLineEdit *errbox, double res, double err, QString unit="")
{
    QString result{}, error{};
    if (fabs(res)<0.001)
        result = QString::number(res,'e',3);
    else if (fabs(res)<0.01)
        result = QString::number(res,'f',6);
    else if (fabs(res)<0.1)
        result = QString::number(res,'f',5);
    else if (fabs(res)<1)
        result = QString::number(res,'f',4);
    else if (fabs(res)<10)
        result = QString::number(res,'f',3);
    else if (fabs(res)<100)
        result = QString::number(res,'f',2);
    else if (fabs(res)<1000)
        result = QString::number(res,'f',1);
    else if (fabs(res)<10000)
        result = QString::number(res,'f',0);
    else
        result = QString::number(res,'e',3);

    if (err<0.01)
        error = QString::number(err,'e',2);
    else if (err<0.1)
        error = QString::number(err,'f',4);
    else if (err<1)
        error = QString::number(err,'f',3);
    else if (err<10)
        error = QString::number(err,'f',2);
    else if (err<100)
        error = QString::number(err,'f',1);
    else if (err<1000)
        error = QString::number(err,'f',0);
    else
        error = QString::number(err,'e',2);

     if (unit != "") {
         valbox->setText(result);
         errbox->setText(error);
     }
     else {
         valbox->setText(result);
         errbox->setText(error);
     }
}


Show_result::~Show_result()
{
    delete ui;
}

void Show_result::showResultOneSite() {
    ui->tabWidget->removeTab(1);
    ui->tabWidget->removeTab(1);
    ui->tabWidget->removeTab(1);
    ui->tabWidget->removeTab(1);
    ui->chi2->setText(QString::number(data->chi2));
    ui->name->setText(data->name);
    ui->datapoints->setText(QString::number(data->responceVector.size()));
    writeLine(ui->kd1, ui->kd_error_1, data->a[2], data->kd_error, data->ligand_unit);
    writeLine(ui->start_res, ui->P_res_error_1, data->a[0], data->dyda[0]);
    writeLine(ui->end_res, ui->PL_res_error_1, data->a[1], data->dyda[1]);

    ui->model->setText(data->model);
}

void Show_result::showResultTwoSite() {
    ui->tabWidget->removeTab(0);
    ui->tabWidget->removeTab(1);
    ui->tabWidget->removeTab(1);
    ui->tabWidget->removeTab(1);
    ui->chi2_2->setText(QString::number(data->chi2));
    ui->name_2->setText(data->name);
    ui->datapoints_2->setText(QString::number(data->responceVector.size()));
    writeLine(ui->kd1_2, ui->kd1_error_2, data->a[3], data->kd1_error, data->ligand_unit);
    writeLine(ui->kd2_2, ui->kd2_error_2, data->a[4], data->kd2_error, data->ligand_unit);
    writeLine(ui->start_res_2, ui->P_res_error_2,data->a[0], data->dyda[0]);
    writeLine(ui->inter_res, ui->PL_res_error_2,data->a[1], data->dyda[1]);
    writeLine(ui->end_res_2, ui->PL2_res_error_2,data->a[2], data->dyda[2]);
    ui->model_2->setText(data->model);
}

void Show_result::showResultFourSite() {
    ui->tabWidget->removeTab(0);
    ui->tabWidget->removeTab(0);
    ui->tabWidget->removeTab(1);
    ui->tabWidget->removeTab(1);
    ui->scrollArea_5->setStyleSheet("background-color:transparent white;"); // there is a linear bg-gradient w/o this line IDK why.
    ui->chi2_four->setText(QString::number(data->chi2));
    ui->name_four->setText(data->name);
    ui->datapoints_four->setText(QString::number(data->responceVector.size()));
    ui->model_four->setText(data->model);
    writeLine(ui->kd1_four, ui->kd1_error_four,data->a[5], data->kd1_error, data->ligand_unit);
    writeLine(ui->kd2_four, ui->kd2_error_four,data->a[6], data->kd2_error, data->ligand_unit);
    writeLine(ui->kd3_four, ui->kd3_error_four,data->a[7], data->kd3_error, data->ligand_unit);
    writeLine(ui->kd4_four, ui->kd4_error_four,data->a[8], data->kd4_error, data->ligand_unit);
    writeLine(ui->P_res_4, ui->P_res_error_four, data->a[0], data->dyda[0]);
    writeLine(ui->PL_res_4, ui->PL_res_error_four, data->a[1], data->dyda[1]);
    writeLine(ui->PL2_res_4, ui->PL2_res_error_four, data->a[2], data->dyda[2]);
    writeLine(ui->PL3_res_4, ui->PL3_res_error_four, data->a[3], data->dyda[3]);
    writeLine(ui->PL4_res_4, ui->PL4_res_error_four, data->a[4], data->dyda[4]);
}

void Show_result::showResultComp() {
    ui->tabWidget->removeTab(0);
    ui->tabWidget->removeTab(0);
    ui->tabWidget->removeTab(0);
    ui->tabWidget->removeTab(1);
    ui->chi2_3->setText(QString::number(data->chi2));
    ui->name_3->setText(data->name);
    ui->datapoints_3->setText(QString::number(data->responceVector.size()));
    writeLine(ui->kd1_3, ui->kd1_error_3, data->a[2], data->kd1_error, data->ligand_unit);
    writeLine(ui->kd2_3, ui->kd2_error_3, data->a[3], data->kd2_error, data->ligand_unit);
    writeLine(ui->kdc_3, ui->kdc_error_3, data->a[4], data->kdc_error, data->ligand_unit);
    writeLine(ui->start_res_3, ui->C_res_error_3, data->a[0], data->dyda[0]);
    writeLine(ui->end_res_3, ui->CL_res_error_3, data->a[1], data->dyda[1]);
    ui->model_3->setText(data->model);
}

void Show_result::showResultCPMG() {
    ui->tabWidget->removeTab(0);
    ui->tabWidget->removeTab(0);
    ui->tabWidget->removeTab(0);
    ui->tabWidget->removeTab(0);
    ui->chi2_4->setText(QString::number(data->chi2));
    ui->name_4->setText(data->name);
    ui->model_4->setText(data->model);
    ui->datapoints_4->setText(QString::number(data->n_cpmgVector.size()));
    ui->kd1_4->setText(QString::number(data->kd_cpmg,'f',2));
    ui->r20->setText(QString::number(data->a[0],'f',2));
    ui->kex->setText(QString::number(data->a[1],'f',2));
    ui->pb->setText(QString::number(data->a[2],'f',2));
    ui->dw->setText(QString::number(data->a[3],'f',2));

    ui->kd_error_cpmg->setText(QString::number(data->kd_error,'f',2));
    ui->r20_error_cpmg->setText(QString::number(data->dyda[0],'f',2));
    ui->kex_error_cpmg->setText(QString::number(data->dyda[1],'f',2));
    ui->pb_error_cpmg->setText(QString::number(data->dyda[2],'f',2));
    ui->dw_error_cpmg->setText(QString::number(data->dyda[3],'f',2));
}
