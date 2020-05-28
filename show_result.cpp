#include "show_result.h"
#include "ui_show_result.h"

Show_result::Show_result(QWidget *parent, Data* &d) :
    QDialog(parent),
    ui(new Ui::Show_result)
{
    ui->setupUi(this);
    data=d;
    if(data->num_bind_site==1){
        ui->tabWidget->removeTab(1);
        ui->tabWidget->removeTab(1);
        ui->tabWidget->removeTab(1);
        ui->chi2->setText(QString::number(data->chi2));
        ui->name->setText(data->name);
        ui->datapoints->setText(QString::number(data->responceVector.size()));
        ui->kd1->setText(QString::number(data->a[2],'f',2)  + " ± " + QString::number(data->kd_error,'f',2) + " " + data->ligand_unit);
        ui->start_res->setText(QString::number(data->a[0],'f',2));
        ui->end_res->setText(QString::number(data->a[1],'f',2));
        ui->model->setText(data->model);
    }
    else if(data->num_bind_site==2){
        ui->tabWidget->removeTab(0);
        ui->tabWidget->removeTab(1);
        ui->tabWidget->removeTab(1);
        ui->chi2_2->setText(QString::number(data->chi2));
        ui->name_2->setText(data->name);
        ui->datapoints_2->setText(QString::number(data->responceVector.size()));
        ui->kd1_2->setText(QString::number(data->a[3],'f',2) + " ± " + QString::number(data->kd1_error,'f',2) + " " + data->ligand_unit);
        ui->kd2_2->setText(QString::number(data->a[4],'f',2) + " ± " + QString::number(data->kd2_error,'f',2) + " " + data->ligand_unit);
        ui->start_res_2->setText(QString::number(data->a[0],'f',2));
        ui->inter_res->setText(QString::number(data->a[1],'f',2));
        ui->end_res_2->setText(QString::number(data->a[2],'f',2));
        ui->model_2->setText(data->model);
    }
    else if(data->num_bind_site==3){
        ui->tabWidget->removeTab(0);
        ui->tabWidget->removeTab(0);
        ui->tabWidget->removeTab(1);
        ui->chi2_3->setText(QString::number(data->chi2));
        ui->name_3->setText(data->name);
        ui->datapoints_3->setText(QString::number(data->responceVector.size()));
        ui->kd1_3->setText(QString::number(data->a[2],'f',2) + " ± " + QString::number(data->kd1_error,'f',2) + " " + data->ligand_unit);
        ui->kd2_3->setText(QString::number(data->a[3],'f',2) + " ± " + QString::number(data->kd2_error,'f',2) + " " + data->ligand_unit);
        ui->kdc_3->setText(QString::number(data->a[4],'f',2) + " ± " + QString::number(data->kdc_error,'f',2) + " " + data->ligand_unit);
        ui->start_res_3->setText(QString::number(data->a[0]));
        ui->end_res_3->setText(QString::number(data->a[1]));
        ui->model_3->setText(data->model);
    }
    else if(data->num_bind_site==4){
        ui->tabWidget->removeTab(0);
        ui->tabWidget->removeTab(0);
        ui->tabWidget->removeTab(0);
        ui->chi2_4->setText(QString::number(data->chi2));
        ui->name_4->setText(data->name);
        ui->model_4->setText(data->model);
        ui->datapoints_4->setText(QString::number(data->n_cpmgVector.size()));
        ui->kd1_4->setText(QString::number(data->kd_cpmg,'f',2)  + " ± " + QString::number(data->kd_error,'f',2) + " " + data->ligand_unit);
        ui->r20->setText(QString::number(data->a[0],'f',2));
        ui->kex->setText(QString::number(data->a[1],'f',2));
        ui->dw->setText(QString::number(data->a[3],'f',2));
        ui->pb->setText(QString::number(data->a[2],'f',2));
    }
}

Show_result::~Show_result()
{
    delete ui;
}
