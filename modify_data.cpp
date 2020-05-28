#include "modify_data.h"
#include "ui_modify_data.h"
#include <QPixmap>

Modify_data::Modify_data(QWidget *parent, QVector<QCPScatterStyle::ScatterShape> s, Data *&d) :
    QDialog(parent),
    ui(new Ui::Modify_data), shapes{s}, data{d}
{
    ui->setupUi(this);
#ifdef Q_OS_MAC
    ui->label_4->setFont(QFont("Helvetica",16));
#endif
    ui->size_curve->setValue(data->pen_curve.width());
    ui->size_data->setValue(int(data->style.size()));
    int m=data->style.shape();
    if (m==5)
        ui->shape_box->setCurrentText("ssDisc");
    else if(m==6)
        ui->shape_box->setCurrentText("ssSquare");
    else if(m==7)
        ui->shape_box->setCurrentText("ssDiamond");
    else if(m==8)
        ui->shape_box->setCurrentText("ssStar");
    else if(m==2)
        ui->shape_box->setCurrentText("ssCross");
    else if(m==3)
        ui->shape_box->setCurrentText("ssPlus");
    else if(m==4)
        ui->shape_box->setCurrentText("ssCircle");
    else if(m==9)
        ui->shape_box->setCurrentText("ssTriangle");
    else if(m==10)
        ui->shape_box->setCurrentText("ssTriangleInverted");
    else if(m==11)
        ui->shape_box->setCurrentText("ssCrossSquare");
    else if(m==12)
        ui->shape_box->setCurrentText("ssPlusSquare");
    else if(m==13)
        ui->shape_box->setCurrentText("ssCrossCircle");
    else if(m==14)
        ui->shape_box->setCurrentText("ssPlusCircle");
    else if(m==15)
        ui->shape_box->setCurrentText("ssPeace");

    QColor temp=data->pen.color();
    QString string=temp.name();
    if (string=="#000000")
        ui->color_box->setCurrentText("Black");
    else if(string=="#ff0000")
        ui->color_box->setCurrentText("Red");
    else if(string=="#0000ff")
        ui->color_box->setCurrentText("Blue");
    else if(string=="#008000")
        ui->color_box->setCurrentText("Green");
    else if(string=="#00008b")
        ui->color_box->setCurrentText("Darkblue");
    else if(string=="#ff00ff")
        ui->color_box->setCurrentText("Magenta");
    else if(string=="#808080")
        ui->color_box->setCurrentText("Grey");
    else if(string=="#00ffff")
        ui->color_box->setCurrentText("Cyan");
    else if(string=="#8b008b")
        ui->color_box->setCurrentText("Darkmagenta");

    temp=data->pen_curve.color();
    string=temp.name();
    if (string=="#000000")
        ui->color_curve_box->setCurrentText("Black");
    else if(string=="#ff0000")
        ui->color_curve_box->setCurrentText("Red");
    else if(string=="#0000ff")
        ui->color_curve_box->setCurrentText("Blue");
    else if(string=="#008000")
        ui->color_curve_box->setCurrentText("Green");
    else if(string=="#00008b")
        ui->color_curve_box->setCurrentText("Darkblue");
    else if(string=="#ff00ff")
        ui->color_curve_box->setCurrentText("Magenta");
    else if(string=="#808080")
        ui->color_curve_box->setCurrentText("Grey");
    else if(string=="#00ffff")
        ui->color_curve_box->setCurrentText("Cyan");
    else if(string=="#8b008b")
        ui->color_curve_box->setCurrentText("Darkmagenta");
}

Modify_data::~Modify_data()
{
    delete ui;
}

void Modify_data::on_okButton_clicked()
{
    data->style.setShape(shapes[ui->shape_box->currentIndex()]);
    data->pen.setColor(ui->color_box->currentText());
    data->pen_curve.setColor(ui->color_curve_box->currentText());
    data->pen_curve.setWidth(ui->size_curve->value());
    data->style.setSize(ui->size_data->value());
    hide();
}

void Modify_data::on_cancelButton_clicked()
{
    hide();
}
