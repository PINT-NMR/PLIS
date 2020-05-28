#include "modify_text.h"
#include "ui_modify_text.h"

Modify_text::Modify_text(QWidget *parent, QCPAxis *axis, QFont &x, QFont &y, QString &xlab, QString &ylab) :
    QDialog(parent), //axis_local(axis),
    ui(new Ui::Modify_text)
{
    ui->setupUi(this);
    mode=0;
    axis_local=axis;
    xlabel_font=&x;
    ylabel_font=&y;
    xlabel=&xlab;
    ylabel=&ylab;
#ifdef Q_OS_MAC
    ui->label->setFont(QFont("Arial",14));
#endif
}

Modify_text::Modify_text(QWidget *parent, QCPTextElement* &title):
    QDialog (parent), ui(new Ui::Modify_text)
{
    ui->setupUi(this);
    title_local=title;
    mode=1;
#ifdef Q_OS_MAC
    ui->label->setFont(QFont("Arial",14));
#endif
}
Modify_text::~Modify_text()
{
    delete ui;
}

void Modify_text::on_font_clicked(){
    hide();
    bool ok{};
    if (mode==0){
        QFont font = QFontDialog::getFont(&ok, axis_local->labelFont(), this, "Plis - Choose Font");
        if (ok){
            if (axis_local->axisType()==QCPAxis::atBottom)
                *xlabel_font=font;
            else if(axis_local->axisType()==QCPAxis::atLeft)
                *ylabel_font=font;
        }
    }
    else if(mode==1){
        QFont font = QFontDialog::getFont(&ok, title_local->font() , this, "Plis - Choose Font");
        if (ok){
            title_local->setFont(font);
        }
    }
}

void Modify_text::on_label_text_clicked()
{
    hide();
    bool ok{};
    if(mode==0){
        QString label = QInputDialog::getText(this, "Plis", "New axis label:", QLineEdit::Normal, axis_local->label(), &ok);
        if (ok && label!="")
        {
            if (axis_local->axisType()==QCPAxis::atBottom)
                *xlabel=label;
            else if(axis_local->axisType()==QCPAxis::atLeft)
                *ylabel=label;
        }
    }
    else if(mode==1){
        QString newTitle = QInputDialog::getText(this, "Plis", "New plot title:", QLineEdit::Normal, title_local->text(), &ok);
        if (ok && newTitle!="")
            title_local->setText(newTitle);
    }
}

void Modify_text::on_both_clicked()
{
    on_font_clicked();
    on_label_text_clicked();
}
