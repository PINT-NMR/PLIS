#include "save_open.h"

void save_project(QMainWindow *this1, const QString &x, const QString &y,
                  const QVector<Data*> &dataSets, const int &current_dataSet, const bool &is_disconnected,
                  int &shape, QCPTextElement *&title, QFont &xaxis_font, QFont &yaxis_font, QString &mode)
{
    std::string saveFilename = QFileDialog::getSaveFileName(this1,"Save as", "Choose a filename", "PROJ(*.proj)").toStdString();
    if(saveFilename=="")
        return;
    std::ofstream file;
    file.open(saveFilename.c_str(), std::ios::out);
    file << "x-axis: " << x.toStdString() << "\n";
    file << "x-axis_font: " << xaxis_font.family().toStdString() << " " << xaxis_font.pointSize() << "\n";
    file << "y-axis: " << y.toStdString() << "\n";
    file << "y-axis_font: " << yaxis_font.family().toStdString() << " " << yaxis_font.pointSize() << "\n";
    file << "DataSets: " << dataSets.size()<< "\n";
    file << "current_dataset: " << current_dataSet<< "\n";
    file << "is_disconnected: " << is_disconnected << "\n";
    file << "mode: " << mode.toStdString() << "\n";
    file << "shape: " << shape << "\n";
    file << "title: " << title->text().toStdString() << "\n";
    QString s=title->font().toString();
    s.replace(',',' ');
    file << "title_font: " << s.toStdString() << "\n";
    for (int i=0; i<dataSets.size(); i++){
        file << "dataSet " << i << "\n";
        file << "name: " << dataSets[i]->name.toStdString() << "\n";
        if(mode=="standard"){
            for  (int j=0; j<dataSets[i]->responceVector.size();j++){
                file << dataSets[i]->responceVector[j] << " " << dataSets[i]->concVector[j] << " "
                     << dataSets[i]->volumeVector[j] << " " << dataSets[i]->protein_conc_vector[j] << " " << dataSets[i]->error_vector[j]<< "\n";
            }
        }
        else if(mode=="cpmg"){
            for  (int j=0; j<dataSets[i]->n_cpmgVector.size();j++){
                file << dataSets[i]->n_cpmgVector[j] << " " << dataSets[i]->R2effVector[j] << " "
                     << dataSets[i]->dyVector[j] << "\n";
            }
        }
        file << "a: ";
        for (unsigned int j=0;j<dataSets[i]->a.size();j++)
            file << dataSets[i]->a[j] << " ";
        file << "\n";
        file << "dyda: ";
        for (unsigned int j=0;j<dataSets[i]->dyda.size();j++)
            file << dataSets[i]->dyda[j] << " ";
        file << " \n";
        file << "protein_conc: " << dataSets[i]->protein_conc << "\n";
        if(dataSets[i]->comp_vector.size()>0)
            file << "comp_conc: " << dataSets[i]->comp_vector[0] << "\n";
        file << "kd_cpmg: " << dataSets[i]->kd_cpmg << "\n";
        file << "ligand_cpmg: " << dataSets[i]->ligand_cpmg << "\n";
        file << "ligand_unit: " << dataSets[i]->ligand_unit.toStdString() << "\n";
        file << "unit: " << dataSets[i]->unit.toStdString() << "\n";
        file << "data_visible: " <<  dataSets[i]->data_visible << "\n";
        file << "curve_visible: " << dataSets[i]->curve_visible << "\n";
        file << "has_curve: " << dataSets[i]->has_curve << "\n";
        file << "num_bind_site: " << dataSets[i]->num_bind_site << "\n";
        file << "kd_error: " << dataSets[i]->kd_error << "\n";
        file << "kd1_error: " << dataSets[i]->kd1_error << "\n";
        file << "kd2_error: " << dataSets[i]->kd2_error << "\n";
        file << "kdc_error: " << dataSets[i]->kdc_error << "\n";
        file << "model: " << dataSets[i]->model.toStdString() << "\n";
        file << "chi2: " << dataSets[i]->chi2 << "\n";
        file << "pen: " << dataSets[i]->pen.color().name().toStdString() << "\n";
        file << "pen_curve: " << dataSets[i]->pen_curve.color().name().toStdString() << "\n";
        file << "style: " << dataSets[i]->style.shape() << " " << dataSets[i]->style.size() << "\n";
        if (dataSets[i]->has_curve==true){
            file << "calc_data: " << "\n";
            for (int k=0; k<dataSets[i]->calc_data[0].size();k++){
                file << dataSets[i]->calc_data[0][k] << " " << dataSets[i]->calc_data[1][k] << "\n";
            }
        }
    }
    file.close();

}

int open_project(QMainWindow *this1, QString &x, QString &y, QVector<Data *> &dataSets, int &current_dataSet,
                  bool &is_disconnected, QVector<QCPGraph*> &graphSets, int &shape, QCPTextElement *&title, QFont &xaxis_font,
                  QFont &yaxis_font, QString &mode)
{
    std::string openFilename = QFileDialog::getOpenFileName(this1, "Open", "A file", "PROJ(*.proj)").toStdString();
    if (openFilename=="")
        return 1;
    std::ifstream file;
    file.open(openFilename.c_str(), std::ios::in);
    dataSets.clear();
    graphSets.clear();
    x="";
    y="";
    int current_dataset_local{0};
    if(file.is_open())
    {
        std::string line{};
        std::vector<std::string> word{};
        bool stop=false;
        while (stop==false){
            getline(file >> std::ws, line);
            if (file.eof())
                stop=true;
            line2words(line, word);
            if (word[0]=="x-axis:")
                for (unsigned int i=1; i<word.size();i++){
                    QString temp;
                    temp=temp.fromStdString(word[i]);
                    x+=temp + " ";
                }
            else if (word[0]=="x-axis_font:"){
                xaxis_font=QFont(QString::fromStdString(word[1]),toInt(word[2]));
            }
            else if (word[0]=="y-axis:"){
                for (unsigned int i=1; i<word.size();i++){
                    QString temp;
                    temp=temp.fromStdString(word[i]);
                    y+=temp+" ";
                }
            }

            else if(word[0]=="y-axis_font:")
                yaxis_font=QFont(QString::fromStdString(word[1]),toInt(word[2]));

            else if(word[0]=="DataSets:"){
                double temp;
                toDouble(word[1],temp);
                for (int i=0; i<temp ;i++){
                    dataSets.push_back(new Data);
                    graphSets.push_back(nullptr);
                    dataSets[i]->responceVector.clear();
                    dataSets[i]->concVector.clear();
                    dataSets[i]->volumeVector.clear();
                    dataSets[i]->protein_conc_vector.clear();
                    dataSets[i]->error_vector.clear();
                }
            }
            else if(word[0]=="current_dataset:")
                current_dataSet=toInt(word[1]);

            else if(word[0]=="is_disconnected:")
                is_disconnected=toInt(word[1]);

            else if(word[0]=="dataSet")
                current_dataset_local=toInt(word[1]);

            else if(word[0]=="name:")
                for (unsigned int i=1; i<word.size();i++){
                    if(i==word.size()-1)
                        dataSets[current_dataset_local]->name+=QString::fromStdString(word[i]);
                    else
                        dataSets[current_dataset_local]->name+=QString::fromStdString(word[i])+" ";
                }
            else if(word[0]=="mode:")
                for (unsigned int i=1; i<word.size();i++){
                    if(i==word.size()-1)
                        mode+=QString::fromStdString(word[i]);
                    else
                        mode+=QString::fromStdString(word[i])+" ";
                }

            else if(word[0]=="shape:")
                shape=toInt(word[1]);

            else if(word[0]=="title:"){
                QString s{};
                for (unsigned int j=1;j<word.size();j++)
                    s+=QString::fromStdString(word[j]) + " ";
                title->setText(s);
            }
            else if(word[0]=="title_font:")
                title->setFont(QFont(QString::fromStdString(word[1]),toInt(word[2])));

            else if(word[0]=="a:"){
                dataSets[current_dataset_local]->a.clear();
                for (unsigned int i=1;i<word.size();i++){
                    dataSets[current_dataset_local]->a.push_back(0);
                    toDouble(word[i],dataSets[current_dataset_local]->a[i-1]);
                }
            }

            else if(word[0]=="dyda:"){
                dataSets[current_dataset_local]->dyda.clear();
                for (unsigned int i=1;i<word.size();i++){
                    dataSets[current_dataset_local]->dyda.push_back(0);
                    toDouble(word[i],dataSets[current_dataset_local]->dyda[i-1]);
                }
            }

            else if(word[0]=="protein_conc:")
                toDouble(word[1],dataSets[current_dataset_local]->protein_conc);

            else if(word[0]=="comp_conc:"){
                double i{};
                toDouble(word[1],i);
                dataSets[current_dataset_local]->update_comp_vector(i);
            }
            else if(word[0]=="ligand_cpmg:")
                toDouble(word[1],dataSets[current_dataset_local]->ligand_cpmg);

            else if(word[0]=="ligand_unit:")
                dataSets[current_dataset_local]->ligand_unit=QString::fromStdString(word[1]);

            else if(word[0]=="unit:")
                dataSets[current_dataset_local]->unit=QString::fromStdString(word[1]);

            else if(word[0]=="data_visible:")
                dataSets[current_dataset_local]->data_visible=toInt(word[1]);

            else if(word[0]=="curve_visible:")
                dataSets[current_dataset_local]->curve_visible=toInt(word[1]);

            else if(word[0]=="has_curve:")
                dataSets[current_dataset_local]->has_curve=toInt(word[1]);

            else if(word[0]=="kd_error:")
                toDouble(word[1],dataSets[current_dataset_local]->kd_error);
            else if(word[0]=="chi2:")
                toDouble(word[1],dataSets[current_dataset_local]->chi2);

            else if(word[0]=="num_bind_site:")
                dataSets[current_dataset_local]->num_bind_site=toInt(word[1]);

            else if(word[0]=="pen:")
                dataSets[current_dataset_local]->pen.setColor(QString::fromStdString(word[1]));
            else if(word[0]=="pen_curve:")
                dataSets[current_dataset_local]->pen_curve.setColor(QString::fromStdString(word[1]));

            else if(word[0]=="style:") {
                void setStyles(QVector<Data *> &dataSets, int current_dataset_local, const std::string &_style, const std::string &_size);
                setStyles(dataSets, current_dataset_local, word[1], word[2]);
            }
            else if(word[0]=="model:"){
                QString s{};
                for (unsigned int j=1;j<word.size();j++)
                    s+=QString::fromStdString(word[j])+ " ";
                dataSets[current_dataset_local]->model=s;
            }

            else if(word[0]=="kd_cpmg:")
                toDouble(word[1], dataSets[current_dataset_local]->kd_cpmg);
            else if(word[0]=="kd2_error:")
                toDouble(word[1], dataSets[current_dataset_local]->kd2_error);
            else if(word[0]=="kd1_error:")
                toDouble(word[1], dataSets[current_dataset_local]->kd1_error);
            else if(word[0]=="kdc_error:")
                toDouble(word[1],dataSets[current_dataset_local]->kdc_error);
            else if(word[0]=="calc_data:"){
                dataSets[current_dataset_local]->calc_data[0].clear();
                dataSets[current_dataset_local]->calc_data[1].clear();
                if(mode=="standard")
                    for (int i=0;i<1000; i++){
                        getline(file >> std::ws, line);
                        toDouble(line, dataSets[current_dataset_local], true);
                    }
                else if(mode=="cpmg")
                    for(int i=0;i<dataSets[current_dataset_local]->n_cpmgVector.size();i++){
                        getline(file >> std::ws, line);
                        toDouble(line, dataSets[current_dataset_local], true);
                    }
                if (dataSets.size()==current_dataset_local+1)
                    break;
            }
            else{
                if (mode=="standard")
                    toDouble(line, dataSets[current_dataset_local]);
                else if(mode=="cpmg")
                    toDouble_cpmg(line, dataSets[current_dataset_local]);
            }
        }
        file.close();
    }
    else
        return 1;
    return 0;
}

void setStyles(QVector<Data *> &dataSets, int current_dataset_local, const std::string &_style, const std::string &_size)
{
    int m=toInt(_style);
    if (m==5)
        dataSets[current_dataset_local]->style.setShape(QCPScatterStyle::ssDisc);
    else if(m==6)
        dataSets[current_dataset_local]->style.setShape(QCPScatterStyle::ssSquare);
    else if(m==7)
        dataSets[current_dataset_local]->style.setShape(QCPScatterStyle::ssDiamond);
    else if(m==8)
        dataSets[current_dataset_local]->style.setShape(QCPScatterStyle::ssStar);
    else if(m==2)
        dataSets[current_dataset_local]->style.setShape(QCPScatterStyle::ssCross);
    else if(m==3)
        dataSets[current_dataset_local]->style.setShape(QCPScatterStyle::ssPlus);
    else if(m==4)
        dataSets[current_dataset_local]->style.setShape(QCPScatterStyle::ssCircle);
    else if(m==9)
        dataSets[current_dataset_local]->style.setShape(QCPScatterStyle::ssTriangle);
    else if(m==10)
        dataSets[current_dataset_local]->style.setShape(QCPScatterStyle::ssTriangleInverted);
    else if(m==11)
        dataSets[current_dataset_local]->style.setShape(QCPScatterStyle::ssCrossSquare);
    else if(m==12)
        dataSets[current_dataset_local]->style.setShape(QCPScatterStyle::ssPlusSquare);
    else if(m==13)
        dataSets[current_dataset_local]->style.setShape(QCPScatterStyle::ssCrossCircle);
    else if(m==14)
        dataSets[current_dataset_local]->style.setShape(QCPScatterStyle::ssPlusCircle);
    else if(m==15)
        dataSets[current_dataset_local]->style.setShape(QCPScatterStyle::ssPeace);
    else if(m==17)
        dataSets[current_dataset_local]->style.setShape(QCPScatterStyle::ssCustom);

    dataSets[current_dataset_local]->style.setSize(toInt(_size));
}

void saveImage(Plis *plis, QCustomPlot *plot_area, QVector<Data*> dataSets, int current_dataSet)
{
    QString saveFilename = QFileDialog::getSaveFileName(plis, "Save as", dataSets[current_dataSet]->name, "PNG(*.png);; TIFF(*.tiff *.tif);; "
                                                                                 "JPEG(*.jpg *.jpeg);; PDF(*.pdf)");
    if (saveFilename=="")
        return;
    QPixmap pixmap(plot_area->size());
    plot_area->render(&pixmap, QPoint(), QRegion());
    if (!saveFilename.contains(".pdf"))
        pixmap.save(saveFilename);
    else
        plot_area->savePdf(saveFilename, true);
}


void line2words(const std::string& line, std::vector<std::string> &word)
{
    word.resize(0);
    std::string sub{};
    std::istringstream iss(line);
    while(iss >> sub)
        word.push_back(sub);
}

void toDouble(std::string line, Data* &data, bool curve){
    std::istringstream iss{line};
    if (curve==false){
        double responce, conc, volume, prot_conc, error;
        iss >> responce >> conc >> volume >> prot_conc >> error;
        data->responceVector.push_back(responce);
        data->concVector.push_back(conc);
        data->volumeVector.push_back(volume);
        data->protein_conc_vector.push_back(prot_conc);
        data->error_vector.push_back(error);
    }
    else if(curve==true){
        double responce, conc;
        iss >> conc >> responce;
        data->calc_data[0].push_back(conc);
        data->calc_data[1].push_back(responce);
    }
}

void toDouble_cpmg(std::string line, Data* &data){
    std::istringstream iss{line};
    double n{}, r{},dy{};
    iss >> n >> r >> dy;
    data->pushBackCPMGdata(n, r, dy);
}

void toDouble(const std::string& word, double &temp){
    std::istringstream iss{word};
    iss >> temp;
}

int toInt(const std::string& word){
    std::istringstream iss{word};
    int i{};
    iss >> i;
    return i;
}

