#include "data.h"
#include "plis.h"

// 200612
// - took absolute values of Kd:s when setting results

Data::Data()
{
    concVector.push_back(0.0);
    responceVector.push_back(0.0);
    volumeVector.push_back(0.0);
    protein_conc_vector.push_back(0.0);
    error_vector.push_back(1.e-4);
    unit="uM";
    ligand_unit="uM";
    model="No fit made";
    data_visible=true;
    curve_visible=false;
    has_curve=false;
}

void Data::read_Data(std::string &f){
    std::ifstream file;
    file.open(f.c_str(), std::ios::in);
    responceVector.clear();
    concVector.clear();
    volumeVector.clear();
    protein_conc_vector.clear();
    error_vector.clear();
    if(file.is_open())
    {
        std::string temp{};
        bool stop=false;
        while (stop==false){
            getline(file >> std::ws, temp);
            if (file.eof())
                stop=true;
            std::istringstream iss{temp};
            double response{}, conc{}, volume{};
            iss >> response >> conc >> volume;
            responceVector.push_back(response);
            concVector.push_back(conc);
            volumeVector.push_back(volume);
            protein_conc_vector.push_back(0);
            error_vector.push_back(1.e-4);
            //y_calc.push_back(0);
        }
        file.close();
    }
    else{//öppna felmeddelande går inte att köra Qmessagebox då ui inte nås från data

    }
}

void Data::change_unit_of_data(QString &arg1, QString data){
    if (data!="protein"){
    if (arg1=="uM"){
        if (ligand_unit=="mM")
            change_data(1000, data); //we have uM and want mM
        else if(ligand_unit=="nM")
            change_data(0.001, data);
        else
            std::cout<<"Something is wrong with changing unit of data, tries to convert from uM to uM" << std::endl;
    }
    else if (arg1=="mM"){
        if (ligand_unit=="uM")
            change_data(1000, data); //we have mM and want uM
        else if(ligand_unit=="nM")
            change_data(1000000, data);
        else
            std::cout<<"Something is wrong with changing unit of data, tries to convert from mM to mM" << std::endl;
    }
    else if (arg1=="nM"){
        if (ligand_unit=="uM")
            change_data(0.001, data); //we have nM and want uM
        else if(ligand_unit=="mM")
            change_data(0.000001, data);
        else
            std::cout<<"Something is wrong with changing unit of data, tries to convert from nM to nM" << std::endl;
    }
    }
    else{ // HÄR STÄMMER INTE OMVANDLINGARNA
        if (arg1=="uM"){
            if (unit=="mM")
                change_data(1000, data);
            else if(unit=="nM")
                change_data(0.001, data);
            else
                std::cout<<"Something is wrong with changing protein conc when fitting a curve" << std::endl;
        }
        else if (arg1=="mM"){
            if (unit=="uM")
                change_data(1000, data);
            else if(unit=="nM")
                change_data(1000000, data);
            else
                std::cout<<"Something is wrong with changing protein conc when fitting a curve" << std::endl;
        }
        else if (arg1=="nM"){
            if (unit=="uM")
                change_data(0.001, data);
            else if(unit=="mM")
                change_data(0.000001, data);
            else
                std::cout<<"Something is wrong with changing protein conc when fitting a curve" << std::endl;
        }

    }
}

void Data::change_data(double ratio, QString type){
    if (type=="conc")
        for (int i=0; i<concVector.size();i++)
            concVector[i]*=ratio;
    else if(type=="responce")
        for (int i=0; i<responceVector.size();i++)
            responceVector[i]*=ratio;
    else if(type=="protein")
        for (int i=0; i<protein_conc_vector.size();i++)
            protein_conc_vector[i]*=ratio;
    else if(type=="calc")
        for (int i=0; i<calc_data[0].size();i++)
            calc_data[0][i]*=ratio;
}

void Data::update_protein_conc(int row){
    if (row!=(-1))
        protein_conc_vector[row]=(protein_conc * volumeVector[0])/volumeVector[row];
    else{
        for (int i=0;i<protein_conc_vector.size();i++){
            protein_conc_vector[i]=(protein_conc * volumeVector[0])/volumeVector[i];
        }
    }
}

void Data::update_comp_vector(double start){
    comp_vector.clear();
    for (int i=0;i<protein_conc_vector.size();i++){
        comp_vector.push_back((start * volumeVector[0])/volumeVector[i]);
    }
}

void Data::clearVectors() {
    responceVector.clear();
    concVector.clear();
    volumeVector.clear();
}

void Data::pushBackData(double response, double conc, double volume, double protConc, double err) {
    responceVector.push_back(response);
    concVector.push_back(conc);
    volumeVector.push_back(volume);
    protein_conc_vector.push_back(protConc);
    error_vector.push_back(err);
}

void Data::insertData(int i, double response, double conc, double volume, double protConc, double err) {
    responceVector.insert(i, response);
    concVector.insert(i, conc);
    volumeVector.insert(i, volume);
    protein_conc_vector.insert(i, protConc);
    error_vector.insert(i, err);
}

void Data::removeData(int i) {
    responceVector.remove(i);
    concVector.remove(i);
    volumeVector.remove(i);
    protein_conc_vector.remove(i);
    error_vector.remove(i);
}

void Data::pushBackCPMGdata(double n_cpmg, double r2eff, double dr2eff) {
    n_cpmgVector.push_back(n_cpmg);
    R2effVector.push_back(r2eff);
    dyVector.push_back(dr2eff);
}

void Data::insertCPMGdata(int i, double n_cpmg, double r2eff, double dr2eff) {
    n_cpmgVector.insert(i,n_cpmg);
    R2effVector.insert(i,r2eff);
    dyVector.insert(i,dr2eff);
}

void Data::removeCPMGdata(int i) {
    n_cpmgVector.remove(i);
    R2effVector.remove(i);
    dyVector.remove(i);
}

void Data::setModel_a_dydaOneSite() {
    a = {responceVector[0], responceVector.last(), -1., concVector[0], protein_conc};
    dyda = std::vector<double>(a.size(), 0.);
    model="One bind site";
    num_bind_site=1;
}

void Data::setModel_a_dydaTwoSite() {
    a = {responceVector[0], -1., -1., -1, -1., concVector[0], protein_conc};
    dyda = std::vector<double>(a.size(), 0.);
    model="Two bind site";
    num_bind_site=2;
}

void Data::setModel_a_dydaFourSite() {
    a = {responceVector[0], -1., -1., -1, -1., -1., -1., -1, -1., concVector[0], protein_conc};
    dyda = std::vector<double>(a.size(), 0.);
    model="Four bind site";
    num_bind_site=5;  // note that this is just an index
}

void Data::setModel_a_dydaCompTwoSite(double comp_conc) {
    a = {responceVector[0], -1., -1., -1, -1., concVector[0], protein_conc, comp_conc};
    dyda = std::vector<double>(a.size(), 0.);
    model="Competitive binding";
    num_bind_site=3;
}

void Data::setModel_a_dydaCPMG()
{
    a = {0., 0., 0., 0., ligand_cpmg};
    dyda = {0., 0., 0., 0., 0.};
    model="CPMG-model";
    num_bind_site=4;
}

void Data::setPenColor(int i) {
    switch (i%5) {
        case 0: pen.setColor("black");  break;
        case 1: pen.setColor("red"); break;
        case 2: pen.setColor("green"); break;
        case 3: pen.setColor("blue"); break;
        case 4: pen.setColor("magenta"); break;
    }
}

void Data::setModelIfNoFit()
{
    model="No fit made";
    num_bind_site=0;
}

void Data::addCPMGresults(const std::vector<double> &fitted_a, double fitted_chi2, const std::vector<std::vector<double>> &temp_vector)
{
    a=fitted_a;
    chi2=fitted_chi2;
    kd_cpmg=CPMG_kd(a);
    calc_data[0] = QVector<double>::fromStdVector(temp_vector[0]);
    calc_data[1] = QVector<double>::fromStdVector(temp_vector[1]);
    curve_visible=true;
    has_curve=true;
}

void Data::removeFittedCurve()
{
    has_curve=false;
    curve_visible=false;
    setModelIfNoFit();
    a[2]=0;
}

void Data::setResult(const std::vector<double> &_a, double _chi2, std::vector<double> &_calc_data0, std::vector<double> &_calc_data1)
{
    a=_a;
    // We are fitting absolute values of Kd:s. N.B. CPMG already taken care of
    if (model=="One bind site")
        a[2]=fabs(a[2]);
    else if (model == "Two bind site") {
        a[3]=fabs(a[3]);
        a[4]=fabs(a[4]);
    }
    else if (model == "Four bind site") {
        a[5]=fabs(a[5]);
        a[6]=fabs(a[6]);
        a[7]=fabs(a[7]);
        a[8]=fabs(a[8]);
    }
    else if (model == "Competitive binding") {
        a[2]=fabs(a[2]);
        a[3]=fabs(a[3]);
        a[4]=fabs(a[4]);
    }
    chi2=_chi2;
    calc_data[0] = QVector<double>::fromStdVector(_calc_data0);
    calc_data[1] = QVector<double>::fromStdVector(_calc_data1);
    curve_visible=true;
    has_curve=true;
}

void Data::dataAndCurveVisible(bool dataVis, bool curveVis)
{
    data_visible=dataVis;
    curve_visible=curveVis;
}
