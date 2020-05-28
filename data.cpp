#include "data.h"
#include "plis.h"

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
