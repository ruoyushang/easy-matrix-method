
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <utility>
#include <math.h> 

#include <algorithm>
#include <functional>

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TChain.h"
#include "TBranch.h"
#include "TRandom.h"
#include "TCut.h"

//#include "../../EventDisplay/VEvndispRunParameter.h"
#include "/home/rshang/MatrixDecompositionMethod/EventDisplay/VEvndispRunParameter.h"

#include <complex>
#include "/home/rshang/MatrixDecompositionMethod/Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
#include "/home/rshang/MatrixDecompositionMethod/Eigen/eigen-eigen-323c052e1731/Eigen/StdVector"
using namespace Eigen;

#include "GetRunList.h"

char target[50] = "";
double SizeSecondMax_Cut = 0.;

double map_center_x = 0.;
double map_center_y = 0.;
double mean_tele_point_ra = 0.;
double mean_tele_point_dec = 0.;
double run_tele_point_ra = 0.;
double run_tele_point_dec = 0.;

string SMI_INPUT;
string SMI_OUTPUT;
string SMI_DIR;
string SMI_AUX;

// EventDisplay variables
double TelElevation = 0;
double TelAzimuth = 0;
double TelRAJ2000 = 0;
double TelDecJ2000 = 0;
double ErecS = 0;
double EChi2S = 0;
double MSCW = 0;
int NImages = 0;
double Xcore = 0.;
double Ycore = 0.;
double SizeSecondMax = 0;
double MSCL = 0;
double Time = 0;
int MJD = 0;
UInt_t MJD_UInt_t = 0;
double Shower_Ze = 0;
double Shower_Az = 0;
double EmissionHeight = 0;
double Xoff = 0;
double Yoff = 0;
double Xoff_derot = 0;
double Yoff_derot = 0;
double R2off = 0;
double Phioff = 0;
double theta2 = 0;
double ra_sky = 0;
double dec_sky = 0;
double exposure_hours = 0.;
double exposure_hours_usable = 0.;

vector<TH2D> Hist_OffData_MSCLW;
vector<TH2D> Hist_OffData_SR_XYoff;
vector<TH2D> Hist_OffData_CR_XYoff;
vector<TH2D> Hist_OffData_Ratio_XYoff;
vector<TH2D> Hist_OnData_MSCLW;
vector<TH2D> Hist_OnData_SR_Skymap;
vector<TH2D> Hist_OnData_CR_Skymap;
vector<TH2D> Hist_OnData_SR_Skymap_Sum;
vector<TH2D> Hist_OnData_CR_Skymap_Sum;
vector<vector<TH2D>> Hist_OnData_SR_Skymap_per_sample;
vector<vector<TH2D>> Hist_OnData_CR_Skymap_per_sample;

pair<double,double> GetRunElevAzim(int int_run_number)
{
    double TelElevation_avg = 0.;
    double TelAzimuth_avg = 0.;
    
    char run_number[50];
    sprintf(run_number, "%i", int_run_number);
    string file_name;
    file_name = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");

    if (!gSystem->AccessPathName(file_name.c_str()))
    {
        char run_number[50];
        sprintf(run_number, "%i", int(int_run_number));
        TFile*  input_file = TFile::Open(file_name.c_str());
        TTree* pointing_tree = nullptr;
        pointing_tree = (TTree*) input_file->Get(TString("run_"+string(run_number)+"/stereo/pointingDataReduced"));
        pointing_tree->SetBranchStatus("*",0);
        pointing_tree->SetBranchStatus("TelElevation",1);
        pointing_tree->SetBranchStatus("TelAzimuth",1);
        pointing_tree->SetBranchAddress("TelElevation",&TelElevation);
        pointing_tree->SetBranchAddress("TelAzimuth",&TelAzimuth);
        double total_entries = (double)pointing_tree->GetEntries();
        pointing_tree->GetEntry(int(total_entries/2.));
        TelElevation_avg = TelElevation;
        TelAzimuth_avg = TelAzimuth;
        input_file->Close();
        //std::cout << "root file elev = " << TelElevation_avg << " azim = " << TelAzimuth_avg << std::endl;
    }

    return std::make_pair(TelElevation_avg,TelAzimuth_avg);
}

void RemoveNonExistingRuns(vector<int>& input_list)
{
    vector<int> new_list;
    for (int run=0;run<input_list.size();run++)
    {
        char run_number[50];
        sprintf(run_number, "%i", input_list.at(run));
        string file_name;
        file_name = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");
        if (gSystem->AccessPathName(file_name.c_str()))
        {
            std::cout << "Run " << run_number << " does not exist!!!" << std::endl;
            //auto it = std::find(input_list.begin(), input_list.end(), run);
            //input_list.erase(it);
            std::cout << "Run " << run_number << " removed." << std::endl;
        }
        else
        {
            new_list.push_back(input_list.at(run));
        }
    }
    input_list.clear();
    input_list = new_list;
}

void SortingList(vector<int>* list, vector<double>* list_pointing)
{
    int temp_list;
    double temp_list_pointing;
    for (int run1=0;run1<list->size()-1;run1++)
    {
        for (int run2=run1+1;run2<list->size();run2++)
        {
            if (list_pointing->at(run1)>list_pointing->at(run2))
            {
                temp_list = list->at(run1);
                temp_list_pointing = list_pointing->at(run1);
                list->at(run1) = list->at(run2);
                list_pointing->at(run1) = list_pointing->at(run2);
                list->at(run2) = temp_list;
                list_pointing->at(run2) = temp_list_pointing;
            }
        }
    }
}

vector<double> GetRunElevationList(vector<int> list)
{
    vector<double> list_elev;
    for (int run=0;run<list.size();run++)
    {
        double run_elev = GetRunElevAzim(int(list.at(run))).first;
        double run_azim = GetRunElevAzim(int(list.at(run))).second;
        if (run_azim>180.) run_azim = 360.-run_azim;
        //list_elev.push_back(run_elev);
        list_elev.push_back(run_azim);
    }
    return list_elev;
}

bool PointingSelection(int int_run_number)
{

    char run_number[50];
    sprintf(run_number, "%i", int_run_number);
    string filename;
    filename = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");

    TelElevation = GetRunElevAzim(int_run_number).first;
    TelAzimuth = GetRunElevAzim(int_run_number).second;
    if (TelElevation<Elev_cut_lower) 
    {
        //input_file->Close();
        std::cout << int_run_number << ", elev = " << TelElevation << ", pointing rejected." << std::endl;
        return false;
    }
    if (TelElevation>Elev_cut_upper)
    {
        //input_file->Close();
        std::cout << int_run_number << ", elev = " << TelElevation << ", pointing rejected." << std::endl;
        return false;
    }
    return true;
}

void SetEventDisplayDL3TreeBranch(TTree* Data_tree)
{
    Data_tree->SetBranchStatus("*",0);
    Data_tree->SetBranchStatus("Xoff",1);
    Data_tree->SetBranchStatus("Yoff",1);
    Data_tree->SetBranchStatus("Xderot",1);
    Data_tree->SetBranchStatus("Yderot",1);
    Data_tree->SetBranchStatus("Energy",1);
    Data_tree->SetBranchStatus("Energy_Err",1);
    Data_tree->SetBranchStatus("MSCW",1);
    Data_tree->SetBranchStatus("MSCL",1);
    Data_tree->SetBranchStatus("NImages",1);
    Data_tree->SetBranchStatus("XCore",1);
    Data_tree->SetBranchStatus("YCore",1);
    Data_tree->SetBranchStatus("EmissionHeight",1);
    Data_tree->SetBranchStatus("timeOfDay",1);
    Data_tree->SetBranchStatus("MJD",1);
    Data_tree->SetBranchAddress("Xoff",&Xoff);
    Data_tree->SetBranchAddress("Yoff",&Yoff);
    Data_tree->SetBranchAddress("Xderot",&Xoff_derot);
    Data_tree->SetBranchAddress("Yderot",&Yoff_derot);
    Data_tree->SetBranchAddress("Energy",&ErecS);
    Data_tree->SetBranchAddress("Energy_Err",&EChi2S);
    Data_tree->SetBranchAddress("MSCW",&MSCW);
    Data_tree->SetBranchAddress("MSCL",&MSCL);
    Data_tree->SetBranchAddress("NImages",&NImages);
    Data_tree->SetBranchAddress("XCore",&Xcore);
    Data_tree->SetBranchAddress("YCore",&Ycore);
    Data_tree->SetBranchAddress("EmissionHeight",&EmissionHeight);
    Data_tree->SetBranchAddress("timeOfDay",&Time);
    Data_tree->SetBranchAddress("MJD",&MJD);
}

vector<pair<double,double>> GetRunTimecuts(int run_number)
{
    string line;
    char delimiter = ' ';
    string acc_runnumber = "";
    string acc_timecuts = "";
    vector<pair<double,double>> timecuts;
    int nth_delimiter = 0;
    std::string::size_type sz;

    // run veritas_db_query.py to get this file
    ifstream myfile (SMI_AUX+"/timecuts_allruns.txt");
    if (myfile.is_open())
    {
        while ( getline(myfile,line) )
        {
            acc_runnumber = "";
            nth_delimiter = 0;
            for(int i = 0; i < line.size(); i++)
            {
                if(line[i] == delimiter)
                {
                    nth_delimiter += 1;
                }
                else if (nth_delimiter==0)
                {
                    acc_runnumber += line[i];
                }
            }
            if (std::stoi(acc_runnumber,nullptr,10)==run_number)
            {
                //std::cout << "find run " << run_number << std::endl;
                nth_delimiter = 0;
                acc_timecuts = "";
                for(int i = 0; i < line.size(); i++)
                {
                    if (line[i] != delimiter)
                    {
                        acc_timecuts += line[i];
                    }
                    else
                    {
                        size_t pos = 0;
                        string time_delimiter = "/";
                        if ((pos = acc_timecuts.find(time_delimiter)) != std::string::npos)
                        {
                            string timecut_start = acc_timecuts.substr(0, pos);
                            acc_timecuts.erase(0, pos + time_delimiter.length());
                            string timecut_end = acc_timecuts;
                            //std::cout << "timecut_start " << timecut_start << std::endl;
                            //std::cout << "timecut_end " << timecut_end << std::endl;
                            double val_timecut_start = std::stod(timecut_start,&sz);
                            double val_timecut_end = std::stod(timecut_end,&sz);
                            timecuts.push_back(std::make_pair(val_timecut_start,val_timecut_end));
                        }
                        acc_timecuts = "";
                        nth_delimiter += 1;
                    }
                }
                size_t pos = 0;
                string time_delimiter = "/";
                if ((pos = acc_timecuts.find(time_delimiter)) != std::string::npos)
                {
                    string timecut_start = acc_timecuts.substr(0, pos);
                    acc_timecuts.erase(0, pos + time_delimiter.length());
                    string timecut_end = acc_timecuts;
                    //std::cout << "timecut_start " << timecut_start << std::endl;
                    //std::cout << "timecut_end " << timecut_end << std::endl;
                    double val_timecut_start = std::stod(timecut_start,&sz);
                    double val_timecut_end = std::stod(timecut_end,&sz);
                    timecuts.push_back(std::make_pair(val_timecut_start,val_timecut_end));
                }
                break;
            }
        }
        myfile.close();
    }
    else std::cout << "Unable to open file timecuts_allruns.txt" << std::endl; 

    return timecuts;
}

double GetRunUsableTime(string file_name,int run_number)
{
    string line;
    char delimiter = ' ';
    string acc_runnumber = "";
    string acc_time = "";
    int nth_delimiter = 0;
    std::string::size_type sz;
    double usable_time = 0.;

    // run veritas_db_query.py to get this file
    ifstream myfile (SMI_AUX+"/usable_time_allruns.txt");
    if (myfile.is_open())
    {
        while ( getline(myfile,line) )
        {
            acc_runnumber = "";
            nth_delimiter = 0;
            for(int i = 0; i < line.size(); i++)
            {
                if(line[i] == delimiter)
                {
                    nth_delimiter += 1;
                }
                else if (nth_delimiter==0)
                {
                    acc_runnumber += line[i];
                }
            }
            if (std::stoi(acc_runnumber,nullptr,10)==run_number)
            {
                //std::cout << "find run " << run_number << std::endl;
                nth_delimiter = 0;
                acc_time = "";
                for(int i = 0; i < line.size(); i++)
                {
                    if (line[i] != delimiter)
                    {
                        acc_time += line[i];
                    }
                    else
                    {
                        acc_time = "";
                        nth_delimiter += 1;
                    }
                }
                usable_time = std::stod(acc_time,&sz);
                //size_t pos = 0;
                //string time_delimiter = ":";
                //if ((pos = acc_time.find(time_delimiter)) != std::string::npos)
                //{
                //    pos = acc_time.find(time_delimiter);
                //    string time_hour = acc_time.substr(0, pos);
                //    acc_time.erase(0, pos + time_delimiter.length());
                //    pos = acc_time.find(time_delimiter);
                //    string time_minute = acc_time.substr(0, pos);
                //    acc_time.erase(0, pos + time_delimiter.length());
                //    string time_second = acc_time;
                //    usable_time = std::stod(time_hour,&sz)*60.*60.;
                //    usable_time += std::stod(time_minute,&sz)*60.;
                //    usable_time += std::stod(time_second,&sz);
                //}
                break;
            }
        }
        myfile.close();
    }
    else std::cout << "Unable to open file usable_time_allruns.txt" << std::endl; 

    if (usable_time==0.)
    {
        char run_number_char[50];
        sprintf(run_number_char, "%i", run_number);
        TFile*  input_file = TFile::Open(file_name.c_str());
        TTree* pointing_tree = nullptr;
        pointing_tree = (TTree*) input_file->Get(TString("run_"+string(run_number_char)+"/stereo/pointingDataReduced"));
        if (!pointing_tree)
        {
            return 0.;
        }
        pointing_tree->SetBranchStatus("*",0);
        pointing_tree->SetBranchStatus("Time",1);
        pointing_tree->SetBranchAddress("Time",&Time);
        double total_entries = (double)pointing_tree->GetEntries();
        pointing_tree->GetEntry(0);
        int time_start = Time;
        pointing_tree->GetEntry(total_entries-1);
        int time_end = Time;
        input_file->Close();
        usable_time = time_end-time_start;
    }
    return usable_time;
}

bool SelectNImages()
{
    if (NImages<min_NImages) return false;
    if (EmissionHeight<EmissionHeight_cut) return false;
    if (Xcore*Xcore+Ycore*Ycore>max_Rcore*max_Rcore) return false;
    if (pow(R2off,0.5)>1.6) return false;
    return true;
}

bool ApplyTimeCuts(double event_time, vector<pair<double,double>> timecut)
{
    bool pass_cut = true;
    for (int cut=0;cut<timecut.size();cut++)
    {
        double timecut_start = timecut.at(cut).first;
        double timecut_end = timecut.at(cut).second;
        //std::cout << "event_time " << event_time << ", timecut_start " << timecut_start << ", timecut_end " << timecut_end << std::endl;
        if (event_time>timecut_start && event_time<timecut_end)
        {
            pass_cut = false;
        }
    }
    return pass_cut;
}

pair<double,double> GetRunRaDec(string file_name, int run)
{

    double TelRAJ2000_tmp = 0.;
    double TelDecJ2000_tmp = 0.;

    char run_number[50];
    sprintf(run_number, "%i", int(run));
    TFile*  input_file = TFile::Open(file_name.c_str());
    TTree* pointing_tree = nullptr;
    pointing_tree = (TTree*) input_file->Get(TString("run_"+string(run_number)+"/stereo/pointingDataReduced"));
    if (!pointing_tree)
    {
        return std::make_pair(0.,0.);
    }
    pointing_tree->SetBranchStatus("*",0);
    pointing_tree->SetBranchStatus("TelRAJ2000",1);
    pointing_tree->SetBranchStatus("TelDecJ2000",1);
    pointing_tree->SetBranchAddress("TelRAJ2000",&TelRAJ2000);
    pointing_tree->SetBranchAddress("TelDecJ2000",&TelDecJ2000);
    double total_entries = (double)pointing_tree->GetEntries();
    pointing_tree->GetEntry(int(total_entries/2.));
    TelRAJ2000_tmp = TelRAJ2000*180./M_PI;
    TelDecJ2000_tmp = TelDecJ2000*180./M_PI;
    input_file->Close();

    return std::make_pair(TelRAJ2000_tmp,TelDecJ2000_tmp);
}

void TrainingRunAnalysis(int int_run_number)
{

    char run_number[50];
    sprintf(run_number, "%i", int_run_number);
    string filename;
    filename = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");

    TFile*  input_file = TFile::Open(filename.c_str());
    TString root_file = "run_"+string(run_number)+"/stereo/DL3EventTree";
    TTree* Data_tree = (TTree*) input_file->Get(root_file);
    if (!Data_tree)
    {
        std::cout << "TTree does not exist: " << root_file << std::endl;
        return;
    }
    SetEventDisplayDL3TreeBranch(Data_tree);
    //std::cout << "Get time cuts for run " << run_number << std::endl;
    vector<pair<double,double>> timecut_thisrun = GetRunTimecuts(int_run_number);
    //std::cout << "Get elev. and azim. for run " << run_number << std::endl;
    double tele_elev = GetRunElevAzim(int_run_number).first;

    //std::cout << "Get telescope pointing RA and Dec for run " << run_number << std::endl;
    std::pair<double,double> tele_point_ra_dec = GetRunRaDec(filename,int_run_number);
    run_tele_point_ra = tele_point_ra_dec.first;
    run_tele_point_dec = tele_point_ra_dec.second;

    Data_tree->GetEntry(0);
    double time_0 = Time;
    Data_tree->GetEntry(Data_tree->GetEntries()-1);
    double time_1 = Time;
    exposure_hours += (time_1-time_0)/3600.;
    //std::cout << "Get usable time for run " << run_number << std::endl;
    double exposure_thisrun = GetRunUsableTime(filename,int_run_number)/3600.;
    exposure_hours_usable += exposure_thisrun;

    double MSCW_plot_upper = gamma_hadron_dim_ratio_w*(MSCW_cut_blind-(-1.*MSCW_cut_blind))+MSCW_cut_blind;
    double MSCL_plot_upper = gamma_hadron_dim_ratio_l*(MSCL_cut_blind-(-1.*MSCL_cut_blind))+MSCL_cut_blind;
    double MSCW_plot_lower = -MSCW_cut_blind;
    double MSCL_plot_lower = -MSCL_cut_blind;

    TH1D Hist_ErecS = TH1D("Hist_ErecS","",N_energy_bins,energy_bins);
    for (int entry=0;entry<Data_tree->GetEntries();entry++) 
    {
        ErecS = 0;
        EChi2S = 0;
        NImages = 0;
        Xcore = 0;
        Ycore = 0;
        SizeSecondMax = 0;
        MSCW = 0;
        MSCL = 0;
        R2off = 0;
        Data_tree->GetEntry(entry);
        R2off = Xoff*Xoff+Yoff*Yoff;
        Phioff = atan2(Yoff,Xoff)+M_PI;
        ra_sky = tele_point_ra_dec.first+Xoff_derot;
        dec_sky = tele_point_ra_dec.second+Yoff_derot;
        theta2 = pow(ra_sky-mean_tele_point_ra,2)+pow(dec_sky-mean_tele_point_dec,2);
        int energy_idx = Hist_ErecS.FindBin(ErecS*1000.)-1;
        if (!SelectNImages()) continue;
        if (!ApplyTimeCuts(Time-time_0, timecut_thisrun)) continue;
        if (energy_idx<0) continue;
        if (energy_idx>=N_energy_bins) continue;
        if (MSCL>MSCL_plot_upper) continue;
        if (MSCL<MSCL_plot_lower) continue;
        if (MSCW>MSCW_plot_upper) continue;
        if (MSCW<MSCW_plot_lower) continue;

        Hist_OffData_MSCLW.at(energy_idx).Fill(MSCL,MSCW);
        if (MSCL<MSCL_cut_blind && MSCW<MSCW_cut_blind)
        {
            Hist_OffData_SR_XYoff.at(energy_idx).Fill(Xoff,Yoff);
        }
        else
        {
            Hist_OffData_CR_XYoff.at(energy_idx).Fill(Xoff,Yoff);
        }

    }
    input_file->Close();

}

void SingleRunAnalysis(int int_run_number)
{

    char run_number[50];
    sprintf(run_number, "%i", int_run_number);
    string filename;
    filename = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");

    TFile*  input_file = TFile::Open(filename.c_str());
    TString root_file = "run_"+string(run_number)+"/stereo/DL3EventTree";
    TTree* Data_tree = (TTree*) input_file->Get(root_file);
    if (!Data_tree)
    {
        std::cout << "TTree does not exist: " << root_file << std::endl;
        return;
    }
    SetEventDisplayDL3TreeBranch(Data_tree);
    //std::cout << "Get time cuts for run " << run_number << std::endl;
    vector<pair<double,double>> timecut_thisrun = GetRunTimecuts(int_run_number);
    //std::cout << "Get elev. and azim. for run " << run_number << std::endl;
    double tele_elev = GetRunElevAzim(int_run_number).first;

    //std::cout << "Get telescope pointing RA and Dec for run " << run_number << std::endl;
    std::pair<double,double> tele_point_ra_dec = GetRunRaDec(filename,int_run_number);
    run_tele_point_ra = tele_point_ra_dec.first;
    run_tele_point_dec = tele_point_ra_dec.second;

    Data_tree->GetEntry(0);
    double time_0 = Time;
    Data_tree->GetEntry(Data_tree->GetEntries()-1);
    double time_1 = Time;
    exposure_hours += (time_1-time_0)/3600.;
    //std::cout << "Get usable time for run " << run_number << std::endl;
    double exposure_thisrun = GetRunUsableTime(filename,int_run_number)/3600.;
    exposure_hours_usable += exposure_thisrun;

    double MSCW_plot_upper = gamma_hadron_dim_ratio_w*(MSCW_cut_blind-(-1.*MSCW_cut_blind))+MSCW_cut_blind;
    double MSCL_plot_upper = gamma_hadron_dim_ratio_l*(MSCL_cut_blind-(-1.*MSCL_cut_blind))+MSCL_cut_blind;
    double MSCW_plot_lower = -MSCW_cut_blind;
    double MSCL_plot_lower = -MSCL_cut_blind;

    TH1D Hist_ErecS = TH1D("Hist_ErecS","",N_energy_bins,energy_bins);
    TH1D Hist_Xoff = TH1D("Hist_Xoff","",N_Xoff_bins,-2.,2.);
    for (int entry=0;entry<Data_tree->GetEntries();entry++) 
    {
        ErecS = 0;
        EChi2S = 0;
        NImages = 0;
        Xcore = 0;
        Ycore = 0;
        SizeSecondMax = 0;
        MSCW = 0;
        MSCL = 0;
        R2off = 0;
        Data_tree->GetEntry(entry);
        R2off = Xoff*Xoff+Yoff*Yoff;
        Phioff = atan2(Yoff,Xoff)+M_PI;
        ra_sky = tele_point_ra_dec.first+Xoff_derot;
        dec_sky = tele_point_ra_dec.second+Yoff_derot;
        theta2 = pow(ra_sky-mean_tele_point_ra,2)+pow(dec_sky-mean_tele_point_dec,2);
        int energy_idx = Hist_ErecS.FindBin(ErecS*1000.)-1;
        int idx_xoff = Hist_Xoff.FindBin(Xoff)-1;
        int idx_yoff = Hist_Xoff.FindBin(Yoff)-1;
        if (!SelectNImages()) continue;
        if (!ApplyTimeCuts(Time-time_0, timecut_thisrun)) continue;
        if (energy_idx<0) continue;
        if (energy_idx>=N_energy_bins) continue;
        if (MSCL>MSCL_plot_upper) continue;
        if (MSCL<MSCL_plot_lower) continue;
        if (MSCW>MSCW_plot_upper) continue;
        if (MSCW<MSCW_plot_lower) continue;

        double map_right_edge = map_center_x+Skymap_size_x;
        double map_left_edge = map_center_x-Skymap_size_x;
        double map_top_edge = map_center_y+Skymap_size_y;
        double map_bottom_edge = map_center_y-Skymap_size_y;
        if (ra_sky>map_right_edge) continue;
        if (ra_sky<map_left_edge) continue;
        if (dec_sky>map_top_edge) continue;
        if (dec_sky<map_bottom_edge) continue;

        Hist_OnData_MSCLW.at(energy_idx).Fill(MSCL,MSCW);
        if (MSCL<MSCL_cut_blind && MSCW<MSCW_cut_blind)
        {
            Hist_OnData_SR_Skymap.at(energy_idx).Fill(ra_sky,dec_sky);
        }
        else
        {
            int xoff_idx = Hist_OffData_Ratio_XYoff.at(energy_idx).GetXaxis()->FindBin(Xoff);
            int yoff_idx = Hist_OffData_Ratio_XYoff.at(energy_idx).GetYaxis()->FindBin(Yoff);
            double weight = Hist_OffData_Ratio_XYoff.at(energy_idx).GetBinContent(xoff_idx,yoff_idx);
            Hist_OnData_CR_Skymap.at(energy_idx).Fill(ra_sky,dec_sky,weight);
        }

    }
    input_file->Close();

}

void FillHistograms(string target_data, bool isON, bool doImposter)
{

    sprintf(target, "%s", target_data.c_str());

    SMI_INPUT = string(std::getenv("SMI_INPUT"));
    SMI_OUTPUT = string(std::getenv("SMI_OUTPUT"));
    SMI_DIR = string(std::getenv("SMI_DIR"));
    SMI_AUX = string(std::getenv("SMI_AUX"));

    TString ONOFF_tag;
    if (isON) 
    {
        source_theta_cut = 0.;
        ONOFF_tag = "ON";
    }
    else
    {
        ONOFF_tag = "OFF";
    }

    SizeSecondMax_Cut = 600.;
    if (TString(target).Contains("V4")) SizeSecondMax_Cut = 400.;
    if (TString(target).Contains("V5")) SizeSecondMax_Cut = 400.;

    std::istringstream ss_input;
    ss_input.str(target);
    std::string source_strip("");
    std::string source_off_strip("");
    std::string segment;
    std::vector<std::string> seglist;
    while(std::getline(ss_input, segment, '_'))
    {
        seglist.push_back(segment);
    }
    for (int s=0;s<seglist.size()-1;s++)
    {
        if (s>0)
        {
            source_strip.append("_");
            source_off_strip.append("_");
        }
        source_strip.append(seglist.at(s));
        source_off_strip.append(seglist.at(s));
    }
    source_off_strip.append("_OFF");

    pair<double,double> source_ra_dec = GetSourceRaDec(TString(target));
    mean_tele_point_ra = source_ra_dec.first;
    mean_tele_point_dec = source_ra_dec.second;
    map_center_x = mean_tele_point_ra;
    map_center_y = mean_tele_point_dec;

    vector<int> Data_runlist_init = GetRunList(source_strip);
    std::cout << "Data_runlist_init.size() = " << Data_runlist_init.size() << std::endl;
    RemoveNonExistingRuns(Data_runlist_init);
    std::cout << "Data_runlist.size() = " << Data_runlist_init.size() << std::endl;

    std::cout << "Getting runlist elevation..." << std::endl;
    vector<double> Data_runlist_init_elev = GetRunElevationList(Data_runlist_init);
    if (Data_runlist_init.size()==0)
    {
        return;
    }
    std::cout << "Sorting On list by elevation..." << std::endl;
    SortingList(&Data_runlist_init,&Data_runlist_init_elev);

    std::cout <<__LINE__ << std::endl;

    double MSCW_plot_upper = gamma_hadron_dim_ratio_w*(MSCW_cut_blind-(-1.*MSCW_cut_blind))+MSCW_cut_blind;
    double MSCL_plot_upper = gamma_hadron_dim_ratio_l*(MSCL_cut_blind-(-1.*MSCL_cut_blind))+MSCL_cut_blind;
    double MSCW_plot_lower = -MSCW_cut_blind;
    double MSCL_plot_lower = -MSCL_cut_blind;
    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        Hist_OnData_MSCLW.push_back(TH2D("Hist_OnData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",mtx_dim_l,MSCL_plot_lower,MSCL_plot_upper,mtx_dim_w,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OffData_MSCLW.push_back(TH2D("Hist_OffData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",mtx_dim_l,MSCL_plot_lower,MSCL_plot_upper,mtx_dim_w,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnData_SR_Skymap.push_back(TH2D("Hist_OnData_SR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_CR_Skymap.push_back(TH2D("Hist_OnData_CR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_SR_Skymap_Sum.push_back(TH2D("Hist_OnData_SR_Skymap_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_CR_Skymap_Sum.push_back(TH2D("Hist_OnData_CR_Skymap_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OffData_SR_XYoff.push_back(TH2D("Hist_OffData_SR_XYoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,-2.,2.,20,-2.,2.));
        Hist_OffData_CR_XYoff.push_back(TH2D("Hist_OffData_CR_XYoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,-2.,2.,20,-2.,2.));
        Hist_OffData_Ratio_XYoff.push_back(TH2D("Hist_OffData_Ratio_XYoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",20,-2.,2.,20,-2.,2.));
    }

    int n_training_samples = 0;
    int nbins_unblind = 0;
    for (int binx=0;binx<mtx_dim_l;binx++)
    {
        for (int biny=0;biny<mtx_dim_w;biny++)
        {
            double bin_center_x = Hist_OnData_MSCLW.at(0).GetXaxis()->GetBinCenter(binx+1);
            double bin_center_y = Hist_OnData_MSCLW.at(0).GetYaxis()->GetBinCenter(biny+1);
            if (bin_center_x>MSCL_cut_blind || bin_center_y>MSCW_cut_blind)
            {
                nbins_unblind += 1;
            }
        }
    }

    std::cout << __LINE__ << std::endl;

    vector<vector<double>> unblinded_elements_per_sample;
    vector<double> unblinded_elements_per_energy;
    vector<double> blinded_element_per_sample;
    double blinded_element_per_energy = 0.;

    int n_samples = 0;

    for (int run=0;run<Data_runlist_init.size();run++)
    {

        if (!PointingSelection(int(Data_runlist_init[run])))
        {
            continue;
        }

        std::cout << "Processing " << run << "/" << Data_runlist_init.size() << "..." << std::endl;
        std::cout << "Getting training data..." << std::endl;
        vector<int> OffData_runlist_init = GetOffRunList(source_strip,int(Data_runlist_init[run]),-1);
        RemoveNonExistingRuns(OffData_runlist_init);

        if (OffData_runlist_init.size()<mtx_dim_w*mtx_dim_l)
        {
            std::cout << "Insufficient training data..." << std::endl;
            continue;
        }

        vector<vector<vector<double>>> off_unblinded_elements;
        vector<vector<double>> off_unblinded_elements_per_sample;
        vector<double> off_unblinded_elements_per_energy;
        vector<vector<double>> off_blinded_element;
        vector<double> off_blinded_element_per_sample;
        double off_blinded_element_per_energy = 0.;

        for (int offrun=0;offrun<OffData_runlist_init.size();offrun++)
        {

            if (!PointingSelection(int(OffData_runlist_init[offrun])))
            {
                continue;
            }

            TrainingRunAnalysis(int(OffData_runlist_init[offrun]));

            off_unblinded_elements_per_sample.clear();
            off_blinded_element_per_sample.clear();
            for (int e=0;e<N_energy_bins;e++) 
            {
                off_unblinded_elements_per_energy.clear();
                off_blinded_element_per_energy = 0.;
                for (int binx=0;binx<mtx_dim_l;binx++)
                {
                    for (int biny=0;biny<mtx_dim_w;biny++)
                    {
                        double bin_center_x = Hist_OffData_MSCLW.at(e).GetXaxis()->GetBinCenter(binx+1);
                        double bin_center_y = Hist_OffData_MSCLW.at(e).GetYaxis()->GetBinCenter(biny+1);
                        double bin_content = Hist_OffData_MSCLW.at(e).GetBinContent(binx+1,biny+1);
                        if (bin_center_x>MSCL_cut_blind || bin_center_y>MSCW_cut_blind)
                        {
                            off_unblinded_elements_per_energy.push_back(bin_content);
                        }
                        else
                        {
                            off_blinded_element_per_energy += bin_content;
                        }
                    }
                }
                off_unblinded_elements_per_sample.push_back(off_unblinded_elements_per_energy);
                off_blinded_element_per_sample.push_back(off_blinded_element_per_energy);
            }
            off_unblinded_elements.push_back(off_unblinded_elements_per_sample);
            off_blinded_element.push_back(off_blinded_element_per_sample);

            for (int e=0;e<N_energy_bins;e++) 
            {
                Hist_OffData_MSCLW.at(e).Reset();
            }
        }

        for (int e=0;e<N_energy_bins;e++) 
        {
            Hist_OffData_Ratio_XYoff.at(e).Add(&Hist_OffData_SR_XYoff.at(e));
            Hist_OffData_Ratio_XYoff.at(e).Divide(&Hist_OffData_CR_XYoff.at(e));
            Hist_OffData_SR_XYoff.at(e).Reset();
            Hist_OffData_CR_XYoff.at(e).Reset();
        }

        std::cout << "Using training data to learn conversion..." << std::endl;

        n_training_samples = off_blinded_element.size();
        std::cout << "n_training_samples = " << n_training_samples << std::endl;
        std::cout << "nbins_unblind = " << nbins_unblind << std::endl;

        // simple ratio method
        vector<double> convert_unblind_to_blind_ratio;
        for (int e=0;e<N_energy_bins;e++) 
        {
            double CR_count = 0.;
            double SR_count = 0.;
            for (int sample=0;sample<n_training_samples;sample++)
            {
                SR_count += off_blinded_element.at(sample).at(e);
                for (int bin=0;bin<nbins_unblind;bin++)
                {
                    CR_count += off_unblinded_elements.at(sample).at(e).at(bin);
                }
            }
            double SR_CR_ratio = 0.;
            if (CR_count>0.)
            {
                SR_CR_ratio = SR_count/CR_count;
            }
            convert_unblind_to_blind_ratio.push_back(SR_CR_ratio);
        }
        

        // linear regression method
        // A * X = B, solve X
        vector<VectorXcd> convert_unblind_to_blind_regression;
        for (int e=0;e<N_energy_bins;e++) 
        {
            VectorXcd vtr_B = VectorXcd::Zero(n_training_samples);
            MatrixXcd mtx_A = MatrixXcd::Zero(n_training_samples,nbins_unblind);
            VectorXcd vtr_X = VectorXcd::Zero(nbins_unblind);

            for (int sample=0;sample<n_training_samples;sample++)
            {
                vtr_B(sample) = off_blinded_element.at(sample).at(e);
                for (int bin=0;bin<nbins_unblind;bin++)
                {
                    mtx_A(sample,bin) = off_unblinded_elements.at(sample).at(e).at(bin);
                }
            }

            BDCSVD<MatrixXcd> bdc_svd(mtx_A, ComputeThinU | ComputeThinV);
            vtr_X = bdc_svd.solve(vtr_B);
            convert_unblind_to_blind_regression.push_back(vtr_X);
        }

        std::cout << "Getting ON data..." << std::endl;

        SingleRunAnalysis(int(Data_runlist_init[run]));
        n_samples += 1;
        unblinded_elements_per_sample.clear();
        blinded_element_per_sample.clear();
        for (int e=0;e<N_energy_bins;e++) 
        {
            unblinded_elements_per_energy.clear();
            blinded_element_per_energy = 0.;
            for (int binx=0;binx<mtx_dim_l;binx++)
            {
                for (int biny=0;biny<mtx_dim_w;biny++)
                {
                    double bin_center_x = Hist_OnData_MSCLW.at(e).GetXaxis()->GetBinCenter(binx+1);
                    double bin_center_y = Hist_OnData_MSCLW.at(e).GetYaxis()->GetBinCenter(biny+1);
                    double bin_content = Hist_OnData_MSCLW.at(e).GetBinContent(binx+1,biny+1);
                    if (bin_center_x>MSCL_cut_blind || bin_center_y>MSCW_cut_blind)
                    {
                        unblinded_elements_per_energy.push_back(bin_content);
                    }
                    else
                    {
                        blinded_element_per_energy += bin_content;
                    }
                }
            }
            unblinded_elements_per_sample.push_back(unblinded_elements_per_energy);
            blinded_element_per_sample.push_back(blinded_element_per_energy);
        }

        std::cout << "Applying conversion to ON data..." << std::endl;


        // A * X = B
        for (int e=0;e<N_energy_bins;e++) 
        {
            // linear regression method
            MatrixXcd mtx_A = MatrixXcd::Zero(1,nbins_unblind);

            for (int bin=0;bin<nbins_unblind;bin++)
            {
                mtx_A(0,bin) = unblinded_elements_per_sample.at(e).at(bin);
            }

            VectorXcd vtr_predict = mtx_A*convert_unblind_to_blind_regression.at(e);
            VectorXcd vtr_truth = VectorXcd::Zero(1); 
            vtr_truth(0) = blinded_element_per_sample.at(e);

            double SR_predict_regression = 0.;
            double total_truth = 0.;
            SR_predict_regression = vtr_predict(0).real();
            total_truth = vtr_truth(0).real();

            double CR_count_map = Hist_OnData_CR_Skymap.at(e).Integral();
            if (CR_count_map>0.)
            {
                Hist_OnData_CR_Skymap.at(e).Scale(SR_predict_regression/CR_count_map);
            }

            // simple ratio method
            double CR_count = 0.;
            for (int bin=0;bin<nbins_unblind;bin++)
            {
                CR_count += unblinded_elements_per_sample.at(e).at(bin);
            }
            double SR_predict_ratio = convert_unblind_to_blind_ratio.at(e)*CR_count; 
            std::cout << "=============================================" << std::endl;
            std::cout << "Energy = " << energy_bins[e] << std::endl;
            std::cout << "Truth = " << total_truth << std::endl;
            std::cout << "Regression predict = " << SR_predict_regression << std::endl;
            std::cout << "Ratio predict = " << SR_predict_ratio << std::endl;
        }


        for (int e=0;e<N_energy_bins;e++) 
        {
            Hist_OnData_SR_Skymap_Sum.at(e).Add(&Hist_OnData_SR_Skymap.at(e));
            Hist_OnData_CR_Skymap_Sum.at(e).Add(&Hist_OnData_CR_Skymap.at(e));
        }

        for (int e=0;e<N_energy_bins;e++) 
        {
            Hist_OnData_MSCLW.at(e).Reset();
            Hist_OnData_SR_Skymap.at(e).Reset();
            Hist_OnData_CR_Skymap.at(e).Reset();
        }

    }


    TFile OutputFile(TString(SMI_OUTPUT)+"/Netflix_"+TString(target)+".root","recreate");
    for (int e=0;e<N_energy_bins;e++) 
    {
        std::cout << "=============================================" << std::endl;
        std::cout << "Energy = " << energy_bins[e] << std::endl;
        std::cout << "Hist_OnData_SR_Skymap_Sum.at(e).Integral() = " << Hist_OnData_SR_Skymap_Sum.at(e).Integral() << std::endl;
        std::cout << "Hist_OnData_CR_Skymap_Sum.at(e).Integral() = " << Hist_OnData_CR_Skymap_Sum.at(e).Integral() << std::endl;
        Hist_OnData_SR_Skymap_Sum.at(e).Write();
        Hist_OnData_CR_Skymap_Sum.at(e).Write();
    }
    OutputFile.Close();

}
