
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
//#include "/home/rshang/MatrixDecompositionMethod/EventDisplay/VEvndispRunParameter.h"

#include <complex>
#include "/home/rshang/MatrixDecompositionMethod/Eigen/eigen-eigen-323c052e1731/Eigen/Dense"
#include "/home/rshang/MatrixDecompositionMethod/Eigen/eigen-eigen-323c052e1731/Eigen/StdVector"
using namespace Eigen;

#include "GetRunList.h"
#include "ResetPublicVariables.C"

char target[50] = "";
double SizeSecondMax_Cut = 0.;

double map_center_x = 0.;
double map_center_y = 0.;
double mean_tele_point_ra = 0.;
double mean_tele_point_dec = 0.;

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
double MVA = 0;
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
double complementary_ra_sky = 0;
double complementary_dec_sky = 0;
double exposure_hours = 0.;
double exposure_hours_usable = 0.;
double total_cr_count = 0.;
double off_cr_count = 0.;
double avg_diff_nsb = 0.;
double avg_diff_el = 0.;
double avg_diff_az = 0.;
double NSB_mean = 0.;
double Elev_mean = 0.;
double Azim_mean = 0.;

vector<double> effective_area;
vector<double> data_count;
vector<double> ratio_bkgd_count;
vector<double> regression_bkgd_count;
vector<double> init_perturbation_bkgd_count;
vector<double> perturbation_bkgd_count;
vector<double> combined_bkgd_count;
vector<double> t00_truth;
vector<double> t01_truth;
vector<double> t10_truth;
vector<double> t11_truth;
vector<double> t00_recon;
vector<double> t01_recon;
vector<double> t10_recon;
vector<double> t11_recon;

int binx_blind_upper_global;
int biny_blind_upper_global;
int binx_blind_lower_global;
int biny_blind_lower_global;

TH2D Hist_Data_AreaTime_Skymap;
TH2D Hist_Data_Norm_Skymap;
TH2D Hist_Data_Elev_Skymap;
TH2D Hist_Data_Azim_Skymap;
TH2D Hist_Data_NSB_Skymap;

vector<TH2D> Hist_OnData_MSCLW;
vector<TH2D> Hist_OffData_MSCLW;
vector<TH2D> Hist_OnData_MSCLW_Fine;
vector<TH2D> Hist_OffData_MSCLW_Fine;
vector<TH2D> Hist_OnBkgd_MSCLW_Fine;
vector<TH2D> Hist_OffData_MSCLW_Sum;
vector<TH2D> Hist_OffData_MSCLW_Fine_Sum;
vector<TH2D> Hist_OnData_SR_XYoff;
vector<TH2D> Hist_OnData_CR_XYoff;
vector<TH2D> Hist_OffData_SR_XYoff;
vector<TH2D> Hist_OffData_CR_XYoff;
vector<TH2D> Hist_OffData_Ratio_XYoff;
vector<TH2D> Hist_OnData_SR_Skymap;
vector<TH2D> Hist_OnData_SR_Skymap_Mask;
vector<TH2D> Hist_OnData_CR_Skymap_Mask;
vector<TH2D> Hist_OnData_CR_Skymap_FoV;
vector<TH2D> Hist_OnData_CR_Skymap_Ratio;
vector<TH2D> Hist_OnData_CR_Skymap_Regression;
vector<TH2D> Hist_OnData_CR_Skymap_Init_Perturbation;
vector<TH2D> Hist_OnData_CR_Skymap_Perturbation;
vector<TH2D> Hist_OnData_Expo_Skymap_Sum;
vector<TH2D> Hist_OnData_SR_Skymap_Sum;
vector<TH2D> Hist_OnData_CR_Skymap_FoV_Sum;
vector<TH2D> Hist_OnData_CR_Skymap_Ratio_Sum;
vector<TH2D> Hist_OnData_CR_Skymap_Regression_Sum;
vector<TH2D> Hist_OnData_CR_Skymap_Init_Perturbation_Sum;
vector<TH2D> Hist_OnData_CR_Skymap_Perturbation_Sum;
vector<TH2D> Hist_OnData_CR_Skymap_Combined_Sum;
vector<TH2D> Hist_OnData_MSCLW_Fine_Sum;
vector<TH2D> Hist_OnBkgd_MSCLW_Fine_Sum;
vector<TH2D> Hist_OnInit_MSCLW_Fine_Sum;

vector<double> CR_off_count_unweighted;
vector<double> CR_on_count_unweighted;
vector<double> CR_on_count_weighted;

vector<vector<double>> GammaSource_Data;
vector<vector<double>> BrightStars_Data;

bool isNaN(double x) {
  return x != x;
}
TObject* getEffAreaHistogram(string file_name, int runnumber, double offset)
{

  TFile* fAnasumDataFile = TFile::Open(file_name.c_str());

  double iSlizeY = -9999;
  //iSlizeY = offset;
  string dirname = "energyHistograms";
  string hisname = "herecEffectiveArea_on";
  if( !fAnasumDataFile )
  {
    return 0;
  }
  
  char dx[600];
  if( runnumber > 1 )
  {
    sprintf( dx, "run_%d/stereo/%s", runnumber, dirname.c_str() );
  }
  else
  {
    if( runnumber == 0 )
    {
      sprintf( dx, "total/stereo/%s", dirname.c_str() );
    }
    else if( runnumber == 1 )
    {
      sprintf( dx, "total_%d/stereo/%s", runnumber, dirname.c_str() );
    }
    else
    {
      sprintf( dx, "total_%d/stereo/%s", -1 * runnumber, dirname.c_str() );
    }
  }
  
  fAnasumDataFile->cd( dx );
  TDirectory* iDir = gDirectory;
  if( !iDir )
  {
    return 0;
  }
  
  TObject* h = ( TObject* )iDir->Get( hisname.c_str() );
  
  if( h && iSlizeY < -9998. )
  {
    return h->Clone();
  }
  else if( h )
  {
    string iClassName = h->ClassName();
    if( iClassName.find( "TH2" ) != string::npos )
    {
      TH2* i_h2 = ( TH2* )h;
      string iN = hisname + "px";
      TH1* i_h = ( TH1* )i_h2->ProjectionX( iN.c_str(), i_h2->GetYaxis()->FindBin( iSlizeY ), i_h2->GetYaxis()->FindBin( iSlizeY ) );
      return i_h->Clone();
    }
  }

  //fAnasumDataFile->Close();
  //delete fAnasumDataFile;
  
  
  return 0;
}

void GetGammaSources()
{
    std::ifstream astro_file(SMI_AUX+"/TeVCat_RaDec.txt");
    std::string line;
    // Read one line at a time into the variable line:
    while(std::getline(astro_file, line))
    {
        if (line.empty()) continue;
        std::vector<double>   lineData;
        std::stringstream  lineStream(line);
        double value;
        // Read an integer at a time from the line
        while(lineStream >> value)
        {
            // Add the integers from a line to a 1D array (vector)
            lineData.push_back(value);
        }
        if (lineData.size()!=2) continue;
        double star_ra = lineData.at(0);
        double star_dec = lineData.at(1);
        GammaSource_Data.push_back(lineData);
    }
    std::cout << "I found " << GammaSource_Data.size() << " gamma-ray sources." << std::endl;
}

void GetBrightStars()
{
    std::ifstream astro_file(SMI_AUX+"/Hipparcos_MAG8_1997.dat");
    std::string line;
    // Read one line at a time into the variable line:
    while(std::getline(astro_file, line))
    {
        if (line.find("#")!=std::string::npos) continue;
        if (line.find("*")!=std::string::npos) continue;
        if (line.empty()) continue;
        std::vector<double>   lineData;
        std::stringstream  lineStream(line);
        double value;
        // Read an integer at a time from the line
        while(lineStream >> value)
        {
            // Add the integers from a line to a 1D array (vector)
            lineData.push_back(value);
        }
        if (lineData.size()!=5) continue;
        double star_ra = lineData.at(0);
        double star_dec = lineData.at(1);
        double star_brightness = lineData.at(3)+lineData.at(4);
        if (pow((mean_tele_point_ra-star_ra)*(mean_tele_point_ra-star_ra)+(mean_tele_point_dec-star_dec)*(mean_tele_point_dec-star_dec),0.5)>Skymap_size_x) continue;
        if (star_brightness<brightness_cut)
        {
            BrightStars_Data.push_back(lineData);
        }
    }
    std::cout << "I found " << BrightStars_Data.size() << " bright stars" << std::endl;
}

bool CoincideWithGammaSources(double ra, double dec, double radius_cut)
{
    bool isCoincident = false;
    for (int star=0;star<GammaSource_Data.size();star++)
    {
        double star_ra = GammaSource_Data.at(star).at(0);
        double star_dec = GammaSource_Data.at(star).at(1);
        double radius = pow((ra-star_ra)*(ra-star_ra)+(dec-star_dec)*(dec-star_dec),0.5);
        if (radius>radius_cut) continue;
        isCoincident = true;
    }
    return isCoincident;
}

bool CoincideWithBrightStars(double ra, double dec, double evt_energy)
{
    bool isCoincident = false;
    for (int star=0;star<BrightStars_Data.size();star++)
    {
        double star_ra = BrightStars_Data.at(star).at(0);
        double star_dec = BrightStars_Data.at(star).at(1);
        double star_brightness = BrightStars_Data.at(star).at(3);
        double radius = pow((ra-star_ra)*(ra-star_ra)+(dec-star_dec)*(dec-star_dec),0.5);
        if (star_brightness>brightness_cut) continue;
        if (evt_energy<398.)
        {
            if (radius>bright_star_radius_cut) continue;
        }
        else
        {
            if (radius>bright_star_radius_cut) continue;
        }
        isCoincident = true;
    }
    return isCoincident;
}
bool FoV(double evt_ra, double evt_dec, bool isON) {

    double x = evt_ra-mean_tele_point_ra;
    double y = evt_dec-mean_tele_point_dec;
    if (source_theta_cut>(pow(x*x+y*y,0.5))) return false;
    if (!isON)
    {
        if (CoincideWithGammaSources(evt_ra,evt_dec,0.25)) return false;
    }

    return true;
}
bool FoV_Mask(double evt_ra, double evt_dec, TString target_name) {

    double x = evt_ra-mean_tele_point_ra;
    double y = evt_dec-mean_tele_point_dec;
    if ((pow(x*x+y*y,0.5))<1.0) return true;
    if (CoincideWithGammaSources(evt_ra,evt_dec,0.25)) return true;
    if (target_name.Contains("SS433"))
    {
        if (pow(pow(evt_ra-288.0833333,2)+pow(evt_dec-4.9166667,2),0.5)<0.3) return true; // SS 433 SNR
        if (pow(pow(evt_ra-288.404,2)+pow(evt_dec-4.930,2),0.5)<0.3) return true; // SS 433 e1
        if (pow(pow(evt_ra-287.654,2)+pow(evt_dec-5.037,2),0.5)<0.3) return true; // SS 433 w1
        if (pow(pow(evt_ra-287.05,2)+pow(evt_dec-6.39,2),0.5)<1.2) return true; // MGRO J1908+06
    }

    return false;
}

double GetRunPedestalVar(int run_number)
{
    string line;
    char delimiter = ' ';
    string acc_runnumber = "";
    string acc_version = "";
    string acc_nsb = "";
    int nth_line = 0;
    int nth_delimiter = 0;
    std::string::size_type sz;
    double NSB = 0.;

    ifstream myfile (SMI_AUX+"/NSB_allruns.txt");
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
                acc_nsb = "";
                for(int i = 0; i < line.size(); i++)
                {
                    if (line[i] != delimiter)
                    {
                        acc_nsb += line[i];
                    }
                    else
                    {
                        acc_nsb = "";
                        nth_delimiter += 1;
                    }
                }
                NSB = std::stod(acc_nsb,&sz);
                break;
            }
        }
        myfile.close();
    }
    else std::cout << "Unable to open file NSB_allruns.txt" << std::endl; 

    return NSB;
}

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
        //delete input_file;
        //std::cout << "root file elev = " << TelElevation_avg << " azim = " << TelAzimuth_avg << std::endl;
    }

    return std::make_pair(TelElevation_avg,TelAzimuth_avg);
}

vector<int> RemoveNonExistingRuns(vector<int>& input_list)
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
            std::cout << "Run " << run_number << " does exist." << std::endl;
            new_list.push_back(input_list.at(run));
        }
    }

    return new_list;
}

void SortingList(vector<int>* list, vector<double>* list_pointing)
{
    int temp_list;
    double temp_list_pointing;
    for (int run1=0;run1<list->size()-1;run1++)
    {
        for (int run2=run1+1;run2<list->size();run2++)
        {
            if (list_pointing->at(run1)<list_pointing->at(run2))
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
        //if (run_azim>180.) run_azim = 360.-run_azim;
        list_elev.push_back(run_elev);
        //list_elev.push_back(run_azim);
    }
    return list_elev;
}

bool PointingSelection(int int_run_number)
{

    char run_number[50];
    sprintf(run_number, "%i", int_run_number);
    string filename;
    filename = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");

    if(gSystem->AccessPathName(filename.c_str()))
    {
        std::cout << "TFile does not exist: " << filename << std::endl;
        return false;
    }

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
    Data_tree->SetBranchStatus("MVA",1);
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
    Data_tree->SetBranchAddress("MVA",&MVA);
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
        usable_time = time_end-time_start;
        input_file->Close();
        //delete input_file;
    }
    return usable_time;
}

bool SelectNImages()
{
    if (NImages<min_NImages) return false;
    if (EmissionHeight<min_EmissionHeight_cut) return false;
    if (EmissionHeight>max_EmissionHeight_cut) return false;
    if (EChi2S>max_Eerr) return false;
    if (Xcore*Xcore+Ycore*Ycore>max_Rcore*max_Rcore) return false;
    if (Xcore*Xcore+Ycore*Ycore<min_Rcore*min_Rcore) return false;
    if (pow(R2off,0.5)>max_Roff) return false;
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
    //delete input_file;

    return std::make_pair(TelRAJ2000_tmp,TelDecJ2000_tmp);
}

double CountCosmicRayEvents(int int_run_number, int input_xoff_idx, int input_yoff_idx)
{

    char run_number[50];
    sprintf(run_number, "%i", int_run_number);
    string filename;
    filename = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");

    if(gSystem->AccessPathName(filename.c_str()))
    {
        std::cout << "TFile does not exist: " << filename << std::endl;
        return 0.;
    }

    TFile*  input_file = TFile::Open(filename.c_str());
    TString root_file = "run_"+string(run_number)+"/stereo/DL3EventTree";
    TTree* Data_tree = (TTree*) input_file->Get(root_file);
    if (!Data_tree)
    {
        std::cout << "TTree does not exist: " << root_file << std::endl;
        return 0.;
    }
    SetEventDisplayDL3TreeBranch(Data_tree);
    vector<pair<double,double>> timecut_thisrun = GetRunTimecuts(int_run_number);

    Data_tree->GetEntry(0);
    double time_0 = Time;
    Data_tree->GetEntry(Data_tree->GetEntries()-1);
    double time_1 = Time;

    double MSCW_bin_size = (MSCW_upper_blind-MSCW_lower_blind)/double(mtx_dim_w_fine);
    double MSCL_bin_size = (MSCL_upper_blind-MSCL_lower_blind)/double(mtx_dim_l_fine);
    double MSCW_plot_upper = MSCW_upper_blind+n_extra_upper_bins*MSCW_bin_size;
    double MSCL_plot_upper = MSCL_upper_blind+n_extra_upper_bins*MSCL_bin_size;
    double MSCW_plot_lower = MSCW_lower_blind-n_extra_lower_bins*MSCW_bin_size;
    double MSCL_plot_lower = MSCL_lower_blind-n_extra_lower_bins*MSCL_bin_size;

    TH1D Hist_ErecS = TH1D("Hist_ErecS","",N_energy_bins,energy_bins);
    TH1D Hist_Xoff = TH1D("Hist_Xoff","",N_Xoff_bins,-max_Roff,max_Roff);
    TH1D Hist_Yoff = TH1D("Hist_Yoff","",N_Yoff_bins,-max_Roff,max_Roff);
    double cosmic_ray_events = 0.;
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
        MVA = 0;
        R2off = 0;
        Data_tree->GetEntry(entry);
        R2off = Xoff*Xoff+Yoff*Yoff;
        Phioff = atan2(Yoff,Xoff)+M_PI;

        int energy_idx = Hist_ErecS.FindBin(ErecS*1000.)-1;
        if (!SelectNImages()) continue;
        if (!ApplyTimeCuts(Time-time_0, timecut_thisrun)) continue;
        if (energy_idx<0) continue;
        if (energy_idx>=N_energy_bins) continue;
        int idx_xoff = Hist_Xoff.FindBin(Xoff)-1;
        int idx_yoff = Hist_Yoff.FindBin(Yoff)-1;
        if (input_xoff_idx!=idx_xoff) continue;
        if (input_yoff_idx!=idx_yoff) continue;


        if (MSCL>MSCL_plot_upper) continue;
        if (MSCL<MSCL_plot_lower) continue;
        if (MSCW>MSCW_plot_upper) continue;
        if (MSCW<MSCW_plot_lower) continue;

        if (MSCL<MSCL_upper_blind && MSCW<MSCW_upper_blind && MSCL>MSCL_lower_blind && MSCW>MSCW_lower_blind)
        {
            continue;
        }
        else
        {
            cosmic_ray_events += 1.;
        }

    }
    input_file->Close();

    return cosmic_ray_events;

}

void TrainingRunAnalysis(int int_run_number, int int_run_number_on, int input_xoff_idx, int input_yoff_idx, double evt_weight)
{

    char run_number[50];
    sprintf(run_number, "%i", int_run_number);
    string filename;
    filename = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");

    if(gSystem->AccessPathName(filename.c_str()))
    {
        std::cout << "TFile does not exist: " << filename << std::endl;
        return;
    }

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
    double off_tele_elev = GetRunElevAzim(int_run_number).first;
    double off_tele_azim = GetRunElevAzim(int_run_number).second;
    if (off_tele_azim>180.) off_tele_azim = 360.-off_tele_azim;
    double off_NSB = GetRunPedestalVar(int_run_number);
    double on_tele_elev = GetRunElevAzim(int_run_number_on).first;
    double on_tele_azim = GetRunElevAzim(int_run_number_on).second;
    if (on_tele_azim>180.) on_tele_azim = 360.-on_tele_azim;
    double on_NSB = GetRunPedestalVar(int_run_number_on);

    //std::cout << "Get telescope pointing RA and Dec for run " << run_number << std::endl;
    std::pair<double,double> tele_point_ra_dec = GetRunRaDec(filename,int_run_number);

    char run_number_on[50];
    sprintf(run_number_on, "%i", int_run_number_on);
    string filename_on;
    filename_on = TString(SMI_INPUT+"/"+string(run_number_on)+".anasum.root");
    std::pair<double,double> on_tele_point_ra_dec = GetRunRaDec(filename_on,int_run_number_on);

    Data_tree->GetEntry(0);
    double time_0 = Time;
    Data_tree->GetEntry(Data_tree->GetEntries()-1);
    double time_1 = Time;
    exposure_hours += (time_1-time_0)/3600.;
    //std::cout << "Get usable time for run " << run_number << std::endl;
    double exposure_thisrun = GetRunUsableTime(filename,int_run_number)/3600.;

    double MSCW_bin_size = (MSCW_upper_blind-MSCW_lower_blind)/double(mtx_dim_w_fine);
    double MSCL_bin_size = (MSCL_upper_blind-MSCL_lower_blind)/double(mtx_dim_l_fine);
    double MSCW_plot_upper = MSCW_upper_blind+n_extra_upper_bins*MSCW_bin_size;
    double MSCL_plot_upper = MSCL_upper_blind+n_extra_upper_bins*MSCL_bin_size;
    double MSCW_plot_lower = MSCW_lower_blind-n_extra_lower_bins*MSCW_bin_size;
    double MSCL_plot_lower = MSCL_lower_blind-n_extra_lower_bins*MSCL_bin_size;

    TH1D Hist_ErecS = TH1D("Hist_ErecS","",N_energy_bins,energy_bins);
    TH1D Hist_Xoff = TH1D("Hist_Xoff","",N_Xoff_bins,-max_Roff,max_Roff);
    TH1D Hist_Yoff = TH1D("Hist_Yoff","",N_Yoff_bins,-max_Roff,max_Roff);
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
        MVA = 0;
        R2off = 0;
        Data_tree->GetEntry(entry);
        R2off = Xoff*Xoff+Yoff*Yoff;
        Phioff = atan2(Yoff,Xoff)+M_PI;
        ra_sky = tele_point_ra_dec.first+Xoff_derot;
        dec_sky = tele_point_ra_dec.second+Yoff_derot;
        double on_ra_sky = on_tele_point_ra_dec.first+Xoff_derot;
        double on_dec_sky = on_tele_point_ra_dec.second+Yoff_derot;
        theta2 = pow(ra_sky-mean_tele_point_ra,2)+pow(dec_sky-mean_tele_point_dec,2);
        int energy_idx = Hist_ErecS.FindBin(ErecS*1000.)-1;
        if (!SelectNImages()) continue;
        if (!ApplyTimeCuts(Time-time_0, timecut_thisrun)) continue;
        if (energy_idx<0) continue;
        if (energy_idx>=N_energy_bins) continue;
        int idx_xoff = Hist_Xoff.FindBin(Xoff)-1;
        int idx_yoff = Hist_Yoff.FindBin(Yoff)-1;
        if (input_xoff_idx!=idx_xoff) continue;
        if (input_yoff_idx!=idx_yoff) continue;

        if (!FoV(ra_sky, dec_sky, false)) continue;

        //MSCL = -1.*2.5*(MVA-0.25);

        Hist_OffData_MSCLW.at(energy_idx).Fill(MSCL,MSCW,evt_weight);
        Hist_OffData_MSCLW_Fine.at(energy_idx).Fill(MSCL,MSCW,evt_weight);

        if (!(MSCL<MSCL_upper_blind && MSCW<MSCW_upper_blind && MSCL>MSCL_lower_blind && MSCW>MSCW_lower_blind))
        {
            CR_off_count_unweighted.at(energy_idx) += 1.;
        }

        if (MSCL>MSCL_plot_upper) continue;
        if (MSCL<MSCL_plot_lower) continue;
        if (MSCW>MSCW_plot_upper) continue;
        if (MSCW<MSCW_plot_lower) continue;

        if (MSCL<MSCL_upper_blind && MSCW<MSCW_upper_blind && MSCL>MSCL_lower_blind && MSCW>MSCW_lower_blind)
        {
            Hist_OffData_SR_XYoff.at(energy_idx).Fill(Xoff,Yoff,evt_weight);
            //Hist_OffData_SR_XYoff.at(energy_idx).Fill(on_ra_sky,on_dec_sky);
        }
        else
        {
            Hist_OffData_CR_XYoff.at(energy_idx).Fill(Xoff,Yoff,evt_weight);
            //Hist_OffData_CR_XYoff.at(energy_idx).Fill(on_ra_sky,on_dec_sky);
            off_cr_count += 1.*evt_weight;
            avg_diff_nsb += (off_NSB-on_NSB)*evt_weight;
            avg_diff_el += (off_tele_elev-on_tele_elev)*evt_weight;
            avg_diff_az += (off_tele_azim-on_tele_azim)*evt_weight;
        }

    }
    input_file->Close();
    //delete input_file;

}

double default_effective_area(double input_log_energy_tev)
{
    double log_energy_tev[23] = {-0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45};
    double effective_area[22] = {15000., 15000., 15000., 20793.1, 48022, 62902.8, 75192.9, 79701.3, 84345, 86425.2, 86141.3, 81266.6, 82327.6, 81233.2, 78587.7, 74493.6, 72086.3, 70820.1, 70299.2, 67367.5, 67767.3, 64167.6};
    TH1D Hist_EffArea = TH1D("Hist_EffArea","",22,log_energy_tev);
    for (int bin=0;bin<Hist_EffArea.GetNbinsX();bin++)
    {
        Hist_EffArea.SetBinContent(bin,effective_area[bin]);
    }
    return Hist_EffArea.GetBinContent( Hist_EffArea.FindBin(input_log_energy_tev) );
}

void SingleRunAnalysis(int int_run_number, int int_run_number_real, int input_xoff_idx, int input_yoff_idx, bool isImposter)
{

    char run_number[50];
    sprintf(run_number, "%i", int_run_number);
    string filename;
    filename = TString(SMI_INPUT+"/"+string(run_number)+".anasum.root");
    
    if(gSystem->AccessPathName(filename.c_str()))
    {
        std::cout << "TFile does not exist: " << filename << std::endl;
        return;
    }

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
    double tele_azim = GetRunElevAzim(int_run_number).second;
    if (tele_azim>180.) tele_azim = 360.-tele_azim;
    double NSB = GetRunPedestalVar(int_run_number);

    NSB_mean += NSB;
    Elev_mean += tele_elev;
    Azim_mean += tele_azim;

    //std::cout << "Get telescope pointing RA and Dec for run " << run_number << std::endl;
    std::pair<double,double> tele_point_ra_dec = GetRunRaDec(filename,int_run_number);

    char run_number_real[50];
    sprintf(run_number_real, "%i", int_run_number_real);
    string filename_real;
    filename_real = TString(SMI_INPUT+"/"+string(run_number_real)+".anasum.root");
    std::pair<double,double> real_tele_point_ra_dec = GetRunRaDec(filename_real,int_run_number_real);

    //TH1* i_hEffAreaP = ( TH1* )getEffAreaHistogram(filename, int_run_number, 0.5);
    //for (int e=0;e<N_energy_bins;e++) 
    //{
    //    //double eff_area = i_hEffAreaP->GetBinContent( i_hEffAreaP->FindBin( log10(0.5*(energy_bins[e]+energy_bins[e+1])/1000.)));
    //    double eff_area = default_effective_area( log10(0.5*(energy_bins[e]+energy_bins[e+1])/1000.) );
    //    effective_area.at(e) += eff_area;
    //}

    Data_tree->GetEntry(0);
    double time_0 = Time;
    Data_tree->GetEntry(Data_tree->GetEntries()-1);
    double time_1 = Time;
    exposure_hours += (time_1-time_0)/3600.;
    //std::cout << "Get usable time for run " << run_number << std::endl;
    double exposure_thisrun = GetRunUsableTime(filename,int_run_number)/3600.;
    exposure_hours_usable += exposure_thisrun;

    double MSCW_bin_size = (MSCW_upper_blind-MSCW_lower_blind)/double(mtx_dim_w_fine);
    double MSCL_bin_size = (MSCL_upper_blind-MSCL_lower_blind)/double(mtx_dim_l_fine);
    double MSCW_plot_upper = MSCW_upper_blind+n_extra_upper_bins*MSCW_bin_size;
    double MSCL_plot_upper = MSCL_upper_blind+n_extra_upper_bins*MSCL_bin_size;
    double MSCW_plot_lower = MSCW_lower_blind-n_extra_lower_bins*MSCW_bin_size;
    double MSCL_plot_lower = MSCL_lower_blind-n_extra_lower_bins*MSCL_bin_size;

    TH1D Hist_ErecS = TH1D("Hist_ErecS","",N_energy_bins,energy_bins);
    TH1D Hist_Xoff = TH1D("Hist_Xoff","",N_Xoff_bins,-max_Roff,max_Roff);
    TH1D Hist_Yoff = TH1D("Hist_Yoff","",N_Yoff_bins,-max_Roff,max_Roff);
    TH1D Hist_Expo_Roff = TH1D("Hist_Expo_Roff","",4,0.,2.);
    TH2D Hist_SingleRun_AreaTime_Skymap = TH2D("Hist_SingleRun_AreaTime_Skymap","",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_x,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y);
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
        MVA = 0.;
        R2off = 0;
        Data_tree->GetEntry(entry);
        R2off = Xoff*Xoff+Yoff*Yoff;
        Phioff = atan2(Yoff,Xoff)+M_PI;
        ra_sky = real_tele_point_ra_dec.first+Xoff_derot;
        dec_sky = real_tele_point_ra_dec.second+Yoff_derot;
        complementary_ra_sky = real_tele_point_ra_dec.first-Xoff_derot;
        complementary_dec_sky = real_tele_point_ra_dec.second-Yoff_derot;
        theta2 = pow(ra_sky-mean_tele_point_ra,2)+pow(dec_sky-mean_tele_point_dec,2);
        int energy_idx = Hist_ErecS.FindBin(ErecS*1000.)-1;
        if (!SelectNImages()) continue;
        if (!ApplyTimeCuts(Time-time_0, timecut_thisrun)) continue;
        if (energy_idx<0) continue;
        if (energy_idx>=N_energy_bins) continue;
        int idx_xoff = Hist_Xoff.FindBin(Xoff)-1;
        int idx_yoff = Hist_Yoff.FindBin(Yoff)-1;
        if (input_xoff_idx!=idx_xoff) continue;
        if (input_yoff_idx!=idx_yoff) continue;

        //MSCL = -1.*2.5*(MVA-0.25);

        double map_right_edge = map_center_x+Skymap_size_x;
        double map_left_edge = map_center_x-Skymap_size_x;
        double map_top_edge = map_center_y+Skymap_size_y;
        double map_bottom_edge = map_center_y-Skymap_size_y;
        //if (ra_sky>map_right_edge) continue;
        //if (ra_sky<map_left_edge) continue;
        //if (dec_sky>map_top_edge) continue;
        //if (dec_sky<map_bottom_edge) continue;

        if (!FoV(ra_sky, dec_sky, true)) continue;
        if (isImposter)
        {
            double imposter_ra_sky = tele_point_ra_dec.first+Xoff_derot;
            double imposter_dec_sky = tele_point_ra_dec.second+Yoff_derot;
            if (!FoV(imposter_ra_sky, imposter_dec_sky, false)) continue;
        }

        Hist_OnData_MSCLW.at(energy_idx).Fill(MSCL,MSCW);
        Hist_OnData_MSCLW_Fine.at(energy_idx).Fill(MSCL,MSCW);

        if (!(MSCL<MSCL_upper_blind && MSCW<MSCW_upper_blind && MSCL>MSCL_lower_blind && MSCW>MSCW_lower_blind))
        {
            CR_on_count_unweighted.at(energy_idx) += 1.;
            Hist_OnData_Expo_Skymap_Sum.at(energy_idx).Fill(ra_sky,dec_sky);
        }

        if (MSCL>MSCL_plot_upper) continue;
        if (MSCL<MSCL_plot_lower) continue;
        if (MSCW>MSCW_plot_upper) continue;
        if (MSCW<MSCW_plot_lower) continue;

        //double evt_eff_area = i_hEffAreaP->GetBinContent( i_hEffAreaP->FindBin( log10(ErecS)));
        double evt_eff_area = 0.;

        if (MSCL<MSCL_upper_blind && MSCW<MSCW_upper_blind && MSCL>MSCL_lower_blind && MSCW>MSCW_lower_blind)
        {
            Hist_OnData_SR_XYoff.at(energy_idx).Fill(Xoff,Yoff);
            Hist_OnData_SR_Skymap.at(energy_idx).Fill(ra_sky,dec_sky);
            if (!FoV_Mask(ra_sky, dec_sky, TString(target)))
            {
                Hist_OnData_SR_Skymap_Mask.at(energy_idx).Fill(ra_sky,dec_sky);
            }
        }
        else
        {
            Hist_OnData_CR_XYoff.at(energy_idx).Fill(Xoff,Yoff);
            int local_xoff_idx = Hist_OffData_Ratio_XYoff.at(energy_idx).GetXaxis()->FindBin(Xoff);
            int local_yoff_idx = Hist_OffData_Ratio_XYoff.at(energy_idx).GetYaxis()->FindBin(Yoff);
            //int local_xoff_idx = Hist_OffData_Ratio_XYoff.at(energy_idx).GetXaxis()->FindBin(ra_sky);
            //int local_yoff_idx = Hist_OffData_Ratio_XYoff.at(energy_idx).GetYaxis()->FindBin(dec_sky);
            double weight = Hist_OffData_Ratio_XYoff.at(energy_idx).GetBinContent(local_xoff_idx,local_yoff_idx);

            if (CoincideWithBrightStars(complementary_ra_sky,complementary_dec_sky,ErecS*1000.))
            {
                Hist_OnData_CR_Skymap_FoV.at(energy_idx).Fill(complementary_ra_sky,complementary_dec_sky,weight*0.5);
                Hist_OnData_CR_Skymap_Ratio.at(energy_idx).Fill(complementary_ra_sky,complementary_dec_sky,weight*0.5);
                Hist_OnData_CR_Skymap_Regression.at(energy_idx).Fill(complementary_ra_sky,complementary_dec_sky,weight*0.5);
                Hist_OnData_CR_Skymap_Init_Perturbation.at(energy_idx).Fill(complementary_ra_sky,complementary_dec_sky,weight*0.5);
                Hist_OnData_CR_Skymap_Perturbation.at(energy_idx).Fill(complementary_ra_sky,complementary_dec_sky,weight*0.5);
                if (!FoV_Mask(ra_sky, dec_sky, TString(target)))
                {
                    Hist_OnData_CR_Skymap_Mask.at(energy_idx).Fill(complementary_ra_sky,complementary_dec_sky,weight*0.5);
                }
            }
            if (CoincideWithBrightStars(ra_sky,dec_sky,ErecS*1000.))
            {
                Hist_OnData_CR_Skymap_FoV.at(energy_idx).Fill(ra_sky,dec_sky,weight*0.5);
                Hist_OnData_CR_Skymap_Ratio.at(energy_idx).Fill(ra_sky,dec_sky,weight*0.5);
                Hist_OnData_CR_Skymap_Regression.at(energy_idx).Fill(ra_sky,dec_sky,weight*0.5);
                Hist_OnData_CR_Skymap_Init_Perturbation.at(energy_idx).Fill(ra_sky,dec_sky,weight*0.5);
                Hist_OnData_CR_Skymap_Perturbation.at(energy_idx).Fill(ra_sky,dec_sky,weight*0.5);
                if (!FoV_Mask(ra_sky, dec_sky, TString(target)))
                {
                    Hist_OnData_CR_Skymap_Mask.at(energy_idx).Fill(ra_sky,dec_sky,weight*0.5);
                }
            }
            else
            {
                Hist_OnData_CR_Skymap_FoV.at(energy_idx).Fill(ra_sky,dec_sky,weight);
                Hist_OnData_CR_Skymap_Ratio.at(energy_idx).Fill(ra_sky,dec_sky,weight);
                Hist_OnData_CR_Skymap_Regression.at(energy_idx).Fill(ra_sky,dec_sky,weight);
                Hist_OnData_CR_Skymap_Init_Perturbation.at(energy_idx).Fill(ra_sky,dec_sky,weight);
                Hist_OnData_CR_Skymap_Perturbation.at(energy_idx).Fill(ra_sky,dec_sky,weight);
                if (!FoV_Mask(ra_sky, dec_sky, TString(target)))
                {
                    Hist_OnData_CR_Skymap_Mask.at(energy_idx).Fill(ra_sky,dec_sky,weight);
                }
            }
            Hist_Expo_Roff.Fill(pow(R2off,0.5),weight);
            //Hist_SingleRun_AreaTime_Skymap.Fill(ra_sky,dec_sky,evt_eff_area*exposure_thisrun*weight);
            Hist_SingleRun_AreaTime_Skymap.Fill(ra_sky,dec_sky,weight);
            CR_on_count_weighted.at(energy_idx) += weight;

            Hist_Data_Norm_Skymap.Fill(ra_sky,dec_sky);
            Hist_Data_Elev_Skymap.Fill(ra_sky,dec_sky,tele_elev);
            Hist_Data_Azim_Skymap.Fill(ra_sky,dec_sky,tele_azim);
            Hist_Data_NSB_Skymap.Fill(ra_sky,dec_sky,NSB);
            total_cr_count += 1.;
        }

    }
    double expo_scaling = 0.;
    if (Hist_Expo_Roff.GetBinContent(1)>0.)
    {
        expo_scaling = exposure_thisrun/Hist_Expo_Roff.GetBinContent(1);
    }
    Hist_Data_AreaTime_Skymap.Add(&Hist_SingleRun_AreaTime_Skymap,expo_scaling);
    input_file->Close();
    //delete input_file;

}

pair<MatrixXcd,VectorXcd> RemoveEmptyRows(MatrixXcd mtx_input,VectorXcd vtr_input)
{
    int n_empty_rows = 0;
    vector<bool> empty_rows;
    for (int row=0;row<mtx_input.rows();row++)
    {
        bool row_is_empty = true;
        for (int col=0;col<mtx_input.cols();col++)
        {
            if (mtx_input(row,col)!=0.)
            {
                row_is_empty = false;
            }
        }
        empty_rows.push_back(row_is_empty);
        if (row_is_empty)
        {
            n_empty_rows += 1;
        }
    }
    //std::cout << "n_empty_rows = " << n_empty_rows << std::endl;
    MatrixXcd mtx_output = MatrixXcd::Zero(mtx_input.rows()-n_empty_rows,mtx_input.cols());
    VectorXcd vtr_output = VectorXcd::Zero(mtx_input.rows()-n_empty_rows);
    int new_row = -1;
    for (int row=0;row<mtx_input.rows();row++)
    {
        if (empty_rows.at(row)) continue;
        new_row += 1;
        vtr_output(new_row) = vtr_input(row);
        for (int col=0;col<mtx_input.cols();col++)
        {
            mtx_output(new_row,col) = mtx_input(row,col);
        }
    }
    return std::make_pair(mtx_output,vtr_output);
}
VectorXcd SolutionWithConstraints(MatrixXcd mtx_big, MatrixXcd mtx_constraints_input, VectorXcd vtr_delta, VectorXcd vtr_constraints_delta_input)
{

    MatrixXcd mtx_constraints = RemoveEmptyRows(mtx_constraints_input,vtr_constraints_delta_input).first;
    VectorXcd vtr_constraints_delta = RemoveEmptyRows(mtx_constraints_input,vtr_constraints_delta_input).second;
    if (mtx_constraints.rows()==0)
    {
        BDCSVD<MatrixXcd> svd(mtx_big, ComputeThinU | ComputeThinV);
        VectorXcd vtr_vari_big = VectorXcd::Zero(mtx_big.cols());
        vtr_vari_big = svd.solve(vtr_delta);
        return vtr_vari_big;
    }

    MatrixXcd BTB = mtx_big.transpose()*mtx_big;
    VectorXcd BTD = mtx_big.transpose()*vtr_delta;

    MatrixXcd mtx_Bigger = MatrixXcd::Zero(BTB.rows()+mtx_constraints.rows(),BTB.cols()+mtx_constraints.rows());
    mtx_Bigger.block(0,0,BTB.rows(),BTB.cols()) = 2.*BTB;
    mtx_Bigger.block(BTB.rows(),0,mtx_constraints.rows(),mtx_constraints.cols()) = mtx_constraints;
    mtx_Bigger.block(0,BTB.cols(),mtx_constraints.cols(),mtx_constraints.rows()) = mtx_constraints.transpose();

    VectorXcd vtr_bigger_delta = VectorXcd::Zero(BTB.rows()+mtx_constraints.rows());
    vtr_bigger_delta.segment(0,BTB.cols()) = 2.*BTD;
    vtr_bigger_delta.segment(BTB.cols(),vtr_constraints_delta.size()) = vtr_constraints_delta;

    ComplexEigenSolver<MatrixXcd> eigensolver_bigger = ComplexEigenSolver<MatrixXcd>(mtx_Bigger);

    VectorXcd vtr_vari_bigger = VectorXcd::Zero(BTB.cols()+mtx_constraints.rows());
    BDCSVD<MatrixXcd> svd(mtx_Bigger, ComputeThinU | ComputeThinV);
    vtr_vari_bigger = svd.solve(vtr_bigger_delta);

    return vtr_vari_bigger.segment(0,BTB.cols());

}


MatrixXcd MatrixPerturbationMethod(int e_idx, MatrixXcd mtx_t_input, MatrixXcd mtx_base_input, MatrixXcd mtx_init_input, MatrixXcd mtx_data_input, int max_rank, int var_rank, bool diagonal_rank, bool use_mtx_t_input, double log_t_input_weight, bool isBlind, int return_type)
{

    MatrixXcd mtx_output = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    if (mtx_init_input.cwiseAbs().sum()<1000. || mtx_data_input.cwiseAbs().sum()<1000. || mtx_base_input.cwiseAbs().sum()<1000.)
    {
        return mtx_output;
    }
    //if (max_rank==1)
    //{
    //    return mtx_output;
    //}

    JacobiSVD<MatrixXd> svd_data(mtx_data_input.real(), ComputeFullU | ComputeFullV);
    MatrixXd mtx_U_data = svd_data.matrixU();
    MatrixXd mtx_V_data = svd_data.matrixV();
    MatrixXd mtx_S_data = MatrixXd::Zero(mtx_data_input.rows(),mtx_data_input.cols());
    for (int entry=0;entry<svd_data.singularValues().size();entry++)
    {
        mtx_S_data(entry,entry) = svd_data.singularValues()(entry);
    }

    JacobiSVD<MatrixXd> svd_init(mtx_init_input.real(), ComputeFullU | ComputeFullV);
    MatrixXd mtx_U_init = svd_init.matrixU();
    MatrixXd mtx_V_init = svd_init.matrixV();
    MatrixXd mtx_S_init = MatrixXd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    for (int entry=0;entry<svd_init.singularValues().size();entry++)
    {
        mtx_S_init(entry,entry) = svd_init.singularValues()(entry);
    }

    JacobiSVD<MatrixXd> svd_base(mtx_base_input.real(), ComputeFullU | ComputeFullV);
    MatrixXd mtx_U_base = svd_base.matrixU();
    MatrixXd mtx_V_base = svd_base.matrixV();
    MatrixXd mtx_S_base = MatrixXd::Zero(mtx_base_input.rows(),mtx_base_input.cols());
    for (int entry=0;entry<svd_base.singularValues().size();entry++)
    {
        mtx_S_base(entry,entry) = svd_base.singularValues()(entry);
    }

    //entry_size = 1;
    //for (int entry=0;entry<svd_init.singularValues().size()-1;entry++)
    //{
    //    entry_size = entry+1;
    //    if (mtx_S_init(entry,entry)/mtx_S_init(entry+1,entry+1)<5.)
    //    {
    //        break;
    //    }
    //}
    //std::cout << "entry_size = " << entry_size << std::endl;

    int size_k = mtx_init_input.cols();
    int size_n = mtx_init_input.cols();
    int length_tkn = size_k*size_n;
    int regularization_size = 0;
    regularization_size = length_tkn;
    //VectorXcd vtr_Delta = VectorXcd::Zero(mtx_init_input.rows()*mtx_init_input.cols());
    //MatrixXcd mtx_A = MatrixXcd::Zero(mtx_init_input.rows()*mtx_init_input.cols(),length_tkn);
    //MatrixXcd mtx_W = MatrixXcd::Zero(mtx_init_input.rows()*mtx_init_input.cols(),mtx_init_input.rows()*mtx_init_input.cols());
    VectorXcd vtr_Delta = VectorXcd::Zero(mtx_init_input.rows()*mtx_init_input.cols()+regularization_size);
    MatrixXcd mtx_A = MatrixXcd::Zero(mtx_init_input.rows()*mtx_init_input.cols()+regularization_size,length_tkn);
    MatrixXcd mtx_W = MatrixXcd::Zero(mtx_init_input.rows()*mtx_init_input.cols()+regularization_size,mtx_init_input.rows()*mtx_init_input.cols()+regularization_size);
    MatrixXcd mtx_Constraint = MatrixXcd::Zero(regularization_size,length_tkn);
    VectorXcd vtr_Constraint_Delta = VectorXcd::Zero(regularization_size);

    double available_bins = 0.;
    double avg_weight = 0.;
    for (int idx_i=0;idx_i<mtx_init_input.rows();idx_i++)
    {
        for (int idx_j=0;idx_j<mtx_init_input.rows();idx_j++)
        {
            int idx_u = idx_j*mtx_init_input.rows() + idx_i;
            double sigma_data = max(1.,pow(mtx_init_input(idx_i,idx_j).real(),0.5));
            double stat_weight = 1./pow(sigma_data*sigma_data,0.5);
            double weight = 1.;
            if (use_stat_err_weight) 
            {
                weight = stat_weight;
            }
            mtx_W(idx_u,idx_u) = weight;
            vtr_Delta(idx_u) = (mtx_data_input-mtx_init_input)(idx_i,idx_j);
            if (isBlind)
            {
                if (idx_i<binx_blind_upper_global && idx_j<biny_blind_upper_global && idx_i>=binx_blind_lower_global && idx_j>=biny_blind_lower_global)
                {
                    mtx_W(idx_u,idx_u) = 0.;
                }
            }
            available_bins += 1.;
            avg_weight += weight*weight;
            for (int idx_k=0;idx_k<size_k;idx_k++)
            {
                int kth_entry = idx_k+1;
                for (int idx_n=0;idx_n<size_n;idx_n++)
                {
                    int nth_entry = idx_n+1;
                    if (kth_entry==nth_entry && !diagonal_rank) continue;
                    if (kth_entry==nth_entry && kth_entry>var_rank) continue;
                    if (kth_entry>max_rank && nth_entry>max_rank) continue;
                    if (kth_entry>max_rank) continue;
                    if (nth_entry>max_rank) continue;
                    //if (kth_entry+nth_entry>max_rank+1) continue;
                    int idx_v = idx_k*size_n + idx_n;
                    mtx_A(idx_u,idx_v) = mtx_U_base(idx_i,idx_k)*mtx_V_base(idx_j,idx_n);
                }
            }
        }
    }
    avg_weight = pow(avg_weight/available_bins,0.5);

    
    //if (isBlind && diagonal_rank)
    //{
    //    int idx_k1 = 1;
    //    int idx_n1 = 1;
    //    int idx_k2 = 0;
    //    int idx_n2 = 0;
    //    int idx_k3 = 0;
    //    int idx_n3 = 1;
    //    int idx_k4 = 1;
    //    int idx_n4 = 0;
    //    int idx_v1 = idx_k1*size_n + idx_n1;
    //    int idx_v2 = idx_k2*size_n + idx_n2;
    //    int idx_v3 = idx_k3*size_n + idx_n3;
    //    int idx_v4 = idx_k4*size_n + idx_n4;
    //    int idx_u1 = idx_v1;
    //    double x00 = coefficient_x00[e_idx]; 
    //    double x01 = coefficient_x01[e_idx]; 
    //    double x10 = coefficient_x10[e_idx]; 
    //    mtx_Constraint(idx_u1,idx_v1) = 1.;
    //    mtx_Constraint(idx_u1,idx_v2) = -x00;
    //    mtx_Constraint(idx_u1,idx_v3) = -x01;
    //    mtx_Constraint(idx_u1,idx_v4) = -x10;
    //}

    if (use_mtx_t_input)
    {
        for (int idx_k=0;idx_k<size_k;idx_k++)
        {
            int kth_entry = idx_k+1;
            for (int idx_n=0;idx_n<size_n;idx_n++)
            {
                int nth_entry = idx_n+1;
                int idx_v = idx_k*size_n + idx_n;
                int idx_u = idx_v + mtx_init_input.rows()*mtx_init_input.cols();
                if (abs(mtx_t_input(idx_k,idx_n))<1.)
                {
                    if (kth_entry==1 && nth_entry==1)
                    {
                        mtx_W(idx_u,idx_u) = pow(10.,-3.)*avg_weight;
                        mtx_A(idx_u,idx_v) = 1.;
                        vtr_Delta(idx_u) = 0.;
                    }
                    if (kth_entry==2 && nth_entry==2 && var_rank>1)
                    {
                        mtx_W(idx_u,idx_u) = pow(10.,log_t_input_weight)*avg_weight;
                        mtx_A(idx_u,idx_v) = 1.;
                        double x00 = coefficient_x00[e_idx]; 
                        vtr_Delta(idx_u)  = x00*mtx_t_input(0,0);
                        double x01 = coefficient_x01[e_idx]; 
                        vtr_Delta(idx_u) += x01*mtx_t_input(0,1);
                        double x10 = coefficient_x10[e_idx]; 
                        vtr_Delta(idx_u) += x10*mtx_t_input(1,0);
                    }
                }
                else
                {
                    mtx_W(idx_u,idx_u) = pow(10.,log_t_input_weight)*avg_weight;
                    mtx_A(idx_u,idx_v) = 1.;
                    vtr_Delta(idx_u) = mtx_t_input(idx_k,idx_n);
                }
            }
        }
    }

    VectorXcd vtr_t = VectorXcd::Zero(length_tkn);

    VectorXcd vtr_Delta_weighted = mtx_A.transpose()*mtx_W*vtr_Delta;
    MatrixXcd mtx_A_weighted = mtx_A.transpose()*mtx_W*mtx_A; 
    BDCSVD<MatrixXcd> bdc_svd(mtx_A_weighted, ComputeThinU | ComputeThinV);


    if (return_type==0)
    {
        for (int idx_k=0;idx_k<size_k;idx_k++)
        {
            int kth_entry = idx_k+1;
            for (int idx_n=0;idx_n<size_n;idx_n++)
            {
                int nth_entry = idx_n+1;
                int idx_kn = idx_k*size_n + idx_n;
                vtr_t(idx_kn) = mtx_t_input(idx_k,idx_n);
            }
        }
    }
    else
    {
        vtr_t = SolutionWithConstraints(mtx_A_weighted, mtx_Constraint, vtr_Delta_weighted, vtr_Constraint_Delta);
        //vtr_t = bdc_svd.solve(vtr_Delta_weighted);
    }


    MatrixXcd mtx_t = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_C = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_D = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_E = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_CDE = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    for (int idx_k=0;idx_k<size_k;idx_k++)
    {
        int kth_entry = idx_k+1;
        std::complex<double> sigma_k = mtx_S_base(idx_k,idx_k);
        for (int idx_n=0;idx_n<size_n;idx_n++)
        {
            int nth_entry = idx_n+1;
            std::complex<double> sigma_n = mtx_S_base(idx_n,idx_n);
            int idx_kn = idx_k*size_n + idx_n;
            int idx_nk = idx_n*size_k + idx_k;
            mtx_t(idx_k,idx_n) = vtr_t(idx_kn); 
            if (kth_entry!=nth_entry)
            {
                mtx_C(idx_k,idx_n) = (vtr_t(idx_kn)/sigma_k+vtr_t(idx_nk)/sigma_n)/(sigma_n/sigma_k-sigma_k/sigma_n); 
                mtx_D(idx_k,idx_n) = (vtr_t(idx_kn)/sigma_n+vtr_t(idx_nk)/sigma_k)/(sigma_n/sigma_k-sigma_k/sigma_n); 
                if (idx_k>idx_n)
                {
                    mtx_CDE(idx_k,idx_n) = mtx_C(idx_k,idx_n);
                    mtx_CDE(idx_n,idx_k) = mtx_D(idx_k,idx_n);
                }
            }
            else
            {
                mtx_E(idx_k,idx_n) = vtr_t(idx_kn); 
                if (kth_entry<=max_rank && nth_entry<=max_rank)
                {
                    mtx_CDE(idx_k,idx_n) = mtx_E(idx_k,idx_n)/abs(mtx_S_base(idx_k,idx_n));
                }
            }
        }
    }

    MatrixXcd mtx_U_vari = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_V_vari = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    MatrixXcd mtx_S_vari = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    for (int idx_i=0;idx_i<mtx_init_input.cols();idx_i++)
    {
        for (int idx_k=0;idx_k<size_k;idx_k++)
        {
            int kth_entry = idx_k+1;
            for (int idx_n=0;idx_n<size_n;idx_n++)
            {
                int nth_entry = idx_n+1;
                mtx_S_vari(idx_k,idx_k) = mtx_E(idx_k,idx_k);
                std::complex<double> coeff_u = mtx_C(idx_k,idx_n);
                mtx_U_vari(idx_i,idx_n) += coeff_u*mtx_U_base(idx_i,idx_k);
                std::complex<double> coeff_v = mtx_D(idx_k,idx_n);
                mtx_V_vari(idx_i,idx_n) += coeff_v*mtx_V_base(idx_i,idx_k);
            }
        }
    }

    mtx_output = mtx_init_input;
    MatrixXcd mtx_vari = MatrixXcd::Zero(mtx_init_input.rows(),mtx_init_input.cols());
    mtx_vari += mtx_U_base*mtx_S_base*mtx_V_vari.transpose();
    mtx_vari += mtx_U_vari*mtx_S_base*mtx_V_base.transpose();
    mtx_vari += mtx_U_base*mtx_S_vari*mtx_V_base.transpose();
    mtx_output += mtx_vari;

    if (return_type==0)
    {
        return mtx_output;
    }
    else
    {
        return mtx_t;
    }

}

void fill2DHistogram(TH2D* hist,MatrixXcd mtx)
{
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            hist->SetBinContent(binx+1,biny+1,mtx(binx,biny).real());
        }
    }
}
MatrixXcd fillMatrix(TH2D* hist)
{
    MatrixXcd matrix(hist->GetNbinsX(),hist->GetNbinsY());
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            matrix(binx,biny) = hist->GetBinContent(binx+1,biny+1);
        }
    }
    return matrix;
}

double GetLeastSquareScale(TH2D* hist_on, TH2D* hist_off)
{
    // x_on = a * x_off
    // (x_off)T x_on = a * (x_off)T x_off
    // a = (x_off)T x_on / ( (x_off)T x_off )
    double XoffTXon = 0.;
    double XoffTXoff = 0.;
    for (int binx=0;binx<hist_on->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist_on->GetNbinsY();biny++)
        {
            double on_count = hist_on->GetBinContent(binx+1,biny+1);
            double off_count = hist_off->GetBinContent(binx+1,biny+1);
            double cell_center_x = hist_on->GetXaxis()->GetBinCenter(binx+1);
            double cell_center_y = hist_on->GetYaxis()->GetBinCenter(biny+1);
            double stat_error = 1./pow(max(on_count,1.),0.5);
            if (cell_center_x<MSCL_upper_blind && cell_center_y<MSCW_upper_blind && cell_center_x>MSCL_lower_blind && cell_center_y>MSCW_lower_blind)
            {
                continue;
            }
            if (cell_center_x>2.*MSCL_upper_blind || cell_center_y>2.*MSCW_upper_blind)
            {
                continue;
            }
            XoffTXon += on_count*stat_error*off_count;
            XoffTXoff += off_count*stat_error*off_count;
        }
    }
    double scale = 0.;
    if (XoffTXoff>0.)
    {
        scale = XoffTXon/XoffTXoff;
    }
    return scale;
}

double GetCRcounts(TH2D* hist)
{
    double total_count = 0.;
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            double local_count = hist->GetBinContent(binx+1,biny+1);
            double cell_center_x = hist->GetXaxis()->GetBinCenter(binx+1);
            double cell_center_y = hist->GetYaxis()->GetBinCenter(biny+1);
            if (cell_center_x<MSCL_upper_blind && cell_center_y<MSCW_upper_blind && cell_center_x>MSCL_lower_blind && cell_center_y>MSCW_lower_blind)
            {
                local_count = 0.;
            }
            if (cell_center_x>2.*MSCL_upper_blind || cell_center_y>2.*MSCW_upper_blind)
            {
                local_count = 0.;
            }
            total_count += local_count;
        }
    }
    return total_count;
}
double GetSRcounts(TH2D* hist)
{
    double total_count = 0.;
    for (int binx=0;binx<hist->GetNbinsX();binx++)
    {
        for (int biny=0;biny<hist->GetNbinsY();biny++)
        {
            double local_count = hist->GetBinContent(binx+1,biny+1);
            double cell_center_x = hist->GetXaxis()->GetBinCenter(binx+1);
            double cell_center_y = hist->GetYaxis()->GetBinCenter(biny+1);
            if (!(cell_center_x<MSCL_upper_blind && cell_center_y<MSCW_upper_blind && cell_center_x>MSCL_lower_blind && cell_center_y>MSCW_lower_blind))
            {
                local_count = 0.;
            }
            total_count += local_count;
        }
    }
    return total_count;
}


void FillHistograms(string target_data, bool isON, int doImposter)
{

    sprintf(target, "%s", target_data.c_str());
    ResetPublicVariables(TString(target));

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
        }
        source_strip.append(seglist.at(s));
    }

    pair<double,double> source_ra_dec = GetSourceRaDec(TString(target));
    mean_tele_point_ra = source_ra_dec.first;
    mean_tele_point_dec = source_ra_dec.second;
    map_center_x = mean_tele_point_ra;
    map_center_y = mean_tele_point_dec;

    vector<int> Data_runlist_init = GetRunList(source_strip);
    std::cout << "Data_runlist_init.size() = " << Data_runlist_init.size() << std::endl;
    Data_runlist_init = RemoveNonExistingRuns(Data_runlist_init);
    std::cout << "Data_runlist.size() = " << Data_runlist_init.size() << std::endl;

    std::cout << "Getting runlist elevation..." << std::endl;
    vector<double> Data_runlist_init_elev = GetRunElevationList(Data_runlist_init);
    if (Data_runlist_init.size()==0)
    {
        return;
    }
    std::cout << "Sorting On list by elevation..." << std::endl;
    SortingList(&Data_runlist_init,&Data_runlist_init_elev);

    GetGammaSources();
    GetBrightStars();

    std::cout <<__LINE__ << std::endl;

    double MSCW_bin_size = (MSCW_upper_blind-MSCW_lower_blind)/double(mtx_dim_w_fine);
    double MSCL_bin_size = (MSCL_upper_blind-MSCL_lower_blind)/double(mtx_dim_l_fine);
    double MSCW_plot_upper = MSCW_upper_blind+n_extra_upper_bins*MSCW_bin_size;
    double MSCL_plot_upper = MSCL_upper_blind+n_extra_upper_bins*MSCL_bin_size;
    double MSCW_plot_lower = MSCW_lower_blind-n_extra_lower_bins*MSCW_bin_size;
    double MSCL_plot_lower = MSCL_lower_blind-n_extra_lower_bins*MSCL_bin_size;


    Hist_Data_AreaTime_Skymap = TH2D("Hist_Data_AreaTime_Skymap","",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_x,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y);
    Hist_Data_Norm_Skymap = TH2D("Hist_Data_Norm_Skymap","",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_x,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y);
    Hist_Data_Elev_Skymap = TH2D("Hist_Data_Elev_Skymap","",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_x,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y);
    Hist_Data_Azim_Skymap = TH2D("Hist_Data_Azim_Skymap","",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_x,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y);
    Hist_Data_NSB_Skymap = TH2D("Hist_Data_NSB_Skymap","",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_x,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y);

    for (int e=0;e<N_energy_bins;e++) 
    {
        char e_low[50];
        sprintf(e_low, "%i", int(energy_bins[e]));
        char e_up[50];
        sprintf(e_up, "%i", int(energy_bins[e+1]));
        Hist_OnData_MSCLW.push_back(TH2D("Hist_OnData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",mtx_dim_l,MSCL_plot_lower,MSCL_plot_upper,mtx_dim_w,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OffData_MSCLW.push_back(TH2D("Hist_OffData_MSCLW_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",mtx_dim_l,MSCL_plot_lower,MSCL_plot_upper,mtx_dim_w,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OffData_MSCLW_Sum.push_back(TH2D("Hist_OffData_MSCLW_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",mtx_dim_l,MSCL_plot_lower,MSCL_plot_upper,mtx_dim_w,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnData_MSCLW_Fine.push_back(TH2D("Hist_OnData_MSCLW_Fine_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",mtx_dim_l_fine+n_extra_lower_bins+n_extra_upper_bins,MSCL_plot_lower,MSCL_plot_upper,mtx_dim_w_fine+n_extra_lower_bins+n_extra_upper_bins,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OffData_MSCLW_Fine.push_back(TH2D("Hist_OffData_MSCLW_Fine_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",mtx_dim_l_fine+n_extra_lower_bins+n_extra_upper_bins,MSCL_plot_lower,MSCL_plot_upper,mtx_dim_w_fine+n_extra_lower_bins+n_extra_upper_bins,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnBkgd_MSCLW_Fine.push_back(TH2D("Hist_OnBkgd_MSCLW_Fine_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",mtx_dim_l_fine+n_extra_lower_bins+n_extra_upper_bins,MSCL_plot_lower,MSCL_plot_upper,mtx_dim_w_fine+n_extra_lower_bins+n_extra_upper_bins,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OffData_MSCLW_Fine_Sum.push_back(TH2D("Hist_OffData_MSCLW_Fine_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",mtx_dim_l_fine+n_extra_lower_bins+n_extra_upper_bins,MSCL_plot_lower,MSCL_plot_upper,mtx_dim_w_fine+n_extra_lower_bins+n_extra_upper_bins,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnData_SR_Skymap.push_back(TH2D("Hist_OnData_SR_Skymap_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_SR_Skymap_Mask.push_back(TH2D("Hist_OnData_SR_Skymap_Mask_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_CR_Skymap_Mask.push_back(TH2D("Hist_OnData_CR_Skymap_Mask_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_CR_Skymap_FoV.push_back(TH2D("Hist_OnData_CR_Skymap_FoV_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_CR_Skymap_Ratio.push_back(TH2D("Hist_OnData_CR_Skymap_Ratio_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_CR_Skymap_Regression.push_back(TH2D("Hist_OnData_CR_Skymap_Regression_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_CR_Skymap_Init_Perturbation.push_back(TH2D("Hist_OnData_CR_Skymap_Init_Perturbation_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_CR_Skymap_Perturbation.push_back(TH2D("Hist_OnData_CR_Skymap_Perturbation_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_Expo_Skymap_Sum.push_back(TH2D("Hist_OnData_Expo_Skymap_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",ExpoMap_nbins,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,ExpoMap_nbins,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_SR_Skymap_Sum.push_back(TH2D("Hist_OnData_SR_Skymap_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_CR_Skymap_FoV_Sum.push_back(TH2D("Hist_OnData_CR_Skymap_FoV_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_CR_Skymap_Ratio_Sum.push_back(TH2D("Hist_OnData_CR_Skymap_Ratio_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_CR_Skymap_Regression_Sum.push_back(TH2D("Hist_OnData_CR_Skymap_Regression_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_CR_Skymap_Init_Perturbation_Sum.push_back(TH2D("Hist_OnData_CR_Skymap_Init_Perturbation_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_CR_Skymap_Perturbation_Sum.push_back(TH2D("Hist_OnData_CR_Skymap_Perturbation_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));
        Hist_OnData_CR_Skymap_Combined_Sum.push_back(TH2D("Hist_OnData_CR_Skymap_Combined_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",Skymap_nbins_x,map_center_x-Skymap_size_x,map_center_x+Skymap_size_x,Skymap_nbins_y,map_center_y-Skymap_size_y,map_center_y+Skymap_size_y));

        Hist_OnData_MSCLW_Fine_Sum.push_back(TH2D("Hist_OnData_MSCLW_Fine_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",mtx_dim_l_fine+n_extra_lower_bins+n_extra_upper_bins,MSCL_plot_lower,MSCL_plot_upper,mtx_dim_w_fine+n_extra_lower_bins+n_extra_upper_bins,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnBkgd_MSCLW_Fine_Sum.push_back(TH2D("Hist_OnBkgd_MSCLW_Fine_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",mtx_dim_l_fine+n_extra_lower_bins+n_extra_upper_bins,MSCL_plot_lower,MSCL_plot_upper,mtx_dim_w_fine+n_extra_lower_bins+n_extra_upper_bins,MSCW_plot_lower,MSCW_plot_upper));
        Hist_OnInit_MSCLW_Fine_Sum.push_back(TH2D("Hist_OnInit_MSCLW_Fine_Sum_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",mtx_dim_l_fine+n_extra_lower_bins+n_extra_upper_bins,MSCL_plot_lower,MSCL_plot_upper,mtx_dim_w_fine+n_extra_lower_bins+n_extra_upper_bins,MSCW_plot_lower,MSCW_plot_upper));

        int XYoff_bins = cr_correction_xyoff_bins[e];
        double map_x_lo = -2.;
        double map_x_up = 2.;
        double map_y_lo = -2.;
        double map_y_up = 2.;
        //double map_x_lo = map_center_x-Skymap_size_x;
        //double map_x_up = map_center_x+Skymap_size_x;
        //double map_y_lo = map_center_y-Skymap_size_y;
        //double map_y_up = map_center_y+Skymap_size_y;
        Hist_OnData_SR_XYoff.push_back(TH2D("Hist_OnData_SR_XYoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",XYoff_bins,map_x_lo,map_x_up,XYoff_bins,map_y_lo,map_y_up));
        Hist_OnData_CR_XYoff.push_back(TH2D("Hist_OnData_CR_XYoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",XYoff_bins,map_x_lo,map_x_up,XYoff_bins,map_y_lo,map_y_up));
        Hist_OffData_SR_XYoff.push_back(TH2D("Hist_OffData_SR_XYoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",XYoff_bins,map_x_lo,map_x_up,XYoff_bins,map_y_lo,map_y_up));
        Hist_OffData_CR_XYoff.push_back(TH2D("Hist_OffData_CR_XYoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",XYoff_bins,map_x_lo,map_x_up,XYoff_bins,map_y_lo,map_y_up));
        Hist_OffData_Ratio_XYoff.push_back(TH2D("Hist_OffData_Ratio_XYoff_ErecS"+TString(e_low)+TString("to")+TString(e_up),"",XYoff_bins,map_x_lo,map_x_up,XYoff_bins,map_y_lo,map_y_up));

        CR_on_count_unweighted.push_back(0.);
        CR_off_count_unweighted.push_back(0.);
        CR_on_count_weighted.push_back(0.);
    }

    binx_blind_upper_global = Hist_OnData_MSCLW_Fine.at(0).GetXaxis()->FindBin(MSCL_upper_blind+0.01)-1;
    biny_blind_upper_global = Hist_OnData_MSCLW_Fine.at(0).GetYaxis()->FindBin(MSCW_upper_blind+0.01)-1;
    binx_blind_lower_global = Hist_OnData_MSCLW_Fine.at(0).GetXaxis()->FindBin(MSCL_lower_blind+0.01)-1;
    biny_blind_lower_global = Hist_OnData_MSCLW_Fine.at(0).GetYaxis()->FindBin(MSCW_lower_blind+0.01)-1;

    int n_training_samples = 0;
    int nbins_unblind = 0;
    for (int binx=0;binx<mtx_dim_l;binx++)
    {
        for (int biny=0;biny<mtx_dim_w;biny++)
        {
            double bin_center_x = Hist_OnData_MSCLW.at(0).GetXaxis()->GetBinCenter(binx+1);
            double bin_center_y = Hist_OnData_MSCLW.at(0).GetYaxis()->GetBinCenter(biny+1);
            if (bin_center_x>MSCL_upper_blind || bin_center_y>MSCW_upper_blind || bin_center_x<MSCL_lower_blind || bin_center_y<MSCW_lower_blind)
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

    vector<vector<vector<double>>> off_unblinded_elements;
    vector<vector<double>> off_unblinded_elements_per_sample;
    vector<double> off_unblinded_elements_per_energy;
    vector<vector<double>> off_blinded_element;
    vector<double> off_blinded_element_per_sample;
    double off_blinded_element_per_energy = 0.;

    int n_samples = 0;
    exposure_hours_usable = 0.;
    total_cr_count = 0.;
    off_cr_count = 0.;
    avg_diff_nsb = 0.;
    avg_diff_el = 0.;
    avg_diff_az = 0.;
    char group_tag[50] = "";
    int group_idx = 0;

    for (int e=0;e<N_energy_bins;e++) 
    {
        effective_area.push_back(0.);
    }

    for (int xoff_idx=0;xoff_idx<N_Xoff_bins;xoff_idx++)
    {
        for (int yoff_idx=0;yoff_idx<N_Yoff_bins;yoff_idx++)
        {
            group_idx = 0;
            for (int run=0;run<Data_runlist_init.size();run++)
            {

                bool analyze_this_run = true;
                if (!PointingSelection(int(Data_runlist_init[run])))
                {
                    analyze_this_run = false;
                }

                std::cout << "X = " << xoff_idx << ", Y = " << yoff_idx << std::endl;
                std::cout << "Processing " << run << "/" << Data_runlist_init.size() << "..." << std::endl;

                std::cout << "Getting training data..." << std::endl;

                int ON_RunID = int(Data_runlist_init[run]);
                if (doImposter>0)
                {
                    int ImposterRunID = GetImposterRunID(source_strip, int(Data_runlist_init[run]), doImposter);
                    ON_RunID = ImposterRunID;
                    if (ImposterRunID==0)
                    {
                        std::cout << "Cannot find imposter data for run " << int(Data_runlist_init[run]) << std::endl;
                        analyze_this_run = false;
                    }
                    char imposter_run_number[50];
                    sprintf(imposter_run_number, "%i", ImposterRunID);
                    string imposter_file_name;
                    imposter_file_name = TString(SMI_INPUT+"/"+string(imposter_run_number)+".anasum.root");
                    if (gSystem->AccessPathName(imposter_file_name.c_str()))
                    {
                        std::cout << "Imposter Run " << imposter_run_number << " does not exist!!!" << std::endl;
                        analyze_this_run = false;
                    }
                }
                vector<int> OffData_runlist_init;
                if (doImposter>0)
                {
                    OffData_runlist_init = GetImposterOffRunList(source_strip,ON_RunID);
                }
                else
                {
                    OffData_runlist_init = GetOffRunList(source_strip,ON_RunID);
                }
                std::cout << "Checking any non-existing files..." << std::endl;
                OffData_runlist_init = RemoveNonExistingRuns(OffData_runlist_init);

                std::cout << "Run " << ON_RunID << " has " << OffData_runlist_init.size() << " OFF runs." << std::endl;
                if (OffData_runlist_init.size()<1)
                {
                    std::cout << "Insufficient training data..." << std::endl;
                    analyze_this_run = false;
                }

                if (analyze_this_run)
                {

                    n_samples += 1;
                    for (int offrun=0;offrun<OffData_runlist_init.size();offrun++)
                    {

                        if (!PointingSelection(int(OffData_runlist_init[offrun])))
                        {
                            continue;
                        }
                        for (int offrun2=0;offrun2<offrun;offrun2++)
                        {
                            if (int(OffData_runlist_init[offrun])==int(OffData_runlist_init[offrun2]))
                            {
                                std::cout << "OFF run " << int(OffData_runlist_init[offrun]) << " has been used. Remove it from training sample." << std::endl;
                                continue;
                            }
                        }

                        double on_run_cr_count = CountCosmicRayEvents(ON_RunID, xoff_idx, yoff_idx);
                        double off_run_cr_count = CountCosmicRayEvents(int(OffData_runlist_init[offrun]), xoff_idx, yoff_idx);
                        double off_evt_weight = 0.;
                        if (off_run_cr_count>0.)
                        {
                            off_evt_weight = on_run_cr_count/off_run_cr_count;
                        }

                        TrainingRunAnalysis(int(OffData_runlist_init[offrun]), ON_RunID, xoff_idx, yoff_idx, off_evt_weight);

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
                                    if (bin_center_x>MSCL_upper_blind || bin_center_y>MSCW_upper_blind || bin_center_x<MSCL_lower_blind || bin_center_y<MSCW_lower_blind)
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
                            Hist_OffData_MSCLW_Sum.at(e).Add(&Hist_OffData_MSCLW.at(e));
                            Hist_OffData_MSCLW_Fine_Sum.at(e).Add(&Hist_OffData_MSCLW_Fine.at(e));
                            Hist_OffData_MSCLW_Fine.at(e).Reset();
                            Hist_OffData_MSCLW.at(e).Reset();
                        }
                    }

                    for (int e=0;e<N_energy_bins;e++) 
                    {
                        for (int binx=0;binx<Hist_OffData_SR_XYoff.at(e).GetNbinsX();binx++)
                        {
                            for (int biny=0;biny<Hist_OffData_SR_XYoff.at(e).GetNbinsY();biny++)
                            {
                                double sr_count = Hist_OffData_SR_XYoff.at(e).GetBinContent(binx+1,biny+1);
                                double cr_count = Hist_OffData_CR_XYoff.at(e).GetBinContent(binx+1,biny+1);
                                double srcr_ratio = 0.;
                                if (cr_count>0)
                                {
                                    srcr_ratio = sr_count/cr_count;
                                }
                                Hist_OffData_Ratio_XYoff.at(e).SetBinContent(binx+1,biny+1,srcr_ratio);
                            }
                        }
                    }

                    std::cout << "Getting ON data..." << std::endl;
                    if (doImposter>0)
                    {
                        SingleRunAnalysis(ON_RunID, int(Data_runlist_init[run]), xoff_idx, yoff_idx, true);
                    }
                    else
                    {
                        SingleRunAnalysis(ON_RunID, int(Data_runlist_init[run]), xoff_idx, yoff_idx, false);
                    }

                    for (int e=0;e<N_energy_bins;e++) 
                    {
                        Hist_OnData_SR_XYoff.at(e).Reset();
                        Hist_OnData_CR_XYoff.at(e).Reset();
                        Hist_OffData_SR_XYoff.at(e).Reset();
                        Hist_OffData_CR_XYoff.at(e).Reset();
                    }

                }

                //if (exposure_hours_usable>=expo_hour_per_analysis || run==Data_runlist_init.size()-1)
                if (total_cr_count>min_CR_count || run==Data_runlist_init.size()-1)
                {

                    std::cout << "Using training data to learn conversion..." << std::endl;

                    n_training_samples = off_blinded_element.size();
                    std::cout << "n_training_samples = " << n_training_samples << std::endl;
                    std::cout << "nbins_unblind = " << nbins_unblind << std::endl;

                    if (n_training_samples<nbins_unblind || n_training_samples<n_samples)
                    {
                        std::cout << "Cannot make prediction! Insufficient training data." << std::endl;
                        continue;
                    }

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
                        if (nbins_unblind>n_training_samples) continue;
                        VectorXcd vtr_B = VectorXcd::Zero(n_training_samples);
                        MatrixXcd mtx_W = MatrixXcd::Zero(n_training_samples,n_training_samples);
                        MatrixXcd mtx_A = MatrixXcd::Zero(n_training_samples,nbins_unblind);
                        VectorXcd vtr_X = VectorXcd::Zero(nbins_unblind);

                        for (int sample=0;sample<n_training_samples;sample++)
                        {
                            vtr_B(sample) = off_blinded_element.at(sample).at(e);
                            if (off_blinded_element.at(sample).at(e)>0.)
                            {
                                double weight = 1.;
                                if (use_stat_err_weight) 
                                {
                                    weight = 1./pow(max(1.,off_blinded_element.at(sample).at(e)),0.5);
                                }
                                mtx_W(sample,sample) = weight;
                            }
                            for (int bin=0;bin<nbins_unblind;bin++)
                            {
                                mtx_A(sample,bin) = off_unblinded_elements.at(sample).at(e).at(bin);
                            }
                        }

                        VectorXcd vtr_B_weighted = mtx_A.transpose()*mtx_W*vtr_B;
                        MatrixXcd mtx_A_weighted = mtx_A.transpose()*mtx_W*mtx_A; 
                        BDCSVD<MatrixXcd> bdc_svd(mtx_A_weighted, ComputeThinU | ComputeThinV);
                        vtr_X = bdc_svd.solve(vtr_B_weighted);
                        convert_unblind_to_blind_regression.push_back(vtr_X);
                    }

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
                                if (bin_center_x>MSCL_upper_blind || bin_center_y>MSCW_upper_blind || bin_center_x<MSCL_lower_blind || bin_center_y<MSCW_lower_blind)
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


                    for (int e=0;e<N_energy_bins;e++) 
                    {

                        double CR_count_map = 0.;

                        // FoV method
                        double SR_count_map = Hist_OnData_SR_Skymap_Mask.at(e).Integral();
                        CR_count_map = Hist_OnData_CR_Skymap_Mask.at(e).Integral();
                        if (CR_count_map>0.)
                        {
                            Hist_OnData_CR_Skymap_FoV.at(e).Scale(SR_count_map/CR_count_map);
                        }
                        
                        // simple ratio method
                        double total_on_cr_count = GetCRcounts(&Hist_OnData_MSCLW_Fine.at(e));
                        double total_off_cr_count = GetCRcounts(&Hist_OffData_MSCLW_Fine_Sum.at(e));
                        double total_off_sr_count = GetSRcounts(&Hist_OffData_MSCLW_Fine_Sum.at(e));
                        double least_square_scale = GetLeastSquareScale(&Hist_OnData_MSCLW_Fine.at(e), &Hist_OffData_MSCLW_Fine_Sum.at(e));
                        TH2D Hist_OffData_MSCLW_CR_scaled = TH2D("Hist_OffData_MSCLW_CR_scaled","",mtx_dim_l_fine+n_extra_lower_bins+n_extra_upper_bins,MSCL_plot_lower,MSCL_plot_upper,mtx_dim_w_fine+n_extra_lower_bins+n_extra_upper_bins,MSCW_plot_lower,MSCW_plot_upper);
                        //if (total_off_cr_count>0.)
                        //{
                        //    Hist_OffData_MSCLW_CR_scaled.Add(&Hist_OffData_MSCLW_Fine_Sum.at(e),total_on_cr_count/total_off_cr_count);
                        //}
                        Hist_OffData_MSCLW_CR_scaled.Add(&Hist_OffData_MSCLW_Fine_Sum.at(e),least_square_scale);

                        //CR_count_map = Hist_OnData_CR_Skymap_Ratio.at(e).Integral();
                        double SR_predict_ratio = Hist_OffData_MSCLW_CR_scaled.Integral(binx_blind_lower_global+1,binx_blind_upper_global,biny_blind_lower_global+1,biny_blind_upper_global);
                        CR_count_map = CR_on_count_weighted.at(e);
                        if (CR_count_map>0.)
                        {
                            Hist_OnData_CR_Skymap_Ratio.at(e).Scale(SR_predict_ratio/CR_count_map);
                        }

                        // linear regression method
                        // A * X = B
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

                        //CR_count_map = Hist_OnData_CR_Skymap_Regression.at(e).Integral();
                        CR_count_map = CR_on_count_weighted.at(e);
                        if (CR_count_map>0.)
                        {
                            Hist_OnData_CR_Skymap_Regression.at(e).Scale(SR_predict_regression/CR_count_map);
                        }

                        // perturbation method
                        MatrixXcd mtx_off_data_cr_scaled = fillMatrix(&Hist_OffData_MSCLW_CR_scaled);
                        MatrixXcd mtx_on_data = fillMatrix(&Hist_OnData_MSCLW_Fine.at(e));

                        std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
                        std::cout << "Energy = " << energy_bins[e] << std::endl;
                        MatrixXcd mtx_on_bkgd = mtx_off_data_cr_scaled;
                        fill2DHistogram(&Hist_OnBkgd_MSCLW_Fine.at(e),mtx_on_bkgd);
                        double total_data_sr_count = GetSRcounts(&Hist_OnData_MSCLW_Fine.at(e));
                        double total_bkgd_sr_count = GetSRcounts(&Hist_OnBkgd_MSCLW_Fine.at(e));
                        if (isNaN(total_bkgd_sr_count) || total_bkgd_sr_count==0.)
                        {
                            fill2DHistogram(&Hist_OnBkgd_MSCLW_Fine.at(e),mtx_off_data_cr_scaled);
                            total_bkgd_sr_count = GetSRcounts(&Hist_OnBkgd_MSCLW_Fine.at(e));
                        }
                        std::cout << "total_data_sr_count = " << total_data_sr_count << ", total_bkgd_sr_count = " << total_bkgd_sr_count << std::endl;

                        //int max_rank = matrix_rank[e];
                        int max_rank = matrix_rank;
                        double log_t_mtx_weight = 3.0;

                        if (total_data_sr_count==0)
                        {
                            Hist_OnData_MSCLW_Fine.at(e).Reset();
                            mtx_on_data = fillMatrix(&Hist_OnData_MSCLW_Fine.at(e));
                        }

                        MatrixXcd mtx_on_t = mtx_off_data_cr_scaled;
                        for (int idx_i=0;idx_i<mtx_on_t.rows();idx_i++)
                        {
                            for (int idx_j=0;idx_j<mtx_on_t.rows();idx_j++)
                            {
                                mtx_on_t(idx_i,idx_j) = 0.;
                            }
                        }

                        bool use_t_input = false;
                        bool use_diagonal = true;
                        bool is_blind = false;

                        max_rank = 2;
                        mtx_on_t = MatrixPerturbationMethod(e,mtx_on_t,mtx_off_data_cr_scaled,mtx_off_data_cr_scaled,mtx_on_data,max_rank,max_rank,use_diagonal,use_t_input,1.0,is_blind,1);
                        std::cout << "mtx_on_t (unblind):" << std::endl;
                        std::cout << mtx_on_t.block(0,0,3,3).real() << std::endl;
                        mtx_on_bkgd = MatrixPerturbationMethod(e,mtx_on_t,mtx_off_data_cr_scaled,mtx_off_data_cr_scaled,mtx_on_data,max_rank,max_rank,use_diagonal,use_t_input,1.0,is_blind,0);
                        fill2DHistogram(&Hist_OnBkgd_MSCLW_Fine.at(e),mtx_on_bkgd);
                        total_data_sr_count = GetSRcounts(&Hist_OnData_MSCLW_Fine.at(e));
                        total_bkgd_sr_count = GetSRcounts(&Hist_OnBkgd_MSCLW_Fine.at(e));
                        if (isNaN(total_bkgd_sr_count) || total_bkgd_sr_count==0.)
                        {
                            fill2DHistogram(&Hist_OnBkgd_MSCLW_Fine.at(e),mtx_off_data_cr_scaled);
                            total_bkgd_sr_count = GetSRcounts(&Hist_OnBkgd_MSCLW_Fine.at(e));
                        }
                        std::cout << "total_data_sr_count = " << total_data_sr_count << ", total_bkgd_sr_count = " << total_bkgd_sr_count << std::endl;
                        t00_truth.push_back(mtx_on_t(0,0).real());
                        t01_truth.push_back(mtx_on_t(0,1).real());
                        t10_truth.push_back(mtx_on_t(1,0).real());
                        t11_truth.push_back(mtx_on_t(1,1).real());

                        is_blind = true;
                        max_rank = 2;
                        if (max_rank>=1)
                        {

                            use_t_input = false;
                            use_diagonal = false;
                            log_t_mtx_weight = 3.0;

                            mtx_on_t = MatrixPerturbationMethod(e,mtx_on_t,mtx_off_data_cr_scaled,mtx_off_data_cr_scaled,mtx_on_data,max_rank,max_rank,use_diagonal,use_t_input,log_t_mtx_weight,is_blind,1);
                            std::cout << "mtx_on_t (blind, off-diagonal):" << std::endl;
                            std::cout << mtx_on_t.block(0,0,3,3).real() << std::endl;
                            use_t_input = true;
                            mtx_on_bkgd = MatrixPerturbationMethod(e,mtx_on_t,mtx_off_data_cr_scaled,mtx_off_data_cr_scaled,mtx_on_data,max_rank,max_rank,use_diagonal,use_t_input,3.0,is_blind,0);
                            fill2DHistogram(&Hist_OnBkgd_MSCLW_Fine.at(e),mtx_on_bkgd);
                            total_data_sr_count = GetSRcounts(&Hist_OnData_MSCLW_Fine.at(e));
                            total_bkgd_sr_count = GetSRcounts(&Hist_OnBkgd_MSCLW_Fine.at(e));
                            if (isNaN(total_bkgd_sr_count) || total_bkgd_sr_count==0.)
                            {
                                fill2DHistogram(&Hist_OnBkgd_MSCLW_Fine.at(e),mtx_off_data_cr_scaled);
                                total_bkgd_sr_count = GetSRcounts(&Hist_OnBkgd_MSCLW_Fine.at(e));
                            }
                            std::cout << "total_data_sr_count = " << total_data_sr_count << ", total_bkgd_sr_count = " << total_bkgd_sr_count << std::endl;
                        }

                        if (max_rank==1 || max_rank==0)
                        {
                            fill2DHistogram(&Hist_OnBkgd_MSCLW_Fine.at(e),mtx_off_data_cr_scaled);
                        }
                        double SR_predict_init_perturbation = Hist_OnBkgd_MSCLW_Fine.at(e).Integral(binx_blind_lower_global+1,binx_blind_upper_global,biny_blind_lower_global+1,biny_blind_upper_global);
                        CR_count_map = CR_on_count_weighted.at(e);
                        if (CR_count_map>0.)
                        {
                            Hist_OnData_CR_Skymap_Init_Perturbation.at(e).Scale(SR_predict_init_perturbation/CR_count_map);
                        }

                        is_blind = true;
                        max_rank = 2;
                        if (max_rank>=1)
                        {
                            use_t_input = true;
                            use_diagonal = true;
                            log_t_mtx_weight = 3.0;

                            for (int rank=1;rank<=2;rank++)
                            {
                                mtx_on_t = MatrixPerturbationMethod(e,mtx_on_t,mtx_off_data_cr_scaled,mtx_off_data_cr_scaled,mtx_on_data,max_rank,rank,use_diagonal,use_t_input,log_t_mtx_weight,is_blind,1);
                                std::cout << "mtx_on_t (blind, diagonal):" << std::endl;
                                std::cout << mtx_on_t.block(0,0,3,3).real() << std::endl;
                                use_t_input = true;
                            }
                            //use_t_input = false;
                            //mtx_on_t = MatrixPerturbationMethod(e,mtx_on_t,mtx_off_data_cr_scaled,mtx_off_data_cr_scaled,mtx_on_data,max_rank,max_rank,use_diagonal,use_t_input,log_t_mtx_weight,is_blind,1);
                            //std::cout << "mtx_on_t (blind, diagonal):" << std::endl;
                            //std::cout << mtx_on_t.block(0,0,3,3).real() << std::endl;
                            //use_t_input = true;

                            mtx_on_bkgd = MatrixPerturbationMethod(e,mtx_on_t,mtx_off_data_cr_scaled,mtx_off_data_cr_scaled,mtx_on_data,max_rank,max_rank,use_diagonal,use_t_input,3.0,is_blind,0);

                            fill2DHistogram(&Hist_OnBkgd_MSCLW_Fine.at(e),mtx_on_bkgd);
                            total_data_sr_count = GetSRcounts(&Hist_OnData_MSCLW_Fine.at(e));
                            total_bkgd_sr_count = GetSRcounts(&Hist_OnBkgd_MSCLW_Fine.at(e));
                            if (isNaN(total_bkgd_sr_count) || total_bkgd_sr_count==0.)
                            {
                                fill2DHistogram(&Hist_OnBkgd_MSCLW_Fine.at(e),mtx_off_data_cr_scaled);
                                total_bkgd_sr_count = GetSRcounts(&Hist_OnBkgd_MSCLW_Fine.at(e));
                            }
                            std::cout << "total_data_sr_count = " << total_data_sr_count << ", total_bkgd_sr_count = " << total_bkgd_sr_count << std::endl;

                        }
                        t00_recon.push_back(mtx_on_t(0,0).real());
                        t01_recon.push_back(mtx_on_t(0,1).real());
                        t10_recon.push_back(mtx_on_t(1,0).real());
                        t11_recon.push_back(mtx_on_t(1,1).real());

                        if (max_rank==0)
                        {
                            fill2DHistogram(&Hist_OnBkgd_MSCLW_Fine.at(e),mtx_off_data_cr_scaled);
                        }
                        double SR_predict_perturbation = Hist_OnBkgd_MSCLW_Fine.at(e).Integral(binx_blind_lower_global+1,binx_blind_upper_global,biny_blind_lower_global+1,biny_blind_upper_global);
                        CR_count_map = CR_on_count_weighted.at(e);
                        if (CR_count_map>0.)
                        {
                            Hist_OnData_CR_Skymap_Perturbation.at(e).Scale(SR_predict_perturbation/CR_count_map);
                        }

                        Hist_OnData_MSCLW_Fine_Sum.at(e).Add(&Hist_OnData_MSCLW_Fine.at(e));
                        Hist_OnBkgd_MSCLW_Fine_Sum.at(e).Add(&Hist_OnBkgd_MSCLW_Fine.at(e));
                        Hist_OnInit_MSCLW_Fine_Sum.at(e).Add(&Hist_OffData_MSCLW_CR_scaled);

                    }


                    for (int e=0;e<N_energy_bins;e++) 
                    {
                        Hist_OnData_SR_Skymap_Sum.at(e).Add(&Hist_OnData_SR_Skymap.at(e));
                        Hist_OnData_CR_Skymap_FoV_Sum.at(e).Add(&Hist_OnData_CR_Skymap_FoV.at(e));
                        Hist_OnData_CR_Skymap_Ratio_Sum.at(e).Add(&Hist_OnData_CR_Skymap_Ratio.at(e));
                        Hist_OnData_CR_Skymap_Regression_Sum.at(e).Add(&Hist_OnData_CR_Skymap_Regression.at(e));
                        Hist_OnData_CR_Skymap_Init_Perturbation_Sum.at(e).Add(&Hist_OnData_CR_Skymap_Init_Perturbation.at(e));
                        Hist_OnData_CR_Skymap_Perturbation_Sum.at(e).Add(&Hist_OnData_CR_Skymap_Perturbation.at(e));

                        Hist_OnData_CR_Skymap_Combined_Sum.at(e).Add(&Hist_OnData_CR_Skymap_Regression.at(e),0.5);
                        Hist_OnData_CR_Skymap_Combined_Sum.at(e).Add(&Hist_OnData_CR_Skymap_Perturbation.at(e),0.5);
                    }

                    off_unblinded_elements.clear();
                    off_unblinded_elements_per_sample.clear();
                    off_unblinded_elements_per_energy.clear();
                    off_blinded_element.clear();
                    off_blinded_element_per_sample.clear();
                    off_blinded_element_per_energy = 0.;

                    for (int e=0;e<N_energy_bins;e++) 
                    {
                        Hist_OffData_MSCLW_Fine_Sum.at(e).Reset();
                        Hist_OffData_MSCLW_Sum.at(e).Reset();
                        Hist_OnData_MSCLW.at(e).Reset();
                        Hist_OnData_MSCLW_Fine.at(e).Reset();
                        Hist_OnBkgd_MSCLW_Fine.at(e).Reset();
                        Hist_OnData_SR_Skymap.at(e).Reset();
                        Hist_OnData_SR_Skymap_Mask.at(e).Reset();
                        Hist_OnData_CR_Skymap_Mask.at(e).Reset();
                        Hist_OnData_CR_Skymap_FoV.at(e).Reset();
                        Hist_OnData_CR_Skymap_Ratio.at(e).Reset();
                        Hist_OnData_CR_Skymap_Regression.at(e).Reset();
                        Hist_OnData_CR_Skymap_Init_Perturbation.at(e).Reset();
                        Hist_OnData_CR_Skymap_Perturbation.at(e).Reset();
                        CR_on_count_unweighted.at(e) = 0.;
                        CR_off_count_unweighted.at(e) = 0.;
                        CR_on_count_weighted.at(e) = 0.;
                    }

                    NSB_mean = NSB_mean/double(n_samples);
                    Elev_mean = Elev_mean/double(n_samples);
                    Azim_mean = Azim_mean/double(n_samples);
                    for (int e=0;e<N_energy_bins;e++) 
                    {
                        effective_area.at(e) = effective_area.at(e)/double(n_samples);
                    }

                    Hist_Data_Elev_Skymap.Divide(&Hist_Data_Norm_Skymap);
                    Hist_Data_Azim_Skymap.Divide(&Hist_Data_Norm_Skymap);
                    Hist_Data_NSB_Skymap.Divide(&Hist_Data_Norm_Skymap);

                    sprintf(group_tag, "_G%d_X%d_Y%d", group_idx, xoff_idx, yoff_idx);
                    TFile OutputFile(TString(SMI_OUTPUT)+"/Netflix_"+TString(target)+group_tag+".root","recreate");
                    group_idx += 1;

                    Hist_Data_AreaTime_Skymap.Write();
                    Hist_Data_Norm_Skymap.Write();
                    Hist_Data_Elev_Skymap.Write();
                    Hist_Data_Azim_Skymap.Write();
                    Hist_Data_NSB_Skymap.Write();

                    for (int e=0;e<N_energy_bins;e++) 
                    {
                        double truth = Hist_OnData_SR_Skymap_Sum.at(e).Integral();
                        double ratio_predict = Hist_OnData_CR_Skymap_Ratio_Sum.at(e).Integral();
                        double regression_predict = Hist_OnData_CR_Skymap_Regression_Sum.at(e).Integral();
                        double init_perturbation_predict = Hist_OnData_CR_Skymap_Init_Perturbation_Sum.at(e).Integral();
                        double perturbation_predict = Hist_OnData_CR_Skymap_Perturbation_Sum.at(e).Integral();
                        double combined_predict = Hist_OnData_CR_Skymap_Combined_Sum.at(e).Integral();
                        data_count.push_back(truth);
                        ratio_bkgd_count.push_back(ratio_predict);
                        regression_bkgd_count.push_back(regression_predict);
                        init_perturbation_bkgd_count.push_back(init_perturbation_predict);
                        perturbation_bkgd_count.push_back(perturbation_predict);
                        combined_bkgd_count.push_back(combined_predict);

                        std::cout << "=============================================" << std::endl;
                        std::cout << "Energy = " << energy_bins[e] << std::endl;
                        std::cout << "Effective area = " << effective_area.at(e) << std::endl;
                        std::cout << 
                            "truth                = " << truth
                            << std::endl;
                        std::cout << 
                            "ratio_predict        = " << ratio_predict 
                            << ", " << -(ratio_predict-truth)/truth*100. << " %"
                            << ", " << -(ratio_predict-truth)/pow(truth,0.5) << " sigma"
                            << std::endl;
                        std::cout << 
                            "init_perturbation_predict = " << init_perturbation_predict 
                            << ", " << -(init_perturbation_predict-truth)/truth*100. << " %"
                            << ", " << -(init_perturbation_predict-truth)/pow(truth,0.5) << " sigma"
                            << std::endl;
                        std::cout << 
                            "perturbation_predict = " << perturbation_predict 
                            << ", " << -(perturbation_predict-truth)/truth*100. << " %"
                            << ", " << -(perturbation_predict-truth)/pow(truth,0.5) << " sigma"
                            << std::endl;
                        std::cout << 
                            "regression_predict   = " << regression_predict 
                            << ", " << -(regression_predict-truth)/truth*100. << " %"
                            << ", " << -(regression_predict-truth)/pow(truth,0.5) << " sigma"
                            << std::endl;
                        std::cout << 
                            "combined_predict   = " << regression_predict 
                            << ", " << -(combined_predict-truth)/truth*100. << " %"
                            << ", " << -(combined_predict-truth)/pow(truth,0.5) << " sigma"
                            << std::endl;
                        Hist_OnData_Expo_Skymap_Sum.at(e).Write();
                        Hist_OnData_SR_Skymap_Sum.at(e).Write();
                        Hist_OnData_CR_Skymap_FoV_Sum.at(e).Write();
                        Hist_OnData_CR_Skymap_Ratio_Sum.at(e).Write();
                        Hist_OnData_CR_Skymap_Regression_Sum.at(e).Write();
                        Hist_OnData_CR_Skymap_Init_Perturbation_Sum.at(e).Write();
                        Hist_OnData_CR_Skymap_Perturbation_Sum.at(e).Write();
                        Hist_OnData_CR_Skymap_Combined_Sum.at(e).Write();
                        Hist_OnData_MSCLW_Fine_Sum.at(e).Write();
                        Hist_OnBkgd_MSCLW_Fine_Sum.at(e).Write();
                        Hist_OnInit_MSCLW_Fine_Sum.at(e).Write();
                    }

                    avg_diff_nsb = avg_diff_nsb/off_cr_count;
                    avg_diff_el = avg_diff_el/off_cr_count;
                    avg_diff_az = avg_diff_az/off_cr_count;

                    TTree InfoTree("InfoTree","info tree");
                    InfoTree.Branch("total_cr_count",&total_cr_count,"total_cr_count/D");
                    InfoTree.Branch("avg_diff_nsb",&avg_diff_nsb,"avg_diff_nsb/D");
                    InfoTree.Branch("avg_diff_el",&avg_diff_el,"avg_diff_el/D");
                    InfoTree.Branch("avg_diff_az",&avg_diff_az,"avg_diff_az/D");
                    InfoTree.Branch("exposure_hours",&exposure_hours_usable,"exposure_hours/D");
                    InfoTree.Branch("NSB_mean",&NSB_mean,"NSB_mean/D");
                    InfoTree.Branch("Elev_mean",&Elev_mean,"Elev_mean/D");
                    InfoTree.Branch("Azim_mean",&Azim_mean,"Azim_mean/D");
                    InfoTree.Branch("effective_area","std::vector<double>",&effective_area);
                    InfoTree.Branch("data_count","std::vector<double>",&data_count);
                    InfoTree.Branch("ratio_bkgd_count","std::vector<double>",&ratio_bkgd_count);
                    InfoTree.Branch("regression_bkgd_count","std::vector<double>",&regression_bkgd_count);
                    InfoTree.Branch("init_perturbation_bkgd_count","std::vector<double>",&init_perturbation_bkgd_count);
                    InfoTree.Branch("perturbation_bkgd_count","std::vector<double>",&perturbation_bkgd_count);
                    InfoTree.Branch("combined_bkgd_count","std::vector<double>",&combined_bkgd_count);
                    InfoTree.Branch("t00_truth","std::vector<double>",&t00_truth);
                    InfoTree.Branch("t01_truth","std::vector<double>",&t01_truth);
                    InfoTree.Branch("t10_truth","std::vector<double>",&t10_truth);
                    InfoTree.Branch("t11_truth","std::vector<double>",&t11_truth);
                    InfoTree.Branch("t00_recon","std::vector<double>",&t00_recon);
                    InfoTree.Branch("t01_recon","std::vector<double>",&t01_recon);
                    InfoTree.Branch("t10_recon","std::vector<double>",&t10_recon);
                    InfoTree.Branch("t11_recon","std::vector<double>",&t11_recon);
                    InfoTree.Fill();
                    InfoTree.Write();

                    exposure_hours_usable = 0.;
                    total_cr_count = 0.;
                    off_cr_count = 0.;
                    avg_diff_nsb = 0.;
                    avg_diff_el = 0.;
                    avg_diff_az = 0.;
                    n_samples = 0;
                    effective_area.clear();
                    NSB_mean = 0.;
                    Elev_mean = 0.;
                    Azim_mean = 0.;
                    for (int e=0;e<N_energy_bins;e++) 
                    {
                        effective_area.push_back(0.);
                    }
                    t00_truth.clear();
                    t01_truth.clear();
                    t10_truth.clear();
                    t11_truth.clear();
                    t00_recon.clear();
                    t01_recon.clear();
                    t10_recon.clear();
                    t11_recon.clear();
                    data_count.clear();
                    ratio_bkgd_count.clear();
                    regression_bkgd_count.clear();
                    init_perturbation_bkgd_count.clear();
                    perturbation_bkgd_count.clear();
                    combined_bkgd_count.clear();

                    Hist_Data_AreaTime_Skymap.Reset();
                    Hist_Data_Norm_Skymap.Reset();
                    Hist_Data_Elev_Skymap.Reset();
                    Hist_Data_Azim_Skymap.Reset();
                    Hist_Data_NSB_Skymap.Reset();
                    for (int e=0;e<N_energy_bins;e++) 
                    {
                        Hist_OnData_Expo_Skymap_Sum.at(e).Reset();
                        Hist_OnData_SR_Skymap_Sum.at(e).Reset();
                        Hist_OnData_CR_Skymap_FoV_Sum.at(e).Reset();
                        Hist_OnData_CR_Skymap_Ratio_Sum.at(e).Reset();
                        Hist_OnData_CR_Skymap_Regression_Sum.at(e).Reset();
                        Hist_OnData_CR_Skymap_Init_Perturbation_Sum.at(e).Reset();
                        Hist_OnData_CR_Skymap_Perturbation_Sum.at(e).Reset();
                        Hist_OnData_CR_Skymap_Combined_Sum.at(e).Reset();
                        Hist_OnData_MSCLW_Fine_Sum.at(e).Reset();
                        Hist_OnBkgd_MSCLW_Fine_Sum.at(e).Reset();
                        Hist_OnInit_MSCLW_Fine_Sum.at(e).Reset();
                    }

                }

            }
        }
    }

}
