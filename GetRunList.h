
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include "NetflixParameters.h"

string RUNLIST_DIR = string(std::getenv("SMI_RUNLIST"));

pair<double,double> GetSourceRaDec(TString source_name)
{
    double Source_RA = 0.;
    double Source_Dec = 0.;
    if (source_name.Contains("1ES1959_p650"))
    {
            Source_RA = 300.00;
            Source_Dec = 65.15;
    }
    if (source_name.Contains("PSR_J2229_p6114"))
    {
            Source_RA = 337.27;
            Source_Dec = 61.23;
    }
    if (source_name.Contains("PSR_B1937_p21"))
    {
            Source_RA = 295.45;
            Source_Dec = 21.44;
    }
    if (source_name.Contains("G079_p0p3"))
    {
            Source_RA = 308.12;
            Source_Dec = 40.33;
    }
    if (source_name.Contains("2FHL_J0431_p5553e"))
    {
            Source_RA = 67.81;
            Source_Dec = 55.89;
    }
    if (source_name.Contains("PSR_J2238_p5903"))
    {
            Source_RA = 339.6079;
            Source_Dec = 59.0487;
    }
    if (source_name.Contains("PSR_J1928_p1746"))
    {
            Source_RA = 292.15;
            Source_Dec = 17.78;
    }
    if (source_name.Contains("PSR_J0631_p1036"))
    {
            Source_RA = 97.9;
            Source_Dec = 10.6;
    }
    if (source_name.Contains("Cisne_HS_2013"))
    {
            Source_RA = 304.08;
            Source_Dec = 36.13;
    }
    if (source_name.Contains("SNR_G150_p4"))
    {
            Source_RA = 67.809;
            Source_Dec = 55.889;
    }
    if (source_name.Contains("Mrk421"))
    {
            Source_RA = 166.0791667;
            Source_Dec = 38.1947222;
    }
    if (source_name.Contains("Tycho"))
    {
            Source_RA = 6.28;
            Source_Dec = 64.17;
    }
    if (source_name.Contains("CTA1"))
    {
            Source_RA = 1.608;
            Source_Dec = 72.983;
    }
    if (source_name.Contains("Cas_A"))
    {
            Source_RA = 350.8075;
            Source_Dec = 58.8072;
    }
    if (source_name.Contains("PSR_J1856_p0245"))
    {
            Source_RA = 284.21;
            Source_Dec = 2.76;
    }
    if (source_name.Contains("LHAASO_J0622_p3754"))
    {
            Source_RA = 95.50;
            Source_Dec = 37.90;
    }
    if (source_name.Contains("PSR_J0359_p5414"))
    {
            Source_RA = 59.858;
            Source_Dec = 54.248;
    }
    if (source_name.Contains("LHAASO_J1956_p2845"))
    {
            Source_RA = 299.05;
            Source_Dec = 28.75;
    }
    if (source_name.Contains("PSR_J2032_p4127"))
    {
            Source_RA = 308.05;
            Source_Dec = 41.46;
    }
    if (source_name.Contains("CTB109"))
    {
            Source_RA = 345.28;
            Source_Dec = 58.88;
    }
    if (source_name.Contains("SNR_G067p6_p0p9"))
    {
            Source_RA = 299.44;
            Source_Dec = 30.88;
    }
    if (source_name.Contains("2HWC_J1953_p294"))
    {
            Source_RA = 298.26;
            Source_Dec = 29.48;
    }
    if (source_name.Contains("SgrA"))
    {
            Source_RA = 266.41;
            Source_Dec = -29.00;
    }
    if (source_name.Contains("PSR_J2021_p3651"))
    {
            Source_RA = 305.27-1.0;
            Source_Dec = 36.85;
    }
    if (source_name.Contains("Geminga"))
    {
            Source_RA = 98.476;
            Source_Dec = 17.770;
    }
    if (source_name.Contains("SNR_G189_p03"))
    {
            Source_RA = 94.213;
            Source_Dec = 22.503;
    }
    if (source_name.Contains("PSR_J2021_p4026"))
    {
            Source_RA = 305.37;
            Source_Dec = 40.45;
    }
    if (source_name.Contains("SS433"))
    {
            Source_RA = 288.0833333;
            Source_Dec = 4.9166667;
    }
    if (source_name.Contains("PSR_J1907_p0602"))
    {
            //Source_RA = 286.975;
            //Source_Dec = 6.03777777778;
            Source_RA = 286.975;
            Source_Dec = 6.33777777778;
    }
    if (source_name.Contains("M82"))
    {
            Source_RA = 148.970;
            Source_Dec = 69.679;
    }
    if (source_name.Contains("1ES0647"))
    {
            Source_RA = 102.694;
            Source_Dec = 25.050;
    }
    if (source_name.Contains("Sky_RA55Dec53"))
    {
            Source_RA = 55.34;
            Source_Dec = 52.97;
    }
    if (source_name.Contains("Crab"))
    {
            Source_RA = 83.633;
            Source_Dec = 22.014;
    }
    if (source_name.Contains("PKS1424"))
    {
            Source_RA = 216.750;
            Source_Dec = 23.783;
    }
    if (source_name.Contains("1ES0229"))
    {
            Source_RA = 38.222;
            Source_Dec = 20.273;
    }
    if (source_name.Contains("H1426"))
    {
            Source_RA = 217.136;
            Source_Dec = 42.673;
    }
    if (source_name.Contains("3C264"))
    {
            Source_RA = 176.271;
            Source_Dec = 19.606;
    }
    if (source_name.Contains("OJ287"))
    {
            Source_RA = 133.705;
            Source_Dec = 20.100;
    }
    if (source_name.Contains("PG1553"))
    {
            Source_RA = 238.936;
            Source_Dec = 11.195;
    }
    if (source_name.Contains("Segue1"))
    {
            Source_RA = 151.767;
            Source_Dec = 16.082;
    }
    if (source_name.Contains("1ES1011"))
    {
            Source_RA = 153.767;
            Source_Dec = 49.434;
    }
    if (source_name.Contains("NGC1275"))
    {
            Source_RA = 49.950;
            Source_Dec = 41.512;
    }
    if (source_name.Contains("BLLac"))
    {
            Source_RA = 330.680416667;
            Source_Dec = 42.2777777778;
    }
    if (source_name.Contains("UrsaMajorII"))
    {
            Source_RA = 132.875;
            Source_Dec = 63.13;
    }
    if (source_name.Contains("UrsaMinor"))
    {
            Source_RA = 227.2854167;
            Source_Dec = 67.2225000;
    }
    if (source_name.Contains("RGB_J0710_p591"))
    {
            Source_RA = 107.6100000;
            Source_Dec = 59.1500000;
    }
    if (source_name.Contains("3C273"))
    {
            Source_RA = 187.277915345;
            Source_Dec = 2.05238856846;
    }
    if (source_name.Contains("1ES0502"))
    {
            Source_RA = 76.9839421535;
            Source_Dec = 67.6234172932;
    }
    if (source_name.Contains("Draco"))
    {
            Source_RA = 260.059729167;
            Source_Dec = 57.9212194444;
    }
    if (source_name.Contains("1ES0414"))
    {
            Source_RA = 64.2206666667;
            Source_Dec = 1.089;
    }
    return std::make_pair(Source_RA,Source_Dec);
}

int GetImposterIDFromFile(string source, int onrun_number, int imposter_idx)
{

    // imposter_idx starting from 1
    string line;
    char delimiter = ' ';
    string txt_on_runnumber = "";
    string txt_off_runnumber = "";
    int nth_line = 0;
    int nth_delimiter = 0;
    std::string::size_type sz;

    string SMI_DIR;
    SMI_DIR = string(std::getenv("SMI_DIR"));
    ifstream myfile (RUNLIST_DIR+"/ImposterList_"+source+".txt");
    if (myfile.is_open())
    {
        int current_imposter_idx = 1;
        int current_on_runnumber = 0;
        while ( getline(myfile,line) )
        {
            txt_on_runnumber = "";
            txt_off_runnumber = "";
            nth_delimiter = 0;
            for(int i = 0; i < line.size(); i++)
            {
                if(line[i] == delimiter)
                {
                    nth_delimiter += 1;
                }
                else if (nth_delimiter==1)
                {
                    txt_off_runnumber += line[i];
                }
                else if (nth_delimiter==0)
                {
                    txt_on_runnumber += line[i];
                }
                if (i==line.size()-1)
                {
                    double this_on_runnumber = std::stod(txt_on_runnumber,&sz);
                    double this_off_runnumber = std::stod(txt_off_runnumber,&sz);
                    if (this_on_runnumber!=current_on_runnumber)
                    {
                        current_on_runnumber = this_on_runnumber;
                        current_imposter_idx = 1;
                    }
                    else
                    {
                        current_imposter_idx += 1;
                    }
                    if (onrun_number==int(this_on_runnumber))
                    {
                        if (current_imposter_idx==imposter_idx)
                        {
                            return int(this_off_runnumber);
                        }
                    }
                }
            }
            nth_line += 1;
        }
        myfile.close();
    }
    else cout << "Unable to open PairList file" << endl; 

    return 0;
}

vector<int> GetImposterPairListFromFile(string source, int onrun_number)
{

    string line;
    char delimiter = ' ';
    string txt_on_runnumber = "";
    string txt_off_runnumber = "";
    int nth_line = 0;
    int nth_delimiter = 0;
    std::string::size_type sz;
    vector<int> offrun_list;

    string SMI_DIR;
    SMI_DIR = string(std::getenv("SMI_DIR"));
    ifstream myfile (RUNLIST_DIR+"/ImposterPairList_"+source+".txt");
    if (myfile.is_open())
    {
        int current_on_runnumber = 0;
        while ( getline(myfile,line) )
        {
            txt_on_runnumber = "";
            txt_off_runnumber = "";
            nth_delimiter = 0;
            for(int i = 0; i < line.size(); i++)
            {
                if(line[i] == delimiter)
                {
                    nth_delimiter += 1;
                }
                else if (nth_delimiter==1)
                {
                    txt_off_runnumber += line[i];
                }
                else if (nth_delimiter==0)
                {
                    txt_on_runnumber += line[i];
                }
                if (i==line.size()-1)
                {
                    double this_on_runnumber = std::stod(txt_on_runnumber,&sz);
                    double this_off_runnumber = std::stod(txt_off_runnumber,&sz);
                    if (this_on_runnumber!=current_on_runnumber)
                    {
                        current_on_runnumber = this_on_runnumber;
                    }
                    if (onrun_number==int(this_on_runnumber))
                    {
                        offrun_list.push_back(int(this_off_runnumber));
                    }
                }
            }
            nth_line += 1;
        }
        myfile.close();
    }
    else cout << "Unable to open PairList file" << endl; 

    return offrun_list;
}

vector<int> GetPairListFromFile(string source, int onrun_number)
{

    string line;
    char delimiter = ' ';
    string txt_on_runnumber = "";
    string txt_off_runnumber = "";
    int nth_line = 0;
    int nth_delimiter = 0;
    std::string::size_type sz;
    vector<int> offrun_list;

    string SMI_DIR;
    SMI_DIR = string(std::getenv("SMI_DIR"));
    ifstream myfile (RUNLIST_DIR+"/PairList_"+source+".txt");
    if (myfile.is_open())
    {
        int current_on_runnumber = 0;
        while ( getline(myfile,line) )
        {
            txt_on_runnumber = "";
            txt_off_runnumber = "";
            nth_delimiter = 0;
            for(int i = 0; i < line.size(); i++)
            {
                if(line[i] == delimiter)
                {
                    nth_delimiter += 1;
                }
                else if (nth_delimiter==1)
                {
                    txt_off_runnumber += line[i];
                }
                else if (nth_delimiter==0)
                {
                    txt_on_runnumber += line[i];
                }
                if (i==line.size()-1)
                {
                    double this_on_runnumber = std::stod(txt_on_runnumber,&sz);
                    double this_off_runnumber = std::stod(txt_off_runnumber,&sz);
                    if (this_on_runnumber!=current_on_runnumber)
                    {
                        current_on_runnumber = this_on_runnumber;
                    }
                    if (onrun_number==int(this_on_runnumber))
                    {
                        offrun_list.push_back(int(this_off_runnumber));
                    }
                }
            }
            nth_line += 1;
        }
        myfile.close();
    }
    else cout << "Unable to open PairList file" << endl; 

    return offrun_list;
}

vector<int> GetRunListFromFile(string source)
{

    string line;
    char delimiter = ' ';
    string txt_on_runnumber = "";
    int nth_line = 0;
    std::string::size_type sz;
    vector<int> onrun_list;

    string SMI_DIR;
    SMI_DIR = string(std::getenv("SMI_DIR"));
    ifstream myfile (RUNLIST_DIR+"/RunList_"+source+".txt");
    std::cout << "Read file " << RUNLIST_DIR+"/RunList_"+source+".txt" << std::endl;
    if (myfile.is_open())
    {
        while ( getline(myfile,line) )
        {
            txt_on_runnumber = "";
            for(int i = 0; i < line.size(); i++)
            {
                txt_on_runnumber += line[i];
                if (i==line.size()-1)
                {
                    double this_on_runnumber = std::stod(txt_on_runnumber,&sz);
                    onrun_list.push_back(int(this_on_runnumber));
                }
            }
            nth_line += 1;
        }
        myfile.close();
    }
    else cout << "Unable to open RunList file" << endl; 

    return onrun_list;
}

vector<int> GetOffRunList(string source, int onrun_number) {

    std::cout << "GetPairList source = " << source << std::endl;

    vector<int> list;
    vector<int> list_temp;

    list_temp = GetPairListFromFile(source, onrun_number);
    list.insert(list.end(), list_temp.begin(), list_temp.end());

    return list;

}

vector<int> GetImposterOffRunList(string source, int onrun_number) {

    std::cout << "GetPairList source = " << source << std::endl;

    vector<int> list;
    vector<int> list_temp;

    list_temp = GetImposterPairListFromFile(source, onrun_number);
    list.insert(list.end(), list_temp.begin(), list_temp.end());

    return list;

}

int GetImposterRunID(string source, int onrun_number, int imposter_idx) {

    return GetImposterIDFromFile(source, onrun_number, imposter_idx);

}

vector<int> GetRunList(string source) {

    std::cout << "GetRunList source = " << source << std::endl;

    vector<int> list;
    vector<int> list_temp;

    list_temp = GetRunListFromFile(source);
    list.insert(list.end(), list_temp.begin(), list_temp.end());

    return list;

}
