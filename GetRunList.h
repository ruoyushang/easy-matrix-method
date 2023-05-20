
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include "NetflixParameters.h"

pair<double,double> GetSourceRaDec(TString source_name)
{
    double Source_RA = 0.;
    double Source_Dec = 0.;
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
    return std::make_pair(Source_RA,Source_Dec);
}

vector<int> GetPairListFromFile(string source, int onrun_number, int imposter_idx)
{

    // imposter_idx starting from 0
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
    ifstream myfile (SMI_DIR+"/runlist/PairList_"+source+".txt");
    if (myfile.is_open())
    {
        int current_imposter_idx = 0;
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
                        current_imposter_idx = 0;
                    }
                    else
                    {
                        current_imposter_idx += 1;
                    }
                    if (onrun_number==int(this_on_runnumber))
                    {
                        if (imposter_idx==-1)
                        {
                            offrun_list.push_back(int(this_off_runnumber));
                        }
                        else if (current_imposter_idx==imposter_idx)
                        {
                            offrun_list.push_back(int(this_off_runnumber));
                        }
                    }
                }
            }
            nth_line += 1;
        }
        myfile.close();
    }
    else cout << "Unable to open file"; 

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
    ifstream myfile (SMI_DIR+"/runlist/RunList_"+source+".txt");
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
    else cout << "Unable to open file"; 

    return onrun_list;
}

vector<int> GetOffRunList(string source, int onrun_number, int imposter_idx) {

    // imposter_idx starting from 0
    std::cout << "GetPairList source = " << source << std::endl;

    vector<int> list;
    vector<int> list_temp;

    list_temp = GetPairListFromFile(source, onrun_number,imposter_idx);
    list.insert(list.end(), list_temp.begin(), list_temp.end());

    return list;

}

vector<int> GetRunList(string source) {

    std::cout << "GetRunList source = " << source << std::endl;

    vector<int> list;
    vector<int> list_temp;

    list_temp = GetRunListFromFile(source);
    list.insert(list.end(), list_temp.begin(), list_temp.end());

    return list;

}
