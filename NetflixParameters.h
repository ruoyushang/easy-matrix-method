
const int N_energy_bins = 11;
double energy_bins[N_energy_bins+1] = {100.,200.,251.,316.,398.,501.,794.,1259.,1995.,3162.,5011.,7943.};
int matrix_rank[N_energy_bins] = {1,2,2,2,2,2,2,1,1,1,1};
int cr_correction_xyoff_bins[N_energy_bins] = {20,20,20,10,10,3,3,1,1,1,1};

double method_ratio_mean[N_energy_bins] =        {1.000,0.009,0.009,0.013,0.014,0.016,0.009,0.003,0.011,0.100,-0.061};
double method_regression_mean[N_energy_bins] =   {1.000,-0.011,-0.002,0.002,0.011,0.015,-0.003,-0.019,-0.083,0.052,0.127};
double method_pertrubation_mean[N_energy_bins] = {1.000,0.007,0.015,0.020,0.030,0.031,0.041,0.003,0.011,0.100,-0.061};
double method_ratio_rms[N_energy_bins] =         {1.000,0.076,0.066,0.063,0.059,0.059,0.074,0.122,0.196,0.399,0.402};
double method_regression_rms[N_energy_bins] =    {1.000,0.053,0.069,0.065,0.061,0.062,0.214,0.156,0.241,0.564,0.540};
double method_pertrubation_rms[N_energy_bins] =  {1.000,0.037,0.043,0.052,0.059,0.061,0.096,0.122,0.196,0.399,0.402};

const int N_Xoff_bins = 1;
const int N_Yoff_bins = 1;

double source_theta_cut = 0.3;
double Elev_cut_lower = 30.;
double Elev_cut_upper = 90.;

double expo_hour_per_analysis = 5.;
//double min_CR_count = 50000.;
double min_CR_count = 100000.;
//double min_CR_count = 200000.;

double MSCW_cut_blind = 0.5;
double MSCL_cut_blind = 0.7;
int n_extra_lower_bins = 1;
int mtx_dim_w = 2;
int mtx_dim_l = 2;
int mtx_dim_w_fine = 12;
int mtx_dim_l_fine = 12;
double gamma_hadron_dim_ratio_w = 1.;
double gamma_hadron_dim_ratio_l = 1.;
//int mtx_dim_w = 3;
//int mtx_dim_l = 3;
//int mtx_dim_w_fine = 12;
//int mtx_dim_l_fine = 12;
//double gamma_hadron_dim_ratio_w = 2.;
//double gamma_hadron_dim_ratio_l = 2.;

bool use_stat_err_weight = true;
//bool use_stat_err_weight = false;

int min_NImages = 3;
double max_Roff = 1.6;
double max_EmissionHeight_cut = 20.;
double min_EmissionHeight_cut = 6.;
double max_Rcore = 250.;
double min_Rcore = 50.;

double Skymap_size_x = 2.5;
double Skymap_size_y = 2.5;
int Skymap_nbins_x = 100;
int Skymap_nbins_y = 100;

double brightness_cut = 6.0;
double bright_star_radius_cut = 0.25;
