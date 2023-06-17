
const int N_energy_bins = 11;
double energy_bins[N_energy_bins+1] = {100.,200.,251.,316.,398.,501.,794.,1259.,1995.,3162.,5011.,7943.};
int matrix_rank[N_energy_bins] = {1,2,2,2,2,2,1,1,1,1,1};
int cr_correction_xyoff_bins[N_energy_bins] = {20,20,20,10,10,3,3,1,1,1,1};

double method_ratio_mean[N_energy_bins] =        {1.000,0.009,0.008,0.013,0.014,0.017,0.010,0.004,0.013,0.102,-0.061};
double method_regression_mean[N_energy_bins] =   {1.000,-0.011,-0.002,0.002,0.014,0.015,-0.004,-0.018,-0.080,0.054,0.127};
double method_pertrubation_mean[N_energy_bins] = {1.000,0.002,0.013,0.015,0.022,0.017,0.010,0.004,0.013,0.102,-0.061};
double method_ratio_rms[N_energy_bins] =         {1.000,0.076,0.066,0.063,0.060,0.059,0.074,0.121,0.197,0.399,0.401};
double method_regression_rms[N_energy_bins] =    {1.000,0.053,0.069,0.064,0.078,0.063,0.216,0.155,0.242,0.563,0.539};
double method_pertrubation_rms[N_energy_bins] =  {1.000,0.032,0.043,0.052,0.062,0.059,0.074,0.121,0.197,0.399,0.401};

const int N_Xoff_bins = 1;
const int N_Yoff_bins = 1;

double source_theta_cut = 0.3;
double Elev_cut_lower = 30.;
double Elev_cut_upper = 90.;

double expo_hour_per_analysis = 5.;
//double min_CR_count = 50000.;
//double min_CR_count = 100000.;
double min_CR_count = 200000.;

double MSCW_lower_blind = -0.4;
double MSCL_lower_blind = -0.6;
double MSCW_upper_blind = 0.5;
double MSCL_upper_blind = 0.6;
int n_extra_lower_bins = 1;
int n_extra_upper_bins = 6;
int mtx_dim_w_fine = 6;
int mtx_dim_l_fine = 6;
int mtx_dim_w = 2;
int mtx_dim_l = 2;

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
