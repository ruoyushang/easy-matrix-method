
const int N_energy_bins = 11;
double energy_bins[N_energy_bins+1] = {100.,200.,251.,316.,398.,501.,794.,1259.,1995.,3162.,5011.,7943.};
int matrix_rank[N_energy_bins] = {1,2,2,2,2,2,1,1,1,1,1};
//int matrix_rank[N_energy_bins] = {1,3,3,3,3,3,1,1,1,1,1};
int cr_correction_xyoff_bins[N_energy_bins] = {20,20,20,10,10,3,3,1,1,1,1};

double method_ratio_mean[N_energy_bins] =        {1.000,0.000,0.007,0.010,0.007,0.011,0.012,0.007,0.001,0.083,0.005};
double method_pertrubation_mean[N_energy_bins] = {1.000,-0.003,0.007,0.008,0.007,0.008,0.012,0.007,0.001,0.083,0.005};
double method_regression_mean[N_energy_bins] =   {1.000,-0.014,0.002,-0.001,0.006,0.011,-0.001,-0.029,-0.101,0.012,0.238};
double method_ratio_rms[N_energy_bins] =         {1.000,0.063,0.064,0.054,0.051,0.050,0.071,0.104,0.162,0.321,0.347};
double method_pertrubation_rms[N_energy_bins] =  {1.000,0.027,0.042,0.045,0.047,0.049,0.071,0.104,0.162,0.321,0.347};
double method_regression_rms[N_energy_bins] =    {1.000,0.040,0.080,0.074,0.050,0.056,0.090,0.108,0.197,0.288,0.567};

double log_coeff_weight = -1.;

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
int n_extra_lower_bins = 2;
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
