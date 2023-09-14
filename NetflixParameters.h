
const int N_energy_bins = 12;
double energy_bins[N_energy_bins+1] = {100.,159.,200.,251.,316.,398.,501.,794.,1259.,1995.,3162.,5011.,7943.};
int cr_correction_xyoff_bins[N_energy_bins] = {1,12,12,12,12,6,6,1,1,1,1,1};

double coefficient_x00[N_energy_bins] = {0.343,0.185,0.168,0.363,0.216,0.178,0.048,0.071,-0.010,-0.090,-0.557,0.000};
double coefficient_x01[N_energy_bins] = {0.290,0.106,0.036,0.043,0.016,-0.091,-0.103,-0.237,-0.184,-0.542,0.039,0.000};
double coefficient_x10[N_energy_bins] = {-0.059,-0.228,-0.151,0.001,-0.096,-0.078,-0.002,0.042,0.134,-0.052,0.071,0.000};

//const int N_energy_bins = 8;
//double energy_bins[N_energy_bins+1] = {100.,167.,300.,538.,965.,1732.,3107.,5574.,10000.};
//int cr_correction_xyoff_bins[N_energy_bins] = {1,12,12,6,3,3,1,1};

double log_coeff_weight = -1.0;

const int N_Xoff_bins = 1;
const int N_Yoff_bins = 1;

double source_theta_cut = 0.3;
//double Elev_cut_lower = 30.;
double Elev_cut_lower = 40.;
double Elev_cut_upper = 90.;

double expo_hour_per_analysis = 5.;
//double min_CR_count = 5000.;
//double min_CR_count = 50000.;
//double min_CR_count = 100000.;
double min_CR_count = 200000.;
//double min_CR_count = 1000000000.;

int min_NImages = 3;
//double MSCW_lower_blind = -0.7;
//double MSCL_lower_blind = -0.7;
//double MSCW_upper_blind = 0.6;
//double MSCL_upper_blind = 0.6;
double MSCW_lower_blind = -0.5;
double MSCL_lower_blind = -0.7;
double MSCW_upper_blind = 0.7;
double MSCL_upper_blind = 0.5;

int matrix_rank = 2;
int n_extra_lower_bins = 0;
int n_extra_upper_bins = 4;
int mtx_dim_w_fine = 4;
int mtx_dim_l_fine = 4;
int mtx_dim_w = 2;
int mtx_dim_l = 2;

bool use_stat_err_weight = true;
//bool use_stat_err_weight = false;

//double max_Roff = 1.6;
//double max_EmissionHeight_cut = 20.;
//double min_EmissionHeight_cut = 6.;
//double max_Rcore = 250.;
//double min_Rcore = 50.;
//double max_Eerr = 1.0;
double max_Roff = 1.7;
double max_EmissionHeight_cut = 20.;
double min_EmissionHeight_cut = 6.;
double max_Rcore = 400.;
double min_Rcore = 0.;
double max_Eerr = 1.5;

double Skymap_size_x = 2.5;
double Skymap_size_y = 2.5;
int Skymap_nbins_x = 100;
int Skymap_nbins_y = 100;
int ExpoMap_nbins = 20;

double brightness_cut = 6.0;
double bright_star_radius_cut = 0.25;
