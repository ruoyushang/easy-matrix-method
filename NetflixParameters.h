
const int N_energy_bins = 12;
double energy_bins[N_energy_bins+1] = {100.,159.,200.,251.,316.,398.,501.,794.,1259.,1995.,3162.,5011.,7943.};
//int cr_correction_xyoff_bins[N_energy_bins] = {1,12,12,12,12,6,6,1,1,1,1,1};
int cr_correction_xyoff_bins[N_energy_bins] = {1,6,6,6,6,6,6,1,1,1,1,1};

double coefficient_t11xt00_sza[N_energy_bins] = {0.409,0.425,0.287,0.340,0.307,0.482,0.375,0.215,0.211,0.303,0.358,0.000};
double coefficient_t11xt01_sza[N_energy_bins] = {0.191,0.075,-0.018,-0.073,-0.081,-0.231,-0.346,-0.428,-0.406,-0.224,0.071,0.000};
double coefficient_t11xt10_sza[N_energy_bins] = {0.087,-0.112,-0.030,0.102,0.198,0.133,0.110,0.105,0.186,0.421,0.576,0.000};
double coefficient_t11xt00_lza[N_energy_bins] = {0.713,0.471,0.123,0.292,0.265,0.224,0.465,0.389,0.321,0.221,0.239,0.000};
double coefficient_t11xt01_lza[N_energy_bins] = {-0.207,0.083,0.065,0.068,0.030,-0.022,-0.227,-0.404,-0.496,-0.368,-0.338,0.000};
double coefficient_t11xt10_lza[N_energy_bins] = {0.164,-0.165,-0.229,-0.142,-0.116,-0.026,-0.010,-0.012,0.089,0.159,0.339,0.000};

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
//double min_CR_count = 100000.;
double min_CR_count = 200000.;

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
int n_extra_lower_bins = 2;
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
