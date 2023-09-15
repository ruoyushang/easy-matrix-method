
const int N_energy_bins = 12;
double energy_bins[N_energy_bins+1] = {100.,159.,200.,251.,316.,398.,501.,794.,1259.,1995.,3162.,5011.,7943.};
//int cr_correction_xyoff_bins[N_energy_bins] = {1,12,12,12,12,6,6,1,1,1,1,1};
int cr_correction_xyoff_bins[N_energy_bins] = {1,6,6,6,6,6,6,1,1,1,1,1};

double coefficient_t00xt01[N_energy_bins] = {-0.322,0.040,0.167,0.197,0.200,0.262,0.384,0.469,0.642,0.523,0.776,0.000};
double coefficient_t00xt10[N_energy_bins] = {-0.268,-0.281,-0.343,-0.284,-0.102,0.014,0.062,0.066,0.034,-0.049,-0.005,0.000};
double coefficient_t11xt00[N_energy_bins] = {0.388,0.245,0.183,0.207,0.246,0.307,0.398,0.278,0.235,0.312,0.233,0.000};
double coefficient_t11xt01[N_energy_bins] = {0.287,0.080,0.025,0.036,-0.009,-0.109,-0.287,-0.412,-0.310,-0.261,-0.215,0.000};
double coefficient_t11xt10[N_energy_bins] = {-0.057,-0.157,-0.118,-0.080,-0.084,-0.094,-0.036,-0.004,0.073,0.232,0.445,0.000};

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
