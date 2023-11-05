
const int N_energy_bins = 12;
double energy_bins[N_energy_bins+1] = {100.,159.,200.,251.,316.,398.,501.,794.,1259.,1995.,3162.,5011.,7943.};
//int cr_correction_xyoff_bins[N_energy_bins] = {1,12,12,12,12,6,6,1,1,1,1,1};
int cr_correction_xyoff_bins[N_energy_bins] = {1,6,6,6,6,6,6,1,1,1,1,1};

double coefficient_t11xt01_sza[N_energy_bins] = {-2.3e-02,2.6e-02,-1.9e-02,-5.0e-02,-2.8e-02,-1.0e-01,-3.4e-01,-5.5e-01,0.0e+00,0.0e+00,0.0e+00,0.0e+00};
double coefficient_t11xt10_sza[N_energy_bins] = {-1.4e-02,-2.6e-01,-2.0e-01,8.3e-02,2.3e-01,1.3e-01,5.4e-02,-1.0e-01,0.0e+00,0.0e+00,0.0e+00,0.0e+00};
double coefficient_t11xt01_lza[N_energy_bins] = {-2.5e-01,-5.5e-02,-5.4e-02,-4.9e-02,-4.9e-02,-9.5e-02,-2.0e-01,-3.3e-01,-4.8e-01,0.0e+00,0.0e+00,0.0e+00};
double coefficient_t11xt10_lza[N_energy_bins] = {-6.3e-01,-3.2e-01,-2.9e-01,-1.3e-01,-3.6e-02,-6.5e-03,6.4e-02,6.9e-02,-6.9e-02,0.0e+00,0.0e+00,0.0e+00};

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
