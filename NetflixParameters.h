
const int N_energy_bins = 12;
double energy_bins[N_energy_bins+1] = {100.,159.,200.,251.,316.,398.,501.,794.,1259.,1995.,3162.,5011.,7943.};
//int cr_correction_xyoff_bins[N_energy_bins] = {1,12,12,12,12,6,6,1,1,1,1,1};
int cr_correction_xyoff_bins[N_energy_bins] = {1,6,6,6,6,6,6,1,1,1,1,1};

double coefficient_t11xt00_lza[N_energy_bins] = {0.313,0.594,0.339,0.393,0.341,0.554,0.466,0.334,0.405,0.095,0.136,0.000};
double coefficient_t11xt01_lza[N_energy_bins] = {0.394,0.224,0.106,0.117,0.052,-0.128,-0.281,-0.453,-0.549,-0.342,-0.256,0.000};
double coefficient_t11xt10_lza[N_energy_bins] = {-0.445,-0.249,-0.239,-0.146,-0.143,-0.119,-0.074,-0.082,0.038,0.183,0.483,0.000};
double coefficient_t11xt00_sza[N_energy_bins] = {0.443,0.324,0.235,0.345,0.410,0.370,0.298,0.311,0.214,0.465,0.040,0.000};
double coefficient_t11xt01_sza[N_energy_bins] = {0.189,0.076,0.038,-0.103,-0.157,-0.277,-0.380,-0.508,-0.260,-0.655,-0.349,0.000};
double coefficient_t11xt10_sza[N_energy_bins] = {0.136,-0.064,-0.096,-0.007,0.082,0.058,0.043,0.078,0.141,0.425,0.008,0.000};

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
