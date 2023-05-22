
const int N_energy_bins = 11;
double energy_bins[N_energy_bins+1] = {100.,200.,251.,316.,398.,501.,794.,1259.,1995.,3162.,5011.,7943.};
int matrix_rank[N_energy_bins] = {1,2,2,2,2,2,2,1,1,1,1};

const int N_Xoff_bins = 1;

double source_theta_cut = 0.3;
double Elev_cut_lower = 45.;
double Elev_cut_upper = 90.;

int n_samples_per_analysis = 10;

//int mtx_dim_w = 2;
//int mtx_dim_l = 2;
int mtx_dim_w = 4;
int mtx_dim_l = 4;
int mtx_dim_w_fine = 12;
int mtx_dim_l_fine = 12;
double MSCW_cut_blind = 0.5;
double MSCL_cut_blind = 0.7;
double gamma_hadron_dim_ratio_w = 1.;
double gamma_hadron_dim_ratio_l = 1.;

int min_NImages = 3;
double EmissionHeight_cut = 6.;
double max_Rcore = 250.;

double Skymap_size_x = 2.5;
double Skymap_size_y = 2.5;
int Skymap_nbins_x = 100;
int Skymap_nbins_y = 100;
