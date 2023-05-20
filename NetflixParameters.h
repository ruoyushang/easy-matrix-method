
const int N_energy_bins = 11;
double energy_bins[N_energy_bins+1] = {100.,200.,251.,316.,398.,501.,794.,1259.,1995.,3162.,5011.,7943.};

const int N_Xoff_bins = 10;

double source_theta_cut = 0.2;
double Elev_cut_lower = 55.;
double Elev_cut_upper = 90.;

int mtx_dim_w = 4;
int mtx_dim_l = 4;
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
