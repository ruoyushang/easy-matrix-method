

void ResetPublicVariables(TString target_name)
{

    string MY_OUTPUT = string(std::getenv("SMI_OUTPUT"));
    std::cout << "MY_OUTPUT = " << MY_OUTPUT << std::endl;

    for (int e=0;e<N_energy_bins;e++) 
    {
        method_ratio_mean[e] = 0.;
        method_regression_mean[e] = 0.;
        method_pertrubation_mean[e] = 0.;
        method_ratio_rms[e] = 1.;
        method_regression_rms[e] = 1.;
        method_pertrubation_rms[e] = 1.;
    }

    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_weight_log_m0p0")
    {
        log_coeff_weight = 0.0;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_weight_log_m0p5")
    {
        log_coeff_weight = -0.5;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_weight_log_m1p0")
    {
        log_coeff_weight = -1.;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_weight_log_m1p5")
    {
        log_coeff_weight = -1.5;
    }
    if (MY_OUTPUT=="/gamma_raid/userspace/rshang/SMI_output/output_weight_log_m2p0")
    {
        log_coeff_weight = -2.;
    }

}
