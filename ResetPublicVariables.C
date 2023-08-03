

void ResetPublicVariables(TString target_name)
{

    string MY_OUTPUT = string(std::getenv("SMI_OUTPUT"));
    std::cout << "MY_OUTPUT = " << MY_OUTPUT << std::endl;

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
