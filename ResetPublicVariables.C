
//void ResetCoefficients(double coefficient_t11xt00[], double coefficient_t11xt01[], double coefficient_t11xt10[])
//{
//    for (int ebin=0;ebin<N_energy_bins;ebin++)
//    {
//        coefficient_t11xt00_incl[ebin] = coefficient_t11xt00[ebin];
//        coefficient_t11xt01_incl[ebin] = coefficient_t11xt01[ebin];
//        coefficient_t11xt10_incl[ebin] = coefficient_t11xt10[ebin];
//    }
//}

void ResetPublicVariables(TString source_name)
{

    string MY_OUTPUT = string(std::getenv("SMI_OUTPUT"));
    std::cout << "MY_OUTPUT = " << MY_OUTPUT << std::endl;

}
