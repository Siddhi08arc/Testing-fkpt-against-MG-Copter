#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <cmath>
#include <algorithm>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <iomanip>
#include <Copter/Cosmology.h>
#include <Copter/LinearPS.h>
#include <Copter/SPT.h>
#include <Copter/GrowthFunction.h>
#include <Copter/SpecialFunctions.h>
#include <Copter/NoWigglePS.h>
#include <Copter/BeyondLCDM.h>

using namespace std;

std::string model = "FR5";
std::string types = "HALO";

vector<double> compute_chi_cpp(
    const string eft_file,           // EFTofLSS_FR5_z1_HALO.dat
    const string ctrl_file,          // CTRL_FR5_Z1_HALO.dat
    const char* invcov_file,         // multidata_invcov.dat
    const double kmax                // kmax
){
    auto loadFileSkip37 = [](const string& file, int ncols) {
        ifstream fin(file);
        if (!fin.is_open()) {
            cerr << "Error: Cannot open file " << file << endl;
            exit(1);
        }
        vector<vector<double>> cols(ncols);
        string line;
        int line_num = 0;

        while (getline(fin, line)) {
            ++line_num;
            istringstream ss(line);
            for (int c = 0; c < ncols; c++) {
                double x;
                ss >> x;
                cols[c].push_back(x);
            }
        }
        return cols;
    };

    // Loading EFTofLSS file (10 columns)
    // k, p1, p2, p3, p4, p5, p6, p7, p8, p9

    auto eft = loadFileSkip37(eft_file, 10);

    vector<double>& k  = eft[0];
    vector<double>& p1 = eft[1];  // EFTofLSS exact P0 part
    vector<double>& p2 = eft[2];  // EFTofLSS exact P2 part
    vector<double>& p3 = eft[3];  // EFTofLSS exact P4 part
    vector<double>& p4 = eft[4];  // EFTofLSS counterterms P0
    vector<double>& p5 = eft[5];  // EFTofLSS counterterms P2
    vector<double>& p6 = eft[6];  // EFTofLSS counterterms P4
    vector<double>& p7 = eft[7];  // model P0 part
    vector<double>& p8 = eft[8];  // model P2 part
    vector<double>& p9 = eft[9];  // model P4 part

    int N = k.size();

    // Loading CTRL file (4 columns)
    // k, ctrlp0, ctrlp2, ctrlp4

    auto ctrl = loadFileSkip37(ctrl_file, 4);

    vector<double>& ctrlp0 = ctrl[1];
    vector<double>& ctrlp2 = ctrl[2];
    vector<double>& ctrlp4 = ctrl[3];

    // DATA multipoles (EXACT)

    vector<double> P0_data(N), P2_data(N), P4_data(N);
    for(int i = 0; i < N; i++) {
        P0_data[i] = p1[i] + p4[i];
        P2_data[i] = p2[i] + p5[i];
        P4_data[i] = p3[i] + p6[i];
    }

    // MODEL multipoles

    vector<double> P0_model(N), P2_model(N), P4_model(N);
    for (int i = 0; i < N; i++) {
        P0_model[i] = p7[i] + ctrlp0[i];
        P2_model[i] = p8[i] + ctrlp2[i];
        P4_model[i] = p9[i] + ctrlp4[i];
    }

    // Loading inverse covariance matrix
    // Format: k  iCov00  iCov02  iCov04  iCov22  iCov24  iCov44

    ifstream fin_cov(invcov_file);
    if (!fin_cov.is_open()) {
        cerr << "Error: Cannot open file multidata_invcov.dat" << endl;
        return {};
    }

    // iCov[bin][row][col] - 3x3 matrix per k-bin

    vector<vector<vector<double>>> iCov(N, vector<vector<double>>(3, vector<double>(3, 0.0)));

    for (int i = 0; i < N; i++) {
        double kval, c00, c02, c04, c22, c24, c44;
        fin_cov >> kval >> c00 >> c02 >> c04 >> c22 >> c24 >> c44;

        iCov[i][0][0] = c00;
        iCov[i][0][1] = c02;
        iCov[i][0][2] = c04;
        iCov[i][1][1] = c22;
        iCov[i][1][2] = c24;
        iCov[i][2][2] = c44;
        // Fill symmetric entries
        iCov[i][1][0] = iCov[i][0][1];
        iCov[i][2][0] = iCov[i][0][2];
        iCov[i][2][1] = iCov[i][1][2];
    }
    fin_cov.close();

    // Map k values to indices
   
    map<double, int> kbins;
    for (int i = 0; i < N; i++) {
        kbins[k[i]] = i;
    }

    double kmin = 0.01;
    vector<double> chi_append;
        vector<int> indices;
        for (auto& kv : kbins) {
            if (kv.first >= kmin && kv.first <= kmax) {
                indices.push_back(kv.second);
            }
        }

        int Nbins = indices.size();
        int Nparams = 3;
        int Ndof = 2 * Nbins - Nparams;  // 2 multipoles per bin

        double chi2_total = 0.0;

        for (int i : indices) {

            // Residuals: [P0_data - P0_model, P2_data - P2_model]
            double resid0 = P0_data[i] - P0_model[i];
            double resid2 = P2_data[i] - P2_model[i];

            // Extract 2x2 submatrix: rows/cols 0 and 2 (P0 and P4 indices in original)
            // But looking at Python code: it uses iCov[i, 0, 0], iCov[i, 0, 2], iCov[i, 2, 0], iCov[i, 2, 2] (ignore, this was for my reference)
            // This corresponds to P0-P0, P0-P4, P4-P0, P4-P4 covariance
            double Cinv00 = iCov[i][0][0];
            double Cinv01 = iCov[i][0][2];
            double Cinv10 = iCov[i][2][0];
            double Cinv11 = iCov[i][2][2];

            // chi2 = resid.T @ Cinv @ resid
            // = [r0, r2] @ [[c00, c01], [c10, c11]] @ [r0, r2]^T
            double t0 = Cinv00 * resid0 + Cinv01 * resid2;
            double t1 = Cinv10 * resid0 + Cinv11 * resid2;
            chi2_total += resid0 * t0 + resid2 * t1;
        }

        chi_append.push_back(chi2_total / Ndof);
    // }

    return chi_append;
}

void compute_spectrum(double rsdpars[], 
                    Cosmology &C, LinearPS &P_l, SPT &spt, IOW &iow,
                    double pars[], double extpars[], double bias[],
                    int mymodel) 
{

    std::string temp = "Final_CTRL_" + model + "_Z1_" + types + ".dat";
    const char* output = temp.c_str();
    int Nk = 30;
    double kmin = 0.01;
    double kmax = 0.3;
    double err = 1e-3;

    FILE* fp = fopen(output, "w");

    for(int i = 0; i < Nk; i++) {

        double k = kmin * exp(i * log(kmax / kmin) / (Nk - 1));
        double p4 = spt.PRSD_mg(5, 1, pars, extpars, rsdpars, bias, k, err, mymodel);
        double p5 = spt.PRSD_mg(5, 2, pars, extpars, rsdpars, bias, k, err, mymodel);
        double p6 = spt.PRSD_mg(5, 3, pars, extpars, rsdpars, bias, k, err, mymodel);

        fprintf(fp, "%e %e %e %e\n", k, p4, p5, p6);
    }

    fclose(fp);
}

struct Range {
    double start;
    double end;
    double step;
};

//ALL THE CHANGES FOR EACH MODEL

int main(int argc, char* argv[]) {
    int mymodel = 2;                                          // 1 for GR, 2 for F4, 2 for F5, 3 for DGP
    double myz = 1;                                           // Unchanged as I am working at z=1
    double omega0 = 0.3132; 
    if (model == "DGP") {
        mymodel = 3;
    }
    double mgpar = 1e-4;                                     //0 for GR, 1e-4 for F4, 1e-5 for F5, 1 for DGP
    if (model == "FR5") {
        mgpar = 1e-5;
    }
    else if (model == "FR4") {
        mgpar = 1e-4;
    }

    double rsd_model = 4;

    // Lagrangian bias params (See 1607.03149 for example) OR Q-bias params (see PRSD_mg in src/SPT.cpp).
    double bias[8];

    if (types == "HALO") {
        // Halos 
      bias[0] = 1.5;              // b1      
      bias[1] = 0.7;              // b2
      bias[2] = 0.3;              // bG2
      bias[3] = -0.7;             // bG3
      // Shot noise terms
      bias[4] = 1/5e-3;           // N
      bias[5] = 5.;               // epsilon_0
      bias[6] = 10.;              // epsilon_2
      bias[7] = 1/bias[4];        // nbar
    } else if (types == "CDM") {
        // CDM values 
      bias[0] = 1.;               // b1
      bias[1] = 0.;               // b2
      bias[2] = 0.;               // bG2
      bias[3] = 0.;               // bG3
      bias[4] = 0.;               // N
      bias[5] = 0.;               // epsilon_0
      bias[6] = 0.;               // epsilon_2
      bias[7] = 0.001;            // nbar
    }

    double err = 1e-3;

    const char* cstr = "transfers/planck";

    real z = 0;
    real epsrel = 1e-3;

    Cosmology C(cstr);
    LinearPS P_l(C, z);
    SPT spt(C, P_l, epsrel);
    IOW iow;

    double pars[10];
    pars[0] = 1./(1.+myz);
    pars[1] = omega0;
    pars[2] = 0.;
    pars[7] = 1e-4;
    pars[8] = 1e-3;
    pars[9] = 1e-3;

    double extpars[maxpars];
    extpars[0] = mgpar;

    double rsdpars[4];

    iow.initnorm(pars, extpars, mymodel);
    rsdpars[0] = spt.sigmav_init(pars, extpars, mymodel);

    rsdpars[1] = 6.3;         //blue
    rsdpars[2] = 25.7;        //yellow
    rsdpars[3] = 39.2;        //green

    iow.initn_lin(pars, extpars, 0.00001, mymodel);
    fl_spt = -G1_nk/F1_nk;
    Dl_spt = F1_nk;

    double sigveds = spt.sigmav_init_linear();
    printf("%e %e \n", rsdpars[0], sigveds);

    // Acceptable χ² range and convergence condition
    auto chi_acceptable = [](double x) {
        return (x >= 0.0 && x <= 4.0);
    };

    auto chi_close_to_one = [](double x) {
        return std::fabs(x - 1.0) < 0.10;   
    };

    // Each counter list is a vector of Range structs {start, end, step}
    // i
    std::vector<Range> counter_0 = {
        {4, 10, 0.8},{5, 7, 0.5}, {5, 7, 0.5}, {6, 7, 0.1}, {5, 8, 0.8}, {6, 10, 0.1}, {6, 7, 0.1}
    };
    // j
    std::vector<Range> counter_1 = {
        {10, 30, 2}, {10, 21, 0.5}, {10, 21, 0.5}, {10, 40, 0.8}, {10, 35, 0.5}, {10, 40, 0.7}, {20, 23, 0.1}
    };
    // k
    std::vector<Range> counter_2 = {
        {10, 70, 2}, {10, 41, 1}, {10, 56, 1}, {10, 56, 1}, {10, 50, 1}, {10, 60, 1}, {30, 40, 0.5}
    };

    // kmax list
    std::vector<double> kmax_list = {1.044009e-01, 1.173923e-01, 1.320002e-01, 1.484260e-01, 1.668957e-01, 1.876638e-01, 2.110161e-01};

    std::string filename = "Chi_values_" + model + types + "_EFTofLSS.dat";

    std::ofstream file(filename, std::ios::app);
    std::string eftfile = "Final_" + model + "_z1_" + types + "_EFTofLSS.dat";
    std::string ctrlfile = "Final_CTRL_" + model + "_Z1_" + types + ".dat";

    // Loop using index (ki)

    for (size_t ki = 5; ki < kmax_list.size(); ki++) {

        double kmax = kmax_list[ki];

        for (double i = counter_0[ki].start; i <= counter_0[ki].end; i += counter_0[ki].step) {
            for (double j = counter_1[ki].start; j <= counter_1[ki].end; j += counter_1[ki].step) {
                for (double k = counter_2[ki].start; k <= counter_2[ki].end; k += counter_2[ki].step) {

    for (double kmax : kmax_list) {
        for (float i=2; i <= 12; i+=1) {
            for (float j=12; j <=40; j+=1) {
                for (float k=15; k <= 50; k+=1) {
                    rsdpars[1] = i;
                    rsdpars[2] = j;
                    rsdpars[3] = k;

                    // Producing spectrum
                    compute_spectrum(rsdpars, C, P_l, spt, iow, pars, extpars, bias, mymodel);
                    std::vector<double> chi = compute_chi_cpp(
                        eftfile,       // Data file
                        ctrlfile,           // Control file from above
                        "multidata_invcov.dat",
                        kmax
                    );

                    // Uncomment to see the chi2 values being processed for each nuisance parameters

                    std::cout << "\nIteration "
                    << "  rsdpars = "
                    << rsdpars[1] << ", "
                    << rsdpars[2] << ", "
                    << rsdpars[3] << "\n";
                    std::cout << "chi² values = ";
                    std::cout << chi[0] << " ";
                    std::cout << chi[1] << " ";
                    std::cout << chi[2] << " ";
                    std::cout << chi[3] << " ";
                    std::cout << chi[4] << " ";
                    std::cout << chi[5] << " ";
                    std::cout << chi[6] << " ";
                    std::cout << chi[7] << " ";
                    std::cout << chi[8] << " ";
                    std::cout << chi[9] << " ";
                    std::cout << "\n";
                    if ((chi[0] >= -1 && chi[1] >= -1 && chi[2] >= -1 && chi[3] >= -1 && chi[4] >= -1 && chi[5] >= -1 && chi[6] >= -1 && chi[7] >= -1 && chi[8] >= -1 ) 
                        && (chi[0] <= 1 && chi[1] <= 1 && chi[2] <= 1 && chi[3] <= 3 && chi[4] <= 3 && chi[5] <= 10 && chi[6] <= 100 && chi[7] <= 500 && chi[8] <= 1000 ) ) {
                    if ((chi[0] >= -1 && chi[1] >= -1 && chi[2] >= -1 && chi[3] >= -1 && chi[4] >= -1 && chi[5] >= -1 && chi[6] >= -1 && chi[7] >= -1 && chi[8] >= -1 ) 
                        && (chi[0] <= 1 && chi[1] <= 1 && chi[2] <= 1 && chi[3] <= 1 && chi[4] <= 1 && chi[5] <= 1.5 && chi[6] <= 8 && chi[7] <= 52 && chi[8] <= 192 ) ) {
                    if ((chi[0] >= -1 && chi[1] >= -1 && chi[2]) 
                        && (chi[0] <= 1.2 && chi[1] <= 3 && chi[2] <= 3) ) {
                            && chi[1] >= -1 && chi[2] >= -1 && chi[3] >= -1 && chi[4] >= -1
                             && chi[1] <= 3 && chi[2] <= 3 && chi[3] <= 100 && chi[4] <= 30000
                        std::cout << "caught one!!!!!!";
                        std::cout << "\n";
                    }
                    if ((chi[0] <= 20) ) {
                        file << " " << model << " " << types << " counters: ";
                        file << i << " " << j << " " << k << " chi: ";
                        for (double x : chi) file << x << " ";
                        file << " kmax: " << kmax;
                        file << "\n";
                    }
                    // Uncomment to generate graph at the same time
                    // std::system("python3 errorgraph.py");
                }
            }
        }
    
    }
    file.close();

    return 0;
}
