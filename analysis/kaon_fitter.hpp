// Methods to interface K → 3π decay data with fitters
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef KAON_FITTER_HPP
#define KAON_FITTER_HPP

#include "amplitudes/kaon.hpp"
#include "isobars/pseudoscalar.hpp"

namespace iterateKT { namespace kaon
{
    
    // Specify the fitter interface
    struct fit
    {
        // Static identifiers for data_set types
        static const int kAll = 0, kWidth = 1, kDalitz = 2;

        // The offset we use for the derivatives of dpars
        static constexpr double derivative_h = 1E-6;

        // Whether or not we are going to use the full set of subtractions allowed by 
        // the froissart bound asymptotics
        static const bool FULL_SUBTRACTIONS = true;

        // String letting us know what is being fit
        static std::string data_type(int i)
        {
            switch (i)
            {
                case kAll:     return "Width & {g, h, k}";
                case kWidth:   return "Width";
                case kDalitz:  return "{g, h, k}";
                default:       return "ERROR!";
            };
        };

        // Function being minimized, sum of chi2s of individual data sets and observables
        static double fcn(std::vector<data_set> & data_vector, amplitude to_fit)
        {
            double chi2_tot = 0;
            for (auto data : data_vector) chi2_tot += chi2(data, to_fit);
            return chi2_tot;
        };

        // Indidivudal chi2 from a single data set
        static double chi2(const data_set & data, amplitude to_fit)
        {
            int type = data._type;
            bool n = (data._type == kAll); // Whether to skip the first slot (width) in the data
            to_fit->set_option(data._option);

            double chi2 = 0;
            // χ² from Γ        
            if (type == kAll || type == kWidth)
            {
                double gam_th = to_fit->width();
                double gam_ex = data._z[0], dgam_ex = data._dz[0];
                chi2 += norm((gam_th - gam_ex)/dgam_ex);
            };
            // χ² from g h k
            if (type == kAll || type == kDalitz)
            {
                double s0 = (data._option == option::P_ppm) ? (norm(M_KAON_PM) + 3*norm(M_PION_PM))/3
                                                            : (norm(M_KAON_PM) + 2*norm(M_PION_0) + norm(M_PION_PM))/3;
                auto dpars  = to_fit->get_dalitz_parameters(derivative_h, s0, {M_PION_PM*M_PION_PM, M_PION_PM*M_PION_PM});
                std::array<double,3> ghk = {dpars[0], dpars[1], dpars[3]};
                for (int i = 0; i < 3; i++) chi2 += norm((ghk[i]-data._z[i+n])/data._dz[i+n]);
            };
            return chi2;
        };

        // We only fit the real parts of the parameters while the imaginary parts are
        // given by requiring Taylor invariants have vanishing imaginary parts
        static std::vector<complex> process_fitter_parameters(std::vector<complex> in_pars, amplitude amp)
        {
            double eps = 1E-7, r = amp->get_kinematics()->s0();
            
            //------------------------------------------------------------------------
            // First we fix the imaginary parts of the M's and N's (total 3π I=1)
            isobar F0 = amp->get_isobar(id::I1_S0);
            isobar F1 = amp->get_isobar(id::I1_P1);
            isobar F2 = amp->get_isobar(id::I1_S2);
            std::array<isobar,3> F = {F0, F1, F2};
            
            // Coefficients of taylor expansion up to cubic
            std::array<std::vector<complex>,3> A, B, C;
            
            uint dim = 3 + FULL_SUBTRACTIONS;
            // First index is isospin, second is basis function ID
            for (uint i = 0; i < 3; i++)
            {
                for (uint n = 0; n < dim; n++)
                {
                    A[i].push_back(F[i]->basis_function(n, 0));
                    B[i].push_back(F[i]->basis_derivative<1>(n, 0, eps));
                    C[i].push_back(F[i]->basis_derivative<2>(n, 0, eps)/2.);
                };
            };

            // Construct the nxn matrix of Taylor invariants
            // First index is which invariant, and second is basis_id
            TArrayD reT_data(dim*dim), imT_data(dim*dim);
            for (uint n = 0; n < dim; n++)
            {
                std::array<complex,4> T;
                T[0] =   A[0][n] + r*B[0][n] + 4*(A[2][n]+r*B[2][n])/3;
                T[1] = 3*A[1][n] +   B[0][n] - 5*B[2][n]/3 + 9*r*(B[1][n] + 2*r*C[1][n]);
                T[2] =   C[2][n] +   B[1][n] + 3*r*C[1][n];
                T[3] = 3*C[0][n] + 4*C[2][n];
               
                for (int j = 0; j < dim; j++)
                {
                    reT_data[dim*j+n] = -real(T[j]);  imT_data[dim*j+n] = imag(T[j]);
                };
            };
            TMatrixD reT(dim,dim,reT_data.GetArray()), imT(dim,dim,imT_data.GetArray());

            // Now we actually solve the matrix equation relating reMu and imMu
            TMatrixD M = reT.Invert()*imT;

            std::vector<Double_t> rePars_data(dim);
            for (int i = 0; i < dim; i++) rePars_data[i] = real(in_pars[i]);
            TVectorD rePars(dim, rePars_data.data()), imPars = M*rePars;
            
            // Assemble together output vector
            std::vector<complex> out_pars;
            for (int i = 0; i < dim; i++) out_pars.push_back(rePars[i]+I*imPars[i]);

            if (!FULL_SUBTRACTIONS) 
            { out_pars.push_back(in_pars.back()); return out_pars; };

            //------------------------------------------------------------------------
            // Now do the same for the H's 

            // These inhabit basis functions 9 & 10
            std::array<isobar,3> H;
            H[0] = nullptr;
            H[1] = amp->get_isobar(id::I2_P1); H[2] = amp->get_isobar(id::I2_S2); 
            for (uint i = 1; i <= 2; i++)
            {
                for (uint n = 0; n <= 1; n++)
                {
                    A[i][n] = H[i]->basis_function     (n+4, 0);
                    B[i][n] = H[i]->basis_derivative<1>(n+4, 0, eps);
                    C[i][n] = H[i]->basis_derivative<2>(n+4, 0, eps)/2.;
                };
            };

            // Construct the 2x2 matrix of Taylor invariants
            // First index is which invariant, and second is basis_id
            TArrayD reTp_data(2*2), imTp_data(2*2);
            for (uint n = 0; n <= 1; n++)
            {
                // See second rows of Eq. 6.5 in [1]
                std::array<complex,2> Tp;
                Tp[0] = 3*A[1][n] - B[2][n] + 9*r*(B[1][n] + 2*C[1][n]);
                Tp[1] = 3*B[1][n] + C[2][n] + 9*r* C[1][n];
                               
                for (int j = 0; j <= 1; j++)
                {
                    reTp_data[2*j+n] = -real(Tp[j]);  imTp_data[2*j+n] = imag(Tp[j]);
                };
            };
            TMatrixD reTp(2,2, reTp_data.GetArray()), imTp(2,2, imTp_data.GetArray());
            TMatrixD Mp = reTp.Invert()*imTp;

            Double_t reNup_data[2];
            for (int i = 0; i <= 1; i++) reNup_data[i] = real(in_pars[i+4]);
            TVectorD reNup(2, reNup_data), imNup = Mp*reNup;

            out_pars.push_back(reNup[0] + I*imNup[0]);
            out_pars.push_back(reNup[1] + I*imNup[1]);
            return out_pars;
        };
    };
}; /* namespace iterateKT */ }; /* namespace kaon_decay */
#endif