/*
 * This code is distributed under the license: MIT License
 * Originally written by Zhenrong WANG
 * mailto: zhenrongwang@live.com
 */

#ifndef REGION1_H_
#define REGION1_H_

#include "if97_general.h"

double gm(double pi, double tau);
double gmpi(double pi, double tau);
double gmpipi(double pi, double tau);
double gmtau(double pi, double tau);
double gmtautau(double pi, double tau);
double gmpitau(double pi, double tau);
void steam_prop_calc_r1(double pressure, double temp, steam_prop* p_prop);
double tempini_ph_r1(double pressure, double spec_enth);
void calcFT_dFT_ph(double* ft, double* dft, double pres, double spec_enth, double tau);
int temp_ph_r1(double* temp, double pres, double spe_enth);
double tempini_ps_r1(double pressure, double spec_entr);
void calcFT_dFT_ps(double* ft, double* dft, double pres, double spec_entr, double tau);
int temp_ps_r1(double* temp, double pres, double spe_entr);
double presini_hs_r1(double spec_enth, double spec_entr);
double tempini_hs_r1(double spec_enth, double spec_entr);
void calcFG_dFG_hs(double* f, double* g, double*dfdp, double* dfdt, double* dgdp, double* dgdt, double spec_enth, double spec_entr, double pi, double tau);
int pt_hs_r1(double* pres, double* temp, double spe_enth, double spe_entr);
void calcFT_dFT_pr(double* ft, double* dft, double pressure, double dens, double T);
int temp_pr_r1(double* temp, double pressure, double dens);
void calcFT_dFT_pu(double* ft, double* dft, double pressure, double spe_ener, double T);
int temp_pu_r1(double* temp, double pressure, double spe_ener);
void calcFPi_dFPi_tr(double* fpi, double* dfpi, double Temp, double dens, double pi);
int pres_tr_r1(double* pres, double temp, double dens);
void calcFPi_dFPi_tu(double* fpi, double* dfpi, double temp, double spe_ener, double pi);
int pres_tu_r1(double* pres, double temp, double spe_ener);
void calcFPi_dFPi_th(double* fpi, double* dfpi, double temp, double spe_enth, double pi);
int pres_th_r1(double* pres, double temp, double spe_enth);
void calcFPi_dFPi_ts(double* fpi, double* dfpi, double temp, double spe_entr, double pi);
int pres_ts_r1(double* pres, double temp, double spe_entr);

#endif