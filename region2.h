/*
 * This code is distributed under the license: MIT License
 * Originally written by Zhenrong WANG
 * mailto: zhenrongwang@live.com
 */

#ifndef REGION2_H_
#define REGION2_H_

#include "if97_general.h"

double gm0(double pi, double tau);
double gm0pi(double pi);
double gm0pipi(double pi);
double gm0tau(double tau);
double gm0tautau(double tau);
double gmr(double pi, double tau);
double gmrpi(double pi, double tau);
double gmrpipi(double pi, double tau);
double gmrtau(double pi, double tau);
double gmrtautau(double pi, double tau);
double gmrpitau(double pi, double tau);
void steam_prop_calc_r2(double pressure, double temp, steam_prop* p_prop);
int sub_region(double pressure, double spec_enth);
int temp_ph_r2(double* temp, double pressure, double spec_enth);
int sub_region_s(double pressure, double spe_entr);
int temp_ps_r2(double* temp, double pressure, double spe_entr);
int sub_region_hs(double spe_enth, double spe_entr);
int pres_hs_r2(double* pressure, double spe_enth, double spe_entr);
void calcFT_dFT_prr2(double* ft, double* dft, double pressure, double dens, double tau);
int temp_pr_r2(double* temp, double pressure, double dens);
void calcFT_dFT_pur2(double* ft, double* dft, double pressure, double spe_ener, double tau);
int temp_pu_r2(double* temp, double pressure, double spe_ener);
void calcFPi_dFPi_trr2(double* fp, double* dfp, double temp, double dens, double pi);
int pres_tr_r2(double* pres, double temp, double dens);
void calcFPi_dFPi_tur2(double* fp, double* dfp, double temp, double spe_ener, double pi);
int pres_tu_r2(double* pres, double temp, double spe_ener);
void calcFPi_dFPi_thr2(double* fp, double* dfp, double temp, double spe_enth, double pi);
int pres_th_r2(double* pres, double temp, double spe_enth);
void calcFPi_dFPi_tsr2(double* fp, double* dfp, double temp, double spe_entr, double pi);
int pres_ts_r2(double* pres, double temp, double spe_entr);
int ptini_hs_r2(double* presini, double* tempini, double spe_enth, double spe_entr);
void calcFG_dFG_hs_r2(double* f, double* g, double*dfdp, double* dfdt, double* dgdp, double* dgdt, double spec_enth, double spec_entr, double pi, double tau);
int pt_hs_r2(double* pres, double* temp, double spe_enth, double spe_entr);

#endif