/*
 * This code is distributed under the license: MIT License
 * Originally written by Zhenrong WANG
 * mailto: zhenrongwang@live.com
 */

#ifndef REGION5_H_
#define REGION5_H_

#include "if97_general.h"

double gm05(double pi, double tau);
double gmr5(double pi, double tau);
double gm0pi5(double pi);
double gm0pipi5(double pi);
double gm0tau5(double tau);
double gm0tautau5(double tau);
double gmrpi5(double pi, double tau);
double gmrpipi5(double pi, double tau);
double gmrtau5(double pi, double tau);
double gmrtautau5(double pi, double tau);
double gmrpitau5(double pi, double tau);
void steam_prop_calc_r5(double pressure, double temp, steam_prop* p_prop);
void calcFT_dFT_prr5(double* ft, double* dft, double pressure, double dens, double tau);
int temp_pr_r5(double *temp, double pressure, double dens);
void calcFT_dFT_pur5(double* ft, double* dft, double pressure, double spe_ener, double tau);
int temp_pu_r5(double* temp, double pressure, double spe_ener);
void calcFT_dFT_phr5(double* ft, double* dft, double pressure, double spe_enth, double tau);
int temp_ph_r5(double* temp, double pressure, double spe_enth);
void calcFT_dFT_psr5(double* ft, double* dft, double pressure, double spe_entr, double tau);
int temp_ps_r5(double* temp, double pressure, double spe_entr);
void calcFP_dFP_trr5(double* fp, double* dfp, double temp, double dens, double pi);
int pres_tr_r5(double* pres, double temp, double dens);
void calcFP_dFP_tur5(double* fp, double* dfp, double temp, double spe_ener, double pi);
int pres_tu_r5(double* pres, double temp, double spe_ener);
void calcFP_dFP_thr5(double* fp, double* dfp, double temp, double spe_enth, double pi);
int pres_th_r5(double* pres, double temp, double spe_enth);
void calcFP_dFP_tsr5(double* fp, double* dfp, double temp, double spe_entr, double pi);
int pres_ts_r5(double* pres, double temp, double spe_entr);
void calcFG_dFG_hs_r5(double* f, double* g, double*dfdp, double* dfdt, double* dgdp, double* dgdt, double spec_enth, double spec_entr, double pi, double tau);
int pt_hs_r5(double* pres, double* temp, double spe_enth, double spe_entr);

#endif