/*
 * This code is distributed under the license: MIT License
 * Originally written by Zhenrong WANG
 * mailto: zhenrongwang@live.com
 */

#ifndef REGION3_H_
#define REGION3_H_

#include "if97_general.h"

double t3abp(double pres);
double t3cdp(double pres);
double t3efp(double pres);
double t3ghp(double pres);
double t3ijp(double pres);
double t3jkp(double pres);
double t3mnp(double pres);
double t3opp(double pres);
double t3qup(double pres);
double t3rxp(double pres);
double tsatp(double pres);
double psatt(double temp);
double t3uvp(double pres);
double t3wxp(double pres);
int pt_subregion_calc1(double pres, double temp);
int pt_subregion_calc2(double pres, double temp);
double calc_vpt(int sr, double pres, double temp);
double spe_vol_pt(double pres, double temp);
double phi(double delta, double tau);
double phid(double delta, double tau);
double phidd(double delta, double tau);
double phit(double delta, double tau);
double phitt(double delta, double tau);
double phidt(double delta, double tau);
void r3prop_calc_tr(double temp, double dens, steam_prop* p_prop);
void calcFD_dFD_pt(double* fd, double* dfd, double pres, double temp, double del);
int hp_dens_pt(double* dens_final, double dens_ini, double pres, double temp);
void calcFT_dFT_prr3(double* ft, double* dft, double pres, double dens, double tau);
int temp_pr_r3(double* temp, double pres, double dens);
void calcFG_DFG_ph(double* f, double* g, double* dfdd, double* dgdd, double* dfdt, double* dgdt, double pres, double spe_enth, double del, double tau);
int tr_ph_r3(double* temp, double* dens, double pres, double spe_enth);
void calcFG_DFG_ps(double* f, double* g, double* dfdd, double* dgdd, double* dfdt, double* dgdt, double pres, double spe_entr, double del, double tau);
int tr_ps_r3(double* temp, double* dens, double pres, double spe_entr);
void calcFG_DFG_pu(double* f, double* g, double* dfdd, double* dgdd, double* dfdt, double* dgdt, double pres, double spe_ener, double del, double tau);
int tr_pu_r3(double* temp, double* dens, double pres, double spe_ener);
void calcFG_DFG_hs(double* f, double* g, double* dfdd, double* dgdd, double* dfdt, double* dgdt, double spe_enth, double spe_entr, double del, double tau);
int pres_tr_r3(double* pres, double temp, double dens);
int tr_hs_r3(double* temp, double* dens, double spe_enth, double spe_entr);
void calcFd_DFd_th(double* fd, double* dfd, double temp, double spe_enth, double del);
int dens_th_r3(double* dens, double temp, double spe_enth);
int pres_th_r3(double* pres, double temp, double spe_enth);
int pt_hs_r3(double* pres, double* temp, double spe_enth, double spe_entr);
void calcFd_DFd_tu(double* fd, double* dfd, double temp, double spe_ener, double del);
int dens_tu_r3(double* dens, double temp, double spe_ener);
int pres_tu_r3(double* pres, double temp, double spe_ener);
void calcFd_DFd_ts(double* fd, double* dfd, double temp, double spe_entr, double del);
int dens_ts_r3(double* dens, double temp, double spe_entr);
int pres_ts_r3(double* pres, double temp, double spe_entr);
int steam_prop_calc_r3(double pressure, double temp, steam_prop* p_prop);

#endif