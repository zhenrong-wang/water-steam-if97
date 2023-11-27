/*
 * This code is distributed under the license: MIT License
 * Originally written by Zhenrong WANG
 * mailto: zhenrongwang@live.com
 */

//units: pressure-Pa, temp-K, spe_vol-m3/kg, dens-kg/m3, spe_ener-J/kg, spe_entr-J/kgK, spe_enth-J/kg, spe_h-J/kgK

#ifndef STEAM_PROPERTY_CALC_H_
#define STEAM_PROPERTY_CALC_H_

void print_prop(steam_prop *p_prop);
int hs_find_tsat(double *t, double spe_enth, double spe_entr);
//normal version 
int steam_prop_calc_pt(double pressure, double temp, steam_prop *p_prop, steam_prop *p_prop_bkup);
int steam_prop_calc_pr(double pressure, double dens, steam_prop* p_prop);
int steam_prop_calc_pu(double pressure, double spe_ener, steam_prop* p_prop);
int steam_prop_calc_ph(double pressure, double spe_enth, steam_prop* p_prop);
int steam_prop_calc_ps(double pressure, double spe_entr, steam_prop* p_prop);
int steam_prop_calc_tr(double temp, double dens, steam_prop* p_prop);
int steam_prop_calc_tu(double temp, double spe_ener, steam_prop* p_prop);
int steam_prop_calc_th(double temp, double spe_enth, steam_prop* p_prop);
int steam_prop_calc_ts(double temp, double spe_entr, steam_prop* p_prop);
int steam_prop_calc_hs(double spe_enth, double spe_entr, steam_prop* p_prop);
int steam_prop_calc_px(double pres, double vf, steam_prop *p_prop);
int steam_prop_calc_tx(double temp, double vf, steam_prop *p_prop);
double steam_visc_calc(double temp, double dens);
int calc_index_lmd2(double dens, double ref_dens);
double steam_thcond_calc(double temp, double dens, double visc, double cp, double cv, double drdp, double pressure);
//special for THE PROGRAM
int steam_prop_calc(double pressure, double temp, double *dens,double *cp, double* cv, double* spe_enth, double* drdp);

#endif