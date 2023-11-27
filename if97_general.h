/*
 * This code is distributed under the license: MIT License
 * Originally written by Zhenrong WANG
 * mailto: zhenrongwang@live.com
 */

//units: pressure-Pa, temp-K, spe_vol-m3/kg, dens-kg/m3, spe_ener-J/kg, spe_entr-J/kgK, spe_enth-J/kg, spe_h-J/kgK

#ifndef IF97_GENERAL_H_
#define IF97_GENERAL_H_

#define MAX_PRES 1e8
#define MAX_TEMP 2273.15
#define INTER_TEMP1 623.15
#define INTER_TEMP2 1073.15
#define MIN_TEMP 273.15
#define INTER_PRES 5e7
#define MAX_ITER_NUMSR2 200
#define MAX_ITER_NUMS 100
#define MAX_ITER_TIMES 900
#define ERR_FOR_DENS 1e-12
#define ZERO 1e-9
#define MAX_ITER_TIMES_R5 200

#define T_C 647.096
#define P_C 22.064e6
#define r_c 322

#define GAS_CONST_STEAM 461.526

typedef struct{
	double pressure;
	double temp;
	double spe_vol;
	double dens;
	double spe_energy;
	double spe_entr;
	double spe_enth;
	double spe_h_v;
	double spe_h_p;
	double speed_sound;
	double vapor_fraction;
	double drdp;
}steam_prop;

int NRP(double pres);
int NRT(double temp);
int NRR(double dens);

#endif