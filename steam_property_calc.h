//units: pressure-Pa, temp-K, spe_vol-m3/kg, dens-kg/m3, spe_ener-J/kg, spe_entr-J/kgK, spe_enth-J/kg, spe_h-J/kgK


#ifndef STEAM_PROPERTY_CALC_H_
#define STEAM_PROPERTY_CALC_H_

#define MAX_PRES 1e8
#define MAX_TEMP 2273.15
#define INTER_TEMP1 623.15
#define INTER_TEMP2 1073.15
#define MIN_TEMP 273.15
#define INTER_PRES 5e7

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

int NRP(double pres)
{
	if(isnan(pres)!=0)
	{
		return -1;
	}
	else
	{
		if(pres<1e-10||pres>1e8)
		{
			return -1;
		}
	}
	return 0;
} 

int NRT(double temp)
{
	if(isnan(temp)!=0)
	{
		return -1;
	}
	else
	{
		if(temp<273.15||temp>2273.15)
		{
			return -1;
		}
	}
	return 0;
}

int NRR(double dens)
{
	if(isnan(dens)!=0)
	{
		return -1;
	}
	else
	{
		if(dens<0||dens>1e10)
		{
			return -1;
		}
	}
	return 0;
}

#endif