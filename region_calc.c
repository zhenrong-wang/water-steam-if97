/*
 * This code is distributed under the license: MIT License
 * Originally written by Zhenrong WANG
 * mailto: zhenrongwang@live.com
 */

#include <stdio.h>
#include <math.h>
#include "if97_general.h"
#include "region_calc.h"

static double sat_coeff[10]={
	0.11670521452767e4,-0.72421316703206e6,-0.17073846940092e2,0.12020824702470e5,-0.32325550322333e7,0.14915108613530e2,-0.48232657361591e4,0.40511340542057e6,-0.23855557567849,0.65017534844798e3
};

static double bound2_3[5]={
	0.34805185628969e3,-0.11671859879975e1,0.10192970039326e-2,0.57254459862746e3,0.13918839778870e2
};

double calc_sat_pres(double temp)
{
	double phi;
	double A,B,C;
	double D;
	double p_ratio;
	phi=temp+sat_coeff[8]/(temp-sat_coeff[9]);
	A=phi*phi+sat_coeff[0]*phi+sat_coeff[1];
	B=phi*phi*sat_coeff[2]+phi*sat_coeff[3]+sat_coeff[4];
	C=sat_coeff[5]*phi*phi+sat_coeff[6]*phi+sat_coeff[7];
	D=2*C/((-B)+sqrt(B*B-4*A*C));
	p_ratio=pow(D,4);
	return p_ratio*1e6;
}

double calc_bound23(double temp)
{
	double p_ratio;
	p_ratio=bound2_3[0]+bound2_3[1]*temp+bound2_3[2]*temp*temp;
	return p_ratio*1e6;
}


int region(double pressure, double temp)
{
	if(pressure<0||pressure>MAX_PRES)
	{
		return -1;
	}
	else if(temp<MIN_TEMP||temp>MAX_TEMP)
	{
		return -1;
	}
	else if(temp>INTER_TEMP2&&pressure>INTER_PRES)
	{
		return -1;
	}
	else if(temp>INTER_TEMP2&&pressure<INTER_PRES)
	{
		return 5;
	}
	else if(temp<INTER_TEMP1+1e-8)
	{
		if(fabs(pressure/calc_sat_pres(temp)-1)<1e-6)
		{
			return 4;
		}
		else if(pressure>calc_sat_pres(temp))
		{
			return 1;
		}
		else if(pressure<calc_sat_pres(temp))
		{
			return 2;
		}
	}
	else if(temp>INTER_TEMP1&&temp<INTER_TEMP2)
	{
		if(pressure>calc_bound23(temp))
		{
			return 3;
		}
		else 
		{
			return 2;
		}
	}
}


/*int main()
{
	printf("\n%lf\n",calc_sat_pres(341))
	printf("%d",region(99568,341));
}
/*	double temp2=373.15;
	double temp3=623.15;
//	for(i=0;i<10;i++)
//	{
//		printf("%lf\n",sat_coeff[i]);
//	}
	
	printf("%lf\n%lf\n",calc_sat_pres(temp2),calc_bound23(temp3));
//	for(i=0;i<3;i++)
//	{
//		printf("%lf\n",calc_sat_pres(temp[i]));
//	}
}*/