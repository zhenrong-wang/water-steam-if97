/*
 * This code is distributed under the license: MIT License
 * Originally written by Zhenrong WANG
 * mailto: zhenrongwang@live.com
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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