#include<stdio.h>
#include<math.h>
#include<time.h> 
#include"steam_property_calc.h"
#include"steam_property_calc.c"

int welcome(void)
{
    time_t rtime;
    struct tm* timeinfo=NULL;
    time(&rtime);
    timeinfo=localtime(&rtime);
	printf("/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */\n");
	printf("/*               Water and Steam Properties Calculation.              */\n");
	printf("/*                    VERSION 1.2(beta) LICENSE: MIT                  */\n");
	printf("/*                                                                    */\n");
	printf("/* WANG ZHENRONG (Edison WANG), reserves all rights of this program.  */\n");
	printf("/* Contacts: zhenrong_w@163.com(email).  K495458966(wechat).          */\n");
	printf("/* Main function: calculating properties by given two values:         */\n");
	printf("/*       1-pt, 2-pr, 3-pu, 4-ph, 5-ps, 6-tr, 7-tu, 8-th, 9-ts, 10-hs. */\n"); 
	printf("/*       11-px, 12-tx (x - vapor fraction).                           */\n"); 
	printf("/*       Input File: _input.dat; Output File: _properties.dat.        */\n");
	printf("/* Range1: 273.15K<T<1073.15K     && P<=100MPa;                       */\n");
	printf("/* Range2: 1073.15.15K<T<2273.15K && P<=50MPa;                        */\n");
	printf("/* NOTE: For detailed information, please read the help doc.          */\n");
	printf("/*       Any bugs found, please contact the author.                   */\n");
	printf("/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */\n");
	printf("\n# CURRENT DATE AND TIME: %s\n",asctime(timeinfo));
	return 0;
}

int main()
{
	steam_prop prop,prop_bkup;
	double viscous;
	double thcond;
	FILE* fin;
	FILE* fout;
	int type,flag,flag2,i=0;
	double v1,v2;
	char sep1,sep2,enter;

	welcome();
	
	fin=fopen("_input.dat","r");
	if(fin==NULL)
	{
		printf("\n! FATAL ERROR: input file not found.\n! Program abort.\n! Please press any key to exit.");
		printf("\n@ Any problems found, please contact the author.\n@ Zhenrong Wang, zhenrong_w@163.com, K495458966(wechat).\n@ All rights reserved.\n");		
		fflush(stdin);
		getchar();
		return -1;
	}
	fout=fopen("_properties.dat","w");
	if(fout==NULL)
	{
		printf("\n! FATAL ERROR: cannot creat output file.\n! Program abort.\n! Please press any key to exit.");
		printf("\n@ Any problems found, please contact the author.\n@ Zhenrong Wang, zhenrong_w@163.com, K495458966(wechat).\n@ All rights reserved.\n");
		fclose(fin);
		fflush(stdin);
		getchar();
		return -1;
	}
	
	printf("\n# Please press ENTER to start calculation:");
	fflush(stdin);
	getchar(); 
	
	fprintf(fout,"\tPRES\t\tTEMP\tSPE_VOL\t\tDENS\t\tu\t\th\t\ts\t\tCv\t\tCp\t\tVsound\t\tVisc\t\tLamd\tV_FRAC\n");
	fprintf(fout,"\tMPa\t\tK\tm^3/kg\t\tkg/m^3\t\tkJ/kg\t\tkJ/kg\t\tkJ/(kg*K)\tkJ/(kg*K)\tkJ/(kg*K)\tm/s\t\tPa.s\t\tW/(m*K)\t-\n");
	while(!feof(fin))
	{
		i++;
		fscanf(fin,"%d%c%lf%c%lf%c",&type,&sep1,&v1,&sep2,&v2,&enter);
		if(type>0&&type<15)
		{
			printf("\n%d\t%.8lf\t%.8lf",type,v1,v2);
		}
	
		if(type==1)
			flag=steam_prop_calc_pt(v1,v2,&prop,&prop_bkup);
		else if(type==2)
			flag=steam_prop_calc_pr(v1,v2,&prop);
		else if(type==3)
			flag=steam_prop_calc_pu(v1,v2,&prop);
		else if(type==4)
			flag=steam_prop_calc_ph(v1,v2,&prop);
		else if(type==5)
			flag=steam_prop_calc_ps(v1,v2,&prop);
		else if(type==6)
			flag=steam_prop_calc_tr(v1,v2,&prop);
		else if(type==7)
			flag=steam_prop_calc_tu(v1,v2,&prop);
		else if(type==8)
			flag=steam_prop_calc_th(v1,v2,&prop);
		else if(type==9)
			flag=steam_prop_calc_ts(v1,v2,&prop);
		else if(type==10)
			flag=steam_prop_calc_hs(v1,v2,&prop);
		else if(type==11)
			flag=steam_prop_calc_px(v1,v2,&prop);
		else if(type==12)
			flag=steam_prop_calc_tx(v1,v2,&prop);
		else if(type==13)
			flag=steam_prop_calc_pt(v1,tsatp(v1),&prop,&prop_bkup);
		else if(type==14)
			flag=steam_prop_calc_pt(calc_sat_pres(v1),v1,&prop,&prop_bkup);
		else
		{
			continue;
		}
					
		if(flag==-1)
		{
			printf("\n! WARNING: Calculation error at line # %d of the input file.\n", i);
			fprintf(fout,"%d\t-1\n",i);
			continue;
		}
		else if(flag==1)
		{
			viscous=steam_visc_calc(prop.temp,prop.dens);
			thcond=steam_thcond_calc(prop.temp,prop.dens,viscous,prop.spe_h_p,prop.spe_h_v,prop.drdp,prop.pressure);
			fprintf(fout,"%d\t%.6lf*\t%.2lf*\t%.8lf\t%.6lf\t%8.4lf\t%8.4lf\t%.8lf\t%.6lf\t%.6lf\t%.4lf\t%.8lf\t%.4lf\t%.4lf\n",i,prop.pressure/1e6,prop.temp,prop.spe_vol,prop.dens,prop.spe_energy/1000,prop.spe_enth/1000,prop.spe_entr/1000,prop.spe_h_v/1000,prop.spe_h_p/1000,prop.speed_sound,viscous,thcond,prop.vapor_fraction);
			fprintf(fout,"\t%.6lf**\t%.2lf*\t%.8lf\t%.6lf\t%8.4lf\t%8.4lf\t%.8lf\t%.6lf\t%.6lf\t%.4lf\t%.8lf\t%.4lf\t%.4lf\n",prop_bkup.pressure/1e6,prop_bkup.temp,prop_bkup.spe_vol,prop_bkup.dens,prop_bkup.spe_energy/1000,prop_bkup.spe_enth/1000,prop_bkup.spe_entr/1000,prop_bkup.spe_h_v/1000,prop_bkup.spe_h_p/1000,prop_bkup.speed_sound,viscous,thcond,prop_bkup.vapor_fraction);		
		}
		else
		{
			viscous=steam_visc_calc(prop.temp,prop.dens);
			thcond=steam_thcond_calc(prop.temp,prop.dens,viscous,prop.spe_h_p,prop.spe_h_v,prop.drdp,prop.pressure);
			fprintf(fout,"%d\t%.6lf\t%.2lf\t%.8lf\t%.6lf\t%.4lf\t%.4lf\t%.8lf\t%.6lf\t%.6lf\t%.4lf\t%.8lf\t%.4lf\t%.4lf\n",i,prop.pressure/1e6,prop.temp,prop.spe_vol,prop.dens,prop.spe_energy/1000,prop.spe_enth/1000,prop.spe_entr/1000,prop.spe_h_v/1000,prop.spe_h_p/1000,prop.speed_sound,viscous,thcond,prop.vapor_fraction);
		}	
	}
	
	printf("\n\n# Calculation finished. Please check _properties.dat for results.\n# Press any key to exit.\n\n@ Any problems found, please contact the author.\n@ Zhenrong Wang, zhenrong_w@163.com, K495458966(wechat).\n@ All rights reserved.\n");
	fclose(fin);
	fclose(fout);
	fflush(stdin);
	getchar();
	return 0;
	
}
