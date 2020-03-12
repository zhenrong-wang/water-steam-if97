#include<stdio.h>
#include<math.h>
#include"steam_property_calc.h"
#include"region_calc.c"
#include"region1.c"
#include"region2.c"
#include"region3.c"
#include"region5.c"

static double coeff_visc1[4]={
	1.67752,2.20462,0.6366564,-0.241605
};

static double coeff_visc2[6][7]={
	5.20094e-1,2.22531e-1,-2.81378e-1,1.61913e-1,-3.25372e-2,0,0,
	8.50895e-2,9.99115e-1,-9.06851e-1,2.57399e-1,0,0,0,
	-1.08374,1.88797,-7.72479e-1,0,0,0,0,
	-2.89555e-1,1.26613,-4.89837e-1,0,6.98452e-2,0,-4.35673e-3,
	0,0,-2.57040e-1,0,0,8.72102e-3,0,
	0,1.20573e-1,0,0,0,0,-5.93264e-4,
}; 

static double coeff_thcond1[5]={
	2.443221e-3,1.323095e-2,6.770357e-3,-3.454586e-3,4.096266e-4
};

static double coeff_thcond2[5][6]={
	1.60397357,-0.646013523,0.111443906,0.102997357,-0.0504123634,0.00609859258,
	2.33771842,-2.78843778,1.53616167,-0.463045512,0.0832827019,-0.00719201245,
	2.19650529,-4.54580785,3.55777244,-1.40944978,0.275418278,-0.0205938816,
	-1.21051378,1.60812989,-0.621178141,0.0716373224,0,0,
	-2.7203370,4.57586331,-3.18369245,1.1168348,-0.19268305,0.012913842,
};

static double coeff_thcond3[6][5]={
	6.53786807199516,6.52717759281799,5.35500529896124,1.55225959906681,1.11999926419994,
	-5.61149954923348,-6.30816983387575,-3.96415689925446,0.464621290821181,0.595748562571649,
	3.39624167361325,8.08379285492595,8.91990208918795,8.93237374861479,9.8895256507892,
	-2.27492629730878,-9.82240510197603,-12.033872950579,-11.0321960061126,-10.325505114704,
	10.2631854662709,12.1358413791395,9.19494865194302,6.1678099993336,4.66861294457414,
	1.97815050331519,-5.54349664571295,-2.16866274479712,-0.965458722086812,-0.503243546373828,
};

static double coeff_thcond4[6][5]={
	6.53786807199516,6.52717759281799,5.35500529896124,1.55225959906681,1.11999926419994,
	-5.61149954923348,-6.30816983387575,-3.96415689925446,0.464621290821181,0.595748562571649,
	3.39624167361325,8.08379285492595,8.91990208918795,8.93237374861479,9.8895256507892,
	-2.27492629730878,-9.82240510197603,-12.033872950579,-11.0321960061126,-10.325505114704,
	10.2631854662709,12.1358413791395,9.19494865194302,6.1678099993336,4.66861294457414,
	1.97815050331519,-5.54349664571295,-2.16866274479712,-0.965458722086812,-0.503243546373828,
};

int hs_find_tsat(double *t, double spe_enth, double spe_entr)
{
	double tt1,tt2,ps1,ps2,xs1,xs2,xh1,xh2,ttm,psm,xhm,xsm;
	steam_prop propw,props,propmw,propms;
	int i=0;
	
	tt1=MIN_TEMP;
	tt2=INTER_TEMP1;
	ps1=calc_sat_pres(tt1);
	ps2=calc_sat_pres(tt2);
	steam_prop_calc_r1(ps1,tt1,&propw);
	steam_prop_calc_r2(ps1,tt1,&props);
	xs1=(spe_entr-propw.spe_entr)/(props.spe_entr-propw.spe_entr);
	xh1=(spe_enth-propw.spe_enth)/(props.spe_enth-propw.spe_enth);
//	printf("\n%lf,%lf\n",xh1,xs1);
	steam_prop_calc_r1(ps2,tt2,&propw);
	steam_prop_calc_r2(ps2,tt2,&props);
	xs2=(spe_entr-propw.spe_entr)/(props.spe_entr-propw.spe_entr);
	xh2=(spe_enth-propw.spe_enth)/(props.spe_enth-propw.spe_enth);
//	printf("\n%lf,%lf\n",xh2,xs2);
	if((xh1-xs1)*(xh2-xs2)>0)
	{
		*t=-1;
		return -1;
	}
	do
	{
		ttm=0.5*(tt1+tt2);
//		printf("\n............%lf\n",ttm);
		psm=calc_sat_pres(ttm);
//		printf("\n...........%lf,%lf\n",psm,ttm);
		steam_prop_calc_r1(psm,ttm,&propmw);
		steam_prop_calc_r2(psm,ttm,&propms);
		xsm=(spe_entr-propmw.spe_entr)/(propms.spe_entr-propmw.spe_entr);
		xhm=(spe_enth-propmw.spe_enth)/(propms.spe_enth-propmw.spe_enth);
//			printf("\n...........%lf,%lf\n",xhm,xsm);
		i++;
		if(fabs(xhm-xsm)<1e-6)
		{
			*t=ttm;
			return 0;
		}
		else if((xhm-xsm)*(xh1-xs1)>0)
		{
			tt1=ttm;
		}
		else if((xhm-xsm)*(xh2-xs2)>0)
		{
			tt2=ttm;
		}
	}while(i<200);
	if(i==200)
	{
		*t=-1;
		return -1;
	}
}

//normal version 
int steam_prop_calc_pt(double pressure, double temp, steam_prop *p_prop, steam_prop *p_prop_bkup)
{
	int region_code;
	int calc_flag;
	
	p_prop_bkup->pressure=-1;
	p_prop_bkup->temp=-1;
	p_prop_bkup->spe_vol=-1;
	p_prop_bkup->dens=-1;
	p_prop_bkup->spe_energy=-1;
	p_prop_bkup->spe_entr=-1;
	p_prop_bkup->spe_enth=-1;
	p_prop_bkup->spe_h_v=-1;
	p_prop_bkup->spe_h_p=-1;
	p_prop_bkup->speed_sound=-1;
	p_prop_bkup->drdp=-1;
	
	region_code=region(pressure,temp);
//	printf("!!!!%d\n",region_code); 
	if(region_code==-1)
	{
		printf("\nFATAL ERROR: PRESSURE/TEMP EXCEEDS BOUNDARIES.\n");
		return -1;
	}
	else if(region_code==1)
	{
		steam_prop_calc_r1(pressure,temp,p_prop);
		return 0;
	}
	else if(region_code==2)
	{
		steam_prop_calc_r2(pressure,temp,p_prop);
		return 0;
	}
	
	else if(region_code==4)
	{
		steam_prop_calc_r1(pressure,temp,p_prop);
		steam_prop_calc_r2(pressure,temp,p_prop_bkup);
		return 1;
	} 
	else if(region_code==3)
	{
		calc_flag=steam_prop_calc_r3(pressure,temp,p_prop);
		if(calc_flag==-1)
		{
			printf("\nFATAL ERROR: DENSITY ITERATION NOT CONVERGED.\n");
			return -2;
		}
		else
		{
//			print_prop(p_prop); 
			return 0;
		}
	}
	else if(region_code==5)
	{
		steam_prop_calc_r5(pressure,temp,p_prop);
		return 0;
	}
}

int steam_prop_calc_pr(double pressure, double dens, steam_prop* p_prop)
{
	int rc,flag;
	double temp,temp_sat;
	double pres_h,pres_l,vapor_frac;
	steam_prop prop_w,prop_s;
	pres_l=calc_sat_pres(MIN_TEMP);
	pres_h=calc_sat_pres(INTER_TEMP1);
	
//	printf("\n\n%lf,%lf,%lf\n",pres_h,pres_l);
	if(pressure>pres_l&&pressure<pres_h)
	{
		temp_sat=tsatp(pressure);
		steam_prop_calc_r1(pressure,temp_sat,&prop_w);
		steam_prop_calc_r2(pressure,temp_sat,&prop_s);
//		printf("\n....%lf,%lf,%lf\n",dens,prop_w.dens,prop_s.dens);

		if(dens>prop_s.dens&&dens<prop_w.dens)
		{
//			printf("\n,,,,,,,,,,,,,,,,,,,,,,\n");
			vapor_frac=prop_s.dens*(prop_w.dens-dens)/(dens*(prop_w.dens-prop_s.dens));
			p_prop->pressure=pressure;
			p_prop->temp=temp_sat;
			p_prop->spe_vol=1/dens;
			p_prop->dens=dens;
			p_prop->spe_energy=vapor_frac*prop_s.spe_energy+(1-vapor_frac)*prop_w.spe_energy;
			p_prop->spe_entr=vapor_frac*prop_s.spe_entr+(1-vapor_frac)*prop_w.spe_entr;
			p_prop->spe_enth=vapor_frac*prop_s.spe_enth+(1-vapor_frac)*prop_w.spe_enth;
			p_prop->spe_h_v=vapor_frac*prop_s.spe_h_v+(1-vapor_frac)*prop_w.spe_h_v;
			p_prop->spe_h_p=vapor_frac*prop_s.spe_h_p+(1-vapor_frac)*prop_w.spe_h_p;
			p_prop->speed_sound=vapor_frac*prop_s.speed_sound+(1-vapor_frac)*prop_w.speed_sound;
			p_prop->drdp=vapor_frac*prop_s.drdp+(1-vapor_frac)*prop_w.drdp;
			p_prop->vapor_fraction=vapor_frac;
			//solve vapor fraction and properties;
			
			return 0;
		}
	}
	
//	printf("\n\tPREEEEEEEEEEEEEEEEEEEEEEE\n");										
	flag=temp_pr_r1(&temp,pressure,dens);
	if(flag==0)
	{
		rc=region(pressure,temp);
		steam_prop_calc_r1(pressure,temp,p_prop);
		if(rc==1&&fabs((p_prop->dens-dens)/dens<1e-4))
		{
			p_prop->dens=dens;
			return 0;
		}
	}

	flag=temp_pr_r2(&temp,pressure,dens);
	if(flag==0)
	{
		rc=region(pressure,temp);
		steam_prop_calc_r2(pressure,temp,p_prop);
		if(rc==2&&fabs((p_prop->dens-dens)/dens<1e-4))
		{
			p_prop->dens=dens;
			return 0;
		}
	}
	
	flag=temp_pr_r3(&temp,pressure,dens);
	if(flag==0)
	{
		rc=region(pressure,temp);
		steam_prop_calc_r3(pressure,temp,p_prop);
		if(rc==3&&fabs((p_prop->dens-dens)/dens<1e-4))
		{
			p_prop->dens=dens;
			return 0;
		}
	}
	flag=temp_pr_r5(&temp,pressure,dens);
	if(flag==0)
	{
		rc=region(pressure,temp);
		steam_prop_calc_r5(pressure,temp,p_prop);
		if(rc==5&&fabs((p_prop->dens-dens)/dens<1e-4))
		{
			p_prop->dens=dens;
			return 0;
		}
	}
	return -1;
}

int steam_prop_calc_pu(double pressure, double spe_ener, steam_prop* p_prop)
{
	int rc,flag;
	double temp,r3dens;
	double pres_h,pres_l,temp_sat,vapor_frac;
	steam_prop prop_s,prop_w;
	
	pres_l=calc_sat_pres(MIN_TEMP);
	pres_h=calc_sat_pres(INTER_TEMP1);
	
	if(pressure>pres_l&&pressure<pres_h)
	{
		temp_sat=tsatp(pressure);
		steam_prop_calc_r1(pressure,temp_sat,&prop_w);
		steam_prop_calc_r2(pressure,temp_sat,&prop_s);
		if(spe_ener<prop_s.spe_energy&&spe_ener>prop_w.spe_energy)
		{
			vapor_frac=(spe_ener-prop_w.spe_energy)/(prop_s.spe_energy-prop_w.spe_energy);
			p_prop->pressure=pressure;
			p_prop->temp=temp_sat;
			p_prop->spe_vol=vapor_frac*prop_s.spe_vol+(1-vapor_frac)*prop_w.spe_vol;
			p_prop->dens=1/p_prop->spe_vol;
			p_prop->spe_energy=spe_ener;
			p_prop->spe_entr=vapor_frac*prop_s.spe_entr+(1-vapor_frac)*prop_w.spe_entr;
			p_prop->spe_enth=vapor_frac*prop_s.spe_enth+(1-vapor_frac)*prop_w.spe_enth;
			p_prop->spe_h_v=vapor_frac*prop_s.spe_h_v+(1-vapor_frac)*prop_w.spe_h_v;
			p_prop->spe_h_p=vapor_frac*prop_s.spe_h_p+(1-vapor_frac)*prop_w.spe_h_p;
			p_prop->speed_sound=vapor_frac*prop_s.speed_sound+(1-vapor_frac)*prop_w.speed_sound;
			p_prop->drdp=vapor_frac*prop_s.drdp+(1-vapor_frac)*prop_w.drdp;
			p_prop->vapor_fraction=vapor_frac;
			//solve vapor fraction and properties;
			return 0;
		}
	}
	
	flag=temp_pu_r1(&temp,pressure,spe_ener);
	if(flag==0)
	{
		rc=region(pressure,temp);
		steam_prop_calc_r1(pressure,temp,p_prop);
		if(rc==1&&fabs((p_prop->spe_energy-spe_ener)/spe_ener)<1e-4)
		{
			p_prop->spe_energy=spe_ener;
			return 0;
		}
	}
	
	flag=temp_pu_r2(&temp,pressure,spe_ener);
	if(flag==0)
	{
		rc=region(pressure,temp);
		steam_prop_calc_r2(pressure,temp,p_prop);
		if(rc==2&&fabs((p_prop->spe_energy-spe_ener)/spe_ener)<1e-4)
		{
			p_prop->spe_energy=spe_ener;
			return 0;
		}
	}
	
	flag=tr_pu_r3(&temp,&r3dens,pressure,spe_ener);
	if(flag==0)
	{
		rc=region(pressure,temp);
		steam_prop_calc_r3(pressure,temp,p_prop);
		if(rc==3&&fabs((p_prop->spe_energy-spe_ener)/spe_ener)<1e-4)
		{
			p_prop->pressure=pressure;
			p_prop->spe_energy=spe_ener;
			return 0;
		}
	}
		
	flag=temp_pu_r5(&temp,pressure,spe_ener);
	if(flag==0)
	{
		rc=region(pressure,temp);
		steam_prop_calc_r5(pressure,temp,p_prop);
		if(rc==5&&fabs((p_prop->spe_energy-spe_ener)/spe_ener)<1e-4)
		{
			p_prop->spe_energy=spe_ener;
			return 0;
		}
	}
	return -1;	
}

int steam_prop_calc_ph(double pressure, double spe_enth, steam_prop* p_prop)
{
	int rc,flag;
	double temp,r3dens;
	double pres_h,pres_l,temp_sat,vapor_frac;
	steam_prop prop_s,prop_w;
	
	pres_l=calc_sat_pres(MIN_TEMP);
	pres_h=calc_sat_pres(INTER_TEMP1);
	
	if(pressure>pres_l&&pressure<pres_h)
	{
		temp_sat=tsatp(pressure);
		steam_prop_calc_r1(pressure,temp_sat,&prop_w);
		steam_prop_calc_r2(pressure,temp_sat,&prop_s);
		if(spe_enth<prop_s.spe_enth&&spe_enth>prop_w.spe_enth)
		{
			vapor_frac=(spe_enth-prop_w.spe_enth)/(prop_s.spe_enth-prop_w.spe_enth);
			p_prop->pressure=pressure;
			p_prop->temp=temp_sat;
			p_prop->spe_vol=vapor_frac*prop_s.spe_vol+(1-vapor_frac)*prop_w.spe_vol;
			p_prop->dens=1/p_prop->spe_vol;
			p_prop->spe_energy=vapor_frac*prop_s.spe_energy+(1-vapor_frac)*prop_w.spe_energy;
			p_prop->spe_entr=vapor_frac*prop_s.spe_entr+(1-vapor_frac)*prop_w.spe_entr;
			p_prop->spe_enth=spe_enth;
			p_prop->spe_h_v=vapor_frac*prop_s.spe_h_v+(1-vapor_frac)*prop_w.spe_h_v;
			p_prop->spe_h_p=vapor_frac*prop_s.spe_h_p+(1-vapor_frac)*prop_w.spe_h_p;
			p_prop->speed_sound=vapor_frac*prop_s.speed_sound+(1-vapor_frac)*prop_w.speed_sound;
			p_prop->drdp=vapor_frac*prop_s.drdp+(1-vapor_frac)*prop_w.drdp;
			p_prop->vapor_fraction=vapor_frac;
			//solve vapor fraction and properties;
			return 0;
		}
	}
	
//		printf("\nEEEEEEEEEEEEEEEEE\n\n");
	flag=temp_ph_r1(&temp,pressure,spe_enth);
	if(flag==0)
	{
		rc=region(pressure,temp);
		steam_prop_calc_r1(pressure,temp,p_prop);
		if(rc==1&&fabs((p_prop->spe_enth-spe_enth)/spe_enth)<1e-4)
		{
			p_prop->spe_enth=spe_enth;
			return 0;
		}
	}
	
	flag=temp_ph_r2(&temp,pressure,spe_enth);
	if(flag==0)
	{
		rc=region(pressure,temp);
		steam_prop_calc_r2(pressure,temp,p_prop);
		if(rc==2&&fabs((p_prop->spe_enth-spe_enth)/spe_enth)<1e-4)
		{
			p_prop->spe_enth=spe_enth;
			return 0;
		}
	}
	
	flag=tr_ph_r3(&temp,&r3dens,pressure,spe_enth);
	if(flag==0)
	{
//			printf("\nEEEEEEEEEE%lf,%lf",temp,r3dens);
		rc=region(pressure,temp);
		steam_prop_calc_r3(pressure,temp,p_prop);
		if(rc==3&&fabs((p_prop->spe_enth-spe_enth)/spe_enth)<1e-4)
		{
			p_prop->pressure=pressure;
			p_prop->spe_enth=spe_enth;
			return 0;
		}
	}
	flag=temp_ph_r5(&temp,pressure,spe_enth);
	if(flag==0)
	{
		rc=region(pressure,temp);
		steam_prop_calc_r5(pressure,temp,p_prop);
		if(rc==5&&fabs((p_prop->spe_enth-spe_enth)/spe_enth)<1e-4)
		{
			p_prop->spe_enth=spe_enth;
			return 0;
		}
	}
	return -1;	
}

int steam_prop_calc_ps(double pressure, double spe_entr, steam_prop* p_prop)
{
	int rc,flag;
	double temp,r3dens;
	double pres_h,pres_l,temp_sat,vapor_frac;
	steam_prop prop_s,prop_w;
	
	pres_l=calc_sat_pres(MIN_TEMP);
	pres_h=calc_sat_pres(INTER_TEMP1);
	
	if(pressure>pres_l&&pressure<pres_h)
	{
		temp_sat=tsatp(pressure);
		steam_prop_calc_r1(pressure,temp_sat,&prop_w);
		steam_prop_calc_r2(pressure,temp_sat,&prop_s);
		if(spe_entr<prop_s.spe_entr&&spe_entr>prop_w.spe_entr)
		{
			vapor_frac=(spe_entr-prop_w.spe_entr)/(prop_s.spe_entr-prop_w.spe_entr);
			p_prop->pressure=pressure;
			p_prop->temp=temp_sat;
			p_prop->spe_vol=vapor_frac*prop_s.spe_vol+(1-vapor_frac)*prop_w.spe_vol;
			p_prop->dens=1/p_prop->spe_vol;
			p_prop->spe_energy=vapor_frac*prop_s.spe_energy+(1-vapor_frac)*prop_w.spe_energy;
			p_prop->spe_entr=spe_entr;
			p_prop->spe_enth=vapor_frac*prop_s.spe_enth+(1-vapor_frac)*prop_w.spe_enth;
			p_prop->spe_h_v=vapor_frac*prop_s.spe_h_v+(1-vapor_frac)*prop_w.spe_h_v;
			p_prop->spe_h_p=vapor_frac*prop_s.spe_h_p+(1-vapor_frac)*prop_w.spe_h_p;
			p_prop->speed_sound=vapor_frac*prop_s.speed_sound+(1-vapor_frac)*prop_w.speed_sound;
			p_prop->drdp=vapor_frac*prop_s.drdp+(1-vapor_frac)*prop_w.drdp;
			p_prop->vapor_fraction=vapor_frac;
			//solve vapor fraction and properties;
			return 0;
		}
	}
	
	flag=temp_ps_r1(&temp,pressure,spe_entr);
	if(flag==0)
	{
		rc=region(pressure,temp);
		steam_prop_calc_r1(pressure,temp,p_prop);
		if(rc==1&&fabs((p_prop->spe_entr-spe_entr)/spe_entr)<1e-4)
		{
			p_prop->spe_entr=spe_entr;
			return 0;
		}
	}
	
	flag=temp_ps_r2(&temp,pressure,spe_entr);
	if(flag==0)
	{
		rc=region(pressure,temp);
		steam_prop_calc_r2(pressure,temp,p_prop);
		if(rc==2&&fabs((p_prop->spe_entr-spe_entr)/spe_entr)<1e-4)
		{
			p_prop->spe_entr=spe_entr;
			return 0;
		}
	}
	
	flag=tr_ps_r3(&temp,&r3dens,pressure,spe_entr);
//	printf("\t#########%lf,%lf",temp,1/r3dens);
	if(flag==0)
	{
		rc=region(pressure,temp);
		steam_prop_calc_r3(pressure,temp,p_prop);
		if(rc==3&&fabs((p_prop->spe_entr-spe_entr)/spe_entr)<1e-4)
		{
			p_prop->pressure=pressure;
			p_prop->spe_entr=spe_entr;
			return 0;
		}
	}
	
	flag=temp_ps_r5(&temp,pressure,spe_entr);
	if(flag==0)
	{
		rc=region(pressure,temp);
		steam_prop_calc_r5(pressure,temp,p_prop);
		if(rc==5&&fabs((p_prop->spe_entr-spe_entr)/spe_entr)<1e-4)
		{
			p_prop->spe_entr=spe_entr;
			return 0;
		}
	}
	return -1;	
}

int steam_prop_calc_tr(double temp, double dens, steam_prop* p_prop)
{
	int rc,flag;
	double pres,pres_sat,vapor_frac;
	steam_prop prop_w,prop_s;

	if(temp>INTER_TEMP1)
	{
		flag=pres_tr_r2(&pres,temp,dens);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r2(pres,temp,p_prop);
			if(rc==2&&fabs((p_prop->dens-dens)/dens)<1e-4)
			{
				p_prop->dens=dens;
				return 0;
			}
		}
	
		flag=pres_tr_r3(&pres,temp,dens);
//	printf("\n.............%d,%lf\n",flag,pres);
		if(flag==0)
		{
			rc=region(pres,temp);
//			printf("...........%d\n",rc);
			steam_prop_calc_r3(pres,temp,p_prop);
			if(rc==3&&fabs((p_prop->dens-dens)/dens)<1e-4)
			{
				p_prop->temp=temp;
				p_prop->dens=dens;
				return 0;
			}
		}
		
		flag=pres_tr_r5(&pres,temp,dens);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r5(pres,temp,p_prop);
			if(rc==5&&fabs((p_prop->dens-dens)/dens)<1e-4)
			{
				p_prop->dens=dens;
				return 0;
			}
		}
		return -1;	
	}
	
	else
	{
		flag=pres_tr_r1(&pres,temp,dens);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r1(pres,temp,p_prop);
			if(rc==1&&fabs((p_prop->dens-dens)/dens)<1e-4)
			{
				p_prop->dens=dens;
				return 0;
			}
		}
	
		flag=pres_tr_r2(&pres,temp,dens);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r2(pres,temp,p_prop);
			if(rc==2&&fabs((p_prop->dens-dens)/dens)<1e-4)
			{
				p_prop->dens=dens;
				return 0;
			}
		}
		pres_sat=calc_sat_pres(temp);
		steam_prop_calc_r1(pres_sat,temp,&prop_w);
		steam_prop_calc_r2(pres_sat,temp,&prop_s);
		
		if(dens>prop_s.dens&&dens<prop_w.dens)
		{
			vapor_frac=prop_s.dens*(prop_w.dens-dens)/(dens*(prop_w.dens-prop_s.dens));
			p_prop->pressure=pres_sat;
			p_prop->temp=temp;
			p_prop->spe_vol=1/dens;
			p_prop->dens=dens;
			p_prop->spe_energy=vapor_frac*prop_s.spe_energy+(1-vapor_frac)*prop_w.spe_energy;
			p_prop->spe_entr=vapor_frac*prop_s.spe_entr+(1-vapor_frac)*prop_w.spe_entr;
			p_prop->spe_enth=vapor_frac*prop_s.spe_enth+(1-vapor_frac)*prop_w.spe_enth;
			p_prop->spe_h_v=vapor_frac*prop_s.spe_h_v+(1-vapor_frac)*prop_w.spe_h_v;
			p_prop->spe_h_p=vapor_frac*prop_s.spe_h_p+(1-vapor_frac)*prop_w.spe_h_p;
			p_prop->speed_sound=vapor_frac*prop_s.speed_sound+(1-vapor_frac)*prop_w.speed_sound;
			p_prop->drdp=vapor_frac*prop_s.drdp+(1-vapor_frac)*prop_w.drdp;
			p_prop->vapor_fraction=vapor_frac;
			//solve vapor fraction and properties;
			return 0;
		}
		return -1;
	}
}

int steam_prop_calc_tu(double temp, double spe_ener, steam_prop* p_prop)
{
	int rc,flag;
	double pres,pres_sat,vapor_frac;
	steam_prop prop_w,prop_s;
	
	if(temp>INTER_TEMP1)
	{
		flag=pres_tu_r2(&pres,temp,spe_ener);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r2(pres,temp,p_prop);
			if(rc==2&&fabs((p_prop->spe_energy-spe_ener)/spe_ener)<1e-4)
			{
				p_prop->spe_energy=spe_ener;
				return 0;
			}
		}
	
		flag=pres_tu_r3(&pres,temp,spe_ener);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r3(pres,temp,p_prop);
			if(rc==3&&fabs((p_prop->spe_energy-spe_ener)/spe_ener)<1e-4)
			{
				p_prop->temp=temp;
				p_prop->spe_energy=spe_ener;
				return 0;
			}
		}
		
		flag=pres_tu_r5(&pres,temp,spe_ener);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r5(pres,temp,p_prop);
			if(rc==5&&fabs((p_prop->spe_energy-spe_ener)/spe_ener)<1e-4)
			{
				p_prop->spe_energy=spe_ener;
				return 0;
			}
		}
		return -1;	
	}
	
	else
	{		
		flag=pres_tu_r1(&pres,temp,spe_ener);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r1(pres,temp,p_prop);
			if(rc==1&&fabs((p_prop->spe_energy-spe_ener)/spe_ener)<1e-4)
			{
				p_prop->spe_energy=spe_ener;
				return 0;
			}
		}
		
		flag=pres_tu_r2(&pres,temp,spe_ener);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r2(pres,temp,p_prop);
			if(rc==2&&fabs((p_prop->spe_energy-spe_ener)/spe_ener)<1e-4)
			{
				p_prop->spe_energy=spe_ener;
				return 0;
			}
		}
		
				pres_sat=calc_sat_pres(temp);
		steam_prop_calc_r1(pres_sat,temp,&prop_w);
		steam_prop_calc_r2(pres_sat,temp,&prop_s);
		
		if(spe_ener>prop_w.spe_energy&&spe_ener<prop_s.spe_energy)
		{
			vapor_frac=(spe_ener-prop_w.spe_energy)/(prop_s.spe_energy-prop_w.spe_energy);
			p_prop->pressure=pres_sat;
			p_prop->temp=temp;
			p_prop->spe_vol=vapor_frac*prop_s.spe_vol+(1-vapor_frac)*prop_w.spe_vol;
			p_prop->dens=1/p_prop->spe_vol;
			p_prop->spe_energy=spe_ener;
			p_prop->spe_entr=vapor_frac*prop_s.spe_entr+(1-vapor_frac)*prop_w.spe_entr;
			p_prop->spe_enth=vapor_frac*prop_s.spe_enth+(1-vapor_frac)*prop_w.spe_enth;
			p_prop->spe_h_v=vapor_frac*prop_s.spe_h_v+(1-vapor_frac)*prop_w.spe_h_v;
			p_prop->spe_h_p=vapor_frac*prop_s.spe_h_p+(1-vapor_frac)*prop_w.spe_h_p;
			p_prop->speed_sound=vapor_frac*prop_s.speed_sound+(1-vapor_frac)*prop_w.speed_sound;
			p_prop->drdp=vapor_frac*prop_s.drdp+(1-vapor_frac)*prop_w.drdp;
			p_prop->vapor_fraction=vapor_frac;
			//solve vapor fraction and properties;
			return 0;
		}
		return -1;
	}
}

int steam_prop_calc_th(double temp, double spe_enth, steam_prop* p_prop)
{
	int rc,flag;
	double pres,pres_sat,vapor_frac;
	steam_prop prop_w,prop_s;
	
	if(temp>INTER_TEMP1)
	{
		flag=pres_th_r2(&pres,temp,spe_enth);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r2(pres,temp,p_prop);
			if(rc==2&&fabs((p_prop->spe_enth-spe_enth)/spe_enth)<1e-4)
			{
				p_prop->spe_enth=spe_enth;
				return 0;
			}
		}
	
		flag=pres_th_r3(&pres,temp,spe_enth);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r3(pres,temp,p_prop);
			if(rc==3&&fabs((p_prop->spe_enth-spe_enth)/spe_enth)<1e-4)
			{
				p_prop->temp=temp;
				p_prop->spe_enth=spe_enth;
				return 0;
			}
		}
		
		flag=pres_th_r5(&pres,temp,spe_enth);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r5(pres,temp,p_prop);
			if(rc==5&&fabs((p_prop->spe_enth-spe_enth)/spe_enth)<1e-4)
			{
				p_prop->spe_enth=spe_enth;
				return 0;
			}
		}
		return -1;	
	}
	
	else
	{
		flag=pres_th_r1(&pres,temp,spe_enth);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r1(pres,temp,p_prop);
			if(rc==1&&fabs((p_prop->spe_enth-spe_enth)/spe_enth)<1e-4)
			{
				p_prop->spe_enth=spe_enth;
				return 0;
			}
		}	
		
		flag=pres_th_r2(&pres,temp,spe_enth);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r2(pres,temp,p_prop);
			if(rc==2&&fabs((p_prop->spe_enth-spe_enth)/spe_enth)<1e-4)
			{
				p_prop->spe_enth=spe_enth;
				return 0;
			}
		}
	
		pres_sat=calc_sat_pres(temp);
		steam_prop_calc_r1(pres_sat,temp,&prop_w);
		steam_prop_calc_r2(pres_sat,temp,&prop_s);
		printf("\n\t%lf,%lf\n",prop_w.spe_enth,prop_s.spe_enth);
		if(spe_enth>prop_w.spe_enth&&spe_enth<prop_s.spe_enth)
		{
			printf("\n\nEEEEEEEEEEE\n");
			vapor_frac=(spe_enth-prop_w.spe_enth)/(prop_s.spe_enth-prop_w.spe_enth);
			p_prop->pressure=pres_sat;
			p_prop->temp=temp;
			p_prop->spe_vol=vapor_frac*prop_s.spe_vol+(1-vapor_frac)*prop_w.spe_vol;
			p_prop->dens=1/p_prop->spe_vol;
			p_prop->spe_energy=vapor_frac*prop_s.spe_energy+(1-vapor_frac)*prop_w.spe_energy;
			p_prop->spe_entr=vapor_frac*prop_s.spe_entr+(1-vapor_frac)*prop_w.spe_entr;
			p_prop->spe_enth=spe_enth;
			p_prop->spe_h_v=vapor_frac*prop_s.spe_h_v+(1-vapor_frac)*prop_w.spe_h_v;
			p_prop->spe_h_p=vapor_frac*prop_s.spe_h_p+(1-vapor_frac)*prop_w.spe_h_p;
			p_prop->speed_sound=vapor_frac*prop_s.speed_sound+(1-vapor_frac)*prop_w.speed_sound;
			p_prop->drdp=vapor_frac*prop_s.drdp+(1-vapor_frac)*prop_w.drdp;
			p_prop->vapor_fraction=vapor_frac;
			//solve vapor fraction and properties;
			return 0;
		}
		
		return -1;
	}	
}

int steam_prop_calc_ts(double temp, double spe_entr, steam_prop* p_prop)
{
	int rc,flag;
	double pres,pres_sat,vapor_frac;
	steam_prop prop_w,prop_s;
	
	if(temp>INTER_TEMP1)
	{
		flag=pres_ts_r2(&pres,temp,spe_entr);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r2(pres,temp,p_prop);
			if(rc==2&&fabs((p_prop->spe_entr-spe_entr)/spe_entr)<1e-4)
			{
				p_prop->spe_entr=spe_entr;
				return 0;
			}
		}
	
		flag=pres_ts_r3(&pres,temp,spe_entr);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r3(pres,temp,p_prop);
			if(rc==3&&fabs((p_prop->spe_entr-spe_entr)/spe_entr)<1e-4)
			{
				p_prop->temp=temp;
				p_prop->spe_entr=spe_entr;
				return 0;
			}
		}
		
		flag=pres_ts_r5(&pres,temp,spe_entr);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r5(pres,temp,p_prop);
			if(rc==5&&fabs((p_prop->spe_entr-spe_entr)/spe_entr)<1e-4)
			{
				p_prop->temp=temp;
				p_prop->spe_entr=spe_entr;
				return 0;
			}
		}
		return -1;	
	}
	
	else
	{
		flag=pres_ts_r1(&pres,temp,spe_entr);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r1(pres,temp,p_prop);
			if(rc==1&&fabs((p_prop->spe_entr-spe_entr)/spe_entr)<1e-4)
			{
				p_prop->spe_entr=spe_entr;
				return 0;
			}
		}
		flag=pres_ts_r2(&pres,temp,spe_entr);
		if(flag==0)
		{
			rc=region(pres,temp);
			steam_prop_calc_r2(pres,temp,p_prop);
			if(rc==2&&fabs((p_prop->spe_entr-spe_entr)/spe_entr)<1e-4)
			{
				p_prop->spe_entr=spe_entr;
				return 0;
			}
		}
		pres_sat=calc_sat_pres(temp);
		steam_prop_calc_r1(pres_sat,temp,&prop_w);
		steam_prop_calc_r2(pres_sat,temp,&prop_s);
		if(spe_entr>prop_w.spe_entr&&spe_entr<prop_s.spe_entr)
		{
			vapor_frac=(spe_entr-prop_w.spe_entr)/(prop_s.spe_entr-prop_w.spe_entr);
			p_prop->pressure=pres_sat;
			p_prop->temp=temp;
			p_prop->spe_vol=vapor_frac*prop_s.spe_vol+(1-vapor_frac)*prop_w.spe_vol;
			p_prop->dens=1/p_prop->spe_vol;
			p_prop->spe_energy=vapor_frac*prop_s.spe_energy+(1-vapor_frac)*prop_w.spe_energy;
			p_prop->spe_entr=spe_entr;
			p_prop->spe_enth=vapor_frac*prop_s.spe_enth+(1-vapor_frac)*prop_w.spe_enth;
			p_prop->spe_h_v=vapor_frac*prop_s.spe_h_v+(1-vapor_frac)*prop_w.spe_h_v;
			p_prop->spe_h_p=vapor_frac*prop_s.spe_h_p+(1-vapor_frac)*prop_w.spe_h_p;
			p_prop->speed_sound=vapor_frac*prop_s.speed_sound+(1-vapor_frac)*prop_w.speed_sound;
			p_prop->drdp=vapor_frac*prop_s.drdp+(1-vapor_frac)*prop_w.drdp;
			p_prop->vapor_fraction=vapor_frac;
		}
		
		return -1;
	}
}

int steam_prop_calc_hs(double spe_enth, double spe_entr, steam_prop* p_prop)
{
	int rc,flag;
	double pres,temp,temp_sat,p_sat,vapor_frac;
	steam_prop prop_w,prop_s;
	flag=pt_hs_r1(&pres,&temp,spe_enth,spe_entr);
	if(flag==0)
	{
//		printf("\n....................\n");
		rc=region(pres,temp);
		steam_prop_calc_r1(pres,temp,p_prop);
		if(rc==1&&fabs((p_prop->spe_entr-spe_entr)/spe_entr)<1e-4&&fabs((p_prop->spe_enth-spe_enth)/spe_enth)<1e-4)
		{
			p_prop->spe_entr=spe_entr;
			p_prop->spe_enth=spe_enth;
			return 0;
		}
	}
	
	flag=pt_hs_r2(&pres,&temp,spe_enth,spe_entr);
	if(flag==0)
	{
//		printf("\n....................%lf,%lf\n",pres,temp);
		rc=region(pres,temp);
		steam_prop_calc_r2(pres,temp,p_prop);
//		printf("\n\n%d,%lf,%lf,%lf,%lf,%lf\n",rc,pres,temp,p_prop->spe_enth,p_prop->spe_entr,p_prop->vapor_fraction);
		if(rc==2&&fabs((p_prop->spe_entr-spe_entr)/spe_entr)<1e-4&&fabs((p_prop->spe_enth-spe_enth)/spe_enth)<1e-4)
		{
			p_prop->spe_entr=spe_entr;
			p_prop->spe_enth=spe_enth;
//			printf("\nGREAT!\n");
			return 0;
		}
	}
	
	flag=pt_hs_r3(&pres,&temp,spe_enth,spe_entr);
	if(flag==0)
	{
//		printf("\n....................%lf,%lf\n",pres,temp);
		rc=region(pres,temp);
//		printf("@@@@@@@@@@@@@@@%d\n",rc);
		steam_prop_calc_r3(pres,temp,p_prop);
//		printf("\n\n%d,%lf,%lf,%lf,%lf,%lf\n",rc,pres,temp,p_prop->spe_enth,p_prop->spe_entr,p_prop->vapor_fraction);
		if(rc==3&&fabs((p_prop->spe_entr-spe_entr)/spe_entr)<1e-4&&fabs((p_prop->spe_enth-spe_enth)/spe_enth)<1e-4)
		{
			p_prop->spe_entr=spe_entr;
			p_prop->spe_enth=spe_enth;
			return 0;
		}
	}
	
	flag=pt_hs_r5(&pres,&temp,spe_enth,spe_entr);
	if(flag==0)
	{
//		printf("\n....................%lf,%lf\n",pres,temp);
		rc=region(pres,temp);
//		printf("@@@@@@@@@@@@@@@%d\n",rc);
		steam_prop_calc_r5(pres,temp,p_prop);
//		printf("\n\n%d,%lf,%lf,%lf,%lf,%lf\n",rc,pres,temp,p_prop->spe_enth,p_prop->spe_entr,p_prop->vapor_fraction);
		if(rc==5&&fabs((p_prop->spe_entr-spe_entr)/spe_entr)<1e-4&&fabs((p_prop->spe_enth-spe_enth)/spe_enth)<1e-4)
		{
			p_prop->spe_entr=spe_entr;
			p_prop->spe_enth=spe_enth;
			return 0;
		}
	}
	
	flag=hs_find_tsat(&temp_sat,spe_enth,spe_entr);
//	printf("\n%d,%lf\n",flag,temp_sat);
	if(flag==-1)
	{
		return -1;
	}

	p_sat=calc_sat_pres(temp_sat);
	steam_prop_calc_r1(p_sat,temp_sat,&prop_w);
	steam_prop_calc_r2(p_sat,temp_sat,&prop_s);
	vapor_frac=(spe_enth-prop_w.spe_enth)/(prop_s.spe_enth-prop_w.spe_enth);
	p_prop->pressure=p_sat;
	p_prop->temp=temp_sat;
	p_prop->spe_vol=vapor_frac*prop_s.spe_vol+(1-vapor_frac)*prop_w.spe_vol;
	p_prop->dens=1/p_prop->spe_vol;
	p_prop->spe_energy=vapor_frac*prop_s.spe_energy+(1-vapor_frac)*prop_w.spe_energy;
	p_prop->spe_entr=spe_entr;
	p_prop->spe_enth=spe_enth;
	p_prop->spe_h_v=vapor_frac*prop_s.spe_h_v+(1-vapor_frac)*prop_w.spe_h_v;
	p_prop->spe_h_p=vapor_frac*prop_s.spe_h_p+(1-vapor_frac)*prop_w.spe_h_p;
	p_prop->speed_sound=vapor_frac*prop_s.speed_sound+(1-vapor_frac)*prop_w.speed_sound;
	p_prop->drdp=vapor_frac*prop_s.drdp+(1-vapor_frac)*prop_w.drdp;
	p_prop->vapor_fraction=vapor_frac;
	return 0;	
}

int steam_prop_calc_px(double pres, double vf, steam_prop *p_prop)
{
	double presmn,presmx;
	double temp_sat;
	steam_prop prop_w,prop_s;
	
	if(vf<0||vf>1)
	{
		return -1;
	}
	presmn=calc_sat_pres(MIN_TEMP);
	presmx=calc_sat_pres(INTER_TEMP1);
	if(pres<presmn||pres>presmx)
	{
		return -1;
	}
	temp_sat=tsatp(pres);
	steam_prop_calc_r1(pres,temp_sat,&prop_w);
	steam_prop_calc_r2(pres,temp_sat,&prop_s);
	p_prop->pressure=pres;
	p_prop->temp=temp_sat;
	p_prop->spe_vol=vf*prop_s.spe_vol+(1-vf)*prop_w.spe_vol;
	p_prop->dens=1/p_prop->spe_vol;
	p_prop->spe_energy=vf*prop_s.spe_energy+(1-vf)*prop_w.spe_energy;
	p_prop->spe_entr=vf*prop_s.spe_entr+(1-vf)*prop_w.spe_entr;
	p_prop->spe_enth=vf*prop_s.spe_enth+(1-vf)*prop_w.spe_enth;
	p_prop->spe_h_v=vf*prop_s.spe_h_v+(1-vf)*prop_w.spe_h_v;
	p_prop->spe_h_p=vf*prop_s.spe_h_p+(1-vf)*prop_w.spe_h_p;
	p_prop->speed_sound=vf*prop_s.speed_sound+(1-vf)*prop_w.speed_sound;
	p_prop->drdp=vf*prop_s.drdp+(1-vf)*prop_w.drdp;
	p_prop->vapor_fraction=vf;
	return 0;	
}

int steam_prop_calc_tx(double temp, double vf, steam_prop *p_prop)
{
	double pres_sat;
	steam_prop prop_w,prop_s;
	
	if(vf<0||vf>1)
	{
		return -1;
	}
	if(temp<MIN_TEMP||temp>INTER_TEMP1)
	{
		return -1;
	}
	pres_sat=calc_sat_pres(temp);
	steam_prop_calc_r1(pres_sat,temp,&prop_w);
	steam_prop_calc_r2(pres_sat,temp,&prop_s);
	p_prop->pressure=pres_sat;
	p_prop->temp=temp;
	p_prop->spe_vol=vf*prop_s.spe_vol+(1-vf)*prop_w.spe_vol;
	p_prop->dens=1/p_prop->spe_vol;
	p_prop->spe_energy=vf*prop_s.spe_energy+(1-vf)*prop_w.spe_energy;
	p_prop->spe_entr=vf*prop_s.spe_entr+(1-vf)*prop_w.spe_entr;
	p_prop->spe_enth=vf*prop_s.spe_enth+(1-vf)*prop_w.spe_enth;
	p_prop->spe_h_v=vf*prop_s.spe_h_v+(1-vf)*prop_w.spe_h_v;
	p_prop->spe_h_p=vf*prop_s.spe_h_p+(1-vf)*prop_w.spe_h_p;
	p_prop->speed_sound=vf*prop_s.speed_sound+(1-vf)*prop_w.speed_sound;
	p_prop->drdp=vf*prop_s.drdp+(1-vf)*prop_w.drdp;
	p_prop->vapor_fraction=vf;
	return 0;	
}

double steam_visc_calc(double temp, double dens)
{
	double visc_0,visc_1;
	double sum1,sum2,sum_temp;
	double ref_t=647.096;
	double ref_dens=322;
	double ref_p=22.064e6;
	double ref_visc=1e-6;
	int i,j;
	
	sum1=0;
	for(i=0;i<4;i++)
	{
		sum1+=coeff_visc1[i]/pow(temp/ref_t,i);
	}
	visc_0=100*sqrt(temp/ref_t)/sum1;
	
	sum1=0;
	for(i=0;i<6;i++)
	{
		sum_temp=pow(ref_t/temp-1,i);
		sum2=0;
		for(j=0;j<7;j++)
		{
			sum2+=coeff_visc2[i][j]*pow(dens/ref_dens-1,j);
		}
		sum1=sum1+(sum_temp*sum2);
	}
	visc_1=exp(dens/ref_dens*sum1);
	return visc_0*visc_1*ref_visc;
}

int calc_index_lmd2(double dens, double ref_dens)
{
	double nd_dens=dens/ref_dens;
	if((nd_dens-0.310559006)<1e-8)
	{
		return 0;
	}
	else if((nd_dens-0.776397516)<1e-8)
	{
		return 1;
	}
	else if((nd_dens-1.242236025)<1e-8)
	{
		return 2;
	}
	else if((nd_dens-1.863354037)<1e-8)
	{
		return 3;
	}
	else 
	{
		return 4;
	}
}

double steam_thcond_calc(double temp, double dens, double visc, double cp, double cv, double drdp, double pressure)
{
	double lamd0,lamd1,lamd2;
	double ref_t=647.096;
	double ref_p=22.064e6;
	double ref_dens=322;
	double ref_lamd=1e-3;
	double ref_visc=1e-6;
	double gas_const=0.46151805e3;
	double nd_visc;
	double nd_cp, real_cp;
	double sp_ratio;
	double delta_x;
	double y,zy;
	
	double A=177.8514;
	double qd=2.5e9;
	double v=0.630;
	double gamma=1.239;
	double ypsilon0=0.13e-9;
	double gamma0=0.06;
	double tr=1.5;
	double ksai_t,ksai_tr,ypsilon;
	
	int indexj;
	
	int i,j;
	double sum1,sum2,sum3;
	sum1=0;
	for(i=0;i<5;i++)
	{
		sum1+=coeff_thcond1[i]/pow(temp/ref_t,i);	
	}
	lamd0=sqrt(temp/ref_t)/sum1;
//	printf("\t%lf\n",lamd0);
	
	sum1=0;
	for(i=0;i<5;i++)
	{
		sum2=pow(ref_t/temp-1,i);
		sum3=0;
		for(j=0;j<6;j++)
		{
			sum3+=coeff_thcond2[i][j]*pow(dens/ref_dens-1,j);
		}
		sum1=sum1+(sum2*sum3);
	}
	lamd1=exp(dens/ref_dens*sum1);
	
//	printf("\t%lf\n",lamd1);
//	printf("%lf,%lf,%lf\n",temp/ref_t,dens/ref_dens,visc/ref_visc);
	
	nd_visc=visc/ref_visc;
	nd_cp=cp/gas_const;
	if(nd_cp<0||nd_cp>1e13)
	{
		nd_cp=1e13;
	}
//	printf("\n%lf\n",nd_cp);
	real_cp=nd_cp*gas_const;
//	printf("\n%lf\n",real_cp);
	sp_ratio=real_cp/cv;
//	printf("\t\t%lf\n",sp_ratio);
	ksai_t=drdp*ref_p/ref_dens;
//	printf("\t%lf\n",ksai_t);
	indexj=calc_index_lmd2(dens,ref_dens);
//	printf("\n%d\n",indexj);
	sum1=0;
	for(i=0;i<6;i++)
	{
		sum1+=coeff_thcond4[i][indexj]*pow(dens/ref_dens,i);
	}
	ksai_tr=1/sum1;
//	printf("\t\t%lf\n\n",ksai_tr);
	delta_x=(ksai_t-ksai_tr*tr*ref_t/temp)*dens/ref_dens;
	if(delta_x<0)
	{
		delta_x=0;
	}
	ypsilon=ypsilon0*pow(delta_x/gamma0,v/gamma);
	if(ypsilon<0||ypsilon>1e4)
	{
		ypsilon=1e4;
	}
	y=qd*ypsilon;
	zy=(((1-1/sp_ratio)*atan(y)+y/sp_ratio)-(1-exp(-1/(1/y+y*y/(3*dens*dens/(ref_dens*ref_dens))))))*2/(3.14159265358*y);
	if(y<1.2e-7)
	{
		zy=0;
	}
	lamd2=A*zy*(dens/ref_dens)*nd_cp*(temp/ref_t)/(visc/ref_visc);
//	printf("\t%lf\n",lamd2);
	return (lamd0*lamd1+lamd2)*ref_lamd;
}


//special for THE PROGRAM
int steam_prop_calc(double pressure, double temp, double *dens,double *cp, double* cv, double* spe_enth, double* drdp)
{
	steam_prop p_prop;
	int region_code;
	int calc_flag;
	
	region_code=region(pressure,temp);
	if(region_code==-1)
	{
		printf("\n%lf,%lfFATAL ERROR: PRESSURE/TEMP EXCEEDS BOUNDARIES.\n",pressure,temp);
		return -1;
	}
	else if(region_code==1)
	{
		steam_prop_calc_r1(pressure,temp,&p_prop);
	}
	else if(region_code==2||region_code==4)
	{
		steam_prop_calc_r2(pressure,temp,&p_prop);
	}
	else if(region_code==3)
	{
		calc_flag=steam_prop_calc_r3(pressure,temp,&p_prop);
		if(calc_flag==-1)
		{
			printf("\nFATAL ERROR: DENSITY ITERATION NOT CONVERGED.\n");
			return -2;
		}
	}
	else if(region_code==5)
	{
		steam_prop_calc_r5(pressure,temp,&p_prop);
	}
//	print_prop(&p_prop);
	*dens=p_prop.dens;
	*cp=p_prop.spe_h_p;
	*cv=p_prop.spe_h_v;
	*spe_enth=p_prop.spe_enth;
	*drdp=p_prop.drdp;
	return 0;
}

/*int main()
{
	double h,s,viscous,thcond,p,t;
	steam_prop prop,propb;
	int flag;
	FILE *fout=fopen("ztest2425.dat","w");
	FILE *fout2=fopen("ztestpt.dat","w"); 
	int i=0;
	
	steam_prop_calc_pt(0.0035e6,300,&prop,&propb);
	printf("#%lf,%Lf,%lf,%lf\n",prop.pressure,prop.temp,prop.dens,prop.spe_entr);
/*	for(s=4.5e3;s<10e3;s=s+100)
	{
		for(h=2.4e6;h<2.5e6;h=h+1e4)
		{
			flag=steam_prop_calc_hs(h,s,&prop);
			if(flag==-1)
			{
				i++;
				printf("\n! WARNING: Calculation error.\n! Calculation abort at line #%d of the input file.\n! Please press any key to exit.",i);
				fprintf(fout,"%lf\t%lf,-1\n",h,s);
			}
			else
			{
				viscous=steam_visc_calc(prop.temp,prop.dens);
				thcond=steam_thcond_calc(prop.temp,prop.dens,viscous,prop.spe_h_p,prop.spe_h_v,prop.drdp,prop.pressure);
				fprintf(fout,"%.6lf\t%.2lf\t%.8lf\t%.6lf\t%.4lf\t%.4lf\t%.8lf\t%.6lf\t%.6lf\t%.4lf\t%.8lf\t%.4lf\t%.4lf\n",prop.pressure/1e6,prop.temp,prop.spe_vol,prop.dens,prop.spe_energy/1000,prop.spe_enth/1000,prop.spe_entr/1000,prop.spe_h_v/1000,prop.spe_h_p/1000,prop.speed_sound,viscous,thcond,prop.vapor_fraction);
			}	
		}
	}*/ 
	
/*	for(t=MIN_TEMP+0.01;t<INTER_TEMP1+0.02;t=t+10)
	{
		for(p=40e6;p<100e6+101;p=p+1e5)
		{
			flag=steam_prop_calc_pt(p,t,&prop,&propb);
			if(flag==-1)
			{
				i++;
				printf("\n! WARNING: Calculation error.\n! Calculation abort at line #%d of the input file.\n! Please press any key to exit.",i);
				fprintf(fout2,"%lf\t%lf,-1\n",p,t);
			}
			else
			{
				viscous=steam_visc_calc(prop.temp,prop.dens);
				thcond=steam_thcond_calc(prop.temp,prop.dens,viscous,prop.spe_h_p,prop.spe_h_v,prop.drdp,prop.pressure);
				fprintf(fout2,"%.6lf\t%.2lf\t%.8lf\t%.6lf\t%.4lf\t%.4lf\t%.8lf\t%.6lf\t%.6lf\t%.4lf\t%.8lf\t%.4lf\t%.4lf\n",prop.pressure/1e6,prop.temp,prop.spe_vol,prop.dens,prop.spe_energy/1000,prop.spe_enth/1000,prop.spe_entr/1000,prop.spe_h_v/1000,prop.spe_h_p/1000,prop.speed_sound,viscous,thcond,prop.vapor_fraction);
			}	
		}
	}
	fclose(fout);
	fclose(fout);
	return 0;
}


/*int main()
{
	
	int flag;
	steam_prop prop,prop2;
//	printf(".......................%lf\n",calc_sat_pres(450));
//	flag=steam_prop_calc_px(1e6,0.98,&prop);
//	printf("%lf,%lf,%lf,%lf\n",prop.spe_vol,prop.spe_entr,prop.spe_enth,prop.vapor_fraction);
//		flag=steam_prop_calc_tx(450,0.98,&prop);
//	printf("%lf,%lf,%lf,%lf\n",prop.spe_vol,prop.spe_entr,prop.spe_enth,prop.vapor_fraction);
//	printf("\n;;;;;; %lf\n",calc_bound23(659.961574204728));
//	printf("\n;;;;;; %lf\n",calc_bound23(623.15));
	flag=steam_prop_calc_hs(2600000,5200,&prop);
	
	if(flag==-1)
	{
		printf("\nEEEEEEEEEEEEEEEEEEEEE!\n");
	}
	else
	{
		printf("%lf,%lf,%lf\n",prop.pressure,prop.temp,prop.spe_enth);
	}
}
/*int main()
{
	int i,j;
	steam_prop prop;
	double pres,temp,viscous,thcond;
	FILE* output=fopen("test_output_pt.dat","w");
	fprintf(output,"\tPRES\t\tTEMP\tSPE_VOL\t\tDENS\t\tu\t\th\t\ts\t\tCv\t\tCp\t\tVsound\t\tVisc\t\tLamd\n");
	fprintf(output,"\tMPa\t\tK\tm^3/kg\t\tkg/m^3\t\tkJ/kg\t\tkJ/kg\t\tkJ/(kg*K)\tkJ/(kg*K)\tkJ/(kg*K)\tm/s\t\tPa.s\t\tW/(m*K)\n");

	for(i=0;i<100;i++)
	{
		for(j=0;j<200;j++)
		{
			pres=1+i*1e6;
			temp=273.16+j*(1073.15-273.15)/200;
			steam_prop_calc_pt(pres,temp,&prop);
			viscous=steam_visc_calc(prop.temp,prop.dens);
			thcond=steam_thcond_calc(prop.temp,prop.dens,viscous,prop.spe_h_p,prop.spe_h_v,prop.drdp,prop.pressure);
			fprintf(output,"%d\t%.6lf\t%.2lf\t%.8lf\t%.6lf\t%.4lf\t%.4lf\t%.8lf\t%.6lf\t%.6lf\t%.4lf\t%.8lf\t%.4lf\n",i,prop.pressure/1e6,prop.temp,prop.spe_vol,prop.dens,prop.spe_energy/1000,prop.spe_enth/1000,prop.spe_entr/1000,prop.spe_h_v/1000,prop.spe_h_p/1000,prop.speed_sound,viscous,thcond);
		}
	}
	return 0;
} 


/*int main()
{
	steam_prop prop;
	double viscous;
	
	printf("PRESSURE\tTEMP\tVISCOUS\t\tTHCOND\n");
	printf("Pa\t\tK\tPa.s\t\tW/(m.K)\n");
	
	steam_prop_calc_orig(20e6,620,&prop);
	viscous=steam_visc_calc(620,prop.dens);
	printf("%.2lf\t%.2lf\t%.9lf\t%.9lf\n",20e6,620.0,viscous,steam_thcond_calc(620,prop.dens,viscous,prop.spe_h_p,prop.spe_h_v,prop.drdp,20e6));

	steam_prop_calc_orig(50e6,620,&prop);
	viscous=steam_visc_calc(620,prop.dens);
	printf("%.2lf\t%.2lf\t%.9lf\t%.9lf\n",50e6,620.0,viscous,steam_thcond_calc(620,prop.dens,viscous,prop.spe_h_p,prop.spe_h_v,prop.drdp,50e6));

	steam_prop_calc_orig(0.3e6,650,&prop);
	viscous=steam_visc_calc(650,prop.dens);
	printf("%.2lf\t%.2lf\t%.9lf\t%.9lf\n",0.3e6,650.0,viscous,steam_thcond_calc(650,prop.dens,viscous,prop.spe_h_p,prop.spe_h_v,prop.drdp,0.3e6));

	steam_prop_calc_orig(50e6,800,&prop);
	viscous=steam_visc_calc(800,prop.dens);
	printf("%.2lf\t%.2lf\t%.9lf\t%.9lf\n",50e6,800.0,viscous,steam_thcond_calc(800,prop.dens,viscous,prop.spe_h_p,prop.spe_h_v,prop.drdp,50e6));

	steam_prop_calc_orig(21.984062713427e6,647.35,&prop);
	viscous=steam_visc_calc(647.35,prop.dens);
	printf("%.2lf\t%.2lf\t%.9lf\t%.9lf\n",21.984062713427e6,647.35,viscous,steam_thcond_calc(647.35,prop.dens,viscous,prop.spe_h_p,prop.spe_h_v,prop.drdp,21.984062713427e6));

	steam_prop_calc_orig(22.132160017547e6,647.35,&prop);
	viscous=steam_visc_calc(647.35,prop.dens);
	printf("%.2lf\t%.2lf\t%.9lf\t%.9lf\n",22.132160017547e6,647.35,viscous,steam_thcond_calc(647.35,prop.dens,viscous,prop.spe_h_p,prop.spe_h_v,prop.drdp,22.132160017547e6));

}


/*int main()
{
	double dens,n1,n2;
	int i;
	for(i=1;i<50;i++)
	{
		steam_prop_calc(101325*i,300,&dens,&n1,&n2);
		n1=101325*i;
		printf("%lf,%lf\n",n1,dens);
	}
	

}
	
/*	double temp[4]={
		300,500,700,900
	};
	double pres[9]={
		100,10000,1e6,3e6,5e6,10e6,20e6,50e6,90e6
	};
	steam_prop prop;
	double pres_cond, temp_cond;
	int i,j;
	
	steam_temp_calc_ph(&pres_cond,0.5e6,0.521976332e7);
	printf("%.4lf\n",pres_cond);
	steam_temp_calc_ph(&pres_cond,8e6,0.520609634e7);
	printf("%.4lf\n",pres_cond);
	steam_temp_calc_ph(&pres_cond,8e6,0.658380291e7);
	printf("%.4lf\n",pres_cond);
	
	steam_temp_calc_ph(&pres_cond,0.255837018e8,0.186343019e7);
	printf("%.4lf\n",pres_cond);
//	temp_ph_r3(&pres_cond,0.222930643e8,0.237512401e7);
//	printf("%.4lf\n",pres_cond);
	steam_temp_calc_ph(&pres_cond,0.222930643e8,0.237512401e7);
	printf("%.4lf\n",pres_cond);
	steam_temp_calc_ph(&pres_cond,0.783095639e8,0.225868845e7);
	printf("%.4lf\n",pres_cond);
	
	
	steam_temp_calc_ph(&pres_cond,3e6,4e6);
	printf("%.4lf\n",pres_cond);
	steam_temp_calc_ph(&pres_cond,25e6,3.5e6);
	printf("%.4lf\n",pres_cond);
	steam_temp_calc_ph(&pres_cond,60e6,2.7e6);
	printf("%.4lf\n",pres_cond);
	
	steam_temp_calc_ph(&pres_cond,3e6,5e5);
	printf("%.4lf\n",pres_cond);
	steam_temp_calc_ph(&pres_cond,80e6,5e5);
	printf("%.4lf\n",pres_cond);
	steam_temp_calc_ph(&pres_cond,80e6,15e5);
	printf("%.4lf\n",pres_cond);
/*	FILE *p_file=fopen("_test_SPC.dat","w+");
//	fprintf(p_file,"P\\T\t\t");
//	for(i=0;i<4;i++)
//	{
//		fprintf(p_file,"%.4lf\t",temp[i]);
//	}
//	fprintf(p_file,"\n");
	for(i=0;i<9;i++)
	{
		pres_cond=pres[i];
//		fprintf(p_file,"%.2e\t",pres_cond);
		for(j=0;j<4;j++)
		{
			temp_cond=temp[j];			
			steam_prop_calc(pres_cond,temp_cond,&prop);
			fprintf(p_file,"%.6lf\t%.6lf\n",prop.temp,prop.spe_energy/prop.spe_h_v);
		}
	}
	fclose(p_file);
}*/



/*
*/