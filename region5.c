#include<stdio.h>
#include<math.h>
#include "steam_property_calc.h"
#define MAX_ITER_TIMES_R5 200

static double region5_coeff0[6][3]={
	1,0,-0.13179983674201e2,
	2,1,0.68540841634434e1,
	3,-3,-0.24805148933466e-1,
	4,-2,0.36901534980333,
	5,-1,-0.31161318213925e1,
	6,2,-0.32961626538917
};

static double region5_coeffr[6][4]={
	1,1,1,0.15736404855259e-2,
	2,1,2,0.90153761673944e-3,
	3,1,3,-0.50270077677648e-2,
	4,2,3,0.22440037409485e-5,
	5,2,9,-0.41163275453471e-5,
	6,3,7,0.37919454822955e-7
};

double gm05(double pi, double tau)
{
	int i;
	double ni,ji;
	double res=log(pi);
	for(i=0;i<6;i++)
	{
		ni=region5_coeff0[i][2];
		ji=region5_coeff0[i][1];
		res+=ni*pow(tau,ji);
	}
	return res;
}

double gmr5(double pi, double tau)
{
	int i;
	double ni,ii,ji;
	double res=0;
	for(i=0;i<6;i++)
	{
		ni=region5_coeffr[i][3];
		ji=region5_coeffr[i][2];
		ii=region5_coeffr[i][1];
		res+=ni*pow(pi,ii)*pow(tau,ji);
	}
	return res;
}

double gm0pi5(double pi)
{
	return 1/pi;
}

double gm0pipi5(double pi)
{
	return -1/(pi*pi);
}

double gm0tau5(double tau)
{
	int i;
	double res=0;
	double ni,ji;
	for(i=0;i<6;i++)
	{
		ni=region5_coeff0[i][2];
		ji=region5_coeff0[i][1];
		res+=ni*ji*pow(tau,ji-1);
	}
	return res;
}

double gm0tautau5(double tau)
{
	int i;
	double res=0;
	double ni,ji;
	for(i=0;i<6;i++)
	{
		ni=region5_coeff0[i][2];
		ji=region5_coeff0[i][1];
		res+=ni*ji*(ji-1)*pow(tau,ji-2);
	}
	return res;
}

double gmrpi5(double pi, double tau)
{
	int i;
	double ni,ii,ji;
	double res=0;
	for(i=0;i<6;i++)
	{
		ni=region5_coeffr[i][3];
		ji=region5_coeffr[i][2];
		ii=region5_coeffr[i][1];
		res+=ni*ii*pow(pi,ii-1)*pow(tau,ji);
	}
	return res;
}

double gmrpipi5(double pi, double tau)
{
	int i;
	double ni,ii,ji;
	double res=0;
	for(i=0;i<6;i++)
	{
		ni=region5_coeffr[i][3];
		ji=region5_coeffr[i][2];
		ii=region5_coeffr[i][1];
		res+=ni*ii*(ii-1)*pow(pi,ii-2)*pow(tau,ji);
	}
	return res;
}

double gmrtau5(double pi, double tau)
{
	int i;
	double ni,ii,ji;
	double res=0;
	for(i=0;i<6;i++)
	{
		ni=region5_coeffr[i][3];
		ji=region5_coeffr[i][2];
		ii=region5_coeffr[i][1];
		res+=ni*pow(pi,ii)*ji*pow(tau,ji-1);
	}
	return res;
}

double gmrtautau5(double pi, double tau)
{
	int i;
	double ni,ii,ji;
	double res=0;
	for(i=0;i<6;i++)
	{
		ni=region5_coeffr[i][3];
		ji=region5_coeffr[i][2];
		ii=region5_coeffr[i][1];
		res+=ni*pow(pi,ii)*ji*(ji-1)*pow(tau,ji-2);
	}
	return res;
}

double gmrpitau5(double pi, double tau)
{
	int i;
	double ni,ii,ji;
	double res=0;
	for(i=0;i<6;i++)
	{
		ni=region5_coeffr[i][3];
		ji=region5_coeffr[i][2];
		ii=region5_coeffr[i][1];
		res+=ni*ii*pow(pi,ii-1)*ji*pow(tau,ji-1);
	}
	return res;
}


void steam_prop_calc_r5(double pressure, double temp, steam_prop* p_prop)
{
	double pi,tau;
	
	pi=pressure/1e6;
	tau=1000/temp;
	
	p_prop->pressure=pressure;
	p_prop->temp=temp;
	p_prop->spe_vol=pi*(gm0pi5(pi)+gmrpi5(pi,tau))*GAS_CONST_STEAM*temp/pressure;
	p_prop->dens=1/p_prop->spe_vol;
	p_prop->spe_energy=GAS_CONST_STEAM*temp*(tau*(gm0tau5(tau)+gmrtau5(pi,tau))-pi*(gm0pi5(pi)+gmrpi5(pi,tau)));
	p_prop->spe_entr=GAS_CONST_STEAM*(tau*(gm0tau5(tau)+gmrtau5(pi,tau))-gm05(pi,tau)-gmr5(pi,tau));
	p_prop->spe_enth=GAS_CONST_STEAM*temp*tau*(gm0tau5(tau)+gmrtau5(pi,tau));
	p_prop->spe_h_v=GAS_CONST_STEAM*(-tau*tau*(gm0tautau5(tau)+gmrtautau5(pi,tau))+pow(1+pi*gmrpi5(pi,tau)-tau*pi*gmrpitau5(pi,tau),2)/(1-pi*pi*gmrpipi5(pi,tau)));
	p_prop->spe_h_p=-GAS_CONST_STEAM*tau*tau*(gm0tautau5(tau)+gmrtautau5(pi,tau));
	p_prop->speed_sound=sqrt(GAS_CONST_STEAM*temp*(1+2*pi*gmrpi5(pi,tau)+pi*pi*gmrpi5(pi,tau)*gmrpi5(pi,tau))/(1-pi*pi*gmrpipi5(pi,tau)+pow(1+pi*gmrpi5(pi,tau)-tau*pi*gmrpitau5(pi,tau),2)/(tau*tau*(gm0tautau5(tau)+gmrtautau5(pi,tau)))));
	p_prop->drdp=-p_prop->dens*p_prop->dens*GAS_CONST_STEAM*temp*(gm0pipi5(pi)+gmrpipi5(pi,tau))/(1e12);
	p_prop->vapor_fraction=1;
}

void calcFT_dFT_prr5(double* ft, double* dft, double pressure, double dens, double tau)
{
	double pi;
	pi=pressure/1e6;
	*ft=gm0pi5(pi)+gmrpi5(pi,tau)-1e6*tau/(dens*GAS_CONST_STEAM*1000);
	*dft=gmrpitau5(pi,tau)-1e6/(dens*GAS_CONST_STEAM*1000);
}

int temp_pr_r5(double *temp, double pressure, double dens)
{
	double tau_prev,tau_nxt,ft,dft;
	int iter_times;
	int ti=0;
	double tau_ini[12]={
		1,0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.45
	};
	
	do{
		iter_times=0;
		calcFT_dFT_prr5(&ft,&dft,pressure,dens,tau_ini[ti]);
		tau_nxt=tau_ini[ti]-ft/dft;
		do{
			tau_prev=tau_nxt;
			calcFT_dFT_prr5(&ft,&dft,pressure,dens,tau_prev);
			tau_nxt=tau_prev-ft/dft;
			iter_times++;
		}while(fabs(tau_nxt-tau_prev)>1e-6&&iter_times<MAX_ITER_TIMES_R5);
		ti++;
		if(iter_times==MAX_ITER_TIMES_R5)
		{
			continue;
		}
		else
		{
			*temp=1000/tau_nxt;
			if(NRT(*temp)==0) 
			break;
		}			
	}while(ti<12);
	if(ti==12)
	{
		*temp=-1;
		return -1;
	}
	else
	{
		return 0;
	}	
} 

void calcFT_dFT_pur5(double* ft, double* dft, double pressure, double spe_ener, double tau)
{
	double pi;
	pi=pressure/1e6;
	
	*ft=tau*(gm0tau5(tau)+gmrtau5(pi,tau))-pi*(gm0pi5(pi)+gmrpi5(pi,tau))-spe_ener*tau/(GAS_CONST_STEAM*1000);
	*dft=gm0tau5(tau)+gmrtau5(pi,tau)+tau*(gm0tautau5(tau)+gmrtautau5(pi,tau))-pi*(gmrpitau5(pi,tau))-spe_ener/(GAS_CONST_STEAM*1000);
}

int temp_pu_r5(double* temp, double pressure, double spe_ener)
{
	double tau_prev,tau_nxt,ft,dft;
	int iter_times;
	int ti=0;
	double tau_ini[12]={
		1,0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.45
	};
	do{
		iter_times=0;
		calcFT_dFT_pur5(&ft,&dft,pressure,spe_ener,tau_ini[ti]);
//	printf("##:%lf,%lf\n",ft,dft);
		tau_nxt=tau_ini[ti]-ft/dft;
		do{
			tau_prev=tau_nxt;
			calcFT_dFT_pur5(&ft,&dft,pressure,spe_ener,tau_prev);
			tau_nxt=tau_prev-ft/dft;
			iter_times++;
//		printf("## %d,%lf,%lf\n",iter_times,ft/dft,540/tau_nxt); 
		}while(fabs(tau_nxt-tau_prev)>1e-6&&iter_times<MAX_ITER_TIMES_R5);
		ti++;
		if(iter_times==MAX_ITER_TIMES_R5)
		{
			continue;
		}
		else
		{
			*temp=1000/tau_nxt;
			if(NRT(*temp)==0) 
		 	break;
		}
	}while(ti<12&&NRT(*temp)!=0);
	if(ti==12)
	{
		*temp=-1;
		return -1;
	}
	else
	{
		return 0;
	}
}


void calcFT_dFT_phr5(double* ft, double* dft, double pressure, double spe_enth, double tau)
{
	double pi;
	pi=pressure/1e6;
	*ft=gm0tau5(tau)+gmrtau5(pi,tau)-spe_enth/(GAS_CONST_STEAM*1000);
	*dft=gm0tautau5(tau)+gmrtautau5(pi,tau);
}

int temp_ph_r5(double* temp, double pressure, double spe_enth)
{
	double tau_prev,tau_nxt,ft,dft;
	int iter_times;
	int ti=0;
	double tau_ini[12]={
		1,0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.45
	};
	do{
		iter_times=0;
		calcFT_dFT_phr5(&ft,&dft,pressure,spe_enth,tau_ini[ti]);
//	printf("##:%lf,%lf\n",ft,dft);
		tau_nxt=tau_ini[ti]-ft/dft;
		do{
			tau_prev=tau_nxt;
			calcFT_dFT_phr5(&ft,&dft,pressure,spe_enth,tau_prev);
			tau_nxt=tau_prev-ft/dft;
			iter_times++;
//		printf("## %d,%lf,%lf\n",iter_times,ft/dft,540/tau_nxt); 
		}while(fabs(tau_nxt-tau_prev)>1e-6&&iter_times<MAX_ITER_TIMES_R5);
		ti++;
		if(iter_times==MAX_ITER_TIMES_R5)
		{
			continue;
		}
		else
		{
			*temp=1000/tau_nxt;
			if(NRT(*temp)==0) 
		 	break;
		}
	}while(ti<12&&NRT(*temp)!=0);
	if(ti==12)
	{
		*temp=-1;
		return -1;
	}
	else
	{
		return 0;
	}
}

void calcFT_dFT_psr5(double* ft, double* dft, double pressure, double spe_entr, double tau)
{
	double pi;
	pi=pressure/1e6;
	*ft=tau*(gm0tau5(tau)+gmrtau5(pi,tau))-(gm05(pi,tau)+gmr5(pi,tau))-spe_entr/GAS_CONST_STEAM;
	*dft=tau*(gm0tautau5(tau)+gmrtautau5(pi,tau));
}

int temp_ps_r5(double* temp, double pressure, double spe_entr)
{
	double tau_prev,tau_nxt,ft,dft;
	int iter_times;
	int ti=0;
	double tau_ini[12]={
		1,0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.45
	};
	do{
		iter_times=0;
		calcFT_dFT_psr5(&ft,&dft,pressure,spe_entr,tau_ini[ti]);
//	printf("##:%lf,%lf\n",ft,dft);
		tau_nxt=tau_ini[ti]-ft/dft;
		do{
			tau_prev=tau_nxt;
			calcFT_dFT_psr5(&ft,&dft,pressure,spe_entr,tau_prev);
			tau_nxt=tau_prev-ft/dft;
			iter_times++;
//		printf("## %d,%lf,%lf\n",iter_times,ft/dft,540/tau_nxt); 
		}while(fabs(tau_nxt-tau_prev)>1e-6&&iter_times<MAX_ITER_TIMES_R5);
		ti++;
		if(iter_times==MAX_ITER_TIMES_R5)
		{
			continue;
		}
		else
		{
			*temp=1000/tau_nxt;
			if(NRT(*temp)==0) 
		 	break;
		}
	}while(ti<12&&NRT(*temp)!=0);
	if(ti==12)
	{
		*temp=-1;
		return -1;
	}
	else
	{
		return 0;
	}
}

void calcFP_dFP_trr5(double* fp, double* dfp, double temp, double dens, double pi)
{
	double tau;
	tau=1000/temp;
	*fp=gm0pi5(pi)+gmrpi5(pi,tau)-1e6/(dens*GAS_CONST_STEAM*temp);
	*dfp=gm0pipi5(pi)+gmrpipi5(pi,tau);
}

int pres_tr_r5(double* pres, double temp, double dens)
{
	double pi_prev,pi_nxt,fp,dfp;
	int iter_times;
	int ti=0;
	double pi_ini[18]={
		25,20,30,15,35,10,40,5,45,1,49,0.5,49.5,0.05,0.005,0.0005,49.9,50
	};
	do{
		iter_times=0;
		calcFP_dFP_trr5(&fp,&dfp,temp,dens,pi_ini[ti]);
//	printf("##:%lf,%lf\n",ft,dft);
		pi_nxt=pi_ini[ti]-fp/dfp;
		do{
			pi_prev=pi_nxt;
			calcFP_dFP_trr5(&fp,&dfp,temp,dens,pi_prev);
			pi_nxt=pi_prev-fp/dfp;
			iter_times++;
//		printf("## %d,%lf,%lf\n",iter_times,ft/dft,540/tau_nxt); 
		}while(fabs(pi_nxt-pi_prev)>1e-6&&iter_times<MAX_ITER_TIMES_R5);
		ti++;
		if(iter_times==MAX_ITER_TIMES_R5)
		{
			continue;
		}
		else
		{
			*pres=pi_nxt*1e6;
			if(NRP(*pres)==0) 
		 	break;
		}
	}while(ti<18&&NRP(*pres)!=0);
	if(ti==18)
	{
		*pres=-1;
		return -1;
	}
	else
	{
		return 0;
	}
}

void calcFP_dFP_tur5(double* fp, double* dfp, double temp, double spe_ener, double pi)
{
	double tau;
	tau=1000/temp;
	*fp=tau*(gm0tau5(tau)+gmrtau5(pi,tau))-pi*(gm0pi5(pi)+gmrpi5(pi,tau))-spe_ener*tau/(GAS_CONST_STEAM*1000);
	*dfp=tau*gmrpitau5(pi,tau)-(gm0pi5(pi)+gmrpi5(pi,tau))-pi*(gm0pipi5(pi)+gmrpipi5(pi,tau));
}

int pres_tu_r5(double* pres, double temp, double spe_ener)
{
	double pi_prev,pi_nxt,fp,dfp;
	int iter_times;
	int ti=0;
	double pi_ini[18]={
		25,20,30,15,35,10,40,5,45,1,49,0.5,49.5,0.05,0.005,0.0005,49.9,50
	};
	do{
		iter_times=0;
		calcFP_dFP_tur5(&fp,&dfp,temp,spe_ener,pi_ini[ti]);
//	printf("##:%lf,%lf\n",ft,dft);
		pi_nxt=pi_ini[ti]-fp/dfp;
		do{
			pi_prev=pi_nxt;
			calcFP_dFP_tur5(&fp,&dfp,temp,spe_ener,pi_prev);
			pi_nxt=pi_prev-fp/dfp;
			iter_times++;
//		printf("## %d,%lf,%lf\n",iter_times,ft/dft,540/tau_nxt); 
		}while(fabs(pi_nxt-pi_prev)>1e-6&&iter_times<MAX_ITER_TIMES_R5);
		ti++;
		if(iter_times==MAX_ITER_TIMES_R5)
		{
			continue;
		}
		else
		{
			*pres=pi_nxt*1e6;
			if(NRP(*pres)==0) 
		 	break;
		}
	}while(ti<18&&NRP(*pres)!=0);
	if(ti==18)
	{
		*pres=-1;
		return -1;
	}
	else
	{
		return 0;
	}
}

void calcFP_dFP_thr5(double* fp, double* dfp, double temp, double spe_enth, double pi)
{
	double tau;
	tau=1000/temp;
	*fp=gm0tau5(tau)+gmrtau5(pi,tau)-spe_enth/(GAS_CONST_STEAM*1000);
	*dfp=gmrpitau5(pi,tau);
}

int pres_th_r5(double* pres, double temp, double spe_enth)
{
	double pi_prev,pi_nxt,fp,dfp;
	int iter_times;
	int ti=0;
	double pi_ini[18]={
		25,20,30,15,35,10,40,5,45,1,49,0.5,49.5,0.05,0.005,0.0005,49.9,50
	};
	do{
		iter_times=0;
		calcFP_dFP_thr5(&fp,&dfp,temp,spe_enth,pi_ini[ti]);
//	printf("##:%lf,%lf\n",ft,dft);
		pi_nxt=pi_ini[ti]-fp/dfp;
		do{
			pi_prev=pi_nxt;
			calcFP_dFP_thr5(&fp,&dfp,temp,spe_enth,pi_prev);
			pi_nxt=pi_prev-fp/dfp;
			iter_times++;
//		printf("## %d,%lf,%lf\n",iter_times,ft/dft,540/tau_nxt); 
		}while(fabs(pi_nxt-pi_prev)>1e-6&&iter_times<MAX_ITER_TIMES_R5);
		ti++;
		if(iter_times==MAX_ITER_TIMES_R5)
		{
			continue;
		}
		else
		{
			*pres=pi_nxt*1e6;
			if(NRP(*pres)==0) 
		 	break;
		}
	}while(ti<18&&NRP(*pres)!=0);
	if(ti==18)
	{
		*pres=-1;
		return -1;
	}
	else
	{
		return 0;
	}
}

void calcFP_dFP_tsr5(double* fp, double* dfp, double temp, double spe_entr, double pi)
{
	double tau;
	tau=1000/temp;
	*fp=tau*(gm0tau5(tau)+gmrtau5(pi,tau))-(gm05(pi,tau)+gmr5(pi,tau))-spe_entr/GAS_CONST_STEAM;
	*dfp=tau*gmrpitau5(pi,tau)-gm0pi5(pi)-gmrpi5(pi,tau);
}

int pres_ts_r5(double* pres, double temp, double spe_entr)
{
	double pi_prev,pi_nxt,fp,dfp;
	int iter_times;
	int ti=0;
	double pi_ini[18]={
		25,20,30,15,35,10,40,5,45,1,49,0.5,49.5,0.05,0.005,0.0005,49.9,50
	};
	do{
		iter_times=0;
		calcFP_dFP_tsr5(&fp,&dfp,temp,spe_entr,pi_ini[ti]);
//	printf("##:%lf,%lf\n",ft,dft);
		pi_nxt=pi_ini[ti]-fp/dfp;
		do{
			pi_prev=pi_nxt;
			calcFP_dFP_tsr5(&fp,&dfp,temp,spe_entr,pi_prev);
			pi_nxt=pi_prev-fp/dfp;
			iter_times++;
//		printf("## %d,%lf,%lf\n",iter_times,ft/dft,540/tau_nxt); 
		}while(fabs(pi_nxt-pi_prev)>1e-6&&iter_times<MAX_ITER_TIMES_R5);
		ti++;
		if(iter_times==MAX_ITER_TIMES_R5)
		{
			continue;
		}
		else
		{
			*pres=pi_nxt*1e6;
			if(NRP(*pres)==0) 
		 	break;
		}
	}while(ti<18&&NRP(*pres)!=0);
	if(ti==18)
	{
		*pres=-1;
		return -1;
	}
	else
	{
		return 0;
	}
}

void calcFG_dFG_hs_r5(double* f, double* g, double*dfdp, double* dfdt, double* dgdp, double* dgdt, double spec_enth, double spec_entr, double pi, double tau)
{
	*f=gm0tau5(tau)+gmrtau5(pi,tau)-spec_enth/(GAS_CONST_STEAM*1000);
	*g=tau*(gm0tau5(tau)+gmrtau5(pi,tau))-(gm05(pi,tau)+gmr5(pi,tau))-spec_entr/GAS_CONST_STEAM;
	*dfdp=gmrpitau5(pi,tau);
	*dfdt=gm0tautau5(tau)+gmrtautau5(pi,tau);
	*dgdp=tau*gmrpitau5(pi,tau)-gm0pi5(pi)-gmrpi5(pi,tau);
	*dgdt=tau*gm0tautau5(tau)+tau*gmrtautau5(pi,tau);
}

int pt_hs_r5(double* pres, double* temp, double spe_enth, double spe_entr)
{
	double pi_prev,pi_nxt;
	double tau_prev,tau_nxt,pp,tt;
	int i=0;
	double f,g,dfdt,dgdt,dfdp,dgdp;
	int p_index,t_index;
	int flag=0;
	double dp,dt,maxd;
	steam_prop prop;
	double pi_ini[16]={
		25,20,30,15,35,10,40,5,45,1,49,0.5,49.5,0.05,0.005,0.0005
	};
	double tau_ini[12]={
		1,0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.45
	};
	for(p_index=0;p_index<16;p_index++)
	{
		for(t_index=0;t_index<12;t_index++)
		{
			printf("\t%lf,%lf\n",pi_ini[p_index],tau_ini[t_index]);
			calcFG_dFG_hs_r5(&f,&g,&dfdp,&dfdt,&dgdp,&dgdt,spe_enth,spe_entr,pi_ini[p_index],tau_ini[t_index]);
	//		printf("%\n!!!!!!%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n\n",f,g,dfdp,dfdt,dgdp,dgdt);
			if(fabs(dfdt*dgdp-dfdp*dgdt)<1e-8)
			{
				continue;
			}
			dt=(g*dfdp-f*dgdp)/(dfdt*dgdp-dgdt*dfdp);
			dp=(g*dfdt-f*dgdt)/(dfdp*dgdt-dgdp*dfdt);
			pi_nxt=pi_ini[p_index]+dp;
			tau_nxt=tau_ini[t_index]+dt;
			flag=0;
//	printf("\n%.12lf\t%.12lf\n",dp,dt);
			do
			{
				pi_prev=pi_nxt;
				tau_prev=tau_nxt;
				calcFG_dFG_hs_r5(&f,&g,&dfdp,&dfdt,&dgdp,&dgdt,spe_enth,spe_entr,pi_prev,tau_prev);
				if(fabs(dfdt*dgdp-dfdp*dgdt)<1e-8)
				{
					flag=-1;
					break;
				}
				dt=(g*dfdp-f*dgdp)/(dfdt*dgdp-dgdt*dfdp);
				dp=(g*dfdt-f*dgdt)/(dfdp*dgdt-dgdp*dfdt);
				pi_nxt=pi_prev+dp;
				tau_nxt=tau_prev+dt;
				if(fabs(dt)>fabs(dp))
				{
					maxd=fabs(dt);
				}
				else
				{
					maxd=fabs(dp);
				}
				i++;
			}while(maxd>1e-8&&i<MAX_ITER_TIMES_R5);
			if(i==MAX_ITER_TIMES_R5||flag==-1)
			{
				continue;
			}	
			else
			{
				if(!isnan(pi_nxt)&&!isnan(tau_nxt))
				{
					pp=pi_nxt*1e6;
					tt=1000/tau_nxt;
					steam_prop_calc_r5(pp,tt,&prop);
					if(region(pp,tt)!=5||fabs((spe_enth-prop.spe_enth)/spe_enth)>1e-4||fabs((spe_entr-prop.spe_entr)/spe_entr)>1e-4)
						continue;
					else
					{
						*pres=pp;
						*temp=tt;
						return 0;
					}
				}
			}		
		}
	}
	return -1;
}

/*void print_prop(steam_prop *p_prop)
{
	printf("%.14lf\n%.14lf\n%.14lf\n%.14lf\n%.14lf\n%.14lf\n%.14lf\n%.14lf\n%.14lf\n%.14lf\n",p_prop->pressure,p_prop->temp,p_prop->spe_vol,p_prop->dens,p_prop->spe_energy,p_prop->spe_entr,p_prop->spe_enth,p_prop->spe_h_v,p_prop->spe_h_p,p_prop->speed_sound);
}

int main()
{
	steam_prop prop;
	double p,t;
	steam_prop_calc_r5(0.5e6,1500,&prop);
	print_prop(&prop);
	steam_prop_calc_r5(30e6,1500,&prop);
	print_prop(&prop);
	steam_prop_calc_r5(30e6,2000,&prop);
	print_prop(&prop);
	pres_tr_r5(&p,1500,1/1.3845509);
	printf("\n...%lf\n",p);
		pres_tr_r5(&p,1500,1/0.0230761299);
	printf("\n...%lf\n",p);
		pres_tr_r5(&p,2000,1/0.0311385219);
	printf("\n...%lf\n",p);
	
		pres_tu_r5(&p,1500,4527.4931e3);
	printf("\n...%lf\n",p);
		pres_tu_r5(&p,1500,4474.95214e3);
	printf("\n...%lf\n",p);
		pres_tu_r5(&p,2000,5637.07038e3);
	printf("\n...%lf\n",p);
	
			pres_th_r5(&p,1500,5219.76855e3);
	printf("\n...%lf\n",p);
		pres_th_r5(&p,1500,5167.23514e3);
	printf("\n...%lf\n",p);
		pres_th_r5(&p,2000,6571.22064e3);
	printf("\n...%lf\n",p);
	
	
				pres_ts_r5(&p,1500,9.65408875e3);
	printf("\n...%lf\n",p);
		pres_ts_r5(&p,1500,7.72970133e3);
	printf("\n...%lf\n",p);
		pres_ts_r5(&p,2000,8.53640523e3);
	printf("\n...%lf\n",p);
	
	
	temp_pr_r5(&t,0.5e6,1/1.3845509);
	printf("\n...%lf\n",t);
	temp_pr_r5(&t,30e6,1/0.0230761299);
	printf("\n...%lf\n",t);
		temp_pr_r5(&t,30e6,1/0.0311385219);
	printf("\n...%lf\n",t);
	
		temp_pu_r5(&t,0.5e6,4527.4931e3);
	printf("\n...%lf\n",t);
	temp_pu_r5(&t,30e6,4474.95214e3);
	printf("\n...%lf\n",t);
		temp_pu_r5(&t,30e6,5637.07038e3);
	printf("\n...%lf\n",t);
	
			temp_ph_r5(&t,0.5e6,5219.76855e3);
	printf("\n...%lf\n",t);
	temp_ph_r5(&t,30e6,5167.23514e3);
	printf("\n...%lf\n",t);
		temp_ph_r5(&t,30e6,6571.22064e3);
	printf("\n...%lf\n",t);
	
				temp_ps_r5(&t,0.5e6,9.65408875e3);
	printf("\n...%lf\n",t);
	temp_ps_r5(&t,30e6,7.72970133e3);
	printf("\n...%lf\n",t);
		temp_ps_r5(&t,30e6,8.53640523e3);
	printf("\n...%lf\n",t);
	
	printf("\t%d\n",pt_hs_r5(&p,&t,5219.76855e3,9.65408875e3));
	printf("\n//////%lf,%lf\n",p,t);
		printf("\t%d\n",pt_hs_r5(&p,&t,5167.23514e3,7.72970133e3));
	printf("\n//////%lf,%lf\n",p,t);
		pt_hs_r5(&p,&t,6571.22064e3,8.53640523e3);
	printf("\n//////%lf,%lf\n",p,t);
	
		steam_prop_calc_r5(59204.842478,1505.193671,&prop);
	print_prop(&prop);
}
¡Á/ 


/*int main()
{
	steam_prop prop;
	double pres=30e6;
	double temp=1500;
	steam_prop_calc_r5(pres,temp,&prop);
	print_prop(&prop);
}

*/