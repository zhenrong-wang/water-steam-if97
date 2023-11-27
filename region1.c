/*
 * This code is distributed under the license: MIT License
 * Originally written by Zhenrong WANG
 * mailto: zhenrongwang@live.com
 */

#include <stdio.h>
#include <math.h>
#include "if97_general.h"
#include "region_calc.h"
#include "region1.h"

const double region1_coeff[34][4]={
	1,0,-2,0.14632971213167,
	2,0,-1,-0.84548187169114,
	3,0,0,-0.37563603672040e1,
	4,0,1,0.33855169168385e1,
	5,0,2,-0.95791963387872,
	6,0,3,0.15772038513228,
	7,0,4,-0.16616417199501e-1,
	8,0,5,0.81214629983568e-3,
	9,1,-9,0.28319080123804e-3,
	10,1,-7,-0.60706301565874e-3,
	11,1,-1,-0.18990068218419e-1,
	12,1,0,-0.32529748770505e-1,
	13,1,1,-0.21841717175414e-1,
	14,1,3,-0.52838357969930e-4,
	15,2,-3,-0.47184321073267e-3,
	16,2,0,-0.30001780793026e-3,
	17,2,1,0.47661393906987e-4,
	18,2,3,-0.44141845330846e-5,
	19,2,17,-0.72694996297594e-15,
	20,3,-4,-0.31679644845054e-4,
	21,3,0,-0.28270797985312e-5,
	22,3,6,-0.85205128120103e-9,
	23,4,-5,-0.22425281908000e-5,
	24,4,-2,-0.65171222895601e-6,
	25,4,10,-0.14341729937924e-12,
	26,5,-8,-0.40516996860117e-6,
	27,8,-11,-0.12734301741641e-8,
	28,8,-6,-0.17424871230634e-9,
	29,21,-29,-0.68762131295531e-18,
	30,23,-31,0.14478307828521e-19,
	31,29,-38,0.26335781662795e-22,
	32,30,-39,-0.11947622640071e-22,
	33,31,-40,0.18228094581404e-23,
	34,32,-41,-0.93537087292458e-25,
};

const double region1_coeff_ph[20][4]={
	1,0,0,-0.23872489924521e3,
	2,0,1,0.40421188637945e3,
	3,0,2,0.11349746881718e3,
	4,0,6,-0.58457616048039e1,
	5,0,22,-0.15285482413140e-3,
	6,0,32,-0.10866707695377e-5,
	7,1,0,-0.13391744872602e2,
	8,1,1,0.43211039183559e2,
	9,1,2,-0.54010067170506e2,
	10,1,3,0.30535892203916e2,
	11,1,4,-0.65964749423638e1,
	12,1,10,0.93965400878363e-2,
	13,1,32,0.11573647505340e-6,
	14,2,10,-0.25858641282073e-4,
	15,2,32,-0.40644363084799e-8,
	16,3,10,0.66456186191635e-7,
	17,3,32,0.80670734103027e-10,
	18,4,32,-0.93477771213947e-12,
	19,5,32,0.58265442020601e-14,
	20,6,32,-0.15020185953503e-16,
};

const double region1_coeff_ps[20][4]={
	1,0,0,0.17478268058307e3,
	2,0,1,0.34806930892873e2,
	3,0,2,0.65292584978455e1,
	4,0,3,0.33039981775489,
	5,0,11,-0.19281382923196e-6,
	6,0,31,-0.24909197244573e-22,
	7,1,0,-0.26107636489332,
	8,1,1,0.22592965981586,
	9,1,2,-0.64256463395226e-1,
	10,1,3,0.78876289270526e-2,
	11,1,12,0.35672110607366e-9,
	12,1,31,0.17332496994895e-23,
	13,2,0,0.56608900654837e-3,
	14,2,1,-0.32635483139717e-3,
	15,2,2,0.44778286690632e-4,
	16,2,9,-0.51322156908507e-9,
	17,2,31,-0.42522657042207e-25,
	18,3,10,0.26400441360689e-12,
	19,3,32,0.78124600459723e-28,
	20,4,32,-0.30732199903668e-30,
};

const double region1_coeff_hs[19][4]={
	1,0,0,-0.691997014660582,
	2,0,1,-0.183612548787560e2,
	3,0,2,-0.928332409297335e1,
	4,0,4,0.659639569909906e2,
	5,0,5,-0.162060388912024e2,
	6,0,6,0.450620017338667e3,
	7,0,8,0.854680678224170e3,
	8,0,14,0.607523214001162e4,
	9,1,0,0.326487682621856e2,
	10,1,1,-0.269408844582931e2,
	11,1,4,-0.319947848334300e3,
	12,1,6,-0.928354307043320e3,
	13,2,0,0.303634537455249e2,
	14,2,1,-0.650540422444146e2,
	15,2,10,-0.430991316516130e4,
	16,3,4,-0.747512324096068e3,
	17,4,1,0.730000345529245e3,
	18,4,4,0.114284032569021e4,
	19,5,0,-0.436407041874559e3,
};

double gm(double pi, double tau)
{
	int i;
	double res=0;
	double ii,ji,ni;
	for(i=0;i<34;i++)
	{
		ii=region1_coeff[i][1];
		ji=region1_coeff[i][2];
		ni=region1_coeff[i][3];
		res+=ni*pow(7.1-pi,ii)*pow(tau-1.222,ji);
	}
	return res;
}


double gmpi(double pi, double tau)
{
	int i;
	double res=0;
	double ii,ji,ni;
	for(i=0;i<34;i++)
	{
		ii=region1_coeff[i][1];
		ji=region1_coeff[i][2];
		ni=region1_coeff[i][3];
		res+=(-ni)*ii*pow(7.1-pi,ii-1)*pow(tau-1.222,ji);		
	}
	return res;
}

double gmpipi(double pi, double tau)
{
	int i;
	double res=0;
	double ii,ji,ni;
	for(i=0;i<34;i++)
	{
		ii=region1_coeff[i][1];
		ji=region1_coeff[i][2];
		ni=region1_coeff[i][3];
		res+=ni*ii*(ii-1)*pow(7.1-pi,ii-2)*pow(tau-1.222,ji);
	}
	return res;
}

double gmtau(double pi, double tau)
{
	int i;
	double res=0;
	double ii,ji,ni;
	for(i=0;i<34;i++)
	{
		ii=region1_coeff[i][1];
		ji=region1_coeff[i][2];
		ni=region1_coeff[i][3];
		res+=ni*pow(7.1-pi,ii)*ji*pow(tau-1.222,ji-1);		
	}
	return res;
}

double gmtautau(double pi, double tau)
{
	int i;
	double res=0;
	double ii,ji,ni;
	for(i=0;i<34;i++)
	{
		ii=region1_coeff[i][1];
		ji=region1_coeff[i][2];
		ni=region1_coeff[i][3];
		res+=ni*pow(7.1-pi,ii)*ji*(ji-1)*pow(tau-1.222,ji-2);		
	}
	return res;
}

double gmpitau(double pi, double tau)
{
	int i;
	double res=0;
	double ii,ji,ni;
	for(i=0;i<34;i++)
	{
		ii=region1_coeff[i][1];
		ji=region1_coeff[i][2];
		ni=region1_coeff[i][3];
		res+=(-ni)*ii*pow(7.1-pi,ii-1)*ji*pow(tau-1.222,ji-1);		
	}
	return res;
}

void steam_prop_calc_r1(double pressure, double temp, steam_prop* p_prop)
{
	double pi,tau;
	
	pi=pressure/16.53e6;
	tau=1386/temp;
	
	p_prop->pressure=pressure;
	p_prop->temp=temp;
	p_prop->spe_vol=pi*gmpi(pi,tau)*GAS_CONST_STEAM*temp/pressure;
	p_prop->dens=pressure/(GAS_CONST_STEAM*temp*pi*gmpi(pi,tau));
	p_prop->spe_energy=GAS_CONST_STEAM*temp*(tau*gmtau(pi,tau)-pi*gmpi(pi,tau));
	p_prop->spe_entr=GAS_CONST_STEAM*(tau*gmtau(pi,tau)-gm(pi,tau));
	p_prop->spe_enth=GAS_CONST_STEAM*temp*tau*gmtau(pi,tau);
	p_prop->spe_h_v=GAS_CONST_STEAM*(-tau*tau*gmtautau(pi,tau)+(gmpi(pi,tau)-tau*gmpitau(pi,tau))*(gmpi(pi,tau)-tau*gmpitau(pi,tau))/gmpipi(pi,tau));
	p_prop->spe_h_p=GAS_CONST_STEAM*(-tau*tau*gmtautau(pi,tau));
	p_prop->speed_sound=sqrt(GAS_CONST_STEAM*temp*gmpi(pi,tau)*gmpi(pi,tau)/(pow(gmpi(pi,tau)-tau*gmpitau(pi,tau),2)/(tau*tau*gmtautau(pi,tau))-gmpipi(pi,tau)));
	p_prop->drdp=-p_prop->dens*p_prop->dens*GAS_CONST_STEAM*temp*gmpipi(pi,tau)/(16.53*16.53*1e12);
	p_prop->vapor_fraction=0;
}

double tempini_ph_r1(double pressure, double spec_enth)
{
	double pi, eta, ni, ii, ji;
	int i;
	double theta=0;
	
	pi=pressure/1e6;
	eta=spec_enth/2.5e6;
	for(i=0;i<20;i++)
	{
		ni=region1_coeff_ph[i][3];
		ii=region1_coeff_ph[i][1];
		ji=region1_coeff_ph[i][2];
		theta+=ni*pow(pi,ii)*pow(eta+1,ji);
	}
	return theta;
}

void calcFT_dFT_ph(double* ft, double* dft, double pres, double spec_enth, double tau)
{
	double pi;
	pi=pres/16.53e6;
	*ft=tau*gmtau(pi,tau)-spec_enth*tau/(GAS_CONST_STEAM*1386);
	*dft=tau*gmtautau(pi,tau)+gmtau(pi,tau)-spec_enth/(GAS_CONST_STEAM*1386);
}
int temp_ph_r1(double* temp, double pres, double spe_enth)
{
	int i=0;
	double temp_ini,tau_ini,tau_prev,tau_nxt,ft,dft;
	temp_ini=tempini_ph_r1(pres,spe_enth);
	tau_ini=1386/temp_ini;
	calcFT_dFT_ph(&ft,&dft,pres,spe_enth,tau_ini);
	tau_nxt=tau_ini-ft/dft;
	do
	{
		tau_prev=tau_nxt;
		calcFT_dFT_ph(&ft,&dft,pres,spe_enth,tau_prev);
		tau_nxt=tau_prev-ft/dft;
		i++;
	}while(fabs(tau_prev-tau_nxt)>1e-8&&i<MAX_ITER_NUMS);
	if(i==MAX_ITER_NUMS)
	{
		*temp=temp_ini;
		return 1;
	}
	else
	{
		*temp=1386/tau_nxt;
		return 0;
	}
}

double tempini_ps_r1(double pressure, double spec_entr)
{
	double pi,sigma,ni,ii,ji;
	int i;
	double theta=0;
	
	pi=pressure/1e6;
	sigma=spec_entr/1e3;
	for(i=0;i<20;i++)
	{
		ni=region1_coeff_ps[i][3];
		ii=region1_coeff_ps[i][1];
		ji=region1_coeff_ps[i][2];
		theta+=ni*pow(pi,ii)*pow(sigma+2,ji);
	}
	return theta;
}

void calcFT_dFT_ps(double* ft, double* dft, double pres, double spec_entr, double tau)
{
	double pi;
	pi=pres/16.53e6;
	*ft=tau*gmtau(pi,tau)-gm(pi,tau)-spec_entr/GAS_CONST_STEAM;
	*dft=tau*gmtautau(pi,tau)+gmtau(pi,tau)-gmtau(pi,tau);
}

int temp_ps_r1(double* temp, double pres, double spe_entr)
{
	int i=0;
	double temp_ini,tau_ini,tau_prev,tau_nxt,ft,dft;
	temp_ini=tempini_ps_r1(pres,spe_entr);
	tau_ini=1386/temp_ini;
	calcFT_dFT_ps(&ft,&dft,pres,spe_entr,tau_ini);
	tau_nxt=tau_ini-ft/dft;
	do
	{
		tau_prev=tau_nxt;
		calcFT_dFT_ps(&ft,&dft,pres,spe_entr,tau_prev);
		tau_nxt=tau_prev-ft/dft;
		i++;
	}while(fabs(tau_prev-tau_nxt)>1e-8&&i<MAX_ITER_NUMS);
	if(i==MAX_ITER_NUMS)
	{
		*temp=temp_ini;
		return 1;
	}
	else
	{
		*temp=1386/tau_nxt;
		return 0;
	}
}

double presini_hs_r1(double spec_enth, double spec_entr)
{
	double eta,sigma,ni,ii,ji;
	int i;
	double pi=0;
	
	eta=spec_enth/3400e3;
	sigma=spec_entr/7.6e3;
	for(i=0;i<19;i++)
	{
		ni=region1_coeff_hs[i][3];
		ii=region1_coeff_hs[i][1];
		ji=region1_coeff_hs[i][2];
		pi+=ni*pow(eta+0.05,ii)*pow(sigma+0.05,ji);
	}
	return pi*1e8;
}
double tempini_hs_r1(double spec_enth, double spec_entr)
{
	return tempini_ph_r1(presini_hs_r1(spec_enth,spec_entr),spec_enth);
}

void calcFG_dFG_hs(double* f, double* g, double*dfdp, double* dfdt, double* dgdp, double* dgdt, double spec_enth, double spec_entr, double pi, double tau)
{
	*f=tau*gmtau(pi,tau)-spec_enth*tau/(GAS_CONST_STEAM*1386);
	*g=tau*gmtau(pi,tau)-gm(pi,tau)-spec_entr/GAS_CONST_STEAM;
	*dfdp=tau*gmpitau(pi,tau);
	*dfdt=gmtau(pi,tau)+tau*gmtautau(pi,tau)-spec_enth/(GAS_CONST_STEAM*1386);
	*dgdp=tau*gmpitau(pi,tau)-gmpi(pi,tau);
	*dgdt=tau*gmtautau(pi,tau);
}

int pt_hs_r1(double* pres, double* temp, double spe_enth, double spe_entr)
{
	double pres_ini,pi_ini,pi_prev,pi_nxt;
	double temp_ini,tau_ini,tau_prev,tau_nxt;
	int i=0;
	double f,g,dfdt,dgdt,dfdp,dgdp;
	double dp,dt,maxd;
	
	temp_ini=tempini_hs_r1(spe_enth,spe_entr);
	pres_ini=presini_hs_r1(spe_enth,spe_entr);
	
//	printf("###%.12lf\n###%.12lf\n",temp_ini,pres_ini);
	pi_ini=pres_ini/16.53e6;
	tau_ini=1386/temp_ini;
	
	calcFG_dFG_hs(&f,&g,&dfdp,&dfdt,&dgdp,&dgdt,spe_enth,spe_entr,pi_ini,tau_ini);
//		printf("%\n!!!!!!%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n\n",f,g,dfdp,dfdt,dgdp,dgdt);
	if(fabs(dfdt*dgdp-dfdp*dgdt)<1e-8)
	{
		*pres=pres_ini;
		*temp=temp_ini;
		return -1;
	}
	dt=(g*dfdp-f*dgdp)/(dfdt*dgdp-dgdt*dfdp);
	dp=(g*dfdt-f*dgdt)/(dfdp*dgdt-dgdp*dfdt);
	pi_nxt=pi_ini+dp;
	tau_nxt=tau_ini+dt;
//	printf("\n%.12lf\t%.12lf\n",dp,dt);
	do
	{
		pi_prev=pi_nxt;
		tau_prev=tau_nxt;
		calcFG_dFG_hs(&f,&g,&dfdp,&dfdt,&dgdp,&dgdt,spe_enth,spe_entr,pi_prev,tau_prev);
		if(fabs(dfdt*dgdp-dfdp*dgdt)<1e-8)
		{
			*pres=pres_ini;
			*temp=temp_ini;
			return -1;
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
	}while(maxd>1e-8&&i<MAX_ITER_NUMS);
	if(i==MAX_ITER_NUMS)
	{
		*pres=pres_ini;
		*temp=temp_ini;
		return -2;
	}	
	else
	{
		*pres=pi_nxt*16.53e6;
		*temp=1386/tau_nxt;
		return 0;
	}	
}

void calcFT_dFT_pr(double* ft, double* dft, double pressure, double dens, double T)
{
	double ni,ii,ji,pi,res,res1,dres,dres1;
	int i;
	pi=pressure/16.53e6;
	res=0;
	dres=0;
	for(i=0;i<34;i++)
	{
		ni=region1_coeff[i][3];
		ii=region1_coeff[i][1];
		ji=region1_coeff[i][2];	
		res1=(-ni)*ii*pow(7.1-pi,ii-1)*T*pow(1386/T-1.222,ji);
		res+=res1;
		dres1=(-ni)*ii*pow(7.1-pi,ii-1)*(pow(1386/T-1.222,ji)-ji*pow(1386/T-1.222,ji-1)*1386/T);
		dres+=dres1;
	}
	res1=dens*GAS_CONST_STEAM;
	res*=res1;
	res-=16.53e6;
	dres*=res1;
	*ft=res;
	*dft=dres;
}

int temp_pr_r1(double* temp, double pressure, double dens)
{
	double temp_prev,temp_nxt,ft,dft;
	int ti=0;
	int iter_times;
	double temp_ini[7]={
		450,400,500,350,550,300,600,
	};
	do{
		iter_times=0;
		calcFT_dFT_pr(&ft,&dft,pressure,dens,temp_ini[ti]);
		temp_nxt=temp_ini[ti]-ft/dft;
		do{
			temp_prev=temp_nxt;
			calcFT_dFT_pr(&ft,&dft,pressure,dens,temp_prev);
			temp_nxt=temp_prev-ft/dft;
			iter_times++;
		}while(fabs(temp_nxt-temp_prev)>1e-6&&iter_times<MAX_ITER_NUMS);
		ti++;
		if(iter_times==MAX_ITER_NUMS)
		{
			continue;
		}
		
		else
		{
			*temp=temp_nxt;
			if(NRT(*temp)==0)
			break;
		}
	}while(ti<7);
	if(ti==7)
	{
		*temp=-1;
		return -1;
	}
	else
	{
//		printf("\n\n%d\n\n",ti);
		return 0;
	}
} 

void calcFT_dFT_pu(double* ft, double* dft, double pressure, double spe_ener, double T)
{
	double ni,ii,ji,pi,res,res1,dres,dres1,dres2;
	int i;
	pi=pressure/16.53e6;
	res=0;
	dres=0;
	for(i=0;i<34;i++)
	{
		ni=region1_coeff[i][3];
		ii=region1_coeff[i][1];
		ji=region1_coeff[i][2];
		res1=GAS_CONST_STEAM*1386*ni*pow(7.1-pi,ii)*ji*pow(1386/T-1.222,ji-1)+GAS_CONST_STEAM*pi*ni*ii*pow(7.1-pi,ii-1)*pow(1386/T-1.222,ji)*T;
		res+=res1;
		dres1=-GAS_CONST_STEAM*1386*ni*pow(7.1-pi,ii)*ji*(ji-1)*pow(1386/T-1.222,ji-2)*1386/(T*T);
		dres2=GAS_CONST_STEAM*pi*ni*ii*pow(7.1-pi,ii)*(pow(1386/T-1.222,ji)-ji*pow(1386/T-1.222,ji-1)*1386/T); 
		dres+=(dres1+dres2);  
	}
	*ft=res-spe_ener;
	*dft=dres; 
}

int temp_pu_r1(double* temp, double pressure, double spe_ener)
{
	double temp_prev,temp_nxt,ft,dft;
	int iter_times;
	int ti=0;
	double temp_ini[7]={
		450,400,500,350,550,300,600,
	};
	do{
		iter_times=0;
		calcFT_dFT_pu(&ft,&dft,pressure,spe_ener,temp_ini[ti]);
		temp_nxt=temp_ini[ti]-ft/dft;
		do{
			temp_prev=temp_nxt;
			calcFT_dFT_pu(&ft,&dft,pressure,spe_ener,temp_prev);
			temp_nxt=temp_prev-ft/dft;
			iter_times++;
		}while(fabs(temp_nxt-temp_prev)>1e-6&&iter_times<MAX_ITER_NUMS);
		ti++;
		if(iter_times==MAX_ITER_NUMS)
		{
			continue;
		}
		else
		{
			*temp=temp_nxt;
			if(NRT(*temp)==0)
			break;
		}	
	}while(ti<7);
	if(ti==7)
	{
		*temp=-1;
		return -1;
	}
	else
	{
		return 0;
	}
}

void calcFPi_dFPi_tr(double* fpi, double* dfpi, double Temp, double dens, double pi)
{
	double ni,ii,ji,tau,res,res1,dres,dres1;
	int i;
	
	tau=1386/Temp;
	res=0;
	dres=0;
	for(i=0;i<34;i++)
	{
		ni=region1_coeff[i][3];
		ii=region1_coeff[i][1];
		ji=region1_coeff[i][2];
		res1=(-ni)*ii*pow(7.1-pi,ii-1)*pow(tau-1.222,ji);
		res+=res1;
		dres1=ni*ii*(ii-1)*pow(tau-1.222,ji)*pow(7.1-pi,ii-2);
		dres+=dres1;
	}
	res-=(16.53e6/(dens*Temp*GAS_CONST_STEAM));
	*fpi=res;
	*dfpi=dres;	
}

int pres_tr_r1(double* pres, double temp, double dens)
{
	double pi_prev,pi_nxt,fpi,dfpi;
	int iter_times;
	int ti=0;
	double pi_ini[11]={
		1,2,1e-1,1e-2,3,4,1e-7,1e-6,1e-5,1e-4,1e-3,	
	};
	do{
		calcFPi_dFPi_tr(&fpi,&dfpi,temp,dens,pi_ini[ti]);
		pi_nxt=pi_ini[ti]-fpi/dfpi;
		do{
			pi_prev=pi_nxt;
			calcFPi_dFPi_tr(&fpi,&dfpi,temp,dens,pi_prev);
			pi_nxt=pi_prev-fpi/dfpi;
			iter_times++;
		}while(fabs(pi_nxt-pi_prev)>1e-8&&iter_times<MAX_ITER_NUMS);
		ti++;
		if(iter_times==MAX_ITER_NUMS)
		{
			continue;
		}
		else
		{
			*pres=pi_nxt*16.53e6;
			if(NRP(*pres)==0)
			break;
		}		
	}while(ti<11);
	if(ti==11)
	{
		*pres=-1;
		return -1;
	}	
	else
	{
		return 0;
	}
}

void calcFPi_dFPi_tu(double* fpi, double* dfpi, double temp, double spe_ener, double pi)
{
	double ni,ii,ji,tau,res,dres,res1,dres1,res2,dres2;
	int i;
	
	tau=1386/temp;
		
	res=0;
	dres=0;
	
	for(i=0;i<34;i++)
	{
		ni=region1_coeff[i][3];
		ii=region1_coeff[i][1];
		ji=region1_coeff[i][2];
		res1=tau*ni*pow(7.1-pi,ii)*ji*pow(tau-1.222,ji-1);
		res2=pi*ni*ii*pow(7.1-pi,ii-1)*pow(tau-1.222,ji);
		res+=(res1+res2);
		dres1=-tau*ni*ji*pow(tau-1.222,ji-1)*ii*pow(7.1-pi,ii-1);
		dres2=ni*ii*pow(tau-1.222,ji)*(pow(7.1-pi,ii-1)-pi*(ii-1)*pow(7.1-pi,ii-2));
		dres+=(dres1+dres2);	
	}
	res-=(spe_ener/(GAS_CONST_STEAM*temp));
	*fpi=res;
	*dfpi=dres;	
} 

int pres_tu_r1(double* pres, double temp, double spe_ener)
{
	double pi_prev,pi_nxt,fpi,dfpi;
	int iter_times;
	int ti=0;
	double pi_ini[11]={
		1,2,1e-1,1e-2,3,4,1e-7,1e-6,1e-5,1e-4,1e-3,	
	};
	
	do{
		iter_times=0;
		calcFPi_dFPi_tu(&fpi,&dfpi,temp,spe_ener,pi_ini[ti]);
		pi_nxt=pi_ini[ti]-fpi/dfpi;
		do{
			pi_prev=pi_nxt;
			calcFPi_dFPi_tu(&fpi,&dfpi,temp,spe_ener,pi_prev);
			pi_nxt=pi_prev-fpi/dfpi;
			iter_times++;
//		printf("#%lf\n",pi_nxt);
		}while(fabs(pi_nxt-pi_prev)>1e-8&&iter_times<MAX_ITER_NUMS);
		ti++;
		if(iter_times==MAX_ITER_NUMS)
		{
			continue;
		}
		else
		{
			*pres=pi_nxt*16.53e6;
			if(NRP(*pres)==0)
			break;
		}
	}while(ti<11);
	if(ti==11)
	{
		*pres=-1;
		return -1;
	} 
	else
	{
		return 0;
	}
}

void calcFPi_dFPi_th(double* fpi, double* dfpi, double temp, double spe_enth, double pi)
{
	double ni,ii,ji,tau,res,dres,res1,dres1;
	int i;
	
	tau=1386/temp;
		
	res=0;
	dres=0;
	for(i=0;i<34;i++)
	{
		ni=region1_coeff[i][3];
		ii=region1_coeff[i][1];
		ji=region1_coeff[i][2];
		res1=ni*pow(7.1-pi,ii)*ji*pow(tau-1.222,ji-1);
		res+=res1;
		dres1=-ni*ji*pow(tau-1.222,ji)*ii*pow(7.1-pi,ii-1);
		dres+=dres1;
	}
	res=res*tau-spe_enth/(GAS_CONST_STEAM*temp);
	dres*=tau;
	*fpi=res;
	*dfpi=dres;
}

int pres_th_r1(double* pres, double temp, double spe_enth)
{
	double pi_prev,pi_nxt,fpi,dfpi;
	int iter_times;
	int ti=0;
	double pi_ini[11]={
		1,2,1e-1,1e-2,3,4,1e-7,1e-6,1e-5,1e-4,1e-3,	
	};
	do{
		iter_times=0;
		calcFPi_dFPi_th(&fpi,&dfpi,temp,spe_enth,pi_ini[ti]);
		pi_nxt=pi_ini[ti]-fpi/dfpi;
		do{
			pi_prev=pi_nxt;
			calcFPi_dFPi_th(&fpi,&dfpi,temp,spe_enth,pi_prev);
			pi_nxt=pi_prev-fpi/dfpi;
			iter_times++;
//		printf("#%lf\n",pi_nxt);
		}while(fabs(pi_nxt-pi_prev)>1e-8&&iter_times<MAX_ITER_NUMS);
		ti++;
		if(iter_times==MAX_ITER_NUMS)
		{
			continue;
		}
		else
		{
			*pres=pi_nxt*16.53e6;
			if(NRP(*pres)==0)
			break;
		}
	}while(ti<11);
	if(ti==11)
	{
		*pres=-1;
		return -1;
	}	
	else
	{
		return 0;
	}
}

void calcFPi_dFPi_ts(double* fpi, double* dfpi, double temp, double spe_entr, double pi)
{
	double ni,ii,ji,tau,res,dres,res1,dres1,res2,dres2;
	int i;
	
	tau=1386/temp;
		
	res=0;
	dres=0;
	for(i=0;i<34;i++)
	{
		ni=region1_coeff[i][3];
		ii=region1_coeff[i][1];
		ji=region1_coeff[i][2];
		res1=tau*ni*pow(7.1-pi,ii)*ji*pow(tau-1.222,ji-1);
		res2=ni*pow(7.1-pi,ii)*pow(tau-1.222,ji);
		res+=(res1-res2);
		dres1=-tau*ni*ji*pow(tau-1.222,ji-1)*ii*pow(7.1-pi,ii-1);
		dres2=ni*pow(tau-1.222,ji)*ii*pow(7.1-pi,ii-1);
		dres+=(dres1+dres2);
	}
	res-=(spe_entr/GAS_CONST_STEAM);
	*fpi=res;
	*dfpi=dres;
} 

int pres_ts_r1(double* pres, double temp, double spe_entr)
{
	double pi_prev,pi_nxt,fpi,dfpi;
	int iter_times;
	int ti=0;
	double pi_ini[11]={
		1,2,1e-1,1e-2,3,4,1e-7,1e-6,1e-5,1e-4,1e-3,	
	};
	
	do{
		iter_times=0;
		calcFPi_dFPi_ts(&fpi,&dfpi,temp,spe_entr,pi_ini[ti]);
		pi_nxt=pi_ini[ti]-fpi/dfpi;
		do{
			pi_prev=pi_nxt;
			calcFPi_dFPi_ts(&fpi,&dfpi,temp,spe_entr,pi_prev);
			pi_nxt=pi_prev-fpi/dfpi;
			iter_times++;
//		printf("#%lf\n",pi_nxt);
		}while(fabs(pi_nxt-pi_prev)>1e-8&&iter_times<MAX_ITER_NUMS);
		ti++;
		if(iter_times==MAX_ITER_NUMS)
		{
			continue;
		}
		else
		{
			*pres=pi_nxt*16.53e6;
			if(NRP(*pres)==0)
			break;
		}	
	}while(ti<11);
	if(ti==11)
	{
		*pres=-1;
		return -1;
	}
	else
	{
		return 0;
	}
}