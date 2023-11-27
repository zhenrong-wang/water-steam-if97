/*
 * This code is distributed under the license: MIT License
 * Originally written by Zhenrong WANG
 * mailto: zhenrongwang@live.com
 */

#include <stdio.h>
#include <math.h>
#include "if97_general.h"
#include "region_calc.h"
#include "region2.h"

static double region2_coeff0[9][3]={
	1,0,-0.96927686500217e1,
	2,1,0.10086655968018e2,
	3,-5,-0.56087911283020e-2,
	4,-4,0.71452738081455e-1,
	5,-3,-0.40710498223928,
	6,-2,0.14240819171444e1,
	7,-1,-0.43839511319450e1,
	8,2,-0.28408632460772,
	9,3,0.21268463753307e-1,
};

static double region2_coeffr[43][4]={
	1,1,0,-0.17731742473213e-2,
	2,1,1,-0.17834862292358e-1,
	3,1,2,-0.45996013696365e-1,
	4,1,3,-0.57581259083432e-1,
	5,1,6,-0.50325278727930e-1,
	6,2,1,-0.33032641670203e-4,
	7,2,2,-0.18948987516315e-3,
	8,2,4,-0.39392777243355e-2,
	9,2,7,-0.43797295650573e-1,
	10,2,36,-0.26674547914087e-4,
	11,3,0,0.20481737692309e-7,
	12,3,1,0.43870667284435e-6,
	13,3,3,-0.32277677238570e-4,
	14,3,6,-0.15033924542148e-2,
	15,3,35,-0.40668253562649e-1,
	16,4,1,-0.78847309559367e-9,
	17,4,2,0.12790717852285e-7,
	18,4,3,0.48225372718507e-6,
	19,5,7,0.22922076337661e-5,
	20,6,3,-0.16714766451061e-10,
	21,6,16,-0.21171472321355e-2,
	22,6,35,-0.23895741934104e2,
	23,7,0,-0.59059564324270e-17,
	24,7,11,-0.12621808899101e-5,
	25,7,25,-0.38946842435739e-1,
	26,8,8,0.11256211360459e-10,
	27,8,36,-0.82311340897998e1,
	28,9,13,0.19809712802088e-7,
	29,10,4,0.10406965210174e-18,
	30,10,10,-0.10234747095929e-12,
	31,10,14,-0.10018179379511e-8,
	32,16,29,-0.80882908646985e-10,
	33,16,50,0.10693031879409,
	34,18,57,-0.33662250574171,
	35,20,20,0.89185845355421e-24,
	36,20,35,0.30629316876232e-12,
	37,20,48,-0.42002467698208e-5,
	38,21,21,-0.59056029685639e-25,
	39,22,53,0.37826947613457e-5,
	40,23,39,-0.12768608934681e-14,
	41,24,26,0.73087610595061e-28,
	42,24,40,0.55414715350778e-16,
	43,24,58,-0.94369707241210e-6,
};

static double region2a_coeff_ph[34][4]={
	1,0,0,0.10898952318288e4,
	2,0,1,0.84951654495535e3,
	3,0,2,-0.10781748091826e3,
	4,0,3,0.33153654801263e2,
	5,0,7,-0.74232016790248e1,
	6,0,20,0.11765048724356e2,
	7,1,0,0.18445749355790e1,
	8,1,1,-0.41792700549624e1,
	9,1,2,0.62478196935812e1,
	10,1,3,-0.17344563108114e2,
	11,1,7,-0.20058176862096e3,
	12,1,9,0.27196065473796e3,
	13,1,11,-0.45511318285818e3,
	14,1,18,0.30919688604755e4,
	15,1,44,0.25226640357872e6,
	16,2,0,-0.61707422868339e-2,
	17,2,2,-0.31078046629583,
	18,2,7,0.11670873077107e2,
	19,2,36,0.12812798404046e9,
	20,2,38,-0.98554909623276e9,
	21,2,40,0.28224546973002e10,
	22,2,42,-0.35948971410703e10,
	23,2,44,0.17227349913197e10,
	24,3,24,-0.13551334240775e5,
	25,3,44,0.12848734664650e8,
	26,4,12,0.13865724283226e1,
	27,4,32,0.23598832556514e6,
	28,4,44,-0.13105236545054e8,
	29,5,32,0.73999835474766e4,
	30,5,36,-0.55196697030060e6,
	31,5,42,0.37154085996233e7,
	32,6,34,0.19127729239660e5,
	33,6,44,-0.41535164835634e6,
	34,7,28,-0.62459855192507e2,
};

static double region2b_coeff_ph[38][4]={
	1,0,0,0.14895041079516e4,
	2,0,1,0.74307798314034e3,
	3,0,2,-0.97708318797837e2,
	4,0,12,0.24742464705674e1,
	5,0,18,-0.63281320016026,
	6,0,24,0.11385952129658e1,
	7,0,28,-0.47811863648625,
	8,0,40,0.85208123431544e-2,
	9,1,0,0.93747147377932,
	10,1,2,0.33593118604916e1,
	11,1,6,0.33809355601454e1,
	12,1,12,0.16844539671904,
	13,1,18,0.73875745236695,
	14,1,24,-0.47128737436186,
	15,1,28,0.15020273139707,
	16,1,40,-0.21764114219750e-2,
	17,2,2,-0.21810755324761e-1,
	18,2,8,-0.10829784403677,
	19,2,18,-0.46333324635812e-1,
	20,2,40,0.71280351959551e-4,
	21,3,1,0.11032831789999e-3,
	22,3,2,0.18955248387902e-3,
	23,3,12,0.30891541160537e-2,
	24,3,24,0.13555504554949e-2,
	25,4,2,0.28640237477456e-6,
	26,4,12,-0.10779857357512e-4,
	27,4,18,-0.76462712454814e-4,
	28,4,24,0.14052392818316e-4,
	29,4,28,-0.31083814331434e-4,
	30,4,40,-0.10302738212103e-5,
	31,5,18,0.28217281635040e-6,
	32,5,24,0.12704902271945e-5,
	33,5,40,0.73803353468292e-7,
	34,6,28,-0.11030139238909e-7,
	35,7,2,-0.81456365207833e-13,
	36,7,28,-0.25180545682962e-10,
	37,9,1,-0.17565233969407e-17,
	38,9,40,0.86934156344163e-14,
};

static double region2c_coeff_ph[23][4]={
	1,-7,0,-0.32368398555242e13,
	2,-7,4,0.73263350902181e13,
	3,-6,0,0.35825089945447e12,
	4,-6,2,-0.58340131851590e12,
	5,-5,0,-0.10783068217470e11,
	6,-5,2,0.20825544563171e11,
	7,-2,0,0.61074783564516e6,
	8,-2,1,0.85977722535580e6,
	9,-1,0,-0.25745723604170e5,
	10,-1,2,0.31081088422714e5,
	11,0,0,0.12082315865936e4,
	12,0,1,0.48219755109255e3,
	13,1,4,0.37966001272486e1,
	14,1,8,-0.10842984880077e2,
	15,2,4,-0.45364172676660e-1,
	16,6,0,0.14559115658698e-12,
	17,6,1,0.11261597407230e-11,
	18,6,4,-0.17804982240686e-10,
	19,6,10,0.12324579690832e-6,
	20,6,12,-0.11606921130984e-5,
	21,6,16,0.27846367088554e-4,
	22,6,20,-0.59270038474176e-3,
	23,6,22,0.12918582991878e-2,
};

static double region2a_coeff_ps[46][4]={
	1,-1.5,-24,-0.39235983861984e6,
	2,-1.5,-23,0.51526573827270e6,
	3,-1.5,-19,0.40482443161048e5,
	4,-1.5,-13,-0.32193790923902e3,
	5,-1.5,-11,0.96961424218694e2,
	6,-1.5,-10,-0.22867846371773e2,
	7,-1.25,-19,-0.44942914124357e6,
	8,-1.25,-15,-0.50118336020166e4,
	9,-1.25,-6,0.35684463560015,
	10,-1.0,-26,0.44235335848190e5,
	11,-1.0,-21,-0.13673388811708e5,
	12,-1.0,-17,0.42163260207864e6,
	13,-1.0,-16,0.22516925837475e5,
	14,-1.0,-9,0.47442144865646e3,
	15,-1.0,-8,-0.14931130797647e3,
	16,-0.75,-15,-0.19781126320452e6,
	17,-0.75,-14,-0.23554399470760e5,
	18,-0.5,-26,-0.19070616302076e5,
	19,-0.5,-13,0.55375669883164e5,
	20,-0.5,-9,0.38293691437363e4,
	21,-0.5,-7,-0.60391860580567e3,
	22,-0.25,-27,0.19363102620331e4,
	23,-0.25,-25,0.42660643698610e4,
	24,-0.25,-11,-0.59780638872718e4,
	25,-0.25,-6,-0.70401463926862e3,
	26,0.25,1,0.33836784107553e3,
	27,0.25,4,0.20862786635187e2,
	28,0.25,8,0.33834172656196e-1,
	29,0.25,11,-0.43124428414893e-4,
	30,0.5,0,0.16653791356412e3,
	31,0.5,1,-0.13986292055898e3,
	32,0.5,5,-0.78849547999872,
	33,0.5,6,0.72132411753872e-1,
	34,0.5,10,-0.59754839398283e-2,
	35,0.5,14,-0.12141358953904e-4,
	36,0.5,16,0.23227096733871e-6,
	37,0.75,0,-0.10538463566194e2,
	38,0.75,4,0.20718925496502e1,
	39,0.75,9,-0.72193155260427e-1,
	40,0.75,17,0.20749887081120e-6,
	41,1,7,-0.18340657911379e-1,
	42,1,18,0.29036272348696e-6,
	43,1.25,3,0.21037527893619,
	44,1.25,15,0.25681239729999e-3,
	45,1.5,5,-0.12799002933781e-1,
	46,1.5,18,-0.82198102652018e-5,
};

const double region2b_coeff_ps[44][4]={
	1,-6,0,0.31687665083497e6,
2,-6,11,0.20864175881858e2,
3,-5,0,-0.39859399803599e6,
4,-5,11,-0.21816058518877e2,
5,-4,0,0.22369785194242e6,
6,-4,1,-0.27841703445817e4,
7,-4,11,0.99207436071480e1,
8,-3,0,-0.75197512299157e5,
9,-3,1,0.29708605951158e4,
10,-3,11,-0.34406878548526e1,
11,-3,12,0.38815564249115,
12,-2,0,0.17511295085750e5,
13,-2,1,-0.14237112854449e4,
14,-2,6,0.10943803364167e1,
15,-2,10,0.89971619308495,
16,-1,0,-0.33759740098958e4,
17,-1,1,0.47162885818355e3,
18,-1,5,-0.19188241993679e1,
19,-1,8,0.41078580492196,
20,-1,9,-0.33465378172097,
21,0,0,0.13870034777505e4,
22,0,1,-0.40663326195838e3,
23,0,2,0.41727347159610e2,
24,0,4,0.21932549434532e1,
25,0,5,-0.10320050009077e1,
26,0,6,0.35882943516703,
27,0,9,0.52511453726066e-2,
28,1,0,0.12838916450705e2,
29,1,1,-0.28642437219381e1,
30,1,2,0.56912683664855,
31,1,3,-0.99962954584931e-1,
32,1,7,-0.32632037778459e-2,
33,1,8,0.23320922576723e-3,
34,2,0,-0.15334809857450,
35,2,1,0.29072288239902e-1,
36,2,5,0.37534702741167e-3,
37,3,0,0.17296691702411e-2,
38,3,1,-0.38556050844504e-3,
39,3,3,-0.35017712292608e-4,
40,4,0,-0.14566393631492e-4,
41,4,1,0.56420857267269e-5,
42,5,0,0.41286150074605e-7,
43,5,1,-0.20684671118824e-7,
44,5,2,0.16409393674725e-8,
};

const double region2c_coeff_ps[30][4]={
	1,-2,0,0.90968501005365e3,
2,-2,1,0.24045667088420e4,
3,-1,0,-0.59162326387130e3,
4,0,0,0.54145404128074e3,
5,0,1,-0.27098308411192e3,
6,0,2,0.97976525097926e3,
7,0,3,-0.46966772959435e3,
8,1,0,0.14399274604723e2,
9,1,1,-0.19104204230429e2,
10,1,3,0.53299167111971e1,
11,1,4,-0.21252975375934e2,
12,2,0,-0.31147334413760,
13,2,1,0.60334840894623,
14,2,2,-0.42764839702509e-1,
15,3,0,0.58185597255259e-2,
16,3,1,-0.14597008284753e-1,
17,3,5,0.56631175631027e-2,
18,4,0,-0.76155864584577e-4,
19,4,1,0.22440342919332e-3,
20,4,4,-0.12561095013413e-4,
21,5,0,0.63323132660934e-6,
22,5,1,-0.20541989675375e-5,
23,5,2,0.36405370390082e-7,
24,6,0,-0.29759897789215e-8,
25,6,1,0.10136618529763e-7,
26,7,0,0.59925719692351e-11,
27,7,1,-0.20677870105164e-10,
28,7,3,-0.20874278181886e-10,
29,7,4,0.10162166825089e-9,
30,7,5,-0.16429828281347e-9,
};

const double region2a_coeff_hs[29][4]={
	1,0,1,-0.182575361923032e-1,
2,0,3,-0.125229548799536,
3,0,6,0.592290437320145,
4,0,16,0.604769706185122e1,
5,0,20,0.238624965444474e3,
6,0,22,-0.298639090222922e3,
7,1,0,0.512250813040750e-1,
8,1,1,-0.437266515606486,
9,1,2,0.413336902999504,
10,1,3,-0.516468254574773e1,
11,1,5,-0.557014838445711e1,
12,1,6,0.128555037824478e2,
13,1,10,0.114144108953290e2,
14,1,16,-0.119504225652714e3,
15,1,20,-0.284777985961560e4,
16,1,22,0.431757846408006e4,
17,2,3,0.112894040802650e1,
18,2,16,0.197409186206319e4,
19,2,20,0.151612444706087e4,
20,3,0,0.141324451421235e-1,
21,3,2,0.585501282219601,
22,3,3,-0.297258075863012e1,
23,3,6,0.594567314847319e1,
24,3,16,-0.623656565798905e4,
25,4,16,0.965986235133332e4,
26,5,3,0.681500934948134e1,
27,5,16,-0.633207286824489e4,
28,6,3,-0.558919224465760e1,
29,7,1,0.400645798472063e-1,
};

const double region2b_coeff_hs[33][4]={
	1,0,0,0.801496989929495e-1,
2,0,1,-0.543862807146111,
3,0,2,0.337455597421283,
4,0,4,0.890555451157450e1,
5,0,8,0.313840736431485e3,
6,1,0,0.797367065977789,
7,1,1,-0.121616973556240e1,
8,1,2,0.872803386937477e1,
9,1,3,-0.169769781757602e2,
10,1,5,-0.186552827328416e3,
11,1,12,0.951159274344237e5,
12,2,1,-0.189168510120494e2,
13,2,6,-0.433407037194840e4,
14,2,18,0.543212633012715e9,
15,3,0,0.144793408386013,
16,3,1,0.128024559637516e3,
17,3,7,-0.672309534071268e5,
18,3,12,0.336972380095287e8,
19,4,1,-0.586634196762720e3,
20,4,16,-0.221403224769889e11,
21,5,1,0.171606668708389e4,
22,5,12,-0.570817595806302e9,
23,6,1,-0.312109693178482e4,
24,6,8,-0.207841384633010e7,
25,6,18,0.305605946157786e13,
26,7,1,0.322157004314333e4,
27,7,16,0.326810259797295e12,
28,8,1,-0.144104158934487e4,
29,8,3,0.410694867802691e3,
30,8,14,0.109077066873024e12,
31,8,18,-0.247964654258893e14,
32,12,10,0.188801906865134e10,
33,14,16,-0.123651009018773e15,
};

const double region2c_coeff_hs[31][4]={
	1,0,0,0.112225607199012,
2,0,1,-0.339005953606712e1,
3,0,2,-0.320503911730094e2,
4,0,3,-0.197597305104900e3,
5,0,4,-0.407693861553446e3,
6,0,8,0.132943775222331e5,
7,1,0,0.170846839774007e1,
8,1,2,0.373694198142245e2,
9,1,5,0.358144365815434e4,
10,1,8,0.423014446424664e6,
11,1,14,-0.751071025760063e9,
12,2,2,0.523446127607898e2,
13,2,3,-0.228351290812417e3,
14,2,7,-0.960652417056937e6,
15,2,10,-0.807059292526074e8,
16,2,18,0.162698017225669e13,
17,3,0,0.772465073604171,
18,3,5,0.463929973837746e5,
19,3,8,-0.137317885134128e8,
20,3,16,0.170470392630512e13,
21,3,18,-0.251104628187308e14,
22,4,18,0.317748830835520e14,
23,5,1,0.538685623675312e2,
24,5,4,-0.553089094625169e5,
25,5,6,-0.102861522421405e7,
26,5,14,0.204249418756234e13,
27,6,8,0.273918446626977e9,
28,6,18,-0.263963146312685e16,
29,10,7,-0.107890854108088e10,
30,12,7,-0.296492620980124e11,
31,16,10,-0.111754907323424e16,
};

double gm0(double pi, double tau)
{
	double ni,ji;
	double res=log(pi);
	int i;
	for(i=0;i<9;i++)
	{
		ni=region2_coeff0[i][2];
		ji=region2_coeff0[i][1];
		res+=(ni*pow(tau,ji));
	}
	return res;
} 

double gm0pi(double pi)
{
	return 1/pi;
}

double gm0pipi(double pi)
{
	return -1/(pi*pi);
}

double gm0tau(double tau)
{
	double ni,ji;
	double res=0;
	int i;
	for(i=0;i<9;i++)
	{
		ni=region2_coeff0[i][2];
		ji=region2_coeff0[i][1];
		res+=(ni*ji*pow(tau,ji-1));
	}
	return res;
}

double gm0tautau(double tau)
{
	double ni,ji;
	double res=0;
	int i;
	for(i=0;i<9;i++)
	{
		ni=region2_coeff0[i][2];
		ji=region2_coeff0[i][1];
		res+=(ni*ji*(ji-1)*pow(tau,ji-2));
	}
	return res;
}

double gmr(double pi, double tau)
{
	int i;
	double res=0;
	double ni,ii,ji;
	for(i=0;i<43;i++)
	{
		ni=region2_coeffr[i][3];
		ii=region2_coeffr[i][1];
		ji=region2_coeffr[i][2];
		res+=(ni*pow(pi,ii)*pow(tau-0.5,ji));
	}
	return res;
}

double gmrpi(double pi, double tau)
{
	int i;
	double res=0;
	double ni,ii,ji;
	for(i=0;i<43;i++)
	{
		ni=region2_coeffr[i][3];
		ii=region2_coeffr[i][1];
		ji=region2_coeffr[i][2];
		res+=(ni*ii*pow(pi,ii-1)*pow(tau-0.5,ji));
	}
	return res;
}

double gmrpipi(double pi, double tau)
{
	int i;
	double res=0;
	double ni,ii,ji;
	for(i=0;i<43;i++)
	{
		ni=region2_coeffr[i][3];
		ii=region2_coeffr[i][1];
		ji=region2_coeffr[i][2];
		res+=(ni*ii*(ii-1)*pow(pi,ii-2)*pow(tau-0.5,ji));
	}
	return res;
}

double gmrtau(double pi, double tau)
{
	int i;
	double res=0;
	double ni,ii,ji;
	for(i=0;i<43;i++)
	{
		ni=region2_coeffr[i][3];
		ii=region2_coeffr[i][1];
		ji=region2_coeffr[i][2];
		res+=(ni*pow(pi,ii)*ji*pow(tau-0.5,ji-1));
	}
	return res;
}

double gmrtautau(double pi, double tau)
{
	int i;
	double res=0;
	double ni,ii,ji;
	for(i=0;i<43;i++)
	{
		ni=region2_coeffr[i][3];
		ii=region2_coeffr[i][1];
		ji=region2_coeffr[i][2];
		res+=(ni*pow(pi,ii)*ji*(ji-1)*pow(tau-0.5,ji-2));
	}
	return res;
}

double gmrpitau(double pi, double tau)
{
	int i;
	double res=0;
	double ni,ii,ji;
	for(i=0;i<43;i++)
	{
		ni=region2_coeffr[i][3];
		ii=region2_coeffr[i][1];
		ji=region2_coeffr[i][2];
		res+=(ni*ii*pow(pi,ii-1)*ji*pow(tau-0.5,ji-1));
	}
	return res;
}




void steam_prop_calc_r2(double pressure, double temp, steam_prop* p_prop)
{
	double pi,tau;

	pi=pressure/1e6;
	tau=540/temp;
	
	p_prop->pressure=pressure;
	p_prop->temp=temp;
	p_prop->spe_vol=pi*(gm0pi(pi)+gmrpi(pi,tau))*GAS_CONST_STEAM*temp/pressure;
	p_prop->dens=pressure/(GAS_CONST_STEAM*temp*(pi*(gm0pi(pi)+gmrpi(pi,tau))));
	p_prop->spe_energy=GAS_CONST_STEAM*temp*(tau*(gm0tau(tau)+gmrtau(pi,tau))-pi*(gm0pi(pi)+gmrpi(pi,tau)));
	p_prop->spe_entr=GAS_CONST_STEAM*(tau*(gm0tau(tau)+gmrtau(pi,tau))-(gm0(pi,tau)+gmr(pi,tau)));
	p_prop->spe_enth=GAS_CONST_STEAM*temp*tau*(gm0tau(tau)+gmrtau(pi,tau));
	p_prop->spe_h_v=GAS_CONST_STEAM*(pow(1+pi*gmrpi(pi,tau)-tau*pi*gmrpitau(pi,tau),2)/(1-pi*pi*gmrpipi(pi,tau))-tau*tau*(gm0tautau(tau)+gmrtautau(pi,tau)));
	p_prop->spe_h_p=-GAS_CONST_STEAM*tau*tau*(gm0tautau(tau)+gmrtautau(pi,tau));
	p_prop->speed_sound=sqrt(GAS_CONST_STEAM*temp*(1+2*pi*gmrpi(pi,tau)+pi*pi*gmrpi(pi,tau)*gmrpi(pi,tau))/(1-pi*pi*gmrpipi(pi,tau)+pow(1+pi*gmrpi(pi,tau)-tau*pi*gmrpitau(pi,tau),2)/(tau*tau*(gm0tautau(tau)+gmrtautau(pi,tau)))));
	p_prop->drdp=-p_prop->dens*p_prop->dens*GAS_CONST_STEAM*temp*(gm0pipi(pi)+gmrpipi(pi,tau))/(1e12);
	p_prop->vapor_fraction=1;
}

int sub_region(double pressure, double spec_enth)
{
	double eta;
	double pi, pi_ref;
	double n1,n2,n3;
	
	n1=0.90584278514723e3;
	n2=-0.67955786399241;
	n3=0.12809002730136e-3;
	
	if(pressure<4e6)
	{
		return 21;
	}
	else
	{
		eta=spec_enth/1000;
		pi=pressure/1e6;
		pi_ref=n1+n2*eta+n3*eta*eta;
		if(pi>pi_ref||pi==pi_ref)
		{
			return 23;
		}
		else
		{
			return 22;
		} 
	}
}

int temp_ph_r2(double* temp, double pressure, double spec_enth)
{
	int sub_reg;
	int i;
	double ni,ii,ji,eta,pi,tau;
	tau=0;
	
	pi=pressure/1e6;
	eta=spec_enth/2e6;
	
	sub_reg=sub_region(pressure,spec_enth);
	if(sub_reg==21)
	{
		for(i=0;i<34;i++)
		{
			ni=region2a_coeff_ph[i][3];
			ii=region2a_coeff_ph[i][1];
			ji=region2a_coeff_ph[i][2];
			tau+=ni*pow(pi,ii)*pow(eta-2.1,ji);
		}
	}
	else if(sub_reg==22)
	{
		for(i=0;i<38;i++)
		{
			ni=region2b_coeff_ph[i][3];
			ii=region2b_coeff_ph[i][1];
			ji=region2b_coeff_ph[i][2];
			tau+=ni*pow(pi-2,ii)*pow(eta-2.6,ji);
		}
	}
	else if(sub_reg==23)
	{
		for(i=0;i<23;i++)
		{
			ni=region2c_coeff_ph[i][3];
			ii=region2c_coeff_ph[i][1];
			ji=region2c_coeff_ph[i][2];
			tau+=ni*pow(pi+25,ii)*pow(eta-1.8,ji);
		}
	}
	*temp=tau;
	return 0;
}

int sub_region_s(double pressure, double spe_entr)
{
	if(pressure<4e6)
	{
		return 21;
	}
	else
	{
		if(spe_entr<5.85e3)
		{
			return 23;
		}
		else
		{
			return 22;
		}
	}
}

int temp_ps_r2(double* temp, double pressure, double spe_entr)
{
	int sub_reg;
	double pi,sigma2a,sigma2b,sigma2c,ni,ii,ji;
	int i;
	double tau=0;
	
	pi=pressure/1e6;
	sigma2a=spe_entr/2e3;
	sigma2b=spe_entr/0.7853e3;
	sigma2c=spe_entr/2.9251e3;
	
	sub_reg=sub_region_s(pressure,spe_entr);
	if(sub_reg==21)
	{
		for(i=0;i<46;i++)
		{
			ni=region2a_coeff_ps[i][3];
			ii=region2a_coeff_ps[i][1];
			ji=region2a_coeff_ps[i][2];
			tau+=ni*pow(pi,ii)*pow(sigma2a-2,ji);
		}
	}
	
	else if(sub_reg==22)
	{
		for(i=0;i<44;i++)
		{
			ni=region2b_coeff_ps[i][3];
			ii=region2b_coeff_ps[i][1];
			ji=region2b_coeff_ps[i][2];
			tau+=ni*pow(pi,ii)*pow(10-sigma2b,ji);
		}
	}
	
	else if(sub_reg==23)
	{
		for(i=0;i<30;i++)
		{
			ni=region2c_coeff_ps[i][3];
			ii=region2c_coeff_ps[i][1];
			ji=region2c_coeff_ps[i][2];
			tau+=ni*pow(pi,ii)*pow(2-sigma2c,ji);
		}
	}
	*temp=tau;
	return 0;
}

int sub_region_hs(double spe_enth, double spe_entr)
{
	double n[4]={-0.349898083432139e4,0.257560716905876e4,-0.421073558227969e3,0.276349063799944e2};
	double h2ab;
	
	h2ab=n[0]+n[1]*spe_entr/1e3+n[2]*spe_entr*spe_entr/1e6+n[3]*spe_entr*spe_entr*spe_entr/1e9;
	
	if(spe_entr<5.85e3)
	{
		return 23;
	}
	else
	{
		if(spe_enth/1e3>h2ab)
		{
			return 22;
		}
		else
		{
			return 21;
		}
	} 
}

int pres_hs_r2(double* pressure, double spe_enth, double spe_entr)
{
	double eta2a,eta2b,eta2c,sig2a,sig2b,sig2c;
	int i;
	int sub_reg;
	double pi;
	double ni,ii,ji;
	
	eta2a=spe_enth/4200e3;
	eta2b=spe_enth/4100e3;
	eta2c=spe_enth/3500e3;
	sig2a=spe_entr/12e3;
	sig2b=spe_entr/7.9e3;
	sig2c=spe_entr/5.9e3;
	
	sub_reg=sub_region_hs(spe_enth,spe_entr);
//	printf("\n\n\\\\\\\\\\\\\\\\\\\\\\\\%d\n",sub_reg);
	if(sub_reg==21)
	{
		pi=0;
		for(i=0;i<29;i++)
		{
			ni=region2a_coeff_hs[i][3];
			ii=region2a_coeff_hs[i][1];
			ji=region2a_coeff_hs[i][2]; 
			pi+=ni*pow(eta2a-0.5,ii)*pow(sig2a-1.2,ji);
		}
		*pressure=pi*pi*pi*pi*4e6;
		return 0;
	}
	else if(sub_reg==22)
	{
		pi=0;
		for(i=0;i<33;i++)
		{
			ni=region2b_coeff_hs[i][3];
			ii=region2b_coeff_hs[i][1];
			ji=region2b_coeff_hs[i][2];
			pi+=ni*pow(eta2b-0.6,ii)*pow(sig2b-1.01,ji);
		}
		*pressure=pi*pi*pi*pi*1e8;
		return 0;
	}
	else if(sub_reg==23)
	{
		pi=0;
		for(i=0;i<31;i++)
		{
			ni=region2c_coeff_hs[i][3];
			ii=region2c_coeff_hs[i][1];
			ji=region2c_coeff_hs[i][2];
			pi+=ni*pow(eta2c-0.7,ii)*pow(sig2c-1.1,ji);
		}
		*pressure=pi*pi*pi*pi*1e8;
//		printf("\n[[[[[[[[[[[[[[[[[[[[[[%lf,%lf\n",pi,*pressure);
		return 0;
	}
	return -1;
}

void calcFT_dFT_prr2(double* ft, double* dft, double pressure, double dens, double tau)
{
	double pi;
	pi=pressure/1e6;
	*ft=gm0pi(pi)+gmrpi(pi,tau)-1e6*tau/(dens*GAS_CONST_STEAM*540);
	*dft=gmrpitau(pi,tau)-1e6/(dens*GAS_CONST_STEAM*540);
}

int temp_pr_r2(double* temp, double pressure, double dens)
{
	double tau_prev,tau_nxt,ft,dft;
	int iter_times;
	int ti=0;
	double tau_ini[10]={
		1.97,1.8,1.5,1.3,1.1,1,0.9,0.7,0.6,0.5,
	};
	
	do{
		iter_times=0;
		calcFT_dFT_prr2(&ft,&dft,pressure,dens,tau_ini[ti]);
		tau_nxt=tau_ini[ti]-ft/dft;
		do{
			tau_prev=tau_nxt;
			calcFT_dFT_prr2(&ft,&dft,pressure,dens,tau_prev);
			tau_nxt=tau_prev-ft/dft;
			iter_times++;
		}while(fabs(tau_nxt-tau_prev)>1e-6&&iter_times<MAX_ITER_NUMSR2);
		ti++;
		if(iter_times==MAX_ITER_NUMSR2)
		{
			continue;
		}
		else
		{
			*temp=540/tau_nxt;
			if(NRT(*temp)==0) 
			break;
		}			
	}while(ti<10);
	if(ti==10)
	{
		*temp=-1;
		return -1;
	}
	else
	{
		return 0;
	}
} 


void calcFT_dFT_pur2(double* ft, double* dft, double pressure, double spe_ener, double tau)
{
	double pi;
	pi=pressure/1e6;
	
	*ft=tau*(gm0tau(tau)+gmrtau(pi,tau))-pi*(gm0pi(pi)+gmrpi(pi,tau))-spe_ener*tau/(GAS_CONST_STEAM*540);
	*dft=gm0tau(tau)+gmrtau(pi,tau)+tau*(gm0tautau(tau)+gmrtautau(pi,tau))-pi*(gmrpitau(pi,tau))-spe_ener/(GAS_CONST_STEAM*540);
}

int temp_pu_r2(double* temp, double pressure, double spe_ener)
{
	double tau_prev,tau_nxt,ft,dft;
	int iter_times;
	int ti=0;
	double tau_ini[10]={
		1.97,1.8,1.5,1.3,1.1,1,0.9,0.7,0.6,0.5,
	};
	do{
		iter_times=0;
		calcFT_dFT_pur2(&ft,&dft,pressure,spe_ener,tau_ini[ti]);
//	printf("##:%lf,%lf\n",ft,dft);
		tau_nxt=tau_ini[ti]-ft/dft;
		do{
			tau_prev=tau_nxt;
			calcFT_dFT_pur2(&ft,&dft,pressure,spe_ener,tau_prev);
			tau_nxt=tau_prev-ft/dft;
			iter_times++;
//		printf("## %d,%lf,%lf\n",iter_times,ft/dft,540/tau_nxt); 
		}while(fabs(tau_nxt-tau_prev)>1e-6&&iter_times<MAX_ITER_NUMSR2);
		ti++;
		if(iter_times==MAX_ITER_NUMSR2)
		{
			continue;
		}
		else
		{
			*temp=540/tau_nxt;
			if(NRT(*temp)==0) 
		 	break;
		}
	}while(ti<10&&NRT(*temp)!=0);
	if(ti==10)
	{
		*temp=-1;
		return -1;
	}
	else
	{
		return 0;
	}
}

void calcFPi_dFPi_trr2(double* fp, double* dfp, double temp, double dens, double pi)
{
	double tau;
	tau=540/temp;
	*fp=dens*GAS_CONST_STEAM*temp*(gm0pi(pi)+gmrpi(pi,tau))-1e6;
	*dfp=dens*GAS_CONST_STEAM*temp*(gm0pipi(pi)+gmrpipi(pi,tau));
}

int pres_tr_r2(double* pres, double temp, double dens)
{
	double pi_prev,pi_nxt,fp,dfp;
	int iter_times;
	int ti=0;
	double pi_ini[13]={
		1,10,1e-1,1e-4,20,50,1e-3,1e-6,80,99,1e-2,1e-5,1e-7,
	};
	
	do{
		iter_times=0;
		calcFPi_dFPi_trr2(&fp,&dfp,temp,dens,pi_ini[ti]);
		pi_nxt=pi_ini[ti]-fp/dfp;
		do{
			pi_prev=pi_nxt;
			calcFPi_dFPi_trr2(&fp,&dfp,temp,dens,pi_prev);
			pi_nxt=pi_prev-fp/dfp;
			iter_times++;
		}while(fabs(pi_nxt-pi_prev)>1e-8&&iter_times<MAX_ITER_NUMSR2);
		ti++;
		if(iter_times==MAX_ITER_NUMSR2)
		{
//			printf("\n!!!!!!!!!!!!!");
			continue;
		}
		else
		{
			*pres=pi_nxt*1e6;
			if(NRP(*pres)==0)
//					printf("\n%d\n",NRP(*pres));
			break;
		}
//		printf("\n%d\n",NRP(*pres));
	}while(ti<13);
	if(ti==13)
	{
		*pres=-1;
		return -1;
	}
	else
	{
		return 0;
	}
}


void calcFPi_dFPi_tur2(double* fp, double* dfp, double temp, double spe_ener, double pi)
{
	double tau;
	tau=540/temp;
	*fp=tau*(gm0tau(tau)+gmrtau(pi,tau))-pi*(gm0pi(pi)+gmrpi(pi,tau))-spe_ener*tau/(GAS_CONST_STEAM*540);
	*dfp=tau*gmrpitau(pi,tau)-pi*(gm0pipi(pi)+gmrpipi(pi,tau))-(gm0pi(pi)+gmrpi(pi,tau));
}

int pres_tu_r2(double* pres, double temp, double spe_ener)
{
	double pi_prev,pi_nxt,fp,dfp;
	int iter_times;
	int ti=0;
	double pi_ini[13]={
		1,10,1e-1,1e-4,20,50,1e-3,1e-6,80,99,1e-2,1e-5,1e-7,
	};
	do{
		iter_times=0;
		calcFPi_dFPi_tur2(&fp,&dfp,temp,spe_ener,pi_ini[ti]);
		pi_nxt=pi_ini[ti]-fp/dfp;
		do{
			pi_prev=pi_nxt;
			calcFPi_dFPi_tur2(&fp,&dfp,temp,spe_ener,pi_prev);
			pi_nxt=pi_prev-fp/dfp;
			iter_times++;
		}while(fabs(pi_nxt-pi_prev)>1e-8&&iter_times<MAX_ITER_NUMSR2);
		ti++;
		if(iter_times==MAX_ITER_NUMSR2)
		{
			continue;
		}
		else
		{
			*pres=pi_nxt*1e6;
			if(NRP(*pres)==0)
			break;
		}	
	}while(ti<13);
	if(ti==13)
	{
		*pres=-1;
		return -1;
	}
	else
	{
		return 0;
	}
}

void calcFPi_dFPi_thr2(double* fp, double* dfp, double temp, double spe_enth, double pi)
{
	double tau;
	tau=540/temp;
	*fp=tau*(gm0tau(tau)+gmrtau(pi,tau))-spe_enth*tau/(GAS_CONST_STEAM*540);
	*dfp=tau*gmrpitau(pi,tau);
}

int pres_th_r2(double* pres, double temp, double spe_enth)
{
	double pi_prev,pi_nxt,fp,dfp;
	int iter_times;
	int ti=0;
	double pi_ini[13]={
		1,10,1e-1,1e-4,20,50,1e-3,1e-6,80,99,1e-2,1e-5,1e-7,
	};
	
	do{
		iter_times=0;
		calcFPi_dFPi_thr2(&fp,&dfp,temp,spe_enth,pi_ini[ti]);
		pi_nxt=pi_ini[ti]-fp/dfp;
		do{
			pi_prev=pi_nxt;
			calcFPi_dFPi_thr2(&fp,&dfp,temp,spe_enth,pi_prev);
			pi_nxt=pi_prev-fp/dfp;
			iter_times++;
		}while(fabs(pi_nxt-pi_prev)>1e-8&&iter_times<MAX_ITER_NUMSR2);
		ti++;
		if(iter_times==MAX_ITER_NUMSR2)
		{
			continue;
		}
		else
		{
			*pres=pi_nxt*1e6;
			if(NRP(*pres)==0)
			break;
		}	
	}while(ti<13);
	if(ti==13)
	{
		*pres=-1;
		return -1;
	}
	else
	{
		return 0;
	}
}

void calcFPi_dFPi_tsr2(double* fp, double* dfp, double temp, double spe_entr, double pi)
{
	double tau;
	tau=540/temp;
	*fp=tau*(gm0tau(tau)+gmrtau(pi,tau))-gm0(pi,tau)-gmr(pi,tau)-spe_entr/GAS_CONST_STEAM;
	*dfp=tau*gmrpitau(pi,tau)-gm0pi(pi)-gmrpi(pi,tau);
}

int pres_ts_r2(double* pres, double temp, double spe_entr)
{
	double pi_prev,pi_nxt,fp,dfp;
	int iter_times;
	int ti=0;
	double pi_ini[13]={
		1,10,1e-1,1e-4,20,50,1e-3,1e-6,80,99,1e-2,1e-5,1e-7,
	};
	
	do{
		iter_times=0;
		calcFPi_dFPi_tsr2(&fp,&dfp,temp,spe_entr,pi_ini[ti]);
//		printf("##fp,dfp:::%lf,%lf\n",fp,dfp);
		pi_nxt=pi_ini[ti]-fp/dfp;
		do{
			pi_prev=pi_nxt;
			calcFPi_dFPi_tsr2(&fp,&dfp,temp,spe_entr,pi_prev);
			pi_nxt=pi_prev-fp/dfp;
			iter_times++;
		}while(fabs(pi_nxt-pi_prev)>1e-8&&iter_times<MAX_ITER_NUMSR2);
		ti++;
		
//		printf("#####%lf\n",pi_nxt);
		if(iter_times==MAX_ITER_NUMSR2)
		{
			continue;
		}
		else
		{
			*pres=pi_nxt*1e6;
			if(NRP(*pres)==0)
			break;
		}		
	}while(ti<13);
	if(ti==13)
	{
		*pres=-1;
		return -1;
	}
	else
	{
		return 0;
	}
}

int ptini_hs_r2(double* presini, double* tempini, double spe_enth, double spe_entr)
{
	double p,t;
	if(pres_hs_r2(&p,spe_enth,spe_entr)==-1)
	{
		return -1;
	}
	if(temp_ph_r2(&t,p,spe_enth)==-1)
	{
		return -1;
	}
//	printf("\n..............%lf,%lf\n",p,t);
	*presini=p;
	*tempini=t;
	return 0;
}

void calcFG_dFG_hs_r2(double* f, double* g, double*dfdp, double* dfdt, double* dgdp, double* dgdt, double spec_enth, double spec_entr, double pi, double tau)
{
	*f=gm0tau(tau)+gmrtau(pi,tau)-spec_enth/(GAS_CONST_STEAM*540);
	*g=tau*(gm0tau(tau)+gmrtau(pi,tau))-(gm0(pi,tau)+gmr(pi,tau))-spec_entr/GAS_CONST_STEAM;
	*dfdp=gmrpitau(pi,tau);
	*dfdt=gm0tautau(tau)+gmrtautau(pi,tau);
	*dgdp=tau*gmrpitau(pi,tau)-gm0pi(pi)-gmrpi(pi,tau);
	*dgdt=tau*gm0tautau(tau)+tau*gmrtautau(pi,tau);
}

int pt_hs_r2(double* pres, double* temp, double spe_enth, double spe_entr)
{
	double pres_ini,pi_ini,pi_prev,pi_nxt;
	double temp_ini,tau_ini,tau_prev,tau_nxt;
	int i=0;
	double f,g,dfdt,dgdt,dfdp,dgdp;
	double dp,dt,maxd;
	
	ptini_hs_r2(&pres_ini,&temp_ini,spe_enth,spe_entr);
	
// 	printf("###%.12lf\n###%.12lf\n",temp_ini,pres_ini);
	pi_ini=pres_ini/1e6;
	tau_ini=540/temp_ini;
	
	calcFG_dFG_hs_r2(&f,&g,&dfdp,&dfdt,&dgdp,&dgdt,spe_enth,spe_entr,pi_ini,tau_ini);
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
		calcFG_dFG_hs_r2(&f,&g,&dfdp,&dfdt,&dgdp,&dgdt,spe_enth,spe_entr,pi_prev,tau_prev);
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
	}while(maxd>1e-8&&i<MAX_ITER_NUMSR2);
	if(i==MAX_ITER_NUMSR2)
	{
		*pres=pres_ini;
		*temp=temp_ini;
		return -2;
	}	
	else
	{
		*pres=pi_nxt*1e6;
		*temp=540/tau_nxt;
//		printf("\n\t......%lf,%lf\n",*pres,*temp);
		return 0;
	}	
}