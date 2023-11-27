/*
 * This code is distributed under the license: MIT License
 * Originally written by Zhenrong WANG
 * mailto: zhenrongwang@live.com
 */

#ifndef REGION_CALC_H_
#define REGION_CALC_H_

double calc_sat_pres(double temp);
double calc_bound23(double temp);
int region(double pressure, double temp);

#endif