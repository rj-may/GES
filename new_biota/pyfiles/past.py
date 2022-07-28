'''
Note that the purpose of this file is to calculate historical emissions with respect to temperature.

This file uses ppm data as a driver for the amount of carbon in the atmosphere
'''
from scipy.integrate import solve_ivp
import numpy as np
import pandas as pd
import params

ka, kd, delta_d, kH, AM, OM, M_upper, delta_a, K_1, K_2, Alk = params.ocean_params
mumix, mudeep, lambdaHat, gamma = params.temp_params
fds, fld, fls, flv, fnd, fns, fnv, frd, frs, fvd, fvs, beta_biota, q10, biota_init_T, npp_0 = params.biota_params

# These two things need to be accessed in multiple places in this file, hence they are global 
Imposed_CO2 = None #np array
ppm_data_global  = None #pd dataframe
LUC = None
Emissions = None
# RF_global = None

ppmtoMol = 1.77e+14 #mult by this for CO2 PPM -> Moles
molCO2toGT =  4.4e-14 #mult by this for moles CO2 -> GTC
ppmtoGT = 7.788 #mult by this for ppm to GTC
molCtoGT = 12/1e15

H_track = np.zeros(10000)
temp_track = np.zeros(10000)


# drive with CO2 means that we will use projected ppm data
def BEAM(t,y, drive_CO2 = True):
	tMix, tDeep, QA, QU, QL, CV, CD, CS = y 
	t_yr = int(np.floor(t))
	if drive_CO2:
		QA = Imposed_CO2[t_yr]
		deltaN = np.interp(t, ppm_data_global["Year"], ppm_data_global["deltaNco2"])
	if drive_CO2 == False:
		currentPPM = QA / AM * 1e6
		k = params.k
		deltaN = k * 3.154e+07 * (np.log(currentPPM / params.ppmCO2_1750))

  #Get pH first
  #Everything is in moles and mole fraction, so we are good:
	a = 1;
	b = K_1 * (1 - QU / Alk);
	c = K_1*K_2 * (1 - 2*QU / Alk);
	#Take H as positive root
	H = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a);
	#	Also track H
	H_track[t_yr] = H

	if t_yr == 2000:
		pf_init_T =  tMix

  #And now our Lambda factor:
	Lambda = 1 + K_1/H + K_1*K_2/H**2;  

  
#   deltaN = np.interp(t, RF_global["Year"], RF_global["SUM"])  * 3.154e+07
	Emit = Emissions[t_yr] #this is really a relic that we don't need because we are imposing
  
	avgT = np.average(temp_track[t_yr-11: t_yr])
	currentPPM = QA / AM * 1e6 

	dQAdt = -ka * QA + ka * kH * QU / (delta_a * Lambda) + Emit + RH_det(CD, avgT) + RH_soil(CS, avgT) - NPP(currentPPM)
  
	dydt = [ (-1 * lambdaHat * tMix - (gamma * (tMix - tDeep)) + deltaN )/ mumix , #change in mixing layer
          (gamma * (tMix - tDeep)) / mudeep, #change in deep layer
          dQAdt , #dQA/dt
          ka * QA - ka * kH * QU / (delta_a * Lambda) - kd * QU + kd * QL / delta_d, #dQU / dt
          kd*QU - kd* QL / delta_d, #dQL/dt
          NPP(currentPPM) * fnv - (CV * (fvd + fvs)) - LUC[t_yr] * flv,   #dCveg
          NPP(currentPPM) * fnd + (CV * fvd) - (CD * fds) - RH_det(CD, avgT) - LUC[t_yr]   * fld,     #dCdet
          NPP(currentPPM) * fns + (CV * fvs)  + (CD * fds) - RH_soil(CS, avgT) - LUC[t_yr] * fls]    #dCsoil

	temp_track[t_yr] = tMix + dydt [0]

	return dydt


def ppmDeltaNco2(ppm): #note this should ppm of CO2
	deltaNco2 = np.zeros(len(ppm))
	k= (5.35 * 3.154e+07)  #k =  J/yr conversion 
#  deltaN[0] = 0
	for i in range(1, len(ppm)):
		deltaNco2[i] = k * np.log(ppm[i] / ppm[0])
	return deltaNco2

def fmn(m, n):
  result = 0.47 * np.log(1 + (2.01e-5 * ((m * n)**0.75)) + (5.31e-15 * m * ((m * n)**1.52)))
  return result
  
def NPP(carbon):
	result = (npp_0 / molCtoGT) * (1 + beta_biota * ( np.log ( carbon /  params.ppmCO2_1750 )))
	return result

def RH_soil(CS, avg_T):
	result = CS * frs * (q10 ** (avg_T/10))
	return result
	
def RH_det(CD, avg_T):
	result = CD * frd * (q10 ** (avg_T/10))
	return result

  	

### This the main function!!!!!!!
  
def solve_past(tspan, inits, ppm_data, carbon_df, bool):
	'''
	This section of code is for calculating the various deltaN values
	'''
	# Adding delta N to the data frame for carbon
	deltaNco2 = ppmDeltaNco2(ppm_data["CO2"])
	
	ppm_data["deltaNco2"] = deltaNco2
	
	#add CH4 and N2O to the dataframe
	alpha_ch4 = 0.036
	alpha_n2o = 0.12
	m0= 720
	n0= 271

	mcol = ppm_data['CH4']
	ncol = ppm_data['N2O']

	ch4array = np.zeros(len(mcol))
	n2oarray = np.zeros(len(ncol))

	count = 0
	for m, n in zip(mcol, ncol): 
		ch4array[count] = (alpha_ch4 * (m**.5 - m0**.5)) - (fmn(m,n0) - fmn(m0,n0))
		n2oarray[count] = (alpha_n2o * (n**.5 - n0**.5)) - (fmn(m0,n) - fmn(m0,n0))
		count += 1
		
	ppm_data['deltaN_ch4'] = ch4array * 3.154e+07   #J/yr conversion
	ppm_data['deltaN_n2o'] = n2oarray * 3.154e+07

	cols = ['deltaNco2', 'deltaN_ch4', 'deltaN_n2o']
	ppm_data["deltaN"] = ppm_data[cols].sum(axis='columns')
		
	
	'''
	Carbon ppm levels This what drives the code
	'''

	#Imposed atmospheric CO2: Convert from ppm to moles 
	global Imposed_CO2 #unfortunately this is probs the best solution
	Imposed_CO2 = np.zeros([10000])
	Imposed_CO2[ppm_data.Year.astype(int)] = ppm_data.CO2 / 1e6 * AM
	
	global ppm_data_global
	ppm_data_global = ppm_data
	
# 	global RF_global
# 	RF_global = RF_df
	
	global LUC
	LUC = np.zeros(10000)
	LUC[carbon_df.Year] = carbon_df["LUC"]  * 1e15 / 12
	
	global Emissions
	Emissions = np.zeros([10000])
	Emissions[carbon_df.Year] = (carbon_df.Fossil + carbon_df.LUC) * 1e15/12
	
	x= solve_ivp(BEAM, t_span=tspan, y0 =inits, args=(bool,), max_step =1)
	
	return x, temp_track