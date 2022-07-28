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
alpha_ch4, alpha_n20, m0, n0 = params.radiative_params
bioCH4 = params.bioCH4  #this is in megatonnes
init_ppm = params.ppmCO2_1750
tao, beta, Cpf_moles, propCH4, propPassive, pf_init_T = params.permafrost_params
fds, fld, fls, flv, fnd, fns, fnv, frd, frs, fvd, fvs, beta_biota, q10, biota_init_T, npp_0 = params.biota_params
k= params.k#for deltaN of CO2


bioCH4ppb = bioCH4 * 1e12/ 16 *  1e9/AM

ppmtoMol = 1.77e+14 #mult by this for CO2 PPM -> Moles



def fmn(m, n):
  result = 0.47 * np.log(1 + (2.01e-5 * ((m * n)**0.75)) + (5.31e-15 * m * ((m * n)**1.52)))
  return result


def deltaN_SOx(SOx):
    # val = 53.8412
    val = 56.2878
    direct = (-0.4 * 3.154e+07) * (SOx/val)
    indirect = (-0.8 * 3.154e+07) * np.log((42 + SOx)/42) * (np.log ((42 + val)/42))**-1
    return (direct + indirect)


def NPP(carbon): #note in ppm 
	result = (npp_0 * 1e15 / 12) * (1 + beta_biota * ( np.log ( carbon /  params.ppmCO2_1750 )))
	return result

def RH_soil(CS, avg_T):
	result = CS * frs * (q10 ** (avg_T/10))
	return result
	
def RH_det(CD, avg_T):
	result = CD * frd * (q10 ** (avg_T/10))
	return result
	


class scenario:
	'''
	Creates a class that returns various emission scenarios
	
	'''
	def __init__(self, tspan, state, temp_array, RCP, custom_ch4 = False, reduction = .9, baseline = 142.0527):
		self.__tspan = tspan
		self.__state = state
		self.__RCP_emit_df = RCP #df
		self.__reduction = reduction #value 0-1
		self.__custom_ch4 = custom_ch4 #boolean
		self.__baseline = baseline
		if self.__custom_ch4 == True:
			new_emit_meth = np.ones(self.__RCP_emit_df["CH4"].size) * baseline
			index_start = 245
			emissions = []
			meth_emit = 330 # megatonnes current
			emissions.append(meth_emit)
			while meth_emit > baseline: #RCP number
			  meth_emit = meth_emit * self.__reduction #reduction by % each year.  
			  emissions.append(meth_emit)

			for i in range(len(emissions)):
			  new_emit_meth[i+index_start] = emissions[i]
		
			self.__RCP_emit_df["CH4"] = new_emit_meth
	
		#carbon emission moles 
		self.__new_emitC = np.zeros([10000])
		self.__new_emitC[self.__RCP_emit_df.Year] = (self.__RCP_emit_df.Fossil + self.__RCP_emit_df.LUC) *  1e15/12 #moles future
		self.__LUC= np.zeros([10000])
		self.__LUC[self.__RCP_emit_df.Year] = self.__RCP_emit_df.LUC  *  1e15/12
		
		#set ppb CH4 by RCP or custom
		self.__emitM = np.zeros(10000)
		self.__emitM[self.__RCP_emit_df.Year] = (self.__RCP_emit_df.CH4) * 1e12/ 16 *  1e9/AM #convert to moles then convert to ppb
		#set ppb N2O by RCP
		self.__emitN = np.zeros(10000)
		self.__emitN[self.__RCP_emit_df.Year]= self.__RCP_emit_df.N2O *  1e12 / 44.013 * 1e9 / AM # convert to moles then convert to ppb
		# Sulf oxide emissions
		self.__emitSOx = np.zeros(10000)
		self.__emitSOx[self.__RCP_emit_df.Year] = self.__RCP_emit_df.SOx 

		self.__H_track = np.zeros(10000)
		self.__H_trackPF = np.zeros(10000)
		
		#permafrost temperature tracker of temp
		self.__temp_trackOG = temp_array
		self.__temp_track = np.copy(self.__temp_trackOG)

		#permafrost tracker of emissions
		self.__emitPF = np.zeros((10000, 2))
		
		# radiative forcing tracker
		self.__rf = np.zeros((10000,5))
		self.get_scenario()
		
	
	def get_RCP(self):
		return self.__RCP_emit_df
	def __getState(self):
		return self.__state
	
	def __get_tspan(self):
		return self.__tspan
		
	def getMethaneEmit(self):
	    return self.__emitM
	def getN2OEmit(self):
	    return self.__emitN
		
	def getH_track(self):
		return self.__H_track
	
	def getRadiativeForcing(self):
		return self.__rf
	
	def getEmitPF(self):
		return self.__emitPF
		
	def get_temp_track(self):
		return self.__temp_track
	def getFinalTemp(self):
	    return self.get_temp_track()[self.__get_tspan()[1]]
	    
	def __PF_extent(self, time):
		extent = (1- beta * ( self.__temp_track[time] - self.__temp_track[2000] ))
		return extent
		
	def __reset_temp(self):
	    self.__temp_track = np.copy(self.__temp_trackOG)
	
	def get_scenarioPF(self):
		self.__reset_temp()
		sol1 = solve_ivp(self.__BEAMfuturePF, t_span=self.__tspan, y0 = self.__state, args=(self.__new_emitC, self.__LUC, self.__emitM, self.__emitN, self.__emitSOx), max_step = .5)
		return sol1

		
	def get_scenario(self):
		self.__reset_temp()
		sol = solve_ivp(self.__BEAMfuture, t_span=self.__tspan, y0 = self.__state[0:10], args=(self.__new_emitC, self.__LUC, self.__emitM, self.__emitN, self.__emitSOx), max_step = .5)
		return sol
		
		
		
	def __BEAMfuture(self, t, y, emissionsC, LUC,  emissionsM, emissionsN, emissionsSOx):

	  tMix, tDeep, QA, QU, QL, CV, CD, CS, ch4, n2o= y 
	  
	  t_yr = int(np.floor(t))

	  #Get pH first
	  #Everything is in moles and mole fraction, so we are good:
	  a = 1;
	  b = K_1 * (1 - QU / Alk);
	  c = K_1*K_2 * (1 - 2*QU / Alk);
	  #Take H as positive root
	  H = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a);
	  #Also track H
	  self.__H_track[t_yr] = H

	  #And now our Lambda factor:
	  Lambda = 1 + K_1/H + K_1*K_2/H**2;  

	  currentppm = QA / ppmtoMol
	  # ppm_array[t_yr] = currentppm
	  m = ch4 #ppb
	  n = n2o #ppb 

	  cERF = k * np.log(currentppm / init_ppm) * 3.154e+07 # J/yr conversion
	  mERF = ((alpha_ch4 * (m**.5 - m0**.5)) - (fmn(m,n0) - fmn(m0,n0)) ) * 3.154e+07
	  nERF = ( (alpha_ch4 * (n**.5 - n0**.5)) - (fmn(m0,n) - fmn(m0,n0)) ) * 3.154e+07 
	  sERF = deltaN_SOx(emissionsSOx[t_yr]) * params.riley_param
	  
	  deltaN = cERF + mERF + nERF + sERF
	  
	  self.__rf[t_yr,:] = np.array([cERF, mERF, nERF, sERF, deltaN])

	  Emit = emissionsC[t_yr];
	  emit_methane = emissionsM[t_yr] # note these are in ppb
	  emit_n20 = emissionsN[t_yr] #note these are in ppb
	  avgT = np.average(self.__temp_track[t_yr-11: t_yr])
	  NPP_val = NPP(currentppm)


	  dydt = [ (-1 * lambdaHat * tMix - (gamma * (tMix - tDeep)) + deltaN )/ mumix , #change in mixing layer
			  (gamma * (tMix - tDeep)) / mudeep, #change in deep layer
			  -ka * QA + ka * kH * QU / (delta_a * Lambda) + Emit + RH_det(CD, avgT) + RH_soil(CS, avgT) - NPP_val, #dQA/dt
			  ka * QA - ka * kH* QU / (delta_a * Lambda) - kd * QU + kd * QL / delta_d, #dQU / dt
			  kd*QU - kd* QL / delta_d,#dQL/dt
			  NPP_val * fnv - (CV * (fvd + fvs)) - LUC[t_yr] * flv,   #dCveg
			  NPP_val * fnd + (CV * fvd) - (CD * fds) - RH_det(CD, avgT) - LUC[t_yr]   * fld,     #dCdet 
			  NPP_val * fns + (CV * fvs)  + (CD * fds) - RH_soil(CS, avgT) - LUC[t_yr] * fls, 
			  emit_methane + bioCH4ppb - ch4/ 12, #decay rate of 1/12 a year methane
			  emit_n20 - n2o/114 ]    #dCsoil
			  
	  self.__temp_track[t_yr] = tMix + dydt[0]
			    
	  return dydt
		
		
	
	def __BEAMfuturePF(self, t, y, emissionsC, LUC, emissionsM, emissionsN, emissionsSOx):
	  tMix, tDeep, QA, QU, QL,  CV, CD, CS, ch4, n2o,  Cpf, Lc, Lc_m= y 
	  
	  t_yr = int(np.floor(t))

	  #Get pH first
	  #Everything is in moles and mole fraction, so we are good:
	  a = 1;
	  b = K_1 * (1 - QU / Alk);
	  c = K_1*K_2 * (1 - 2*QU / Alk);
	  #Take H as positive root
	  H = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a);
	  #Also track H
	  self.__H_track[t_yr] = H

	  #And now our Lambda factor:
	  Lambda = 1 + K_1/H + K_1*K_2/H**2;  

	  currentppm = QA / ppmtoMol
	  # ppm_array[t_yr] = currentppm
	  m = ch4 #ppb
	  n = n2o #ppb 

	  cERF = k * np.log(currentppm / init_ppm) * 3.154e+07 # J/yr conversion
	  mERF = ((alpha_ch4 * (m**.5 - m0**.5)) - (fmn(m,n0) - fmn(m0,n0)) ) * 3.154e+07
	  nERF = ( (alpha_ch4 * (n**.5 - n0**.5)) - (fmn(m0,n) - fmn(m0,n0)) ) * 3.154e+07 
	  
	  sERF = deltaN_SOx(emissionsSOx[t_yr]) * params.riley_param
	  deltaN = cERF + mERF + nERF + sERF; 
	  
	  self.__rf[t_yr,:] = np.array([cERF, mERF, nERF, sERF, deltaN]);
	  
	  Emit = emissionsC[t_yr];
	  emit_methane = emissionsM[t_yr] # note these are in ppb
	  emit_n20 = emissionsN[t_yr] #note these are in ppb
	  
	  #permafrost stuff
	  # shoudl be negative 
	  change_extent = self.__PF_extent(t_yr) - self.__PF_extent(t_yr-1)
	  if change_extent > 0:
	      change_extent = 0
	  dCpfdt =  Cpf * change_extent * (1 -propPassive) * (1 - propCH4) +  (Cpf * change_extent * (1-propPassive) * (propCH4))
	  dLc = Cpf * -1 * change_extent * (1 -propPassive) * (1 - propCH4) - 1 / tao * Lc 
	  dLc_m = Cpf * -1 * change_extent * (1-propPassive) * (propCH4) - 1 / tao * Lc_m
	  # print(PF_extent(t_yr) - PF_extent(t_yr-1))

	  # print(dCpfdt)
	  PF_emit_co2 = 1 /tao * Lc 
	  PF_emit_ch4 = 1/tao * Lc_m  # should be moles of methane 
	  self.__emitPF[t_yr,:] = np.add(self.__emitPF[t_yr - 1, : ], np.array([PF_emit_co2, PF_emit_ch4]) )
	  
	  PF_emit_ch4_ppb = PF_emit_ch4 *  1e9 / AM #convert to parts per billion
	  
	  avgT = np.average(self.__temp_track[t_yr-11: t_yr])
	  NPP_val = NPP(currentppm)
	  
	  dydt = [ (-1 * lambdaHat * tMix - (gamma * (tMix - tDeep)) + deltaN )/ mumix , #change in mixing layer
			  (gamma * (tMix - tDeep)) / mudeep, #change in deep layer
			  -ka * QA + ka * kH * QU / (delta_a * Lambda) + Emit + + RH_det(CD, avgT) + RH_soil(CS, avgT) - NPP_val  +PF_emit_co2, #dQA/dt
			  ka * QA - ka * kH* QU / (delta_a * Lambda) - kd * QU + kd * QL / delta_d, #dQU / dt
			  kd*QU - kd* QL / delta_d,#dQL/dt
			  NPP_val * fnv - (CV * (fvd + fvs)) - LUC[t_yr] * flv,   #dCveg
			  NPP_val * fnd + (CV * fvd) - (CD * fds) - RH_det(CD, avgT) - LUC[t_yr]   * fld,     #dCdet 
			  NPP_val * fns + (CV * fvs)  + (CD * fds) - RH_soil(CS, avgT) - LUC[t_yr] * fls,
			  emit_methane + bioCH4ppb  + PF_emit_ch4_ppb - ch4/ 12, #decay rate of 1/12 a year methane #PPB
			  emit_n20 -  n2o/114, #change  in N2o  PPB
			  dCpfdt, # d Cpf /dt
			  dLc,  # Labile carbon 
			  dLc_m] # Labile carbon as methane 
			  
	  self.__temp_track[t_yr] = tMix + dydt[0]

	  return dydt