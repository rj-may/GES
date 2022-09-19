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


# bioCH4ppb = bioCH4 * 1e12/ 16 *  1e9/AM
bioCH4moles = bioCH4 * 1e12/ 16

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
		self.__fixState()
		self.__RCP_emit_df = RCP #df
		self.__reduction = reduction #value 0-1
		self.__custom_ch4 = custom_ch4 #boolean
		self.__baseline = baseline
		if self.__custom_ch4 == True:
			new_emit_meth = np.ones(self.__RCP_emit_df["CH4"].size) * self.__baseline
			index_start = 245 #corresponds with year 2010
			emissions = []
			meth_emit = 330 # megatonnes current
			while meth_emit > baseline and len(emissions) <= self.__RCP_emit_df["CH4"].size - index_start: #RCP number
			  emissions.append(meth_emit)
			  meth_emit = meth_emit * self.__reduction #reduction by % each year.  

			for i in range(len(emissions)-1):
			  new_emit_meth[i+index_start] = emissions[i]
		
			self.__RCP_emit_df["CH4"] = new_emit_meth
	
		#carbon emission moles 
		self.__new_emitC = np.zeros([10000])
		self.__new_emitC[self.__RCP_emit_df.Year] = (self.__RCP_emit_df.Fossil + self.__RCP_emit_df.LUC) *  1e15/12 #moles future
		self.__LUC= np.zeros([10000])
		self.__LUC[self.__RCP_emit_df.Year] = self.__RCP_emit_df.LUC  *  1e15/12
		
		#set CH4 by RCP or custom
		self.__emitM = np.zeros(10000)
		# self.__emitM[self.__RCP_emit_df.Year] = (self.__RCP_emit_df.CH4) * 1e12/ 16 *  1e9/AM #convert to moles then convert to ppb
		self.__emitM[self.__RCP_emit_df.Year] = (self.__RCP_emit_df.CH4) * 1e12/ 16  # just moles
		#set  N2O by RCP
		self.__emitN = np.zeros(10000)
		# self.__emitN[self.__RCP_emit_df.Year]= self.__RCP_emit_df.N2O *  1e12 / 44.013 * 1e9 / AM # convert to moles then convert to ppb
		self.__emitN[self.__RCP_emit_df.Year]= self.__RCP_emit_df.N2O *  1e12 / 44.013 #just moles

		# Sulf oxide emissions
		self.__emitSOx = np.zeros(10000)
		self.__emitSOx[self.__RCP_emit_df.Year] = self.__RCP_emit_df.SOx 

		self.__H_track = np.zeros(10000)
		# self.__H_trackPF = np.zeros(10000)
		
		#permafrost temperature tracker of temp
		self.__temp_trackOG = temp_array
		self.__temp_track = np.copy(self.__temp_trackOG)

		#permafrost tracker of emissions
		self.__emitPF = np.zeros((10000, 2))
		
		self.__meth_track = np.zeros(10000)
		# radiative forcing tracker
		self.__rf = np.zeros((10000,5))
		
		self.get_scenarioPF()
		
	
	def __fixState(self):
		fixing = np.copy(self.__getState())
		fixing[8] = fixing[8] * AM/1e9 #Converting back to moles
		fixing[9] = fixing[9] * AM/1e9 #converting back to moles
		self.__state= fixing

	def get_RCP(self):
		return self.__RCP_emit_df
	def __getState(self):
		return self.__state
	
	def __get_tspan(self):
		return self.__tspan
		
	def getMethaneEmit(self):
	    return self.__emitM * AM /1e9 * 16 / 1e12
	    
	def getN2OEmit(self):
	    return self.__emitN  / 1e12 * 44.013 / 1e9 * AM
		
	def getCarbonEmit(self):
	    return self.__new_emitC * 12 / 1e15
	    
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

	def getTemp2100(self):
	    return self.get_temp_track()[2100]
	    
	def __PF_extent(self, time):
		extent = (1- beta * ( self.__temp_track[time] - self.__temp_track[2010] ))
		return extent
	
	def __frozen(self, time):
	    cnew = Cpf_moles * self.__PF_extent(time)
	    return cnew
	
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
	  m = ch4 *1e9/AM #ppb
	  n = n2o  *1e9/AM #ppb 

	  cERF = k * np.log(currentppm / init_ppm) * 3.154e+07 # J/yr conversion
	  mERF = ((alpha_ch4 * (m**.5 - m0**.5)) - (fmn(m,n0) - fmn(m0,n0)) ) * 3.154e+07
	  nERF = ( (alpha_ch4 * (n**.5 - n0**.5)) - (fmn(m0,n) - fmn(m0,n0)) ) * 3.154e+07 
	  sERF = deltaN_SOx(emissionsSOx[t_yr]) * params.riley_param
	  
	  deltaN = cERF + mERF + nERF + sERF
	  
	  self.__rf[t_yr,:] = np.array([cERF, mERF, nERF, sERF, deltaN])

	  Emit = emissionsC[t_yr];
	  emit_methane = emissionsM[t_yr] 
	  emit_n20 = emissionsN[t_yr] 
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
			  emit_methane + bioCH4moles - ch4/ 12, #decay rate of 1/12 a year methane
			  emit_n20 - n2o/114 ]    #dCsoil
			  
	  self.__temp_track[t_yr] = tMix + dydt[0]
			    
	  return dydt
		
		
	
	def __BEAMfuturePF(self, t, y, emissionsC, LUC, emissionsM, emissionsN, emissionsSOx):
	  tMix, tDeep, QA, QU, QL,  CV, CD, CS, ch4, n2o,  Cpf, Lc, Lm= y 
	  
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
	  m = ch4 *1e9/AM #ppb
	  n = n2o *1e9/AM #ppb 

	  cERF = k * np.log(currentppm / init_ppm) * 3.154e+07 # J/yr conversion
	  mERF = ((alpha_ch4 * (m**.5 - m0**.5)) - (fmn(m,n0) - fmn(m0,n0)) ) * 3.154e+07
	  nERF = ( (alpha_ch4 * (n**.5 - n0**.5)) - (fmn(m0,n) - fmn(m0,n0)) ) * 3.154e+07 
	  
	  sERF = deltaN_SOx(emissionsSOx[t_yr]) * params.riley_param
	  deltaN = cERF + mERF + nERF + sERF; 
	  
	  self.__rf[t_yr,:] = np.array([cERF, mERF, nERF, sERF, deltaN]);
	  
	  Emit = emissionsC[t_yr];
	  emit_methane = emissionsM[t_yr] 
	  emit_n20 = emissionsN[t_yr]
	  
	  #permafrost stuff
	  # shoudl be negative 
	  B = (self.__frozen(t_yr) - Cpf )
	  if B > 0:
	      B = 0
	  alpha = 1
	  dCpfdt = alpha * B
	  dLc = -B * (1 - propPassive) * (1- propCH4) - 1 / tao * Lc 
	  dLm = -B * (1-propPassive) * (propCH4) - 1 / tao * Lm
	  # print(PF_extent(t_yr) - PF_extent(t_yr-1))

	  self.__meth_track[t_yr] = dLm  + emit_methane * .25 #track new methane so it can be converted to carbon. #abritarily say that 1/4 of methane is new carbon from petroleum 

	  carbon_from_meth = self.__meth_track[t_yr -12]
	#   carbon_from_meth = 0


	  # print(dCpfdt)
	  PF_emit_co2 = 1 /tao * Lc 
	  PF_emit_ch4 = 1/tao * Lm  # should be moles of methane 
	  self.__emitPF[t_yr,:] = np.add(self.__emitPF[t_yr - 1, : ], np.array([PF_emit_co2, PF_emit_ch4]) )
	  
	#   PF_emit_ch4_ppb = PF_emit_ch4 *  1e9 / AM #convert to parts per billion
	  
	  avgT = np.average(self.__temp_track[t_yr-11: t_yr])
	  NPP_val = NPP(currentppm)
	  
	  dydt = [ (-1 * lambdaHat * tMix - (gamma * (tMix - tDeep)) + deltaN )/ mumix , #change in mixing layer
			  (gamma * (tMix - tDeep)) / mudeep, #change in deep layer
			  -ka * QA + ka * kH * QU / (delta_a * Lambda) + Emit + + RH_det(CD, avgT) + RH_soil(CS, avgT) - NPP_val  +PF_emit_co2 + carbon_from_meth, #dQA/dt
			  ka * QA - ka * kH* QU / (delta_a * Lambda) - kd * QU + kd * QL / delta_d, #dQU / dt
			  kd*QU - kd* QL / delta_d,#dQL/dt
			  NPP_val * fnv - (CV * (fvd + fvs)) - LUC[t_yr] * flv,   #dCveg
			  NPP_val * fnd + (CV * fvd) - (CD * fds) - RH_det(CD, avgT) - LUC[t_yr]   * fld,     #dCdet 
			  NPP_val * fns + (CV * fvs)  + (CD * fds) - RH_soil(CS, avgT) - LUC[t_yr] * fls,
			  emit_methane + bioCH4moles  + PF_emit_ch4 - ch4/ 12, #decay rate of 1/12 a year methane 
			  emit_n20 -  n2o/114, 
			  dCpfdt, # d Cpf /dt
			  dLc,  # Labile carbon 
			  dLm] # Labile carbon as methane 
			  
	  self.__temp_track[t_yr] = tMix + dydt[0]

	  return dydt
