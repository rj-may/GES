from scipy.integrate import solve_ivp
import numpy as np
import params
from params import ocean_params

ka, kd, delta_d, kH, AM, OM, M_upper, delta_a, K_1, K_2, Alk = ocean_params
fds, fld, fls, flv, fnd, fns, fnv, frd, frs, fvd, fvs, beta_biota, q10, biota_init_T, npp_0 = params.biota_params


#Define a variable H_track to track pH
#H_track in mole fraction
H_track = np.zeros(10000)

def RH_soil(CS, avg_T):
	result = CS * frs * (q10 ** (avg_T/10))
	return result
	
def RH_det(CD, avg_T):
	result = CD * frd * (q10 ** (avg_T/10))
	return result




def beam_init(t, Q):
  QA, QU, QL , CV, CD, CS = Q
  # print(t)

  t_yr = int(np.floor(t))
  #Get pH first
  #Everything is in moles and mole fraction, so we are good:
  a = 1;
  b = K_1 * (1 - QU / Alk);
  c = K_1*K_2 * (1 - 2*QU / Alk);
  #Take H as positive root
  H = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a);
  #Also track H
  H_track[t_yr] = H
  # print("This")
  H = H_track[t_yr]
  
  #And now our Lambda factor:
  Lambda = 1 + K_1/H + K_1*K_2/H**2;
  avgT = 0
  NPP = (npp_0 * 1e15 / 12)
  dQdt = [0, #dQA/dt assumes no change for carbon levels
		ka * QA - ka * kH* QU / (delta_a * Lambda) - kd * QU + kd * QL / delta_d, #dQU / dt
		kd*QU - kd* QL / delta_d, #
		NPP * fnv - (CV * (fvd + fvs)) - 0 * flv,   #dCveg 
		NPP * fnd + (CV * fvd) - (CD * fds) - RH_det(CD, avgT) - 0   * fld,     #dCdet 
		NPP * fns + (CV * fvs)  + (CD * fds) - RH_soil(CS, avgT) - 0 * fls] #d
  
  return dQdt
  
  
# initial ppm for the data set- pre industrial 
# init_ppm = ppm_data.loc[ppm_data['Year'] == 1750]["CO2"].values[0]



def get_C_init(init_ppm):
	
	init_moles = init_ppm / 1e6 * AM 
	init_GTC_QA = init_moles * 12 / (1e15)
	# init_GTC_QA = init_ppm * 7.788

	# Q0_GTC = np.array([init_GTC_QA, 713, 35658]) #GTC of carbon
	carbon_GTC = np.array([init_GTC_QA, 500, 35400, 550, 55, 1782]) #GTC of carbon

	carbon_moles = carbon_GTC  * (1e15) /12

	max_t= 5000
	sol = solve_ivp(beam_init, t_span=[0,max_t], y0 =carbon_moles, max_step =.5)
	
	
	last_index = len(sol.y[0])-1

	initC = sol.y[:,last_index]
	
	return initC
	
