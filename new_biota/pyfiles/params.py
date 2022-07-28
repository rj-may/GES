#Parameters

ka = 0.2 #yrs^-1
kd = 0.05 #yrs ^-1
delta_d = 50 # M_L / M_U

kH = 1.23 * 10 ** 3 #unitless

AM = 1.77e20
OM = 7.8e22

#Solve for delta_a, given AM, OM, and delta_d
M_upper = OM / (delta_d + 1)
delta_a = M_upper / AM

#DIC dissociation constants. Convert to mole fraction:
K_1 = 8e-7 *  18/1000
K_2 = 4.63e-10 * 18/1000

#And moles of alkalinity
Alk = 767e15 / 12 * 1.02 


ocean_params = [ka, kd, delta_d, kH, AM, OM, M_upper, delta_a, K_1, K_2, Alk]

fds = 0.6     #detritus to soil
fld = 0.01    #LUC  flux from detritus
fls = 0.89    # LUC from soil
flv = 0.10   #LUC from veg
fnd =  0.6  # npp to detritus
fns = 0.05 #npp to soil
fnv = 0.35  # npp to veg
frd = 0.25   #resp to det
frs = 0.02  # resp to soil
fvd = 0.034 #veg to det
fvs = 0.001 # veg to soil
beta_biota = 0.36 
# beta_biota = .5
q10 = 2.45
# q10 = 2
biota_init_T = 13.42
npp_0 = 50 #gtc

biota_params = [fds, fld, fls, flv, fnd, fns, fnv, frd, frs, fvd, fvs, beta_biota, q10, biota_init_T, npp_0]


mumix= 3.154e+08
mudeep= 6.307e+09
lambdaHat= (1.2 * 3.154e+07)
gamma= (1.2 * 3.154e+07)

temp_params = [mumix, mudeep, lambdaHat, gamma]


#add CH4 and N2O Delta N:
alpha_ch4 = 0.036
alpha_n2o = 0.12
m0=720 # 1750 level of ch4 ppb
n0= 271.2 # 1750 level of n2o ppb 

radiative_params = [alpha_ch4, alpha_n2o, m0, n0]

bioCH4 = 300 #Mt CH4
ppmCO2_1750 = 276.8 
k = 5.35 # W / m^2 for CO2 conversion

# riley_param = 0.46576551191432586
riley_param = 0




tao = 70 #mean 70, std 30
beta = 0.172 # std 0.0261
Cpf = 1035 #Gigatonnes of carbon
Cpf_moles = Cpf * 1e15 /12
propCH4 = 2.3/100
propPassive = 40/100

pf_init_T = 0.6936118298467074 #year 2000 temperature change  I took it from the running the model previously. This is the baseline year 
permafrost_params = [tao, beta,Cpf_moles, propCH4, propPassive, pf_init_T]



