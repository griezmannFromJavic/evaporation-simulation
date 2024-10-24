from funkcije import *
import CoolProp.CoolProp as CP

tlak_ekspanzije = 2e5

p_iz = 75e5
t_ul = 273.15 + 260 # cca 285
h_ul = CP.PropsSI('H', 'T', t_ul, 'P', p_iz, rt)
rho_ul = CP.PropsSI('D', 'T', t_ul, 'P', p_iz, rt)
u_ul = brzina_ekspanzije(tlak_ekspanzije, rho_ul)
grad_p = dpdl_trenja(p_iz, h_ul, u_ul * rho_ul, epsilon)
dp = grad_p * l_c
p_ul = p_iz + dp
ts_ul = t_ul
q_ul = 0

p_0 = np.linspace(p_ul, p_iz, n_l)
t_0 = np.zeros(n_l) + t_ul
rho_0 = CP.PropsSI('D', 'P', p_0, 'T', t_0, rt)
h_0 = CP.PropsSI('H', 'P', p_0, 'T', t_0, rt)
u_0 = rho_ul * u_ul / rho_0
ts_0 = t_0 + 10
q_0 = np.zeros(n_l) + 10000

stanje_0 = np.stack( (u_0, p_0, h_0, rho_0, t_0, ts_0, q_0) )

# rubni uvjeti

h_visokotlacno_min = h_ul
p_visokotlacno = p_ul + tlak_ekspanzije 
p_niskotlacno = p_iz
h_visokotlacno_max = h_ul


#h_visokotlacno = np.linspace(h_visokotlacno_max, h_visokotlacno_min, n_t)
h_visokotlacno = np.repeat(h_visokotlacno_min, n_t)
p_visokotlacno = np.repeat(p_visokotlacno, n_t)
p_niskotlacno = np.repeat(p_niskotlacno, n_t)
