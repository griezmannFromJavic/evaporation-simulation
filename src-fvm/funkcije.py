from velicine import *
import CoolProp.CoolProp as CP
import os

odmak_x = 0.001

def manhattan_vektora(vektor_1, vektor_2):
    return np.average(np.absolute((vektor_2 - vektor_1) / vektor_1))

def manhattan_rjesenja(rjesenje_1, rjesenje_2, smanjenje):
    uzorak_1 = rjesenje_1[:, ::smanjenje]
    uzorak_2 = rjesenje_2[:, ::smanjenje] # provjerava se svaka n-ta ćelija n=smanjenje
    n = len(uzorak_1)
    rezidual = np.zeros(n)
    for i in range(n):
        rezidual[i] = manhattan_vektora(uzorak_1[i], uzorak_2[i])
    return rezidual

def eksplicitna(izracunata_velicina, pretpostavka, podrelaksacija):
    return podrelaksacija * izracunata_velicina + (1 - podrelaksacija) * pretpostavka

def zom_brzina(gustoca, gustoca_proslog_trenutka, brzina_na_ulazu, gustoca_na_ulazu, dx, pretpostavka_brzine,
               podrelaksacija):
    n = len(dx)
    a = np.zeros(n)
    a[1:-1] = - gustoca[1:-1]
    a[0] = 0
    a[-1] = - gustoca[-1]
    a = a[1:]
    b = np.zeros(n)
    b[1:-1] = gustoca[2:]
    b[0] = gustoca[1]
    b[-1] = gustoca[-1]
    b = b / podrelaksacija
    c = np.zeros(n)
    c[1:-1] = 0
    c[0] = 0
    c[-1] = 0
    c = c[:-1]
    d = np.zeros(n)
    d[1:-1] = - (gustoca[1:-1] - gustoca_proslog_trenutka[1:-1]) * dx[1:-1] / dt
    d[0] = - (gustoca[0] - gustoca_proslog_trenutka[0]) * dx[0] / dt + brzina_na_ulazu * gustoca_na_ulazu
    d[-1] = - (gustoca[-1] - gustoca_proslog_trenutka[-1]) * dx[-1] / dt
    d = d + (1 - podrelaksacija) * b * pretpostavka_brzine # b = b / podrelaksacija
    return tdma(a, b, c, d)

def zokg_brzina(pretpostavka_brzine, brzina_na_ulazu, tlak, tlak_na_izlazu, entalpija, gustoca, gustoca_na_ulazu, dx,
                podrelaksacija):
    n = len(dx)
    a = np.zeros(n)
    a[1:-1] = - gustoca[1:-1] * (dx[1:-1] / dt + pretpostavka_brzine[:-2])
    a[0] = 0
    a[-1] = - gustoca[-1] * (dx[-1] / dt + pretpostavka_brzine[-2])
    a = a[1:]
    b = np.zeros(n)
    b[1:-1] = gustoca[2:] * (dx[1:-1] / dt + pretpostavka_brzine[1:-1])
    b[0] = gustoca[1] * (dx[0] / dt + pretpostavka_brzine[0])
    b[-1] = gustoca[-1] * (dx[-1] / dt + pretpostavka_brzine[-1])
    b = b / podrelaksacija
    c = np.zeros(n - 1)
    d = np.zeros(n)
    d[1:-1] = - dpdl_trenja(tlak[1:-1], entalpija[1:-1], pretpostavka_brzine[1:-1] * gustoca[1:-1], epsilon) * dx[1:-1] \
              - gustoca[1:-1] * g * sin_beta * dx[1:-1] \
              - (tlak[2:] - tlak[1:-1])
    d[0] = gustoca_na_ulazu * brzina_na_ulazu * (dx[0] / dt  + brzina_na_ulazu) \
           - dpdl_trenja(tlak[0], entalpija[0], pretpostavka_brzine[0] * gustoca[0], epsilon) * dx[0] \
           - gustoca[0] * g * sin_beta * dx[0] \
           - (tlak[1] - tlak[0])
    d[-1] = - dpdl_trenja(tlak[-1], entalpija[-1], pretpostavka_brzine[-1] * gustoca[-1], epsilon) * dx[-1] \
            - gustoca[-1] * g * sin_beta * dx[-1] \
            - (tlak_na_izlazu - tlak[-1])
    d = d + (1 - podrelaksacija) * b * pretpostavka_brzine
    return tdma(a, b, c, d)

def zokg_tlak(brzina, brzina_proslog_trenutka, brzina_na_ulazu, gustoca, gustoca_proslog_trenutka, gustoca_na_ulazu,
              tlak_na_izlazu, entalpija, dx, pretpostavka_tlaka, podrelaksacija):
    n = len(dx)
    a = np.zeros(n)
    a[1:-1] = 0
    a[0] = 0
    a[-1] = 0
    a = a[1:]
    b = np.zeros(n)
    b[1:-1] = 1
    b[0] = 1
    b[-1] = 1
    b = b / podrelaksacija
    c = np.zeros(n)
    c[1:-1] = - 1
    c[0] = - 1
    c[-1] = 0
    c = c[:-1]
    d = np.zeros(n)
    d[1:-1] = (gustoca[1:-1] * brzina[1:-1] - gustoca_proslog_trenutka[1:-1] * brzina_proslog_trenutka[1:-1]) * dx[1:-1] / dt \
              + gustoca[1:-1] * brzina[:-2] ** 2 \
              - gustoca[2:] * brzina[1:-1] ** 2 \
              + dpdl_trenja(pretpostavka_tlaka[1:-1], entalpija[1:-1], brzina[1:-1] * gustoca[1:-1], epsilon) * dx[1:-1] \
              + gustoca[1:-1] * g * sin_beta * dx[1:-1]
    d[0] = + (gustoca[0] * brzina[0] - gustoca_proslog_trenutka[0] * brzina_proslog_trenutka[0]) * dx[0] / dt \
           + gustoca_na_ulazu * brzina_na_ulazu ** 2 \
           - gustoca[1] * brzina[0] ** 2 \
           + dpdl_trenja(pretpostavka_tlaka[0], entalpija[0], brzina[0] * gustoca[0], epsilon) * dx[0] \
           + gustoca[0] * g * sin_beta * dx[0]
    d[-1] = (gustoca[-1] * brzina[-1] - gustoca_proslog_trenutka[-1] * brzina_proslog_trenutka[-1]) * dx[-1] / dt \
            + gustoca[-1] * brzina[-2] ** 2 \
            - gustoca[-1] * brzina[-1] ** 2 \
            + tlak_na_izlazu \
            + dpdl_trenja(pretpostavka_tlaka[0], entalpija[0], brzina[0] * gustoca[0], epsilon) * dx[-1] \
            + gustoca[-1] * g * sin_beta * dx[-1]
    d = d + (1 - podrelaksacija) * b * pretpostavka_tlaka
    return tdma(a, b, c, d)

def zom_zokg_tlak(brzina, pretpostavka_tlaka, tlak_visokotlacno, tlak_niskotlacno, entalpija, gustoca,
                  gustoca_proslog_trenutka, dx):
    n = len(dx)
    a = np.zeros(n)
    a[1:-1] = - gustoca[:-2] / (gustoca[1:-1] * (dx[1:-1] / dt + brzina[:-2]))
    a[0] = 0
    a[-1] = - gustoca[-2] / (gustoca[-1] * (dx[-1] / dt + brzina[-2]))
    a = a[1:]
    b = np.zeros(n)
    b[1:-1] = gustoca[:-2] / (gustoca[1:-1] * (dx[1:-1] / dt + brzina[:-2])) \
              + gustoca[1:-1] / (gustoca[2:] * (dx[1:-1] / dt + brzina[1:-1]))
    b[0] = gustoca[0] / (gustoca[1] * (dx[0] / dt + brzina[0])) \
           + gustoca[0] / (np.sqrt(2 * gustoca[0] * (tlak_visokotlacno - pretpostavka_tlaka[0])))
    b[-1] = gustoca[-2] / (gustoca[-1] * (dx[-1] / dt + brzina[-2])) \
            + gustoca[-1] / (gustoca[-1] * (dx[-1] / dt + brzina[-1]))
    c = np.zeros(n)
    c[1:-1] = - gustoca[1:-1] / (gustoca[2:] * (dx[1:-1] / dt + brzina[1:-1]))
    c[0] = - gustoca[0] / (gustoca[1] * (dx[0] / dt + brzina[0]))
    c[-1] = 0
    c = c[:-1]
    brzina_i_minus_dva = np.append(brzina[0], brzina)[:-3]
    d = np.zeros(n)
    d[1:-1] = - dx[1:-1] / dt * (gustoca[1:-1] - gustoca_proslog_trenutka[1:-1]) \
              + (gustoca[:-2] * brzina_i_minus_dva * (dx[1:-1] / dt + brzina_i_minus_dva)
                 - dpdl_trenja(pretpostavka_tlaka[:-2], entalpija[:-2], brzina[:-2] * gustoca[:-2], epsilon) * dx[1:-1]
                 - gustoca[:-2] * g * sin_beta * dx[1:-1]) \
              / (gustoca[1:-1] * (dx[1:-1] / dt + brzina[:-2])) * gustoca[:-2] \
              - (gustoca[1:-1] * brzina[:-2] * (dx[1:-1] / dt + brzina[:-2])
                 - dpdl_trenja(pretpostavka_tlaka[1:-1], entalpija[1:-1], brzina[1:-1] * gustoca[1:-1], epsilon) * dx[1:-1]
                 - gustoca[1:-1] * g * sin_beta * dx[1:-1]) \
              / (gustoca[2:] * (dx[1:-1] / dt + brzina[1:-1])) * gustoca[1:-1]
    d[0] = dx[0] / dt * (gustoca[0] - gustoca_proslog_trenutka[0]) \
           + np.sqrt(2 * gustoca[0] * (tlak_visokotlacno - pretpostavka_tlaka[0])) \
           + pretpostavka_tlaka[0] * np.sqrt(gustoca[0]) / np.sqrt(2 * (tlak_visokotlacno - pretpostavka_tlaka[0])) \
           - (gustoca[0] * brzina[0] * (dx[0] / dt + brzina[0])
              - dpdl_trenja(pretpostavka_tlaka[0], entalpija[0], brzina[0] * gustoca[0], epsilon) * dx[0]
              - gustoca[0] * g * sin_beta * dx[0]) \
           / (gustoca[1] * (dx[0] / dt + brzina[0])) * gustoca[0]
    d[-1] = gustoca[-1] / (gustoca[-1] * (dx[-1] / dt + brzina[-1])) * tlak_niskotlacno \
            - dx[-1] / dt * (gustoca[-1] - gustoca_proslog_trenutka[-1]) \
            + (gustoca[-2] * brzina[-3] * (dx[-1] / dt + brzina[-3])
               - dpdl_trenja(pretpostavka_tlaka[-2], entalpija[-2], brzina[-2] * gustoca[-2], epsilon) * dx[-1]
               - gustoca[-2] * g * sin_beta * dx[-1]) \
            / (gustoca[-1] * (dx[-1] / dt + brzina[-2])) * gustoca[-2] \
            - (gustoca[-1] * brzina[-2] * (dx[-1] / dt + brzina[-2])
               - dpdl_trenja(pretpostavka_tlaka[-1], entalpija[-1], brzina[-1] * gustoca[-1], epsilon) * dx[-1]
               - gustoca[-1] * g * sin_beta * dx[-1]) \
            / (gustoca[-1] * (dx[-1] / dt + brzina[-1])) * gustoca[-1]
    return tdma(a, b, c, d)

def zoe_entalpija(brzina, brzina_proslog_trenutka, brzina_na_ulazu, gustoca, gustoca_proslog_trenutka, gustoca_na_ulazu,
                  tlak, tlak_proslog_trenutka, entalpija_proslog_trenutka, entalpija_na_ulazu, q_w, dx,
                  pretpostavka_entalpije, podrelaksacija):
    n = len(dx)
    a = np.zeros(n)
    a[1:-1] = - gustoca[:-2] * brzina[:-2]
    a[0] = 0
    a[-1] = - gustoca[-2] * brzina[-2]
    a = a[1:]
    b = np.zeros(n)
    b[1:-1] = gustoca[1:-1] * dx[1:-1] / dt + gustoca[1:-1] * brzina[1:-1]
    b[0] = gustoca[0] * dx[0] / dt + gustoca[0] * brzina[0]
    b[-1] = gustoca[-1] * dx[-1] / dt + gustoca[-1] * brzina[-1]
    b = b / podrelaksacija
    c = np.zeros(n - 1)
    d = np.zeros(n)
    d[1:-1] = + (tlak[1:-1] - tlak_proslog_trenutka[1:-1]) * dx[1:-1] / dt \
              + q_w[1:-1] * dx[1:-1] * d_u * np.pi / a_u \
              - gustoca[1:-1] * brzina[1:-1] ** 2 / 2 * dx[1:-1] / dt \
              + gustoca_proslog_trenutka[1:-1] * (entalpija_proslog_trenutka[1:-1] + brzina_proslog_trenutka[1:-1] ** 2 / 2) * dx[1:-1] / dt \
              - gustoca[1:-1] * brzina[1:-1] ** 3 / 2 \
              + gustoca[:-2] * brzina[:-2] ** 3 / 2 \
              - gustoca[1:-1] * brzina[1:-1] * g * sin_beta * dx[1:-1]
    d[0] = + (tlak[0] - tlak_proslog_trenutka[0]) * dx[0] / dt \
           + q_w[0] * dx[0] * d_u * np.pi / a_u \
           - gustoca[0] * brzina[0] ** 2 / 2 * dx[0] / dt \
           + gustoca_proslog_trenutka[0] * (entalpija_proslog_trenutka[0] + brzina_proslog_trenutka[0] ** 2 / 2) * dx[0] / dt \
           - gustoca[0] * brzina[0] ** 3 / 2 \
           + gustoca_na_ulazu * brzina_na_ulazu ** 3 / 2 \
           + gustoca_na_ulazu * brzina_na_ulazu * entalpija_na_ulazu \
           - gustoca[0] * brzina[0] * g * sin_beta * dx[0]
    d[-1] = + (tlak[-1] - tlak_proslog_trenutka[-1]) * dx[-1] / dt \
            + q_w[-1] * dx[-1] * d_u * np.pi / a_u \
            - gustoca[-1] * brzina[-1] ** 2 / 2 * dx[-1] / dt \
            + gustoca_proslog_trenutka[-1] * (entalpija_proslog_trenutka[-1] + brzina_proslog_trenutka[-1] ** 2 / 2) * dx[-1] / dt \
            - gustoca[-1] * brzina[-1] ** 3 / 2 \
            + gustoca[-2] * brzina[-2] ** 3 / 2 \
            - gustoca[-1] * brzina[-1] * g * sin_beta * dx[-1]
    d = d + (1 - podrelaksacija) * b * pretpostavka_entalpije
    return tdma(a, b, c, d)

def temp_stijenke(temp_stijenke_proslog_trenutka, temperatura_fluida, alpha_u, dx):
    n = len(dx)
    a = np.zeros(n)
    a[1:-1] = - lambda_c / (rho_c * c_c * dx[1:-1])
    a[0] = 0
    a[-1] = lambda_c / (rho_c * c_c * dx[-1])
    a = a[1:]
    b = np.zeros(n)
    b[1:-1] = dx[1:-1] / dt \
              + 2 * lambda_c / (rho_c * c_c * dx[1:-1]) \
              + 4 * d_v * alpha_vanjsko / (rho_c * c_c * (d_v ** 2 - d_u ** 2)) \
              + 4 * d_u * alpha_u[1:-1] / (rho_c * c_c * (d_v ** 2 - d_u ** 2))
    b[0] = dx[0] / dt \
           + lambda_c / (rho_c * c_c * dx[0]) \
           + 4 * d_v * alpha_vanjsko / (rho_c * c_c * (d_v ** 2 - d_u ** 2)) \
           + 4 * d_u * alpha_u[0] / (rho_c * c_c * (d_v ** 2 - d_u ** 2))
    b[-1] = dx[-1] / dt \
            - lambda_c / (rho_c * c_c * dx[-1]) \
            + 4 * d_v * alpha_vanjsko / (rho_c * c_c * (d_v ** 2 - d_u ** 2)) \
            + 4 * d_u * alpha_u[-1] / (rho_c * c_c * (d_v ** 2 - d_u ** 2))
    c = np.zeros(n)
    c[1:-1] = - lambda_c / (rho_c * c_c * dx[1:-1])
    c[0] = - lambda_c / (rho_c * c_c * dx[0])
    c[-1] = 0
    c = c[:-1]
    d = dx / dt * temp_stijenke_proslog_trenutka \
        + 4 * d_v * alpha_vanjsko / (rho_c * c_c * (d_v ** 2 - d_u ** 2)) * temperatura_ogrijevnog_medija \
        + 4 * d_u * alpha_u / (rho_c * c_c * (d_v ** 2 - d_u ** 2)) * temperatura_fluida
    return tdma(a, b, c, d)

def stanje_gustoca(tlak, entalpija):
    return CP.PropsSI('D', 'P', tlak, 'H', entalpija, rt)

def stanje_temperatura(tlak, entalpija):
    return CP.PropsSI('T', 'P', tlak, 'H', entalpija, rt)

def tdma(a, b, c, d):
    n = len(d)
    w = np.zeros(n - 1, float)
    g = np.zeros(n, float)
    p = np.zeros(n, float)

    w[0] = c[0] / b[0]
    g[0] = d[0] / b[0]

    for i in range(1, n - 1):
        w[i] = c[i] / (b[i] - a[i - 1] * w[i - 1])
    for i in range(1, n):
        g[i] = (d[i] - a[i - 1] * g[i - 1]) / (b[i] - a[i - 1] * w[i - 1])
    p[n - 1] = g[n - 1]
    for i in range(n - 1, 0, -1):
        p[i - 1] = g[i - 1] - w[i - 1] * p[i]
    return p

def sigma_func(x, x0, alpha):
    return 1 / (1 + np.exp(-(x - x0) / alpha))

def swamee_jain(re, epsilon):
    f = 0.25 / (np.log10((epsilon / 3.7) + (5.74 / (np.absolute(re) ** 0.9)))) ** 2
    return f

#slican swamee jain, ne koristi se
def moody(re, epsilon):
    f = 0.0055 * (1 + (2e4 * epsilon + 1e6 / re) ** (1 / 3))
    return f

def x_termo(p, h):
    h_fluid = CP.PropsSI('H', 'P', p, 'Q', 0, rt)
    h_para = CP.PropsSI('H', 'P', p, 'Q', 1, rt)
    x_t = (h - h_fluid) / (h_para - h_fluid)
    return x_t

def x_realno(p, h):
    x = np.clip(x_termo(p, h), odmak_x, 1 - odmak_x)
    return x

def brzina_ekspanzije(tlak_ekspanzije, gustoca):
    return np.sqrt(2 * tlak_ekspanzije / gustoca)

def dpdl_bankoff(p, h, m, epsilon):
    x = x_termo(p, h)
    x_r = x_realno(p, h)
    rho_l = CP.PropsSI('D', 'P', p, 'Q', 0, rt)
    rho_g = CP.PropsSI('D', 'P', p, 'Q', 1, rt)
    nu_l = CP.PropsSI('VISCOSITY', 'P', p, 'Q', 0, rt)
    nu_g = CP.PropsSI('VISCOSITY', 'P', p, 'Q', 1, rt)
    re_l = m * d_u / nu_l * (1 - x_r)
    re_g = m * d_u / nu_g * x_r
    A = (swamee_jain(re_l, epsilon) / d_u) * m ** 2 / (2 * rho_l) # wolwerine [13.2.43]
    B = (swamee_jain(re_g, epsilon) / d_u) * m ** 2 / (2 * rho_g)

    omjer_gustoca = rho_g / rho_l
    gamma = (0.71 + 2.35 * omjer_gustoca) / (1 + (1 - x_r) / x_r * omjer_gustoca)
    phi_bf = 1 / (1 - x_r) * (1 - gamma * (1 - omjer_gustoca)) ** (3 / 7) * (1 + x_r * (1 / omjer_gustoca - 1))
    C = A * phi_bf ** (7 / 4)

    alpha_sigma = 0.01
    dpdl_mali_x = (1 - sigma_func(x, 0 + odmak_x, alpha_sigma)) * A + sigma_func(x, 0 + odmak_x, alpha_sigma) * C
    dpdl_veliki_x = (1 - sigma_func(x, 1 - odmak_x, alpha_sigma)) * C + sigma_func(x, 1 - odmak_x, alpha_sigma) * B
    dpdl = (1 - sigma_func(x, 0.5, alpha_sigma)) * dpdl_mali_x + sigma_func(x, 0.5, alpha_sigma) * dpdl_veliki_x
    return dpdl

def dpdl_trenja(p, h, m, epsilon): # Müller-Steinhagen and Heck
    x = x_termo(p, h)
    x_r = x_realno(p, h)
    rho_l = CP.PropsSI('D', 'P', p, 'Q', 0, rt)
    rho_g = CP.PropsSI('D', 'P', p, 'Q', 1, rt)
    nu_l = CP.PropsSI('VISCOSITY', 'P', p, 'Q', 0, rt)
    nu_g = CP.PropsSI('VISCOSITY', 'P', p, 'Q', 1, rt)
    re_l = m * d_u / nu_l * (1 - x_r)
    re_g = m * d_u / nu_g * x_r
    A = (swamee_jain(re_l, epsilon) / d_u) * m ** 2 / (2 * rho_l) # wolwerine [13.2.43]
    B = (swamee_jain(re_g, epsilon) / d_u) * m ** 2 / (2 * rho_g) # wolwerine [13.2.43]
    C = (A + 2 * (B - A) * x_r) * (1 - x_r) ** (1 / 3) + B * x_r ** 3

    alpha_sigma = 0.01
    dpdl_mali_x = (1 - sigma_func(x, 0 + odmak_x, alpha_sigma)) * A + sigma_func(x, 0 + odmak_x, alpha_sigma) * C
    dpdl_veliki_x = (1 - sigma_func(x, 1 - odmak_x, alpha_sigma)) * C + sigma_func(x, 1 - odmak_x, alpha_sigma) * B
    dpdl = (1 - sigma_func(x, 0.5, alpha_sigma)) * dpdl_mali_x + sigma_func(x, 0.5, alpha_sigma) * dpdl_veliki_x
    return dpdl

def dpdl_trenja_regularizirano(p, h, m, epsilon):
    dh = 3e4
    dpdl_reg = (dpdl_trenja(p, h - dh / 2, m, epsilon) + dpdl_trenja(p, h + dh / 2, m, epsilon)) / 2
    return dpdl_reg

def petukhov(f, re, pr):
    return (f / 8) * re * pr / (1.07 + 12.7 * np.sqrt(f / 8) * (pr ** (2 / 3) - 1))

def alpha_unutarnje(p, h, m, q): # Gungor and Winterton (1987)
    x = x_termo(p, h)
    x_r = x_realno(p, h)
    k_l = CP.PropsSI('CONDUCTIVITY', 'P', p, 'Q', 0, rt)
    k_g = CP.PropsSI('CONDUCTIVITY', 'P', p, 'Q', 1, rt)
    pr_l = CP.PropsSI('PRANDTL', 'P', p, 'Q', 0, rt)
    pr_g = CP.PropsSI('PRANDTL', 'P', p, 'Q', 1, rt)
    mu_l = CP.PropsSI('VISCOSITY', 'P', p, 'Q', 0, rt)
    mu_g = CP.PropsSI('VISCOSITY', 'P', p, 'Q', 1, rt)
    rho_l = CP.PropsSI('D', 'P', p, 'Q', 0, rt)
    rho_g = CP.PropsSI('D', 'P', p, 'Q', 1, rt)
    h_lg = CP.PropsSI('H', 'P', p, 'Q', 1, rt) - CP.PropsSI('H', 'P', p, 'Q', 0, rt)
    bo = np.absolute( q / (m * h_lg) )
    e_new = 1 + 3000 * bo ** 0.86 + 1.12 * (x_r / (1 - x_r)) ** 0.75 * (rho_l / rho_g) ** 0.41
    # e_new = 0.3 * e_new #### KOREKCIJA ZA STABILNOST
    re_l = (1 - x_r) * m * d_u / mu_l
    re_g = x_r * m * d_u / mu_g
    f_l = swamee_jain(re_l, epsilon)
    f_g = swamee_jain(re_g, epsilon)

    alpha_l = petukhov(f_l, re_l, pr_l) * (k_l / d_u)
    alpha_g = petukhov(f_g, re_g, pr_g) * (k_g / d_u)
    alpha_tp = e_new * alpha_l

    alpha_sigma = 0.01

    alpha_malo_x = (1 - sigma_func(x, 0 + odmak_x, alpha_sigma)) * alpha_l + sigma_func(x, 0 + odmak_x, alpha_sigma) * alpha_tp
    alpha_veliko_x = (1 - sigma_func(x, 1 - odmak_x, alpha_sigma)) * alpha_tp + sigma_func(x, 1 - odmak_x, alpha_sigma) * alpha_g
    alpha = (1 - sigma_func(x, 0.5, alpha_sigma)) * alpha_malo_x + sigma_func(x, 0.5, alpha_sigma) * alpha_veliko_x
    return alpha

def spec_toplinski_tok(alpha, temperatura_visa, temperatura_niza):
    q = alpha * (temperatura_visa - temperatura_niza)
    return q

def citaj_stanje(trenutak):
    stanje = np.zeros( (broj_velicina, n_l) )
    for i in range( broj_velicina ):
        putanja_direktorija = os.path.join('..', 'results/', str(trenutak))
        putanja_datoteke = putanja_direktorija + '/' + raspored_velicina[i] + '.csv'
        data = np.genfromtxt(putanja_datoteke, delimiter=',')
        stanje[i, :] = data
    return stanje
