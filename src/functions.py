from velicine import *
import CoolProp.CoolProp as CP
from scipy.interpolate import CubicSpline, RectBivariateSpline
import os
import matplotlib.pyplot as plt


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

def x_termo(p, h):
    return (h - entalpija_f(p)) / (entalpija_g(p) - entalpija_f(p))

def x_termo_coolprop(p, h):
    h_f = CP.PropsSI('H', 'P', p, 'Q', 0, rt)
    h_g = CP.PropsSI('H', 'P', p, 'Q', 1, rt)
    return (h - h_f) / (h_g - h_f)

def x_realno(p, h):
    x = np.clip(x_termo(p, h), odmak_x, 1 - odmak_x)
    # x = np.clip(interpolirani_xt(p, h), odmak_x, 1 - odmak_x)
    return x

def brzina_ekspanzije(tlak_ekspanzije, gustoca):
    a = np.sign(tlak_ekspanzije)
    brzina = a * np.sqrt(2 * np.absolute(tlak_ekspanzije) / gustoca)
    # return np.sqrt(2 * tlak_ekspanzije / gustoca)
    return brzina

def dpdl_trenja(p, h, m, epsilon): # Müller-Steinhagen and Heck
    x = x_termo(p, h)
    x_r = x_realno(p, h)
    rho_l = gustoca_f(p)
    rho_g = gustoca_g(p)
    nu_l = viskoznost_f(p)
    nu_g = viskoznost_g(p)
    re_l = m * d_u / nu_l * (1 - x_r)
    re_g = m * d_u / nu_g * x_r
    A = (swamee_jain(re_l, epsilon) / d_u) * m ** 2 / (2 * rho_l) # wolwerine [13.2.43]
    B = (swamee_jain(re_g, epsilon) / d_u) * m ** 2 / (2 * rho_g) # wolwerine [13.2.43]
    C = (A + 2 * (B - A) * x_r) * (1 - x_r) ** (1 / 3) + B * x_r ** 3

    # alpha_sigma = 0.01
    alpha_sigma = 0.015
    dpdl_mali_x = (1 - sigma_func(x, 0 + odmak_x, alpha_sigma)) * A + sigma_func(x, 0 + odmak_x, alpha_sigma) * C
    dpdl_veliki_x = (1 - sigma_func(x, 1 - odmak_x, alpha_sigma)) * C + sigma_func(x, 1 - odmak_x, alpha_sigma) * B
    dpdl = (1 - sigma_func(x, 0.5, alpha_sigma)) * dpdl_mali_x + sigma_func(x, 0.5, alpha_sigma) * dpdl_veliki_x
    return dpdl

def petukhov(f, re, pr):
    return (f / 8) * re * pr / (1.07 + 12.7 * np.sqrt(f / 8) * (pr ** (2 / 3) - 1))

def alpha_unutarnje(p, h, m, q): # Gungor and Winterton (1987)
    x = x_termo(p, h)
    x_r = x_realno(p, h)
    k_l = toplinska_provodnost_f(p)
    k_g = toplinska_provodnost_g(p)
    pr_l = prandtl_f(p)
    pr_g = prandtl_g(p)
    mu_l = viskoznost_f(p)
    mu_g = viskoznost_g(p)
    rho_l = gustoca_f(p)
    rho_g = gustoca_g(p)
    h_lg = entalpija_g(p) - entalpija_f(p)
    bo = np.absolute( q / (m * h_lg) )
    e_new = 1 + 3000 * bo ** 0.86 + 1.12 * (x_r / (1 - x_r)) ** 0.75 * (rho_l / rho_g) ** 0.41
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
        putanja_direktorija = '/home/josip/Desktop/diplomski/kod/rezultati/' + str(trenutak) + '/'
        putanja_datoteke = putanja_direktorija + raspored_velicina[i] + '.csv'
        data = np.genfromtxt(putanja_datoteke, delimiter=',')
        stanje[i, :] = data
    return stanje

def generiraj_matrice_gustoce(p_min, p_max, broj_podjela_p, h_min, h_max, broj_podjela_h):
    tlakovi = np.linspace(p_min, p_max, broj_podjela_p)
    entalpije = np.linspace(h_min, h_max, broj_podjela_h)
    matrica = np.zeros((broj_podjela_p, broj_podjela_h))
    for i in range(broj_podjela_p):
        matrica[i] = stanje_gustoca(tlakovi[i], entalpije)
    return matrica

def generiraj_matrice_temperature(p_min, p_max, broj_podjela_p, h_min, h_max, broj_podjela_h):
    tlakovi = np.linspace(p_min, p_max, broj_podjela_p)
    entalpije = np.linspace(h_min, h_max, broj_podjela_h)
    matrica = np.zeros((broj_podjela_p, broj_podjela_h))
    for i in range(broj_podjela_p):
        matrica[i] = stanje_temperatura(tlakovi[i], entalpije)
    return matrica

def snimi_matrice(p_min, p_max, broj_podjela_p, h_min, h_max, broj_podjela_h):
    putanja_direktorija = '/home/josip/Desktop/diplomski/kod/velicine_stanja/'
    os.makedirs(os.path.dirname(putanja_direktorija), exist_ok=True)
    np.savetxt(putanja_direktorija + 'informacije.txt', [p_min, p_max, broj_podjela_p, h_min, h_max, broj_podjela_h],
               delimiter=",")

    tlakovi = np.linspace(p_min, p_max, broj_podjela_p)
    np.savetxt(putanja_direktorija + 'tlakovi.csv', tlakovi, delimiter=",")

    entalpije = np.linspace(h_min, h_max, broj_podjela_h)
    np.savetxt(putanja_direktorija + 'entalpije.csv', entalpije, delimiter=",")

    toplinska_provodnost_fluida = CP.PropsSI('CONDUCTIVITY', 'P', tlakovi, 'Q', 0, rt)
    np.savetxt(putanja_direktorija + 'toplinska_provodnost_fluida.csv', toplinska_provodnost_fluida, delimiter=",")

    toplinska_provodnost_plina = CP.PropsSI('CONDUCTIVITY', 'P', tlakovi, 'Q', 1, rt)
    np.savetxt(putanja_direktorija + 'toplinska_provodnost_plina.csv', toplinska_provodnost_plina, delimiter=",")

    prandtl_fluida = CP.PropsSI('PRANDTL', 'P', tlakovi, 'Q', 0, rt)
    np.savetxt(putanja_direktorija + 'prandtl_fluida.csv', prandtl_fluida, delimiter=",")

    prandtl_plina = CP.PropsSI('PRANDTL', 'P', tlakovi, 'Q', 1, rt)
    np.savetxt(putanja_direktorija + 'prandtl_plina.csv', prandtl_plina, delimiter=",")

    viskoznost_fluida = CP.PropsSI('VISCOSITY', 'P', tlakovi, 'Q', 0, rt)
    np.savetxt(putanja_direktorija + 'viskoznost_fluida.csv', viskoznost_fluida, delimiter=",")

    viskoznost_plina = CP.PropsSI('VISCOSITY', 'P', tlakovi, 'Q', 1, rt)
    np.savetxt(putanja_direktorija + 'viskoznost_plina.csv', viskoznost_plina, delimiter=",")

    gustoca_fluida = CP.PropsSI('D', 'P', tlakovi, 'Q', 0, rt)
    np.savetxt(putanja_direktorija + 'gustoca_fluida.csv', gustoca_fluida, delimiter=",")

    gustoca_plina = CP.PropsSI('D', 'P', tlakovi, 'Q', 1, rt)
    np.savetxt(putanja_direktorija + 'gustoca_plina.csv', gustoca_plina, delimiter=",")

    entalpija_fluida = CP.PropsSI('H', 'P', tlakovi, 'Q', 0, rt)
    np.savetxt(putanja_direktorija + 'entalpija_fluida.csv', entalpija_fluida, delimiter=",")

    entalpija_plina = CP.PropsSI('H', 'P', tlakovi, 'Q', 1, rt)
    np.savetxt(putanja_direktorija + 'entalpija_plina.csv', entalpija_plina, delimiter=",")

    gustoce = generiraj_matrice_gustoce(p_min, p_max, broj_podjela_p, h_min, h_max, broj_podjela_h)
    np.savetxt(putanja_direktorija + 'gustoce.csv', gustoce, delimiter=",")

    temperature = generiraj_matrice_temperature(p_min, p_max, broj_podjela_p, h_min, h_max, broj_podjela_h)
    np.savetxt(putanja_direktorija + 'temperature.csv', temperature, delimiter=",")
    return 'snimljeno!'

putanja_stanja = '/home/josip/Desktop/diplomski/kod/velicine_stanja/'

tlakovi = np.genfromtxt(putanja_stanja + 'tlakovi.csv', delimiter=',')
entalpije = np.genfromtxt(putanja_stanja + 'entalpije.csv', delimiter=',')

toplinska_provodnost_fluida = np.genfromtxt(putanja_stanja + 'toplinska_provodnost_fluida.csv', delimiter=',')
toplinska_provodnost_plina = np.genfromtxt(putanja_stanja + 'toplinska_provodnost_plina.csv', delimiter=',')
prandtl_fluida = np.genfromtxt(putanja_stanja + 'prandtl_fluida.csv', delimiter=',')
prandtl_plina = np.genfromtxt(putanja_stanja + 'prandtl_plina.csv', delimiter=',')
viskoznost_fluida = np.genfromtxt(putanja_stanja + 'viskoznost_fluida.csv', delimiter=',')
viskoznost_plina = np.genfromtxt(putanja_stanja + 'viskoznost_plina.csv', delimiter=',')
gustoca_fluida = np.genfromtxt(putanja_stanja + 'gustoca_fluida.csv', delimiter=',')
gustoca_plina = np.genfromtxt(putanja_stanja + 'gustoca_plina.csv', delimiter=',')
entalpija_fluida = np.genfromtxt(putanja_stanja + 'entalpija_fluida.csv', delimiter=',')
entalpija_plina = np.genfromtxt(putanja_stanja + 'entalpija_plina.csv', delimiter=',')
matrica_gustoce = np.genfromtxt(putanja_stanja + 'gustoce.csv', delimiter=',')
matrica_temperature = np.genfromtxt(putanja_stanja + 'temperature.csv', delimiter=',')



def toplinska_provodnost_f(p):
    return CubicSpline(tlakovi, toplinska_provodnost_fluida).__call__(p)

def toplinska_provodnost_g(p):
    return CubicSpline(tlakovi, toplinska_provodnost_plina).__call__(p)

def prandtl_f(p):
    return CubicSpline(tlakovi, prandtl_fluida).__call__(p)

def prandtl_g(p):
    return CubicSpline(tlakovi, prandtl_plina).__call__(p)

def viskoznost_f(p):
    return CubicSpline(tlakovi, viskoznost_fluida).__call__(p)

def viskoznost_g(p):
    return CubicSpline(tlakovi, viskoznost_plina).__call__(p)

def gustoca_f(p):
    return CubicSpline(tlakovi, gustoca_fluida).__call__(p)

def gustoca_g(p):
    return CubicSpline(tlakovi, gustoca_plina).__call__(p)

def entalpija_f(p):
    return CubicSpline(tlakovi, entalpija_fluida).__call__(p)

def entalpija_g(p):
    return CubicSpline(tlakovi, entalpija_plina).__call__(p)

def gustoca_interpolirana(p, h):
    return RectBivariateSpline(tlakovi, entalpije, matrica_gustoce).ev(p, h)

def temperatura_interpolirana(p, h):
    return RectBivariateSpline(tlakovi, entalpije, matrica_temperature).ev(p, h)

# im = plt.imshow(matrica_gustoce, cmap='hot')
# plt.colorbar(im, orientation='horizontal')

# 15e5_30e5_4e5_3.5e6_400