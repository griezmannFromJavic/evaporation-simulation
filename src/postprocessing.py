import numpy as np
import CoolProp.CoolProp as CP
from funkcije import alpha_unutarnje, dpdl_trenja, gustoca_interpolirana, temperatura_interpolirana, x_termo_coolprop
from velicine import *
import os
import scipy.ndimage
import matplotlib.pyplot as plt

jedinice = ['[m / s]', '[Pa]', '[J / kg]', '[kg / m^3]', '[K]', '[K]', '[W / m^2]']

def prikaz_velicina_prostorno(trenutak, velicina):
    n = len(rjesenje[trenutak, velicina, :])
    plt.plot(np.linspace(0, l_c, n), rjesenje[trenutak, velicina, :], label=str(vremena[trenutak]))
    plt.legend()
    plt.xlabel('udaljenost od ustrujne povrsine [m]')
    plt.ylabel(raspored_velicina[velicina] + ' ' + jedinice[velicina])
    return None

def prikaz_velicina_u_svim_vremenima(velicina):
    for trenutak in range(len(vremena)):
        plt.plot(np.linspace(0, l_c, n_l), rjesenje[trenutak, velicina, :], label=str(vremena[trenutak]))
    plt.legend()
    plt.xlabel('udaljenost od ustrujne povrsine [m]')
    plt.ylabel(raspored_velicina[velicina] + ' ' + jedinice[velicina])
    return None

def prikaz_velicina_u_vremenskom_rasponu(velicina, pocetni_trenutak, krajnji_trenutak):
    for trenutak in range(pocetni_trenutak, krajnji_trenutak):
        plt.plot(np.linspace(0, l_c, n_l), rjesenje[trenutak, velicina, :], label=str(vremena[trenutak]))
    plt.xlabel('udaljenost od ustrujne povrsine [m]')
    plt.ylabel(raspored_velicina[velicina] + ' ' + jedinice[velicina])
    plt.legend()
    return None

def prikaz_izracunatih_velicina_prostorno(trenutak, velicina):
    plt.plot(np.linspace(0, l_c, n_l), rjesenje[trenutak, :], label=str(vremena[trenutak]))
    plt.legend()
    plt.xlabel('udaljenost od ustrujne povrsine [m]')
    plt.ylabel(raspored_velicina[velicina] + ' ' + jedinice[velicina])
    return None

def prikaz_velicina_vremenski(pocetni_trenutak, zavrsni_trenutak, velicina, lokacija):
    plt.plot(vremena[pocetni_trenutak : zavrsni_trenutak], rjesenje[pocetni_trenutak : zavrsni_trenutak, velicina, lokacija])
    plt.title(raspored_velicina[velicina] + ', ' + str(lokacije_x[lokacija]) + ' m od ulaza')
    plt.xlabel('vrijeme [s]')
    plt.ylabel(raspored_velicina[velicina] + ' ' + jedinice[velicina])
    return None

def citaj_stanja(ime_direktorija):
    putanja_rezultata = '/home/josip/Desktop/diplomski/kod/' + ime_direktorija + '/'
    broj_volumena = len( np.genfromtxt(putanja_rezultata + '0/brzina.csv', delimiter=',') )
    broj_rjesenih_trenutaka = len(os.listdir(putanja_rezultata))
    memorija = np.zeros((broj_rjesenih_trenutaka, broj_velicina, broj_volumena))
    for t in range(broj_rjesenih_trenutaka):
        for i in range(broj_velicina):
            putanja = putanja_rezultata + str(t) + '/' + raspored_velicina[i] + '.csv'
            memorija[t, i, :] = np.genfromtxt(putanja, delimiter=',')
    return memorija

direktorij_rezultata = 'rezultati'

rjesenje = citaj_stanja(direktorij_rezultata)

def citaj_info(file):
    putanja_rezultata = '/home/josip/Desktop/diplomski/kod/' + file +'/'
    broj_rjesenih_trenutaka = len( os.listdir(putanja_rezultata) )
    informacije = np.zeros((broj_rjesenih_trenutaka, 4))
    for t in range(broj_rjesenih_trenutaka):
        putanja = putanja_rezultata + str(t) + '/informacije.csv'
        informacije[t] = np.genfromtxt(putanja, delimiter=',')
    return informacije

informacije = citaj_info(direktorij_rezultata)
vremena = informacije[:, 0]
broj_iteracija = informacije[:, 1]
postignuti_reziduali = informacije[:, 2]
podrelaksacija = informacije[:, 3]

def udio_pare(rjesenje):
    tlak = rjesenje[:, 1, :]
    entalpija = rjesenje[:, 2, :]
    t = len(tlak)
    x = np.zeros(np.shape(tlak))
    for i in range(t):
        x[i] = CP.PropsSI('Q', 'P', tlak[i], 'H', entalpija[i], rt)
    return x

def termo_udio_pare():
    tlak = rjesenje[:, 1, :]
    entalpija = rjesenje[:, 2, :]
    t = len(tlak)
    h_fluid= np.zeros(np.shape(tlak))
    h_para = np.zeros(np.shape(tlak))
    x_termo = np.zeros(np.shape(tlak))
    for i in range(t):
        h_fluid[i] = CP.PropsSI('H', 'P', tlak[i], 'Q', 0, rt)
        h_para[i] = CP.PropsSI('H', 'P', tlak[i], 'Q', 1, rt)
        x_termo[i] = (entalpija[i] - h_fluid[i]) / (h_para[i] - h_fluid[i])
    return x_termo

def temperatura_isparavanja():
    tlak = rjesenje[:, 1, :]
    t = len(tlak)
    ts = np.zeros(np.shape(tlak))
    for i in range(t):
        ts[i] = CP.PropsSI('T', 'P', tlak[i], 'Q', 0, rt)
    return ts

def entalpija_isparavanja():
    tlak = rjesenje[:, 1, :]
    t = len(tlak)
    hs = np.zeros(np.shape(tlak))
    for i in range(t):
        hs[i] = CP.PropsSI('H', 'P', tlak[i], 'Q', 0, rt)
    return hs

def courant(trenutak):
    co = rjesenje[trenutak, 0, :] * (vremena[trenutak] - vremena[trenutak - 1]) / dx
    plt.plot(co)
    return co

def transformacija_rjesenja_na_drugu_mrezu(stari_broj_celija, novi_broj_celija, trenutak):
    info = np.array([0, 0, 0, 0])
    putanja_direktorija = '/home/josip/Desktop/diplomski/kod/rezultati_transformirani_na_' + str(novi_broj_celija) + '_celija/0/'
    os.makedirs(os.path.dirname(putanja_direktorija), exist_ok=True)
    putanja_info = putanja_direktorija + 'informacije.csv'
    np.savetxt(putanja_info, info, delimiter=",")
    for i in range(broj_velicina):
        novo_rjesenje = scipy.ndimage.zoom(rjesenje[trenutak, i, :], novi_broj_celija / stari_broj_celija)
        putanja_velicine = putanja_direktorija + raspored_velicina[i] + '.csv'
        np.savetxt(putanja_velicine, novo_rjesenje, delimiter=",")
    return print('snimljeno!')

def provjera_zoe_stacionarnog_stanja(ulazna_granica, izlazna_granica, trenutak):
    toplinski_tok_predan_fluidu = sum(rjesenje[trenutak, 6, ulazna_granica : izlazna_granica] * d_u * np.pi *
                                                    dx[ulazna_granica : izlazna_granica])
    toplinski_tok_predviden = (d_u ** 2 * np.pi / 4) *\
                              (rjesenje[trenutak, 3, izlazna_granica] * rjesenje[trenutak, 0, izlazna_granica] *
                               (rjesenje[trenutak, 2, izlazna_granica] + rjesenje[trenutak, 0, izlazna_granica] ** 2 / 2) -
                               rjesenje[trenutak, 3, ulazna_granica] * rjesenje[trenutak, 0, ulazna_granica] *
                               (rjesenje[trenutak, 2, ulazna_granica] +  rjesenje[trenutak, 0, ulazna_granica] ** 2 / 2))
    return toplinski_tok_predan_fluidu / toplinski_tok_predviden

def provjera_stacionarnosti_stijenke(trenutak):
    ts = rjesenje[trenutak, 5, :]
    qu = rjesenje[trenutak, 6, :]
    q_vanjsko = (temperatura_ogrijevnog_medija - ts) * alpha_vanjsko * d_v * np.pi * dx
    q_unutarnje = qu * d_u * np.pi * dx
    omjer = q_vanjsko / q_unutarnje
    plt.plot(omjer)
    plt.plot([0, n_l], [1, 1])
    return omjer

def provjera_alphe(p, x_min, x_max, m, q):
    h_fluid = CP.PropsSI('H', 'P', p, 'Q', 0, rt)
    h_para = CP.PropsSI('H', 'P', p, 'Q', 1, rt)
    h_isparavanja = h_para - h_fluid
    h_min = h_fluid + x_min * h_isparavanja
    h_max = h_fluid + x_max * h_isparavanja
    alpha = alpha_unutarnje(p, np.linspace(h_min, h_max, 1000), m, q)
    alpha_0 = alpha_unutarnje(p, h_fluid, m, q)
    omjer = alpha / alpha_0
    plt.plot(np.linspace(x_min, x_max, 1000), omjer, label=str(q))
    plt.legend()
    return None

def prikaz_temperatura(trenutak):
    ts = temperatura_isparavanja()
    plt.plot(np.linspace(0, l_c, n_l), ts[trenutak], label='temperatura isparavanja')
    # plt.plot([0, l_c], [temperatura_ogrijevnog_medija, temperatura_ogrijevnog_medija], label='temperatura ogrijevnog medija')
    plt.plot(np.linspace(0, l_c, n_l), rjesenje[trenutak, 4, :], label='temperatura fluida')
    plt.plot(np.linspace(0, l_c, n_l), rjesenje[trenutak, 5, :], label='temperatura stijenke')
    plt.legend()
    return None

def prikaz_x_termo(trenutak):
    xt = termo_udio_pare()[trenutak]
    plt.plot(np.linspace(0, l_c, n_l), xt)
    plt.plot([0, l_c], [0, 0])
    plt.xlabel('udaljenost od ustrujne površine [m]')
    plt.ylabel('termodinamički sadržaj pare [kg / kg]')
    return xt

def provjera_dpdl_trenja(p, x_min, x_max, m, epsilon):
    h_fluid = CP.PropsSI('H', 'P', p, 'Q', 0, rt)
    h_para = CP.PropsSI('H', 'P', p, 'Q', 1, rt)
    h_isparavanja = h_para - h_fluid
    h_min = h_fluid + x_min * h_isparavanja
    h_max = h_fluid + x_max * h_isparavanja
    dpdl = dpdl_trenja(p, np.linspace(h_min, h_max, 1000), m, epsilon)
    omjer = dpdl / dpdl_trenja(p, h_fluid, m, epsilon)

    plt.plot(np.linspace(x_min, x_max, 1000), omjer, label='cisto')
    plt.legend()
    return None

def provjera_aproksimiranih_gustoca(p, x_min, x_max):
    n = 200
    h_fluid = CP.PropsSI('H', 'P', p, 'Q', 0, rt)
    h_para = CP.PropsSI('H', 'P', p, 'Q', 1, rt)
    h_isparavanja = h_para - h_fluid
    h_min = h_fluid + x_min * h_isparavanja
    h_max = h_fluid + x_max * h_isparavanja
    entalpije = np.linspace(h_min, h_max, n)
    xt = np.linspace(x_min, x_max, n)
    gustoce_stvarno = CP.PropsSI('D', 'P', p, 'H', entalpije, rt)
    gustoce_aproksimacija = gustoca_interpolirana(np.repeat(p, n), entalpije)
    plt.plot(xt, gustoce_stvarno, label='gustoce iz baze podataka')
    plt.plot(xt, gustoce_aproksimacija, label='aproksimirane_gustoce')
    plt.legend()
    return None

def provjera_aproksimiranih_temperatura(p, x_min, x_max):
    n = 200
    h_fluid = CP.PropsSI('H', 'P', p, 'Q', 0, rt)
    h_para = CP.PropsSI('H', 'P', p, 'Q', 1, rt)
    h_isparavanja = h_para - h_fluid
    h_min = h_fluid + x_min * h_isparavanja
    h_max = h_fluid + x_max * h_isparavanja
    entalpije = np.linspace(h_min, h_max, n)
    xt = np.linspace(x_min, x_max, n)
    temperature_stvarno = CP.PropsSI('T', 'P', p, 'H', entalpije, rt)
    temperature_aproksimacija = temperatura_interpolirana(np.repeat(p, n), entalpije)
    plt.plot(xt, temperature_stvarno, label='gustoce iz baze podataka')
    plt.plot(xt, temperature_aproksimacija, label='aproksimirane_gustoce')
    plt.legend()
    return None

def prikaz_lokacija_isparavanja():
    xt = termo_udio_pare()
    umnozak_susjeda = xt[:, 1:] * xt[:, :-1]
    nultocke = np.ones(np.shape(umnozak_susjeda)) * (umnozak_susjeda < 0)
    elementi = np.transpose(np.nonzero(nultocke))[:,1]
    lokacije = [lokacije_x[i] for i in elementi]
    plt.plot(vremena, lokacije)
    return lokacije

def odnos_gustoca(tlak):
    return CP.PropsSI('D', 'P', tlak, 'Q', 0, rt) / CP.PropsSI('D', 'P', tlak, 'Q', 1, rt)

# nedovrseno
def achard_bezdimenzijske(trenutak):
    brzina = rjesenje[trenutak, 0, :]
    tlak = rjesenje[trenutak, 1, :]
    srednji_tlak = np.average(tlak)
    entalpija = rjesenje[trenutak, 2, :]
    gustoca = rjesenje[trenutak, 3, :]
    toplinski_tok = rjesenje[trenutak, 6, :]
    q_0 = np.average(toplinski_tok)

    p_h = d_u * np.pi

    v_f = 1 / CP.PropsSI('D', 'P', srednji_tlak, 'Q', 0, rt)
    v_g = 1 / CP.PropsSI('D', 'P', srednji_tlak, 'Q', 1, rt)
    v_fg = v_g - v_f
    h_f = CP.PropsSI('H', 'P', srednji_tlak, 'Q', 0, rt)
    h_g = CP.PropsSI('H', 'P', srednji_tlak, 'Q', 1, rt)
    h_fg = h_g - h_f
    h_i = rjesenje[trenutak, 2, 0]

    v_0 = (q_0 * p_h * v_f * l_c) / (a_u * (h_f - h_i))

    f = d_u / (2 * gustoca * brzina ** 2) * dpdl_trenja(tlak, entalpija, brzina * gustoca, epsilon)
    f = np.average(f)

    n_sub = (v_fg * (h_f - h_i)) / (v_f * h_fg)
    gamma = f * l_c / (2 * d_u)
    fr = v_0 ** 2 / (g * l_c)
    fr_inverz = 1 / fr
    j = brzina[0] / v_0 # provjeri

    return print('N_SUB = ', n_sub, '\nGAMMA = ', gamma, '\nFROUDE^-1 = ', fr_inverz, '\nj = ', j)