from funkcije import brzina_ekspanzije, zokg_brzina, zom_zokg_tlak, zoe_entalpija, eksplicitna, alpha_unutarnje, \
    temp_stijenke, manhattan_rjesenja, citaj_stanje
from velicine import *
from uvjeti import h_visokotlacno, p_visokotlacno, p_niskotlacno, stanje_0
import os
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP


def vanjska_iteracija_velocity_based(trenutak, pretpostavka, vrijednost_proslog_trenutka, podrelaksacija):
    tlak_na_ulazu = pretpostavka[1, 0]
    brzina_na_ulazu = brzina_ekspanzije(p_visokotlacno[trenutak] - tlak_na_ulazu, pretpostavka[3,0])
    entalpija_na_ulazu = h_visokotlacno[trenutak] - brzina_na_ulazu ** 2 / 2
    gustoca_na_ulazu = CP.PropsSI('D', 'P', tlak_na_ulazu, 'H', entalpija_na_ulazu, rt)

    brzina = zokg_brzina(pretpostavka[0], brzina_na_ulazu, pretpostavka[1], p_niskotlacno[trenutak], pretpostavka[2],
        pretpostavka[3], gustoca_na_ulazu, dx, podrelaksacija[0])

    tlak = zom_zokg_tlak(brzina, pretpostavka[1], p_visokotlacno[trenutak], p_niskotlacno[trenutak], pretpostavka[2],
        pretpostavka[3], vrijednost_proslog_trenutka[3], dx)
    tlak = eksplicitna(tlak, pretpostavka[1], podrelaksacija[1])

    entalpija = zoe_entalpija(brzina, vrijednost_proslog_trenutka[0], brzina_na_ulazu, pretpostavka[3],
        vrijednost_proslog_trenutka[3], gustoca_na_ulazu, tlak, vrijednost_proslog_trenutka[1],
        vrijednost_proslog_trenutka[2], entalpija_na_ulazu, pretpostavka[6], dx, pretpostavka[2],
        podrelaksacija[2])

    gustoca = gustoca_interpolirana(tlak, entalpija)
    gustoca = eksplicitna(gustoca, pretpostavka[3], podrelaksacija[3])

    temperatura = temperatura_interpolirana(tlak, entalpija)
    temperatura = eksplicitna(temperatura, pretpostavka[4], podrelaksacija[4])

    alpha_unutarnje_izracunato = alpha_unutarnje(tlak, entalpija, brzina * gustoca, pretpostavka[6]) # ne zapisuje se

    temperatura_stijenke = temp_stijenke(vrijednost_proslog_trenutka[5], temperatura, alpha_unutarnje_izracunato, dx)
    temperatura_stijenke = eksplicitna(temperatura_stijenke, pretpostavka[5], podrelaksacija[5])

    specificni_toplinski_tok = alpha_unutarnje_izracunato * (temperatura_stijenke - temperatura)
    specificni_toplinski_tok = eksplicitna(specificni_toplinski_tok, pretpostavka[6], podrelaksacija[6])

    rjesenje = np.stack((brzina, tlak, entalpija, gustoca, temperatura, temperatura_stijenke, specificni_toplinski_tok))
    return rjesenje

# krivo, a dobro radi, mozda cak stabilniji od običnog koda
def vanjska_iteracija(trenutak, pretpostavka, vrijednost_proslog_trenutka, podrelaksacija):
    tlak = zom_zokg_tlak(pretpostavka[0], pretpostavka[1], p_visokotlacno[trenutak], p_niskotlacno[trenutak], pretpostavka[2],
        pretpostavka[3], vrijednost_proslog_trenutka[3], dx)
    tlak = eksplicitna(tlak, pretpostavka[1], podrelaksacija[1])

    tlak_na_ulazu = tlak[0]

    brzina_na_ulazu = brzina_ekspanzije(p_visokotlacno[trenutak] - tlak_na_ulazu, pretpostavka[3,0])
    entalpija_na_ulazu = h_visokotlacno[trenutak] - brzina_na_ulazu ** 2 / 2
    gustoca_na_ulazu = CP.PropsSI('D', 'P', tlak_na_ulazu, 'H', entalpija_na_ulazu, rt)

    brzina = zokg_brzina(pretpostavka[0], brzina_na_ulazu, tlak, p_niskotlacno[trenutak], pretpostavka[2],
        pretpostavka[3], gustoca_na_ulazu, dx, podrelaksacija[0])

    tlak = zom_zokg_tlak(brzina, pretpostavka[1], p_visokotlacno[trenutak], p_niskotlacno[trenutak], pretpostavka[2],
        pretpostavka[3], vrijednost_proslog_trenutka[3], dx)
    tlak = eksplicitna(tlak, pretpostavka[1], podrelaksacija[1])

    entalpija = zoe_entalpija(brzina, vrijednost_proslog_trenutka[0], brzina_na_ulazu, pretpostavka[3],
        vrijednost_proslog_trenutka[3], gustoca_na_ulazu, tlak, vrijednost_proslog_trenutka[1],
        vrijednost_proslog_trenutka[2], entalpija_na_ulazu, pretpostavka[6], dx, pretpostavka[2],
        podrelaksacija[2])

    gustoca = CP.PropsSI('D', 'P', tlak, 'H', entalpija, rt)
    gustoca = eksplicitna(gustoca, pretpostavka[3], podrelaksacija[3])

    temperatura = CP.PropsSI('T', 'P', tlak, 'H', entalpija, rt)
    temperatura = eksplicitna(temperatura, pretpostavka[4], podrelaksacija[4])

    alpha_unutarnje_izracunato = alpha_unutarnje(tlak, entalpija, brzina * gustoca, pretpostavka[6]) # ne zapisuje se

    temperatura_stijenke = temp_stijenke(vrijednost_proslog_trenutka[5], temperatura, alpha_unutarnje_izracunato, dx)
    temperatura_stijenke = eksplicitna(temperatura_stijenke, pretpostavka[5], podrelaksacija[5])

    specificni_toplinski_tok = alpha_unutarnje_izracunato * (temperatura_stijenke - temperatura)
    specificni_toplinski_tok = eksplicitna(specificni_toplinski_tok, pretpostavka[6], podrelaksacija[6])

    rjesenje = np.stack((brzina, tlak, entalpija, gustoca, temperatura, temperatura_stijenke, specificni_toplinski_tok))
    return rjesenje

def rjesavac_trenutka(trenutak, vrijednost_proslog_trenutka, podrelaksacija, rezidual_limit):
    global broj_iteracija
    broj_iteracija = 0
    rjesenje_1 = vrijednost_proslog_trenutka
    uvjet = True
    while uvjet:
        rjesenje_2 = vanjska_iteracija(trenutak, rjesenje_1, vrijednost_proslog_trenutka, podrelaksacija)
        rezidual = manhattan_rjesenja(rjesenje_1, rjesenje_2, 1)
        uvjet = any( np.greater(rezidual, rezidual_limit) )
        rjesenje_1 = rjesenje_2
        broj_iteracija = broj_iteracija + 1
        print("broj unutarnjih iteracija: ", broj_iteracija)
        print("rezidual: ", rezidual)
    return rjesenje_1

def rjesavac(pocetni_trenutak, podrelaksacija, rezidual_limit):
    for t in range(pocetni_trenutak, n_t):
        prosli_trenutak = citaj_stanje(t - 1)
        rjesenje_trenutka = rjesavac_trenutka(t, prosli_trenutak, podrelaksacija, rezidual_limit)
        print(t)
        snimi_trenutak(rjesenje_trenutka, t)
    return print('rjeseno')

def snimi_trenutak(stanje, vremenski_trenutak):
    info = np.array([vremena[vremenski_trenutak], broj_iteracija, np.log10(zeljeni_rezidual), imp])
    putanja_direktorija = '../' + 'results/' + str(vremenski_trenutak) + '/'
    os.makedirs(putanja_direktorija, exist_ok=True)
    putanja_info = os.path.join(putanja_direktorija, 'informacije.csv')
    np.savetxt(putanja_info, info, delimiter=",")
    for i in range( broj_velicina ):
        putanja_velicine = putanja_direktorija + raspored_velicina[i] + '.csv'
        np.savetxt(putanja_velicine, stanje[i], delimiter=",")
    return print('\nvremenski trenutak: ', vremenski_trenutak, '\nbroj iteracija: ', broj_iteracija)

zeljeni_rezidual = 5e-4
imp = 0.01
eksp = 0.1

reziduali = np.repeat(zeljeni_rezidual, broj_velicina)

podrelaksacije = [imp, eksp, imp, eksp, eksp, eksp, eksp]
#podrelaksacije[1] = 0.1
#podrelaksacije[3] = 0.02
#podrelaksacije[6] = 0.5

rjesavac(1, podrelaksacije, reziduali)
