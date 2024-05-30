import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize
from math import*

with open('synthese-fra.json','r') as fichier:
    donnes = json.load(fichier)

donnes_nuplets = []

for objet in donnes:
    nuplet = (
        objet["casConfirmes"], #NULL apcr, il faut le prednre en cmpte dans le code pour ne pas avoir d'erreur
        objet["deces"],
        66992159 - objet["deces"]
        )
    donnes_nuplets.append(nuplet)

#for nuplet in donnes_nuplets:
#    print(nuplet)

print(donnes_nuplets[11]) # envoie le nuplets que je veux





# Fonction décrivant les équations différentielles du modèle SEIR avec deux compartiments retirés
def equations_SEIR(y, t, a, b, c, f):
    S, E, I, R1, R2 = y
    N = S + E + I + R1 + R2
    dSdt = -a * S * I / N
    dEdt = a * S * I / N - b * E
    dIdt = b * E - c * I
    dR1dt = f * c * I
    dR2dt = (1 - f) * c * I
    return dSdt, dEdt, dIdt, dR1dt, dR2dt

# Conditions initiales
S0 = 999999
E0 = 1
I0 = 0
R10 = 0
R20 = 0
y0 = S0, E0, I0, R10, R20

# Vecteur de temps
t = np.linspace(0, 200, 200)

# Paramètres du modèle
a = 0.3  # Taux de contact effectif
b = 0.1  # Taux auquel les personnes infectées en phase latente deviennent infectieuses
c = 0.2  # Taux moyen auquel les personnes infectieuses sont isolées
f = 0.8  # Fraction d'individus infectieux comptabilisés parmi les cas confirmés

# Exécution de la simulation
solution = odeint(equations_SEIR, y0, t, args=(a, b, c, f))
S, E, I, R1, R2 = solution.T



def calculate_error(params, S0, E0, I0, R10, R20, donnes):
    a, b, c, f = params
    _, _, I, _, _ = equations_SEIR(y, days, a, b, c, f)
    IR,_,_ = donnes
    error = np.mean((I-IR)**2)
    return error

S0 = 999999
E0 = 1
I0 = 0
R10 = 0
R20 = 0
y0 = S0, E0, I0, R10, R20
a_v = []
b_v = []
c_v = []
f_v = []
day = 200
ini_params = [0.3,0.1,0.2,0.8]

for day in range(1, day + 1):

    result = minimize(calculate_error, ini_params, args=( S0, E0, I0, R10, R20, day, donnes_nuplets[day+1] ))#jour un par

    a_opt, b_opt, c_opt, f_opt = result.x

    a_v.append(a_opt)
    b_v.append(b_opt)
    c_v.append(c_opt)
    f_v.append(f_opt)



    S, R, I, R1, R2 = SEIR_model(y, 2, a_opt, b_opt, c_opt, f_opt)
    
        
    ini_params = [a_opt, b_opt, c_opt, f_opt]

    plt.figure(figsize=(12,8))
    
    plt.subplot(2,2,1)
    plt.plot(S_prediction, label='S (modèle SEIR)')
    plt.plot(external_date, label='Données réelles (S)', linestyle='dashed')
    plt.xlabel('Jours')
    plt.ylabel('Nombre de personnes')
    plt.title('Sains (S)')
    plt.legend()

    plt.subplot(2,2,2)
    plt.plot(S_prediction, label='E (modèle SEIR)')
    plt.plot(external_date, label='Données réelles (E)', linestyle='dashed')
    plt.xlabel('Jours')
    plt.ylabel('Nombre de personnes')
    plt.title('Exposés (E)')
    plt.legend()

    plt.tight_layout()
    plt.show()



# Visualisation des résultats
plt.figure(figsize=(10, 6))
plt.plot(t, S, label='Susceptibles')
plt.plot(t, E, label='Infectés en phase latente')
plt.plot(t, I, label='Infectieux sans protection')
plt.plot(t, R1, label='Retirés confirmés')
plt.plot(t, R2, label='Retirés non confirmés')
plt.xlabel('Temps (jours)')
plt.ylabel('Nombre de personnes')
plt.title('Modèle SEIR avec deux compartiments retirés')
plt.legend()
plt.grid(True)
plt.show()
