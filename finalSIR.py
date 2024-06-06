import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize

# Charger les données
with open('synthese-fra.json', 'r') as fichier:
    donnes = json.load(fichier)

donnes_nuplets = []

for objet in donnes:
    if objet["casConfirmes"] is not None:
        nuplet = (
            objet["casConfirmes"], 
            objet["deces"],
            66992159 - objet["deces"]
        )
        donnes_nuplets.append(nuplet)

# Exemple
print(donnes_nuplets[11])

# modéle SEIR
def equations_SIR(y, t, a, b):
    S, I, R = y
    N = S  + I + R
    dSdt = -a * S * I / N
    dIdt =  a * S * I / N - b * I
    dRdt = b  * I
    return dSdt, dIdt, dRdt

# Conditions initiales
S0 = 999999
I0 = 1
R0 = 0
y0 = S0, I0, R0

# Temps
t = np.linspace(0, 200, 200)

# Paramétres initiaux
a, b = 0.3, 0.1

# Calcul d'erreur
def calculate_error(params, S0, I0, R0, donnes):
    a, b = params
    t = np.arange(len(donnes))
    solution = odeint(equations_SIR, [S0, I0, R0], t, args=(a, b))
    S, I, R = solution.T
    casConfirmes = np.array([x[0] for x in donnes])
    error = np.mean((I - casConfirmes)**2)  # erreur au carré (Euclidienne)
    return error

#Paramétres initiaux
ini_params = [a, b]

# Optimisation paramétres
result = minimize(calculate_error, ini_params, args=(S0, I0, R0, donnes_nuplets[:200]), method='L-BFGS-B')
a_opt, b_opt = result.x

# Simulation avec les paramétres optimisés
solution = odeint(equations_SIR, y0, t, args=(a_opt, b_opt))
S, I, R = solution.T

# affichage des résultats
"""
plt.figure(figsize=(12, 8))
plt.subplot(2, 1, 1)
plt.plot(t, S, label='Susceptible (SIR)')
plt.xlabel('Jours')
plt.ylabel('Nombre d individus')
plt.title('Susceptible (S)')
plt.legend()
"""
"""
plt.subplot(2, 2, 2)
plt.plot(t, E, label='Exposés (SIR model)')
plt.xlabel('Jours')
plt.ylabel('Nombre d individus')
plt.title('Exposés (E)')
plt.legend()"""

plt.subplot(2, 1, 1)
plt.plot(t, I, label='Infectieux (SIR)')
plt.xlabel('Jours')
plt.ylabel('Nombres d individus')
plt.title('Infectieux (I)')
plt.legend()

plt.subplot(2, 1, 1)
plt.plot(t, R, label='Rétablis (SIR)')
plt.xlabel('Jours')
plt.ylabel('Nombres d individus')
plt.title('Infectieux (I)')
plt.legend()

"""
plt.subplot(2, 2, 4)
plt.plot(t, R1, label='cas confirmés effacés(SIR model)')
plt.plot(t, R2, label='cas non confirmés effacés (SIR model)')
plt.xlabel('Jours')
plt.ylabel('Nombre d individus')
plt.title('(R1 et R2)')
plt.legend()
"""
plt.tight_layout()
plt.show()
