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

# Temps
t = np.linspace(0, 200, 200)

# Paramétres initiaux
a, b, c, f = 0.3, 0.1, 0.2, 0.8

# Calcul d'erreur
def calculate_error(params, S0, E0, I0, R10, R20, donnes):
    a, b, c, f = params
    t = np.arange(len(donnes))
    solution = odeint(equations_SEIR, [S0, E0, I0, R10, R20], t, args=(a, b, c, f))
    S, E, I, R1, R2 = solution.T
    casConfirmes = np.array([x[0] for x in donnes])
    error = np.mean((I - casConfirmes)**2)  # Mean squared error
    return error

#Paramétres initiaux
ini_params = [a, b, c, f]

# Optimisation paramétres
result = minimize(calculate_error, ini_params, args=(S0, E0, I0, R10, R20, donnes_nuplets[:200]), method='L-BFGS-B')
a_opt, b_opt, c_opt, f_opt = result.x

# Simulation avec les paramétres optimisés
solution = odeint(equations_SEIR, y0, t, args=(a_opt, b_opt, c_opt, f_opt))
S, E, I, R1, R2 = solution.T

# affichage des résultats
plt.figure(figsize=(12, 8))
plt.subplot(2, 2, 1)
plt.plot(t, S, label='Susceptible (SEIR)')
plt.xlabel('Jours')
plt.ylabel('Nombre d individus')
plt.title('Susceptible (S)')
plt.legend()

plt.subplot(2, 2, 2)
plt.plot(t, E, label='Exposés (SEIR model)')
plt.xlabel('Jours')
plt.ylabel('Nombre d individus')
plt.title('Exposés (E)')
plt.legend()

plt.subplot(2, 2, 3)
plt.plot(t, I, label='Infectieux (SEIR)')
plt.xlabel('Jours')
plt.ylabel('Nombres d individus')
plt.title('Infectieux (I)')
plt.legend()

plt.subplot(2, 2, 4)
plt.plot(t, R1, label='cas confirmés effacés(SEIR model)')
plt.plot(t, R2, label='cas non confirmés effacés (SEIR model)')
plt.xlabel('Jours')
plt.ylabel('Nombre d individus')
plt.title('(R1 et R2)')
plt.legend()

plt.tight_layout()
plt.show()