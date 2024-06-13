import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize
import sys

jours = int(sys.argv[1])

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

# Modèle SEIR
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
S0 = 66992159 - donnes_nuplets[0][0]
E0 = 1
I0 = donnes_nuplets[0][0]
R10 = 0
R20 = 0
y0 = S0, E0, I0, R10, R20

# Temps
t = np.linspace(0, jours, jours)

# Paramètres initiaux
a, b, c, f = 0.3, 0.1, 0.2, 0.2

# Calcul d'erreur
def calculate_error(params, S0, E0, I0, R10, R20, donnes):
    a, b, c, f = params
    t = np.arange(len(donnes))
    solution = odeint(equations_SEIR, [S0, E0, I0, R10, R20], t, args=(a, b, c, f))
    S, E, I, R1, R2 = solution.T
    casConfirmes = np.array([x[0] for x in donnes])
    error = np.mean((I - casConfirmes)**2)  # Mean squared error
    return error

# Paramètres initiaux
ini_params = [a, b, c, f]

# Optimisation des paramètres
result = minimize(calculate_error, ini_params, args=(S0, E0, I0, R10, R20, donnes_nuplets[:jours]), method='L-BFGS-B')
a_opt, b_opt, c_opt, f_opt = result.x

print("a =", a_opt, ", b =", b_opt, ", c =", c_opt, ", f =", f_opt)

# Simulation avec les paramètres optimisés
solution = odeint(equations_SEIR, y0, t, args=(a_opt, b_opt, c_opt, f_opt))
S, E, I, R1, R2 = solution.T

I3 = np.array([x[0] for x in donnes_nuplets[0:jours]])

# Simulation avec les paramètres initiaux pour comparaison
solution2 = odeint(equations_SEIR, y0, t, args=(a, b, c, f))
S1, E1, I1, R11, R21 = solution2.T

# Affichage des résultats
plt.figure(figsize=(12, 8))

# Données réelles
plt.subplot(2, 2, 1)
plt.plot(t, I3, label='Infectieux (données réelles)')
plt.ylabel('Nombre d\'individus')
plt.legend()

# Modèle SEIR optimisé
plt.subplot(2, 2, 2)
plt.plot(t, E, label='Exposés (modèle SEIR)')
plt.plot(t, I, label='Infectieux (modèle SEIR)')
plt.plot(t, R1, label='Cas confirmés récupérés (modèle SEIR)')
plt.plot(t, R2, label='Cas non confirmés récupérés (modèle SEIR)')
plt.xlabel('Jours')
plt.ylabel('Nombre d\'individus')
plt.title('Modèle SEIR optimisé')
plt.legend()

# Modèle SEIR initial
plt.subplot(2, 2, 3)
plt.plot(t, E1, label='Exposés (modèle SEIR)')
plt.plot(t, I1, label='Infectieux (modèle SEIR)')
plt.plot(t, R11, label='Rétablis (modèle SEIR)')
plt.plot(t, R21, label='Cas non confirmés récupérés (modèle SEIR)')
plt.xlabel('Jours')
plt.ylabel('Nombre d\'individus')
plt.title('Modèle SEIR initial')
plt.legend()

plt.tight_layout()
plt.show()
