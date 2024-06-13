import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize
import sys

# Nombre de jours à simuler (limite maximale)
jours_max = 300

# Charger les données et limiter à 300 jours
with open('synthese-fra.json', 'r') as fichier:
    donnees = json.load(fichier)

donnees_nuplets = []
for objet in donnees:
    # Vérifier et remplacer les valeurs null par 0
    cas = objet["casConfirmes"] if objet["casConfirmes"] is not None else 0
    deces = objet["deces"] if objet["deces"] is not None else 0
    hospitalises = objet["hospitalises"] if objet["hospitalises"] is not None else 0
    ventiles = objet["reanimation"] if objet["reanimation"] is not None else 0
    # Ajouter le tuple à donnees_nuplets
    nuplet = (cas, deces, hospitalises, ventiles)
    donnees_nuplets.append(nuplet)

    # Arrêter après 300 jours
    if len(donnees_nuplets) >= jours_max:
        break

# Exemple
print(donnees_nuplets[11])

# Modèle SEIRHVD
def equations_SEIRHVD(y, t, a, b, c, d, e, f, g, h):
    S, E, I, H, V, R, D = y
    N = S + E + I + H + V + R + D
    dSdt = -a * S * I / N
    dEdt = a * S * I / N - b * E
    dIdt = b * E - c * I - d * I - e * I
    dHdt = d * I - f * H - g * H
    dVdt = e * I - h * V
    dRdt = c * I + f * H
    dDdt = g * H + h * V
    return dSdt, dEdt, dIdt, dHdt, dVdt, dRdt, dDdt

# Conditions initiales
S0 = 66992159 - donnees_nuplets[0][0]
E0 = 1
I0 = donnees_nuplets[0][0]
H0 = 0
V0 = 0
R0 = 0
D0 = 0
y0 = S0, E0, I0, H0, V0, R0, D0

# Temps limité à 300 jours
t = np.linspace(0, jours_max - 1, jours_max)

# Paramètres initiaux
a, b, c, d, e, f, g, h = 0.3, 0.1, 0.2, 0.05, 0.01, 0.05, 0.01, 0.01

# Calcul d'erreur
def calculate_error(params, S0, E0, I0, H0, V0, R0, D0, donnees):
    alpha, beta, gamma, delta, eta, zeta, theta, kappa = params
    solution = odeint(equations_SEIRHVD, [S0, E0, I0, H0, V0, R0, D0], t, args=(alpha, beta, gamma, delta, eta, zeta, theta, kappa))
    S, E, I, H, V, R, D = solution.T
    casConfirmes = np.array([x[0] for x in donnees])
    deces = np.array([x[1] for x in donnees])
    hospitalises = np.array([x[2] for x in donnees])
    ventiles = np.array([x[3] for x in donnees])
    error = np.mean((I - casConfirmes)**2 + (D - deces)**2 + (H - hospitalises)**2 + (V - ventiles)**2)
    return error

# Paramètres initiaux
ini_params = [a, b, c, d, e, f, g, h]

# Optimisation des paramètres
result = minimize(calculate_error, ini_params, args=(S0, E0, I0, H0, V0, R0, D0, donnees_nuplets[:jours_max]), method='L-BFGS-B')
a_opt, b_opt, c_opt, d_opt, e_opt, f_opt, g_opt, h_opt = result.x

print("Paramètres Optimisés :")
print(f"a_opt = {a_opt:.4f}, b_opt = {b_opt:.4f}, c_opt = {c_opt:.4f}, d_opt = {d_opt:.4f}, e_opt = {e_opt:.4f}, f_opt = {f_opt:.4f}, g_opt = {g_opt:.4f}, h_opt = {h_opt:.4f}")

# Simulation avec les paramètres optimisés
solution_opt = odeint(equations_SEIRHVD, y0, t, args=(a_opt, b_opt, c_opt, d_opt, e_opt, f_opt, g_opt, h_opt))
S_opt, E_opt, I_opt, H_opt, V_opt, R_opt, D_opt = solution_opt.T

# Affichage des résultats
plt.figure(figsize=(12, 8))

# Données réelles
plt.subplot(2, 2, 1)
plt.plot(t, [x[0] for x in donnees_nuplets], label='Cas Confirmés (réels)')
plt.plot(t, [x[1] for x in donnees_nuplets], label='Décès (réels)')
plt.plot(t, [x[2] for x in donnees_nuplets], label='Hospitalisés (réels)')
plt.plot(t, [x[3] for x in donnees_nuplets], label='En Réanimation (réels)')
plt.xlabel('Jours')
plt.ylabel('Nombre d\'individus')
plt.title('Données Réelles')
plt.legend()

# Modèle avec paramètres optimisés
plt.subplot(2, 2, 2)
plt.plot(t, I_opt, label='Infectieux (paramètres optimisés)')
plt.plot(t, D_opt, label='Décès (paramètres optimisés)')
plt.plot(t, H_opt, label='Hospitalisés (paramètres optimisés)')
plt.plot(t, R_opt, label='Rétablis (paramètres optimisés)')
plt.plot(t, V_opt, label='Ventilés (paramètres optimisés)')
plt.xlabel('Jours')
plt.ylabel('Nombre d\'individus')
plt.title('Modèle avec Paramètres Optimisés')
plt.legend()

plt.tight_layout()
plt.show()
