import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sys

# Paramètres du modèle
beta = float(sys.argv[1])  # Taux de transmission
gamma = float(sys.argv[2])  # Taux de récupération
population = float(sys.argv[3])  # Taille de la population initiale
infecte = float(sys.argv[4])  # nombre initial de personnes infectées
jours = int(sys.argv[5])  # Nombre de jours à simuler
temps_new_susceptibles = int(sys.argv[6])


ratio_infecte = infecte/population
ratio_sain = 1 - ratio_infecte

# Fonction décrivant les équations différentielles du modèle SIR
def equations_SIR(y, t, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I
    dIdt = beta * S * I  - (gamma * I)
    dRdt = gamma * I
    return dSdt, dIdt, dRdt

# Simulation du modèle SIR
def simulation_SIR(beta, gamma, population, infectes_initiaux, jours):
    # Conditions initiales
    S0 = ratio_sain - ratio_infecte
    I0 =  ratio_infecte
    R0 = 0

    # Vecteur de temps
    t = np.linspace(0, jours, jours)

    # Résolution des équations différentielles
    y0 = S0, I0, R0
    solution = odeint(equations_SIR, y0, t, args=(beta, gamma))
    S, I, R = solution.T
    return t, S, I, R

# Exécution de la simulation
t, S, I, R = simulation_SIR(beta, gamma, ratio_sain, ratio_infecte, jours)

# Visualisation des résultats
plt.figure(figsize=(10, 6))
plt.plot(t, S, label='Susceptibles')
plt.plot(t, I, label='Infectés')
plt.plot(t, R, label='Récupérés')
plt.xlabel('Temps (jours)')
plt.ylabel('Pourcentage de personnes')
plt.title('Modèle SIR de propagation d\'une maladie infectieuse')
plt.legend()
plt.grid(True)
plt.show()

