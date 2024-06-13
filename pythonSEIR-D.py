import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize

# Fonction des équations SEIRD
def equations_SEIRD(y, t, beta, gamma, delta, epsilon):
    S, E, I, R, D = y
    N = S + E + I + R + D
    dSdt = -beta * S * I / N
    dEdt = beta * S * I / N - gamma * E
    dIdt = gamma * E - delta * I
    dRdt = delta * I - epsilon * R
    dDdt = epsilon * R
    return dSdt, dEdt, dIdt, dRdt, dDdt

# Chargement des données
donnees = [
    (200, 10),   # Exemple de données (cas confirmés, décès)
    (400, 20),
    # Ajoutez vos données ici ou chargez-les à partir d'un fichier
]

# Convertir les données en numpy arrays
donnees = np.array(donnees)

# Population totale
population = 67000000

# Conditions initiales
S0 = population - donnees[0, 0]  # Susceptibles initiaux
E0 = 1   # Exposés initiaux
I0 = donnees[0, 0]   # Infectieux initiaux
R0 = 0   # Rétablis initiaux
D0 = donnees[0, 1]   # Décès initiaux
y0 = S0, E0, I0, R0, D0

# Temps (300 jours)
jours = 300
t = np.linspace(0, jours, jours)

# Fonction pour calculer l'erreur entre les données réelles et simulées
def calculate_error(params, y0, t, donnees):
    beta, gamma, delta, epsilon = params
    solution = odeint(equations_SEIRD, y0, t, args=(beta, gamma, delta, epsilon))
    I, D = solution[:, 2], solution[:, 4]
    error = np.mean((I - donnees[:, 0])**2 + (D - donnees[:, 1])**2)
    return error

# Paramètres initiaux et optimisation
ini_params = [0.3, 0.1, 0.05, 0.01]  # Paramètres initiaux (à ajuster si nécessaire)
result = minimize(calculate_error, ini_params, args=(y0, t, donnees), method='L-BFGS-B')
beta_opt, gamma_opt, delta_opt, epsilon_opt = result.x

print(f"Paramètres optimisés : beta_opt = {beta_opt:.4f}, gamma_opt = {gamma_opt:.4f}, delta_opt = {delta_opt:.4f}, epsilon_opt = {epsilon_opt:.4f}")

# Simulation avec les paramètres optimisés
solution = odeint(equations_SEIRD, y0, t, args=(beta_opt, gamma_opt, delta_opt, epsilon_opt))
S, E, I, R, D = solution.T

# Tracé des résultats
plt.figure(figsize=(12, 8))

plt.plot(t, donnees[:, 0], 'o-', label='Cas Confirmés (réels)')
plt.plot(t, I, label='Infectieux (modèle SEIRD)')
plt.plot(t, donnees[:, 1], 's-', label='Décès (réels)')
plt.plot(t, D, label='Décès (modèle SEIRD)')

plt.xlabel('Jours')
plt.ylabel('Nombre d\'individus')
plt.title('Modèle SEIRD vs Données Réelles (300 jours)')
plt.legend()

plt.tight_layout()
plt.show()
