import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Fonction décrivant les équations différentielles du modèle SEIR avec deux compartiments retirés
def equations_SEIR_R(y, t, a, b, c, f):
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
solution = odeint(equations_SEIR_R, y0, t, args=(a, b, c, f))
S, E, I, R1, R2 = solution.T

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
