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
a, b = 0.23523446569746287, 0.1

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
result = minimize(calculate_error, ini_params, args=(S0, I0, R0, donnes_nuplets[:200]), tol = 0.01 ,method='L-BFGS-B')
a_opt, b_opt = result.x
print ("a = ",a_opt,", b =",b_opt)

# Simulation avec les paramétres optimisés
solution = odeint(equations_SIR, y0, t, args=(a_opt, b_opt))
S, I, R = solution.T

solution2 = odeint(equations_SIR, y0, t, args=(a,b))
S1, I1 ,R1 = solution2.T

# affichage des résultats


plt.subplot(2, 2, 1)
plt.plot(t, I, label='Infectieux (SIR)')
plt.legend()

plt.subplot(2, 2, 1)
plt.plot(t, R, label='Rétablis (SIR)')
plt.xlabel('Jours')
plt.ylabel('Nombres d individus')
plt.title('Optimisé')
plt.legend()



plt.subplot(2, 2, 2)
plt.plot(t,donnes_nuplets[:200] , label='Infectieux (SIR)')
plt.legend()

plt.subplot(2, 2, 2)
plt.plot(t, R, label='Rétablis (SIR)')
plt.xlabel('Jours')
plt.ylabel('Nombres d individus')
plt.title('Données réelles')
plt.legend()



plt.subplot(2, 2, 3)
plt.plot(t, I1, label='Infectieux (SIR)')
plt.legend()

plt.subplot(2, 2, 3)
plt.plot(t, R1, label='Rétablis (SIR)')
plt.xlabel('Jours')
plt.ylabel('Nombres d individus')
plt.title('Para ini')
plt.legend()



plt.tight_layout()
plt.show()
