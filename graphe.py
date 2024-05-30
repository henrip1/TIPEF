import numpy as np
import matplotlib.pyplot as plt

import json
from math import*


deb = 0
fin = 500


with open('synthese-fra.json','r') as fichier:
    donnes = json.load(fichier)

donnes_nuplet = []

for objet in donnes:
    nuplet = (
        objet["casConfirmes"], #NULL apcr, il faut le prednre en cmpte dans le code pour ne pas avoir d'erreur
        objet["deces"],
        66992159 - objet["deces"] #population francaise en 2018 selon data.gouv déduit des décès
        )
    donnes_nuplet.append(nuplet)

    
def donnees_du_jour (jour):
    return donnes_nuplet[jour] # envoie le nuplets que je veux

def donnees_sur_interv (deb, fin):
    tab = []
    for i in range (deb, fin):
        tab.append(donnees_du_jour(i))
    return tab    

plt.plot(np.linspace(deb, fin, fin-deb), donnees_sur_interv(deb, fin))
plt.show()
