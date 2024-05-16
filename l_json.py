import json
from math import*

with open('synthese-fra.json','r') as fichier:
    donnes = json.load(fichier)

donnes_nuplets = []

for objet in donnes:
    nuplet = (
        objet["casConfirmes"], #NULL apcr, il faut le prednre en cmpte dans le code pour ne pas avoir d'erreur
        objet["deces"]
        )
    donnes_nuplets.append(nuplet)

#for nuplet in donnes_nuplets:
#    print(nuplet)

print(donnes_nuplets[10]) # envoie le nuplets que je veux
