
# --- 1 EXTRACTION DES COORDOONNÉES ATOMIQUES ---

residus = {}

with open("../data/thioredoxin.pdb", "r") as fichier:
    for ligne in fichier:

        if ligne.startswith("ATOM"):

            nom_atome = ligne[12:16].strip()
            nom_residue = ligne[17:20].strip()
            numero_residue = int(ligne[22:26])

            x = float(ligne[30:38])
            y = float(ligne[38:46])
            z = float(ligne[46:54])

            if nom_atome in ["N", "CA", "C", "O", "H"]:

                if numero_residue not in residus:
                    residus[numero_residue] = {
                        "nom": nom_residue,
                        "N": None,
                        "CA": None,
                        "C": None,
                        "O": None,
                        "H" : None
                    }

                residus[numero_residue][nom_atome] = (x, y, z) 

# print(residus)




# --- 2 CALCULE DE DISTANCE + ÉNERGIE ---

from math import sqrt

def distance(coord_atome1, coord_atome2):
	x1, y1, z1 = coord_atome1
	x2, y2, z2 = coord_atome2
	return sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

q1 = 0.42
q2 = 0.20
f = 332

liaisons_H = []

for residu1 in residus:
	for residu2 in residus:
		if (residu1 != residu2 and residus[residu1]["O"] is not None and residus[residu1]["C"] is not None and residus[residu2]["N"] is not None and residus[residu2]["H"] is not None):
		
			r_ON = distance(residus[residu1]["O"], residus[residu2]["N"]) 
			r_CH = distance(residus[residu1]["C"], residus[residu2]["H"])
			r_OH = distance(residus[residu1]["O"], residus[residu2]["H"])
			r_CN = distance(residus[residu1]["C"], residus[residu2]["N"])

			Energie = q1 * q2 * f * ((1 / r_ON) + (1 / r_CH) - (1 / r_OH) - (1 / r_CN))

		#	print("Résidus :", residu1, residu2, "--> E =", Energie)

			if Energie < -0.5:
				liaisons_H.append((residu1, residu2))
				

			#	print("Résidu", residu1, residus[residu1]["nom"], "-->", "Résidu", residu2, residus[residu2]["nom"], "| Énergie =", Energie, "kcal/mol")


print("Nombre total de liaisons hydrogènes :", len(liaisons_H))



# --- 3 ASSIGNATION STRUCTURE SECONDAIRE ---

# A) 4-turn - Hélice alpha

turn4 = []
for residu1, residu2 in liaisons_H:
	if residu2 - residu1 == 4:
		turn4.append((residu1,residu2))
		print("4-turn:", residu1, "-", residu2)


def trouver_helices(turn4):
    helices = []

    if len(turn4) == 0:
        return helices

    debut = turn4[0][0]
    fin = turn4[0][1]
    nombre_turns = 1

    for numero_turn in range(1, len(turn4)):

        if (turn4[numero_turn][0] == turn4[numero_turn - 1][0] + 1
                and turn4[numero_turn][1] == turn4[numero_turn - 1][1] + 1):

            fin = turn4[numero_turn][1]
            nombre_turns += 1

        else: 
            if nombre_turns >= 2:
                helices.append((debut, fin))

            debut = turn4[numero_turn][0]
            fin = turn4[numero_turn][1]
            nombre_turns = 1

    if nombre_turns >= 2:
        helices.append((debut, fin))

    return helices
    
    
helices_alpha = trouver_helices(turn4)

for debut, fin in helices_alpha:
    print("Hélice alpha :", debut, "-", fin)




# B) B-bridge parallèle (pont parallèle) / B-bridge antiparallèle (pont antiparallèle)

bridge_paralleles = []

for i in residus:
    for j in residus:
        if i < j and abs(i - j) >= 3:  
            if ((i - 1, j) in liaisons_H and (j, i + 1) in liaisons_H) \
                    or ((j - 1, i) in liaisons_H and (i, j + 1) in liaisons_H):  
                bridge_paralleles.append((i, j))
                print("Pont parallèle :", i, "-", j)

print(bridge_paralleles)


bridge_antiparalleles = []

for i in residus:
    for j in residus:
        if i < j and abs(i - j) >= 3:
            if ((i, j) in liaisons_H and (j, i) in liaisons_H) \
                    or ((i - 1, j + 1) in liaisons_H and (j - 1, i + 1) in liaisons_H):
                bridge_antiparalleles.append((i, j))
                print("Pont antiparallèle :", i, "-", j)

print(bridge_antiparalleles)



# C) Ladders parallèle / Ladders antiparallèle

ladders = {}
nom_ladders_parallele = "a"

for residu1, residu2 in bridge_paralleles:

    if (residu1 - 1, residu2 - 1) in bridge_paralleles:
        continue

    bridges = [(residu1, residu2)]

    while (residu1 + 1, residu2 + 1) in bridge_paralleles:
        residu1 += 1
        residu2 += 1
        bridges.append((residu1,residu2))
       
    ladders[nom_ladders_parallele] = bridges
    nom_ladders_parallele = chr(ord(nom_ladders_parallele) + 1)
   
   
   
nom_ladders_antiparallele = "A"

for residu1, residu2 in bridge_antiparalleles:

    if (residu1 - 1, residu2 + 1) in bridge_antiparalleles:
        continue

    bridges = [(residu1, residu2)]

    while (residu1 + 1, residu2 - 1) in bridge_antiparalleles:
        residu1 += 1
        residu2 -= 1
        bridges.append((residu1, residu2))

    ladders[nom_ladders_antiparallele] = bridges
    nom_ladders_antiparallele = chr(ord(nom_ladders_antiparallele) + 1)
    
for nom_ladder, bridges in ladders.items():
	print("Ladder", nom_ladder, "-->", bridges)


# print(ladders)


# D) Sheets

# 1e étape
residus_communs = []
noms = list(ladders)

for i in range(len(noms)):
    for j in range(i + 1, len(noms)):

        ladder1 = noms[i]
        ladder2 = noms[j]

        residus1 = set()
        residus2 = set()

        for bridge in ladders[ladder1]:
            residus1.add(bridge[0])
            residus1.add(bridge[1])

        for bridge in ladders[ladder2]:
            residus2.add(bridge[0])
            residus2.add(bridge[1])

        if residus1 & residus2:
        	residus_communs.append((ladder1,ladder2))
        	print("Ladder", ladder1, "Ladder", ladder2, "-->", residus1 & residus2)

# 2e étape
Sheets = {}
nom_sheet = "A"

for ladder1, ladder2 in residus_communs:

    connecte = False

    for sheet in Sheets:

        if ladder1 in Sheets[sheet] or ladder2 in Sheets[sheet]:

            if ladder1 not in Sheets[sheet]:
                Sheets[sheet].append(ladder1)

            if ladder2 not in Sheets[sheet]:
                Sheets[sheet].append(ladder2)

            connecte = True

    if not connecte:
        Sheets[nom_sheet] = [ladder1, ladder2]
        nom_sheet = chr(ord(nom_sheet) + 1)

# print(Sheets)


# --- 5 SUMMARY ---


SUMMARY = {}
for residu in residus:
	SUMMARY[residu] = ""
	
for nom_ladder in ladders:
	if len(ladders[nom_ladder]) == 1:
		lettre = "B"
	else:
		lettre = "E"
		
	for bridge in ladders[nom_ladder]:
		for residu in bridge:
			SUMMARY[residu] = lettre
	
for debut, fin in helices_alpha:
    for residu in range(debut, fin + 1):
        SUMMARY[residu] = "H"


print(SUMMARY)



# --- 6 RÉSULTAT ---

import pandas as pd
from codes_acides_amines import codes_aa

tableau = []

for numero_residu in SUMMARY:
	nom_residu = residus[numero_residu]["nom"]
	code = codes_aa[nom_residu]
	
	tableau.append([numero_residu, code, SUMMARY[numero_residu]])

df = pd.DataFrame(tableau, columns=["RESIDUE", "SEQUENCE" , "SUMMARY"])



