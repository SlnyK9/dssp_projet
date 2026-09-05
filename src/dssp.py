residus = {}    

# 1 EXTRACTION DES COORDOONNÉES ATOMIQUES

with open("thioredoxin.pdb", "r") as fichier:
    for ligne in fichier:

        if ligne.startswith("ATOM"):

            nom_atome = ligne[12:16].strip()
            nom_residue = ligne[17:20].strip()
            nombre_residue = int(ligne[22:26])

            x = float(ligne[30:38])
            y = float(ligne[38:46])
            z = float(ligne[46:54])

            if nom_atome in ["N", "CA", "C", "O", "H"]:

                if nombre_residue not in residus:
                    residus[nombre_residue] = {
                        "nom": nom_residue,
                        "N": None,
                        "CA": None,
                        "C": None,
                        "O": None,
                        "H" : None
                    }

                residus[nombre_residue][nom_atome] = (x, y, z) 

# print(residus)


# 2 VÉRIFICATION PRÉSENCE DE NH + CALCULE DE dist_NH

from math import sqrt

def distance(atome1, atome2):
	x1, y1, z1 = atome1
	x2, y2, z2 = atome2
	return sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)


#print("Vérification des distances N-H")

for i in residus:
    N = residus[i]["N"]
    H = residus[i]["H"]

    if N is not None and H is not None:
        d_NH = distance(N, H)
   #     print(
    #        "Résidu", i,
     #       residus[i]["nom"],
      #      "N-H =", round(d_NH, 3), "Å"
     #   )
   # elif N is not None and H is None:
   #     print("Abscense de NH -->", i, residus[i]["nom"])


# 3 CALCULE DE DISTANCE + ÉNERGIE 

q1 = 0.42
q2 = 0.20
f = 332

liaisons_H = []

for i in residus:
	for j in residus:
		if (i != j 
		and residus[j]["O"]
		and residus[j]["C"]
		and residus[j]["N"]
		and residus[j]["H"] is not None):
		
			r_ON = distance(residus[i]["O"], residus[j]["N"]) # en A
			r_CH = distance(residus[i]["C"], residus[j]["H"])
			r_OH = distance(residus[i]["O"], residus[j]["H"])
			r_CN = distance(residus[i]["C"], residus[j]["N"])
			
		#	print("Résidus :", i, j)
		#	print("r_ON =", r_ON)
		#	print("r_CH =", r_CH)
		#	print("r_OH =", r_OH)
		#	print("r_CN =", r_CN)


			Energie = q1 * q2 * f * ((1 / r_ON) + (1 / r_CH) - (1 / r_OH) - (1 / r_CN))

		#	print("Résidus :", i, j)
		#	print("Énergie :", Energie)

			if Energie < -0.5:
				liaisons_H.append((i, j))

#for i, j, energie in liaisons_H:
#	print("Résidu", i, residus[i]["nom"], "-->", "Résidu", j, residus[j]["nom"],
#        "| Énergie =", energie, "kcal/mol"
   # )

# Dans l'article ordre de grandeur attendu pour une liaison H idéale est : -3 kcal


# 4 ASSIGNATION STRUCTURE SECONDAIRE

# 4-turn => difference de 4 entre résidu i et j 

turn4 = []
for i, j in liaisons_H:
	if j - i == 4:
		turn4.append((i,j))
		print("4-turn:", i, "-", j)
		
helices_alpha = []
for i, j in turn4:	
	if (i + 1, j + 1) in turn4:
		helices_alpha.append((i, j))
		print("Hélice alpha :", i, "-", j)



# B-bridge parallèle (pont parallèle)

bridge_paralleles = []

for i in residus:
    for j in residus:
        if i >= j:
            continue

        if abs(i - j) < 3:
            continue

        if ((i - 1, j) in liaisons_H and (j, i + 1) in liaisons_H) \
        or ((j - 1, i) in liaisons_H and (i, j + 1) in liaisons_H):

            bridge_paralleles.append((i, j))
            print("Pont parallèle :", i, "-", j)

# Ladders parallèle


ladders_paralleles = []

for i, j in bridge_paralleles:

    # on cherche seulement le début d'une ladder
    if (i - 1, j - 1) in bridge_paralleles:
        continue

    debut_i = i
    debut_j = j

    while (i + 1, j + 1) in bridge_paralleles:
        i = i + 1
        j = j + 1

    if i > debut_i:
        ladders_paralleles.append((debut_i, debut_j, i, j))


for ladder in ladders_paralleles:
    print("Ladder parallèle :", ladder)



# B-bridge antiparallèle (pont antiparallèle)

bridge_antiparalleles = []

for i in residus:
    for j in residus:
        if i >= j:
            continue

        if abs(i - j) < 3:
            continue

        if ((i, j) in liaisons_H and (j, i) in liaisons_H) \
        or ((i - 1, j + 1) in liaisons_H and (j - 1, i + 1) in liaisons_H):

            bridge_antiparalleles.append((i, j))
        #    print("Pont antiparallèle :", i, "-", j)


# Ladders antiparallèle

ladders_antiparalleles = []

for i, j in bridge_antiparalleles:

    # on cherche seulement le début d'une ladder
    if (i - 1, j + 1) in bridge_antiparalleles:
        continue

    debut_i = i
    debut_j = j

    while (i + 1, j - 1) in bridge_antiparalleles:
        i = i + 1
        j = j - 1

    if i > debut_i:
        ladders_antiparalleles.append((debut_i, debut_j, i, j))


for ladder in ladders_antiparalleles:
    print("Ladder antiparallèle :", ladder)
    
    



