residus = {}    

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

# 53: {'nom': 'ILE', 'N': (-8.342, -3.574, 2.878), 'CA': (-7.109, -3.762, 3.695), 'C': (-6.017, -2.808, 3.195), 'O': (-5.963, -2.483, 2.021), 'H': (-8.423, -4.008, 2.0)}





# vérifier si résidus sont complets

# for numero, residue in residus.items():
#    if None in [residue["N"], residue["CA"], residue["C"], residue["O"]]:
#        print("Résidu incomplet :", numero)
     


# Calcule de distance   

from math import sqrt

def distance(atome1, atome2):
	x1, y1, z1 = atome1
	x2, y2, z2 = atome2
	return sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)


q1 = 0.42
q2 = 0.20
f = 332

liaisons_H = []

for i in residus:
	for j in residus:
		if i != j and residus[j]["H"] is not None:
		
			r_ON = distance(residus[i]["O"], residus[j]["N"]) # en A
			r_CH = distance(residus[i]["C"], residus[j]["H"])
			r_OH = distance(residus[i]["O"], residus[j]["H"])
			r_CN = distance(residus[i]["C"], residus[j]["N"])
			
			print("Résidus :", i, j)
			print("r_ON =", r_ON)
			print("r_CH =", r_CH)
			print("r_OH =", r_OH)
			print("r_CN =", r_CN)

			print("Résidu i =", i, residus[i]["nom"])
			print("Résidu j =", j, residus[j]["nom"])
			print("O vient du résidu", i)
			print("H vient du résidu", j)
			print("r_OH =", r_OH)


			Energie = q1 * q2 * f * ((1 / r_ON) + (1 / r_CH) - (1 / r_OH) - (1 / r_CN)) # en kcal/mol

			#print("Résidus :", i, j)
			#print("Énergie :", Energie)

			if Energie < -0.5:
				liaisons_H.append((i, j, Energie))


# VÉRIFICATION DES LIAISONS HYDROGÈNES

print("Liaisons hydrogène trouvées :")

for i, j, energie in liaisons_H:
    print(
        "Résidu", i, residus[i]["nom"],
        "-->",
        "Résidu", j, residus[j]["nom"],
        "| Énergie =", round(energie, 3), "kcal/mol"
    )
# Dans l'article ordre de grandeur attendu pour une liaison H idéale est : -3 kcal


print("Vérification des distances N-H")

for i in residus:
    N = residus[i]["N"]
    H = residus[i]["H"]

    if N is not None and H is not None:
        d_NH = distance(N, H)
        print(
            "Résidu", i,
            residus[i]["nom"],
            "N-H =", round(d_NH, 3), "Å"
        )









# 4 ASSIGNATION STRUCTURE SECONDAIRE

# 4-turn => difference de 4 entre résidu i et j 

turn4 = []
for i, j in liaisons_H:
	if j - i == 4:
		turn4.append((i,j))
		print("4-turn:", i, "-", j)


# B-bridge parallèle (pont parallèle)

bridge_paralleles = []

for i in residus:
	for j in residus:
		if i == j:
			continue
		if ((i - 1, j) in liaisons_H and (j , i + 1) in liaisons_H) or ((j - 1, i) in liaisons_H and (i , j + 1) in liaisons_H):
			bridge_paralleles.append((i, j))
			print("Pont parallèle:", i, "-", j) 



# B-bridge antiparallèle (pont antiparallèle)

bridge_antiparalleles = []

for i in residus:
	for j in residus:
		if i >= j:
			continue
		if ((i, j) in liaisons_H and (j, i) in liaisons_H) or ((i - 1, j + 1) in liaisons_H and (j - 1, i + 1) in liaisons_H):
			bridge_antiparalleles.append((i,j))
			print("Pont antiparallèle:", i, "-", j)









