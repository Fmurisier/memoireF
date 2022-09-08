"""
redock GalaxyDock3
Editor : Murisier Frederic
Date : 2022
effectue le test 0 pour GalaxyDock3, docking de chaque recepteur avec son / ses ligands et evaluation si le docking a
bien fonctionne

Python version : 3.8


"""
import os
from os import listdir
from os.path import isfile, join

path = os.getcwd()
terminal = False
# variable environnement
PATH = '/home/rene/bin/GalaxyDock3/'


def liste_file(patern, file_path=''):
    """
    effectue la liste des fichier comprenant un certain motif donne en argument dans le dossier donne en argument si la
    liste des fichiers doit etre faites dans un sous dossier de repertoire dans lequel le script se trouve
    :param patern: motif a rechercher dans le nom des fichiers
    :param file_path: chemin du repertoire souhaite, si vide = repertoire dans lequel le script se trouve
    :return:
    """
    list_file = []

    # fait la liste des fichier dans le repertoire ou se trouve script.py
    fichiers = [f for f in listdir(path + file_path) if isfile(join(path + file_path, f))]

    # pour chaque nom de fichier on separe le nom par '.' pour recuperer et stocker tous les noms de fichier dans une
    # liste
    for ligne in fichiers:
        if patern in ligne:
            name = ligne.split('.')
            list_file.append(name[0])
    list_file.sort()

    return list_file


def dock_galaxy(e):
    """
    effectue le docking avec galaxyDock en allant prealablement chercher les informations de la box 
    :param e: 
    :return: 
    """
    if terminal:
        if e not in listdir(path + '/GalaxyD'):
            os.system('mkdir GalaxyD/' + e)

    file = open('VINA/BoxTxt/' + e.split('_')[0] + '_box.txt', 'r').readlines()
    cx = file[0].split()[-1]
    cy = file[1].split()[-1]
    cz = file[2].split()[-1]
    sx = file[3].split()[-1]
    sy = file[4].split()[-1]
    sz = file[5].split()[-1]

    # galaxy dock
    commande = '/home/rene/bin/GalaxyDock3/script/run_GalaxyDock3.py -d /home/rene/bin/GalaxyDock3 -p RECEPTEUR/' + \
               e.split('_')[0] + '_recepteur.pdb -l LIGAND/MOL2/' + e + '.mol2 -x ' + cx + ' -y ' + cy + ' -z ' + \
               cz + ' -size_x ' + sx + ' -size_y ' + sy + ' -size_z ' + sz + ' --n_proc 1'

    if terminal:
        os.system(commande)
    else:
        print('docking--->' + commande)

    # split le resultat
    commande = 'python2 /home/rene/bin/GalaxyDock3/script/split_mol2.py GD3_cl.mol2 GalaxyD/' + e + '/' + e
    if terminal:
        os.system(commande)
    else:
        print('split---->' + commande)

    # ouverture ficher resultat, suppression trois derniere colonne, ajout de 'rmsd' comme derniere colonne
    result = open('GD3_cl.E.info', 'r').readlines()
    resultfinal = open('GalaxyD/' + e + '_result_galaxy.txt', 'w')
    resultfinal.write(result[0] + result[1][:-15] + 'RMSD\n' + result[0])

    # calcul RMSD
    liste_galaxy = liste_file(e, '/GalaxyD/' + e)
    nbr_model = len(liste_galaxy)
    compteur = 0
    for model in liste_galaxy:
        commande = '/home/rene/bin/DockRMSD/DockRMSD LIGAND/MOL2/' + e + '.mol2 GalaxyD/' + e + '/' + model + \
                   '.mol2 > rescore_result.txt'

        if terminal:
            os.system(commande)
            # else:
            print('rmsd ---> ' + commande)

        # extraction du RMSD
        rescore_lignes = open('rescore_result.txt', 'r').readlines()

        dockrmsd = 'NA'
        for i in rescore_lignes:
            if 'Calculated Docking RMSD:' in i:
                dockrmsd = i.split()[3]

        # ecriture du RMSD dans le fichier
        resultfinal.write(result[3 + compteur][:-27] + dockrmsd + '\n')
        compteur += 1

    resultfinal.close()


def triage():
    """
    cree un nouveau fichier result_galaxy_tried.txt contenant les resultats de galaxyDock trie par energie croissant
    :return: 
    """
    l = liste_file('result_galaxy.txt', '/GalaxyD')
    for s in l:
        print(s)
        file = result = open('GalaxyD/' + s + '.txt', 'r').readlines()
        c = 0
        for i in file:
            c += 1
        print(c)
        co = 'head -3 GalaxyD/' + s + '.txt > GalaxyD/' + s + '_tried.txt'
        commande = 'tail -' + str(c - 3) + ' GalaxyD/' + s + '.txt | sort -n -k 2 >> GalaxyD/' + s + '_tried.txt'

        if terminal:
            os.system(co)
            os.system(commande)
        else:
            print(co)
            print(commande)


def reconstruction(e, model):
    """
    fonction pour la reconstruction des complexes
    :param e: 
    :param model: 
    :return: 
    """
    # reconstruction du complex
    commande = 'obabel -imol2 GalaxyD/' + e + '/' + model + ' -O fichier_temporaire.pdb'
    if terminal:
        os.system(commande)
    else:
        print(commande)


def conversion_mol2():
    """
    conversion en mol de tous les fichiers ligand_smile
    :return: 
    """
    if terminal:
        if 'MOL2' not in listdir(path + '/SMILE'):
            os.system('mkdir SMILE/MOL2')
    liste_smile = liste_file('smile', '/SMILE/')
    for smi in liste_smile:
        commande = '/usr/bin/obabel -ipdb SMILE/' + smi + '.pdb -O SMILE/MOL2/' + smi + '.mol2'
        if terminal:
            os.system(commande)


def dockG():
    """
    effectue le docking pour toutes les structures contenue dans la liste
    :return: 
    """
    liste_gal = liste_file('.', '/LIGAND/MOL2')
    #liste_gal = ['5NZ1_ligand1', '5NZ1_ligand2', '5NZ1_ligand3', '6HNX_ligand1', '6HNZ_ligand1', '6HO0_ligand1', '6HO1_ligand1', '6HO2_ligand1', '6HO3_ligand1', '6HO4_ligand1', '6HO5_ligand1', '6HO6_ligand1', '6HO7_ligand1', '6HO8_ligand1', '6HO9_ligand1', '6HOA_ligand1', '6HOB_ligand1', '6HOC_ligand1', '6HOD_ligand1', '6HOE_ligand1', '6HOF_ligand1', 'BDM15048_ligand1', 'BDM15048_ligand2', 'Cmpd1_ligand1', 'Cmpd1_ligand2', 'Cmpd2_ligand1', 'Cmpd3_ligand1', 'Cmpd3_ligand2', 'Cmpd4_ligand1', 'Cmpd4_ligand2', 'Cmpd5_ligand1', 'Cmpd5_ligand2', 'Cmpd6_ligand1', 'Cmpd6_ligand2', 'Cmpd7_ligand1', 'Cmpd8_ligand1', 'Cmpd9_ligand1']
    #liste_gal = ['3Q0V', '5NZ1']
    n=0
    for i in liste_gal:
        n += 1
        print('######################## Docking Galaxy ' + str(n) + ' sur ' + str(len(liste_gal)) + ' ----> ' + i)
        dock_galaxy(i)


if __name__ == '__main__':
    if terminal:
        if 'GalaxyD' not in listdir(path):
            os.system('mkdir GalaxyD')

    dockG()

    #conversion_mol2()
    # triage()
    # e = '5MXV_ligand1'
    # dock_galaxy(e)
