"""
redock GalaxyDock3
Editor : Murisier Frederic
Date : 2022
redocking for GalaxyDock3 : docking each receptor with his ligand(s) and check if the docking went well

Python version : 3.8


"""
import datetime
import os
from os import listdir
from os.path import isfile, join

path = os.getcwd()
terminal = False
# variable environnement
PATH = '/home/rene/bin/GalaxyDock3/'


def liste_file(pattern, file_path=''):
    """
    Return the list of file with the pattern given in argument present in the directory specified in argument
    the name in the list do not have the extension file tag
    :param pattern:
    :param file_path: path of the directory, if empty then current directory
    :return:
    """
    list_file = []

    # Generate the list of all file in the directory where the script.py is located
    fichiers = [f for f in listdir(path + file_path) if isfile(join(path + file_path, f))]

    # For loop to split each file name by '.' to keep only the name of the files in a list
    for ligne in fichiers:
        if pattern in ligne:
            name = ligne.split('.')
            list_file.append(name[0])
    list_file.sort()

    return list_file


def dock_galaxy(e):
    """
    Dock the structure given in argument automaticaly searching for the box information
    :param e: 
    :return: 
    """

    # Creation of a directory for each structure to store the specific result
    if terminal:
        if e not in listdir(path + '/GalaxyD'):
            os.system('mkdir GalaxyD/' + e)

    # Open the box file
    file = open('VINA/BoxTxt/' + e.split('_')[0] + '_box.txt', 'r').readlines()
    # Extracting each coordinate of the box
    cx = file[0].split()[-1]
    cy = file[1].split()[-1]
    cz = file[2].split()[-1]
    sx = file[3].split()[-1]
    sy = file[4].split()[-1]
    sz = file[5].split()[-1]

    # Commande line for the Docking with Galaxy dock
    commande = '/home/rene/bin/GalaxyDock3/script/run_GalaxyDock3.py -d /home/rene/bin/GalaxyDock3 -p RECEPTEUR/' + \
               e.split('_')[0] + '_recepteur.pdb -l LIGAND/MOL2/' + e + '.mol2 -x ' + cx + ' -y ' + cy + ' -z ' + \
               cz + ' -size_x ' + sx + ' -size_y ' + sy + ' -size_z ' + sz + ' --n_proc 1'

    if terminal:
        os.system(commande)

    print('docking--->' + commande)
    fichier_log.write('docking--->' + commande + '\n')

    # Using a Python script to split the result
    commande = 'python2 /home/rene/bin/GalaxyDock3/script/split_mol2.py GD3_cl.mol2 GalaxyD/' + e + '/' + e
    if terminal:
        os.system(commande)
    print('split---->' + commande)
    fichier_log.write('split--->' + commande + '\n')

    # Openning the result file, suppression of the three last colunm, adding a new colunm (at the end) 'rmsd'
    result = open('GD3_cl.E.info', 'r').readlines()
    resultfinal = open('GalaxyD/' + e + '_result_galaxy.txt', 'w')
    resultfinal.write(result[0] + result[1][:-15] + 'RMSD\n' + result[0])

    # Calcul the RMSD
    liste_galaxy = liste_file(e, '/GalaxyD/' + e)
    nbr_model = len(liste_galaxy)
    compteur = 0
    # For loop to generate the cmd line
    for model in liste_galaxy:
        commande = '/home/rene/bin/DockRMSD/DockRMSD LIGAND/MOL2/' + e + '.mol2 GalaxyD/' + e + '/' + model + \
                   '.mol2 > rescore_result.txt'

        if terminal:
            os.system(commande)
        print('rmsd ---> ' + commande)
        fichier_log.write('rmsd--->' + commande + '\n')

        # RMSD extraction
        rescore_lignes = open('rescore_result.txt', 'r').readlines()

        # Initialisation at 'NA' for the RMSD value in case the calcul of the rmsd went wrong
        dockrmsd = 'NA'
        for i in rescore_lignes:
            if 'Calculated Docking RMSD:' in i:
                dockrmsd = i.split()[3]

        # Writing the RMSD in the final file
        resultfinal.write(result[3 + compteur][:-27] + dockrmsd + '\n')
        compteur += 1

    resultfinal.close()


def triage():
    """
    Create a new file result_galaxy_tried.txt with the sorted result of galaxyDock
    :return: 
    """
    l = liste_file('result_galaxy.txt', '/GalaxyD')
    for s in l:
        print(s)
        fichier_log.write(s)
        file = result = open('GalaxyD/' + s + '.txt', 'r').readlines()
        c = 0
        for i in file:
            c += 1
        print(c)
        fichier_log.write(str(c))
        co = 'head -3 GalaxyD/' + s + '.txt > GalaxyD/' + s + '_tried.txt'
        commande = 'tail -' + str(c - 3) + ' GalaxyD/' + s + '.txt | sort -n -k 2 >> GalaxyD/' + s + '_tried.txt'

        if terminal:
            os.system(co)
            os.system(commande)

        print(co)
        print(commande)


def reconstruction(e, model):
    """
    Function for the reconstruction of the complexes
    :param e: 
    :param model: 
    :return: 
    """
    commande = 'obabel -imol2 GalaxyD/' + e + '/' + model + ' -O fichier_temporaire.pdb'
    if terminal:
        os.system(commande)
    print(commande)
    fichier_log.write(commande)


def conversion_mol2():
    """
    Convert all the ligand_smile file in mol
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
        fichier_log.write(commande)
        print(commande)


def dockG():
    """
    Dock the structure present in the /LIGAND/MOL2 repertory
    :return: 
    """
    liste_gal = liste_file('.', '/LIGAND/MOL2')
    #liste_gal = ['5NZ1_ligand1', '5NZ1_ligand2', '5NZ1_ligand3', '6HNX_ligand1', '6HNZ_ligand1', '6HO0_ligand1', '6HO1_ligand1', '6HO2_ligand1', '6HO3_ligand1', '6HO4_ligand1', '6HO5_ligand1', '6HO6_ligand1', '6HO7_ligand1', '6HO8_ligand1', '6HO9_ligand1', '6HOA_ligand1', '6HOB_ligand1', '6HOC_ligand1', '6HOD_ligand1', '6HOE_ligand1', '6HOF_ligand1', 'BDM15048_ligand1', 'BDM15048_ligand2', 'Cmpd1_ligand1', 'Cmpd1_ligand2', 'Cmpd2_ligand1', 'Cmpd3_ligand1', 'Cmpd3_ligand2', 'Cmpd4_ligand1', 'Cmpd4_ligand2', 'Cmpd5_ligand1', 'Cmpd5_ligand2', 'Cmpd6_ligand1', 'Cmpd6_ligand2', 'Cmpd7_ligand1', 'Cmpd8_ligand1', 'Cmpd9_ligand1']
    #liste_gal = ['3Q0V', '5NZ1']
    n = 0
    for i in liste_gal:
        n += 1

        message = "#" * 24 + " Docking Galaxy " + str(n) + ' on ' + str(len(liste_gal)) + ' ----> ' + str(i) + '\n'
        print(message)
        fichier_log.write(message)
        print()
        dock_galaxy(i)


if __name__ == '__main__':
    fichier_log = open('log_VINA_reDock.txt', 'a')
    fichier_log.write(str(datetime.date.today()) + ' -' * 40)

    if terminal:
        if 'GalaxyD' not in listdir(path):
            os.system('mkdir GalaxyD')

    dockG()

"""    
    message = '' + '\n'
    print(message)
    fichier_log.write(message)"""
    #conversion_mol2()
    # triage()
    # e = '5MXV_ligand1'
    # dock_galaxy(e)
