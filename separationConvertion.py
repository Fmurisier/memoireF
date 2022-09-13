"""
Separation ligand recepteur et conversions lig ref en pdbqt et les ligand pour docking en smile puis pdb puis pdbqt
ensuite conversion de chaque fichier en mol2

Editor : Murisier Frederic

Separation : Prend tous les fichiers qui ont ete transforme et va les separer en un fichier ligand et un fichier
recepteur
Conversion : converti les fichier en de pdb a pdbqt, de pdb a smile, de smile a pdb, de pdbqt a mol2
Date : 2022
cree une liste des fichiers qui ont ete transforme (= les fichiers qui ont l'extension '_transformed.pdb' et va ensuite
cree de nouveaux fichiers separant les proteines et les ligands en creant un nouveau fichier lorsqu'il y a plusieurs
ligands
Resultat : nouveaux fichiers separe : ligand(s) et recepteur pour chaque fichier

Python version : 3.8

"""

import os
from os import listdir
from os.path import isfile, join

path = os.getcwd()
terminal = False
# variable environnement
PATH = '/home/rene/bin/ADFRsuite-1.0/bin/'
liste_soluté = ['GOL', 'HOH', 'MES', 'TRS', 'MG', 'SO4', 'CL', 'PO4', 'AZI', 'ACT', 'K', 'MN', 'SCN', 'ZN', 'BR']


def obtention_liste_pdb_transformed():
    """
    va effectuer la liste des fichiers pdb transforme dans le repertoire transformed
    :return: list
    """
    list_pdb = []
    c = 0

    # fait la liste des fichier dans le repertoire ou se trouve script.py
    fichiers = [f for f in listdir(path + '/transformed') if isfile(join(path + '/transformed', f))]

    # pour chaque nom de fichier on separe le nom par '_' pour chercher uniquement les fichiers comportant le prefix
    # '_transformed.pdb' et on stocke tous les noms de fichier dans une liste
    for ligne in fichiers:
        name = ligne.split('_')

        if name[-1] == 'transformed.pdb':
            list_pdb.append(name[0])
            c += 1

    list_pdb.sort()
    print('Nombre de pdb file transforme dans le dossier : ', c)

    return list_pdb


def separation():
    """
    va appeler une fonction pour lister les fichiers transforme dans le repertoire et va prendre chaque fichier et
    appelera une seconde fonction qui separera le fichier d'origine en deux nouveaux fichier recepteur et ligand(s)
    :return:
    """
    # transf_liste = obtention_liste_pdb_transformed()
    transf_liste = [f for f in listdir(path + '/transformed') if isfile(join(path + '/transformed', f))]

    for element in transf_liste:
        e = element.split('_')
        if os.path.exists(e[0] + '_ligand1.pdb'):
            os.system('rm ' + e[0] + '_lig*')
        write_sep_file(e[0])


def write_sep_file(e):
    """
    recoit le nom d'un fichier qui va etre lu pour ensuite creer de nouveaux fichier : recepteur (contenant uniquement
    les informations du recepteur et ligand (contenant uniquement les informations du ligands), creera plusieurs fichier
    ligand si le fichier d'origine contient plusieurs ligand
    :param e:
    :return:
    """

    fichier = open('transformed/' + e + '_transformed.pdb', 'r')
    fichier_recepteur = open(e + '_recepteur.pdb', 'w')
    lines = fichier.readlines()
    lineone = lines[0].split()
    num_lig = lineone[5]
    lc = 0
    chain = lineone[4]
    first = True
    ligand_ref = ''

    for i in lines:
        line = i.split()
        # empeche de faire une erreur quand la ligne ne contiens qu'une seule colonne
        if len(line) > 4:
            # verifie si la 3eme colonne n'a pas un nom qui fait que la 3 eme et 4eme colonne sont accolee
            if len(line[2]) > 3:
                d = -1
            else:
                d = 0
            # prend uniquement la chaine A du recepteur
            if line[0] == 'ATOM' and line[4 + d] == chain:
                fichier_recepteur.write(i)
            # prend uniquement les ligands et ne prend pas les solutes
            elif line[0] == 'HETATM' and line[3 + d] not in liste_soluté:
                if d == -1:
                    ligand = line[2][3:]
                else:
                    ligand = line[3]
                if first:
                    ligand_ref = ligand
                    first = False
                if ligand == ligand_ref:
                    # incremente le nom du fichier pour le cas ou il y a deux ligands
                    if num_lig != line[5 + d]:
                        chain = line[4 + d]
                        num_lig = line[5 + d]
                        lc += 1
                    elif chain != line[4 + d]:
                        chain = line[4 + d]
                        lc += 1

                    fichier_ligand = open(e + '_ligand' + str(lc) + '.pdb', 'a')
                    fichier_ligand.write(i)
                    fichier_ligand.close()

    fichier_recepteur.close()
    fichier.close()


def conversion_ligand_smile_pdbqt():
    """
    /home/rene/bin/ADFRsuite-1.0/bin/prepare_ligand -l 1T56_ligand1.pdb
    obabel -ipdb 1T56_ligand1.pdb -osmi -O 1T56_ligand1.smi
    acedrg -i 1T56_ligand1.smi -o 1T56_ligand1_smile
    /home/rene/bin/ADFRsuite-1.0/bin/prepare_ligand -l 1T56_ligand1_smile.pdb
    :return:
    """
    c = 0
    fichier_log = open('log_smile', 'w')
    fichiers = liste_file('_ligand1')
    fichier = ['3Q0V', '5MYM', '5NZ1']
    compteur = len(fichiers)

    for ligne in fichiers:
        c += 1
        print('####################################################################################################\n' +
              str(c) + '------------------>' + str(compteur) + '    ' + ligne +
              '\n###################################################################################################\n')

        # transforme le ligand de ref en pdbqt ################# pdb to pdbqt ################
        conversion_pdbqt(ligne + '_ligand1.pdb')

        # commande qui transforme un pdb en smi ################# pdb to smi ################
        if terminal:
            os.system('obabel -ipdb ' + ligne + '_ligand1.pdb -osmi -O ' + ligne + '_ligand1.smi')
        print('obabel -ipdb ' + ligne + '_ligand1.pdb -osmi -O ' + ligne + '_ligand1.smi')

        # commande qui transforme un smi en pdb ################### smi to pdb ##################
        if terminal:
            os.system('acedrg -i ' + ligne + '_ligand1.smi -o ' + ligne + '_ligand1_smile')
            fichiers_smile_pdb = liste_file('_ligand1_smile')
            if ligne not in fichiers_smile_pdb:
                fichier_log.write('\nerror smile : ' + ligne)
        print('acedrg -i ' + ligne + '_ligand1.smi -o ' + ligne + '_ligand1_smile')

        # commande qui transforme un pdb en pdbqt ################# pdb to pdbqt ################
        conversion_pdbqt(ligne + '_ligand1_smile.pdb')

    fichier_log.close()
    if terminal:
        if 'LIGAND' not in listdir(path):
            os.system('mkdir LIGAND')
        if 'SMILE' not in listdir(path):
            os.system('mkdir SMILE')
        if 'smi' not in listdir(path):
            os.system('mkdir smi')
        if 'smileLog' not in listdir(path):
            os.system('mkdir smileLog')
        os.system('mv *_ligand1_smile.pdb SMILE/')
        os.system('mv *.smi smi/')
        os.system('mv *_TMP smileLog/')
        os.system('mv *cif smileLog/')

        if 'PDBQT' not in listdir(path + '/LIGAND'):
            os.system('mkdir LIGAND/PDBQT')
        os.system('mv *.pdbqt LIGAND/PDBQT')
        os.system('mv *_ligand1.pdb LIGAND')
        if 'ligand2_3' not in listdir(path + '/LIGAND'):
            os.system('mkdir LIGAND/ligand2_3')
        os.system('mv *_ligand* LIGAND/ligand2_3')


def conversion_ligand_smile_pdbqt_liste(files):
    """
    /home/rene/bin/ADFRsuite-1.0/bin/prepare_ligand -l 1T56_ligand1.pdb
    obabel -ipdb 1T56_ligand1.pdb -osmi -O 1T56_ligand1.smi
    acedrg -i 1T56_ligand1.smi -o 1T56_ligand1_smile
    /home/rene/bin/ADFRsuite-1.0/bin/prepare_ligand -l 1T56_ligand1_smile.pdb
    :return:
    """
    c = 0
    fichier_log = open('log_smile', 'w')

    compteur = len(files)

    for ligne in files:
        c += 1
        print('####################################################################################################\n' +
              str(c) + '------------------>' + str(compteur) + '    ' + ligne +
              '\n###################################################################################################\n')

        # transforme le ligand de ref en pdbqt ################# pdb to pdbqt ################
        conversion_pdbqt(ligne + '_ligand1.pdb')

        # commande qui transforme un smi en pdb ################### smi to pdb ##################
        if terminal:
            os.system('acedrg -i ' + ligne + '_ligand1.smi -o ' + ligne + '_ligand1_smile')
            fichiers_smile_pdb = liste_file('_ligand1_smile')
            if ligne not in fichiers_smile_pdb:
                fichier_log.write('\nerror smile : ' + ligne)
        print('acedrg -i ' + ligne + '_ligand1.smi -o ' + ligne + '_ligand1_smile')

        # commande qui transforme un pdb en pdbqt ################# pdb to pdbqt ################
        conversion_pdbqt(ligne + '_ligand1_smile.pdb')

    fichier_log.close()
    if terminal:
        os.system('mv *_ligand1_smile.pdb SMILE/')
        os.system('mv *.smi SMILE/smi/')
        os.system('mv *_TMP SMILE/smileLog/')
        os.system('mv *cif SMILE/smileLog/')


def conversion_pdbqt(nomfichier, lig=True):
    """
    convertion d'un fichier pdb fournis en argument en pdbqt
    /home/rene/bin/ADFRsuite-1.0/bin/prepare_ligand -l _ligand1.pdb
    :param lig: boolean
    :param nomfichier:
    :return:
    """
    # commande qui transforme un pdb en pdbqt
    # os.system('/home/rene/bin/ADFRsuite-1.0/bin/prepare_ligand -l ' + nomfichier)
    if lig:
        if terminal:
            os.system(PATH + 'prepare_ligand -l ' + nomfichier)
        print(PATH + 'prepare_ligand -l ' + nomfichier)
    else:
        commande1 = 'python3 ../../../home/rene/bin/pdb2pqr/pdb2pqr --titration-state-method propka --drop-water ' \
                    '--pdb-output ' + nomfichier + '_H.pdb --with-ph 7.4 ' + nomfichier + '.pdb ' + nomfichier + '.pqr'
        if terminal:
            os.system(commande1)
            os.system(PATH + 'prepare_receptor -r ' + nomfichier + '_H.pdb')
        print(commande1)
        print(PATH + 'prepare_receptor -r ' + nomfichier + '_H.pdb')


def conversion_recepteur_pdbqt():
    c = 0
    fichiers = [f for f in listdir(path) if isfile(join(path, f))]

    fichiers = liste_file('_recepteur')
    fichier = ['3Q0V', '5MYM', '5NZ1']
    compteur = len(fichiers)

    for ligne in fichiers:
        c += 1
        print('####################################################################################################\n' +
              str(c) + '------------------>' + str(compteur) + '           ' + ligne +
              '\n###################################################################################################\n')
        conversion_pdbqt(ligne + '_recepteur', False)

    if terminal:
        dossier = 'RECEPTEUR'
        if dossier not in listdir(path):
            os.system('mkdir ' + dossier)
        os.system('mv *_recepteur.pdb RECEPTEUR')
        dossier = 'PDBQT'
        if dossier not in listdir(path + '/RECEPTEUR'):
            os.system('mkdir /RECEPTEUR/' + dossier)
        os.system('mv *_H.pdbqt RECEPTEUR/PDBQT')
        dossier = 'pdbqtLOG'
        if dossier not in listdir(path + '/RECEPTEUR'):
            os.system('mkdir /RECEPTEUR/' + dossier)
        os.system('mv *_recepteur* RECEPTEUR/pdbqtLOG')


def liste_file(patern, file_path=''):
    list_file = []

    # fait la liste des fichier dans le repertoire ou se trouve script.py
    fichiers = [f for f in listdir(path + file_path) if isfile(join(path + file_path, f))]

    # pour chaque nom de fichier on separe le nom par '_' pour chercher uniquement les fichiers comportant le prefix
    # '_transformed.pdb' et on stocke tous les noms de fichier dans une liste
    for ligne in fichiers:
        if patern in ligne:
            name = ligne.split('_')
            list_file.append(name[0])
    list_file.sort()

    return list_file


def mol2():
    '''

    '/usr/bin/obabel -ipdbqt fichier_ligand1.pdbqt -O XX.mol2'
    :return:
    '''
    if terminal:
        dossier = 'MOL2'
        if dossier not in listdir(path + '/LIGAND'):
            os.system('mkdir LIGAND/' + dossier)
    l = os.listdir('LIGAND/PDBQT')
    l2 = os.listdir('LIGAND/ligand2_3')
    l = ['3Q0V_ligand1.pdb', '5MYM_ligand1.pdb', '5NZ1_ligand1.pdb']
    l2 = ['3Q0V_ligand2.pdb', '5MYM_ligand2.pdb', '5NZ1_ligand2.pdb', '5NZ1_ligand3.pdb']

    for h in l:
        if '_smile.pdbqt' not in h:
            if terminal:
                os.system('/usr/bin/obabel -ipdbqt LIGAND/PDBQT/' + h + ' -O LIGAND/MOL2/' + h[:-3] + 'mol2')
            print('/usr/bin/obabel -ipdbqt LIGAND/PDBQT/' + h + ' -O LIGAND/MOL2/' + h[:-3] + 'mol2')

    for h in l2:
        if '.pdb' in h:
            if terminal:
                os.system('/usr/bin/obabel -ipdb LIGAND/ligand2_3/' + h + ' -O LIGAND/ligand2_3/' + h[:-3] + 'pdbqt')
                os.system('/usr/bin/obabel -ipdbqt LIGAND/ligand2_3/' + h[:-3] + 'pdbqt' + ' -O LIGAND/MOL2/' + h[:-3]
                          + 'mol2')
            print('/usr/bin/obabel -ipdb LIGAND/ligand2_3/' + h + ' -O LIGAND/ligand2_3/' + h[:-3] + 'pdbqt')
            print('/usr/bin/obabel -ipdbqt LIGAND/ligand2_3/' + h[:-3] + 'pdbqt' + ' -O LIGAND/MOL2/' + h[:-3] + 'mol2')


if __name__ == '__main__':

    separation()

    conversion_recepteur_pdbqt()

    conversion_ligand_smile_pdbqt()

    mol2()
