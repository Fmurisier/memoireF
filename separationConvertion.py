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
import datetime

path = os.getcwd()
terminal = False
# variable environnement
PATH = '/home/rene/bin/ADFRsuite-1.0/bin/'
REF = '1T56'
LIGABS = []
liste_solute = ['GOL', 'HOH', 'MES', 'TRS', 'MG', 'SO4', 'CL', 'PO4', 'AZI', 'ACT', 'K', 'MN', 'SCN', 'ZN', 'BR']


def separation():
    """
    List the transformed file present in the repertory then for each transformed file split it in receptor and ligand(s)
    files
    :return:
    """
    transf_liste = [f for f in listdir(path + '/transformed') if isfile(join(path + '/transformed', f))]

    for element in transf_liste:
        e = element.split('_')
        if os.path.exists(e[0] + '_ligand1.pdb'):
            # delete previous lig file
            os.system('rm ' + e[0] + '_lig*')
        # give only the name of the file in argument, removing the '_transformed.pdb'
        write_sep_file(e[0])


def write_sep_file(e):
    """
    Take in argument the file transformed aimed to be split in receptor and ligand(s) files.
    :param e:
    :return:
    """
    fichier = open('transformed/' + e + '_transformed.pdb', 'r')
    # fichier = open('../Donnee_memoire/transformed_box/' + e + '_transformed.pdb', 'r')
    fichier_recepteur = open(e + '_recepteur.pdb', 'w')
    lines = fichier.readlines()
    lineone = lines[0].split()
    num_lig = lineone[5]
    lc = 0  # set the ligand number
    chain = 'A'  # set the chain
    first = True
    ligand_ref = ''
    error = False
    errorlig = False
    res_alt = []
    lig_alt = []
    num_res = ecriture_box(True)[0]

    for i in lines:
        line = i.split()
        # avoid error when the line contain only one column
        if len(line) > 4:
            # check if the third column didn't have a name tag too long merging with the next column
            if len(line[2]) > 3:
                d = -1
            else:
                d = 0
            # taking only the chain A
            if line[0] == 'ATOM' and line[4 + d] == chain:
                fichier_recepteur.write(i)
                if float(i.split()[-3]) < 1:
                    if i.split()[-7] in num_res and i.split()[-7] not in res_alt:
                        res_alt.append(i.split()[-7])
                        error = True
            # keep only the ligands, which is not in the list of solute and only from the chain A
            elif line[0] == 'HETATM' and line[3 + d] not in liste_solute and line[4 + d] == 'A':
                if d == -1:
                    ligand = line[2][3:]
                else:
                    ligand = line[3]
                if first:
                    ligand_ref = ligand
                    first = False
                if ligand == ligand_ref:
                    # incremente the name of the file in case of multiple ligands
                    if num_lig != line[5 + d]:
                        chain = line[4 + d]
                        num_lig = line[5 + d]
                        lc += 1
                    elif chain != line[4 + d]:
                        chain = line[4 + d]
                        lc += 1
                    if float(i.split()[-3]) < 1 and i.split()[-7] not in lig_alt:
                        lig_alt.append(i.split()[-7])
                        errorlig = True

                    fichier_ligand = open(e + '_ligand' + str(lc) + '.pdb', 'a')
                    fichier_ligand.write(i)
                    fichier_ligand.close()
    m = ''
    # safety warning if there is no ligand detected for a structure
    if first:
        LIGABS.append(e)
        m = 'No ligand detected ! ' + e
        print(m)
        fichier_log.write(m + '\n')
    if errorlig:
        m = 'alternatif LIGAND conformation for : ' + e
        print(m)
        fichier_log.write(m + '\n')
        # print('residu from the box having alternative conformation : ')
        print(lig_alt)
    if error:
        m = 'alternatif RECEPTOR conformation for : ' + e
        print(m)
        fichier_log.write(m + '\n')
        # print('residu from the box having alternative conformation : ')
        print(res_alt)
    fichier_recepteur.close()
    fichier.close()
    print(m)
    fichier_log.write('END SEPARATION' + '\n' * 2)


def lecture_residu():
    """
    extraction of the name and number of the different residus of the box
    :return:
    """
    resi_file = open('residus.txt', 'r').readlines()
    residus = ''
    if len(resi_file) > 1:
        for e in resi_file:
            residus += e[:-1]
    else:
        residus += resi_file[0]

    if 'residues' not in residus:
        residus = 'residues A:' + residus
    if residus[-1] == '\n':
        residus = residus[:-1]
    print(residus)
    fichier_log.write('Residus = ' + residus + '\n' *2)
    return residus


def check_error_grid(file):
    """
    Check if there is any error during the creation of the box. If yes, error message with the information in a
    fileerror.txt file
    :param file: fichier a verifier
    :return:
    """
    restart = False
    fichier = open(file, 'r')
    for l in fichier:
        if 'ERROR' in l:
            print('\nWARNING PROBLEM agdr crashes....\n' + l)
            fichier_log.write('\nWARNING PROBLEM agdr crashes....\n' + l + '\n')
            restart = True
    if not restart:
        print('\nNo problem detected in the grid creation\n')
        fichier_log.write('\nNo problem detected in the grid creation\n' + '\n')
    return restart


def grids(fichier):
    """
    Create the target file (= grid or box) for the docking
    an exemple of commande line :
    /home/rene/bin/ADFRsuite-1.0/bin/agfr -b residues A:LEU87,LEU90,MET102,TRP103,ILE107,PHE110,PHE114,TRP138,MET142,
    TRP145,TYR148,THR149,VAL152,ASN176,ASN179,GLU180,LEU183,PHE184,TRP207 -r RECEPTEUR/PDBQT/1T56_recepteur_H.pdbqt -o
    1T56
    :param fichier: name of the structure for whitch the grid need to be created
    :return:
    """
    RESIDUES = lecture_residu()
    commande = PATH + 'agfr -b ' + RESIDUES + ' -r ' + fichier + '.pdb -o ' + fichier

    if terminal:
        os.system(commande)
    else:
        print(commande)
    fichier_log.write('\n\nCommande line :\n' + commande + '\n')


def verification_ligand_box():
    """
    The function verify if the ligand are in the box, if not print an error message
    :return:
    """
    grids(REF)
    fichier_box = open('1T56.log', 'r').readlines()
    center = ''
    length = ''
    for l in fichier_box:
        l = l[:-1]
        if 'Box ' in l:
            e = l.split(' ')
            for i in range(e.count('')):
                e.remove('')
            if e[1] == 'center:':
                center = e
            elif e[1] == 'length:':
                length = e
    ligands = liste_file_complet('ligand')
    print(center)
    print(length)
    fichier_log.write('center = ' + center + '\n')
    fichier_log.write('length = ' + length + '\n')

    xmin = float(center[2]) - float(length[2]) / 2
    xmax = float(center[2]) + float(length[2]) / 2
    ymin = float(center[3]) - float(length[3]) / 2
    ymax = float(center[3]) + float(length[3]) / 2
    zmin = float(center[4]) - float(length[4]) / 2
    zmax = float(center[4]) + float(length[4]) / 2
    print('x ', xmin, xmax)
    print('y ', ymin, ymax)
    print('z ', zmin, zmax)
    fichier_log.write('x ' + str(xmin) + '---' + str(xmax) + '\n')
    fichier_log.write('y ' + str(ymin) + '---' + str(ymax) + '\n')
    fichier_log.write('z ' + str(zmin) + '---' + str(zmax) + '\n')
    for i in ligands:
        # lig = open('../Donnee_memoire/ligand/' + i).readlines()
        lig = open(i).readlines()
        e = lig[0].split(' ')
        for h in range(e.count('')):
            e.remove('')
        e = e[:-4]
        z = float(e.pop())
        y = float(e.pop())
        x = float(e.pop())
        in_box = False
        if not xmin < x < xmax:
            in_box = True
        if not ymin < y < ymax:
            in_box = True
        if not zmin < z < zmax:
            in_box = True
        if in_box:
            print('WARNING ! ligand ' + i + ' is not in the box')
            fichier_log.write('\n\nWARNING !\n ligand ' + i + ' is not in the box\n\n')
    print('Verification done')
    fichier_log.write('Verification done\n')
    if LIGABS:
        print('ligand absent from Chain A: ' + '  '.join(LIGABS))
        fichier_log.write('WARNING !\nligand absent from Chain A: ' + '  '.join(LIGABS) + '\n')


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
        print('#' * 100 + '\n' + str(c) + '-' * 17 + str(compteur) + ' ' * 4 + ligne + '\n' + '#' * 100 + '\n')

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
    # fichier = ['3Q0V', '5MYM', '5NZ1']
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


def liste_file_complet(patern, file_path=''):
    """
    MODIFIER POUR QU'IL DONNE DIRECT LE NOM COMPLET
    :param patern:
    :param file_path:
    :return:
    """
    list_file = []

    # fait la liste des fichier dans le repertoire ou se trouve script.py
    fichiers = [f for f in listdir(path + file_path) if isfile(join(path + file_path, f))]

    # pour chaque nom de fichier on separe le nom par '_' pour chercher uniquement les fichiers comportant le prefix
    # '_transformed.pdb' et on stocke tous les noms de fichier dans une liste
    for ligne in fichiers:
        if patern in ligne:
            list_file.append(ligne)
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


def ecriture_box(liste=False):
    dossier = 'residus.txt'
    if dossier in listdir(path):
        resi_file = open('residus.txt', 'r').readlines()
        resi_liste = []
        resi_code3 = []
        for e in resi_file:
            if '\n' in e:
                e = e[:-1]
                new_line = e.split(',')
                for e in new_line:
                    if ':' not in e and e != '':
                        resi_code3.append(e)
                        resi_liste.append(e[3:])
                    elif e != '':
                        resi_code3.append(e.split(':')[-1])
                        resi_liste.append(e.split(':')[-1][3:])
        if not liste:
            print('The box is made of ' + str(len(resi_liste)) + \
                  ' residus, if it is not the case then please check that the file is filled correctly (' \
                  ' exemple : \'residues A:LEU87,LEU90,MET102,TRP103,ILE107,... \')')

        commande = '_poche, resi ' + '+'.join(resi_liste) + ' and model '
        if liste:
            commande = resi_liste
        return commande, resi_code3
    else:
        print('error file residus.txt don\'t exist ! ')
        return 'error', ''


if __name__ == '__main__':
    fichier_log = open('log_pipeline.txt', 'a')
    fichier_log.write(str(datetime.date.today()) + ' -' * 40)
    fichier_log.write('START SEPARATION')

    separation()

    verification_ligand_box()

    # conversion_recepteur_pdbqt()

    # conversion_ligand_smile_pdbqt()

    # mol2()

    fichier_log.close()

    '''
    listeX = ['1T56', '5NIZ', '5NIM']

    for e in listeX:
        write_sep_file(e)
    
    '''
