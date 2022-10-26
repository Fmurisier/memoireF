"""
Redock autodock FR
Editor : Murisier Frederic

Date : mars 2022

redocking for AutodockFR : docking each receptor with his ligand(s) and check if the docking went well

Python version : 3.8
"""
import datetime
import os
from os import listdir
from os.path import isfile, join

path = os.getcwd()
terminal = False
# variable environnement
PATH = '/home/rene/bin/ADFRsuite-1.0/bin/'


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
    fichier_log.write('Residus = ' + residus + '\n' * 2)
    return residus


def liste_file(pattern, file_path=''):
    """
    Return the list of file with the pattern given in argument present in the directory specified in argument
    the name in the list do not have the extension file tag
    :param pattern:
    :param file_path:
    :return:
    """
    list_file = []

    # list of file present in the directory specified in argument
    fichiers = [f for f in listdir(path + file_path) if isfile(join(path + file_path, f))]

    # split the name of each file by '_' to select only the file name with the pattern given in argument
    for ligne in fichiers:
        if pattern in ligne:
            name = ligne.split('_')
            list_file.append(name[0])
    list_file.sort()

    return list_file


def obtention_liste_pdbqt_H():
    """
    Return the list of the pdb file transformed in the current directory
    :return: list
    """
    list_pdbqt = []
    c = 0

    # list of the pdb file transformed in the current directory
    fichiers = [f for f in listdir(path) if isfile(join(path, f))]

    # for each file name split by '_' to find only the file with '_transformed.pdb' and then store the names in a list
    for ligne in fichiers:
        name = ligne.split('_')

        if name[-1] == 'H.pdbqt':
            list_pdbqt.append(name[0])
            c += 1

    list_pdbqt.sort()
    message = 'Number of pdb file in the directory : ' + str(c) + '\n'
    print(message)
    fichier_log.write(message)

    return list_pdbqt


def grids(fichier):
    """
    Create the target file (= grid or box) for the docking
    an exemple of commande line :
    /home/rene/bin/ADFRsuite-1.0/bin/agfr -b residues A:LEU87,LEU90,MET102,TRP103,ILE107,PHE110,PHE114,TRP138,MET142,
    TRP145,TYR148,THR149,VAL152,ASN176,ASN179,GLU180,LEU183,PHE184,TRP207 -r RECEPTEUR/PDBQT/1T56_recepteur_H.pdbqt -o
    1T56
    loop to minimise the faillure when creating the box
    :param fichier: name of the receptor file from which we need to create the grid
    :return:
    """
    boucle = True
    commande = PATH + 'agfr -b ' + RESIDUES + ' -r RECEPTEUR/PDBQT/' + fichier + '_recepteur_H.pdbqt -o target/' + \
               fichier

    while boucle:
        if terminal:
            os.system(commande)
        print(commande)
        fichier_log.write('\n\nCommande line :\n' + commande + '\n')
        boucle = check_error_grid(fichier + '.log')


def check_error_grid(file):
    """
    Check if there is any error during the creation of the box. If yes, error message with the information in a
    fileerror.txt file
    :param file: file to check
    :return:
    """
    restart = False
    fichier = open('target/' + file, 'r')
    fichierlog = open('fileERROR_Grid.txt', 'a')
    for l in fichier:
        if 'ERROR' in l:
            message = '\nWARNING PROBLEM agdr crashes....\n' + l + '\n'
            print(message)
            fichierlog.write(message)
            fichier_log.write(message)
            restart = True
    if not restart:
        message = '\n\n\nNo problem detected in the grid creation\n\n\n'
        print(message)
        fichier_log.write(message)
    fichierlog.close()
    return restart


def grid_for_all():
    """
    Create the grid for each structure calling the grid function for every file
    :return:
    """
    liste_grid = obtention_liste_pdbqt_H()
    for element in liste_grid:
        grids(element)


def dock(file):
    """
    Run the docking for the structure in argument
    example of commande:
    /home/rene/bin/ADFRsuite-1.0/bin/adfr -l LIGAND/PDBQT/1T56_ligand1_smile.pdbqt -t target/1T56.trg -r
    LIGAND/PDBQT/1T56_ligand1.pdbqt -p 150 -e 5000000 -n 100 -c 6 -O -o 1T56_result -J dock
    :param file:
    :return:
    """
    commande = PATH + 'adfr -l LIGAND/PDBQT/' + file + '_ligand1_smile.pdbqt -t target/' + file + \
               '.trg -r LIGAND/PDBQT/' + file + '_ligand1.pdbqt -p 150 -e 5000000 -n 100 -c 6 -O -o AUTODOCKFR/' + \
               file + '_result -J dock'
    if terminal:
        os.system(commande)
    print(commande)
    fichier_log.write(commande + '\n')


def grid_dock_liste(start=0, end=82):
    """
    Create the grid and docking for each structur in the liste of all receptor from indice start to end
    the aim of this parameter if to allows the running of some part of the list
    :param start:
    :param end:
    :return:
    """
    liste = liste_file('H.pdbqt', '/RECEPTEUR/PDBQT/')
    if end == 82:
        end = len(liste)
    c = 0
    sep = '#' * 110 + '\n'
    for element in range(start, end):
        if terminal:
            c += 1
            info_message = sep + str(c) + '-' * 18 + str(end - start) + ' ' * 4 + liste[element] + '\n' + sep
            print(info_message)
            fichier_log.write(info_message)
            grids(liste[element])
            dock(liste[element])

        print(str(element) + ' --> ' + liste[element])
        fichier_log.write(str(element) + ' --> ' + liste[element] + '\n')


def check_summary():
    """
    Check the length of the file in the result directory. If the file is too small (< 3000 octet) the docking failled.
    :return:
    """
    error = open('taille_error.txt', 'w')
    if terminal:
        os.system('ls -l AUTODOCKFR/*summary* > taille.txt')

    taille = open('taille.txt', 'r')
    fichier_log.write('\n\n\nChecking if the docking failled...\n')
    for ligne in taille:
        sep = ligne.split(' ')
        d = 0
        if sep[4] == '':
            d = 1
        if int(sep[4 + d]) < 3000:
            name = sep[-1].split('_')
            error.write(name[0] + '-' * 8 + sep[4 + d] + '-' * 8 + sep[-1])
            fichier_log.write(name[0] + '-' * 8 + sep[4 + d] + '-' * 8 + sep[-1])
    taille.close()
    fichier_log.write('\n\n...Checking Done\n\n')


def ttt_resultat(name, file_result):
    """
    file + _result_summary.dlg
    :return:
    """

    file = open('AUTODOCKFR/' + name + '_result_summary.dlg', 'r')

    p = False
    t = False
    table = []
    best_score = 100
    top = ''
    for ligne in file:
        if 'runs failed' in ligne:
            p = True
        if 'Clustering information' in ligne:
            p = False
        if p:
            if '-----' in ligne:
                t = True
            if t:
                table.append(ligne.split())
    table.pop(0)
    for e in table:
        if e:
            if float(e[3]) < best_score:
                best_score = float(e[3])
                top = e[0]
    file_result.write(name + ' top : ' + top + '  score :  ' + str(best_score) + ' \n')
    print(name + '---> top : ' + top + '  score :  ' + str(best_score))


def ttt_resultat_all():
    """
        ### one by one
        n = 8
    liste = liste_file('summary', '/AUTODOCKFR/')
    resultat_adfr = open('resultatADFR.txt', 'w')
    print(liste[n])
    ttt_resultat(liste[n])
    resultat_adfr.close()
    :return:
    """
    lise_error = ['3Q0V', '5MYM', '5NZ1']

    liste = liste_file('summary', '/AUTODOCKFR/')
    resultat_adfr = open('resultatADFR.txt', 'w')

    for e in liste:
        if e not in lise_error:
            ttt_resultat(e, resultat_adfr)
    resultat_adfr.close()


def check_result():
    """
    classement des resultats par top
    :return:
    """
    result = open('resultatADFR.txt', 'r')

    top1 = []
    top3 = []
    top10 = []
    top_else = []
    fail = []

    for e in result:
        cut = e.split(' ')
        if float(cut[-2]) > 2:
            fail.append(e)
        elif float(cut[3]) == 1:
            top1.append(e)
        elif float(cut[3]) <= 3:
            top3.append(e)
        elif float(cut[3]) <= 10:
            top10.append(e)
        else:
            top_else.append(e)

    names = ['top1', 'top3', 'top10', 'else', 'fail']  # nom des barres

    values = [len(top1), len(top3), len(top10), len(top_else), len(fail)]

    plt.bar(names, values)
    plt.show()  # Tracer

    data1 = len(top1)
    data2 = len(top3)
    data3 = len(top10)
    data4 = len(top_else)
    data5 = len(fail)

    year = ["AutodockFR"]

    plt.figure(figsize=(10, 10))
    plt.bar(year, data4, color="red", label="else --> " + str(data4))
    plt.bar(year, data3, color="orange", bottom=np.array(data4), label="top10 --> " + str(data3))
    plt.bar(year, data2, color="yellow", bottom=np.array(data4) + np.array(data3), label="top3 --> " + str(data2))
    plt.bar(year, data1, color="green", bottom=np.array(data4) + np.array(data3) + np.array(data2),
            label="top1 --> " + str(data1))

    plt.legend(loc="lower left", bbox_to_anchor=(0.8, 1.0))
    plt.show()


if __name__ == '__main__':
    fichier_log = open('log_ADFR_reDock.txt', 'a')
    fichier_log.write(str(datetime.date.today()) + ' -' * 40)

    RESIDUES = lecture_residu()
    if terminal:
        ldoc = ['target', 'AUTODOCKFR']
        for dossier in ldoc:
            if dossier not in listdir(path):
                os.system('mkdir ' + dossier)

    liste = ['5MYM']
    # grid_for_all()
    # grids('model_alphaFold')
    # dock('1T56')
    for e in liste:
        grids(e)
        dock(e)
    # grid_dock_liste(0, 25)
    check_summary()
    # ttt_resultat_all()
    # check_result()
