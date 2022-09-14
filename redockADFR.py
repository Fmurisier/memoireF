"""
Redock autodock FR
Editor : Murisier Frederic

Date : 2022

effectue le test 0 pour AutodockFR, docking de chaque recepteur avec son / ses ligands et evaluation si le docking a
bien fonctionne

Python version : 3.8
"""
import os
from os import listdir
from os.path import isfile, join

path = os.getcwd()
terminal = False
# variable environnement
PATH = '/home/rene/bin/ADFRsuite-1.0/bin/'



def lecture_residu():
    resi_file = open('residus.txt', 'r').readlines()
    residus = ''
    if len(resi_file) > 1:
        for e in resi_file:
            residus += e[:-1]
    else:
        residus += resi_file[0]

    if 'residues' not in residus:
        residus = 'residues A:' + residus
    print(residus)
    return residus


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
            name = ligne.split('_')
            list_file.append(name[0])
    list_file.sort()

    return list_file


def obtention_liste_pdbqt_H():
    """
    va effectuer la liste des fichiers pdb transforme dans le repertoire dans lequel le script.py se trouve
    :return: list
    """
    list_pdbqt = []
    c = 0

    # fait la liste des fichier dans le repertoire ou se trouve script.py
    fichiers = [f for f in listdir(path) if isfile(join(path, f))]

    # pour chaque nom de fichier on separe le nom par '_' pour chercher uniquement les fichiers comportant le prefix
    # '_transformed.pdb' et on stocke tous les noms de fichier dans une liste
    for ligne in fichiers:
        name = ligne.split('_')

        if name[-1] == 'H.pdbqt':
            list_pdbqt.append(name[0])
            c += 1

    list_pdbqt.sort()
    print('Nombre de pdb file transforme dans le dossier : ', c)

    return list_pdbqt


def grids(fichier):
    """
    effectue la creation du fichier target = la grille/ box pour le docking a partir du fichier recepteur
    un exemple de la commande qui doit etre effectue
    /home/rene/bin/ADFRsuite-1.0/bin/agfr -b residues A:LEU87,LEU90,MET102,TRP103,ILE107,PHE110,PHE114,TRP138,MET142,
    TRP145,TYR148,THR149,VAL152,ASN176,ASN179,GLU180,LEU183,PHE184,TRP207 -r RECEPTEUR/PDBQT/1T56_recepteur_H.pdbqt -o
    1T56
    fonction qui contien une boucle pour essayer de minimiser les echecs de creation de grille, en cas d'echec de
    creation de grille relancement de la ligne de commande pour effectuer de nouveau la creation de la grille
    :param fichier: nom de la structure recepteur a partir de laquelle il faut cree la grille
    :return:
    """
    boucle = True
    while boucle:
        if terminal:
            os.system(PATH + 'agfr -b ' + RESIDUES + ' -r RECEPTEUR/PDBQT/' + fichier + '_recepteur_H.pdbqt -o target/'
                      + fichier)
        else:
            print(PATH + 'agfr -b ' + RESIDUES + ' -r RECEPTEUR/PDBQT/' + fichier + '_recepteur_H.pdbqt -o target/' +
                  fichier)
        boucle = check_error_grid(fichier + '.log')


def check_error_grid(file):
    """
    verifie s'il y a eu des erreurs lors de la creation des gilles, si oui affichage d'un message d'erreur sur le
    terminal indiquant que agfr a crashe et sur quel structure. enregistrement egalement de ce message d'erreur dans un
    fichier fileerror.txt pour une consultation rapide des erreurs potentielles
    :param file: fichier a verifier
    :return:
    """
    restart = False
    fichier = open('target/' + file, 'r')
    fichierlog = open('fileERROR.txt', 'a')
    for l in fichier:
        if 'ERROR' in l:
            print('\nWARNING PROBLEM agdr crashes....\n')
            print(l)
            fichierlog.write('\nWARNING PROBLEM agdr crashes....\n' + l)
            restart = True
    if not restart:
        print('\nNo problem detected in the grid creation\n')
    fichierlog.close()
    return restart


def grid_for_all():
    """
    effectue la creation de grille pour le docking pour toutes les structures
    :return:
    """
    liste = obtention_liste_pdbqt_H()
    for e in liste:
        grids(e)


def dock(file):
    """
    effectue le docking pour la structure donnee en argument.
    exemple de commande:
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
    else:
        print(commande)


def grid_dock_liste(start=0,end=82):
    """
    effectue la creation de la grille et du docking pour chaque structure de la liste d'indice arg1 jusqu'a l'indice
    arg2

    0,20
    20,40
    40,60
    60,82
    :param start:
    :param end:
    :return:
    """
    liste = liste_file('H.pdbqt', '/RECEPTEUR/PDBQT/')
    c=0
    for element in range(start, end):
        if terminal:
            c += 1
            print(
                '##################################################################################################\n' +
                str(c) + '------------------>' + str(end - start) + '    ' + liste[element] +
                '\n#################################################################################################\n')
            grids(liste[element])
            dock(liste[element])
        else:
            print(str(element) + ' --> ' + liste[element])


def check_summary():
    """
    verifie la taille des fichiers resultat, dans le cas ou la taille du fichier est inferieur a 3000 octets cela
    signifie que le docking a echoue ou qu'il n'a pas pu s'effectuer.
    :return:
    """

    error = open('taille_error.txt', 'w')
    if terminal:
        os.system('ls -l AUTODOCKFR/*summary* > taille.txt')

    taille = open('taille.txt', 'r')
    for ligne in taille:
        sep = ligne.split(' ')
        d = 0
        if sep[4] == '':
            d = 1
        if int(sep[4 + d]) < 3000:
            name = sep[-1].split('_')
            error.write(name[0] + '-------->' + sep[4 + d] + '-------->' + sep[-1])
    taille.close()


def ttt_resultat(name,file_result):
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
    plt.bar(year, data1, color="green", bottom=np.array(data4) + np.array(data3) + np.array(data2), label="top1 --> "  + str(data1))

    plt.legend(loc="lower left", bbox_to_anchor=(0.8, 1.0))
    plt.show()


if __name__ == '__main__':
    RESIDUES = lecture_residu()
    if terminal:
        ldoc = ['target', 'AUTODOCKFR']
        for dossier in ldoc:
            if dossier not in listdir(path):
                os.system('mkdir ' + dossier)

    liste = [ '5MYM']
    # grid_for_all()
    #grids('model_alphaFold')
    # dock('1T56')
    for e in liste:
        grids(e)
        dock(e)
    #grid_dock_liste(0, 25)
    check_summary()
    #ttt_resultat_all()
    #check_result()
