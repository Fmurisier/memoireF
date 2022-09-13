"""
Superposition pymol
Editor : Murisier Frederic

superpose toutes les structures par rapport a la structure 1T56

fait la liste de tous les fichiers pdb compresse ou pas dans le repertoire dans lequel le script.py se trouve, cree un
fichier script.pml contenant les lignes commande necessaire a la superposition de toutes les structures par rapport a la
structure de reference.
fait executer via le terminal le scipt pour pymol resultant la creation de fichier des structures superposee
Date : 2022
cree une liste de tous les fichiers pdb ou pdb.gz contenu dans le dossier et qui n'ont pas deja ete transforme
cree un script pymol qui sera lance a la fin de l'algorithme pour superposer toutes les structures par rapport a la
structure de reference (1T56)
Resultat : nouveaux fichier transforme

Python version : 3.8
"""

import sys
from os import listdir
import os
from os.path import isfile, join
import pathlib

path = os.getcwd()
terminal = False


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


def ecriture_pymol_all(liste_pdb, ref):
    """
    recoit en argument une liste de fichier et une structure de reference, la fonction va ensuite ecrire le script pymol
    dans un fichier = les lignes de commande pymol effectuant les superpositions tout en enregistrant les nouveaux
    fichiers superpose dans de nouveaux fichier
    :param ref:
    :param liste_pdb:
    :return:
    """
    fichier = open('scriptPymol.pml', 'w')
    fichier.write('output = open("rmsd_result.txt", "w")\n')
    for e in liste_pdb:
        fichier.write('load ' + ref + '.pdb.gz\n')
        fichier.write('load ' + e + '.pdb.gz\n')
        fichier.write('super ' + e + '//A//CA, ' + ref + '//A//CA\n')
        fichier.write('save ' + e + '_transformed.pdb, ' + e + '\n')

        fichier.write('data = cmd.super("' + ref + '", "' + e + '")\n')
        fichier.write('output.write("' + e + '=")\n')
        fichier.write('output.write(" %f\\n" % data[0])\n')

    fichier.write('output.close()\n')
    fichier.write('print("END")\n quit')
    fichier.close()


def ecriture_box():

    dossier = 'residus.txt'
    if dossier not in listdir(path):
        print('error file residus.txt don\'t exist ! ')
    resi_file = open('residus.txt', 'r').readlines()
    resi_liste = []
    for e in resi_file:
        if '\n' in e:
            e = e[:-1]
            new_line = e.split('+')
            for e in new_line:
                if e.isdigit():
                    resi_liste.append(e)
    print(resi_liste)
    print('La box est composee de ' + str(len(resi_liste)) + ' residus, si ce n\'est pas le cas verifiez que le fichier'
                                                             ' residu soie ecrit correctement ( chaque residu doit etre'
                                                             ' separe par un +)')

    #commande = 'select ' + name + '_poche, resi ' + '+'.join(resi_liste) + ' and model ' + name
    commande = '_poche, resi ' + '+'.join(resi_liste) + ' and model '
    print(commande)



def check_file(lfichier):
    """
    regarde ce que l'utilisateur a entre comme structure de reference et verifie si celle ci est presente dans le
    dossier par defaut 1T56 si le fichier est inexistant alors message d erreur
    :param lfichier:
    :return:
    """
    r = False

    x = sys.argv
    x.insert(1, '1T56')
    print(x)
    x = x[-1]
    if x in lfichier:
        r = True

    return r, x


def check_rmsd():
    """
    verifie si le rmsd est plus haut que 1, si c'est le cas affiche un message d'erreur
    :return:
    """
    problem = False
    fichier = open('rmsd_result.txt', 'r')
    for l in fichier:
        ligne = l.split('=')
        # on converti le string en float et on ne doit pas prendre le dernier caratere qui est le passage a la ligne
        # checking if the rmsd is higher than 1
        # if it is print a warning message with the name of the structure and the rmsd value
        if float(ligne[1][:-1]) > 1:
            problem = True
            print('\nWARNING ' + str(l))
    if not problem:
        print('\nNo RMSD problem detected\n')


def superpose_all():
    """
    prend l'argument donne par l'utilisateur lors de l'entree de la ligne de commande, check si la structure est
    presente dans le fichier, si non affiche un message d'erreur
    si oui lance le programme pour l'ecriture du script pymol
    :return:
    """
    liste_pdb = liste_file('pdb')
    for i in range(len(liste_pdb)):
        liste_pdb[i] = liste_pdb[i].split('.')[0]
    present, ref_structure = check_file(liste_pdb)
    if present:
        ecriture_pymol_all(liste_pdb, ref_structure)
        if terminal:
            os.system('pymol -cp scriptPymol.pml')
    else:
        print('Fichier absent du repertoire entrez la commande : python superposition.py NomFichierValide')

    check_rmsd()


def superpose_liste(liste_super):
    """
    recoit une liste et va effectuer la supperposition de chaque structure de la liste sur la structure de reference
    entre par l'utilisateur, ici par default 1T56
    :param liste_super:
    :return:
    """
    ecriture_pymol_all(liste_super, '1T56')
    if terminal:
        os.system('pymol -cp scriptPymol.pml')


if __name__ == '__main__':
    if 'transformed' not in listdir(path):
        os.system('mkdir transformed')


    superpose_all()
    #os.system('mv *_transformed.pdb transformed/')


