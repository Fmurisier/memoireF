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


def obtention_liste_pdb():
    """
    va effectuer la liste des fichiers pdb et pdb compresse dans le repertoire dans lequel le script.py se trouve
    :return : list
    """
    list_pdb = []
    list_ext = ['transformed', 'ligand', 'recepteur']
    c = 0
    # fait la liste des fichier dans le repertoire ou se trouve script.py
    fichiers = [f for f in listdir(path) if isfile(join(path, f))]

    # pour chaque nom de fichier on separe le nom et son/ses extension(s)
    for ligne in fichiers:
        name = ligne.split('.')

        # on va selectionner uniquement les fichiers pdb et les fichiers pdb compresse
        if name[-2] == 'pdb' or name[-1] == 'pdb':
            # on separe le nom du fichier pour verifier s'il n'y a pas deja des fichiers qui ont ete transforme et
            # eviter de les resuperposer une deuxieme fois et cree des doublons
            name0 = name[0].split('_')

            add = True
            for e in list_ext:
                for i in name0:
                    if e in i:
                        add = False

            # ajoute chaque nom de fichier dans la liste
            if add:
                list_pdb.append(name[0])

    # Trie la liste, pas forcement necessaire mais utile pour nous car alphabetiquement notre structure de ref est 1ere
    list_pdb.sort()
    # print('Nombre de pdb file non transforme dans le dossier : ', len(list_pdb))

    return list_pdb


def ecriture_pymol(name):
    """
    ecris un fichier contenant le script pour pymol pour superposer la structure fournie en argument avec 1T56
    :param name:
    :return:
    """
    ref = '1T56'
    fichier = open('scriptPymol.pml', 'w')

    fichier.write('load ' + ref + '.pdb.gz\n')
    fichier.write('load ' + name + '.pdb.gz\n')
    # on ajoute //A pour que la superposition s'effectue sur la chaine A au cas ou le ligand aurais plusieurs chaines
    fichier.write('super ' + name + '//A//CA, ' + ref + '//A//CA\n')
    # fichier.write('save ' + name + '_transformed.pdb, ' + name + '\n')
    fichier.write('save ../file_prep/' + name + '_transformed.pdb, ' + name + '\n')
    fichier.close()


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
    liste_pdb = obtention_liste_pdb()
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
    os.system('mv *_transformed.pdb transformed/')


