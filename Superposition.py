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
dossier_path = '/../Donnee_memoire/'
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


def ecriture_pymol_all(liste_pdb, ref, box):
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
        name = e.split('.')[0]
        ref_name = ref.split('.')[0]
        fichier.write('load ' + ref + '\n')
        fichier.write('load ' + e + '\n')
        # selection des poches a superposer, d'abord celle de reference puis celle a superposer
        fichier.write('select ' + ref_name + box + ref_name + '\n')
        fichier.write('select ' + name + box + name + '\n')
        # superposition des deux poches
        fichier.write('super ' + name + '_poche////CA, ' + ref_name + '_poche////CA\n')
        fichier.write('save ' + name + '_transformed.pdb, ' + name + '\n')

        fichier.write('data = cmd.super("' + name + '_poche", "' + ref_name + '_poche")\n')
        fichier.write('output.write("' + name + '=")\n')
        fichier.write('output.write(" %f\\n" % data[0])\n')

    fichier.write('output.close()\n')
    fichier.write('print("END")\n quit')
    fichier.close()


def ecriture_box():
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
        print('The box is made of ' + str(len(resi_liste)) + \
              ' residus, if it is not the case then please check that the file is filled correctly (' \
              ' exemple : \'residues A:LEU87,LEU90,MET102,TRP103,ILE107,... \')')

        commande = '_poche, resi ' + '+'.join(resi_liste) + ' and model '
        return commande, resi_code3
    else:
        print('error file residus.txt don\'t exist ! ')
        return 'error', ''


def check_file():
    """
    regarde ce que l'utilisateur a entre comme structure de reference et verifie si celle ci est presente dans le
    dossier par defaut 1T56 si le fichier est inexistant alors message d erreur
    :param lfichier:
    :return:
    """
    r = False
    b = 'error'
    x = sys.argv
    x.insert(1, '1T56.pdb.gz')
    x = x[-1]
    check_box = False

    if x in [f for f in listdir(path + dossier_path) if isfile(join(path + dossier_path, f))]:
        r = True
        if x[-4:] != '.pdb':
            decompression_pdb(x)
        b, resi3 = ecriture_box()
        check_box = lecture_ref_file(x[:-3], resi3)


        if b == 'error' and check_box:
            r = False

    return r, x, b


def lecture_ref_file(ref, liste_code3):
    #ref_file = open('../Donnee_memoire/model_alphaFold.pdb', 'r').readlines()
    ref_file = open(dossier_path + ref, 'r').readlines()
    dico_res = {}
    for e in ref_file:
        if 'ATOM' in e:
            e = e.split(' ')
            for i in range(e.count('')):
                e.remove('')
            if e[5] not in dico_res:
                dico_res[e[5]] = e[3]

    box_ok = True
    for e in liste_code3:
        code = e[:3]
        resnum = e[3:]
        if dico_res.get(resnum) != code:
            box_ok = False
    if not box_ok:
        print('WARNING !!! Error the residus of the box are not present in the structure,'
              ' please check again the residus')
    return box_ok


def decompression_pdb(file):
    fichier = open('scriptPymol.pml', 'w')
    fichier.write('load ' + file + '\n')
    fichier.write('save ' + file[:-3] + '\n')
    fichier.write('print("END")\n quit')
    fichier.close()
    os.system('pymol -cp scriptPymol.pml')


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
    present, ref_structure, b = check_file()
    if present:
        ecriture_pymol_all(liste_file('pdb'), ref_structure, b)
        if terminal:
            os.system('pymol -cp scriptPymol.pml')
            check_rmsd()
            os.system('mv *_transformed.pdb transformed/')
    else:
        print('Fichier absent du repertoire entrez la commande : python superposition.py NomFichierValide')


def superpose_liste(file_prob):
    present, ref_structure, b = check_file()
    ecriture_pymol_all(file_prob, ref_structure, b)
    if terminal:
        os.system('pymol -cp scriptPymol.pml')
        check_rmsd()
        os.system('mv *_transformed.pdb transformed/')

if __name__ == '__main__':
    if 'transformed' not in listdir(path):
        os.system('mkdir transformed')

    superpose_all()
    file_prob = []
    superpose_liste(file_prob)

