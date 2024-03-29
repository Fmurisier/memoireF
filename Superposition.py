
###################################################################
# Description: Superposition pymol                                #
# # Superpose every structure to the structure reference given in #
# # argument listing every (compressed) pdb file present in the   #
# # repertory then create a pymol script with the instruction     #
# # to superpose every structure to the reference structure then  #
# # execute the script                                            #
# Python version : 3.8                                            #
# Date: 2022/02                                                   #
# Author: Frederic MURISIER                                       #
###################################################################


import sys
from os import listdir
import os
from os.path import isfile, join
import datetime
import pathlib

path = os.getcwd()
dossier_path = '/../Donnee_memoire/'
terminal = False


def liste_file(patern, file_path=''):
    """
    search in the directory given in argument the file with the pattern given in argument and create a list of file
    :param patern:
    :param file_path:
    :return:
    """
    list_file = []
    fichiers = [f for f in listdir(path + file_path) if isfile(join(path + file_path, f))]

    # each name is splited by the symbol '_' to find only the file 'transformed
    for ligne in fichiers:
        if patern in ligne:
            name = ligne.split('_')
            list_file.append(name[0])
    list_file.sort()
    return list_file


def ecriture_pymol_all(liste_pdb, ref, box):
    """
    write a pymol script with the command line to superpose each structure from the list given in argument to the
    reference structure also given in argument. the box argument is the box where the structure are superposed
    :param ref:
    :param liste_pdb:
    :return:
    """
    fichier = open('scriptPymol.pml', 'w')
    fichier.write('output = open("rmsd_result.txt", "w")\n')
    fichier_log.write('WRITING PYMOL SCRIPT\n\n')
    for e in liste_pdb:
        name = e.split('.')[0]
        ref_name = ref.split('.')[0]
        fichier.write('load ' + ref + '\n')
        fichier.write('load ' + e + '\n')
        # selection of the residue to be superposed
        fichier.write('select ' + ref_name + box + ref_name + '\n')
        fichier.write('select ' + name + box + name + '\n')
        # superposition
        fichier.write('super ' + name + '_poche////CA, ' + ref_name + '_poche////CA\n')
        fichier.write('save ' + name + '_transformed.pdb, ' + name + '\n')

        fichier.write('data = cmd.super("' + name + '_poche", "' + ref_name + '_poche")\n')
        fichier.write('output.write("' + name + '=")\n')
        fichier.write('output.write(" %f\\n" % data[0])\n')

    fichier.write('output.close()\n')
    fichier.write('print("END")\n quit')
    fichier.close()
    fichier_log.write('PYMOL SCRIPT DONE\n\n')


def ecriture_box(liste=False):
    """
    read the residus.txt file, extract the name of each residus and create the box
    :param liste:
    :return:
    """
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
            information = 'The box is made of ' + str(len(resi_liste)) + \
                  ' residus, if it is not the case then please check that the file is filled correctly (' \
                  ' exemple : \'residues A:LEU87,LEU90,MET102,TRP103,ILE107,... \')'
            print(information)
            fichier_log.write(information + '\n\n')

        commande = '_poche, resi ' + '+'.join(resi_liste) + ' and model '
        if liste:
            commande = resi_liste
        return commande, resi_code3
    else:
        error_message = 'error file residus.txt don\'t exist ! '
        print(error_message)
        fichier_log.write(error_message + '\n\n')
        return 'error', ''


def check_file():
    """
    verification of the argument given by the user: if the reference structure exist or not
    if it exists check if the file is compressed or not. if it is then the function will decompress it automatically
    :param:
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
            fichier_log.write('Decompression ref file\n\n')
            decompression_pdb(x)
        b, resi3 = ecriture_box()
        check_box = lecture_ref_file(x[:-3], resi3)

        if b == 'error' and check_box:
            r = False

    return r, x, b


def lecture_ref_file(ref, liste_code3):
    """
    Check if the residus are in the reference structure
    :param ref:
    :param liste_code3:
    :return:
    """
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
        warning_message = 'WARNING !!! Error the residus of the box are not present in the structure, ' \
                          'please check again the residus'
        print(warning_message)
        fichier_log.write(warning_message+ '\n\n')

    return box_ok


def decompression_pdb(file):
    """
    decompress the file given in argument
    :param file:
    :return:
    """
    fichier = open('scriptPymol.pml', 'w')
    fichier.write('load ' + file + '\n')
    fichier.write('save ' + file[:-3] + '\n')
    fichier.write('print("END")\n quit')
    fichier.close()
    os.system('pymol -cp scriptPymol.pml')
    fichier_log.write('Decompression of the reference file\n\n')


def check_rmsd():
    """
    Check if the rmsd of the superposition is higher than 1. If it is then give a warning
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
            fichier_log.write('\nWARNING ' + str(l)+ '\n\n')
    if not problem:
        print('\nNo RMSD problem detected\n')
        fichier_log.write('\nNo RMSD problem detected\n\n')


def superpose_all():
    """
    Check if the reference structure exists in the directory if yes then executes the function writing the pymol file
    :return:
    """
    fichier_log.write('START CHECK SUPERPOSITION\n\n')
    liste_pdb = liste_file('pdb')
    for i in range(len(liste_pdb)):
        liste_pdb[i] = liste_pdb[i].split('.')[0]
    present, ref_structure, b = check_file()
    if present:
        fichier_log.write('CHECK SUPERPOSITION OK\n\n')
        fichier_log.write('START SUPERPOSITION\n\n')
        ecriture_pymol_all(liste_file('pdb'), ref_structure, b)
        if terminal:
            os.system('pymol -cp scriptPymol.pml')
            check_rmsd()
            os.system('mv *_transformed.pdb transformed/')
    else:
        error_mess = 'File absent from directory, try : python superposition.py ValidFileName\n\n'
        print(error_mess)
        fichier_log.write(error_mess)


def superpose_liste(file_p):
    """
    superpose a list of structure (usualy the ones that failed)
    :param file_p:
    :return:
    """
    present, ref_structure, b = check_file()
    ecriture_pymol_all(file_p, ref_structure, b)
    if terminal:
        os.system('pymol -cp scriptPymol.pml')
        check_rmsd()
        # os.system('mv *_transformed.pdb transformed/')


if __name__ == '__main__':

    fichier_log = open('log_pipeline.txt', 'a')
    fichier_log.write(str(datetime.date.today()) + ' -' * 20 + '\n\n')
    fichier_log.write('START \n\n')
    if 'transformed' not in listdir(path):
        os.system('mkdir transformed')

    superpose_all()
    file_prob = []
    superpose_liste(file_prob)

    fichier_log.close()

