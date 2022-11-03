"""
Redock Vina
Editor : Murisier Frederic
Date : mars 2022

redocking for Vina : docking each receptor with his ligand(s) and check if the docking went well

Python version : 3.8
"""
import os
from os import listdir
from os.path import isfile, join

path = os.getcwd()
terminal = False
# variable environnement
PATH = '/home/rene/bin/AutoDockVina1.2/vina_1.2.2_linux_x86_64'


def createBoxTxt(f):
    """
    creation de la box = fichier texte contenant les coordonnees de la box
    :param f: nom de la structure pour laquelle il faut cree la box
    :return:
    """
    fichier = open('target/' + f + '.log', 'r')
    fichier_box = open('VINA/BoxTxt/' + f + '_box.txt', 'w')
    lines = fichier.readlines()
    fichier.close()
    for i in range(0, 7):
        if 'center' in lines[i]:
            center = lines[i].split()[-3:]
            fichier_box.write('center_x = ' + center.pop(0) + '\n')
            fichier_box.write('center_y = ' + center.pop(0) + '\n')
            fichier_box.write('center_z = ' + center.pop(0) + '\n')
        elif 'length' in lines[i]:
            size = lines[i].split()[-3:]
            fichier_box.write('size_x = ' + size.pop(0) + '\n')
            fichier_box.write('size_y = ' + size.pop(0) + '\n')
            fichier_box.write('size_z = ' + size.pop(0))
    fichier_box.close()


def vina(element):
    """
    effectue le docking pour la structure donnee en argument, en prenant les parametres par default, le fichier
    recepteur, les coordonnees de la box et le nom du fichier resultat
    :param element:
    :return:
    """
    commande = '/home/rene/bin/AutoDockVina1.2/vina_1.2.2_linux_x86_64 --receptor RECEPTEUR/PDBQT/' + element + \
               '_recepteur_H.pdbqt --ligand LIGAND/PDBQT/' + element + '_ligand1.pdbqt --config VINA/BoxTxt/' + element\
               + '_box.txt --exhaustiveness=132 --cpu 10  --num_modes 10 --out VINA/RESULT/' + element +\
               '_ligand_vina_out.pdbqt'
    if terminal:
        os.system(commande)
    else:
        print(commande)


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


def result_vina(file):
    """
    effectue la separation des dix modeles (nombre choisit en parametre) contenu dans le fichier resultat de vina en dix
    fichier pour faciliter la comparaison
    :param file: 
    :return: 
    """
    dossier = 'model_' + file
    if terminal:
        if dossier not in listdir(path + '/VINA/RESULT/'):
            os.system('mkdir VINA/RESULT/model_' + file)

    fichier = open('VINA/RESULT/' + file + '_ligand_vina_out.pdbqt', 'r')
    c = 1
    model = open('VINA/RESULT/model_' + file + '/' + file + 'result_' + str(c) + '.pdbqt', 'w')
    for ligne in fichier:
        if 'MODEL' in ligne:
            commande = '/usr/bin/obabel -ipdbqt VINA/RESULT/model_' + file + '/' + file + 'result_' + str(c) +\
                       '.pdbqt -O ' + 'VINA/RESULT/model_' + file + '/' + file + '_model_' + str(c) + '.mol2'
            if terminal:
                os.system(commande)
            else:
                print(commande)
            model.close()
            model = open('VINA/RESULT/model_' + file + '/' + file + 'result_' + str(c) + '.pdbqt', 'w')
            c += 1
        model.write(ligne)
    model.close()
    fichier.close()
    
    # conversion du fichier resultat en mol2 en prevision du rescoring
    c -=1
    commande = '/usr/bin/obabel -ipdbqt VINA/RESULT/model_' + file + '/' + file + 'result_' + str(c) + '.pdbqt -O ' +\
               'VINA/RESULT/model_' + file + '/' + file + '_model_' + str(c) + '.mol2'
    if terminal:
        os.system(commande)
    else:
        print(commande)

    return c


def ecriturelig1(structure, m, lig):
    """
    creation du fichier comparaison et ajout du calcul rmsd par dock rmsd pour le ligand 1 de la structure entree en 
    argument 
    :param structure: 
    :param m: 
    :param lig: 
    :return: 
    """
    # initialisation du document contenant la comparaison
    resultfile = open('VINA/RESULT/comparaison_' + structure + '.txt', 'w')
    resultfile.write('Nombre de models = ' + str(m) + '\n')
    m += 1
    resultfile.write('Model | affinite | Docks RMSD |\n')
    resultfile.write('      | Kcal/mol | Angstrom   |\n')
    # conversion et remplissage du doc
    for n in range(1, m):
        # conversion du model en mol2, pas besoin de le faire pour le ligand car il a deja ete transforme avant
        commande = '/usr/bin/obabel -ipdbqt VINA/RESULT/model_' + structure + '/' + structure + 'result_' + str(n) + \
                   '.pdbqt -O ' + 'VINA/RESULT/model_' + structure + '/' + structure + '_model_' + str(n) + '.mol2'
        if terminal:
            os.system(commande)
        print(commande)

        # dockrmsd : on rescore le rmsd; et on stocke le resultat dans un fichier temporaire (sera ecrase a la prochaine
        # boucle) nome rescore_result.txt
        commande = '/home/rene/bin/DockRMSD/DockRMSD LIGAND/MOL2/' + lig + '.mol2 VINA/RESULT/model_' + structure + '/'\
                   + structure + '_model_' + str(n) + '.mol2 > rescore_result.txt'
        if terminal:
            os.system(commande)
        print(commande)

        # extraction du rescoring
        dockrmsd='NA'
        rescore_file = open('rescore_result.txt', 'r')
        rescore_lignes = rescore_file.readlines()
        for i in rescore_lignes:
            if 'Calculated Docking RMSD:' in i:
                dockrmsd = i.split()[3]
        rescore_file.close()
        print(dockrmsd)
        e = 6
        if n > 9:
            e = 5
        # ecriture dans le doc
        model = open('VINA/RESULT/model_' + structure + '/' + structure + 'result_' + str(n) + '.pdbqt', 'r')
        lignes = model.readlines()
        resultfile.write(str(n) + ' ' * e + lignes[1].split()[3] + ' ' * 6 + str(dockrmsd) + '\n')

    resultfile.close()


def ecritureligs(structure, m, lig):
    """
    ajout du calcul rmsd par dock rmsd pour le ligand 2 ou 3 de la structure entree en argument au fichier comparaison
    :param structure: 
    :param m: 
    :param lig: 
    :return: 
    """
    resultfile = open('VINA/RESULT/comparaison_' + structure + '.txt', 'r')
    old = resultfile.readlines()
    resultfile.close()
    resultfile = open('VINA/RESULT/comparaison_' + structure + '.txt', 'w')
    print(old)
    n=0
    for l in old:
        if l == old[0]:
            resultfile.write(l)
            n = 0
        elif l == old[1]:
            resultfile.write(l[:-1] + ' Docks RMSD lig' + str(n+2) + ' |\n')
        elif l == old[2]:
            resultfile.write(l[:-1] + ' Angstrom        |\n')
        else:
            n += 1

            # dockrmsd : on rescore le rmsd avec le ligand 2 ou 3 ; et on stocke le resultat dans un fichier temporaire
            # (sera ecrase a la prochaine boucle) nome rescore_result.txt
            commande = '/home/rene/bin/DockRMSD/DockRMSD LIGAND/MOL2/' + lig + '.mol2 VINA/RESULT/model_' + structure +\
                       '/' + structure + '_model_' + str(n) + '.mol2 > rescore_result.txt'
            if terminal:
                os.system(commande)
            print(commande)

            # extraction du rescoring
            dockrmsd = 'NA'
            rescore_file = open('rescore_result.txt', 'r')
            rescore_lignes = rescore_file.readlines()
            for i in rescore_lignes:
                if 'Calculated Docking RMSD:' in i:
                    dockrmsd = i.split()[3]
            print(dockrmsd)
            resultfile.write(l[:-1] + ' ' * 12 + dockrmsd + '\n')

    resultfile.close()


def ecriture_energie_galaxy(structure):
    """
    ajout dans le fichier comparaison de du rescoring effectuee avec le champ de force ge galaxyDock
    :param structure: 
    :return: 
    """
    boxfile = open('VINA/BoxTxt/' + structure + '_box.txt', 'r')
    box = boxfile.readlines()
    boxfile.close()
    resultfile = open('VINA/RESULT/model_' + structure + '/comparaison_' + structure + '.txt', 'r')
    old = resultfile.readlines()
    resultfile.close()
    resultfile = open('VINA/RESULT/model_' + structure + '/comparaison_' + structure + '.txt', 'w')
    print(old)
    n=0
    for l in old:
        if l == old[0]:
            resultfile.write(l)
            n = 0
        elif l == old[1]:
            resultfile.write(l[:-1] + ' GalaxyDock |\n')
        elif l == old[2]:
            resultfile.write(l[:-1] + ' Kcal/mol   |\n')
        else:
            n += 1
            # galaxyDock : on rescore l'energie' avec le ligand 2 ou 3 ; et on obtient deux fichier de resultat
            if terminal:
                os.system('/home/rene/bin/GalaxyDock3/script/calc_energy.py -d /home/rene/bin/GalaxyDock3 -p RECEPTEUR/'
                          + structure + '_recepteur.pdb -l VINA/RESULT/model_' + structure + '/' + structure + '_model_'
                          + str(n) + '.mol2 -x ' + box[0].split()[-1] + ' -y ' + box[1].split()[-1] + ' -z ' +
                          box[2].split()[-1] + ' -size_x ' + box[3].split()[-1] + ' -size_y ' + box[4].split()[-1] +
                          ' -size_z ' + box[5].split()[-1])
            print('model numero' + str(n))

            # extraction du resuultat du rescoring
            fileresult = open('calc_energy.log', 'r')
            energy = fileresult.readlines()
            fileresult.close()
            galaxy_energy = energy[0].split()[-1]
            print(galaxy_energy)
            resultfile.write(l[:-1] + ' ' * 8 + galaxy_energy + '\n')

    resultfile.close()


def comparaison(structure, m):
    """
    creation et remplissage du fichier comparaison pour la structure donnee en argument 
    :param structure: 
    :param m: 
    :return: 
    """
    stop = True
    lig = '_ligand'
    c = 1
    while stop:
        name = structure + lig + str(c)
        if c == 1:
            print(name)
            ecriturelig1(structure, m, name)
        elif name in listeMol2:
            ecritureligs(structure, m, name)
        else:
            stop = False
        c += 1


def galaxy():
    """
    fonction qui va calculer l'energie avec le champ de force de galaxydock
    resultfile = open('VINA/RESULT/model_1T56/comparaison_1T56.txt', 'r')
    old = resultfile.readlines()
    resultfile.close()
    resultfile = open('VINA/RESULT/model_1T56/comparaison_1T56.txt', 'w')
    print(old)
    for l in old:
        if l == old[0]:
            resultfile.write(l)
        elif l == old[1]:
            resultfile.write(l[:-1] + ' Docks RMSD lig2 |\n')
        elif l == old[2]:
            resultfile.write(l[:-1] + ' Angstrom        |\n')
        else:
            resultfile.write(l[:-1] + ' ' * 12 + 'xxx\n')

    resultfile.close()

    :return:
    """
    structure = '1T56'
    n = 1
    boxfile = open('VINA/BoxTxt/' + structure + '_box.txt', 'r')
    box = boxfile.readlines()

    print('/home/rene/bin/GalaxyDock3/script/calc_energy.py -d /home/rene/bin/GalaxyDock3 -p '
          'RECEPTEUR/1T56_recepteur.pdb -l VINA/RESULT/model_1T56/1T56_model_1.mol2 -x 26.441 -y 11.537 -z 14.814 '
          '-size_x 76 -size_y 86 -size_z 72')
    print('/home/rene/bin/GalaxyDock3/script/calc_energy.py -d /home/rene/bin/GalaxyDock3 -p RECEPTEUR/' + structure +
          '_recepteur.pdb -l VINA/RESULT/model_' + structure + '/' + structure + '_model_' + str(n) +
          '.mol2 -x ' + box[0].split()[-1] + ' -y ' + box[1].split()[-1] + ' -z ' + box[2].split()[-1] + ' -size_x ' +
          box[3].split()[-1] + ' -size_y ' + box[4].split()[-1] + ' -size_z ' + box[5].split()[-1])

    boxfile.close()


if __name__ == '__main__':
    if terminal:
        os.system('mkdir VINA')
        os.system('mkdir VINA/BoxTxt')
        os.system('mkdir VINA/RESULT')
    listeMol2 = liste_file('.mol2', '/LIGAND/MOL2/')
    l = liste_file('.log', '/target')
    #l = ['3Q0V', '5MYM', '5NZ1']

    modelliste = []

    for i in l:
        createBoxTxt(i)
        #vina(i)
        #models = result_vina(i)
        #comparaison(i, models)
        #ecriture_energie_galaxy(i)



    # comparaison('5MYN', 10)
    # ecriture_energie_galaxy(i)
