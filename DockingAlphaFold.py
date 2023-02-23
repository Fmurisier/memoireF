"""
AlphaFold
Editor : Murisier Frederic
Date : 04/2022

Python version : 3.8
"""
import datetime
import os
from os import listdir
from os.path import isfile, join

path = os.getcwd()
terminal = True
# variable environnement
PATH = '/home/rene/bin/'
RESIDUES = 'residues A:LEU87,LEU90,MET102,TRP103,ILE107,PHE110,PHE114,TRP138,MET142,TRP145,TYR148,THR149,VAL152,' \
           'ASN176,ASN179,GLU180,LEU183,PHE184,TRP207'


def liste_file(patern, file_path=''):
    list_file = []

    # fait la liste des fichier dans le repertoire ou se trouve script.py
    fichiers = [f for f in listdir(path + file_path) if isfile(join(path + file_path, f))]

    # pour chaque nom de fichier on separe le nom par '_' pour chercher uniquement les fichiers comportant le prefix
    # '_transformed.pdb' et on stocke tous les noms de fichier dans une liste
    for ligne in fichiers:
        if patern in ligne:
            name = ligne.split('.')
            list_file.append(name[0])
    list_file.sort()

    return list_file


def run_lig():
    lig = liste_file('.', '/LIGAND/MOL2')


def ecriture_pymol(e, ref):
    """
    recoit en argument une liste de fichier et une structure de reference, la fonction va ensuite ecrire le script pymol
    dans un fichier = les lignes de commande pymol effectuant les superposition tout en enregistrant les nouveaux
    fichiers superpose dans de nouveaux fichier
    :param ref:
    :param liste_pdb:
    :return:
    """

    fichier = open('scriptPymol.pml', 'w')
    fichier.write('output = open("rmsd_result.txt", "w")\n')

    fichier.write('load ' + ref + '.pdb.gz\n')
    fichier.write('load ' + e + '.pdb\n')
    fichier.write('super ' + e + '//A//CA, ' + ref + '//A//CA\n')
    fichier.write('save ' + e + '_transformed.pdb, ' + e + '\n')

    fichier.write('data = cmd.super("' + ref + '", "' + e + '")\n')
    fichier.write('output.write("' + e + '=")\n')
    fichier.write('output.write(" %f\\n" % data[0])\n')

    fichier.write('output.close()\n')
    fichier.write('print("END")\n quit')
    fichier.close()

    if terminal:
        os.system('pymol -cp scriptPymol.pml')


def conversion_recepteur_pdbqt(nomfichier):
    commande1 = 'python3 ../../../home/rene/bin/pdb2pqr/pdb2pqr --titration-state-method propka --drop-water ' \
                '--pdb-output ' + nomfichier + '_H.pdb --with-ph 7.4 ' + nomfichier + '_transformed.pdb ' + nomfichier \
                + '_H.pqr'
    commande2 = '/home/rene/bin/ADFRsuite-1.0/bin/prepare_receptor -r ' + nomfichier + '_H.pdb'

    if terminal:
        os.system(commande1)
        os.system(commande2)
    else:
        print(commande1)
        print(commande2)


    dossier = 'RECEPTEUR'
    if terminal:
        if dossier not in listdir(path + '/AlphaFold'):
            os.system('mkdir AlphaFold/' + dossier)
        os.system('mv model_alphaFold* AlphaFold/' + dossier)
        os.system('cp AlphaFold/model_alphaFold.pdb ~/pdb1')


def grids(fichier):
    """
    /home/rene/bin/ADFRsuite-1.0/bin/agfr -b residues A:LEU87,LEU90,MET102,TRP103,ILE107,PHE110,PHE114,TRP138,MET142,
    TRP145,TYR148,THR149,VAL152,ASN176,ASN179,GLU180,LEU183,PHE184,TRP207 -r RECEPTEUR/PDBQT/1T56_recepteur_H.pdbqt -o
    1T56
    :param fichier:
    :return:
    """
    dossier = 'target'
    if terminal:
        if dossier not in listdir(path + '/AlphaFold'):
            os.system('mkdir AlphaFold/' + dossier)

    boucle = True
    while boucle:
        commande = '/home/rene/bin/ADFRsuite-1.0/bin/agfr -b ' + RESIDUES + ' -r AlphaFold/RECEPTEUR/' + fichier + \
                   '_H.pdbqt -o AlphaFold/target/' + fichier
        if terminal:
            os.system(commande)
        #else:
            print(commande)
        boucle = check_error_grid(fichier + '.log')


def check_error_grid(file):
    """

    :param file:
    :return:
    """
    restart = False
    fichier = open('AlphaFold/target/' + file, 'r')
    fichierlog = open('fileERROR.txt', 'a')
    for l in fichier:
        if 'ERROR' in l:
            print('\nWARNING PROBLEM agdr crashes....\n')
            print(l)
            fichierlog.write('\nWARNING PROBLEM agdr crashes....\n' + l)
            restart = True
    if not restart:
        print('\nNo problem detected in the grid creation with model alphaFold\n')
    fichierlog.close()
    return restart


def dockadfr(ligand):
    """
    /home/rene/bin/ADFRsuite-1.0/bin/adfr -l LIGAND/PDBQT/1T56_ligand1_smile.pdbqt -t target/1T56.trg -r
    LIGAND/PDBQT/1T56_ligand1.pdbqt -p 150 -e 5000000 -n 100 -c 6 -O -o 1T56_result -J dock

    :param ligand:
    :return:
    """
    dossier = 'AUTODOCKFR'
    if terminal:
        if dossier not in listdir(path + '/AlphaFold'):
            os.system('mkdir AlphaFold/' + dossier)

    commande = '/home/rene/bin/ADFRsuite-1.0/bin/adfr -l LIGAND/PDBQT/' + ligand + '_smile.pdbqt -t AlphaFold/target/'\
               'model_alphaFold.trg -r LIGAND/PDBQT/' + ligand + '.pdbqt -p 150 -e 5000000 -n 100 -c 6 -O -o ' \
               'AlphaFold/AUTODOCKFR/' + ligand + '_result -J dock'
    if terminal:
        os.system(commande)
    else:
        print(commande)


def check_summary():
    error = open('AlphaFold/AUTODOCKFR/taille_error.txt', 'w')
    if terminal:
        os.system('ls -l AlphaFold/AUTODOCKFR/*summary* > AlphaFold/AUTODOCKFR/taille.txt')

    taille = open('AlphaFold/AUTODOCKFR/taille.txt', 'r')
    for ligne in taille:
        sep = ligne.split(' ')
        d = 0
        if sep[4] == '':
            d = 1
        # print(sep[4 + d] + '-------->' + sep[-1])
        if int(sep[4 + d]) < 3000:
            name = sep[-1].split('_')
            error.write(name[0] + '-------->' + sep[4 + d] + '-------->' + sep[-1])
            # dock(name[0])
    taille.close()
    error.close()


def createBoxTxt(f):
    dossier = 'VINA'
    if terminal:
        if dossier not in listdir(path + '/AlphaFold'):
            os.system('mkdir AlphaFold/' + dossier)
    fichier = open('AlphaFold/target/' + f + '.log', 'r')
    fichier_box = open('AlphaFold/VINA/' + f + '_box.txt', 'w')
    lines = fichier.readlines()
    fichier.close()
    # print(lines[0:7])
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


def vina(element, ref):
    commande = '/home/rene/bin/AutoDockVina1.2/vina_1.2.2_linux_x86_64 --receptor ' \
               'AlphaFold/RECEPTEUR/model_alphaFold_H.pdbqt --ligand LIGAND/PDBQT/' + element + \
               '.pdbqt --config AlphaFold/VINA/' + ref + '_box.txt --exhaustiveness=132 --cpu 10  --num_modes 10 ' \
                                                         '--out AlphaFold/VINA/' + element + '_vina_out.pdbqt'
    if terminal:
        os.system(commande)
    else:
        print(commande)


def result_vina(file):
    print('####################### result vina #################################')
    dossier = 'RESULT'
    if terminal:
        if dossier not in listdir(path + '/AlphaFold/VINA'):
            os.system('mkdir AlphaFold/VINA/' + dossier)
    dossier = 'model_' + file
    if terminal:
        if dossier not in listdir(path + '/AlphaFold/VINA/RESULT/'):
            os.system('mkdir AlphaFold/VINA/RESULT/model_' + file)
    fichier = open('AlphaFold/VINA/' + file + '_vina_out.pdbqt', 'r')
    c = 1
    model = open('AlphaFold/VINA/RESULT/model_' + file + '/' + file + '_result_' + str(c) + '.pdbqt', 'w')
    for ligne in fichier:
        if 'MODEL' in ligne:
            commande = '/usr/bin/obabel -ipdbqt AlphaFold/VINA/RESULT/model_' + file + '/' + file + '_result_' + str(c)\
                    + '.pdbqt -O ' + 'AlphaFold/VINA/RESULT/model_' + file + '/' + file + '_model_' + str(c) + '.mol2'
            if terminal:
                os.system(commande)
            else:
                print(commande)

            model.close()
            model = open('AlphaFold/VINA/RESULT/model_' + file + '/' + file + '_result_' + str(c) + '.pdbqt', 'w')
            c += 1
        model.write(ligne)
    model.close()
    fichier.close()
    # conversion du fichier resultat en mol2 en prevision du rescoring
    c -= 1
    commande = '/usr/bin/obabel -ipdbqt AlphaFold/VINA/RESULT/model_' + file + '/' + file + '_result_' + str(c) + \
               '.pdbqt -O ' + 'AlphaFold/VINA/RESULT/model_' + file + '/' + file + '_model_' + str(c) + '.mol2'
    if terminal:
        os.system(commande)
    else:
        print(commande)

    return c


def ecriturelig1(structure, m):
    # initialisation du document contenant la comparaison
    resultfile = open('AlphaFold/VINA/RESULT/comparaison_' + structure + '.txt', 'w')
    resultfile.write('Nombre de models = ' + str(m) + '\n')
    m += 1
    resultfile.write('Model | affinite | Docks RMSD |\n')
    resultfile.write('      | Kcal/mol | Angstrom   |\n')
    # conversion et remplissage du doc
    for n in range(1, m):
        # conversion du model en mol2, pas besoin de le faire pour le ligand car il a deja ete transforme avant
        commande = '/usr/bin/obabel -ipdbqt AlphaFold/VINA/RESULT/model_' + structure + '/' + structure + '_result_' + \
                   str(n) + '.pdbqt -O ' + 'AlphaFold/VINA/RESULT/model_' + structure + '/' + structure + '_model_' + \
                   str(n) + '.mol2'
        if terminal:
            os.system(commande)
        else:
            print(commande)

        # dockrmsd : on rescore le rmsd; et on stocke le resultat dans un fichier temporaire (sera ecrase a la prochaine
        # boucle) nome rescore_result.txt
        commande = '/home/rene/bin/DockRMSD/DockRMSD LIGAND/MOL2/' + lig + '.mol2 AlphaFold/VINA/RESULT/model_' + \
                   structure + '/' + structure + '_model_' + str(n) + '.mol2 > rescore_result.txt'
        if terminal:
            os.system(commande)
        #else:
            print(commande)

        # extraction du rescoring
        dockrmsd = 'NA'
        rescore_file = open('rescore_result.txt', 'r')
        rescore_lignes = rescore_file.readlines()
        for i in rescore_lignes:
            if 'Calculated Docking RMSD:' in i:
                dockrmsd = i.split()[3]
        rescore_file.close()
        print(structure, dockrmsd)
        e = 6
        if n > 9:
            e = 5
        # ecriture dans le doc
        model = open('AlphaFold/VINA/RESULT/model_' + structure + '/' + structure + '_result_' + str(n) + '.pdbqt', 'r')

        lignes = model.readlines()
        resultfile.write(str(n) + ' ' * e + lignes[1].split()[3] + ' ' * 6 + str(dockrmsd) + '\n')

    resultfile.close()


def ecritureligs(structure, m):
    resultfile = open('AlphaFold/VINA/RESULT/model_' + structure + '/comparaison_' + structure + '.txt', 'r')
    old = resultfile.readlines()
    resultfile.close()
    resultfile = open('AlphaFold/VINA/RESULT/model_' + structure + '/comparaison_' + structure + '.txt', 'w')
    print(old)
    n=0
    for l in old:
        if l == old[0]:
            resultfile.write(l)
            n = 0
        elif l == old[1]:
            resultfile.write(l[:-1] + ' Docks RMSD lig2 |\n')
        elif l == old[2]:
            resultfile.write(l[:-1] + ' Angstrom        |\n')
        else:
            n += 1

            # dockrmsd : on rescore le rmsd avec le ligand 2 ou 3 ; et on stocke le resultat dans un fichier temporaire
            # (sera ecrase a la prochaine boucle) nome rescore_result.txt
            commande = '/home/rene/bin/DockRMSD/DockRMSD LIGAND/MOL2/' + structure + \
                       '.mol2 AlphaFold/VINA/RESULT/model_' + structure + '/' + structure + '_model_' + str(n) + \
                       '.mol2 > rescore_result.txt'
            if terminal:
                os.system(commande)
            else:
                print(commande)

            # extraction du rescoring
            dockrmsd = 'NA'
            rescore_file = open('rescore_result.txt', 'r')
            rescore_lignes = rescore_file.readlines()
            for i in rescore_lignes:
                if 'Calculated Docking RMSD:' in i:
                    dockrmsd = i.split()[3]
            print(structure, dockrmsd)
            resultfile.write(l[:-1] + ' ' * 12 + dockrmsd + '\n')

    resultfile.close()


def comparaison(structure, m):
    stop = True
    c = 1
    while stop:
        name = structure + str(c)
        if c == 1:
            ecriturelig1(structure, m)
        elif name in liste_lig:
            ecritureligs(structure, m)
        else:
            stop = False
        c += 1


def dock_galaxy(e):
    dossier = 'GalaxyD'
    if terminal:
        if dossier not in listdir(path + '/AlphaFold/'):
            os.system('mkdir AlphaFold/' + dossier)
    if terminal:
        if e not in listdir(path + '/AlphaFold/GalaxyD/'):
            os.system('mkdir AlphaFold/GalaxyD/' + e)

    file = open('AlphaFold/VINA/model_alphaFold_box.txt', 'r').readlines()
    cx = file[0].split()[-1]
    cy = file[1].split()[-1]
    cz = file[2].split()[-1]
    sx = file[3].split()[-1]
    sy = file[4].split()[-1]
    sz = file[5].split()[-1]

    # galaxy dock
    commande = '/home/rene/bin/GalaxyDock3/script/run_GalaxyDock3.py -d /home/rene/bin/GalaxyDock3 -p ' + \
               'AlphaFold/RECEPTEUR/model_alphaFold_transformed.pdb -l LIGAND/MOL2/' + e + '.mol2 -x ' + cx + ' -y ' +\
               cy + ' -z ' + cz + ' -size_x ' + sx + ' -size_y ' + sy + ' -size_z ' + sz + ' --n_proc 1'

    if terminal:
        os.system(commande)
    else:
        print('docking--->' + commande)

    # split le resultat
    commande = 'python2 /home/rene/bin/GalaxyDock3/script/split_mol2.py GD3_cl.mol2 AlphaFold/GalaxyD/' + e + '/' + e
    if terminal:
        os.system(commande)
    else:
        print('split---->' + commande)

    # ouverture ficher resultat, suppression trois derniere colonne, ajout de 'rmsd' comme derniere colonne
    result = open('GD3_cl.E.info', 'r').readlines()
    resultfinal = open('AlphaFold/GalaxyD/' + e + '_result_galaxy.txt', 'w')
    resultfinal.write(result[0] + result[1][:-15] + 'RMSD\n' + result[0])

    # calcul RMSD
    liste_galaxy = liste_file(e, '/AlphaFold/GalaxyD/' + e)
    nbr_model = len(liste_galaxy)
    compteur = 0
    for model in liste_galaxy:
        commande = '/home/rene/bin/DockRMSD/DockRMSD LIGAND/MOL2/' + e + '.mol2 AlphaFold/GalaxyD/' + e + '/' + model +\
                   '.mol2 > rescore_result.txt'

        if terminal:
            os.system(commande)
            # else:
            print('rmsd ---> ' + commande)

        # extraction du RMSD
        rescore_lignes = open('rescore_result.txt', 'r').readlines()

        dockrmsd = 'NA'
        for i in rescore_lignes:
            if 'Calculated Docking RMSD:' in i:
                dockrmsd = i.split()[3]

        # ecriture du RMSD dans le fichier
        resultfinal.write(result[3 + compteur][:-27] + dockrmsd + '\n')
        compteur += 1

    resultfinal.close()


def dockG(liste_gal):
    c = 1
    for i in liste_gal:
        print('docking ' + str(c) + ' sur ' + str(len(liste_gal)))
        dock_galaxy(i)
        c += 1


def triage():
    l = liste_file('result_galaxy.txt', '/AlphaFold/GalaxyD')
    for s in l:
        print(s)
        file = result = open('AlphaFold/GalaxyD/' + s + '.txt', 'r').readlines()
        c = 0
        for i in file:
            c += 1
        print(c)
        # tail -14 1T56_ligand1 _result_galaxy.txt | sort -n -k 2
        co = 'head -3 AlphaFold/GalaxyD/' + s + '.txt > AlphaFold/GalaxyD/' + s + '_tried.txt'
        commande = 'tail -' + str(c - 3) + ' AlphaFold/GalaxyD/' + s + '.txt | sort -n -k 2 >> AlphaFold/GalaxyD/' + s\
                   + '_tried.txt'

        if terminal:
            os.system(co)
            os.system(commande)
        else:
            print(co)
            print(commande)


def ecritureligs2_3(structure):
    resultfile = open('AlphaFold/VINA/RESULT/comparaison_' + structure[:-1] + '1.txt', 'r')
    old = resultfile.readlines()
    resultfile.close()
    resultfile = open('AlphaFold/VINA/RESULT/comparaison_' + structure[:-1] + '1.txt', 'w')
    print(old)
    n=0
    for l in old:
        if l == old[0]:
            resultfile.write(l)
            n = 0
        elif l == old[1]:
            resultfile.write(l[:-1] + ' Docks RMSD lig' + structure[-1] + ' |\n')
        elif l == old[2]:
            resultfile.write(l[:-1] + ' Angstrom        |\n')
        else:
            n += 1

            # dockrmsd : on rescore le rmsd avec le ligand 2 ou 3 ; et on stocke le resultat dans un fichier temporaire
            # (sera ecrase a la prochaine boucle) nome rescore_result.txt
            commande = '/home/rene/bin/DockRMSD/DockRMSD LIGAND/MOL2/' + structure + \
                       '.mol2 AlphaFold/VINA/RESULT/model_' + structure[:-1] + '1/' + structure[:-1] + '1_model_' + \
                       str(n) + '.mol2 > rescore_result.txt'
            if terminal:
                os.system(commande)
            else:
                print(commande)

            # extraction du rescoring
            dockrmsd = 'NA'
            rescore_file = open('rescore_result.txt', 'r')
            rescore_lignes = rescore_file.readlines()
            for i in rescore_lignes:
                if 'Calculated Docking RMSD:' in i:
                    dockrmsd = i.split()[3]
            print(structure, dockrmsd)
            resultfile.write(l[:-1] + ' ' * 12 + dockrmsd + '\n')

    resultfile.close()


def ecritureligs2_3_GD(structure):
    resultfile = open('AlphaFold/GalaxyD/' + structure[:-1] + '1_result_galaxy.txt', 'r')
    old = resultfile.readlines()
    resultfile.close()
    resultfile = open('AlphaFold/GalaxyD/' + structure[:-1] + '1_result_galaxy.txt', 'w')
    print(old)
    n=0
    for l in old:
        if l == old[0]:
            resultfile.write(l)
            n = 0
        elif l == old[1]:
            resultfile.write(l[:-1] + ' Docks RMSD lig' + structure[-1] + ' |\n')
        elif l == old[2]:
            resultfile.write(l[:-1] + ' Angstrom        |\n')
        else:
            n += 1
            if n < 10:
                x = '0' + str(n)
            else:
                x = str(n)

            # dockrmsd : on rescore le rmsd avec le ligand 2 ou 3 ; et on stocke le resultat dans un fichier temporaire
            # (sera ecrase a la prochaine boucle) nome rescore_result.txt
            commande = '/home/rene/bin/DockRMSD/DockRMSD LIGAND/MOL2/' + structure + \
                       '.mol2 AlphaFold/GalaxyD/' + structure[:-1] + '1/' + structure[:-1] + '100' + x + \
                       '.mol2 > rescore_result.txt'
            if terminal:
                os.system(commande)
            else:
                print(commande)

            # extraction du rescoring
            dockrmsd = 'NA'
            rescore_file = open('rescore_result.txt', 'r')
            rescore_lignes = rescore_file.readlines()
            for i in rescore_lignes:
                if 'Calculated Docking RMSD:' in i:
                    dockrmsd = i.split()[3]
            print(structure, dockrmsd)
            resultfile.write(l[:-1] + ' ' * 12 + dockrmsd + '\n')

    resultfile.close()



if __name__ == '__main__':
    fichier_log = open('log_VINA_reDock.txt', 'a')
    fichier_log.write(str(datetime.date.today()) + ' -' * 40)

    dossier = 'AlphaFold'
    if terminal:
        if dossier not in listdir(path):
            os.system('mkdir ' + dossier)


    t = False
    alpha = 'model_alphaFold'
    #liste_lig = liste_file('.mol2', '/LIGAND/MOL2')
    #ecriture_pymol(alpha, '1T56')
    #conversion_recepteur_pdbqt(alpha)
    print('####################### Start #################################')



    liste_lig = liste_file('.pdbqt', '/LIGAND/ligand2_3')
    le = ['3Q0V_ligand2', '5MYM_ligand2', '5NZ1_ligand2', '5NZ1_ligand3']
    for i in liste_lig:
        if i not in le:
            print('lig ---------------------------> ' + i)
            ecritureligs2_3_GD(i)
    triage()



    if t:
        alpha = 'model_alphaFold'
        liste_lig = liste_file('.pdbqt', '/LIGAND/MOL2')
        ecriture_pymol(alpha, '1T56')
        conversion_recepteur_pdbqt(alpha)

        ########### adfr #########################
        grids(alpha)
        for lig in liste_lig:
            print('####################### lig ----> ' + lig)
            dockadfr(lig)
        check_summary()

        ########### vina #########################
        createBoxTxt(alpha)
        for lig in liste_lig:
            print('####################### lig ----> ' + lig)
            vina(lig, alpha)
            models = result_vina(lig)
            comparaison(lig, models)

        ########### Galaxy #########################
        dockG(liste_lig)
        triage()
