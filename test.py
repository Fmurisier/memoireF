
import os
from os import listdir
from os.path import isfile, join

path = os.getcwd()
terminal = False

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
            elif line[0] == 'HETATM' and line[3 + d] != 'HOH':
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


def ecriture_box():
    dossier = 'residus.txt'
    if dossier in listdir(path):

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
        print('La box est composee de ' + str(len(resi_liste)) + \
              ' residus, si ce n\'est pas le cas verifiez que le fichier residus.txt soie ecrit correctement (' \
              ' chaque residu doit etre separe par un +)')

        # commande = 'select ' + name + '_poche, resi ' + '+'.join(resi_liste) + ' and model ' + name
        commande = '_poche, resi ' + '+'.join(resi_liste) + ' and model '
        return commande
    else:
        print('error file residus.txt don\'t exist ! ')
        return 'error'


if __name__ == '__main__':

    l = '1T56'
    ecriture_box(l)

    '''
    select 1T56_poche, resi 87+90+102+103+107+110+114+138+142+145+148+149+152+176+179+180+183+184+207 and model 1T56

    select 1U9N_poche, resi 87+90+102+103+107+110+114+138+142+145+148+149+152+176+179+180+183+184+207 and model 1U9N

    super 1T56_poche////CA,1U9N_poche////CA
'''


