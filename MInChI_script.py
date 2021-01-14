#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################
#@author: Nina Sachdev
#@date created: 1/13/21

###########################
import json, rdkit
from rdkit import Chem


def getFile(fileName):
    filepath = fileName
    filename= fileName.split('.')[0]
    return filepath, filename


def createChemDict(fileName):
    filepath, filename = getFile(fileName)

    with open (filepath, 'r') as infile:
        lines = infile.readlines()
        chem_keys =  lines[0].split('\t')  

    chem_dict = []
    for i in range(1, len(lines)): #want to skip the first row of labels
        current_dict = {}
        for k in range(len(chem_keys)):
            current_dict[chem_keys[k]] = lines[i].split('\t')[k]

        chem_dict.append(current_dict)

    with open(filename + '_dict.json', 'w') as outfile:
        json.dump(chem_dict, outfile)

    return chem_dict


def filterInChI(chem_dict):
    filter_InChI = []
    for i in range(len(chem_dict)):
        if len(chem_dict[i]['InChI_Final']) > 0:
            filter_InChI.append(chem_dict[i])

    return filter_InChI

def getInChI(filter_InChI):
    InChI_list = []
    InChI_key_list = []
    for i in range(len(filter_InChI)):
        InChI_list.append(filter_InChI[i]['InChI_Final'])
        InChI_key_list.append(filter_InChI[i]['InChIKey_Final'])

    return InChI_list, InChI_key_list

def getSMILES(filter_InChI):
    SMILES_list = []
    for i in range(len(filter_InChI)):
        SMILES_list.append(filter_InChI[i]['SMILES_Final'])

    return SMILES_list

def getMolfile(SMILES_list):
    molfile_list = []
    for i in range(len(SMILES_list)):
        m = Chem.MolFromSmiles(SMILES_list[i])
        molfile_string = Chem.MolToMolBlock(m)
        molfile_list.append(molfile_string)

    return molfile_list

#eventually need to add quantity and units
#def createMixfileDict(filter_InChI, InChI_list, InChI_key_list, SMILES, molfiles):


#don't forget to remove prefix in front of InChI string!
#def createMInChIString


########################################
if __name__=='__main__':

    mixtureFiles = ['EU_REACH.txt', 'mixtures.txt', 'random.txt'] #want to eventually loop through and compute for all three files


    filepath, filename = getFile('random.txt')

    chem_dict = createChemDict(filepath)
    filter_InChI = filterInChI(chem_dict)
    InChI_list, InChI_key_list = getInChI(filter_InChI)
    SMILES = getSMILES(filter_InChI)
    molfiles = getMolfile(SMILES)

    print(InChI_list)
    print(len(InChI_key_list))
    print(len(SMILES))
    print(len(molfiles))