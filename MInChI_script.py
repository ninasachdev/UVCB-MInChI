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

def getFeatures(filter_InChI):
    names = []
    InChI_list = []
    InChI_key_list = []
    SMILES_list = []
    for i in range(len(filter_InChI)):
        names.append(filter_InChI[i]['Chemical_Name_Final'])
        InChI_list.append(filter_InChI[i]['InChI_Final'])
        InChI_key_list.append(filter_InChI[i]['InChIKey_Final'])
        SMILES_list.append(filter_InChI[i]['SMILES_Final'])

    return names, InChI_list, InChI_key_list, SMILES_list

def getMolfile(SMILES_list):
    molfile_list = []
    for i in range(len(SMILES_list)):
        m = Chem.MolFromSmiles(SMILES_list[i])
        molfile_string = Chem.MolToMolBlock(m)
        molfile_list.append(molfile_string) 

    return molfile_list

#all mixtures or just the ones that have an InChI??

def getConcentration(filter_InChI, names):
    quantities = ['' for i in range(len(names))]
    units = ['' for i in range(len(names))]
    
    for n in range(len(names)):
        name = names[n]
        for l in range(len(name)):
            if name[l] == '(':
                if name[l+2] == ':' and name[l+4] == ')':
                    quantities[n] = name[l:l+5]
                    units[n] = 'vp' #MInChI units for ratio
        
    return quantities, units




'''
Note that ratio units are handled as a special case,
since the Mixfle stores them in [numerator, denominator] form, while in the MInChI notation, only the numerators are listed, and the denominator is implied
'''
#technically don't need the units parameter?

def createMixfileDict(filename, InChI_list, InChI_key_list, names, molfiles, SMILES, quantities, units):
    mixfiles = []
    for i in range(len(InChI_list)):
        mixfile_dict = {}
        mixfile_dict['mixfileVersion'] = 0.01
        mixfile_dict['name'] = names[i]
        mixfile_dict['molfile'] = molfiles[i]
        mixfile_dict['inchi'] = InChI_list[i]
        mixfile_dict['inchiKey'] = InChI_key_list[i]
        mixfile_dict['smiles'] = SMILES[i]
        if len(quantities[i]) > 0:
            ratio = '[' + str(quantities[i][1]) + ', ' + str(quantities[i][3] + ']')
            mixfile_dict['ratio'] = ratio

        outfile = open('mixfiles/' + filename + 'mixture' + str(i) + '.mixfile', 'w')
        outfile.write('{"mixfileVersion":' + str(mixfile_dict['mixfileVersion']) + ',\n')
        outfile.write('"name":' + str(mixfile_dict['name']) + ',\n')
        outfile.write('"molfile":' + str(mixfile_dict['molfile']) + ',\n')
        outfile.write('"inchi":' + str(mixfile_dict['inchi']) + ',\n')
        outfile.write('"inchiKey":' + str(mixfile_dict['inchiKey']) + ',\n')
        outfile.write('"smiles":' + str(mixfile_dict['smiles']) + ',\n')

        if len(quantities[i]) > 0:
            ratio = quantities[i]
            mixfile_dict['ratio'] = ratio
            outfile.write('"ratio":' + '[' + str(quantities[i][1]) + ', ' + str(quantities[i][3]) + ']' + ',\n')

        outfile.write('}')
        outfile.close()

        mixfiles.append(mixfile_dict)

    return mixfiles

def createMInChIString(mixfiles):
    minchis = {}
    for mixfile in mixfiles:
        header = 'MInChI=0.00.1S'
        identifier = mixfile['inchi'][9:]
        indexing = 'n1'
        concentration = ''
        if 'ratio' in mixfile.keys():
            concentration = mixfile['ratio'][1] + ':' + mixfile['ratio'][3] + 'vp'

        minchi = header + '/' + identifier + '/' + indexing + '/' + concentration
        minchis[mixfile['name']] = minchi
    
    return minchis


########################################
if __name__=='__main__':

    mixtureFiles = ['EU_REACH.txt', 'mixtures.txt', 'random.txt'] #want to eventually loop through and compute for all three files


    filepath, filename = getFile('EU_REACH.txt')

    chem_dict = createChemDict(filepath)
    filter_InChI = filterInChI(chem_dict)
    names, InChI_list, InChI_key_list, SMILES = getFeatures(filter_InChI)
    molfiles = getMolfile(SMILES)

    #print(names)
    #print(InChI_list[2])
    #print(len(InChI_key_list))
    #print(len(SMILES))
    #print(molfiles)


    #m = Chem.MolFromSmiles(SMILES[3])
    #print(rdkit.Chem.inchi.MolToInchi(m))

    quantities, units = getConcentration(filter_InChI, names)
    #print(quantities)
    #print(len(units))

    mixfiles = createMixfileDict(filename, InChI_list, InChI_key_list, names, molfiles, SMILES, quantities, units)

    print(createMInChIString(mixfiles))