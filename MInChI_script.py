#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
@author: Nina Sachdev
@date created: 1/13/21 (modified from my original JavaScript code created on 1/5/21)

'''
import json, rdkit, os, shutil
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

    infile.close()

    if os.path.isdir(filename + '_MInChI'):
        shutil.rmtree(filename + '_MInChI')
        
    if os.path.isdir(filename + '_MInChI/mixtures'):
        shutil.rmtree(filename + '_MInChI/mixtures')

    os.mkdir(filename + '_MInChI')
    os.chdir(filename + '_MInChI')

    with open(filename + '_dict.json', 'w') as outfile:
        json.dump(chem_dict, outfile)

    return chem_dict, chem_keys


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

def createMixfileDict(filename, InChI_list, InChI_key_list, names, molfiles, SMILES, quantities, units):
    os.mkdir('mixfiles')

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


        '''
        outfile = open('mixfiles/' + filename + '_mixture' + str(i) + '.mixfile', 'w')
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

        '''

        with open('mixfiles/' + filename + '_mixture' + str(i) + '.mixfile', 'w') as outfile:
            json.dump(mixfile_dict, outfile)

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

def addMInChIFile(minchis, chem_dict, chem_keys, filename):

    outfile = open(filename + '_MInChI.txt', 'w')

    chem_keys[-1] = 'Iodo'
    chem_keys.append('MInChI\n')
    outfile.write('\t'.join(chem_keys))

    for i in range(len(chem_dict)):
        name = chem_dict[i]['Chemical_Name_Final']
        elts = list(chem_dict[i].values())
        elts[-1] = elts[-1].strip('\n')
        if name in minchis:
            elts.append(minchis[name])
        else:
            elts.append('')

        outfile.write('\t'.join(elts) + '\n')


def summary_stats(chem_dict, minchis, filename):
    stats = {}
    for i in range(len(chem_dict)):
        keys = []
        for key, val in chem_dict[i].items():
            if len(chem_dict[i][key]) > 0 and (key == 'InChI_Final' or key == 'SMILES_Final'):
                keys.append(key)

        name = chem_dict[i]['Chemical_Name_Final']
        for l in range(len(name)):
            if name[l] == '(':
                if name[l+2] == ':' and name[l+4] == ')':
                    keys.append('Ratio_Concentration')


        stats[name] = keys

    inchi_smiles_ratio = inchi_smiles = inchi_ratio = smiles_ratio = inchi = smiles = ratio = 0

    for s in stats:
        if stats[s] == ['InChI_Final', 'SMILES_Final', 'Ratio_Concentration']:
            inchi_smiles_ratio += 1
        elif stats[s] == ['InChI_Final', 'SMILES_Final']:
            inchi_smiles += 1
        elif stats[s] == ['InChI_Final', 'Ratio_Concentration']:
            inchi_ratio += 1
        elif stats[s] == ['SMILES_Final', 'Ratio_Concentration']:
            smiles_ratio += 1
        elif stats[s] == ['InChI_Final']:
            inchi += 1
        elif stats[s] == ['SMILES_Final']:
            smiles += 1
        elif stats[s] == ['Ratio_Concentration']:
            ratio += 1

    outfile = open(filename + 'summary_stats.txt', 'w')
    labels = ['Num_MInChIs', 'InChI_SMILES_Ratio', 'InChI_SMILES', 'InChI_Ratio', 'SMILES_Ratio', 'InChI', 'SMILES', 'Ratio']
    stats_list = [len(minchis), inchi_smiles_ratio, inchi_smiles, inchi_ratio, smiles_ratio, inchi, smiles, ratio]

    for i in range(len(labels)):
        outfile.write(labels[i] + ': ' + str(stats_list[i]) + ' UVCB substances')
        outfile.write('\n')
    

def main(currentFile):
    filepath, filename = getFile(currentFile)

    chem_dict, chem_keys = createChemDict(filepath)

    filter_InChI = filterInChI(chem_dict)

    names, InChI_list, InChI_key_list, SMILES = getFeatures(filter_InChI)
    
    molfiles = getMolfile(SMILES)

    quantities, units = getConcentration(filter_InChI, names)

    mixfiles = createMixfileDict(filename, InChI_list, InChI_key_list, names, molfiles, SMILES, quantities, units)

    minchis = createMInChIString(mixfiles)

    addMInChIFile(minchis, chem_dict, chem_keys, filename)

    summary_stats(chem_dict, minchis, filename)


########################################
if __name__=='__main__':

    mixtureFiles = ['EU_REACH.txt', 'mixtures.txt', 'random.txt']

    mainDir = os.getcwd()
    
    for mixtureFile in mixtureFiles:
        os.chdir(mainDir)
        main(mixtureFile)

    
    

  
