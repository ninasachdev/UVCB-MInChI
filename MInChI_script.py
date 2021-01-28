#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
@author: Nina Sachdev
@date created: 1/13/21 (modified from my original JavaScript code created on 1/5/21)

'''
import json, rdkit, os, shutil
from rdkit import Chem


'''
@param fileName: name of the file containing UVCB data
@return filepath: the path of the file (including .txt)
@return filename: the name of the file (not including .txt)
'''
def getFile(fileName):
    filepath = fileName
    filename= fileName.split('.')[0]
    return filepath, filename

'''
Creates a list of dictionaries to hold all UVCB data
Each dictionary corresponds to one row (substance) in the file
Each key in the dictionary corresponds to a column in the file
Each value in the dictionary corresponds to the entry associated with that column in the file
Saves chem_dict as a JSON file
@param fileName: name of the file containing UVCB data (not including .txt)
@return chem_dict: list of dictionaries containing all UVCB data
@return chem_keys: list of column names in the file
'''
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

'''
Filters out all substances that do not have an InChI 
@param chem_dict: list of dictionaries containing all UVCB data
@return filter_InChI: list of dictionaries containing all substances that have an InChI
'''
def filterInChI(chem_dict):
    filter_InChI = []
    for i in range(len(chem_dict)):
        if len(chem_dict[i]['InChI_Final']) > 0:
            filter_InChI.append(chem_dict[i])

    return filter_InChI

'''
Extracts features needed to create the Mixfile and MInChI string
@param filter_InChI: list of dictionaries containing all substances that have an InChI
@return names: list of chemical names
@return InChI_list: list of InChI values
@return InChI_key_list: list of InChI keys
@return SMILES_list: list of SMILES values
'''
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

'''
Uses Chem module from rdkit library to convert a list of SMILES to their corresponding molfiles
@param SMILES_list: list of SMILES values
@return molfile_list: list of molfiles associated with the SMILES list
'''
def getMolfile(SMILES_list):
    molfile_list = []
    for i in range(len(SMILES_list)):
        m = Chem.MolFromSmiles(SMILES_list[i])
        molfile_string = Chem.MolToMolBlock(m)
        molfile_list.append(molfile_string) 

    return molfile_list


'''
Finds the concentration ratio for each substance, if it exists
Matches the corresponding unit (vp) to the ratio 
@param filter_InChI: list of dictionaries containing all substances that have an InChI
@param names: chemical names of all substances in filter_InChI
@return quantities: list of concentration ratios in filter_InChI, empty string if no ratio exists
@return units: list of units associated with concentration ratios, empty string if no ratio exists
'''
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
Constructs a Mixfile for each UVCB substance that has an InChI 
Saves the Mixfile in JSON format in a new /mixfiles/ directory
@param filename: name of the file containing UVCB data (not including .txt)
@param InChI_list: list of InChI values 
@param InChI_key_list: list of InChI keys
@param names: list of chemical names
@param molfiles: list of molfiles
@param SMILES: list of SMILES values
@param quantities: list of concentration ratios
@return mixfiles: list of dictionaries, where each dictionary corresponds to a substance's Mixfile
'''
def createMixfileDict(filename, InChI_list, InChI_key_list, names, molfiles, SMILES, quantities):
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

        with open('mixfiles/' + filename + '_mixture' + str(i) + '.mixfile', 'w') as outfile:
            json.dump(mixfile_dict, outfile)

        mixfiles.append(mixfile_dict)

    return mixfiles

'''
Builds the MInChI string from the substance's corresponding Mixfile
@param mixfiles: list of dictionaries, where each dictionary corresponds to a substance's Mixfile
@return minchis: a dictionary of MinChI strings, where each key is the chemical name and each value is the MInChI string
'''
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


'''
Adds the MInChI strings back into the original UVCB data
Creates a new MInChI column
Saves the updated data in a new .txt file 
@param minchis: dictionary of MInChI strings
@param chem_dict: list of dictionaries containing all UVCB data
@param chem_keys: list of column names in the original file
@param filename: name of file containing original UVCB data (not including .txt)
'''
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

'''
Computes various summary statistics for a given UVCB dataset
Finds how many substances had values for InChI, SMILES, concentration ratio, and a combination of the three
Writes summary statistics to a new file
@param chem_dict: list of dictionaries containing all UVCB data
@param minchis: dictionary of MInChI strings
@param filename: name of file containing all UVCB data (not including .txt)
'''
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
    
'''
Executes all functions in this script
For each substance in the UVCB dataset that contains an InChI, generates its Mixfile and MInChI string
Adds Mixfile to /mixfiles/ directory
Adds all MInChI strings to the dataset
Computes summary statistics for the dataset
@param currentFile: file path of the UVCB dataset (including .txt)
'''
def main(currentFile):
    filepath, filename = getFile(currentFile)

    chem_dict, chem_keys = createChemDict(filepath)

    filter_InChI = filterInChI(chem_dict)

    names, InChI_list, InChI_key_list, SMILES = getFeatures(filter_InChI)
    
    molfiles = getMolfile(SMILES)

    quantities, units = getConcentration(filter_InChI, names)

    mixfiles = createMixfileDict(filename, InChI_list, InChI_key_list, names, molfiles, SMILES, quantities)

    minchis = createMInChIString(mixfiles)

    addMInChIFile(minchis, chem_dict, chem_keys, filename)

    summary_stats(chem_dict, minchis, filename)


########################################
if __name__=='__main__':

    #List of UVCB datasets
    mixtureFiles = ['EU_REACH.txt', 'mixtures.txt', 'random.txt']

    #Current directory (should also contain UVCB datasets) 
    mainDir = os.getcwd()
    
    #Run the main function on the list of UVCB datasets
    for mixtureFile in mixtureFiles:
        os.chdir(mainDir)
        main(mixtureFile)

    
    

  
