
var mixtureFiles = ['EU_REACH.txt', 'mixtures.txt', 'random.txt'];

var fileArray;

var fs = require('fs'); 

fileArray = JSON.parse(fs.readFileSync('EU_REACHArray.json')); //eventually need to add for loop to include all text files

//console.log(fileArray[0][2]);
//to create Mixfile: need name, molfile, quantity, units, inchi, inchiKey

//create a list of dictionaries
var chem_dict = []
var chem_keys = fileArray[0]
for (i = 1; i < filerray.length; i++) { //want to skip the first row of labels
    current_dict = {}
    for (k=0; k < chem_keys.length; k++) {
        current_dict[chem_keys[k]] = fileArray[i][k];
    }
    chem_dict.push(current_dict);
}

//console.log(chem_dict[50]);

filter_InChI = [];
for (i = 0; i < chem_dict.length; i++) {
    if (chem_dict[i]['InChI_Final'].length > 0) {
        filter_InChI.push(chem_dict[i]);
    }
}

console.log(filter_InChI.length); //there are only 26 entries in EU_REACH that have an InChI value