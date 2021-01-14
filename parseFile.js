var filepath;
var filename;
function getFile(fileName) {
    filepath = fileName;
    filename = fileName.split('.')[0];
}

getFile('random.txt'); //change file name accordingly

const lineReader = require('line-reader'); //load line-reader module installed from npm
var rawFileContent = [];
lineReader.eachLine(filepath, function(line, last) {
    rawFileContent.push(line)
    if(last) {
        return false;
    }
}); 

var fileContent;
var fileArray = [];
function getFileContent() {
    setTimeout(function () {
        fileContent = rawFileContent;
        for (i=0; i < fileContent.length; i++) {
            fileArray.push(fileContent[i].split('\t'));
            var arrayToJSON = JSON.stringify(fileArray);
            var fs = require('fs');
            fs.writeFile(filename + 'Array.json', arrayToJSON, function(err) {
                if (err) {
                    console.log(err);
                }
            })
        } 
    }
    , 8000); //might need to increase timeout length if there are errors in the resulting JSON file
}


getFileContent();
