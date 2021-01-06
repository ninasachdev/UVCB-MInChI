var filepath;
function getFile(fileName) {
    filepath = fileName;
}



getFile('EU_REACH.txt');

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
            fs.writeFile('fileArray.json', arrayToJSON, function(err) {
                if (err) {
                    console.log(err);
                }
            })
        } 
    }
    , 2000);
}


getFileContent();
