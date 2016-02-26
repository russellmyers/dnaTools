/**
 * Created by RussellM on 12/08/2015.
 */

function addObserver(subject, property, callbackHandler) {
    subject['_' + property] = subject[property];
    Object.defineProperty(subject, property, {
        //Return the default value of the property
        //("value" automatically gives you the property's current value)
        get: function () {
            /// return value;
            return subject['_' + property];
        },

        //Set the property with a new value
        set: function (newValue) {
            //Assign the new value
            var oldVal = subject['_' + property];
            //value = newValue;
            subject['_' + property] = newValue;
            callbackHandler(subject,property,oldVal); //(oldVal);
            //Bind the observer's changeHandler to the subject
        },
        //Set the default parameters for how this property can be accessed and changed.
        //You probably don't need to change these unless you want to lock down the
        //property values to prevent your program from changing them
        enumerable: true,
        configurable: true,
        writeable: true
    });
}

function removeObserver(subject, property) {
    //Delete the changeHandler
   // delete subject.changeHandler;

    //Reset the getter and setter
    Object.defineProperty(subject, property, {
        get: function () {
            /// return value;
            return subject['_' + property];
        },
        set: function (newValue) {
            //value = newValue;
            subject['_' + property] = newValue;
        },
        enumerable: true,
        configurable: true,
        writeable: true
    });
}

function transpose(a) {
    //transpose array cols to rows

    return Object.keys(a[0]).map(
        function (c) { return a.map(function (r) { return r[c]; }); }
    );
}

//Random Routines

function normaliseProbDist(probDist) {
    //normalises probability distribution so adds up to 1
    var tot = probDist.reduce(function(a,b) {
        return a + b;

    });
    if (tot == 1.0) {
        return probDist;
    }

    var normed = probDist.map(function(el) {
       return el * 1.0 / tot;
    });

    return normed;

}

function getRandomFromProbDist(probDist) {

    var normed = normaliseProbDist(probDist);

    var cum = 0.0;
    var cumul = normed.map(function(el) {
        cum += el;
        return cum;

    });

    var r = getRandomArbitrary(0,1);

    for (var i = 0;i < cumul.length;++i) {
        if (r <= cumul[i]) {
            return i;
        }
    }
    return -1;
}




function getRandomArbitrary(min, max) {
    // Returns a random number between min (inclusive) and max (exclusive)
    return Math.random() * (max - min) + min;
}

function getRandomInt(min,max) {
    //returns a random int between min and max, both inclusive
    return Math.floor(Math.random() * (max - min + 1)) + min;

}

function alleleRandValueAlloc(numPerGene,initRange) {
    var alleleList = [];
    for (var i = 0;(i <numPerGene); ++i) {
        alleleList.push(getRandomArbitrary(initRange * -1,initRange));
    }
    return alleleList;
}


Date.prototype.yyyymmdd = function() {
    var yyyy = this.getFullYear().toString();
    var mm = (this.getMonth()+1).toString(); // getMonth() is zero-based
    var dd  = this.getDate().toString();
    return yyyy + (mm[1]?mm:"0"+mm[0]) + (dd[1]?dd:"0"+dd[0]); // padding
};


function glowEnded(el) {
  //  setTimeout(function() {
        el.style.animation = null;
        console.log('glow ended: ' + el.innerHTML);
        el.offsetWidth = el.offsetWidth;


  //  },10);

}

function glowAnimation(el,anCSS,backTime) {
    el.addEventListener("animationend", glowEnded(el),false);
    setTimeout(function() {
        el.style.animation = anCSS;

    },10);

    /*
    setTimeout(function() {
        el.style.animation = '';
     },backTime);
     */


}

function getKeys(obj){
    var keys = [];
    for (var key in obj) {
        if (obj.hasOwnProperty(key)) { keys[keys.length] = key; }
    }
    return keys;
}

function setUpFileReader(fileInput,fileExtVal,readCompleteCallback,fileDisplayArea) {

    if (window.File && window.FileReader && window.FileList && window.Blob) {
        // Great success! All the File APIs are supported.

    } else {
        alert('The File APIs are not fully supported in this browser.');
    }

    fileInput.addEventListener('change', function (e) {
        var file = fileInput.files[0];

        var fileExt = file.name.split('.').pop();

        //var textType = /text.*/;
        //var csvType = /csv.*/;
        //if ((file.type.match(textType)) || (file.type.match(csvType))) {
        if (fileExt == 'json') {// (file.type.match('application/json')) {

            var reader = new FileReader();

            reader.onload = function (e) {

                fileContents = reader.result;
                if (fileDisplayArea) {
                    fileDisplayArea.innerText = fileContents;
                }

                readCompleteCallback(fileContents);
            };
            reader.readAsText(file);
        } else {
            fileDisplayArea.innerText = "File not supported!"
        }
    });
}
function createStudDropDown(ddId,inclNone) {

    inclNone = inclNone || false;

    var studDropDown = document.createElement("select");
    studDropDown.id = ddId;

    var optGrp = document.createElement("optgroup");
    studDropDown.appendChild(optGrp);

    var option;



    if (inclNone) {
        option = document.createElement("option");
        option.text = "  ";
        //studDropDown.add(option);
        optGrp.appendChild(option);
    }
    gPool.genes.forEach(function(g) {
        option = document.createElement("option");
        option.text = g.name;
        //studDropDown.add(option);
        optGrp.appendChild(option);

    });

    return studDropDown;

}

function createStudOptGrp(inclNone) {
    var option;

    var optGrp = document.createElement("optgroup");


    if (inclNone) {
        option = document.createElement("option");
        option.text = "  ";
        //studDropDown.add(option);
        optGrp.appendChild(option);
    }
    gPool.genes.forEach(function(g) {
        option = document.createElement("option");
        option.text = g.name;
        //studDropDown.add(option);
        optGrp.appendChild(option);

    });

    return optGrp;

}

function updateStudDropDown(cr,inclNone) {

    inclNone = inclNone || false;

    cr.removeChild(cr.childNodes[0]);

    var optGrp = document.createElement("optgroup");

    cr.appendChild(optGrp);


    var option;



    if (inclNone) {
        option = document.createElement("option");
        option.text = "  ";
        //studDropDown.add(option);
        optGrp.appendChild(option);
    }
    gPool.genes.forEach(function(g) {
        option = document.createElement("option");
        option.text = g.name;
        //studDropDown.add(option);
        optGrp.appendChild(option);

    });


}


function stripTags(inTag) {
    var startedTag = false;
    var outTag = '';


    if (inTag) {

    }
    else {
        return ['','none'];
    }


    for (var i = 0; i < inTag.length; ++i) {
        var lastOne = false;
        if (inTag[i] == '<') {
            startedTag = true;
        }
        if (inTag[i] == '>') {
            startedTag = false;
            lastOne = true;
        }
        if (startedTag) {

        }
        else {
            if (lastOne) {

            }
            else {
                outTag += inTag[i];
            }
        }

    }

    var howGood;
    //console.log('intag: ' + inTag);
    if (inTag.search('red') != -1) {
        howGood = 'excl';

    }
    else if (inTag.search('green') != -1) {
        howGood = 'good';
    }
    else  {
        howGood = 'nogood';
    }
    return [outTag,howGood];
}

function updateSpecialTable(table,tableData) {
   // console.log('xxx   start updatespec table');



   // table.rows.forEach(function(row,i) {
      //  row.cells.forEach(function(cell,j) {
    for (var i = 0;i < table.rows.length; ++i) {
      //  console.log('update spec: ' + i);
        if (i == 0) continue; //header

        var row = table.rows[i];
        for (var j = 0; j < row.cells.length; ++j) {
           // console.log('i: ' + i + ' j: ' + j + ' row cells: ' + row.cells.length);
            var cell = row.cells[j];
            if (j == 0)  {
                cell.childNodes[0].value = tableData[i-1][j];

                //cell.innerHTML = tableData[i-1][j];
            }
            else {

                var sel = cell.children[0];
               // console.log(cell.id);
               // var cr = createStudDropDown(cName + '_' + r + '_' + cNum,true);
                //cr.value = 3;
                //var optGrp =  createStudOptGrp(true);
                /* restore this to update if students change
                updateStudDropDown(sel,true);
                */
                var stripped = stripTags(tableData[i-1][j]);
                //cr.selectedIndex = stripTags(cellData)[0] - 1;
                sel.selectedIndex = stripped[0];
                sel.className = stripped[1];
            }
        }
    }
    //var gg = 0;

    //console.log('xxx   end update spec table');


}

function studNameChanged() {
    console.log('changed name of ' + this.value + ' id: ' + this.id);
    var spl = this.id.split('_');
    var g = gPool.geneWithIndex(spl[1]);
    g.name = this.value;

    updateDisplay();
    changeDisplayedIndiv();


}

function createSpecialTable(tableData,headerData,tabType) {

    tabType = tabType || 'fixed';
    headerData = headerData || null;

    var table = document.createElement('table')
        , tableBody = document.createElement('tbody');
    table.style.tableLayout = tabType;

    var row = document.createElement('tr');
    if (headerData) {
        headerData.forEach(function(head) {
            var cell = document.createElement('td');

            // cell.appendChild(document.createTextNode(cellData));
            cell.innerHTML = head;
            row.appendChild(cell);

        });
        tableBody.appendChild(row);
    }

    tableData.forEach(function(rowData,r) {
        var row = document.createElement('tr');

        rowData.forEach(function(cellData,c) {
            var cell = document.createElement('td');
            //cell.className = 'studDropCell';
            // cell.appendChild(document.createTextNode(cellData));
            if (c > 0) {
                var cNum = c - 1;
                var cName = 'studPref';

                if (cNum > StudentGene.maxPrefs - 1) { //must be an excl
                    cNum -= StudentGene.maxPrefs;
                    cName = 'studExcl';
                }
                var cr = createStudDropDown(cName + '_' + r + '_' + cNum,true);
                //cr.value = 3;
                if (!cellData) {
                   //alert('not cell data r: ' + r + ' c: ' + c);
                }
                var stripped = stripTags(cellData);
                //cr.selectedIndex = stripTags(cellData)[0] - 1;
                cr.selectedIndex = stripped[0];
                cr.className = stripped[1];
                cr.onchange = function() {
                    //alert('changed' + this.id);
                    var spl = this.id.split('_');

                    var g = gPool.geneWithIndex(spl[1]);
                    if (spl[0] == 'studPref') {
                        g.setPrefWithIndex(spl[2],'' + this.selectedIndex);

                    }
                    else {
                        g.setExclWithIndex(spl[2],'' + this.selectedIndex);
                    }

                    updateDisplay();
                    changeDisplayedIndiv();

                };
                cell.appendChild(cr);
                row.appendChild(cell);


            }
            else {
                /*
                var nameDiv = document.createElement('div');
                nameDiv.appendChild(cell);
                nameDiv.width = 50;
                nameDiv.overflow = 'hidden';
                //cell.appendChild(nameDiv);
                cell.innerHTML = cellData;
                cell.contentEditable = true;
                cell.overflow = 'hidden';
                row.appendChild(nameDiv);
                */
                var nameInp = document.createElement('input');
                nameInp.id = 'nameInput_' + r;
                nameInp.className = 'nameInputClass';
                nameInp.value = cellData;
                nameInp.onchange = studNameChanged;
                cell.appendChild(nameInp);
                //nameInp.appendChild(cell);
                //cell.innerHTML = cellData;
                //cell.contentEditable = true;
                row.appendChild(cell);

            }


        });

        tableBody.appendChild(row);
    });

    table.appendChild(tableBody);

    // document.body.appendChild(table);
    return table;
}

function createTable(tableData,headerData,editableCols,startEditableCol,keyPressCallback,onBlurCallback,onFocusCallback,tabIdPrefix,tabType) {
    tabType = tabType || 'fixed';
    editableCols = editableCols || 0;
    startEditableCol = startEditableCol || null;
    keyPressCallback = keyPressCallback || null;
    onBlurCallback = onBlurCallback || null;
    onFocusCallback = onFocusCallback || null;
    tabIdPrefix = tabIdPrefix || '';
    headerData = headerData || null;

    var table = document.createElement('table')
        , tableBody = document.createElement('tbody');
    table.style.tableLayout = tabType;

    var row = document.createElement('tr');
    row.className += ' createdTableHeaderRow';
    if (headerData) {
        headerData.forEach(function(head) {
            var cell = document.createElement('th');
            // cell.appendChild(document.createTextNode(cellData));
            cell.innerHTML = head;
            row.appendChild(cell);

        });
        tableBody.appendChild(row);
    }

    tableData.forEach(function(rowData,rowNum) {
        var row = document.createElement('tr');
        if (rowNum % 2 == 0) {
            row.className += ' createdTableEvenRow';
        }
        else {
            row.className += ' createdTableOddRow';
        }

        rowData.forEach(function(cellData,colNum) {
            var cell = document.createElement('td');
            cell.id =  tabIdPrefix + rowNum + 'c' + colNum;
            if ((colNum >= startEditableCol)
                && (colNum < (startEditableCol + editableCols))) {
                cell.contentEditable = true;
                cell.onkeypress = function(event) {
                    keyPressCallback(event);
                };
                cell.onblur = function(event) {
                    onBlurCallback(event);
                };
                cell.onfocus = function(event) {
                    onFocusCallback(event);
                };


            }
            // cell.appendChild(document.createTextNode(cellData));
            cell.innerHTML = cellData;
            row.appendChild(cell);
        });

        tableBody.appendChild(row);
    });

    table.appendChild(tableBody);

    // document.body.appendChild(table);
    return table;
}



function checkEnter(event)
{
    var keynum;
    var keychar;
    var enttest;

    if (window.event)
    {
        keynum = event.keyCode;
    }
    else
    {
        if (event.which)
        {
            keynum = event.which;
        }
    }

    keychar = String.fromCharCode(keynum);
    enttest = "\r";

    if (enttest == keychar)
    {
        event.srcElement.blur();
        return false;
    }
}
//Tab stuff
function tab_click(x,onComplete){
    var myParams = Params.getInstance();
    var tabActive=myParams.tabActive;
    var tabs=document.getElementById('tabs').getElementsByTagName('A');
   // var tabs_data=document.getElementById('tabs_data').getElementsByTagName('fieldset');
    var tabs_data=document.getElementById('tabs_data').getElementsByTagName('section');

    if(x > -1 && x < tabs.length && x < tabs_data.length) {
        tabs[tabActive].setAttribute('class', '');
        tabs_data[tabActive].style.display = 'none';
        // tabs_data[tabActive].style.opacity = 0;
        //tabs_data[tabActive].style.width = 0;
        //tabs_data[tabActive].style.height = 0;


        tabActive = x;
        tabs[tabActive].setAttribute('class', 'active');
        tabs_data[tabActive].style.display = 'block';
        //tabs_data[tabActive].style.opacity = 1;
        //tabs_data[tabActive].style.width = 300;
        // tabs_data[tabActive].style.height = 300;


        myParams.tabActive = tabActive;
        if (onComplete) {
           onComplete(tabActive);
        }

    } return false;
}
/*
document.onkeydown=function (evnt){
    if(evnt.keyCode==37 || evnt.keyCode==39)
        tab_click(tabActive+evnt.keyCode-38);
}
*/
//End of tab stuff

//File creation stuff
makeTextFile = function (text) {
    var data = new Blob([text], {type: 'text/plain'});

    // If we are replacing a previously generated file we need to
    // manually revoke the object URL to avoid memory leaks.
    if (textFile !== null) {
        window.URL.revokeObjectURL(textFile);
    }

    var textFile = window.URL.createObjectURL(data);

    return textFile;
};
