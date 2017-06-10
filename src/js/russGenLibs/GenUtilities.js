/**
 * @file General utility routines
 * @author RussellM on 12/08/2015
 */

/**
 * @module General Utilities module
 */

/* 19/4/16: Added ArrayMax and ArrayMin routines
 21/4/16: Added squishString
 4/3/17: Added binary search
 9/5/17: Added file jsdoc comment
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

function arrayMax(arr) {
    return arr.reduce(function (p, v) {
        return ( p > v ? p : v );
    });
}

function arrayMin(arr) {
    return arr.reduce(function(p,v) {
        return (p < v  ? p : v);
    });
}

function squishString(str,max) {
    var etc = '...';
    if (str.length > max) {
        return str.substring(0,max) + etc;
    }
    else {
        return str;
    }

}

//added 12/7/16
function getSelectedRadioEl(name) {
    //input: radio group name
    // output: element selected

    var rads = document.getElementsByName(name);

    var selected;

    for (var i = 0;i < rads.length;++i) {
        if (rads[i].checked) {
            selected = rads[i];
            break;

        }
    }
    return selected;
    
}


//added 4/3/17
function binarySearch(el,ar) {
    //assumes sorted array
    
    var N = ar.length;
    
    var bottom = 0;
    var top = N -1;
    var curr;
    
    var done = false;
    
    while (!done) {
        if (top == bottom) {
            if (el == ar[bottom]) {
                return bottom;
            }
            else {
                return -1;
            }
        }
        curr = bottom +  Math.floor((top - bottom) / 2);
        if (el == ar[curr]) {
            return curr + 1; // 1 based index
        }
        if (el > ar[curr]) {
            bottom = curr + 1;
            
        }
        else {
            top = curr;
        }
        
    }
    
    

}

//added 4/3/17
function insertionSort(ar,numSwapsFlag) {
    //if numSwapsFlag is true, also returns number of swaps required
    
    
    if (numSwapsFlag == null) {
        numSwapsFlag = false;
    }
    
    var totSwaps = 0;
    var swapsThisPass = 0;
    
    var done = false;
    
    while (!done) {
        swapsThisPass = 0;
        for (var i = 1;i < ar.length;++i) {
            if (ar[i] < ar[i-1]) {
                ++swapsThisPass;
                var tmp = ar[i];
                ar[i] = ar[i-1];
                ar[i-1] = tmp;
                
            }
        }
        totSwaps+= swapsThisPass;
        if (swapsThisPass == 0) {
            done = true;
        }
    }
    
    return numSwapsFlag ? [ar,totSwaps] : ar;
}

//added: 4/3/17

function majorityElement(ar) {

    var sortedAr = insertionSort(ar);

    //if there is a majority element, it must cross the middle
    var potential = sortedAr[Math.floor(sortedAr.length / 2) - 1];

    var st = -1;
    var end = -1;
    for (var i = 0;i < sortedAr.length;++i) {
        if (sortedAr[i] == potential) {
            if (st == -1) {
                st = i;
            }
        }
        else {
            if (st == -1) {

            }
            else {
                if (end == -1) {
                    end = i - 1;
                }
            }
        }

    }

    if (sortedAr[sortedAr.length -1] == potential) {
        end = sortedAr.length - 1;
    }

    var majElFound = false;
    if ((end - st + 1) > (sortedAr.length / 2) ) {
        majElFound = true;
        
    }
    

    return [majElFound,potential,st,end];

}

//added 4/3/17
function mergeSorted(ar1,ar2) {
    //assumes sorted arrays

    /* recursive:
    if (ar1.length == 0) {
        return ar2;
    }
    else if (ar2.length == 0) {
        return ar1;
    }

    var firstEl;
    if (ar1[0] < ar2[0]) {
        firstEl = ar1.shift();
    }
    else {
        firstEl = ar2.shift();
    }
    return [firstEl].concat(mergeSorted(ar1,ar2));
    */

    var leftProcessed = 0;
    var rightProcessed = 0;

    var done = false;

    var merged = [];

    while (!done) {
        if ((ar1.length == 0) && (ar2.length == 0)) {
            done = true;
            break;
        }
        else if (ar1.length == 0) {
            merged = merged.concat(ar2);
            done = true;
            break;
        }
        else if (ar2.length == 0) {
            merged = merged.concat(ar1);
            done = true;
            break;
        }

        if (ar1[0] < ar2[0]) {
            merged.push(ar1.shift());
        }
        else {
            merged.push(ar2.shift());
        }
    }

    return merged;

}

//added 5/3/17:
function mergeSort(ar) {

    /*
    if (ar.length == 1) {
        return ar;
    }
    
    var halfPoint = Math.ceil(ar.length / 2);
    var leftSide = ar.splice(0,halfPoint);
    
    var leftSorted = mergeSort(leftSide);
    var rightSorted = mergeSort(ar);
    return mergeSorted(leftSorted,rightSorted);
*/
    var activeArrays = [];
    ar.forEach(function(el) {
        activeArrays.push([el]);
    });
    while (activeArrays.length > 1) {
        if (activeArrays.length < 20) {
            var jj = 1;
        }
        var left = activeArrays.shift();
        var right = activeArrays.shift();
        var merged = mergeSorted(left,right);
        activeArrays.push(merged);

    }
    return activeArrays[0];
}

//added: 5/3/17:
function countInversions(ar) {
    var invs = 0;
    for (var i = 0;i < ar.length;++i) {
        for (var j = i;j < ar.length;++j) {
            if (ar[i] > ar[j]) {
                ++invs;
                
            }
        }
    }
    return invs;
    
}

//added 7/3/17:
function Heap(ar) {
    this.inputArray = ar;

    this.heapAr = [];
    
    this.swapElements = function(ind1,ind2) {
        var  tmp = this.heapAr[ind1];
        this.heapAr[ind1] = this.heapAr[ind2];
        this.heapAr[ind2] = tmp;
        
    };

    this.orderChildren = function(ch1,ch2) {
      //default largest val first, ie max heap
        return (this.heapAr[ch1] >  this.heapAr[ch2]) ? [ch1,ch2] : [ch2,ch1];
    };
    
    this.childrenInds = function(parentInd) {
        var childInd1 = (parentInd + 1) * 2 - 1; //zero based
        var childInd2 = (parentInd + 1) * 2;
        if ((childInd1 >= this.heapAr.length) && (childInd2 >= this.heapAr.length)) {
            //no children
            return [-1,-1];
        }
        else if (childInd1 >= this.heapAr.length) {
            return [childInd2,-1];
        }
        else if (childInd2 >= this.heapAr.length) {
            return [childInd1,-1];
        }
        else {
            return this.orderChildren(childInd1,childInd2);
            /*
            if (this.heapAr[childInd2] > this.heapAr[childInd1]) {
                return [childInd2,childInd1];
            }
            else {
                return [childInd1,childInd2];
            }
            */
        }
        
        
    };
    
  

    this.siftNeedToSwap = function(curInd,childInd) {
        //default max heap
        return this.heapAr[curInd] < this.heapAr[childInd];
    };


    this.siftDown = function() {
        var svdThis = this;
        var done = false;
        var curInd = 0;
        while (!done) {
            if (curInd >=this.heapAr.length) {
                done = true;
                break;
            }
            var childInds = this.childrenInds(curInd);
            if (childInds[0] == -1) {
                done = true;
                break;
            }
            else {
               // if (svdThis.heapAr[curInd] < this.heapAr[childInds[0]]) {
                if (this.siftNeedToSwap(curInd,childInds[0])) {
                    svdThis.swapElements(curInd, childInds[0]);
                    curInd = childInds[0];
                }
                else {
                    done = true;
                    break;
                }

            }
            

        }



    };

    this.addToSorted = function(sorted,el) {
        sorted.unshift(el);
    };

    this.sort = function(num) {
        //if num is populated, only return the first num elements

        if (num == null) {
            num = -1; //return all
        }

        var sorted = [];
        
        var svdHeap = this.heapAr.map(function(el) {
           return el; 
        });

        while (this.heapAr.length > 1) {
            this.swapElements(0, this.heapAr.length - 1);
            this.addToSorted(sorted,this.heapAr.pop());
            this.siftDown();
            if (num == -1) {

            }
            else if (sorted.length >= num) {
                break;

            }
        }

        if (num == -1) {
            sorted.unshift(this.heapAr.pop());
        }

        this.heapAr = svdHeap;

        return sorted;

    };

    this.bubbleNeedToSwap = function(curInd) {
        return this.heapAr[curInd] > this.parentVal(curInd);

    };

    this.bubbleUp = function() {
        var svdThis = this;
        var done = false;
        var curInd = svdThis.heapAr.length - 1;
        while (!done) {
            if (curInd <=0) {
                done = true;
                break;
            }
           // if (svdThis.heapAr[curInd] > svdThis.parentVal(curInd)) {
            if (svdThis.bubbleNeedToSwap(curInd)) {
                svdThis.swapElements(curInd,svdThis.parentInd(curInd));
                curInd = svdThis.parentInd(curInd);

            }
            else {
                done = true;
                break;
            }


        }

    };

    this.init = function() {

        var svdThis = this;
        this.inputArray.forEach(function(el) {
            svdThis.heapAr.push(el);
            svdThis.bubbleUp();


        });

    }
        
    this.parentInd = function(ind) {
            return Math.floor((ind+1)/2) - 1;
            
    };
        
    this.parentVal = function(ind) {
            return this.heapAr[this.parentInd(ind)];
            
    };
        
    
    

}

function MinHeap(ar) {

    Heap.apply(this, [ar]);

    this.orderChildren = function(ch1,ch2) {
        //small first
        return (this.heapAr[ch1] <  this.heapAr[ch2]) ? [ch1,ch2] : [ch2,ch1];
    };

    this.siftNeedToSwap = function(curInd,childInd) {
        return this.heapAr[curInd] > this.heapAr[childInd];
    };
    
    this.bubbleNeedToSwap = function(curInd) {
        return this.heapAr[curInd] < this.parentVal(curInd);
    };

    this.addToSorted = function(sorted,el) {
        sorted.push(el);
    };

}
//added 9/3/17:

function twoWayPartition(ar) {
    var first = ar[0];
    var higher = 0;
    for (var i = ar.length -1;i >=0;--i) {
        if (i % 1000 == 0) {
            var tst = 1;
        }
        if (ar[i] > first) {
            ++higher;
            //push to end
            ar.push(ar.splice(i,1)[0]);
            
        }
    }

}

//added 9/3/17:

function threeWayPartition(ar,v) {
    if (v == null) {
        v = ar[0];
    }

    var numVs = 0;

    var higher = 0;
    for (var i = ar.length -1;i >=0;--i) {
        if (i % 1000 == 0) {
            var tst = 1;
        }
        if (ar[i] > v) {
            ++higher;
            //push to end
            ar.push(ar.splice(i,1)[0]);

        }
        else if (ar[i] == v) {
            ++numVs;
            ar.splice(i,1);

        }
    }



    for (var i =0;i < numVs;++i) {
        ar.splice(ar.length - higher,0,v);
    }




     return [ar.length - higher - numVs,numVs,higher];


    
}
//added 9/3/17

function threeWayPartitionNonInPlace(ar,v) {
    if (v == null) {
        v = ar[0];
    }

    var numVs = 0;

    var smallAr = [];
    var largeAr = [];

    var higher = 0;
    for (var i = 0;i < ar.length;++i) {
        if (i % 1000 == 0) {
            var tst = 1;
        }
        if (ar[i] > v) {
            ++higher;
            //push to end
            largeAr.push(ar[i]);


        }
        else if (ar[i] == v) {
            ++numVs;


        }
        else {
            smallAr.push(ar[i]);
        }
    }

    var smallCount = smallAr.length;
    var largeCount = largeAr.length;

    ar.length = 0; // clear in place

    for (var i = 0;i < smallAr.length;++i) {
        ar.push(smallAr[i]);
    }

    for (var i =0;i < numVs;++i) {
        ar.push(v);
    }

    for (var i = 0;i < largeAr.length;++i) {
        ar.push(largeAr[i]);
    }

    ar = ar.concat(largeAr);




    return [smallCount,numVs,largeCount];



}

//added 11/3/17:
function quickSort(ar) {

    if (ar.length <= 1) {
        return ar;
    }

    var r = getRandomInt(0,ar.length-1);
    var v = ar[r];
    var inds = threeWayPartition(ar,v);

    var leftSorted = quickSort(ar.slice(0,inds[0]));
    var mid = ar.slice(inds[0],inds[0] + inds[1]);
    var rightSorted = quickSort(ar.slice(inds[0] + inds[1]),ar.length);
    
    var sorted = leftSorted.concat(mid);
    sorted = sorted.concat(rightSorted);
    
    return sorted;


}



//added 9/3/17:
function kthSmallest(ar,k) {
    if ((k > ar.length) || (k <= 0)) {
        return -1;
    }
    
    var r = getRandomInt(0,ar.length-1);
    var v = ar[r];
    var inds = threeWayPartition(ar,v);

    if (k <= inds[0]) {
        ar.splice(inds[0],ar.length-inds[0]);
        return kthSmallest(ar,k);

    }
    else if (k <= (inds[0] + inds[1])) {
        return v;
    }
    else {
        ar.splice(0,inds[0] + inds[1]);
        return kthSmallest(ar,k - (inds[0] + inds[1]));

    }


    
    
}

//added 12/3/17:
function deltaTurnpike(inAr) {
    var delta = [];
    for (var i = 0;i < inAr.length;++i) {
        for (var j = 0;j < inAr.length;++j) {
            //if (inAr[j] - inAr[i] > 0) {
            if (inAr[j] - inAr[i] == 1264) {
                var ggg = 1;
            }
                delta.push(inAr[j] - inAr[i]);
           // }
        }
    }
    delta.sort(function(a, b) {
        return a - b;
    });

    /*
    for (var i = 0;i < inAr.length;++i) {
        delta.unshift(0);
    }
    */
    
    return delta;
    
}

//added: 13/3/17:

function findBacktrack(turns,inputAr) {
    var currInd = turns.length - 1;
    for (var i = currInd;i >=0;--i) {
        var currTurn = turns[i];
        var currPossibles = currTurn[2];
        var currPossibleNums = currTurn[1];
        if ((currPossibles[0]) && (currPossibles[1])) { //had two choices
            currPossibles[0] = false;
            inputAr.pop();
            inputAr.push(currPossibleNums[1]);//now try second one
            break;
        }
        else { //only one possible, no other choice, keep trying backwards
            turns.pop();
            inputAr.pop();
            
        }
    }
    
    
}

function turnpikeBacktrack(turn) {

    var lenInputAr = Math.sqrt(turn.length);

    var done = false;

    var inputAr = [];
    var numDone = 0;

    turn = turn.filter(function(el) {
       return (el > 0);
    });

    turns = [];

    var newTurn = turn.map(function(el) {
        return el;
    });

    inputAr = [0];
    var absHighest = newTurn.pop();
    inputAr.push(absHighest);

    var turns = [];
    turns.push([newTurn,[absHighest]]);

    while (!done) {
        newTurn = newTurn.map(function(el) {
            return el;
        });

        for (var ii = 0;ii <inputAr.length - 1;++ii) {
            var diff = Math.abs(inputAr[inputAr.length-1] - inputAr[ii]);
            var ind = newTurn.indexOf(diff);
            if (ind == -1) {

            }
            else {
                newTurn.splice(ind,1);
            }
        }


        var highestInd = newTurn.length - 1;
        var highest = newTurn[highestInd];
        var other = absHighest - highest;
        var otherInd = newTurn.indexOf(other);

        var possibleNums = [other,highest];

        var possibles = [true,true];

        var impossible = false;

        if (otherInd == -1) {
            impossible = true;
            findBacktrack(turns,inputAr);
            newTurn = turns[turns.length-1][0].map(function(el) {
                return el;
            });
            
        }

        else {



            possibleNums.forEach(function (possEl,i) {
                var possible = true;
                inputAr.forEach(function (inEl) {
                    var ind = newTurn.indexOf(Math.abs(inEl - possEl));
                    if (ind == -1) {
                        possibles[i] = false;
                    }
                });
            });

            if ((possibles[0] == false) && (possibles[1] == false)) {
                impossible = true;
                findBacktrack(turns,inputAr);
                newTurn = turns[turns.length-1][0].map(function(el) {
                    return el;
                });
            }
            else {
                //newTurn.pop();
                //newTurn.splice(otherInd, 1);
                turns.push([newTurn,possibleNums,possibles]);
                for (var i = 0;i < possibleNums.length;++i) {
                    if (possibles[i]) {
                        inputAr.push(possibleNums[i]);

                        break;
                    }
                }
                if (inputAr.length >= lenInputAr) {
                    done = true;
                }
            }
        }


    }

    inputAr.sort(function(a,b) {
        return a-b;
    });
    
    return inputAr;

}

//added: 13/3/17

function turnpikeRecursive(turn) {

    if (turn.length == 3) {
        return [turn[0],turn[2]];
    }

    var newTurn = turn.filter(function(el) {
        return (el >= 0);
    });

    //var lessTurn =


    
}

//added: 11/3/17
function turnpike(turn) {

    var lenInputAr = Math.sqrt(turn.length);
    var inputArs = [];

    turn = turn.filter(function(el) {
        return el >= 0;
    });

    var inputAr = [];

   // turn.shift(); // remove zero diff to zero

  //  var highest = turn.pop();
    var highest = turn[turn.length - 1];

    inputAr.push(highest);

    inputAr.push(0);

    inputArs.push(inputAr);

   // turn.shift(); //remove zero diff to highest

    /*
    countDict = {};

       for (var i = 0;i < turn.length;++i) {
        if (turn[i] in countDict) {
            ++countDict[turn[i]];

        }

        else {
            countDict[turn[i]] = 1;
        }



    }

    var mustBeUnknownDiffs = [];

    var possibleNums = [];

    for (var key in countDict) {

        if (highest - parseInt(key) == parseInt(key) ) { //equal
            if (countDict[parseInt(key)] % 2 == 1) {
                //odd
                mustBeUnknownDiffs.push(parseInt(key));
            }


        }
        if (highest - parseInt(key) in countDict) {
            var numThisOne = countDict[key];
            var numComplement = countDict[highest - parseInt(key)];
            for (var i = 0; i < (numThisOne - numComplement); ++i) {
                mustBeUnknownDiffs.push(parseInt(key));

            }
        }

        else {
            for (var i = 0; i < countDict[key]; ++i) {
                mustBeUnknownDiffs.push(parseInt(key));
            }

        }
    }

    var mustBeUnknownDiffsWorking = mustBeUnknownDiffs.map(function(el) {
        return el;
    })

    */

    var done = false;

    var count = 0;

    var workingTurn;

    while (!done) {

        var newInputArs = [];
        inputArs.forEach(function(inputAr,inInd) {
            // if ((inputArs.length < 100) || (inInd % 5000 == 0)) {
            var delta = deltaTurnpike(inputAr);
            workingTurn = turn.map(function (el) {
                return el;

            });

            delta.forEach(function (el) {
                var ind = workingTurn.indexOf(el);
                workingTurn.splice(ind, 1);

            });

            var countDict = {};

            for (var i = 0; i < workingTurn.length; ++i) {
                if (workingTurn[i] in countDict) {
                    ++countDict[workingTurn[i]];

                }

                else {
                    countDict[workingTurn[i]] = 1;
                }


            }


            var lastEl = -1;
            var consistentCount = 0;

            for (var i = 0;i < workingTurn.length;++i) {

                el = workingTurn[i];// workingTurn.forEach(function (el, wInd) {
            // if ((workingTurn.length < 1000) || (wInd % 1000 == 0)) {

                if ((el == lastEl) || (inputAr.indexOf(el) > -1) || (el <= inputAr[inputAr.length - 1])) { // only put in greater nus

                }
                else {
                    lastEl = el;
                    var consistent = true;
                    inputAr.forEach(function (inEl) {
                        var diff = Math.abs(el - inEl);
                        var ind = workingTurn.indexOf(diff);
                        if (ind == -1) {
                            consistent = false;
                        }
                        else {
                            if ((diff == el) && (inEl > 0)) {
                                if (countDict[el] > 1) {

                                }
                                else {
                                    consistent = false;
                                }

                            }
                        }

                    });
                    if (consistent) {

                        var newAr = inputAr.map(function (inAr) {
                            return inAr;
                        });
                        newAr.push(el);
                        newInputArs.push(newAr);
                    }

                }
                // }

                if (newInputArs.length > 5000) {
                    break;
                }
            }//});

           // }


        });

        inputArs = newInputArs;

        ++count;
        console.log('candidates: ' +   inputArs.length +  ' workingturn: ' + workingTurn.length + ' count: ' + count + ' inparlen: ' + inputArs[0].length + ' targ: ' + lenInputAr + 'inputar0: ' + inputArs[0]);
        /*
        if (count > 100) {
            done = true;
        }
        */
        if (inputArs[0].length >= lenInputAr) {
            done = true;
        }
    }

    var possibleNums = turn.filter(function(el) {
        var ind = mustBeUnknownDiffsWorking.indexOf(el);
        if (ind == -1) {
            return true;
        }
        mustBeUnknownDiffsWorking.splice(ind,1);



    });

    var next = possibleNums[0];
    inputAr.push(next);
    inputAr.forEach(function(el) {
       var diff = Math.abs(el - next);
       var ind = turn.indexOf(diff);
       turn.splice(ind,1);
    });


    var smallestPoss = 1;

    var lastEl = -1;

    var nextPossibles = [];
    turn.forEach(function(el,i) {
        if (el == lastEl) {

        }
        else {
            lastEl = el;
            if (el >= smallestPoss) {
                var foundMatch = true;
                inputAr.forEach(function(inEl) {
                    var diff = Math.abs(inEl - el);
                    var ind = turn.indexOf(diff);
                    if ((ind == -1)) {
                        foundMatch = false;
                    }

                });
                if (foundMatch) {
                    nextPossibles.push(el);
                }

            }
        }

    });

    countDict = {};

    for (var i = 0;i < turn.length;++i) {
        if (turn[i] in countDict) {
            ++countDict[turn[i]];

        }

        else {
            countDict[turn[i]] = 1;
        }



    }




    turn.forEach(function(el) {

    });


}


//added: 4/3/17
function arrayToString(ar,sep) {
    var str = ''
    if (sep == null) {
        sep = ' ';
    }
    
    ar.forEach(function(el,i) {
        if (i == ar.length - 1) {
            str += el;
        }
        else {
            str += el + sep;
        }
        
    });
    
    return str;
    
}