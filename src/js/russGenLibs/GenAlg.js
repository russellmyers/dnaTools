/**
 * Created by RussellM on 29/07/2015.
 */

//Singleton
function Params() {
    // if (Params.instance instanceof 'object') {
    //     return Params.instance;
    // }
    this._indivToDisplay = 1;
    this.indivToDisplay = 1;
    this.bestIndivToDisplay = -1;
    this.mutRate = 0.25;
    this.popSize = 50;
    this.crossoverRate = 1.0;
    this.numToGenAtOnce = 50;
    this.eliteFlag = true;
    this.catastropheGens = 500;
    this.isDirty = false;
    this.classTolerance = 1;

    this.tabActive = 0; //tabs stuff

    this.numGensNoProgress = 0; //number of gens where high score not increased

    this.trialClass = null;
    this.trialBest = null;

    this.localSavedExists = false; //file saved in local storage

    Params.instance = this;


}

Params.getInstance = function() {

    if (typeof Params.instance === 'object') {
        return Params.instance;
    }
    return new Params();


};

function Gene(num) {
    this.alleles = {};
    this.num = num;
    this.nextAlleleNum = num*1000000 + 1;
    this.info = 'Gen';

}

Gene.prototype.addAllele = function(alleleVal) {
    this.alleles[this.nextAlleleNum] = alleleVal;
    ++this.nextAlleleNum;
    return this.nextAlleleNum - 1;

};

Gene.prototype.numAlleles = function()    {
    var size = 0, key;
    for (key in this.alleles) {

        if (this.alleles.hasOwnProperty(key)) size++;
    }
    return size;

};

Gene.prototype.alleleKey = function(alleleVal) {
    //var keys = getKeys(this.alleles);
    for (var key in this.alleles) {
        if (this.alleles.hasOwnProperty(key)) {
            if (this.alleles[key] == alleleVal) {
                return key;
            }
        }

    }
    return -1;

};

Gene.prototype.toSerial = function() {
    var out = '';
    out+= 'Gene: ' + this.num;

    return out;

};
Gene.prototype.printGene = function() {
    var ret = 'gene: ';
    ret += this.num;


    for (var g in this.alleles) {


        if (this.alleles.hasOwnProperty(g)) {
            ret += ' ';
            ret += this.alleles[g];

        }

    }
    ret += ' ';
    ret +=this.info;

    return ret;
};


function StudentGene(num)  {

    StudentGene.maxPrefs = 5;
    StudentGene.maxExcls = 3;

    Gene.call(this,num);
    this.name = 'None';
    this.prefs = ['None'];
    this.excls = ['None'];

}

StudentGene.prototype = new Gene(null);

StudentGene.prototype.loadStudentInfo = function(info)   {
    //this.studentInfo = info;
    var  spl = info.split(',');
    this.name = spl[0];//[1];
    this.prefs = [spl[1]  || '',spl[2] || '',spl[3] || '',spl[4] || '',spl[5] || ''];
    this.excls = [spl[6] || '',spl[7] || '',spl[8] ? spl[8].trim() : ''];


};

StudentGene.prototype.toSerial = function() {
    var out = '';
    out+= this.name;
    out+=',';
    this.prefs.forEach(function(el) {
        out+=el + ',';

    });
    this.excls.forEach(function(el) {
        out+=el + ',';
    });

    return out;

};

StudentGene.prototype.printGene = function() {
    var notUsed = Gene.prototype.printGene.call(this);

    var ret = ' Name: ' + this.name;
    ret += ' prefs: ' + this.prefs;
    ret+=  ' Excls: ' +  this.excls;
    return ret;
};

StudentGene.prototype.markText = function(color,text)    {
    return '<mark class = ' + '"' + color + '"> ' +  text + '</mark>';

};


StudentGene.prototype.prefWithIndex = function(ind) {
    return this.prefs[ind];
};

StudentGene.prototype.exclWithIndex = function(ind) {
    return this.excls[ind];
};

StudentGene.prototype.setPrefWithIndex = function(ind,val) {
    if (val == 0) {
        this.prefs[ind] = '';
    }
    else {
        this.prefs[ind] = val;
    }

};

StudentGene.prototype.setExclWithIndex = function(ind,val) {
    if (val == 0) {
        this.excls[ind] = '';
    }
    else {
        this.excls[ind] = val;
    }


};

StudentGene.prototype.printGeneWithHits =function(prefHits,exclHits)    {
    var unused = Gene.prototype.printGene.call(this);

    var genePrintAr = [];

    var ret = ' Name: ' + this.name;
    genePrintAr.push(this.name);

    ret += ' Prefs: ';
    var i;
    for (i = 0; i < this.prefs.length; ++i)  {
        var curPref = '';

        if (prefHits && (prefHits[i] == 1))    {
            curPref=this.markText('green',this.prefs[i]);
        }
        else if (prefHits && (prefHits[i] == -1)) {
            curPref=this.markText('error',this.prefs[i]);
        }
        else {

            curPref= this.prefs[i];
        }
        ret+=curPref + ',';
        genePrintAr.push(curPref);

    }

    //ret += '<mark class = "green"> '  + ' prefs: ' + this.prefs;
    ret+=  ' Excls: ' ;
    for (i = 0; i < this.excls.length; ++i)  {
        var curExcl = '';
        if (exclHits && (exclHits[i]) == 1)    {
            curExcl=this.markText('red',this.excls[i]);
        }
        else if (exclHits && (exclHits[i] == -1)) {
            curExcl = this.markText('error', this.excls[i]);

        }
        else if (exclHits && (exclHits[i] == 0)) {
            curExcl=this.markText('green',this.excls[i]);
        }
        else    {
            //curExcl=this.markText('green',this.excls[i]);
            curExcl= this.excls[i];
        }
        ret+= curExcl + ',';
        genePrintAr.push(curExcl);

    }
    //console.log(ret);
    return [ret,genePrintAr];

};


function GenePool(alleleValueAlloc,numGenesToCreate,allelesPerGene,initRange)   {
    this.genes = [];
    this.allelesPerGene = allelesPerGene;
    this.initRange = initRange;
    this.alleleValueAllocRoutine = alleleValueAlloc;



    var g;
    for (var i = 0;i < numGenesToCreate; ++i) {
        //g = new Gene(i+1);
        g = new StudentGene(i + 1, 'No stud info');
        /*
         var alList = alleleValueAlloc(allelesPerGene,initRange);
         for (var j = 0;j < alList.length; ++j ) {
         g.addAllele(alList[j]);
         }
         */
        //this.addGene(g);
        if ((this.allelesPerGene > 0) && this.alleleValueAllocRoutine) {

            this.addAllelesToGenes(this.allelesPerGene);
        }

        this.genes.push(g);
    }

    GenePool.instance = this;

/*
    this.addInfo = function(geneNum,info)   {
        this.gene[geneNum].info = info;
    };
*/

/*
    this.printGenes = function(n)   {
        var ret = '';
        for (var i = 0; i < n; ++i)   {
            //ret += this.genes[i].printGeneWithHits([1,1,0,1,0],[1,0,0])//printGene();
            ret += this.genes[i].printGene();
            ret += '<br>';
        }
        return ret;

    };
    */

}

GenePool.getInstance = function() {

    if (typeof GenePool.instance === 'object') {
        return GenePool.instance;
    }
    return new GenePool();


};
GenePool.prototype.printGenes = function(n)  {
    var ret = '';
    for (var i = 0; i < n; ++i)   {
        //ret += this.genes[i].printGeneWithHits([1,1,0,1,0],[1,0,0])//printGene();
        ret += this.genes[i].printGene();
        ret += '<br>';
    }
    return ret;

};

GenePool.prototype.clearAllAllelesFromGenes = function() {
    this.allelesPerGene = 0;
    this.genes.forEach(function(g) {
        g.alleles = [];
    });


};
GenePool.prototype.addAllelesToGenes = function(numAls) {

    this.allelesPerGene = numAls;
    var alList;
    var that = this;
    this.genes.forEach(function(g) {
        alList = that.alleleValueAllocRoutine(that.allelesPerGene,that.initRange);
        for (var j = 0;j < alList.length; ++j ) {
            g.addAllele(alList[j]);
        }
    });

};

GenePool.prototype.addGene = function(g)  {
    this.genes.push(g);
};

GenePool.prototype.numGenes = function() {
    return this.genes.length;

};

GenePool.prototype.maxNumAlleles = function() {
    var b = 0;
    for (var i = 0;i < this.numGenes(); ++i)    {

        if (this.genes[i].numAlleles() > b) {
            b = this.genes[i].numAlleles();
        }

    }
    return b;
};

GenePool.prototype.geneWithNum = function(num)    {
    for (var i = 0; i < this.numGenes(); ++i) {
        if (this.genes[i].num == num) {
            return this.genes[i];
        }
    }
    return null;

};

GenePool.prototype.getNum = function(ind) {
    return this.genes[ind].num;
};

GenePool.prototype.geneWithIndex = function(ind)  {
    return this.genes[ind];
};

GenePool.prototype.getIndex = function(num)   {
    //num = num.trim();
    for (var i = 0;i < this.numGenes(); ++i) {
        if (this.genes[i].num == num)   {
            return i;
        }

    }
    //console.log('not found: ' + num);
    return -1;
};
/*
 this.addInfo = function(geneNum,info)   {
 this.gene[geneNum].info = info;
 };
 */

GenePool.prototype.toSerial = function() {
    out = '';
    this.genes.forEach(function(g) {
        out+= g.toSerial() + '\r\n';
    });
    return out;

};

GenePool.prototype.printGenePool = function(numToPrint) {
    var n = numToPrint || this.numGenes();
    if (n > this.numGenes())    {
        n = this.numGenes();
    }
    var ret = 'Num Genes: ';
    ret+= this.numGenes();
    ret += '<br> Max num alleles per gene:  ';
    ret+= this.maxNumAlleles();
    ret += '<br>';

    ret += this.printGenes(n);
    return ret;

};
/*
 this.printGenes = function(n)   {
 var ret = '';
 for (var i = 0; i < n; ++i)   {
 //ret += this.genes[i].printGeneWithHits([1,1,0,1,0],[1,0,0])//printGene();
 ret += this.genes[i].printGene();
 ret += '<br>';
 }
 return ret;

 };
 */

GenePool.prototype.alleleVal = function(geneNum,alleleId) {
    //var tmp = this.genes[geneNum].alleles[alleleId];
    return this.genes[geneNum].alleles[alleleId];
};


function StudentGenePool(alleleValueAlloc,allelesPerGene,studentArray) {



    GenePool.call(this, alleleValueAlloc, 0, allelesPerGene, 50);

    this.printGenes = function(n)   {
        var ret = '';
        for (var i = 0; i < n; ++i)   {
            //ret += this.genes[i].printGeneWithHits([1,1,0,1,0],[1,0,0])//printGene();
            ret += this.genes[i].printGene();
            ret += '<br>';
        }
        return ret;
    };


    this.printStudentGenes = function(n,hitArray)   {
        var ret = '';
        var studentGenesAr = [];
        for (var i = 0; i < n; ++i)   {
           // ret += this.genes[i].printGeneWithHits([1,1,0,1,0],[1,0,0])//printGene();
            if (hitArray) {
                ret += this.genes[i].printGeneWithHits(hitArray[i][0], hitArray[i][1])[0];
                studentGenesAr.push(this.genes[i].printGeneWithHits(hitArray[i][0], hitArray[i][1])[1]);
            }
            else {
                ret += this.genes[i].printGeneWithHits(null,null)[0];
                studentGenesAr.push(this.genes[i].printGeneWithHits(null,null)[1]);

            }
            //ret += this.genes[i].printGene();
            ret += '<br>';
        }
        return [ret,studentGenesAr];
    };

    this.printStudentGenePool = function(numToPrint,bestIndiv) {
       // console.log('xxx start print gp');
        var n = numToPrint || this.numGenes();

        var hitArray = null;
        var bestId = null;
        var prefTotsArray = null;
        var prefNumsArray = null;
        var c = null;

        //console.log('xxx   start bestindiv stuff');
        if (bestIndiv) {
            hitArray = bestIndiv.checkHits();
            bestId = bestIndiv.num;
            prefTotsArray = bestIndiv.checkPrefsTotals();
            prefNumsArray = bestIndiv.checkNumPrefs();
        }
       // console.log('xxx   end bestindiv stuff');

        if (n > this.numGenes())    {
            n = this.numGenes();
        }
        var ret = 'Num Genes: ';
        ret+= this.numGenes();
        ret += ' Max num alleles per gene:  ';
        ret+= this.maxNumAlleles();
        ret += '<br>';
        ret+= 'BestIndiv: ' + (bestId ? bestId : '-') + ' Cost: ' + (c ? c : '-');
        ret+='<br>';
        ret+= 'Prefs Met: ' + (prefTotsArray ? prefTotsArray[0] : '-') + ' Excls not met: ' + (prefTotsArray ? prefTotsArray[1] : '-');
        ret+='<br>';
        ret+='At Least Prefs Met: ' + (prefNumsArray ? prefNumsArray[0] : '-') + ' At Least Excls not met: ' + (prefNumsArray ? prefNumsArray[1] : '-');
        ret+='<br>';

        //console.log('xxx   end ret stuff');
        var r =  this.printStudentGenes(n,hitArray);
        //ret += r[0];  this prints out simple  text version. Now returning array instead to be formatted as html
        var studGenesAr = r[1];

       // console.log('xxx end print gp');
        return [ret,studGenesAr];

    };

    this.loadStudentsAsGenes = function(alleleValueAlloc,allelesPerGene,studentArray)   {
        var ret = '';
        for (var i = 0;i < studentArray.length; ++i)    {
            var studentStr = studentArray[i];
            var spl = studentStr.split(',');
            //console.log('spl: ' + spl);
            //alert(spl);
            // alert(spl[0]);
            if (spl[0].length == 0) {
                continue;
            }
            var studG= new StudentGene(i+1);  //spl[0]);//gPool.geneWithNum(spl[0]);
            if (studG) studG.loadStudentInfo(studentArray[i]);
            var alList = alleleValueAlloc(this.allelesPerGene,this.initRange);
            for (var j = 0;j < alList.length; ++j ) {
                studG.addAllele(alList[j]);
            }
            //this.addGene(g);
            this.genes.push(studG);

            //if (g)  g.info = studentStr;
            // gPool.addInfo(spl[0],studentStr);
            //if (g) alert(g.info);
            ret+=studentArray[i];
            ret+='<BR>';
        }
    };

    if (studentArray) {
        this.loadStudentsAsGenes(alleleValueAlloc, allelesPerGene, studentArray);
    }
}

StudentGenePool.prototype = new GenePool(null,null,null,null);

function Indiv(num,numInPop,parent1,parent2) {

    this.genome = [];
    this.ancestors = [];
    this.rankScore = 0;
    this.rankCumScore = 0;
    this.numInPop = numInPop;
    this.genNum = null; //only used for best population (ie which has mixed generations - so diff indivs in a pop may have diff generations

   // this.classTolerance = 1;

    //this.singletonParams = Params.getInstance();

    // now uses myParam.  this.crossOverProb = 1.0;//1.0;//1.0; //1.0;// 1.0; //0.5

    if (num) {
        this.num = num;
    } else {
        this.num =  ++Indiv.lastIndivNum;//++lastIndivNum;
    }

    this.parent1 = parent1 || null;
    this.parent2 = parent2 || null;

    if ((this.parent1) && (this.parent2)) {
        this.genome = this.parent1.genome.slice(0);
        this.crossOver();

    }
    else if (this.parent1) {
        this.genome = this.parent1.genome.slice(0);
    }
    else if (this.parent2)   {
        this.genome = this.parent2.genome.slice(0);
    }
    else {
        if (num) { // do nothing - this is a copy ie num provided but not parents provided
        }
        else {
            this.genesis();
        }
    }
    if (this.parent1) {

        this.ancestors = this.parent1.ancestors.slice(0); //patriarchal line
        var num1,num2;
        if (this.parent1)   {
            num1 = this.parent1.num;
        }
        else    {
            num1 = null;
        }
        if (this.parent2)   {
            num2 = this.parent2.num;
        }
        else    {
            num2 = null;
        }
        this.ancestors.push([num1,num2]);

    }





}
Indiv.lastIndivNum = 0;

Indiv.prototype.genesis = function () {

    var done = false;

    while (!done) {

        this.genome = [];
        gPool.genes.forEach(function (g) {
            //console.log('for eaching: ' + g);
            if (g) {

            }
            else alert('no g');
            var n = g.numAlleles();
            var r = Math.floor(Math.random() * n);
            var keys = getKeys(g.alleles);
            this.genome.push(keys[r]);

        }, this);


        var c = this.getAlleleValCounts()[1];

        var good = true;

        var tol = Params.getInstance().classTolerance;
        if (tol == 0) {
            if  (gPool.numGenes() % gPool.allelesPerGene != 0) {  //impossible to have even classes
                tol = 1;
            }
        }

        for (var i = 0; i < c.length; ++i) {
            if ((c[i] < ((gPool.numGenes() * 1.0 / gPool.allelesPerGene) -  tol)) // Params.getInstance().classTolerance)) //this.classTolerance))
                || (c[i] > ((gPool.numGenes() * 1.0 / gPool.allelesPerGene) +  tol))) { //Params.getInstance().classTolerance))) { //this.classTolerance))) {
                good = false;
                break;
                // totFitness -=50;
            }

        }



       // if ((good) || ((Params.getInstance().classTolerance == 0) && (gPool.numGenes() % gPool.allelesPerGene != 0))) { //if class tolerance is zero, don't try for even classes if impossible
        if (good) {
            done = true;
        }
    }



    /*
     for (var i = 0;i < gPool.numGenes();++i) {

     var g = gPool.genes[i];
     if (g) {
     } else alert('no g');
     var n = g.numAlleles();
     var r = Math.floor(Math.random() * n);
     var keys = getKeys(gPool.genes[i].alleles);
     this.genome.push(keys[r]);
     }
     */

};

Indiv.prototype.toSerial = function() {
    ret = '';
    if (gPool) {
        this.genome.forEach(function(allele,i) {
            ret+=gPool.alleleVal(i,allele);
        });
    }
    else {
        ret +='No gPool';
    }
    return ret;

};


Indiv.prototype.printIndiv = function(hamDist)   {
    hamDist = hamDist || null;


    var ret = '<p id = "indiv' + this.numInPop + '">';
    if (this.genNum != null) {
        ret+= ' Gen: ' + String("0000" + this.genNum).slice(-5);
    }
    ret+=' ' + String("0000" + (this.numInPop + 1)).slice(-3);
    ret += ' ';

    this.genome.forEach(function(allele,i) {
        ret+=gPool.alleleVal(i,allele);
    });
    /*
     for (var i=0;i <this.genome.length;++i) {
     //ret +=gPool.genes[i].alleleVal(this.genome[i]);
     ret += gPool.alleleVal(i,this.genome[i]);


     }
     */
    ret += '       ' + this.getAlleleValCounts()[0];

    ret += ' Fit: ' + this.cost();

    if (hamDist) {
        ret += ' HD: ' + hamDist.toFixed(2);
    }

    ret += '</p>';

    return ret;

};

Indiv.prototype.toTable = function() {
    var allocTable = [];

    var classIndexes = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'J':9,'K':10,'L':11,'M':12,'N':13,'O':14,'P':15};
    var classArray = [];
    var aC = this.getAlleleValCounts();
    var acArray = aC[1];
    var numAlleleVals = acArray.length;
    var biggestClass = 0;
    var i;
    for (i = 0;i < numAlleleVals;++i) {
        if (acArray[i] > biggestClass) {
            biggestClass = acArray[i];
        }
    }

    for (i = 0;i < gPool.allelesPerGene; ++i) {
        allocTable.push(0);
    }


    var row;
    var j;
    for (i = 0;i < biggestClass; ++i) {
        row = [];  //new Array(gPool.allelesPerGene);//(numAlleleVals);
        for (j = 0; j < gPool.allelesPerGene; ++j) { //numAlleleVals; ++j) {
            row.push('Bob');
        }
        classArray.push(row);

    }

    var alVal;
    this.genome.forEach(function(allele,i) {
        alVal = gPool.alleleVal(i,allele);
        var alInd = classIndexes[alVal];
        var alAllocInd = allocTable[alInd];
        classArray[alAllocInd][alInd] = i + '';//  gPool.geneWithIndex(i).name; // i;
        allocTable[alInd]+=1;
    });

    return classArray;


};

Indiv.prototype.getAlleleCounts = function()  {
    var aCount = {};
    //var i;
    this.genome.forEach(function(allele) {
        aCount[allele] +=1;
    });
    /*
     for (i = 0;i < this.genome.length; ++i) {
     aCount[this.genome[i]] +=1;
     }
     */
    var ret = '';
    for (var c in aCount) {
        if (aCount.hasOwnProperty(c)) {
            ret += 'allele: ' + c + ' count: ' + aCount[c];
        }

    }

    return ret;

};

Indiv.prototype.getAlleleValCounts = function()    {

    var aCount = {};
    //var i;
    var val;

    this.genome.forEach(function(allele,i) {
        val = gPool.alleleVal(i, allele);
        if (val in aCount) {
            aCount[val] += 1;
        }
        else {
            aCount[val] = 1;
        }
    } );
    /*
     for (i = 0;i < this.genome.length; ++i) {
     val = gPool.alleleVal(i,this.genome[i]);
     if (val in aCount)  {
     aCount[val] +=1;
     }
     else    {
     aCount[val] = 1;
     }

     }
     */
    var ret = '';
    var keys = getKeys(aCount);

    var aCountArray = [];

    keys.forEach(function(k)    {
        ret += ' ' + k + ': ' + aCount[k];
        aCountArray.push(aCount[k]);

    });
    /*
     for (i = 0;i < keys.length; ++i )   {
     ret += ' ' + keys[i] + ': ' + aCount[keys[i]];
     aCountArray.push(aCount[keys[i]]);
     }
     */


    return [ret,aCountArray];

};

Indiv.prototype.redVector = function(a,b) {

    var ha = a[0];
    var hb = b[0];
    var ea = a[1];
    var eb = b[1];


    ha = ha.map(function(el) {
        if (el == -1) {
            return 0;
        }
        else {
            return el;
        }

    });
    hb = hb.map(function(el) {
        if (el == -1) {
            return 0;
        }
        else {
            return el;
        }

    });
    ea = ea.map(function(el) {
        if (el == -1) {
            return 0;
        }
        else {
            return el;
        }

    });
    eb = eb.map(function(el) {
        if (el == -1) {
            return 0;
        }
        else {
            return el;
        }

    });

    var out = [];
    if (ha) {
        out = ha.map(function (item, ind) {
            return item + hb[ind];

        });
    }
    else {

        out = hb;
    }
    var out2 = [];
    if (ea) {
        out2 = ea.map(function (item, ind) {
            return item + eb[ind];

        });
    }
    else {

        out2 = eb;
    }


    return [out,out2];
};

Indiv.prototype.checkPrefsTotals = function() {
   // console.log('xxx  Indiv to disp start check pref totals');
    var hitsArray = this.checkHits();
    var prefsTots =hitsArray.reduce(this.redVector);
    //var exclsTots =hitsArray.reduce(this.redVector);

    //console.log('xxx  Indiv to disp end check pref totals');
    return prefsTots;


};

Indiv.prototype.checkHits = function() {
    var hitsArray = [];
    var h;

    this.genome.forEach(function(allele,i) {
        h = this.alleleHits(i);
        hitsArray.push(h);
    },this);
    /*
     for (var i = 0; i < this.genome.length; ++i) {

     h = this.alleleHits(i);
     hitsArray.push(h);

     }
     */
    return hitsArray;
};

Indiv.prototype.checkNumPrefs = function() {

    var hitsArray = this.checkHits();
    var numPrefsArray = hitsArray.map(function (el) {
        var numPrefs = el[0].reduce(function (e1, e2) {
            if (e1 == -1) e1 = 0;
            if (e2 == -1) e2 = 0;
            return e1 + e2;
        });
        var numExcls = el[1].reduce(function (e1, e2) {
            if (e1 == -1) e1 = 0;
            if (e2 == -1) e2 = 0;
            return e1 + e2;
        });
        return [numPrefs, numExcls];


    });


    var numPrefsHash = {};
    var numExclsHash = {};


    var i;
    for (i = 0; i < StudentGene.maxPrefs + 1; ++i) {

        numPrefsHash[i] = 0;
    }
    for (i = 0; i < StudentGene.maxExcls + 1; ++i) {

        numExclsHash[i] = 0;
    }




    numPrefsArray.forEach(function(numPref) {
        if (numPrefsHash.hasOwnProperty(numPref[0])) {
            numPrefsHash[numPref[0]] +=1;
        }
        else {
            numPrefsHash[numPref[0]] = 1;
        }
        if (numExclsHash.hasOwnProperty(numPref[1])) {
            numExclsHash[numPref[1]] +=1;
        }
        else {
            numExclsHash[numPref[1]] = 1;
        }

    });

    var prefsHashKeys = getKeys(numPrefsHash);
    var exclsHashKeys = getKeys(numExclsHash);
    //var prefsCum = 0;
    //var exclsCum = 0;

    prefsHashKeys.sort(function(a, b) {
        return a - b;
    });
    exclsHashKeys.sort(function(a, b) {
        return a - b;
    });

    var atLeastPrefs = prefsHashKeys.map(function(el,ind) {
        var tot = 0;
        if (ind == 0) {
            tot = numPrefsHash[0];
        }
        else {
            for (var i = ind; i < prefsHashKeys.length; ++i) {
                tot += numPrefsHash[i];


            }

        }
        return tot;
    });
    var atLeastExcls = exclsHashKeys.map(function(el,ind) {
        var tot = 0;
        if (ind == 0) {
            tot = numExclsHash[0];
        }
        else {
            for (var i = ind; i < exclsHashKeys.length; ++i) {
                tot += numExclsHash[i];


            }

        }
        return tot;
    });






    return [atLeastPrefs,atLeastExcls];


};

Indiv.prototype.alleleHits = function(ind)    {

    var g = gPool.geneWithIndex(ind);

    var alPrefHits; // = new Array(g.prefs.length);
    var i;

    alPrefHits = g.prefs.map(function() {
        return 0;
    });
    /*
     for (i = 0;i < g.prefs.length; ++i) {
     alPrefHits[i] = 0;
     }
     */


    var alExclHits = g.excls.map(function() {
        return 0;
    });
    /*
     var alExclHits = new Array(g.excls.length);
     for (i = 0;i < g.excls.length; ++i) {
     alExclHits[i] = 0;
     }
     */


    //console.log('ind: ' + ind);
    var indOther;
    for (i = 0; i < g.prefs.length;++i) {
        indOther = gPool.getIndex(g.prefs[i]);
        if (indOther == -1) {
            //console.log('not found in allele hits: ' + g.prefs[i]);
            alPrefHits[i] = -1;
            continue;
        }
        //console.log('  ind other: ' + indOther);
        var otherG;
        otherG = gPool.geneWithIndex(indOther);
        //console.log('genome: ' + this.genome[ind] + ' other genome: ' + this.genome[indOther]);
        var val = gPool.alleleVal(ind,this.genome[ind]);
        var valOther = gPool.alleleVal(indOther,this.genome[indOther]);
        //console.log('genome val: ' + val + ' other genome val: ' + valOther);
        if (val == valOther)  {
            alPrefHits[i] = 1;

        }

    }
    //console.log('cost of: ' + g.name + ' ' + alCost);


    for (i = 0; i < g.excls.length;++i) {
        //console.log('excl: ' + i);
        indOther = gPool.getIndex(g.excls[i]);
        if (indOther == -1) {
            //console.log('not found in allele excl hits: ' + g.excls[i]);
            alExclHits[i] = -1;
            continue;
        }
        //console.log('excl index other: ' + indOther);
        otherG = gPool.geneWithIndex(indOther);
        val = gPool.alleleVal(ind,this.genome[ind]);
        valOther = gPool.alleleVal(indOther,this.genome[indOther]);
        //console.log('genome excl val: ' + val + ' other genome excl val: ' + valOther);
        if (val == valOther)  {
            alExclHits[i] = 1;

        }

    }
    //console.log('cost of excl: ' + g.name + ' ' + alExclCost);
    //console.log ('   alcost: ' + (alCost - alExclCost));
    return [alPrefHits,alExclHits];


};



Indiv.prototype.cost = function()  {

    var totFitness = 0;
    //console.log('working out cost of: ' + this.genome);
    var i;
    for (i = 0;i < this.genome.length;++i)   {

        //console.log('val of: ' + this.genome[i] + ' val: ' + gPool.alleleVal(i,this.genome[i]));

    }


    var ser = this.toSerial();
    for (i = 0;i < this.genome.length;++i)   {


       // totFitness += this.alleleCost(i);
        totFitness += this.alleleQuickCost(i,ser); //much faster!

    }

    //console.log('Tot cost for indiv: ' + i + ' ' + totFitness);

    var c = this.getAlleleValCounts()[1];

    var tol = Params.getInstance().classTolerance;
    if (tol == 0) {
        if  (gPool.numGenes() % gPool.allelesPerGene != 0) {  //impossible to have even classes
            tol = 1;
        }
    }

    for (i = 0;i < c.length; ++i) {
        if ((c[i] < ((gPool.numGenes() * 1.0 / gPool.allelesPerGene) -  tol)) // Params.getInstance().classTolerance)) // this.classTolerance))
            ||     (c[i] > ((gPool.numGenes() * 1.0 / gPool.allelesPerGene) +  tol ))) //Params.getInstance().classTolerance))) // this.classTolerance)))
        {

            // totFitness -=50;
        }
        else {
            totFitness +=50;
        }
    }

    if (totFitness <= 0)    {
        totFitness = 1; //arbitrarily assign 1 as lowest possible fitness
    }
    return totFitness;
};

Indiv.prototype.alleleQuickCost = function(ind,ser)   {
    var alCost = 0;

    var i;
  //  var indOther;
    var g = gPool.geneWithIndex(ind);

    for (i = 0; i < g.prefs.length;++i) {


        if (g.prefs[i].length == 0) {
            //console.log('not found: ' + g.prefs[i]);
            continue;
        }
        //console.log('  ind other: ' + indOther);

        if (ser[ind] == ser[g.prefs[i] - 1])  {
            alCost+=(g.prefs.length - i); // was 1 each

        }

    }
    if (alCost == 0) {
        //alCost-=10;
    } //penalise if no prefs at all
    else {
        alCost+=12; //10
    }
    //console.log('cost of: ' + g.name + ' ' + alCost);

    var alExclCost= 0;
    for (i = 0; i < g.excls.length;++i) {

        if (g.excls[i].length == 0) {
            //console.log('not found: ' + g.excls[i]);
            continue;
        }
        //console.log('excl index other: ' + indOther);

        //console.log('genome excl val: ' + val + ' other genome excl val: ' + valOther);
        if (ser[ind] == ser[g.excls[i] - 1]) {
        }
        else {
            alExclCost+=10;

        }

    }
    //console.log('cost of excl: ' + g.name + ' ' + alExclCost);
    //console.log ('   alcost: ' + (alCost - alExclCost));
    return  alCost + alExclCost; //alCost - alExclCost;


};

Indiv.prototype.alleleCost = function(ind)   {
    var alCost = 0;
    var g = gPool.geneWithIndex(ind);
    //console.log('ind: ' + ind);
    var i;
    var indOther;
    for (i = 0; i < g.prefs.length;++i) {
        indOther = gPool.getIndex(g.prefs[i]);
        if (indOther == -1) {
            //console.log('not found: ' + g.prefs[i]);
            continue;
        }
        //console.log('  ind other: ' + indOther);
        var otherG = gPool.geneWithIndex(indOther);
        //console.log('genome: ' + this.genome[ind] + ' other genome: ' + this.genome[indOther]);
        var val = gPool.alleleVal(ind,this.genome[ind]);
        var valOther = gPool.alleleVal(indOther,this.genome[indOther]);
        //console.log('genome val: ' + val + ' other genome val: ' + valOther);
        if (val == valOther)  {
            alCost+=(g.prefs.length - i); // was 1 each

        }

    }
    if (alCost == 0) {
        //alCost-=10;
    } //penalise if no prefs at all
    else {
        alCost+=12; //10
    }
    //console.log('cost of: ' + g.name + ' ' + alCost);

    var alExclCost= 0;
    for (i = 0; i < g.excls.length;++i) {
        //console.log('excl: ' + i);
        indOther = gPool.getIndex(g.excls[i]);
        if (indOther == -1) {
            //console.log('not found: ' + g.excls[i]);
            continue;
        }
        //console.log('excl index other: ' + indOther);
        otherG = gPool.geneWithIndex(indOther);
        val = gPool.alleleVal(ind,this.genome[ind]);
        valOther = gPool.alleleVal(indOther,this.genome[indOther]);
        //console.log('genome excl val: ' + val + ' other genome excl val: ' + valOther);
        if (val == valOther) {
        }
        else {
            alExclCost+=10;

        }

    }
    //console.log('cost of excl: ' + g.name + ' ' + alExclCost);
    //console.log ('   alcost: ' + (alCost - alExclCost));
    return  alCost + alExclCost; //alCost - alExclCost;


};



Indiv.prototype.hamDist = function(other) {
    var dist = 0;
    this.genome.forEach(function(el,i) {
        if (el === other.genome[i]) {

        }
        else {
            ++dist;
        }

    });
    return dist;

};

Indiv.prototype.sameAs = function(compareIndiv) {

    var doneChar = '-';
    //var notSame = false;

    var aArray = this.genome.map(function (el,i) {
        //var gNum  = gPool.getNum(i);
        return gPool.alleleVal(i,el);

    });

    var bArray =  compareIndiv.genome.map(function (el,i) {
        //var gNum  = gPool.getNum(i);
        return gPool.alleleVal(i,el);

    });

    for (var i = 0; i < gPool.maxNumAlleles() - 1; ++i) {
        var startChar = '';
        var translateChar = '';
        var j;
        for (j = 0; j < this.genome.length; ++j) {

            if (aArray[j] == doneChar) {

            }
            else {
                startChar = aArray[j];
                translateChar = bArray[j];
                break;

            }

        }
        for (j = 0; j < this.genome.length; ++j) {
            if (aArray[j] == startChar) {
                if (bArray[j] == translateChar) {
                    aArray[j] = doneChar;
                }
                else {
                    return false;
                }
            }
        }

    }

    return true;

};

Indiv.prototype.copyIndiv = function()    {

    var newIndiv = new Indiv(this.num,this.numInPop);
    newIndiv.genome = this.genome.slice(0);
    newIndiv.ancestors = this.genome.slice(0);
    newIndiv.crossOverProb = this.crossOverProb;
    newIndiv.genNum = this.genNum;
    //newIndiv.classTolerance = this.classTolerance;
    return newIndiv;


};

Indiv.prototype.crossOver = function() {
    var r = Math.random();
    var par = Params.getInstance();
    if (r < par.crossoverRate) {  // this.singletonParams.crossoverRate) { //this.crossOverProb) {
        var crossPoint = Math.floor(Math.random() * (this.genome.length - 1)) + 1; //was 15 15); //was 30
        this.genome = this.parent1.genome.slice(0,crossPoint).concat(this.parent2.genome.slice(crossPoint));
    }

};


/*
var dummyFuncHandle;
//function dummyFunc(i,tot)    {
function dummyFunc()    {
    var i = 20;
    var tot = 100;
    //setTimeout(dummyFunc,50);
    document.getElementById("progress").innerHTML +=  'Loading..' + (i * 100.0 / tot) + '%';
    console.log('in dummy func');
   // dummyFunc();


}
*/
//var pos = 0;
//var globInterval;





function Population(callBackFun,genNum,initSize,prevPop) {

    this.firstName = 'Blah';
    this.lastName = 'Bloo';
    this.mutRate = 0.25;//0.25

    //this.singletonParams = Params.getInstance();
    var par = Params.getInstance();

    if (prevPop) {
        this.prevPop = prevPop;
        this.genNum = this.prevPop.genNum + 1;
        this.mutRate = this.prevPop.mutRate;
        //this.elite = this.prevPop.elite;
        this.elite = par.eliteFlag; //this.singletonParams.eliteFlag;
        this.round = this.prevPop.round;
    }
    else    {
        this.prevPop = null;
        this.genNum = genNum;
        this.mutRate = par.mutRate; //this.singletonParams.mutRate;
        this.elite = par.eliteFlag; //this.singletonParams.eliteFlag;
        this.round = 0;
    }

    this.indivs = [];
    this.initSize = initSize;

   // this.genNum = genNum;

    //this.elite =  false; //true; //true; //true; //true; //true; //true; //true;

    this.callBackFun = callBackFun;

    var catastropheThreshold = 500; //500
   //if (this.genNum % 2000 == 0) {  //(this.genNum == 0) {
    var par = Params.getInstance();

    if (this.callBackFun) {
        if ((this.genNum == 0) || (par.numGensNoProgress >= par.catastropheGens )) { // (this.singletonParams.numGensNoProgress >= catastropheThreshold )) {
            ++this.round;
            this.genesis();
        }
        else {
            this.nextGeneration();
        }
    }

/*
    this.callBackFun = function() {
        console.log(' in inline cb');
        console.log('DONE');
        console.log('in cb pop size: = ' + myPop.popSize() + ' indivs len: ' + myPop.indivs.length);
        var fred = gPool.printGenePool();
        var bob = myPop.printPopulation();
        //alert(myPop.indivs[0].getAlleleValCounts());
        document.getElementById("gp").innerHTML =  fred;
        var popEl = document.getElementById("pop");
        popEl.innerHTML =  bob;
        //document.getElementById()
        popEl.scrollTop = 10;


        var best = myPop.getBestCost();
        var pickIndiv = myPop.indivs[best[1]];
        bestCostArray.push(best[0]);
        //var chk = pickStud.checkHits();
        fred = gPool.printStudentGenePool(null,pickIndiv);
        document.getElementById("gp").innerHTML =  fred;
        //console.log('check hits: ' + chk[0]);
        // console.log('check excls: ' + chk[1]);
        //console.log('best: ' + best);

        document.getElementById("popToDisp").value = best[1]+1;

        var all = myPop.getAllCosts();

        document.getElementById("Gen").innerHTML = myPop.genNum;
        document.getElementById("Best").innerHTML = myPop.getBestCost()[0];

        graphBestCosts();


        overallGen +=1;
        //pos = 0;

    };
    */


    /*
   //one technique (not sure though - because doesn't wait to do next one)
    this.doABitMoreTest = function()    {
        console.log('in do a bit more test');
        var interval = this.initSize / 10;
        var tot = this.initSize;
        var targetThisCall = pos + interval;
        console.log('targ: ' + targetThisCall);
        for (var i = pos;i < targetThisCall;++i) {
            var indii = new Indiv(0); //(i + 1);

            this.indivs.push(indii);


        }
        pos += interval;
        document.getElementById("progress").innerHTML =  'Loading..' + (pos * 100.0 / tot) + '%';
        if (pos >= this.initSize)   {
            clearInterval(globInterval);
            this.callBackFun();
        }

    };
    */

/* one technique - but may not finish one iter before next one
            globInterval = setInterval(function() {
                that.doABitMoreTest()},100);
*/

        /*
        for (var i = 0; i < this.initSize; ++i) {
            if (i % interval == 0)  {
                //window.setTimeout(dummyFunc,50);
                //dummyFunc();
                //window.setTimeout(dummyFunc(i),50);
                  //window.requestAnimationFrame(function()   {
                  //    dummyFunc(i,tot);
                  //})
                  //dummyFunc);
            }





            var indii = new Indiv(0); //(i + 1);

            this.indivs.push(indii);
        }
        */

}

//Test only:
Object.defineProperty(Population.prototype, 'fullName', {
    get: function() {
        return this.firstName + ' ' + this.lastName;
    },
    set: function(name) {
        var words = name.split(' ');
        this.firstName = words[0] || '';
        this.lastName = words[1] || '';
        //alert('first: ' + this.firstName + ' last: ' + this.lastName + ' full: ' + this.fullName);
    }
});

Population.CreateFromJSON = function(j) { //static method!
    var newPop = new Population();
    for (var prop in j) {
        if (j.hasOwnProperty(prop))  {
            newPop[prop] = j[prop];

        }
    }
    return newPop;

}

Population.prototype.genesis = function() {

    //var intervalHandle = setInterval(function() {
    //    dummyFunc()
    //},100);

    // var interval = this.initSize / 10;
    //  var tot = this.initSize;


    this.iterationF = function(currPos,interval,tot)    {
        //console.log('in iterF pos: ' + currPos + ' interval: ' + interval +  ' pos+inter: ' + currPos + interval);
        //console.log(' indivs size: ' + this.indivs.length);
        var that = this;
        var targetThisCall = currPos + interval;
        //console.log('targ: ' + targetThisCall);
        for (var i = currPos;i < targetThisCall;++i) {
            var indii = new Indiv(0,i); //(i + 1);

            that.indivs.push(indii);


        }
        currPos += interval;

        //should move this out of here. Shouldn't referenc dom objects, as Population is a model
        if (typeof document !== 'undefined') {
            document.getElementById("progress").innerHTML = 'Loading..' + (currPos * 100.0 / tot) + '%';
        }
        if (currPos < tot)  {
            setTimeout(function() {
                that.iterationF(currPos,interval,tot)
            },50);

        }
        else    {
            var par = Params.getInstance();
          //  if (this.singletonParams.trialClass) {
            if (par.trialClass) {

                var trialGenome = [];

                for (i = 0; i <  par.trialClass.length; ++i) { //this.singletonParams.trialClass.length; ++i) {
                    trialGenome = [];
                    var currTrial =  par.trialClass[i];  //this.singletonParams.trialClass[i];
                    for (var j = 0;j < currTrial.length; ++j) {
                        trialGenome.push(gPool.geneWithIndex(j).alleleKey(currTrial[j]));
                    }
                    this.indivs[i].genome = trialGenome;

                }
               // trialGenome.forEach(function(el,ind) {
                //    that.indivs[ind].genome = el;

               // });
                //that.indivs[0].genome = trialGenome;

            }
            setTimeout(this.callBackFun,50);
            //this.callBackFun();
        }

    };

    //var that = this;
    if (this.initSize > 0) {


        if (this.genNum > 0) { //ie re-genesising after long stagnation. Dp directly, not in background
            for (var i = 0;i < this.initSize;++i) {
                var indii = new Indiv(0, i); //(i + 1);
                this.indivs.push(indii);
            }
            return;

        }

        //console.log('about to iter. pos: ' + 0 + ' interval: ' + this.initSuze/10 + ' inits: ' + this.initSize);

        //should move this out of here. Shouldn't referenc dom objects, as Population is a model
        if (typeof document !== 'undefined') {
             document.getElementById("progress").style.visibility = 'visible';
        }
        var stepSize = Math.floor(this.initSize/10);
        if (stepSize == 0) {
            stepSize = this.initSize;
        }
        this.iterationF(0,stepSize,this.initSize);
        /* one technique - but may not finish one iter before next one
         globInterval = setInterval(function() {
         that.doABitMoreTest()},100);
         */
        //do {
        //  console.log('setting timeout. pos: ' + pos);

        // } while (pos < this.initSize);
    }
    /*
     for (var i = 0; i < this.initSize; ++i) {
     if (i % interval == 0)  {
     //window.setTimeout(dummyFunc,50);
     //dummyFunc();
     //window.setTimeout(dummyFunc(i),50);
     //window.requestAnimationFrame(function()   {
     //    dummyFunc(i,tot);
     //})
     //dummyFunc);
     }





     var indii = new Indiv(0); //(i + 1);

     this.indivs.push(indii);
     }
     */

    //clearInterval(intervalHandle);

};


Population.prototype.popSize = function()  {
    return this.indivs.length;

};

Population.prototype.mutate = function()    {
    //var randArray = ['A','B','C','D'];
    var startMutate = 0;
    if (this.elite) startMutate+=1;

    for (var i = startMutate; i < this.popSize();++i) {
        // var maybe = Math.floor(Math.random() * 8); //was 15 15); //was 30
        var maybe = Math.random();
        if (maybe <  this.mutRate) { //> 5)  {//was 13) { //was 28
            var rGene = Math.floor(Math.random() * gPool.numGenes());
            var g = gPool.genes[rGene];
            var rAllele = Math.floor(Math.random() * gPool.allelesPerGene);
            var keys = getKeys(g.alleles);
            this.indivs[i].genome[rGene] = keys[rAllele];  //gPool.alleles[rGene][rAllele];//randArray[rAllele];
        }
    }
};

Population.prototype.chooseRoulette = function(allCosts)    {

    var cumArray = allCosts[2];
    var cumTot = allCosts[2][allCosts[2].length-1];
    var rNum = Math.floor(Math.random() * cumTot);
    for (var i = 0;i < cumArray.length; ++i) {
        if (rNum < cumArray[i]) {
            return i;
        }
    }

};

Population.prototype.nextGeneration = function() {
    var numInPop = 0;

    var rouletteNumber;
    rouletteNumber = this.prevPop.popSize();
    if (this.elite) {
        rouletteNumber -= 1; //Reserve one place in pop for best individual
    }
    this.initSize = this.prevPop.initSize;

    var indii;

    var allCosts = this.prevPop.getAllCosts();


    if (this.elite) {
        var bestPrev = this.prevPop.getBestCost();
        indii = new Indiv(0,numInPop,this.prevPop.indivs[bestPrev[1]]);
        this.indivs.push(indii);
        ++numInPop;
    }


    for (var i = 0;i < rouletteNumber;++i)  {
        // indii = new Indiv(0,numInPop,this.prevPop.indivs[this.chooseRoulette(allCosts)]);
        indii = new Indiv(0,numInPop,this.prevPop.indivs[this.chooseRoulette(allCosts)],this.prevPop.indivs[this.chooseRoulette(allCosts)]);
        this.indivs.push(indii);
        ++numInPop;

    }
    this.mutate();

};

/*
if (this.initSize > 0) {  //will be zero in copy population
    if (this.genNum) {
        this.nextGeneration();

    }
    else {

        this.genesis();
    }
}
*/



Population.prototype.getAllCosts = function()   {
    var allCostsIndexes = [];
    var allCostsCosts = [];
    var allCostsCumCosts = [];
    var cumCost = 0;

    this.indivs.forEach(function(indii,i) {
        var c = indii.cost();
        cumCost +=c;
        allCostsIndexes.push(i);
        allCostsCosts.push(c);
        allCostsCumCosts.push(cumCost);

    });
    /*
     for (var i =0;i < this.popSize(); ++i)  {
     var c = this.indivs[i].cost();
     cumCost +=c;
     allCostsIndexes.push(i);
     allCostsCosts.push(c);
     allCostsCumCosts.push(cumCost);

     }
     */
    return [allCostsIndexes,allCostsCosts,allCostsCumCosts];

};


/*
 this.dummyFunc = function(i) {

 document.getElementById("progress").innerHTML =  'Loading..' + (i * 100.0 / this.initSize);
 return;
 }
 */

Population.prototype.getBestCost = function()   {
    var bestCost = 0;
    var bestInd =  0;

    this.indivs.forEach(function(indii,i)   {
        var c = indii.cost();
        if (c > bestCost)   {
            bestCost = c;
            bestInd = i;
        }

    });
    /*
     for (var i =0;i < this.popSize(); ++i)  {
     var c = this.indivs[i].cost();
     if (c > bestCost)   {
     bestCost = c;
     bestInd = i;
     }

     }
     */
    return ([bestCost,bestInd]);
};

Population.prototype.getHamDists = function() {
    var that = this;
    var hamDists  = [];
    this.indivs.forEach(function (indiv,i) {
        var indivHamDist = 0;
        for (var j = 0;j < that.popSize();++j) {
            if (i == j) {

            }
            else {
                indivHamDist += indiv.hamDist(that.indivs[j]);

            }

        }
        hamDists.push(indivHamDist / that.popSize());

    });
    return hamDists;
};

Population.prototype.copyPopulation = function()    {
    //this.prevPop =
    var newPop = new Population(this.callBackFun,0,0);
    newPop.initSize = this.initSize;
    newPop.genNum = this.genNum;
    newPop.mutRate = this.mutRate;
    newPop.indivs = [];
    newPop.elite = this.elite;
    newPop.round = this.round;
  //  newPop.singletonParams = this.singletonParams;
    var newIndiv;
    for (var i = 0;i < this.popSize();++i)  {
        newIndiv = this.indivs[i].copyIndiv();
        newPop.indivs.push(newIndiv);
    }
    return newPop;
};

Population.prototype.toSerial = function(isBestPop) {
    var out = '';
    this.indivs.forEach(function(el) {
        out+= isBestPop ? 'Best,' : 'Class,';
        out+= el.toSerial() + '\r\n';
    });
    return out;

};
Population.prototype.printPopulation = function (hamDists) {
    //console.log('xxxStart Print Pop');
    hamDists = hamDists || null;

    //var ret = 'Generation: ' + this.genNum;
    var ret = '';

    this.indivs.forEach(function(indii,i) {
        if (hamDists) {

            ret += indii.printIndiv(hamDists[i]);
        }
        else {
            ret += indii.printIndiv();
        }
        ret += ' ';
    });
    /*
     for (var i = 0; i < this.indivs.length; ++i) {
     ret += this.indivs[i].printIndiv();
     ret += ' ';

     }
     */
   // console.log('xxxEnd print pop');
    return ret;

};
/*
function Student(studentString) {
    spl = studentString.split(',');
    this.studentId = spl[0];
    this.name = spl[1];
    this.prefs = [spl[2],spl[3],spl[4],spl[5],spl[6]];
    this.excls = [spl[7],spl[8],spl[9]];

}


    this.indivs = [];
*/
