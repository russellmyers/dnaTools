/**
 * Created by RussellM on 16/12/2015.
 */

/*
var c_Bases = ['A','C','G','T'];
var c_NumBases = 4;
*/

//global used in tree recursion
var mfk = [];

var allMFK = []; //used in ltclump
var clumps = []; //used in ltclump

var freqArray = []; //used in ltclump
var bestFreqArray = []; //used in ltclump

var revComplAlreadyDone = [];

var kMersDone = {}


var mismatchSeqs = [];

var nodesVisited = 0; //debugging

var startTime;
/*
//not importing from bioutilities for some reason

function allKmers(k) {

    if (k < 1) {
        return []; //just in case called with zero
    }
    else if (k == 1) {
        return c_Bases;
    }
    else {
        var lower = allKmers(k-1);
        var current = [];

        c_Bases.forEach(function(b) {


            lower.forEach(function(el) {
                current.push(b + el);

            });

        });
        return current;
    }

}
*/


/*
function reverseComplement(dna) {
//Bioutility

    var revCompl = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'a': 't',
        'c': 'g',
        'g': 'c',
        't': 'a'

    };

    var rev = '';

    for (var i = dna.length-1;i >= 0;--i) {
        rev += revCompl[dna[i]];
    }

    return rev;
}


function hamDist(dna1,dna2) {
//Bioutility
    //assumes equal length, otherwise returns -1

    var dist = 0;

    if (dna1 === dna2) {
        return 0;
    }
    else if (dna1.length == dna2.length ) {
        for (var i = 0; i < dna1.length; ++i) {
            if (dna1[i] === dna2[i]) {

            }
            else {
                ++dist;
            }

        }
        return dist;

    }
    else {
        return -1;
    }
}


function kMersWithMaxDist(kmer,maxDist) {
//Bioutility

    var allDist = [];
    for (var ii = 0;ii <= maxDist;++ii) {
        var ret = kMersWithDist(kmer, ii);
        allDist = allDist.concat(ret);
    }

    return allDist;


}

function kMersWithDist(kmer,dist) {
//Bioutility

    var otherBases = {'A':['C','G','T'],'C':['A','G','T'],'G':['A','C','T'],'T':['A','C','G']};

    //magic next 3 lines :-)
    if (dist == 0) {
        return [kmer];
    }

    if (kmer.length == 1) {
        if (dist == 0) {
            return [kmer];
        }
        else if (dist == 1) {
            return otherBases[kmer];
        }
        else {
            return [];
        }


    }
    else {

        if (kmer.length < dist) {
            return [];
        }
        else {
            var firstBase = kmer.substring(0,1);
            var rest = kmer.substring(1);
            var others = otherBases[firstBase];

            var otherFirst = [];
            var fixedFirst = [];

            others.forEach(function(c) {
                var otherFirstOne = kMersWithDist(rest,dist-1);
                otherFirstOne = otherFirstOne.map(function(el) {
                    return c + el;
                });
                otherFirst = otherFirst.concat(otherFirstOne);
            });

            if (kmer.length > dist) {
                //can keep first char the same
                fixedFirst = kMersWithDist(rest, dist);
                fixedFirst = fixedFirst.map(function (el) {
                    return firstBase + el;
                });
            }

            return fixedFirst.concat(otherFirst);

        }
    }



}
*/

function findKmerInDnaStrings(dnaStrings,kmer,d,inclRevCompl) {

    var foundInAll = true;

/*
    dnaStrings.forEach(function(dnaStr) {
        var countAr =  countKMer(kmer,dnaStr,inclRevCompl,d);
        if ((countAr[0].length + countAr[1].length) > 0) {

        }
        else {
            foundInAll = false;
        }


    });

    return foundInAll;
*/
    return dnaStrings.map(function(dnaStr) {
        return countKMer(kmer,dnaStr,inclRevCompl,d);

    });

}
function motifSearch(dnaStrings,k,d,inclRevCompl,progCallback) {

    var motifs = {};

    var progThreshold = 1;

    dnaStrings.forEach(function(dnaStr,strNum) {


        if (progCallback) {
                if ((strNum % progThreshold == 0) || strNum == dnaStrings.length - 1) {
                    progCallback('Process string',strNum,dnaStrings.length - 1);
                }
        }
        for (var i = 0; i < dnaStr.length - k + 1;++i) {
            if (progCallback) {
                if ((i % progThreshold == 0) || i == dnaStr.length - k) {
                    progCallback('String: ' + strNum + ' Base: ',i,dnaStr.length - k);
                }
            }
            var kmer = dnaStr.substring(i,i+k);
            var allDist = kMersWithMaxDist(kmer,d);
            allDist.forEach(function(neighbour) {
                var foundInAll = true;
                var foundArray = findKmerInDnaStrings(dnaStrings,neighbour,d,inclRevCompl);
                foundArray.forEach(function(el) {
                    if (el[0].length + el[1].length < 1) {
                        foundInAll = false;
                    }


                });
                if (foundInAll) {
                    motifs[neighbour] = foundArray;
                    //motifs.push(neighbour);

                }
            });


        }


    });

    var newAr = [];

    for (var key in motifs) {
        if (motifs.hasOwnProperty(key)) {
            var locs = [motifs[key],key];
            //var locs = [motifs[key][0][0],motifs[key][0][1],key];
            //var entry = [[0],[0],key];
            newAr.push(locs);




        }

    }
    return newAr;

}


function processBrute(dnaStrings,k,d,laplace,progCallback) {

    var kmers = [];
    for (var i = 0;i < dnaStrings[0].length - k + 1; ++i) {
        var kmer = dnaStrings[0].substring(i,i+k);
        var restKmers = processBrute(dnaStrings.slice(1),k,d,laplace,progCallback);
    }

}


function bruteMotifSearch(dnaStrings,k,d,laplace,progCallback) {

    var motifs = {};

    var progThreshold = 100;

    var iAr = dnaStrings.map(function(el) {
        return 0;
    });

    var limit = dnaStrings[0].length - k;

    var last = iAr.length - 1;

    var finished = false;
    var currPos;
    var shuffling;

    var bestScore = 99999;
    var bestInds = [];
    var bestConsensus = '';

    var totLoops = Math.pow(dnaStrings[0].length - k + 1,dnaStrings.length);

    var thisLoop = 0;

    var kmers = iAr.map(function(el,i) {
        return dnaStrings[i].substring(el,el+k);

    });


    while (!finished) {

        if (progCallback) {
            if ((thisLoop % progThreshold == 0) || (thisLoop == totLoops - 1)) {
                progCallback('Processing',thisLoop + 1,totLoops);
                //('Processing',thisLoop +1, totLoops);
            }
        }

        var res = scoreMotifs(kmers,laplace);
        var score = res[1];
        var consensus = res[0];
        if (score < bestScore) {
            bestScore = score;
            bestInds =  iAr.map(function(el) {
                return el;
            });
            bestConsensus = consensus;
        }

        currPos = last;
        shuffling = true;
        while ((shuffling) && (currPos >= 0)) {
             ++iAr[currPos];
             if (iAr[currPos] > limit) {
                    iAr[currPos] = 0;
                    kmers[currPos] = dnaStrings[currPos].substring(iAr[currPos],iAr[currPos] + k);
                    --currPos;
                    if (currPos < 0) {
                        finished = true;
                    }

             }
             else {
                  kmers[currPos] = dnaStrings[currPos].substring(iAr[currPos],iAr[currPos] + k);
                  shuffling = false;
             }
        }

        ++thisLoop;
    }


    var res = [
        [
            [

            ],
            bestConsensus
        ]
    ];
    bestInds.forEach(function(el) {
        res[0][0].push([[el],[]]);
    });
    /*
     dnaStrings.forEach(function(el) {
     res[0][0].push([[0],[]]);

     });
     */
    return [res,bestScore];
}


function medianString(dnaStrings,k,d,inclRevCompl,progCallback) {

    var motifs = {};

    var progThreshold = 100;

    var r = reverseComplement('AAA');


    //var allPatterns = allKmers(k);

    var minDist = 999999;
    var minPosns;
    var median = '';

    // allPatterns.forEach(function(pat) {
    var limit = Math.pow(4, k);
    for (var i = 0; i < limit; ++i) {
        if (progCallback) {
            if ((i % progThreshold == 0) || (i == limit-1)) {
                progCallback('Processing pattern',i, limit - 1);
            }
        }
        var pat = indToKmer(i,k);
        var res = patternSequencesDist(pat, dnaStrings, true);
        if (res[0] < minDist) {
            minDist = res[0];
            median = pat;
            minPosns = res[1];

        }
    }
    //});

    /*
    var newAr = [];

    for (var key in motifs) {
        if (motifs.hasOwnProperty(key)) {
            var locs = [motifs[key],key];
            //var locs = [motifs[key][0][0],motifs[key][0][1],key];
            //var entry = [[0],[0],key];
            newAr.push(locs);




        }

    }
    */
    //return newAr;

    var res = [
        [
            [

            ],
            median
        ]
    ];
    minPosns.forEach(function(el) {
        res[0][0].push([[el],[]]);
    });
    /*
    dnaStrings.forEach(function(el) {
        res[0][0].push([[0],[]]);

    });
    */
    return [res,minDist];
}

function greedyMotif(dnaStrings,k,d,laplace,progCallback) {

    var motifs = {};

    var progThreshold = 1;


    var minDist = 999999;
    var minPosns;
    var greedy = '';


    var first = dnaStrings[0];

    var bestScore = 999999;
    var bestConsensus = '';
    var bestMotifs = [];
    var bestIndexes = [];

    for (var i = 0;i < first.length - k + 1;++i) {
        if (progCallback) {
            if ((i % progThreshold == 0) || i == first.length - k) {
                progCallback('Process string',i,first.length - k);
            }
        }
        var kmer = first.substring(i,i+k);
        var kmersSoFar = [kmer];
        var kmerIndsSoFar = [i];
        var profMatrix = scoreMotifs(kmersSoFar,laplace)[5];

        for (var j = 0;j < dnaStrings.length - 1;  ++j) {
            var mostProb = profileMostProbable(profMatrix,dnaStrings[j + 1],k);
            kmersSoFar.push(mostProb[0]);
            kmerIndsSoFar.push(mostProb[2]);
            profMatrix = scoreMotifs(kmersSoFar,laplace)[5];

        }
        var res = scoreMotifs(kmersSoFar,laplace);
        var consensus = res[0];
        var score = res[1];
        if (score < bestScore) {
            bestScore = score;
            bestConsensus = consensus;
            bestMotifs = kmersSoFar;
            bestIndexes = kmerIndsSoFar;
        }


    }

    var b = bestMotifs;

    var res = [
        [
            [

            ],
            bestConsensus
        ]
    ];

    /*
    minPosns = dnaStrings.map(function(el) {
        return 0;

    });
    */

    bestIndexes.forEach(function(el) {
        res[0][0].push([[el],[]]);
    });


    return [res,bestScore];
}


function randomMotif(dnaStrings,k,numIters,d,laplace,progCallback) {

    var motifs = {};

    var progThreshold = 1;


    var minDist = 999999;
    var minPosns;
    var greedy = '';

 //   var first = dnaStrings[0];


    var bestConsensus = '';
    var bestMotifs = [];
    var bestIndexes = [];


    var bestScore = 99999;
    var entropyOfBestScore = 99999;
    var bestKmers = [];
    var bestKmerInds = [];

    var kmers = [];
    var kmerInds = [];


    for (var i = 0; i < numIters; ++i) {
        if (progCallback) {
            if ((i % progThreshold == 0) || i == numIters - 1) {
                progCallback('Iteration', i, numIters - 1);
            }
        }

        kmers = [];
        kmerInds = [];

        dnaStrings.forEach(function(el) {
            var r = getRandomInt(0, el.length - k);
            kmers.push(el.substring(r,r+k));
            kmerInds.push(r);
        });

        var lastScore = 99999;
        var res = scoreMotifs(kmers,laplace);
        var score = res[1];
        var entropyScore = res[2];
        var profMatrix = res[5];
        var consensus = res[0];

        while (score < lastScore) {

            lastScore = score;

            var res = motifsMostProbable(profMatrix, dnaStrings, k);
            kmers = res.map(function (el) {
                return el[0];
            });
            kmerInds = res.map(function(el) {
                return el[2];

            });
            res = scoreMotifs(kmers,laplace);
            consensus = res[0];
            score = res[1];
            entropyScore = res[2];
            profMatrix = res[5];


        }

        if (score < bestScore) {
            bestScore = score;
            bestKmers = kmers;
            bestKmerInds = kmerInds;
            bestConsensus = consensus;
            entropyOfBestScore = entropyScore;
        }



    }

    var res = [
        [
            [

            ],
            bestConsensus
        ]
    ];

    /*
     minPosns = dnaStrings.map(function(el) {
     return 0;

     });
     */

    bestKmerInds.forEach(function(el) {
        res[0][0].push([[el],[]]);
    });


    return [res,bestScore,entropyOfBestScore];
}

function gibbsSampler(dnaStrings,k,numIters,d,laplace,progCallback) {

    var motifs = {};

    var progThreshold = 1;


    var minDist = 999999;
    var minPosns;
    var greedy = '';

    //   var first = dnaStrings[0];


    var bestConsensus = '';
    var bestMotifs = [];
    var bestIndexes = [];


    var bestScore = 99999;
    var bestKmers = [];
    var bestKmerInds = [];

    var kmers = [];
    var kmerInds = [];

    var numStarts = 2000;

    for (var s = 0;s < numStarts; ++s ) {
        if (progCallback) {
            if ((s % progThreshold == 0) || s == numStarts - 1) {
                progCallback('Starts', s, numStarts - 1);
            }
        }

        kmers = [];
        kmerInds = [];

        dnaStrings.forEach(function (el) {
            var r = getRandomInt(0, el.length - k);
            kmers.push(el.substring(r, r + k));
            kmerInds.push(r);
        });

        var res = scoreMotifs(kmers, laplace);
        var score = res[1];
        var profMatrix = res[5];
        var consensus = res[0];

        if (score < bestScore) {
            bestScore = score;
            bestConsensus = consensus;
            bestKmers = kmers;
            bestKmerInds = kmerInds;
        }

        for (var i = 0; i < numIters; ++i) {

            var randRow = getRandomInt(0, dnaStrings.length - 1);
            kmers.splice(randRow, 1);
            kmerInds.splice(randRow, 1);
            var res = scoreMotifs(kmers, laplace);
            var profMatrix = res[5];
            var res = profileProbDist(profMatrix, dnaStrings[randRow], k);
            var kmer = res[0];
            var ind = res[2];
            kmers.splice(randRow, 0, kmer);
            kmerInds.splice(randRow, 0, ind);
            var res = scoreMotifs(kmers, laplace);
            consensus = res[0];
            score = res[1];
            if (score < bestScore) {
                bestScore = score;
                bestKmers = kmers;
                bestKmerInds = kmerInds;
                bestConsensus = consensus;
            }


        }
    }
    /*
        for (var i = 0; i < numIters; ++i) {
            if (progCallback) {
                if ((i % progThreshold == 0) || i == numIters - 1) {
                    progCallback('Iteration', i, numIters - 1);
                }
            }

            kmers = [];
            kmerInds = [];

            dnaStrings.forEach(function(el) {
                var r = getRandomInt(0, el.length - k);
                kmers.push(el.substring(r,r+k));
                kmerInds.push(r);
            });

            var lastScore = 99999;
            var res = scoreMotifs(kmers,laplace);
            var score = res[1];
            var profMatrix = res[5];
            var consensus = res[0];

            while (score < lastScore) {

                lastScore = score;

                var res = motifsMostProbable(profMatrix, dnaStrings, k);
                kmers = res.map(function (el) {
                    return el[0];
                });
                kmerInds = res.map(function(el) {
                    return el[2];

                });
                res = scoreMotifs(kmers,laplace);
                consensus = res[0];
                score = res[1];
                profMatrix = res[5];


            }

            if (score < bestScore) {
                bestScore = score;
                bestKmers = kmers;
                bestKmerInds = kmerInds;
                bestConsensus = consensus;
            }



        }
        */

    var res = [
        [
            [

            ],
            bestConsensus
        ]
    ];

    /*
     minPosns = dnaStrings.map(function(el) {
     return 0;

     });
     */

    bestKmerInds.forEach(function(el) {
        res[0][0].push([[el],[]]);
    });


    return [res,bestScore];
}


function countKMer(kmer,dna,inclRevCompl,maxMismatch,progCallback) {

    if (kmer === 'GCACACAGAC') {
        console.log('aha, processing: ' + kmer);
    }
    var progThreshold = 10000;

    var k = kmer.length;
    indexArray = [];
    var revIndexArray = [];

    var revComplKmer;

    maxMismatch = maxMismatch || 0;

    inclRevCompl = inclRevCompl || false;

    if (inclRevCompl) {
        revComplKmer = reverseComplement(kmer);

    }
    var count = 0;

    for (var j = 0; j < dna.length - k + 1; ++j) {
        if (progCallback) {
            if (j % progThreshold == 0) {
                progCallback(j, dna.length - k + 1);
            }
        }
        var tester = dna.substring(j, j + k);


        /*
         if ((tester === kmer) ||
         ( inclRevCompl && (tester === revComplKmer)))
         {
         ++count;
         indexArray.push(j);
         //console.log('found: ' + kmer);

         }
         */
        // if (tester === kmer) {
        if (hamDist(tester, kmer) <= maxMismatch) {
            if (kmer === 'acctgggatc') {
                console.log('gatcc...');
            }
            indexArray.push(j);
        }
        // else if (tester === revComplKmer) {
        if (inclRevCompl) {
            /*
            if (revComplKmer === kmer) {

            }

            else {
            */
                if (hamDist(tester, revComplKmer) <= maxMismatch) {
                    if (kmer === 'acctgggatc') {
                        console.log('gatcc... rev');
                    }
                    if (tester === 'gatcccaggt') {
                        console.log('gatcc..ggt rev');
                    }
                    revIndexArray.push(j);
                }

            /*
            }
            */
        }

    }


    if (kmer === 'acctgggatc') {
        console.log('ind arr: ' + indexArray);
        console.log('rev: ' + revIndexArray);
    }

    if (totIndArrayLen([indexArray,revIndexArray]) >= 4) {
        //console.log('in count. at least 4. kmer is: ' + kmer);
    }
    return [indexArray,revIndexArray,kmer];

}

function checkEqual(ar1,ar2) {
    if (ar1.length != ar2.length) {
        return false;
    }

    for (var i = 0; i < ar1.length; ++i) {
        if (ar1[i] === ar2[i]) {

        }
        else {
            return false;
        }
    }

    return true;


}

function cleanUpIdentical(mostArrays) {
    var delIndexes = [];
    var sortedArrays = mostArrays.map(function(el) {
        var newEl = el[0].concat(el[1]); //combine indexes and revindexes
        newEl.sort();
        return [newEl,el[2]];
    });
    for (var i = 0;i < sortedArrays.length; ++i) {
        for (var j = i+1;j < sortedArrays.length;++j) {
            //if (sortedArrays[i] === sortedArrays[j]) {
            if (checkEqual(sortedArrays[i][0],sortedArrays[j][0])) {
                if (sortedArrays[i][1] === sortedArrays[j][1]) { //also check if string equals
                    delIndexes.push(j);
                }

            }
        }
    }

   var newArrays = mostArrays.filter(function(el,i) {
       if (delIndexes.indexOf(i) > -1) {
           return false;

       }
       else {
           return true;
       }

   });
    return newArrays;

}

function mostFrequentKMersSort(k,dna,inclRevCompl,maxMismatch,ltClumpThreshold,debug,progCallback) {

    var progThreshold = 100;


    var freqArray = [];
    /*
    for (var i = 0;i <= Math.pow(c_NumBases,k) - 1;++i) {
        freqArray[i] = [];
    }
*/
    var kMersDone = [];
    var mostCount = 0;

    mostArray = [];

    var kMersDone = {}; //trying object instead of array

    maxMismatch = maxMismatch || 0;

    inclRevCompl = inclRevCompl || false;

    debug = debug || false;


    for (var i = 0;i < dna.length - k + 1; ++i) {
        if (progCallback) {
            if ((i % progThreshold == 0) || (i == dna.length - k)) {
                progCallback('Build Freq',i,dna.length - k);
            }
        }

        var count = 0;
        var kmer = dna.substring(i,i+k);

        var allDist = [];
        for (var ii = 0;ii <= maxMismatch;++ii) {
            var ret = kMersWithDist(kmer, ii);
            allDist = allDist.concat(ret);
        }

        /*
         if (inclRevCompl) {
         for (var ii = 0;ii <= maxMismatch;++ii) {
         var ret = kMersWithDist(reverseComplement(kmer), ii);
         allDist = allDist.concat(ret);
         }

         }
         */

        for (var j = 0;j < allDist.length;++j) {

            var ind = kMerToInd(allDist[j]);
            if (inclRevCompl) {
                var rInd = kMerToInd(reverseComplement(allDist[j]));
                freqArray.push([Math.min(ind,rInd),i,ind]);
                if (ind === rInd) { //palindrome
                    freqArray.push([Math.min(ind,rInd),i,-1]);

                }
            }
            else {
                freqArray.push([ind,i,ind]);
            }
            //freqArray[ind].push(i);


        }





    }

    if (progCallback) {

        progCallback('Sorting',1,1);

    }

    freqArray.sort(function(a,b) {
         return a[0] - b[0];
    });

    var countArray = [];

    var prev = -1;

    var currCount = 0;
    var currPosns = [[],[]];

    freqArray.forEach(function(el,i) {
        if (progCallback) {
            if ((i % progThreshold == 0) || (i == freqArray.length - 1)) {
                progCallback('Process count',i,freqArray.length - 1);
            }
        }
        if (el[0] != prev) {
            prev = el[0];
            currCount = 1;
            if (el[0] === el[2]) {
                currPosns = [[el[1]],[]];
            }
            else {
                currPosns = [[],[el[1]]];
            }

        }
        else {
            ++currCount;
            if (el[0] == el[2]) {
                var c = currPosns[0].concat(el[1]);
                currPosns = [c,currPosns[1]];
            }
            else {
                var c = currPosns[1].concat(el[1]);
                currPosns = [currPosns[0],c];
            }
        }
        countArray.push([currCount,currPosns]);


    });

    var thresh;
    if (ltClumpThreshold > 0 ) {
        thresh = ltClumpThreshold;

    }
    else {
        thresh = -1;
        countArray.forEach(function(el,i) {
            if (el[0] > thresh) {
                thresh = el[0];
            }

        });
        //max = Array.max(countArray);
    }



    var mostLen = -1;

    var newFreqArray = [];

    countArray.forEach(function(el,i) {
        if (progCallback) {
            if ((i % progThreshold == 0) || (i == countArray.length - 1)) {
                progCallback('Consol count',i,countArray.length - 1);
            }
        }
        if (
            (el[0] >= thresh)
            && ( (i == countArray.length - 1) || (countArray[i+1][0] <= el[0]))) { //ensure if ltclumpthresh entered, only best is counted
            if (inclRevCompl) {
                var km = indToKmer(freqArray[i][0], k);
                if (reverseComplement(km) === km) { //palindrome
                    //  if (freqArray[i][0] === freqArray[i][2]) {//palindrome
                    newFreqArray.push([el[1][0], el[1][1], indToKmer(freqArray[i][0], k)]);
                }
                else {

                    if ((el[1][0].length > 0) && (el[1][1].length > 0)) {
                        km = dna.substring(el[1][0][0],el[1][0][0] + k);
                        //newFreqArray.push([el[1][0], el[1][1], indToKmer(freqArray[i][0], k)]);
                        newFreqArray.push([el[1][0], el[1][1],km]);
                        km = dna.substring(el[1][1][0],el[1][1][0] + k);

                       // newFreqArray.push([el[1][1], el[1][0], indToKmer(freqArray[i][2], k)]);
                        newFreqArray.push([el[1][1], el[1][0], km]);
                    }
                    else if (el[1][0].length > 0) {
                        km = dna.substring(el[1][0][0],el[1][0][0] + k);
                        newFreqArray.push([el[1][0],[],km ]);
                    }
                    else {
                        km = dna.substring(el[1][1][0],el[1][1][0] + k);
                        newFreqArray.push([el[1][1],[], km]);

                    }

                }


            }
            else {
                newFreqArray.push([el[1][0], el[1][1], indToKmer(freqArray[i][0], k)]);

            }

        }

    });



    //console.log('most array: ' + mostArray);

    //console.log('most arrr: ' + mostArray);
    /*
     var highest = -1;
     var mostArray = [];
     freqArray.forEach(function(el,i) {
     if ( el > highest) {
     mostArray = [];
     mostArray.push()

     }


     });
     for (var i = 0;i < freqArray.length;++i) {
     */



    // return cleaned;

    var cleaned = cleanUpIdentical(newFreqArray);

    //return newFreqArray;

    if (debug) {
        return [cleaned,[]];
    }
    else {
        return cleaned;
    }


}

function buildBestFreqArray(freqArray,progCallback) {

    var progThreshold = 1;

    var bestFreqArray = freqArray.map(function(el,i) {
        //console.log('i: ' + i);
        if (progCallback) {
            if ((i % progThreshold == 0) || i == freqArray.length - 1) {
                progCallback('Build Best Freq',i,freqArray.length - 1);
            }
        }

        var newEl0,newE1;
       if (el.length == 0) {
           return [];
       }
       if (el[0]) {
           newEl0 = el[0].map(function (e) {
               return e;

           }) ;

       }
        else {
           newEl0 = [];
       }


        if (el[1]) {
            newEl1 = el[1].map(function (e) {
                return e;

            }) ;


        }
        else {
            newEl1 = [];
        }

        return[newEl0,newEl1];

    });

    return bestFreqArray;

}



function deleteKmerFromFreqArray(kmer,pos,freqArray,inclRevCompl) {
    var ind = kMerToInd(kmer);


    var checkPos = freqArray[ind][0][0];

    if (checkPos == pos) {
        freqArray[ind][0].splice(0,1);
    }

    if (inclRevCompl) {
        var revInd = kMerToInd(reverseComplement(kmer));
        var checkRevPos = freqArray[revInd][1][0];
        if (checkRevPos == pos) {
            freqArray[revInd][1].splice(0, 1);
        }
    }

    return freqArray;



}

function addKmerToFreqArray(kmer,pos,freqArray,bestFreqArray,inclRevCompl) {

    var ind = kMerToInd(kmer);

    if (!freqArray[ind]) freqArray[ind] = [];

    if (!freqArray[ind][0]) {
        freqArray[ind][0] = [];
    }
    freqArray[ind][0].push(pos);

    //if ( (freqArray[ind][0].length + freqArray[ind][1].length) > (bestFreqArray[ind][0].length + bestFreqArray[ind][1].length) ) {
    if (totIndArrayLen(freqArray[ind]) >= totIndArrayLen(bestFreqArray[ind])) {
        var newEl0;
        if (freqArray[ind][0]) {

            newEl0 = freqArray[ind][0].map(function (el) {
                return el;

            });
        }
        else {
            newEl0 = [];
        }

        var newEl1;
        if (freqArray[ind][1]) {


            newEl1 = freqArray[ind][1].map(function (el) {
                return el;

            });
        }
        else {
            newEl1 = [];
        }
        bestFreqArray[ind] = [newEl0,newEl1];
    }


    if (inclRevCompl) {
        var revInd = kMerToInd(reverseComplement(kmer));
        if (!freqArray[revInd]) freqArray[revInd] = [];

        if (!freqArray[revInd][1]) {
            freqArray[revInd][1] = [];
        }
        freqArray[revInd][1].push(pos);

       // if ( (freqArray[revInd][0].length + freqArray[revInd][1].length) > (bestFreqArray[revInd][0].length + bestFreqArray[revInd][1].length) ) {
       if (totIndArrayLen(freqArray[revInd]) >= totIndArrayLen(bestFreqArray[revInd])) {
           var newEl0;
           if (freqArray[revInd][0]) {


               newEl0 = freqArray[revInd][0].map(function (el) {
                   return el;

               });
           }
           else {
               newEl0 = [];
           }

           var newEl1;
           if (freqArray[revInd][1]) {


               newEl1 = freqArray[revInd][1].map(function (el) {
                   return el;

               });
           }
           else {
               newEl1 = [];

           }
           bestFreqArray[revInd] = [newEl0,newEl1];
        }

    }

    return [freqArray,bestFreqArray];


}




function buildFreqArray(k,dna,inclRevCompl,maxMismatch,ltClumpThreshold,debug,progCallback) {

    var progThreshold = 100;

    var freqArray = [];
    //freqArray = [<array of posns>,<array of rev compl posns>] - one for each kmer hash index, eg AA = 0, AC = 1 etc

    /*
    for (var i = 0;i <= Math.pow(c_NumBases,k) - 1;++i) {
        if (progCallback) {
            if (i % progThreshold == 0) {
                progCallback('Build Freq 1',i,Math.pow(c_NumBases,k) - 1);
            }
        }
        //freqArray[i] = [[],[]];

        freqArray[i] = [];

    }
    */

    progThreshold = 100;
    for (var i = 0;i < dna.length - k + 1; ++i) {
        if (progCallback) {
            if ((i % progThreshold == 0) || i == dna.length - k) {
                progCallback('Build Freq',i,dna.length - k);
            }
        }

        var kmer = dna.substring(i,i+k);

        /*
        var allDist = [];
        for (var ii = 0;ii <= maxMismatch;++ii) {
            var ret = kMersWithDist(kmer, ii);
            allDist = allDist.concat(ret);
        }
        */
        var allDist = kMersWithMaxDist(kmer,maxMismatch);


        for (var j = 0;j < allDist.length;++j) {

            var ind = kMerToInd(allDist[j]);
            if (!freqArray[ind]) {
                freqArray[ind] = [];
            }
            if (!freqArray[ind][0]) {
                freqArray[ind][0] = [];
            }
            freqArray[ind][0].push(i);


            if (inclRevCompl) {
                var revKM = reverseComplement(allDist[j]);
                var revInd = kMerToInd(revKM);
                if (!freqArray[revInd]) freqArray[revInd] = [];


                if (!freqArray[revInd][1]) {
                    freqArray[revInd][1] = [];
                }
                freqArray[revInd][1].push(i);
            }

        }

    }

    return freqArray;


}


function buildFreqObject(k,dna,inclRevCompl,maxMismatch,ltClumpThreshold,debug,progCallback) {

    var progThreshold = 10;

    var freqObjectKeyCount = 0;

    var freqObject = {};
    //freqArray = [<array of posns>,<array of rev compl posns>] - one for each kmer hash index, eg AA = 0, AC = 1 etc

    /*
     for (var i = 0;i <= Math.pow(c_NumBases,k) - 1;++i) {
     if (progCallback) {
     if (i % progThreshold == 0) {
     progCallback('Build Freq 1',i,Math.pow(c_NumBases,k) - 1);
     }
     }
     //freqArray[i] = [[],[]];

     freqArray[i] = [];

     }
     */

    progThreshold = 10;
    for (var i = 0;i < dna.length - k + 1; ++i) {
        if (progCallback) {
            if ((i % progThreshold == 0) || i == dna.length - k) {
                progCallback('Build Freq',i,dna.length - k);
            }
        }

        var kmer = dna.substring(i,i+k);

        /*
         var allDist = [];
         for (var ii = 0;ii <= maxMismatch;++ii) {
         var ret = kMersWithDist(kmer, ii);
         allDist = allDist.concat(ret);
         }
         */
        var allDist = kMersWithMaxDist(kmer,maxMismatch);


        for (var j = 0;j < allDist.length;++j) {

            var ind = kMerToInd(allDist[j]);

            if (!freqObject[ind]) {
                freqObject[ind] = {};
                ++freqObjectKeyCount;
            }
            if (!freqObject[ind].posns) {
                freqObject[ind].posns =  [];
            }
            freqObject[ind].posns.push(i);


            if (inclRevCompl) {
                var revKM = reverseComplement(allDist[j]);
                var revInd = kMerToInd(revKM);
                if (!freqObject[revInd]) {
                    freqObject[revInd] = {};
                    ++freqObjectKeyCount;
                }


                if (!freqObject[revInd].revPosns) {
                    freqObject[revInd].revPosns = [];
                }
                freqObject[revInd].revPosns.push(i);
            }

        }

    }

    return [freqObject,freqObjectKeyCount];

}



function mostFrequentKMersFreqArray(k,dna,inclRevCompl,maxMismatch,ltClumpThreshold,debug,progCallback) {

    var progThreshold = 100;

    maxMismatch = maxMismatch || 0;

    inclRevCompl = inclRevCompl || false;

    debug = debug || false;

    var freqArray = buildFreqArray(k, dna, inclRevCompl, maxMismatch, ltClumpThreshold, debug, progCallback);
    // var bestFreqArray = buildBestFreqArray(freqArray);


    /*
     var freqArray = [];

     for (var i = 0;i <= Math.pow(c_NumBases,k) - 1;++i) {
     freqArray[i] = [];
     }

     var kMersDone = [];
     var mostCount = 0;

     mostArray = [];

     var kMersDone = {}; //trying object instead of array



     for (var i = 0;i < dna.length - k + 1; ++i) {
     if (progCallback) {
     if (i % progThreshold == 0) {
     progCallback(i,dna.length - k + 1);
     }
     }

     var count = 0;
     var kmer = dna.substring(i,i+k);

     var allDist = [];
     for (var ii = 0;ii <= maxMismatch;++ii) {
     var ret = kMersWithDist(kmer, ii);
     allDist = allDist.concat(ret);
     }


     for (var j = 0;j < allDist.length;++j) {

     var ind = kMerToInd(allDist[j]);
     freqArray[ind].push(i);

     }





     }

     */
    var newFreqArray = [];


    var mostLen = -1;

    freqArray.forEach(function (el, i) {

        if (progCallback) {
            if ((i % progThreshold == 0) || i == freqArray.length - 1) {
                progCallback('Consol most freq', i, freqArray.length - 1);
            }
        }

        if (el) {

            var km = indToKmer(i, k);
            /*


             if (inclRevCompl) {
             var revKM = reverseComplement(km);
             rEl = freqArray[kMerToInd(revKM)];
             }
             else {
             rEl = [];
             }
             */
            var totLen = totIndArrayLen(el);
            // var totLen = el[0].length + el[1].length;


            if (ltClumpThreshold > 0) {
                if (totLen >= ltClumpThreshold) {
                    if (!el[0]) el[0] = [];
                    if (!el[1]) el[1] = [];

                    newFreqArray.push([el[0], el[1], km]);
                }
            }
            else {
                if (totLen > mostLen) {
                    newFreqArray = [];
                    if (!el[0]) el[0] = [];
                    if (!el[1]) el[1] = [];
                    newFreqArray.push([el[0], el[1], km]);
                    mostLen = totLen;
                }
                else if (totLen == mostLen) {
                    if (!el[0]) el[0] = [];
                    if (!el[1]) el[1] = [];
                    newFreqArray.push([el[0], el[1], km]);

                }
                else {

                }
            }
        }

    });



    //console.log('most array: ' + mostArray);

    //console.log('most arrr: ' + mostArray);
    /*
     var highest = -1;
     var mostArray = [];
     freqArray.forEach(function(el,i) {
     if ( el > highest) {
     mostArray = [];
     mostArray.push()

     }


     });
     for (var i = 0;i < freqArray.length;++i) {
     */


    // var cleaned = cleanUpIdentical(mostArray);

    // return cleaned;
    if (debug) {
        var freqCountArray = freqArray.map(function(el) {

            // return el[0].length + el[1].length;
            return (totIndArrayLen(el));

        });
        return [newFreqArray,freqCountArray];
    }
    else {
        return newFreqArray;
    }

}



function mostFrequentKMersObject(k,dna,inclRevCompl,maxMismatch,ltClumpThreshold,debug,progCallback) {

    var progThreshold = 100;

    maxMismatch = maxMismatch || 0;

    inclRevCompl = inclRevCompl || false;

    debug = debug || false;

    var res = buildFreqObject(k, dna, inclRevCompl, maxMismatch, ltClumpThreshold, debug, progCallback);
    var freqObj = res[0];
    var freqCount = res[1];

  //  var freqArray = buildFreqArray(k, dna, inclRevCompl, maxMismatch, ltClumpThreshold, debug, progCallback);
    // var bestFreqArray = buildBestFreqArray(freqArray);


    /*
     var freqArray = [];

     for (var i = 0;i <= Math.pow(c_NumBases,k) - 1;++i) {
     freqArray[i] = [];
     }

     var kMersDone = [];
     var mostCount = 0;

     mostArray = [];

     var kMersDone = {}; //trying object instead of array



     for (var i = 0;i < dna.length - k + 1; ++i) {
     if (progCallback) {
     if (i % progThreshold == 0) {
     progCallback(i,dna.length - k + 1);
     }
     }

     var count = 0;
     var kmer = dna.substring(i,i+k);

     var allDist = [];
     for (var ii = 0;ii <= maxMismatch;++ii) {
     var ret = kMersWithDist(kmer, ii);
     allDist = allDist.concat(ret);
     }


     for (var j = 0;j < allDist.length;++j) {

     var ind = kMerToInd(allDist[j]);
     freqArray[ind].push(i);

     }





     }

     */
    var newFreqArray = [];


    var mostLen = -1;

    var i = 0;

    for (var key in freqObj) {
        if (freqObj.hasOwnProperty(key)) {


            if (progCallback) {
                if ((i % progThreshold == 0) || i == freqCount - 1) {
                    progCallback('Consol most freq obj', i, freqCount - 1);
                }
            }
            ++i;

            var km = indToKmer(key, k);



            var lenPos = freqObj[key].posns ? freqObj[key].posns.length : 0;
            var lenRevPos = freqObj[key].revPosns ? freqObj[key].revPosns.length : 0;
            var totLen = lenPos + lenRevPos;

            if (ltClumpThreshold > 0) {
                if (totLen >= ltClumpThreshold) {
                    if (!freqObj[key].posns) freqObj[key].posns = [];
                    if (!freqObj[key].revPosns) freqObj[key].revPosns = [];

                    newFreqArray.push([freqObj[key].posns,freqObj[key].revPosns, km]);
                }
            }
            else {
                if (totLen > mostLen) {
                    newFreqArray = [];
                    if (!freqObj[key].posns) freqObj[key].posns = [];
                    if (!freqObj[key].revPosns) freqObj[key].revPosns = [];

                    newFreqArray.push([freqObj[key].posns,freqObj[key].revPosns, km]);
                    mostLen = totLen;
                }
                else if (totLen == mostLen) {
                    if (!freqObj[key].posns) freqObj[key].posns = [];
                    if (!freqObj[key].revPosns) freqObj[key].revPosns = [];
                    newFreqArray.push([freqObj[key].posns,freqObj[key].revPosns, km]);

                }
                else {

                }
            }
        }

        }




    //console.log('most array: ' + mostArray);

    //console.log('most arrr: ' + mostArray);
    /*
     var highest = -1;
     var mostArray = [];
     freqArray.forEach(function(el,i) {
     if ( el > highest) {
     mostArray = [];
     mostArray.push()

     }


     });
     for (var i = 0;i < freqArray.length;++i) {
     */


    // var cleaned = cleanUpIdentical(mostArray);

    // return cleaned;
    if (debug) {
        var freqCountArray = freqArray.map(function(el) {

            // return el[0].length + el[1].length;
            return (totIndArrayLen(el));

        });
        return [newFreqArray,freqCountArray];
    }
    else {
        return newFreqArray;
    }

}



function mostFrequentKMers(k,dna,inclRevCompl,maxMismatch,ltClumpThreshold,debug,progCallback) {

    //If ltClumpThreshold is blank: only return most frequent kmers
    //otherwise: return all kmers which occur >= ltClumpThreshold

    var progThreshold = 1;


    var kMersDone = [];
    var mostCount = 0;

    mostArray = [];


    var kMersDone = {}; //trying object instead of array

    debug = debug || false;


    maxMismatch = maxMismatch || 0;

    inclRevCompl = inclRevCompl || false;


    for (var i = 0;i < dna.length - k + 1; ++i) {
        if (progCallback) {
            if (i % progThreshold == 0) {
                progCallback('Processing',i,dna.length - k + 1);
            }
        }

        var count = 0;
        var kmer = dna.substring(i,i+k);


        if (kMersDone[kmer] === undefined) {

        }
        else {
            ++kMersDone[kmer];
            continue;
        }

        /* following is logic to skip if already done. Now replaced with a cleanup at end
        var revComplKmer;
        if (inclRevCompl) {
            revComplKmer = reverseComplement(kmer);
        }

        var foundAlready = false;
        for (var ii=0; ii < kMersDone.length;++ii) {//ii<dna.length; ii++) {
            if (
                (kmer === kMersDone[ii])
                ||
                (inclRevCompl && (revComplKmer === kMersDone[ii]))
            )
            {
                foundAlready = true;
                break;

            }
            else {

            }
        }
        if (foundAlready) {
            continue;
        }
        else {
            kMersDone.push(kmer);
        }

        */

        //console.log('kmer: ' + kmer);

        /*
         indexArray = [];
         for (var j = 0;j < dna.length - k + 1; ++j) {
         var tester = dna.substring(j,j+k);
         if ((tester === kmer) ||
         ( inclRevCompl && (tester === revComplKmer)))
         {
         ++count;
         indexArray.push(j);
         //console.log('found: ' + kmer);

         }
         }
         */


        var allDist = [];
        for (var ii = 0;ii <= maxMismatch;++ii) {
            var ret = kMersWithDist(kmer, ii);
            allDist = allDist.concat(ret);
        }
        if (inclRevCompl) {
            for (var ii = 0;ii <= maxMismatch;++ii) {
                var ret = kMersWithDist(reverseComplement(kmer), ii);
                allDist = allDist.concat(ret);
            }

        }

       for (var j = 0;j < allDist.length;++j) {


           if (kMersDone[allDist[j]] === undefined) {
               kMersDone[allDist[j]] = 1;

           }
           else {
               ++kMersDone[allDist[j]];
               continue;
           }


           var indexArrays = countKMer(allDist[j], dna, inclRevCompl, maxMismatch);

           if (totIndArrayLen(indexArrays) >= 4) {
               // console.log('ahha ' + totIndArrayLen(indexArrays));
           }

           if (ltClumpThreshold > 0) {
               if ((totIndArrayLen(indexArrays) >= ltClumpThreshold)) { //(indexArray.length >= ltClumpThreshold) {
                   mostArray.push(indexArrays);
               }
           }
           else {
               if ((totIndArrayLen(indexArrays) == mostCount)) {
                   mostArray.push(indexArrays);

               }
               else if ((totIndArrayLen(indexArrays) > mostCount)) {
                   mostArray = [];
                   mostArray.push(indexArrays);
                   mostCount = totIndArrayLen(indexArrays);
               }
           }
       }



    }
    //console.log('most array: ' + mostArray);

    //console.log('most arrr: ' + mostArray);
    var cleaned = cleanUpIdentical(mostArray);

    if (debug) {
        return [cleaned,[]];
    }
    else {
        return cleaned;
    }


}

//Tree routines

function Node(letter) {
    this.letter = letter;
    this.sequence = letter;
    this.positions = [];

    this.subNodes = [];

}


function buildTree(dna,k,ltClumpThresh,includeRevCompl,progCallback) {


    var progThreshold = 10000;

    var justAddedPrev = [];
    var justAddedCurr = [];

    var baseNode = new Node('');

    justAddedCurr = [baseNode];
    for (var i = 0;i < dna.length; ++i) {
        if (progCallback) {
            if (i % progThreshold == 0) {
                progCallback('build',i,dna.length);
            }
        }


        justAddedPrev = [];
        justAddedCurr.forEach(function (e) {
            justAddedPrev.push(e);
        });
        justAddedCurr = [baseNode];



        justAddedPrev.forEach(function(ja) {

            if (ja.sequence.length > k) {

            }
            else {

                if (ja.subNodes.length == 0) {
                    var newNode = new Node(dna[i]);
                    newNode.positions = [i];
                    newNode.sequence = ja.sequence + newNode.letter;
                    ja.subNodes.push(newNode);
                    justAddedCurr.push(newNode);

                }
                else {
                    var snFound = false;
                    var snInd = -1;


                    ja.subNodes.forEach(function (sn, snNum) {
                        if (sn.letter === dna[i]) {
                            snFound = true;
                            snInd = snNum;
                        }
                    });
                    if (snFound) {
                        var sn = ja.subNodes[snInd];
                        sn.positions.push(i);
                        justAddedCurr.push(sn);
                    }
                    else {
                        var newNode = new Node(dna[i]);
                        newNode.positions = [i];
                        newNode.sequence = ja.sequence + newNode.letter;
                        ja.subNodes.push(newNode);
                        justAddedCurr.push(newNode);

                    }


                }
            }


        });


    }

    return [baseNode,justAddedCurr];

}


function totIndArrayLen(indArrays) {

    var tot = 0;
    if (!indArrays) {
        return 0;
    }

    var len0 = indArrays[0] ? indArrays[0].length : 0;
    var len1 = indArrays[1] ? indArrays[1].length : 0;

    //return indArrays[0].length + indArrays[1].length;
    return len0 + len1;

    /*
    indArrays.forEach(function(el) {
        tot+= el.length;

    });
    return tot;
    */

    /*
    var ret = indArrays.reduce(function(el1,el2) {
        var tot = 0;
        if (el1) {
            tot+= el1.length;

        }
        if (el2) {
            tot+= el2.length;
        }
        return tot;
    });
    return ret;
    */
}

function traverseTreeOneNode(baseNode,seq,ltClumpThresh,inclRevCompl,maxMismatch) {
    var posnsToUse = [];


    /*
     node.positions.forEach(function(el) {
     posnsToUse.push(el);

     });
     */


    //var indArrays = [indArray];
    var indArrays = [[],[],seq];

   // if (maxMismatch > 0) {
        mismatchSeqs = [];
        findSequenceInTreeWithMismatches(baseNode,seq,maxMismatch,0);
        mismatchSeqs =  mismatchSeqs.map(function (n) {
            return n - seq.length + 1;

        });
        indArrays[0] = indArrays[0].concat(mismatchSeqs);
   // }



    if (inclRevCompl) {
        // console.log('rev compl');
        var revCompl = reverseComplement(seq);


        /*
         if (revCompl === node.sequence) { //palindrome
         indArrays.push([]);
         }

         else {
         */

        var revComplPosns = [];
        //var revComplPosns = findSequenceInTree(baseNode, reverseComplement(node.sequence));
        //console.log('rev compl found in tree');

        mismatchSeqs = [];
        //if (maxMismatch > 0) {

            findSequenceInTreeWithMismatches(baseNode, revCompl, maxMismatch, 0);
            mismatchSeqs = mismatchSeqs.map(function (n) {
                return n - seq.length + 1;

            });
       // }// indArrays[0] = indArrays[0].concat(mismatchSeqs);
        /*
        if (!revComplPosns) {
            indArrays.push(mismatchSeqs);

        }
        else {

            revComplPosns.forEach(function (el) {
                posnsToUse.push(el);

            });
            revComplPosns = revComplPosns.map(function (n) {
                return n - node.sequence.length + 1;

            });
            revComplPosns = revComplPosns.concat(mismatchSeqs);
            indArrays.push(revComplPosns);
            revComplAlreadyDone.push(node.sequence);
        }
        */
        indArrays[1] = mismatchSeqs;
        revComplAlreadyDone.push(seq);
        //console.log('rev compl pushed');
        /*
         }
         */

    }
    /*
    else {
        indArrays.push([]);
    }
    */

    if (ltClumpThresh > 0) {
        if  (totIndArrayLen(indArrays) >= ltClumpThresh) { //(posnsToUse.length >= ltClumpThresh) {
            mfk.push(indArrays);
            /*mfk.push(posnsToUse.map(function (n) {
             return n - node.sequence.length + 1;

             }));
             */

        }
    }
    else {
        if (mfk.length > 0) {
            if  (totIndArrayLen(indArrays) > totIndArrayLen(mfk[0])) { //(posnsToUse.length > mfk[0].length) {
                mfk = [];
                mfk.push(indArrays);
                /*mfk.push(posnsToUse.map(function (n) {
                 return n - node.sequence.length + 1;

                 }));
                 */
            }
            else if (totIndArrayLen(indArrays)  == totIndArrayLen(mfk[0]) ) { //(posnsToUse.length == mfk[0].length) {
                /* mfk.push(posnsToUse.map(function (n) {
                 return n - node.sequence.length + 1;

                 }));
                 */
                mfk.push(indArrays);
            }
            else {

            }

        }
        else {
            /*
             mfk.push(posnsToUse.map(function (n) {
             return n - node.sequence.length + 1;

             }));
             */
            mfk.push(indArrays);

        }
    }


}
function traverseTree(baseNode,node,k,ltClumpThresh,inclRevCompl,maxMismatch,progCallback) {

    //uses globals mkf and revComplAlreadyDone during recursion

    maxMismatch = maxMismatch || 0;

    if (node.sequence.length == k) {
        //console.log('traversing....');
        ++nodesVisited;
        if (kMersDone[node.sequence] === undefined) {

        }
        else {
            ++kMersDone[node.sequence];
            return;
        }

        var allDist = [];
        for (var ii = 0;ii <= maxMismatch;++ii) {
            var ret = kMersWithDist(node.sequence, ii);
            allDist = allDist.concat(ret);
        }
        if (inclRevCompl) {
            for (var ii = 0;ii <= maxMismatch;++ii) {
                var ret = kMersWithDist(reverseComplement(node.sequence), ii);
                allDist = allDist.concat(ret);
            }

        }

        for (var j = 0;j < allDist.length;++j) {


            if (kMersDone[allDist[j]] === undefined) {
                kMersDone[allDist[j]] = 1;

            }
            else {
                ++kMersDone[allDist[j]];
                continue;
            }


            traverseTreeOneNode(baseNode, allDist[j], ltClumpThresh, inclRevCompl, maxMismatch);

        }

        /*


        var posnsToUse = [];

        posnsToUse = node.positions.map(function (n) {
            return n - node.sequence.length + 1;

        });


        var indArray = posnsToUse.map(function(el) {
            return el;
        });
        var indArrays = [indArray];

        if (maxMismatch > 0) {
            mismatchSeqs = [];
            findSequenceInTreeWithMismatches(baseNode,node.sequence,maxMismatch,0);
            mismatchSeqs =  mismatchSeqs.map(function (n) {
                return n - node.sequence.length + 1;

            });
            indArrays[0] = indArrays[0].concat(mismatchSeqs);
        }



        if (inclRevCompl) {
           // console.log('rev compl');
            var revCompl = reverseComplement(node.sequence);



                var revComplPosns = findSequenceInTree(baseNode, reverseComplement(node.sequence));
                //console.log('rev compl found in tree');

                mismatchSeqs = [];
                if (maxMismatch > 0) {

                    findSequenceInTreeWithMismatches(baseNode, reverseComplement(node.sequence), maxMismatch, 0);
                    mismatchSeqs = mismatchSeqs.map(function (n) {
                        return n - node.sequence.length + 1;

                    });
                }// indArrays[0] = indArrays[0].concat(mismatchSeqs);
                if (!revComplPosns) {
                    indArrays.push(mismatchSeqs);

                }
                else {

                    revComplPosns.forEach(function (el) {
                        posnsToUse.push(el);

                    });
                    revComplPosns = revComplPosns.map(function (n) {
                        return n - node.sequence.length + 1;

                    });
                    revComplPosns = revComplPosns.concat(mismatchSeqs);
                    indArrays.push(revComplPosns);
                    revComplAlreadyDone.push(node.sequence);
                }
                //console.log('rev compl pushed');


        }
        else {
            indArrays.push([]);
        }

        if (ltClumpThresh > 0) {
            if  (totIndArrayLen(indArrays) >= ltClumpThresh) { //(posnsToUse.length >= ltClumpThresh) {
                mfk.push(indArrays);


            }
        }
        else {
            if (mfk.length > 0) {
                if  (totIndArrayLen(indArrays) > totIndArrayLen(mfk[0])) { //(posnsToUse.length > mfk[0].length) {
                    mfk = [];
                    mfk.push(indArrays);

                }
                else if (totIndArrayLen(indArrays)  == totIndArrayLen(mfk[0]) ) { //(posnsToUse.length == mfk[0].length) {

                    mfk.push(indArrays);
                }
                else {

                }

            }
            else {

                mfk.push(indArrays);

            }
        }
        */

        return;

    }


    /*
     if ((node.sequence.length == k) && (node.positions.length >= ltClumpThresh)) {
     //mfk.push(node.positions);

     mfk.push(node.positions.map(function (n) {
     return n - node.sequence.length + 1;

     }));
     ///console.log('seq: ' + node.sequence + ' num posns: ' + node.positions.length);
     }
     */



    node.subNodes.forEach(function (sn,ind) {
        if (node.sequence.length == 0) {

            if (progCallback) {
                progCallback('traverse ' + sn.sequence + ' branch',ind + 1,node.subNodes.length);
            }
        }

        traverseTree(baseNode,sn,k,ltClumpThresh,inclRevCompl,maxMismatch,progCallback);

    });





}

function findSequenceInTreeWithMismatches(node,seq,maxMismatch,currMismatch) {


    if (node.sequence.length == seq.length) {
        var d = hamDist(node.sequence,seq);
        if (d == 0) {
            //used to skip this, as already counted. but now count here
            mismatchSeqs = mismatchSeqs.concat(node.positions);
        }
        else if (d <= maxMismatch) {
            mismatchSeqs = mismatchSeqs.concat(node.positions);
        }
        else {

        }
        return;

    }
    else {
        var nextChar = node.sequence.length;
        //var found = false;
        // var foundPos;
        var ret;
        for (var i = 0; i < node.subNodes.length; ++i) {

            if (node.subNodes[i].letter === seq[nextChar]) {

                ret = findSequenceInTreeWithMismatches(node.subNodes[i], seq, maxMismatch, currMismatch);
            }
            else {
                if (currMismatch + 1 <= maxMismatch) {
                    ret = findSequenceInTreeWithMismatches(node.subNodes[i], seq, maxMismatch, currMismatch + 1);
                }
                else {

                }



            }
        }

        return false;

    }

}



function findSequenceInTree(node,seq) {


    if (node.sequence === seq) {
        return node.positions;

    }
    else {
        var nextChar = node.sequence.length;
        //var found = false;
        // var foundPos;
        var ret;
        for (var i = 0; i < node.subNodes.length; ++i) {
            if (node.subNodes[i].letter === seq[nextChar]) {

                ret = findSequenceInTree(node.subNodes[i], seq);
                if (ret) {
                    return ret;
                }
                else {
                    return false;
                }


            }
        }

        return false;

    }

}

function ltClump(dna,l,t,k,curr,useTree,inclRevCompl,maxMismatch,progCallback) {
    // console.log('st timeout');

    //colourDNAProgressNew(dna, [curr], clumps,l);

    if (curr % 1000 == 1) {
        console.log('curr start: ' + curr);
    }

    var progThreshold = 100;

    var limit = dna.length - l;

    //document.getElementById('progress').innerHTML = curr + ' / ' + limit;

    var ltClumpThresh =  t;

    var mf;

    if (progCallback) {
        if (curr % progThreshold == 0) {
            progCallback('Clumps',curr,dna.length - l + 1);
        }

    }

    if (useTree) {
        // console.log('st build');
        var added;
        if (curr == 0) {
            added = buildTree(dna.substring(curr, curr + l),k,ltClumpThresh,inclRevCompl);
            baseNode = added[0];
            justAdded = added[1];
        }
        else {
            var seq = dna.substring(curr-1,curr-1 + k);
            delKmerFromTree(baseNode,seq);
            added = addSingleBaseToTree(justAdded,baseNode,dna.substring((curr+l-1),curr+l),curr+l-1,k);
            baseNode = added[0];
            justAdded = added[1];
        }
        // var baseNode = buildTree(dna.substring(curr, curr + l),k,ltClumpThresh,inclRevCompl);
        mfk = [];

        revComplAlreadyDone = [];

        kMersDone = {};
        mismatchSeqs = [];

        //console.log('st traverse');

        nodesVisited = 0;

        traverseTree(baseNode,baseNode,k,ltClumpThresh,inclRevCompl,maxMismatch);

        if (curr % 1000 == 1) {
            console.log('curr: ' + curr + ' nodes visited: ' + nodesVisited);
        }

        mfk.forEach(function(el) {
            var tmp = dna.substring(el[0][0], el[0][0] + k);
            el.push(tmp);


        });
        mfk  = cleanUpIdentical(mfk);

        mf = mfk;
    }
    else {

        if (curr == 0) {
           freqArray = buildFreqArray(k,dna.substring(0,l),inclRevCompl,maxMismatch,ltClumpThresh,false,progCallback);
           bestFreqArray = buildBestFreqArray(freqArray,progCallback);
            mf = [];

        }
        else {
            var oldestKmer = dna.substring(curr-1,curr-1 + k);
            var allOldDist = kMersWithMaxDist(oldestKmer,maxMismatch);
            allOldDist.forEach(function(el) {
                freqArray = deleteKmerFromFreqArray(el,curr-1,freqArray,inclRevCompl);
                });

            var newKmer = dna.substring(curr+l-k,curr+l);
            var allNewDist = kMersWithMaxDist(newKmer,maxMismatch);
            allNewDist.forEach(function(el) {
                var res = addKmerToFreqArray(el,curr+l-k,freqArray,bestFreqArray,inclRevCompl);
                freqArray = res[0];
                bestFreqArray = res[1];

            });

            mf = [];

        }

        //mf = mostFrequentKMersSort(k, dna.substring(curr, curr + l), inclRevCompl,maxMismatch, ltClumpThresh);

    }

    //console.log('curr: ' + curr + ' ' + mostFrequentKMers(k,dna.substring(curr,curr+l),document.getElementById('includeRevCompl').checked,ltClumpThresh));
    //console.log('curr: ' + curr + ' mf: ' + mf.length + ' mf innards: ' + mf[0] + ' tot mf: ' + mf);
    // console.log('st adjust');
    if (mf.length > 0) {

        if (useTree) {

        }
        else {
            /*
            var adjMf = mf.map(function (el) {
                adjEl0 = el[0].map(function (e) {
                    return e + curr;
                });
                adjEl1 = el[1].map(function (e) {
                    return e + curr;
                });
                return [adjEl0,adjEl1,el[2]];

            });
            //// console.log('curr: ' + curr + ' mf: ' + mf.length);
            mf = adjMf;
            */
        }
    }
    clumps.push(mf.length > 0);

    //console.log('st add to allmfk');

    if (useTree) {
        mf.forEach(function (el) {
            var alreadyThere = false;
            for (var i = 0; i < allMFK.length; ++i) {

              //  if (dna.substring(el[0][0], el[0][0] + k) === dna.substring(allMFK[i][0][0], allMFK[i][0][0] + k)) {
                if (el[2] === allMFK[i][2]) {
                    alreadyThere = true;
                    break;


                }

            }
            if (alreadyThere) {
                //No need to add. Note: this entry may actually have more frequency then prev one, so may need to tweak this
                if ((el[0].length + el[1].length) >= (allMFK[i][0].length + allMFK[i][1].length)) {
                    allMFK[i] = el;
                }
            }
            else {
                allMFK.push(el);
            }


        });
    }
    else {


    }



    if ((curr >= dna.length - l)) {  //|| stop) {
        mfk = [];

        if (useTree) {

        }
        else {
            allMFK = [];
            progThreshold = 10;
            bestFreqArray.forEach(function(el,i) {
                if (progCallback) {
                    if ((i % progThreshold == 0) || i == bestFreqArray.length - 1) {
                        progCallback('consol best',i,bestFreqArray.length - 1);
                    }

                }

                if (totIndArrayLen(el) >= ltClumpThresh) {
               // if ( (el[0].length + el[1].length) >= ltClumpThresh) {
                    var km = indToKmer(i,k);
                   if (!el[0]) {
                       el[0] = [];
                   }
                   if (!el[1]) {
                       el[1] = [];
                   }

                    allMFK.push([el[0],el[1],km]);


                }
            });
        }

        allMFK.forEach(function(el) {
            mfk.push(el);

        });
        //colourDNA(dna,clumps);

    }
    else {
        /*
        setTimeout(function() {
            timeOutLoop(dna,l,k,curr+1,useTree,inclRevCompl);
        },0);
        */

    }

    if (curr % 1000 == 20) {
        console.log('curr end: ' + curr);
    }

}

function delKmerFromTree(node,seq) {
    //deletes one occurence of kmer and its subkmers from tree (eg kmer ACGTG - deletes one A, AC, ACG, ACGT, ACGTG)
    var currNode = node;
    for (var i = 0;i < seq.length;++i) {
        for (var j = 0;j < currNode.subNodes.length;++j) {
            if (currNode.subNodes[j].letter === seq[i]) {
                if (currNode.subNodes[j].positions.length == 1) {
                    currNode.subNodes.splice(j,1);
                    i = seq.length;
                    break;
                }
                currNode = currNode.subNodes[j];
                currNode.positions.shift(); //assumes entry is in first position
                break;
            }
        }
    }

}

function addSingleBaseToTree(justAddedPrev,baseNode,letter,posn,k) {

    var justAddedCurr = [];

    justAddedCurr = [baseNode];

    justAddedPrev.forEach(function(ja) {

        if (ja.sequence.length > k) {

        }
        else {

            if (ja.subNodes.length == 0) {
                var newNode = new Node(letter);
                newNode.positions = [posn];
                newNode.sequence = ja.sequence + newNode.letter;
                ja.subNodes.push(newNode);
                justAddedCurr.push(newNode);

            }
            else {
                var snFound = false;
                var snInd = -1;


                ja.subNodes.forEach(function (sn, snNum) {
                    if (sn.letter === letter) {
                        snFound = true;
                        snInd = snNum;
                    }
                });
                if (snFound) {
                    var sn = ja.subNodes[snInd];
                    sn.positions.push(posn);
                    justAddedCurr.push(sn);
                }
                else {
                    var newNode = new Node(letter);
                    newNode.positions = [posn];
                    newNode.sequence = ja.sequence + newNode.letter;
                    ja.subNodes.push(newNode);
                    justAddedCurr.push(newNode);

                }


            }
        }


    });




    return [baseNode,justAddedCurr];


}


/*
function countProgMostFreqTree(stage,prog,tot) {
    self.postMessage({'task': 'mostFrequentKmers','msgType' : 'prog', 'stage' :stage, 'soFar' : prog, 'total' :tot});

}
*/


function countProgMostFreq(prog,tot) {
    self.postMessage({'task': 'mostFrequentKmers','msgType' : 'prog', 'soFar' : prog, 'total' :tot});

}

function countProgMostFreqWithStages(stage,prog,tot) {
    var now = Date.now();
    var diff = (now - startTime) / 1000;
    var mins = Math.floor(diff / 60);
    var secs = Math.round(diff -  (mins * 60));
    var formattedTime = mins + "'" + secs + '"';

    self.postMessage({'task': 'mostFrequentKmers','msgType' : 'prog', 'stage':stage,'soFar' : prog, 'total' :tot, 'elapsed':formattedTime});

}

function countProgltClump(stage,prog,tot) {
    var now = Date.now();
    var diff = (now - startTime) / 1000;
    var mins = Math.floor(diff / 60);
    var secs = Math.round(diff -  (mins * 60));
    var formattedTime = mins + "'" + secs + '"';

    self.postMessage({'task': 'ltClump','msgType' : 'prog', 'stage':stage,'soFar' : prog, 'total' :tot, 'elapsed':formattedTime});

}

function countProg(prog,tot) {
    self.postMessage({'task': 'searchKmer','msgType' : 'prog', 'soFar' : prog, 'total' :tot});

}

self.addEventListener('message', function(e) {

    console.log('Worker starting. Task: ' + e.data.task);
    importScripts('russGenLibs/BioUtilities.js');
    importScripts('russGenLibs/GenUtilities.js');

    startTime = Date.now();

    switch (e.data.task) {
        case 'searchKmer' :
            //console.log('dist: ' + dist('ACGCGGTCGAGCGCGG','GCGCGGCGCGAGCGCGG'));


            var dna = e.data.dna;
            var kmer = e.data.kmer;
            var inclRevCompl = e.data.inclRevCompl;
            var maxMismatch = e.data.maxMismatch;
            var indArray = countKMer(kmer,dna,inclRevCompl,maxMismatch,countProg);
            var txtArray = '';
            indArray.forEach(function(el) {
                txtArray+=' ' + el;
            });
            self.postMessage({'task': 'searchKmer','msgType' : 'result','txtStuff' : txtArray,'indArray' : indArray});
            break;
        case 'mostFrequentKmers':

            var dna = e.data.dna;
            var k = e.data.k;
            var inclRevCompl = e.data.inclRevCompl;
            var ltClumpThreshold = e.data.ltClumpThreshold;
            var maxMismatch = e.data.maxMismatch;
            var debug = e.data.debug;

            var freqCountArray = [];
            var debugInfo;

            startTime = Date.now();


            if (e.data.method === 'brute') {

                var res;


                switch(e.data.variant) {
                    case 1:
                        res = mostFrequentKMers(k, dna, inclRevCompl, maxMismatch,ltClumpThreshold,debug,countProgMostFreqWithStages);
                        break;
                    case 2:
                        res = mostFrequentKMersFreqArray(k, dna, inclRevCompl, maxMismatch,ltClumpThreshold,debug,countProgMostFreqWithStages);

                        break;
                    case 3:
                        res = mostFrequentKMersSort(k, dna, inclRevCompl, maxMismatch,ltClumpThreshold,debug,countProgMostFreqWithStages);
                        break;
                    case 4:
                        res = mostFrequentKMersObject(k, dna, inclRevCompl, maxMismatch,ltClumpThreshold,debug,countProgMostFreqWithStages);
                        break;

                    default:
                        break;
                }
                //mfk = mostFrequentKMers(k, dna, inclRevCompl, maxMismatch,ltClumpThreshold, countProgMostFreq);
                //var res = mostFrequentKMersSort(k, dna, inclRevCompl, maxMismatch,ltClumpThreshold,debug,countProgMostFreq);
                //var res = mostFrequentKMersFreqArray(k, dna, inclRevCompl, maxMismatch,ltClumpThreshold,debug,countProgMostFreqWithStages);
                if (debug) {
                    mfk = res[0];
                    debugInfo = res[1];
                }
                else {
                    mfk = res;
                    debugInfo = [];
                }

            }
            else {
                var added = buildTree(dna,k,ltClumpThreshold,inclRevCompl,countProgMostFreqWithStages);
                baseNode = added[0];
                //var justAdded = added[1];
                mfk = [];
                revComplAlreadyDone = [];
                kMersDone = {}; //trying object instead of array
                traverseTree(baseNode,baseNode,k,ltClumpThreshold,inclRevCompl,maxMismatch,countProgMostFreqWithStages);
                mfk.forEach(function(el) {
                    var tmp = dna.substring(el[0][0], el[0][0] + k);
                    el.push(tmp);


                });
                mfk  = cleanUpIdentical(mfk);
                debugInfo = [];


            }

            self.postMessage({'task': 'mostFrequentKmers','msgType' : 'result','indArray' : mfk, 'debugResult':debugInfo});//freqCountArray});


            break;
        case 'ltClump':
            var method = e.data.useTree ? 'tree' : 'brute';
            var dna = e.data.dna;
            var k = e.data.k;
            var ltClumpL = e.data.clumpL;
            if (ltClumpL > dna.length) { //clump  cannot be larger than DNA
                ltClumpL = dna.length;
            }

            var ltClumpT = e.data.clumpT;
            var inclRevCompl = e.data.inclRevCompl;
            var maxMismatch = e.data.maxMismatch;
            mfk = [];
            allMFK = [];
            clumps = [];
            revComplAlreadyDone = [];
            kMersDone = {}; //trying object instead of array

            freqArray = [];
            bestFreqArray = [];

            startTime = Date.now();

            for (var curr = 0;curr <= dna.length - ltClumpL;++curr) {
                ltClump(dna, ltClumpL,ltClumpT, k, curr, e.data.useTree, inclRevCompl,maxMismatch,countProgltClump);
            }
            if (method  === 'brute') {
               console.log('worker commencing ltClump brute: ' + e.data);
            }
            else {
                console.log('worker commencing ltClump tree: ' + e.data);
            }

            self.postMessage({'task': 'ltClump','msgType' : 'result','indArray' : mfk, 'clumps' :clumps});



            // self.postMessage({'task': 'ltClump','msgType' : 'result','indArray' : mfk});


            break;
        case 'motifSearch' :
            //console.log('dist: ' + dist('ACGCGGTCGAGCGCGG','GCGCGGCGCGAGCGCGG'));

            var dnaStrings = e.data.dnaStrings;
            dnaStrings = dnaStrings.filter(function(el) {
                if (el.length == 0) {
                    return false;
                }
                else {
                    return true;
                }
            });

            startTime = Date.now();

            var k = e.data.k;
            var inclRevCompl = e.data.inclRevCompl;
            var maxMismatch = e.data.maxMismatch;
            var indArray = motifSearch(dnaStrings,k,maxMismatch,inclRevCompl,countProgMostFreqWithStages);
           /* var indArray = countKMer(kmer,dna,inclRevCompl,maxMismatch,countProg);
            var txtArray = '';
            indArray.forEach(function(el) {
                txtArray+=' ' + el;
            });
            */
            self.postMessage({'task': 'motifSearch','msgType' : 'result','txtStuff' : txtArray,'indArray' : indArray});
            break;

        case 'bruteMotifSearch' :

            var res = bruteMotifSearch(e.data.dnaStrings, e.data.k, e.data.maxMismatch, e.data.laplace,countProgMostFreqWithStages);
            self.postMessage({'task': 'bruteMotifSearch','msgType' : 'result','txtStuff' : 'Score: ' + res[1],'indArray' : res[0]});
            break;


        case 'medianString' :
            //console.log('dist: ' + dist('ACGCGGTCGAGCGCGG','GCGCGGCGCGAGCGCGG'));

            var dnaStrings = e.data.dnaStrings;
            dnaStrings = dnaStrings.filter(function(el) {
                if (el.length == 0) {
                    return false;
                }
                else {
                    return true;
                }
            });

            startTime = Date.now();

            var k = e.data.k;
            var inclRevCompl = e.data.inclRevCompl;
            var maxMismatch = e.data.maxMismatch;
            var res = medianString(dnaStrings,k,maxMismatch,inclRevCompl,countProgMostFreqWithStages);
            var indArray = res[0];
            /* var indArray = countKMer(kmer,dna,inclRevCompl,maxMismatch,countProg);
             var txtArray = '';
             indArray.forEach(function(el) {
             txtArray+=' ' + el;
             });
             */
            var txtArray = 'Score: ' + res[1];
            self.postMessage({'task': 'medianString','msgType' : 'result','txtStuff' : txtArray,'indArray' : indArray});
            break;

        case 'greedyMotif' :
            //console.log('dist: ' + dist('ACGCGGTCGAGCGCGG','GCGCGGCGCGAGCGCGG'));

            var dnaStrings = e.data.dnaStrings;
            dnaStrings = dnaStrings.filter(function(el) {
                if (el.length == 0) {
                    return false;
                }
                else {
                    return true;
                }
            });
            dnaStrings = dnaStrings.map(function(el) {
                return el.trim();

            });

            startTime = Date.now();

            var k = e.data.k;
            var laplace = e.data.laplace;
            var maxMismatch = e.data.maxMismatch;
            var res = greedyMotif(dnaStrings,k,maxMismatch,laplace,countProgMostFreqWithStages);
            var indArray = res[0];
            var txtArray = 'Score: ' + res[1];
            /* var indArray = countKMer(kmer,dna,inclRevCompl,maxMismatch,countProg);
             var txtArray = '';
             indArray.forEach(function(el) {
             txtArray+=' ' + el;
             });
             */
            self.postMessage({'task': 'greedyMotif','msgType' : 'result','txtStuff' : txtArray,'indArray' : indArray});
            break;

        case 'randomMotif' :

            var res = randomMotif(e.data.dnaStrings, e.data.k, e.data.numIters, e.data.maxMismatch, e.data.laplace,countProgMostFreqWithStages);
            self.postMessage({'task': 'randomMotif','msgType' : 'result','txtStuff' : 'Score: ' + res[1] + '\nEntropy: ' + res[2],'indArray' : res[0]});
            break;

        case 'gibbsSampler' :

            var res = gibbsSampler(e.data.dnaStrings, e.data.k, e.data.numIters, e.data.maxMismatch, e.data.laplace,countProgMostFreqWithStages);
            self.postMessage({'task': 'gibbsSampler','msgType' : 'result','txtStuff' : 'Score: ' + res[1],'indArray' : res[0]});
            break;



        default:
            break;
    }

}, false);

/*
//Bioinformatics utilities

function indToKmer(ind,k) {

    var curr = ind;
    var pow;

    var inds = '';

    for (var i = 0;i < k;++i) {
        pow = k - i - 1;
        var sum = Math.floor(curr / Math.pow(c_NumBases,pow));
        var rem = curr % Math.pow(c_NumBases,pow);
        inds+= c_Bases[sum];
        curr = rem;

    }

    return inds;



}

function kMerToInd(kmer) {

    var pow = 0;
    var sum = 0;
    var vals = {'A':0,'C':1,'G':2,'T':3};


    for (var i = kmer.length-1;i >= 0; --i) {

        sum+= vals[kmer[i]] * (Math.pow(c_NumBases,pow));
        ++pow;

    }

    return sum;


}

*/

Array.max = function( array ){
    return Math.max.apply( Math, array );
};

Array.min = function( array ){
    return Math.min.apply( Math, array );
};