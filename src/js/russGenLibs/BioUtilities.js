/**
 * Created by RussellM on 20/01/2016.
 */


//Bioinformatics routines

var c_Bases = ['A','C','G','T'];
var c_NumBases = 4;
var c_BaseInds = {'A':0,'C':1,'G':2,'T':3};


function hamDist(dna1,dna2) {
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

    var allDist = [];
    for (var ii = 0;ii <= maxDist;++ii) {
        var ret = kMersWithDist(kmer, ii);
        allDist = allDist.concat(ret);
    }

    return allDist;


}

function kMersWithDist(kmer,dist) {

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

function gcSkew(dna,numBases,skewAlready) {

    //skewAlready contains skews already processed for part of DNA

    var skewArray = [0];
    var min = 0;
    var max = 0;
    var minPos = -1;
    var maxPos = -1;
    var minPosArray = [];
    var maxPosArray = [];
    var curr = 0;
    var startPos = 0;

    if (skewAlready) {
        skewArray = skewAlready[2];
        min = skewAlready[0];
        max = skewAlready[1];
        minPosArray = skewAlready[3];
        maxPosArray = skewAlready[4];
        curr = skewArray[skewArray.length - 1];
        startPos = skewArray.length - 1;

    }

    var skewCounts = {
        'A': 0,
        'C': -1,
        'G': 1,
        'T': 0,
        'a': 0,
        'c': -1,
        'g': 1,
        't': 0

    };
    // var skewArray = [0];
    //  var curr = 0;
    for (var j = startPos;j < startPos + numBases;++j) {
        if (j >= dna.length) {
            break;
        }
        curr+= skewCounts[dna[j]];
        skewArray.push(curr);
        if (curr < min) {
            min = curr;
            minPosArray = [j+1];
        }
        else if (curr == min) {
            minPosArray.push(j+1);
        }

        if (curr > max) {
            max = curr;
            maxPosArray = [j+1];
        }
        else if (curr == max) {
            maxPosArray.push(j+1);

        }
    }
    return [min,max,skewArray,minPosArray,maxPosArray];
}



function gcCount(dna) {


    var bCount = baseCount(dna);

    return bCount[1] + bCount[2];
}

function baseCount(dna) {

    var baseCount = [0,0,0,0];

    for (var i = 0; i < dna.length;++i) {
        switch (dna[i].toUpperCase()) {
            case 'A':
                ++baseCount[0];
                break;

            case 'C':
                ++baseCount[1];
                break;
            case 'G':
                ++baseCount[2];
                break;

            case 'T':
                ++baseCount[3];
                break;

            default:
                break;

        }

    }
    return baseCount;

}

function reverseComplement(dna) {

    if (!dna) {
        return '';
    }

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

function randomDNA(len) {

    var randDNA = '';

    for (var i = 0;i < len;++i) {
        var r = getRandomInt(0, c_NumBases - 1);
        randDNA += c_Bases[r];
    }


    return randDNA;


}

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

function motifsToColumns(motifs) {
    var len = motifs[0].length;

    motifs.forEach(function(el) {
        if (el.length == len) {

        }
        else  {
            return -1;
        }
    });

    var colArray = [];
    for (var i = 0;i < len;++i) {
        var colStr = '';
        motifs.forEach(function(el) {
            colStr += el.substring(i,i+1);
        });
        colArray.push(colStr);

    }

    return colArray;
}

function scoreMotifs(motifs,laplace) {
    //assumes equal length motifs, otherwise returns -1
    // laplace flag determines whether to apply rule of succession (ie inc each count by 1, to avoid zeros)

    var cols = motifsToColumns(motifs);

    var score = 0;
    cols.forEach(function(col,c) {
        var counts = baseCount(col);
        var m = counts.reduce(function(a,b){
            return Math.max(a,b);

        });

        score += col.length - m;
        /*
        if (laplace) {
            score += (col.length - m) + c_NumBases - 1;
        }
        else {
            score += col.length - m;
        }
        */



    });

    var countMatrix = [];

    c_Bases.forEach(function(b) {
        var entry = [];
        for (var c = 0;c < motifs[0].length;++c) {
            entry.push(0);
        }
        countMatrix.push(entry);
    });

    cols.forEach(function(col,c) {

        var counts = baseCount(col);
        counts.forEach(function(baseCount,r) {
            countMatrix[r][c] = laplace ? baseCount + 1 : baseCount;
        });


    });

    var profileMatrix = countMatrix.map(function(el,r) {
        return el.map(function(col,c) {
            if (laplace) {
                return col * 1.0 / (motifs.length + c_NumBases);
            }
            else {
                return col * 1.0 / motifs.length;
            }

        });

    });

    var entropyMatrix  = profileMatrix.map(function(el,r) {
        return el.map(function(col,c) {
            return (col == 0) ? col : -1.0 * col * Math.log2(col);

        });

    });




    var consensus = '';
    var entropyScore = 0;

    for (var c = 0;c < motifs[0].length;++c) {
        var m = -1;
        var mInd = -1;
        for (var r = 0;r < c_NumBases;++r) {
            if (countMatrix[r][c]  > m) {
                m = countMatrix[r][c];
                mInd = r;
            }
            entropyScore+= entropyMatrix[r][c];


        }
        consensus+= c_Bases[mInd];


    }
   // entropyScore *= -1;

    return [consensus,score,entropyScore,countMatrix,entropyMatrix,profileMatrix];
}

function calcMotifLogo(motifs,laplace) {
    var res = scoreMotifs(motifs,laplace);

    var profileMatrix = res[5];
    var entropyMatrix = res[4];
    var consensus = res[0];

    var colEnts = [];

    for(var c = 0;c < consensus.length;++c) {
        var colEnt = 0;
        for (var r = 0;r < entropyMatrix.length; ++r) {
            colEnt+=entropyMatrix[r][c];

        }
        colEnts.push(colEnt);

    }

    var colHeightFactors = colEnts.map(function(colEnt) {
        return  1.0 - colEnt / 2.0;

    });
    var charHeightMatrix = profileMatrix.map(function(rowEl,r) {
        return rowEl.map(function(colEl,c) {
            return [c_Bases[r],colEl * colHeightFactors[c]];
        });

    });

    var colMatrix = [];
    for (var c = 0;c < consensus.length;++c) {
        var colRows = [];
        for (var r = 0; r < profileMatrix.length;++r) {
            colRows.push(charHeightMatrix[r][c]);
       }
        colRows.sort(function(a,b) {
            if (a[1] == b[1]) {
                return 0;
            }
            else if (a[1] < b[1]) {
                return 1;
            }
            else {
                return -1;
            }

        });
        colMatrix.push(colRows);
    }

    var tr = transpose(colMatrix);

    return tr;

}


function profileMostProbable(profMatrix,dna,k) {
    //returns most probable kmer in dna based on profile matrix,its probability, and its position

    var bestProb = -1;
    var bestInd = -1;

    for (var i = 0;i < dna.length - k +1 ; ++i) {
        var prob = 1.0;
        var kmer = dna.substring(i,i+k);
        for (var c = 0;c < kmer.length;++c) {
            var r = c_BaseInds[kmer[c]];
            prob*= profMatrix[r][c];
        }
        if (prob > bestProb) {
            bestInd = i;
            bestProb = prob;

        }

    }
    return [dna.substring(bestInd,bestInd+k),bestProb,bestInd];

}

function profileProbDist(profMatrix,dna,k) {
    //selects a random kmer from dna with probability distribution based on profile matrix
    //returns the kmer,its probability, and its position



    var probDist = [];

    for (var i = 0;i < dna.length - k +1 ; ++i) {
        var prob = 1.0;
        var kmer = dna.substring(i,i+k);
        for (var c = 0;c < kmer.length;++c) {
            var r = c_BaseInds[kmer[c]];
            prob*= profMatrix[r][c];
        }
        probDist.push(prob);
   }
   var r = getRandomFromProbDist(probDist);

    return [dna.substring(r,r+k),probDist[r],r];

}


function motifsMostProbable(profMatrix,dnaStrings,k) {
    return dnaStrings.map(function(el) {
        return profileMostProbable(profMatrix,el,k);
    });
}

function patternMotifsDist(pattern,motifs) {


    var motifDists = motifs.map(function(el) {
        return hamDist(pattern,el);
     });
    return motifDists.reduce(function(el1,el2) {
        return el1 + el2;
    });
}

function patternDNADist(pattern,dna) {
    var k = pattern.length;

    var min = 999999;
    var minInd = -1;

    for (var i = 0;i < dna.length - k + 1;++i) {
        var hd = hamDist(pattern,dna.substring(i,i+k));
        if (hd  < min) {
            minInd = i;
            min = hd;
       }

    }
    return [min,minInd];

}

function patternSequencesDist(pattern,sequences,returnPosns) {

    var sequencesUpper = sequences.map(function(el) {
        return el.toUpperCase();

    });
    var patternUpper = pattern.toUpperCase();

    var d = 0;
    var dPosns = [];
    sequencesUpper.forEach(function(dna) {
        var res = patternDNADist(patternUpper,dna);
        d+= res[0];
        if (returnPosns) {
            dPosns.push(res[1]);
        }

    });
    /*
    var seqDists = sequencesUpper.map(function(dna) {
        return patternDNADist(patternUpper,dna)[0];

    });
    return seqDists.reduce(function(d1,d2) {
        return d1 + d2;
    });
    */
    if (returnPosns) {
        return [d, dPosns];
    }
    else {
        return d;
    }

}

function motif(pattern,dna) {

    var res = patternDNADist(pattern,dna);
    if (res[1] > -1 ) {
        return dna.substring(res[1],res[1] + pattern.length);
    }
    else {
        return '';
    }


}