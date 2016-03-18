
//Global Variables

var k = 2;
var mfk = [];
var mfMotif = []; //used for motif searches

var sk;

var dnaMaster = '';
var dnaPage = 0;
var dnaPageOffset = 0;

var dnaMasterStrings = []; // used for Motif search

var basesPerPage = 20000;

var clumps = [];
var stop = false;
var myParams = Params.getInstance();

var w = new Worker('js/worker.js');

//Initialisation

initialisePage();
testStuff();


function testStuff() {

    //var d = patternMotifsDist('ACG',['ACT','GGG','TCG']);

    //var t = motif('GATTCTCA', 'GCAAAGACGCTGACCAA');


    //var allK = allKmers(0);

    //var f = kMerToInd('GCTCCTTGGGCGCAGATCT');
    //var g = indToKmer(8026,11);

    //var hd = hamDist('TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC','GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA');
    //document.getElementById("debugText").innerHTML = 'ham: ' + hd;

    //var neighbours = kMersWithMaxDist('ACGT',3);
    //var output = '';
    //neighbours.forEach(function(el) {
    //    output+=el + '\n';
    //});

   // var testMotifs = ['TCGGGGGTTTTT','CCGGTGACTTAC','ACGGGGATTTTC','TTGGGGACTTTT','AAGGGGACTTCC','TTGGGGACTTCC','TCGGGGATTCAT','TCGGGGATTCCT','TAGGGGAACTAC','TCGGGTATAACC'];
   // logoCanvas(testMotifs);


    var rands = [];
    for (var i = 0;i < 100; ++i) {
        rands.push(getRandomFromProbDist([0.05,0.05,0.05,0.7,0.05,0.05]));
    }


    var spaceString = 'ACATTGGAATTATTGAAGTCGGCAAAAGAGCCAGAAGGGGTTTAATATGGGCAGGACTACATGGGAGACGTCGAATGTACAAGTGCGGCTCGCTTACACACTG GGATACGGAATGAGTATAATTGATTTGAGCCGGGCACGGCTCAGGAGCTACAGAGGCTACTTGGTAAGGGGTCTGGTAACTTGTGGAATCGCCCCTGTGGCAG GAATTAGCGAGATGCTCTCCTAAGTTCGTTCCTTAGATTCTGTAATAAGGACTACGAATATGAAGCCCTCGTAGCGATAATACGATCCGCTGATAAGGCACAT GGTCAGATCTCAGGCGTGTTGGCGTTTAGCTGTGGCACCAGCTTTAACTGCGCGAGCATCCACTTGATGGCTGTTACGAATAATTGCGAGCTCATAAACATCA TTCGTTACATTTTCCCGCCAGTCTGGCACTCCAGCCGGAGTGGTTCACTCACTTTTCGGCTGATTTTTGCCCAATACGGGGGAAGGCCCTCTAAAGAACAACT GTTGCAGTCGGGAAATTAGAACGTATTCTTGGTTGACTATCCTAGCGCGGAGTATCACCCAGCGGACTGCCAACTCCTTCCCTGGTATGTCTTTGACTTCTTA GCACCTTGACTGGAACCCGCGGTAACAGGAATAAGAATTTTGAACCAGACGGGAGAAAGGGTCCAGCACTTGTCCATGCATGGGGAGCAAAGGTTAGTTTACT GTAACGACCCCGAGGAGTCAACCGCAGATTACGGGCCCGAGCGGGGATGTTGCCGGGTCGGTGTTGGGGCATTCGCACCTGGTGGGGGTTGGAAGTACCACCA CGCACCGAGTCGTACTATCGTTCTTGACAGACTTGGCTCGACACCTGGACCGTACATTTGCGGCCAATGACCAAGCCTGTAGCACACCTCAAACAGCCGCTCA GACTGGGTGGATGTTCAGTAATTACTCGTATCATGGGCGCTTCGTTGATCCGAAGCCCGAAACGACTCTACGAACCTGTCCAAGAGGGTCCCGGTCGAAACGT GTGTCAAAAGAGTAATTAGTAAGGCAAACTACGCAGGTCCCACGTAAAGTGGCTCTAGCGTCGACGGCCAGAACGGTAGCGACCTTAATTTGGATAAAGAACT AGCCAGTCATGTGCATCCTCGTGAGGATTCTCTCTGCAGTGCGTTGCGAACATTGCGCAAGGCGAGATTCTCCATCCTAAGAGTGTGGTTTATTTCTGCAAGT CTGAGTAATCGGCTTAGGTATTATGTCGGCGTAAGTCCCCTGACTCCAGCCGCCAGAAGTGTCTAACTTGAAGGTTACTCGCCACTGACGGTGATTTGCCGTT TACGCTATAATGATATAGTCACGGCACGAGAGCAAGTACTTTGCGCACAATTGGGCGTATCCCAGTTGATCTACAACCACTGGAACAACTTGCCGACTATCGG AACATAGAAAAGCGTAACCCGCGACGATTGCTGTTATTGAATCAAGTAAAATTCAGATCAGTAGAGCGATTTGCGGTGCATCTTTCAGGACACGTTATTCGGC ATCCCACCTATCCCTCCGTCATGTAATTACAATGAGGAGGCGTATCAATTGTTACAGAAGCTTCCACATATAGTTGTTGTGGTATCTTTTCGTTAGAAAGCCT GTCGAAGATACCAGTAATCATCGATCAGCGATTAGTATGATGTGGAGGGCGTGCATCTCCTCGGAGACACGTTCTATAGACAATGAGCGGAATAGTCAGTTAA GGGCGGAACCCGGAACTAAATTGAGCGTCCGAGGACCTGGTCTTTACAGATGGAGGGGAAGTGGGACTAGTTGATTTTGAGAGGAAACCATGGGAGTTTGGGG CGCGTTCAAATATGCCTCATCGGATCGCGTTTCCCTACAACCACCAGATCAAAACTGCCCAAAGATCAAACGCTTACAGCCTGGCAGCCTTCATTCCTTGATG ATCCCATCTAAGACGCCCGCAATTTCCAGCGACTTGTGTCAGGCTGTACTTTAATGCCATATATGGTGGACAGCCAGCGCACTCATCGGTGATAATCCAAAAA AGCAGAGGCAGAACTGTTGGTCCTTCTATGGGCATACGGGTGTAAGTTTCTCCCGTGGAGCACCTATCCGTCCTACATAACAGGCGAGGCACTAGGGTCTGCT GGTATTCATGTCCTCGCACACATGTAGTACTATGCGATTTTGCAAAATCAAGAGCATTATGCCCGTCTACTCCGTTATCATTTTAGAAACCCTCTGATGTAGG ACTCATTTAACGCCCGGGTCATGAGGAAATTACCTCTACTAGCACTCATCTTGCGAACGACACCTCATAACGATTGAGTTTCTTAAGCCGTCAATAGTCTTTA ATCAAAGGTTCAGTGAAGCTTCTGCACTGAAGCGCACATTTCTTTGCAGGTGACGTTCAAGCCCGCTATGGAATTGATACGGAGATTTCCCCGCCGGTGATGG AGGAAATCAGTGCTTTTCCTCCAATCCCTTGGCCCTTTTGTCGGCTCCGTAAATAACCAAGGTGCAGTTCGAATTCGATGAATAGCCATGCTGACGGACCTTC TGTCGCTGGGCGCAAGGGCACCACATGTTTCTACGTGCTTGGGTGCTAACTTGGTGAATTGGTTAAGGGTGCGCGGACGCTCGTCACCAACCGGACAACAGTA AGAGCTACGACAGCTTAGGCTCCGGAGGATTCGCTTTCCGGTAAAGATCACCTCTAACCGATGCAACGCCCTCCGGTCAAACACGAGTCCAAAATTTTACACC TGGGTTGAGATCCTCAGAGTCCAATGCGCGGGCATCCTCACAGAACAATATTTCTGAGTTATCTATAACGTCCTTCTCATGTTTCAGCCTTGGCCTAATGTTG CTTTTGTAGGGCATGTATAATGCTGACAATGACGAAAGCTTTTTCGGCCCCATTTGAACTCCGGCTTTGCCTTCATTAATTAATTTAGGCCAGCAGTCATTCC AGGCTTTGGAGGCAAGCTGTTCGATAGGGTGGTATTTGTAGGATACGACAAAGCAGGTGGAGACTGATTCCATGAAGCAGCCTGGCTCCAGCCAGACACAGAC GCACGTGCGCAGTAATATGTATCTTCCCTTTGGGTCGCTAGCTGAGCGTCTACTCCGACGGTTGCCCTCTTCCATGTCTAACTCGGAGGCGGATCCTGCAGGA AACGAATTTAGAACGAGCGTTAAGCGGCTCAATAAAATTGAGGGCCAAAACCGACTAGGTCAAATGTGCTGAAGTAGCACATAAACTAAAACATGCGGGGGAA GACGCTAAAACTCGCTTATAGTTGGCTTTACATTGGTCTCCAAAAAAATCGTTCGACTCAGCGCGAGCACCTTTTCACCTGACTAGACGTGTGTAGATTATAA CTATGACACGCCAAGTCGTAAGTGTGCATGAAGCCCATTGGACTTCACGTGCTAGATGGCAAAATCGCAGCAGCCGGTAGCTTGCTCACACCACCCGAATAAT';
    var arr = spaceString.split(' ');
    var filtered = arr.filter(function(el) {
        if (el.length == 0) {
            return false;
        }
        else {
            return true;
        }
    });


   // var cyc = cycloPeptideSequencing('0 113 128 186 241 299 314 427');
    var cyc = cycloPeptideSequencing('0 113 114 128 129 242 242 257 370 371 484');


    //var all = allPossiblePeptidesWithWeight(1024);

    var dna = new DNA('AGTUTGGAGCTGATGTATGGCCAAAACGCTTAGCTAGCCGTG');
    //alert(dna.rnaTranscript());

    var rna = new RNA(dna.rnaTranscript());

    var bob = rna.translate(0);

    var str = bob.toMedString('');

    var str2 = bob.toLongString();

    var w = bob.getIntegerWeight();

    var pep = new Peptide(Peptide.AminoArrFromStr('AVGQMM'));

    var wStr = pep.toWeightString();

    var tryro = new Peptide(Peptide.AminoArrFromStr('VKLFPWFNQY'));

    var tst = new Peptide(Peptide.AminoArrFromStr('NQEL'));

    var tstSpec = tst.linearSpectrum();

    var d2 = patternSequencesDist('CTCTGG',filtered);

    d2 = patternSequencesDist('aaa',
        [
            '',
            '',
            '',
            '',
            '',
            '',
            '',
            '',
            ''
        ],true);

    d2 = patternSequencesDist('aaa',
        [
            '',
            '',
            '',
            '',
            '',
            '',
            '',
            '',
            ''
        ],true);


}

function initialisePage() {

    myParams.tabActive = 5;
    tab_click(4); //tab stuff
    tab_click(3);
    tab_click(2);
    tab_click(1);
    tab_click(0);

    processHTMLParams();
    setUpWorkerListeners();

    document.getElementById('fileInput')
        .addEventListener('change', readSingleFile, false);

    document.getElementById('fileMotifInput')
        .addEventListener('change', readSingleSequencesFile, false);

//  document.getElementById('motifBrute').checked = true;
    motifRadClicked('motifBrute');
    sequencingRadClicked('seqPath');
    sequencingInputRadClicked('seqDNA');

}


function setUpWorkerListeners() {
    w.addEventListener('message', function(e) {

        var resString, numKmer, numKmerVal, tot, dna, km;

        if (e.data.msgType === 'result') {
            console.log('Worker finished. Task: ' + e.data.task);
        }

        if (e.data.msgType === 'prog') {
            outputProgress(e);
        }
        else {

            switch (e.data.task) {
                case 'searchKmer':
                     // console.log('Worker finished');
                    // alert('kmer ind is: ' + e.data.txtStuff);
                    // mfk = [indexArray];
                    mfk = [e.data.indArray];
                    //k = kmer.length;
                    document.getElementById('numKmers').value = 1;
                    numKmerChanged();
                    //colourDNA(dnaMaster,null,document.getElementById("includeRevComplKS").checked);

                    if (document.getElementById('debugKS').checked) {
                        resString = '';

                        var combined = mfk[0][0].concat(mfk[0][1]);
                        combined.sort();
                        combined.forEach(function (el) {
                            resString += el + ' ';
                        });

                        /*
                         mfk[0][0].forEach(function(el) {
                         resString+=el + ' ';
                         });
                         if (mfk[0][1].length == 0) [

                         ]
                         else {
                         resString+='\n';
                         mfk[0][1].forEach(function(el) {
                         resString+=el + ' ';
                         });

                         }
                         */

                        document.getElementById('debugText').value  += '\n' + resString;
                        //alert(debugRes);
                    }

                    break;

                case 'mostFrequentKmers':

                    mfk = e.data.indArray;
                    tot = 0;
                    if (mfk.length == 0) {

                    }
                    else {
                        tot = mfk[0][0].length + mfk[0][1].length;
                    }
                    // mfk[0].forEach(function(el) {
                    //    tot+= el.length;

                    //});

                    document.getElementById('kMerOccurrences').innerHTML = tot;//mfk[0].length;
                    document.getElementById('kMerNumEqualMostFreq').innerHTML = mfk.length;

                    ///console.log('mfk: ' + mfk + ' len: ' + mfk.length);
                    /*
                     mfk.forEach(function (m) {
                     console.log('mfk entry: ' + m);
                     });
                     */

                    numKmer = document.getElementById('numKmers');
                    numKmer.max = mfk.length;
                    if (mfk.length == 0) {
                        numKmer.min = 0;
                        numKmer.value = 0;
                    }
                    else {
                        numKmer.min = 1;
                        numKmer.value = 1;
                    }


                    numKmerVal = parseInt(numKmer.value);
                    dna = dnaMaster; // document.getElementById('dnaInput').value;

                    if (mfk.length == 0) {
                        document.getElementById('kMerMostFreq').innerHTML = 'None';
                        document.getElementById('kMerRevCompl').innerHTML = 'None';

                    }
                    else {
                        km = dna.substring(mfk[numKmerVal - 1][0][0], mfk[numKmerVal - 1][0][0] + k);
                        // var km = dna.substring(mfk[numKmerVal - 1][0], mfk[numKmerVal - 1][0] + k);
                        document.getElementById('kMerMostFreq').innerHTML = dna.substring(mfk[numKmerVal - 1][0][0], mfk[numKmerVal - 1][0][0] + k);
                        document.getElementById('kMerRevCompl').innerHTML = reverseComplement(km);
                    }


                    revComplAlreadyDone = [];
                    colourDNA(dna, null, document.getElementById('includeRevComplMF').checked);

                    if (document.getElementById('debugMF').checked) {
                        resString = '';
                        var debugRes = e.data.debugResult;
                         debugRes.forEach(function (el) {
                            resString += el;
                            resString += ' ';
                        });
                        if (resString.length == 0) {
                            document.getElementById('debugText').value += '\nNone';
                        }
                        else {
                            document.getElementById('debugText').value += '\n' + resString;
                        }
                        //alert(debugRes);
                    }

                  break;

                case 'ltClump':

                    mfk = e.data.indArray;
                    var clumps = e.data.clumps;
                    tot = 0;
                    if (mfk.length == 0) {

                    }
                    else {
                        tot = mfk[0][0].length + mfk[0][1].length;
                    }
                    // mfk[0].forEach(function(el) {
                    //    tot+= el.length;

                    //});

                    document.getElementById('kMerOccurrences').innerHTML = tot;//mfk[0].length;
                    document.getElementById('kMerNumEqualMostFreq').innerHTML = mfk.length;

                    console.log('mfk: ' + mfk + ' len: ' + mfk.length);
                    mfk.forEach(function (m) {
                        console.log('mfk entry: ' + m);
                    });

                    numKmer = document.getElementById('numKmers');
                    numKmer.max = mfk.length;
                    if (mfk.length == 0) {
                        numKmer.min = 0;
                        numKmer.value = 0;
                    }
                    else {
                        numKmer.min = 1;
                        numKmer.value = 1;
                    }


                    numKmerVal = parseInt(numKmer.value);
                    dna = dnaMaster; //document.getElementById('dnaInput').value;

                    if (mfk.length == 0) {
                        document.getElementById('kMerMostFreq').innerHTML = 'None';
                        document.getElementById('kMerRevCompl').innerHTML = 'None';


                    }
                    else {
                        km = dna.substring(mfk[numKmerVal - 1][0][0], mfk[numKmerVal - 1][0][0] + k);
                        // var km = dna.substring(mfk[numKmerVal - 1][0], mfk[numKmerVal - 1][0] + k);
                        document.getElementById('kMerMostFreq').innerHTML = dna.substring(mfk[numKmerVal - 1][0][0], mfk[numKmerVal - 1][0][0] + k);
                        document.getElementById('kMerRevCompl').innerHTML = reverseComplement(km);


                    }


                    revComplAlreadyDone = [];
                    colourDNA(dna, clumps, document.getElementById('includeRevComplLT').checked);

                    break;

                case 'motifSearch':

                     processReturnedMotif(e);
                     break;

                case 'medianString':

                     processReturnedMotif(e);
                     break;

                case 'greedyMotif':

                    processReturnedMotif(e);
                    break;

                case 'randomMotif':

                    processReturnedMotif(e);
                    break;

                case 'gibbsSampler':

                    processReturnedMotif(e);
                    break;

                case 'bruteMotifSearch':

                    processReturnedMotif(e);
                    break;

                default:
                     break;
            }
        }
        }, false);



}


function outputProgress(e) {
    //console.log('most freq kmer prog: ' + e.data.soFar);
    document.getElementById('progress').innerHTML =    e.data.soFar + ' / ' + e.data.total + ' ' + e.data.elapsed;
    if (e.data.stage) {
        document.getElementById('progress').innerHTML = e.data.stage + ' ' + document.getElementById('progress').innerHTML ;
    }

}

function processReturnedMotif(e) {
    mfMotif = e.data.indArray;
    var tot = 0;
    if (mfMotif.length == 0) {

    }
    else {
       // tot = mfMotif[0][0][0].length + mfMotif[0][0][1].length;
        tot = mfMotif[0].length;
    }


    document.getElementById('kMerOccurrences').innerHTML = tot;//mfk[0].length;
    document.getElementById('kMerNumEqualMostFreq').innerHTML = mfMotif.length;


    var numKmer = document.getElementById('numKmers');
    numKmer.max = mfMotif.length;
    if (mfMotif.length == 0) {
        numKmer.min = 0;
        numKmer.value = 0;
    }
    else {
        numKmer.min = 1;
        numKmer.value = 1;
    }


    var numKmerVal = parseInt(numKmer.value);
    var dna = dnaMaster; // document.getElementById('dnaInput').value;

    if (mfMotif.length == 0) {
        document.getElementById('kMerMostFreq').innerHTML = 'None';
        document.getElementById('kMerRevCompl').innerHTML = 'None';

    }
    else {
        //var km = dna.substring(mfk[numKmerVal - 1][0][0], mfk[numKmerVal - 1][0][0] + k);
        var km = mfMotif[0][0][2];
        // var km = dna.substring(mfk[numKmerVal - 1][0], mfk[numKmerVal - 1][0] + k);
        document.getElementById('kMerMostFreq').innerHTML =  km ;//dna.substring(mfk[numKmerVal - 1][0][0], mfk[numKmerVal - 1][0][0] + k);
        document.getElementById('kMerRevCompl').innerHTML = reverseComplement(km);
    }


    revComplAlreadyDone = [];
    colourDNA(dna,null,false);

    if (document.getElementById('debugMS').checked) {

        var txtStuff = e.data.txtStuff;
        document.getElementById('debugText').value = txtStuff;

        var extraInfoAr = [];

        mfMotif[0][0].forEach(function(el,i) {
            var ind = el[0][0];
            var km = dnaMasterStrings[i].substring(ind,ind+k);

            extraInfoAr.push(km);
        });

        logoCanvas(extraInfoAr);
        /*
         mfMotif.forEach(function(el) {
         extraInfoAr.push(el[1]);

         });
         extraInfoAr.sort();
         */
        var extraInfo = '';
        extraInfoAr.forEach(function(el) {
            extraInfo+=el + '\n';
        });
        document.getElementById('moreDetailText').value = extraInfo;
        //alert(debugRes);
    }

}


function readSingleFile(e) {
    var file = e.target.files[0];
    if (!file) {
        return;
    }
    var reader = new FileReader();
    reader.onload = function(e) {
        var contents = e.target.result;

        dnaMaster = cleanContents(contents);
        dnaMasterChanged();
        document.getElementById('dnaInput').value = '';
    };
    reader.readAsText(file);
}

function readSingleSequencesFile(e) {
    var file = e.target.files[0];
    if (!file) {
        return;
    }
    var reader = new FileReader();
    reader.onload = function(e) {
        var contents = e.target.result;


        var parts = contents.split(/\r\n|\r|\n/g);
        parts = parts.filter(function(e) {
            if (e.substring(0,1) === '>') {
                return false;
            }
            else {
                return true;
            }
        });

        var joined = parts.join('\n');


        document.getElementById('dnaStrings').value = joined;
        motifsInput();
    };
    reader.readAsText(file);
    //e.target.files[0] = '';
   // document.getElementById('fileMotifInput').value = '';
}


function clearFileName(e) {
    e.target.value = '';
}


function processHTMLParams()
{
    var urlParams = location.search.substring(1).split("&");
    if (urlParams === '') {
        return;
    }

    urlParams.forEach(function(el) {
       var param = el.split('=')[0];
       var val = el.split('=')[1];
       switch (param) {
           case 'dna':
               document.getElementById('dnaInput').value = val;
               dnaInput();
               break;
           case 'tab':
               tab_click(parseInt(val));
               break;
           default:
               break
       }

    });

    /*
    l = unescape(temp[1]);
    temp = urlParams[1].split("=");
    p = unescape(temp[1]);
    document.getElementById("log").innerHTML = l;
    document.getElementById("pass").innerHTML = p;
    */
}

function displayContents(contents) {
    var element = document.getElementById('file-content');
    element.innerHTML = contents;
}





//tab_click(1);
//tab_click(0);


//colourDNATest('CTGGTTTACTGGTTCAACCCAGCTGGTAACCCAGCGGTAGATATGGCTGCTCTGGTCGCGATCGGCTTTTTGATCTACGGCCCTGTGATGCTGATTGGCCTTTACGCTCTGGAACTGGCTCCGAAGAAAGCCGCCGGTACCGCAGCAGGTCTGACTGGTCTCTTTGGCTACTTAGGTGGTGCTGTGGCAGCTAACGCGATATTGGGCTATACCGTTGACCACTTCGGTTGGGATGGCGGCTTCATGGTCTTGGTTGCCTCTTGTGTACTCTCAGTGCTCTGCTTGATTTACGCTTTCGTTGGCGAACGCGCTCACCATAACGATAAGCTTAAACAAGCGACCATTTAAGAGCGTGTGTCATCAAGCAATAACATCGTCAGTATAGCTTGTCCCCCTCGCCGCGCTCTTACCTGATTAGGAGCGCGGCAACTCTCTCTCCCACCCTAACGTCACGCGCCATTTGACGCCTTAGGGCCGGGAAGAACAAAAAAGGAATTCATGATGCTAAAACCATTCTCGCTTTCTCTGCTGGCTCTGGCCTGCTCAACGTCTTTATTTTCAAGCATTGTTTCTGCAGAACCAATAGTGATTGCTCACCGTGGAGCCTCTGGCTACTTGCCAGAACATACATTAGAAGCCAAAACACTGGCTTATGCGATGAAACCGGATTACATCGAGCAAGATGTGGTGATGACCAAAGACGATCAATTGGTGGTATTACATGATCACTATTTGGATCGCGTCACCGATGTTGCGGAGCGTTTCCCTAACCGCGCACGAGCCGATGGCCGTTATTACGCGATTGACTTTACCTTAGCCGAAATCAAAACCCTGCGTGTCACGGAAGGGTTTGATATTGATGCGCAAGGCAATAAAGTCGCCGGTTTTCCTGATCGTTTCCCTCTTTGGAAAGGGGATTTCACTGTCCCGACTCTCGCAGAAGAAATTGAGCTGATTCAAGGGCTCAATAAAACGCTCGGTTACAACATTGGTATCTACCCTGAAATCAAAGCACCTTGGTTTCACCGCCACGAAGGCAAAGATATTTCTCAAGCCGTCCTCAAAGTGCTAAAGCAGTATGGTTACGACAGCAAAGACGACAAAATCTATCTGCAATGTTTTGACCCTATCGAGCTAAAACGCATTAATGATGAGCTGCTTCCTGCAATGAAGATGGATCTCAATCTAGTTCAGCTGCTGGCGTACACCGACTGGAACGAAACCATGGTTTATCAAGCAGACCAAGCCACACCTTATGACTATGACTGGATGTTTGCCGAAGGCGGCATGGCCAAAGTCGCGCAATACGCCGATGGCATCGGCCCTTGGAAACCTATGTTGGTTGATGATGCTTCCACCAAAGACAACATCATGATTAAGCCGTTGATGAAGCAAGCAAAAGAGGCTGGCCTCGTGGTGCATCCTTACACTTTCCGAGCCGATAAAGGGCGCATTGCACCTTGGGCAGACAACTTCGAAGGGATGTTGGATGTGTTCTACAACCAAGTGAAAGTCGATGGTCTGTTCACTGACTTCCCCGATAAAGCGGTGGCTTTCTTAAATCAATAATGCTTTAAATTGACATTATCGCCAGTCAGCGAATAAGCGAGGGACACCTCGCTTATTCAAAGCTGGCCACAAAAGCACGCAGAGCTTCAAGATCTTTGATGATAAGTCCACGCTCACCTTTTTCAACCAATCCCCTATCCACCAACTCTTTCACCGCTCGGCGATACACTCGGCTCGATGTACCAAAGCGTTCCGCTTCTGGCTCCATCCGTTCAAAGCCACCCAAATTGACCCGAGTTTCATGCTGGAGCAGCAGATCATACGCAATGTTAAAGCTGATTGAATGCATCAAACGCTGCGTGTAAATGTCCAGCGACTCTTGGTAATCCTGCGCTAAAGCAGAAGCAAAAAACAGGCTAAACACAGGTTGCTGTTCAAGGCACTCTTGCAGCTTCTGAGCACAAATGACCATCGCTTGAATCGGTTCTTCAGCCACCACATTCCACTGACAAGGCATAGCGGTAAAAAACTCCATCTCGCCAAACAGGTGATCGTCACACTGCACTTCACCGAGTTGAAAACGGCGTCCATTAGCCGCCAGAATATGCATCGAAACTCGACCACAAGACACCACATACAGTGACTCAACAGGTTGGCCTTGCTGCAACAGCGCCTCACCCGCATCAAACCTTTTTGTCGCGGTTTGACATTGAGCCAACGCCAGACGAAACGCAGGCAGTTGCGCATTCAAATAGCCTGAGAAACGGCCGGTAGTATTCGCGGTAAGTTGCATAGCGTTTTGAAATTCGATTCACATGAGCATAGGCTACCTGTTTTTTGCTATTAGGCGAAGCTCTACCTAGATTCCTCACTCAACTAAGCAATACCGCTCATTGAGTGAGGCCATGTAATGGGTCGTGAATAGGTTGGGACTAATCGAATTTAAGCGCGTTTAATGCATGCAAAGCAGCACGTGCAGCCAGCGTTTCACCCTGCATCATTCTTTGATCAGCATTGACGTATTGATTCACGCCCTCTTGCACATCTTGCTCATACAACACCACATTATTGAGAACGGAAGCGACCTGTAACCCGCGTAAACGGCCAACCGTGAGTAGCGCAGAGGTTTCCATATCCGCCGCAAGAATCCCTTTGCGATGCCAATAACGGCACAGTTCTGCTTCTTCATCGGTATAAAAGCTGTCATGCGAGCGAACAATACCTCGGTGAATGGGTACTGACTGCTCAGCTAAAAAGCGTTGCATCTCCAGCACCAATTCAAAGCTTGAGTAAGCCGGATAAGCGGCGCCAATATATGCCTTTGAACCGCCCTCATCACGCACTGCGCCTTCCACCAAAATCAGTTCACCTAATCCGATTTCTGACTGCATTGCGCCTGCAGAACCTACTCGCACAATGGCTTTTGCACCACTTCGCGCTAACTCTTCCACCGCGATGATTATGGATGGCGCACCAATCCCTGTGCTGCATACGGTCATCGGCTGCCCTTCGAACTCACCACTGAATAAGCGGTATTCGCGGTTTTCCGCCACAAGCTTAGCGT');

/*
for (var i = 0;i < 300;++i) {

    setTimeout(function() {
        colourDNA('AGACTAACCCCATTAACTCGATGGAAGCGATGGAAGAGACTAACCTACCTAATTGCTTAGGTAGACTAACCCGATGGAAGCGATGGAAGCGATGGAAGCGATGGAAGAGACTAACCTACCTAATCGATGGAAGTACCTAATTGCTTAGGTCGATGGAAGCGATGGAAGCCATTAACTAGACTAACCAGACTAACCAGACTAACCCGATGGAAGTACCTAATCCATTAACTTACCTAATCCATTAACTTGCTTAGGTTACCTAATTGCTTAGGTCGATGGAAGTACCTAATCGATGGAAGCGATGGAAGTACCTAATCCATTAACTTACCTAATAGACTAACCAGACTAACCTGCTTAGGTTACCTAATAGACTAACCCCATTAACTCCATTAACTTACCTAATTGCTTAGGTAGACTAACCCCATTAACTAGACTAACCAGACTAACCCGATGGAAGCCATTAACTCCATTAACTCCATTAACTCCATTAACTCGATGGAAGAGACTAACCAGACTAACCTACCTAATCGATGGAAGTGCTTAGGTTACCTAATCCATTAACTCCATTAACTCCATTAACTTGCTTAGGTTGCTTAGGTTACCTAATAGACTAACCTGCTTAGGTAGACTAACCCGATGGAAGTACCTAATTACCTAATTACCTAATAGACTAACCCGATGGAAGCCATTAACTTACCTAATTACCTAATTGCTTAGGTCCATTAACTTACCTAATCGATGGAAGTGCTTAGGTTACCTAATCGATGGAAGCGATGGAAGCGATGGAAGTGCTTAGGTAGACTAACCAGACTAACCAGACTAACCAGACTAACCTGCTTAGGTTGCTTAGGT', i);
    },0);

}
*/



//var baseNode = buildTree('AGACTAACCCCATTAACTCGATGGAAGCGATGGAAGAGACTAACCTACCTAATTGCTTAGGTAGACTAACCCGATGGAAGCGATGGAAGCGATGGAAGCGATGGAAGAGACTAACCTACCTAATCGATGGAAGTACCTAATTGCTTAGGTCGATGGAAGCGATGGAAGCCATTAACTAGACTAACCAGACTAACCAGACTAACCCGATGGAAGTACCTAATCCATTAACTTACCTAATCCATTAACTTGCTTAGGTTACCTAATTGCTTAGGTCGATGGAAGTACCTAATCGATGGAAGCGATGGAAGTACCTAATCCATTAACTTACCTAATAGACTAACCAGACTAACCTGCTTAGGTTACCTAATAGACTAACCCCATTAACTCCATTAACTTACCTAATTGCTTAGGTAGACTAACCCCATTAACTAGACTAACCAGACTAACCCGATGGAAGCCATTAACTCCATTAACTCCATTAACTCCATTAACTCGATGGAAGAGACTAACCAGACTAACCTACCTAATCGATGGAAGTGCTTAGGTTACCTAATCCATTAACTCCATTAACTCCATTAACTTGCTTAGGTTGCTTAGGTTACCTAATAGACTAACCTGCTTAGGTAGACTAACCCGATGGAAGTACCTAATTACCTAATTACCTAATAGACTAACCCGATGGAAGCCATTAACTTACCTAATTACCTAATTGCTTAGGTCCATTAACTTACCTAATCGATGGAAGTGCTTAGGTTACCTAATCGATGGAAGCGATGGAAGCGATGGAAGTGCTTAGGTAGACTAACCAGACTAACCAGACTAACCAGACTAACCTGCTTAGGTTGCTTAGGT');
//('AGACTAACCCCATTAACTCGATGGAAGCGATGGAAGAGACTAACCTACCTAATTGCTTAGGTAGACTAACCCGATGGAAGCG'); //('TATTAT'); //('ACGATTCACTGACT');
//alert('base node: ' + baseNode);
//traverseTree(baseNode);



/*
function countKMer(kmer,dna,inclRevCompl) {
    indexArray = [];
    var k = kmer.length;
    var revComplKmer;

    if (inclRevCompl) {
        revComplKmer = reverseComplement(kmer);

    }
    var count = 0;

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
    return indexArray;

}
*/

/*
function mostFrequentKMers(k,dna,inclRevCompl,ltClumpThreshold) {

    //If ltClumpThreshold is blank: only return most frequent kmers
        //otherwise: return all kmers which occur >= ltClumpThreshold



    var kMersDone = [];
    var mostCount = 0;

    mostArray = [];


    for (var i = 0;i < dna.length - k + 1; ++i) {
        var count = 0;
        var kmer = dna.substring(i,i+k);
        var revComplKmer;

        var foundAlready = false;
        for (var ii=0; ii<dna.length; ii++) {
            if (kmer === kMersDone[ii]) {
                foundAlready = true;

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
        //console.log('kmer: ' + kmer);

        var indexArray = countKMer(kmer,dna,inclRevCompl);

        if (ltClumpThreshold > 0){
            if (indexArray.length >= ltClumpThreshold) {
                mostArray.push(indexArray);
            }
        }
        else {
            if (indexArray.length == mostCount) {
                mostArray.push(indexArray);

            }
            else if (indexArray.length > mostCount) {
                mostArray = [];
                mostArray.push(indexArray);
                mostCount = indexArray.length;
            }
        }



    }
    //console.log('most array: ' + mostArray);




    return mostArray;


}
*/
/*
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
*/

/*
function delKmerFromTree(node,seq) {
    //deletes one occurence of kmer and its subkmers from tree (eg kmer ACGTG - deletes one A, AC, ACG, ACGT, ACGTG)
    var currNode = node;
    for (var i = 0;i < seq.length;++i) {
        for (var j = 0;j < currNode.subNodes.length;++j) {
            if (currNode.subNodes[j].letter === seq[i]) {
                currNode = currNode.subNodes[j];
                currNode.positions.shift(); //assumes entry is in first position
                break;
            }
        }
    }

}
*/

//Actions

function findSequenceInTreeNaive(node,seq) {


    if (node.sequence === seq) {
        return node.positions;

    }
    else {

        var found = false;
        var foundPos = [];
        node.subNodes.forEach(function (sn) {
            var ret = findSequenceInTreeNaive(sn, seq);
            if (ret) {
                found = true;
                ret.forEach(function (el) {
                    foundPos.push(el);

                });
            }

        });
        if (found) {
            return foundPos;
        }
        else {
            return false;
        }
    }


}

/*

function traverseTree(baseNode,node,k,ltClumpThresh,inclRevCompl) {


    if (node.sequence.length == k) {
        var posnsToUse = [];
        node.positions.forEach(function(el) {
            posnsToUse.push(el);

        });
        if (inclRevCompl) {
            var revCompl = reverseComplement(node.sequence);

            if (revCompl === node.sequence) { //palindrome
            }
            else if   (revComplAlreadyDone.indexOf(revCompl) > -1) {
               // alert('already done: ' + node.sequence + ' ' + revCompl);

            }
            else {
                var revComplPosns = findSequenceInTree(baseNode, reverseComplement(node.sequence));
                if (!revComplPosns) {

                }
                else {

                    revComplPosns.forEach(function (el) {
                        posnsToUse.push(el);

                    });
                    revComplAlreadyDone.push(node.sequence);
                }
            }

        }

        if (ltClumpThresh > 0) {
            if (posnsToUse.length >= ltClumpThresh) {
                mfk.push(posnsToUse.map(function (n) {
                    return n - node.sequence.length + 1;

                }));

            }
        }
        else {
            if (mfk.length > 0) {
                if (posnsToUse.length > mfk[0].length) {
                    mfk = [];
                    mfk.push(posnsToUse.map(function (n) {
                        return n - node.sequence.length + 1;

                    }));
                }
                else if (posnsToUse.length == mfk[0].length) {
                    mfk.push(posnsToUse.map(function (n) {
                        return n - node.sequence.length + 1;

                    }));
                }
                else {

                }

            }
            else {
                mfk.push(posnsToUse.map(function (n) {
                    return n - node.sequence.length + 1;

                }));

            }
        }
        return;

    }



    node.subNodes.forEach(function (sn) {
        traverseTree(baseNode,sn,k,ltClumpThresh,inclRevCompl);

    });





}
*/

/*
function Node(letter) {
    this.letter = letter;
    this.sequence = letter;
    this.positions = [];

    this.subNodes = [];

}
*/
/*
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
*/

/*
function buildTree(dna,k,ltClumpThresh,includeRevCompl) {



    var justAddedPrev = [];
    var justAddedCurr = [];

    var baseNode = new Node('');

    justAddedCurr = [baseNode];
    for (var i = 0;i < dna.length; ++i) {

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
*/


//Execute stuff in gui

function skTimeoutLoop(dna,numAtATime) {
    //var basesDone = 0;
    //if (sk) {
    //    basesDone = sk[2].length - 1;
   // }
   // while (basesDone < dna.length) {
   
    sk = gcSkew(dna,numAtATime, sk);


    skewCanvas(sk);

    var basesDone = sk[2].length - 1;
    if (basesDone >= dna.length) {
        //return sk;
         var stats = collectStats(dna,sk);
         document.getElementById('rightTwo').innerHTML = stats;


    }
    else {
        setTimeout(function() {
            skTimeoutLoop(dna,numAtATime);
        },0);

    }


}

function mostFrequentKmersBrute(e,variant) {

    initialiseResults();

    var numKmer = document.getElementById('numKmers');
    numKmer.value = 0;

    var mfThresh = 0;

    document.getElementById('stats').style.display = "block";
    //var inDNA = document.getElementById('dnaInput');

    // document.getElementById('dnaView').innerHTML = inDNA.value;

    //document.getElementById('dnaLength').innerHTML = inDNA.value.length;

    var radio = document.getElementById('mostMost');
    if (radio.checked) {

        mfThresh = 0;
    }
    else {
        mfThresh = parseInt(document.getElementById('mfThresh').value);

    }

    var maxMismatch = parseInt(document.getElementById('maxMismatchMF').value);




    k = parseInt(document.getElementById('kMerLenMF').value);

    w.postMessage({'task' : 'mostFrequentKmers', 'method' : 'brute', 'variant':variant,'dna' : dnaMaster, 'k':k,'inclRevCompl' : document.getElementById('includeRevComplMF').checked, 'maxMismatch':maxMismatch ,'ltClumpThreshold':mfThresh,'debug':document.getElementById('debugMF').checked }); // Start the worker.
    // mfk = mostFrequentKMers(k,inDNA.value,document.getElementById('includeRevComplMF').checked,mfThresh);

    /* moved to worker message listener:
     document.getElementById('kMerOccurrences').innerHTML = mfk[0].length;
     document.getElementById('kMerNumEqualMostFreq').innerHTML = mfk.length;

     console.log('mfk: ' + mfk + ' len: ' + mfk.length);
     mfk.forEach(function(m) {
     console.log('mfk entry: ' + m);
     });

     var numKmer = document.getElementById('numKmers');
     numKmer.max = mfk.length;
     if (mfk.length == 0) {
     numKmer.min = 0;
     }
     else {
     numKmer.min = 1;
     }


     var numKmerVal = parseInt(numKmer.value);
     var dna = document.getElementById('dnaInput').value;

     var km = dna.substring(mfk[numKmerVal-1][0],mfk[numKmerVal-1][0] + k);
     document.getElementById('kMerMostFreq').innerHTML = dna.substring(mfk[numKmerVal-1][0],mfk[numKmerVal-1][0] + k);
     document.getElementById('kMerRevCompl').innerHTML =  reverseComplement(km);



     revComplAlreadyDone = [];
     colourDNA(inDNA.value);
     */


}


function mostFrequentKmersTree() {

    initialiseResults();

    var numKmer = document.getElementById('numKmers');
    numKmer.value = 1;

   // var dnaTest = document.getElementById('dnaInput').value;

    var mfThresh = 0;

    var maxMismatch = parseInt(document.getElementById('maxMismatchMF').value);


    //document.getElementById('dnaView').innerHTML = inDNA.value;

    //document.getElementById('dnaLength').innerHTML = inDNA.value.length;

    var radio = document.getElementById('mostMost');
    if (radio.checked) {

        mfThresh = 0;
    }
    else {
        mfThresh = parseInt(document.getElementById('mfThresh').value);

    }



    document.getElementById('stats').style.display = "block";
    //var inDNA = document.getElementById('dnaInput');
    //document.getElementById('dnaView').innerHTML = inDNA.value;
   // document.getElementById('dnaLength').innerHTML = inDNA.value.length;

    // ltClumpThresh = parseInt(document.getElementById('ltClumpThresh').value);

    k = parseInt(document.getElementById('kMerLenMF').value);

    w.postMessage({'task' : 'mostFrequentKmers', 'method' : 'tree','dna' : dnaMaster, 'k':k,'inclRevCompl' : document.getElementById('includeRevComplMF').checked,'maxMismatch':maxMismatch, 'ltClumpThreshold':mfThresh }); // Start the worker.

/*
    var added = buildTree(dnaTest,k,mfThresh,document.getElementById('includeRevComplMF').checked);
    baseNode = added[0];
    //var justAdded = added[1];
    mfk = [];
    revComplAlreadyDone = [];
    traverseTree(baseNode,baseNode,k,mfThresh,document.getElementById('includeRevComplMF').checked);
*/

    // document.getElementById('kMerOccurrences').innerHTML = mfk[0].length;
    //document.getElementById('kMerNumEqualMostFreq').innerHTML = mfk.length;

    /*
     var numKmer = document.getElementById('numKmers');
     numKmer.max = mfk.length;
     if (mfk.length == 0) {
     numKmer.min = 0;
     }
     else {
     numKmer.min = 1;
     }

     var numKmerVal = parseInt(numKmer.value);
     */
    var dna =  dnaMaster; //document.getElementById('dnaInput').value;

    if (mfk.length > 0) {
        var km = dna.substring(mfk[0][0], mfk[0][0] + k);
        document.getElementById('kMerMostFreq').innerHTML = dna.substring(mfk[0][0], mfk[0][0] + k);
        document.getElementById('kMerRevCompl').innerHTML = reverseComplement(km);
    }



    colourDNA(dnaMaster);



}
/*
function timeOutLoop(dna,l,k,curr,useTree,inclRevCompl) {
    // console.log('st timeout');

    //colourDNAProgressNew(dna, [curr], clumps,l);
    var limit = dna.length - l;
    document.getElementById('progress').innerHTML = curr + ' / ' + limit;

    ltClumpThresh = parseInt(document.getElementById('ltClumpThresh').value);
    var mf;
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
        //console.log('st traverse');
        traverseTree(baseNode,baseNode,k,ltClumpThresh,inclRevCompl);
        mf = mfk;
    }
    else {
        mf = mostFrequentKMers(k, dna.substring(curr, curr + l), inclRevCompl, ltClumpThresh);
    }

    //console.log('curr: ' + curr + ' ' + mostFrequentKMers(k,dna.substring(curr,curr+l),document.getElementById('includeRevCompl').checked,ltClumpThresh));
    //console.log('curr: ' + curr + ' mf: ' + mf.length + ' mf innards: ' + mf[0] + ' tot mf: ' + mf);
    // console.log('st adjust');
    if (mf.length > 0) {

        if (useTree) {

        }
        else {
            var adjMf = mf.map(function (el) {
                adjEl = el.map(function (e) {
                    return e + curr;
                });
                return adjEl;

            });
            //// console.log('curr: ' + curr + ' mf: ' + mf.length);
            mf = adjMf;
        }
    }
    clumps.push(mf.length > 0);

    //console.log('st add to allmfk');
    mf.forEach(function(el) {
        var alreadyThere = false;
        for (var i = 0;i < allMFK.length;++i) {
            if (dna.substring(el[0],el[0] + k) === dna.substring(allMFK[i][0],allMFK[i][0] + k)) {
                alreadyThere = true;
                break;



            }
        }
        if (alreadyThere) {
            //No need to add. Note: this entry may actually have more frequency then prev one, so may need to tweak this
        }
        else {
            allMFK.push(el);
        }


    });



    if ((curr >= dna.length - l) || stop) {
        mfk = [];
        allMFK.forEach(function(el) {
            mfk.push(el);

        });
        colourDNA(dna,clumps);

    }
    else {
        setTimeout(function() {
            timeOutLoop(dna,l,k,curr+1,useTree,inclRevCompl);
        },0);

    }

}
*/

function runLTClump() {
    // var dnaTest = 'CTTGATCAAGACGTCCAACTAATGACAACTCAAAGTTCGCACCCATGTACTTAAAGTGAGGACGAACCGAGAGCGACACCTGATACCAACCTGGATTGCCTTCCACATCCATCACTTCAATGCGCGCAGCACGAAGTGGACGACGGCTACGTACGTCTGCAGGTGGGTTCTCTTGATCAGCAACGTATTGTTTGATCCATGAGTTCAGTTCACGCTCAAGATCTTGACGCTCTTTCCAAGCACCGATCTGCTCACGTTGCAGAACTTTCACATAGTGCGCCAAACGGTTGATGATCATCATGTACGGCAACTGGGTACCCAACTTGTAGTTGGTTTCCGCTTCTTTGCCTTCTTTGGTATTTGGGAAAACCTTAGGTTTTTGAATGGAGTTTGCAGAGAAGAACGCCGCGTTATCACTGCCTTTACGCATAGTAAGAGCAATAAAACCTTCTTCCGCCAGTTCAAACTCTTTACGGTCCGTGATCAGGACTTCGGTTGGGATCTTGCTTTGCAATGCACCCATAGATTCAAAGACATGCACCGGCAGATCTTCAACTGCACCACCACTTTGTGGACCGATAATGTTTGGACACCAGCGATATTTAGCAAAGCTATCCGTCAAACGAGTTGCGAAGGCAAATGCCGTGTTACCCCACAGGTAGTGCTCGTGCGAAGCACTGACGTTTTCCGCATAATTGAACGACTTCACTGGATTTTCGATTGGATCGTAAGGAACACGCAGCAGGAAACGAGGCGCAGTCAAACCAAGATAGCGCGCATCTTCCGATTCACGCAATGAACGCCATTTGGTGTATTTCGGGCTTTCAAATGTCGACTTGAGATCTTTAATGTTAGGCAGTTCTTCAAAAGAATCGATACCAAAGAATTCAGGACCTACGCTTGAAATGAAAGGAGCATGCGCCATGGCACCCAGTGCGCCCATGTATTGCAGCAGCTTCATATCTGGCGTTGAAGGGGTAAACGCATAGTTACCAATGATCGCGCCAACAGGTTCGCCACCAAATTGACCATAACCGGCAGAATAAACGTGCTTGTAAAGACCGGACTGAGCCGTTTCTGGAGCAAACTCGAAATCTTCCAGCAGTTCATCTTTGGTTACGTGAAGGATTTCGACTTTGTTATTTTCACGAAAATCAGTGCGATCCACGAACAGCTTCAAACCGCGCCACGCCGATTCCATCGCTTGAAATTGCGAGTTGTGCAGGATTTCATCCATCTGTGCACTGATTTTCTTATCCAGTTCAACCAACATTTGGTCAACCAGAGATTTGTTGACAGGCTCAGCAGAGTGTTGTGAACCCATAAGATTTTCGATAAACGCTGCAACACCTTTTTTCGCGATGTCGTAACCCTCTTCGCTTGGTGCGATACGGGTTTGCGCCATAATTTCATCAAGAAGGCTGCCTTGAGCAAGCTGTGGCCTTTCCAATACCTTTTCAGTCGTAGACATCATTAAGTTCCTAATGAGTGATTAATTGATTAATAACGTTTGCTTACGCTTGTGGCTCTTCTTGACCACTGAGCAGATTCAGTTCTGCCAACAGTTTTTCTCTCGACTCTTCTGAGTTGAGTAATGACTGTAAACGCTCACGAAATGCGGGAATGTTGCCTAGCGGCCCTTTAAGGGCAACTAACGCTTCACGCAACTCAATCAATTTTTTCAGTTCTGGAACTTGTGATGCCACCGCATCAGGAGCGAAGTCGGCTAAGGATTTGAAATTGAGTTCAACAGGAAGCTCGGCATTCTCATCATCAGTCAGCTTGTTTTTCACCGTGGCGGTGATTTTCAGCTCGCTCTCGCGCATTACGGCTTCAAAGTTGTTCTTAATCTACCGTGACTGTTGCACGCTCTTCCAATGGGGTTTGCTCCGCATGCCCTTTGAAATCACCTACAACTAGGGTTTTGAGTGGTAGCTCAACCTCAGCCTGTGCATCCCCCGTCGCCGGAATATACTTGATATTAATCCGCTCTTTGGGAGCTACACTTCCTTCTTTAGACATATTACGTCTCCAATACCTATGCCAAACGTTGTCAATGAAACAGTTGATTGAAATCATTCAATCAACCAAGAATCTTGAACTTATCCTTAGTGATAACTAAGGGCTCTCATTTCAAAATATGAAATAGAGTCAATCACTTCCGTGACACTGTGAGAAAATAAATCCCACACTAAGCCATGTCACAGAGTGTTGTTGTATTTTGATATATTGATGAAAAATGGGCTTATTATGCAAACTTATCATTAGCAGCTTTTTGCTGCCATGTATCAGTTATCTGGTATAAATATTGTTTGCGGTGATCTAATACCCAATAAATATGAGATGGCATTCTATTCATGGAGACAAAAAAGTAATCGAAAAAAAACCCGAAGCGTGATCTTTTACCCAAAAATATATCCAATAACCATACTTAAATAATACTTAACAGGGAACTAATTATCCCAGTAAGCATCATGGTTGAAAAAAGGCTGATAATTTTCGAGTTGAAGCTCTATTTTCTCTGCAATTTCTTTCTCATGCCCATTTAATTGCATTGAAATATCGAGCAGTTGCTTTTGGCTATTTTGGTAGACATGACTCAAGCGATCGCTGCCTGCCACTTCTAAATATCGATTGAAGTAAGCCTTGGCCTGCTTAAGATCAGGCTCGACAAAACGATAAAACTCACTCTGTCCACTGAGAATACCGCCCATCGCCACCAAAGAAGTGAGATCGCCACGCTGTACTGCCTGCTGGCGCCAATAATAGGCATCTTGGTATTTACCTTTATTTTCAAGAAGTTGAATATAGTTTCTTATCGCAGGGATATAGCCCGTTTGCGCCGCTTTCAGATACGTTTTTCGAGCTTCGGTTGCTCGGCTACCAGGCATGGGATAACTGCCCTCTCCTTGCTCAATCTGTTGCGCCATAAGCATCAGAGCTGCGGGGTGTTCAAATTCAACTGCAATCGTTAAATAATGCTGGGCAAGCTTGGGGTTATCCTTTTGATGGAACAGCGCGAGCTCATAAGCGGCCTGACTCGGTTCACGAGAGCCTAATTCAATCAAAGTATCGTAATAAGTTTGTTGCCATAAACTGCGTTCTGCAGCTCTGAGCCACTCACCGTTTTGATACAACACTCGCATCGCTTGACGGCTTCCACCTTGAGCGGCGAGCAATAAATAGTGCTGGCTTTCAGGCGGGGTTCGAATGGTAGTTCGGTAGTTGGCCACTTCCATCGCGTAAAGATAAGCCGCATTAGGCTCACCAAGGTCTGCCGCGTATTTCAAATATTCGCGAGCCGCTAAGTTCTTAAATTGCGCACGTAGCAAACGACCATTTTCGTAGGCTTGCTGAGCACTCAACTGCTGCATATCGTAGTTTTTTTCTTCGCTAAAACTACGTCCTATAAACAGAAGAGTTACGATGAAAATCAGAGTCTTAACCAATGTTGACACTGCTTGCTCCCTGTAATACGCCACCACAATCGACCGCATCTCCAACCCTTGCTGCGGGTTTTCCATCAATCATCACCGTCCCAGAACCCGCAGCAATCGCTCGGCCATGCGAGGGATGCTTTGGTTTATCGTGAGGGGCCAAGGGGTCCCCAAGCCGAGCAGCTGGAATCCCATCATAACGCACAGTTGCTGAGCCTGCGGTCACTGGCGTTGGCGGGAAACCATCATGATCTGTACCTAAATGCCCTACTACGATTCCGTTACCCATATTCTGTCCTCATCGGTTAGTAAAATGTGTTAAAGCAATTAATGTTCCATACTTTGGGAGGTCATTTTTGCTTATCAACCTCATGACTGGGCAATATATTGCCCAAAAACAACACCTCCCTGCTAGGCGGGCAAGATCTTACCTGAAAGGAAAAATTGTCATAGAAACGTCAGGAAAATATAAAACAATTCATTTTTGCGTTTCTTCTTCATTTAAAGTGAATATTAAAAGTGAATAAAGTTAATAGTGTTATTTCATTATCTCCTAATTGAGAAATCTATTAAATTAATACCCTTTCAGCCTCATAACTTAGTCGCAGAGTTTAATAGGATTAATTAGTTCATTCTTATTAAGCAATTTTTACTTAAGCGAATCCTAGACCATACATTAGAGTCAATAACCCTTTATCACGAACATTTGGCCAGAATTCACTCATATTGATGATATAACAGGCAAATACAACATTTGGAAATGATATGTTTTTGAAATGTAACATCGACCATTAATTAGGGTGAACAAATATCCCTACATTGAACTGTTGTAAGTTGAAAGAGCGATAACCAAACTACTTGGTGTTGCAGAAAGGTGGCAAGTAAGTTCATCATCATGAGCATAAATAAACTATGTGATTGGGGTGAGCGAACGTAACCCATATCGCTGCAACTTCAAGTAGGCACGGTATAGGTCAACAATCACCAATAGAAACATAGACTGACAAGTTTCACCCAAGCAGGTAAGATGCGTCTTTGCTTTACCCACTCTGTCATTTGTATGTCTGCACTGCCCGATTGCTTCACTCGTTTCACTTCTGTCATCGCTCATATCGCGCTACCCGAGCGCTTTACTTTTCCATTTTGTTATCAGCCACACCCACTGTCTGAAATTGCGGCCGATGAGCTACAGCAACATTTGCTCACTCAAACCGATTGGTATCACCCTTTTGGTTTAACGGAAACCGACCCACAAGCTCACGGCAAAATGTTTGGAGTACTGGTTGTGCAGCATCGCTCTGGAACCTTGGGTTATTTAGCGGCGTTCTCTGGTCAACTGGCTGAACAAAATCGGCTGCCGGGATTTGTTCCTCCTGTGTTTGATCGCTTTGCTGATGCGTCTTTTTTCAGAGTGGATGGCGATCAGATTGCGGCCATTAACCAACAAGTGCGAGAGCAAGAAACTGATCCGGAGCTTGCGGAACTCGCTAGTGCCCTTAAGCAATCTCAGTTACAAGCCGAGCATGAGCTGAGCCAGCGGCGTATATTGCACAACAAGCAGCGTCACTTTCGAAAACAGCAGCGCTTACAAGCAGAAGAGCTACCAGAAAAAGAAAAACACACTCTGTTGGCTCGCCTAGCGGAAGAAAGTGTGCAACAAAAGCGTCAGTTACAGTGCCTCAAACAAGAGTGGGAACAGCGTATTGCGGCTTTACAGCATAGGTTAGATCATAAACTCATGCGCATTGAGCAACTGAAACAGCAACGTAAACAGCGCTCTGCGGCGCTGCAGAAAAAACTGTTCTCCGCCTATCGTTTTACCAATATCCGTGGCGTTGAAAAGGATTTGGTTGAGCTATTTTCCGTGACCAAAAATCCGCTTCCACCCGCAGGTTCAGGCGAGTGCGCGGCTCCAAAACTATTACACTATGCGTTCCAACATCAGCTGCGACCCATTGCATTGGCAGAGTTCTGGTGGGGGCGCTCCCCCAAATCGGAAATTCGTCAGCATAAGAAGTTCTATCCTGCATGTCAGAGTAAATGCCAACCGATTTTGGCGCACATGTTAGAAGGGATGCCGCTGGAAGATAATCCACTGCTCAGCAATCCGGCGCAAGGTCAAGACATCACCATTGTGTATCAAGATGACGCCATCGTAGTGGTCAATAAGCCTGCTGAATTTCTCTCTGTCCCCGGTGTGCATGTGCATGATTCTGTGTTGACGCGCTTAAAGGCTCAATTTACCCAAGCAGAAGGGGTATTTGCCCTACATCGTCTCGATATGTCTACTTCTGGACTTTTGGTGTTTGCGCTGACGCGACGTGCTAACAAACAGTTGCAAAAGCAGTTTATTTCTCGTGCCGTCCAAAAGCGCTATGTCGCACTTATCGAAGGAAAATTAAGCGAAACACAAGGTGAGATCCAGTTGCCTTTGTGTGGCGATCTGGATGATAGACCAAGACAAAAAGTGTGCTGGCAGCAAGGTAAGCCTGCACTCACTCATTGGGAAACCGTGCAAGTGGAGCAAAACCGAACTCGTGTTTACCTCTACCCACATACAGGACGTACCCATCAATTGCGAGTACATTGCGCCCATCATCTTGGCTTAAATGCGCCGATTGTCGGTGACGATTTGTATGGCTTGCAAGACAAACGTTTGTTTTTGCATGCCGAACAACTGAGTTTTGCTCATCCCTACACTAAGCAGCCCATGACTTTCCAAGTCGATGCGGATTTTTAACCGAAAAAGAAAAGGGGCCCAATTTGAGCCCCTAACCATTCTTAAGCCACTTCGGTCGTGGTTGGTTTTAGCAACACGCTGATATCTGGTTCCTTCTGCACCACAGGGACATGGCCTTCTTCCATCGCACTTAGGTAACAACCACTCAGTTTGATGGCTAAATACATATCGGCACGCAGCACGATTCCCCCTTCAGCGGTATTGATCCCGGTGAGATCACACCAATGGCGTAAACTGCGCTTACGCCCCGCTAGAATGAGACGGATACCGCGTTTTTTCAGGATGCCATGTAAGTCAGCCAGCATCGCCATCACACTCAGATCTAAGTGGGTAAAACTGGCCACGGCATCAATGATCACACAACCCACTTGGGCGCCTTCACGCTCCGTCTGATCGAGAATACGCCGCTTAAAATAGGGCGCGTTAAAATAGGTCAGCGGTGAATTGAAACGGAAAATCACCATGCCGGGGATCGGTTTTGCTTTCTCTGAACCATCGAGTGTGCGCAGCGTTCCTTCCTCATCCAATCCCATCATTTGATCAGTGGGGCGCATCACCAATTTGAGAAATTGGAATAAACCCAGCAGTACCGCCAAAGTGATCCCTGGGATCACCCCAATCACTAACACTGCGATAAAGGTGATGAGTGCGAGGTAAAACGCATCTTTATCGCGCTTTCTTAAATTCCAAACGCCTTTTAAATCGAGCAGAGACAATGAGGCGATGATCAGCACCACGCCCAAAGCCGCGACAGGGATAAATTGCAAAGGCTGATAAGCAAAAACCGCCACCAGCGCAATAAACAAGGCTGCGATGACCGAAACCAGTTGCGATTTACCACCATTGGCGTCATTGACTGCGGTACGGGAATCGGCACCACTGATCGCAAAGCCTTGCGAAAAGGCGGCTGCGACATTCGCGACACCAAGCGCTCGGAACTCTTTATCCGCATCAATATCGTAGCCATTTTTAGCCGCAAAGCTGCGAGCGGTGAGCATCATACTCACAAAGCTCACCATCGCTAAGTTGAGCGCAGGCATAACGAGTTCACGGCTGATCCCTAAATCAAAAGCCGGAGCTTGGAACTCAGGTAAACCGCCTTGGATAACCCCCACCACTTGTACGCCCACACTTTCAAGATTCAGTGCCCACACCAGCAAGGCCGCGACCATGATTGCAAACATCGCCGCAGGCCATCTCGGCTGCCAGCGCTTGATCACTAAGTAAATCGCTAAAGTGAGCGCACTCAAACCCAAGGTTTGCCAATGCAAGGAGTACAGCAGCTCTGGCGCTTCAACAATCCGTTCAAGCAGATAACGTTTTTCATACTTGAGCCCCAATACTTTGGCAAACTGACCAACAATAATCGTCAGCGCCACCCCATTAAGCAGACCTAACAAAATCGGCCTTGAGAGAAAATCGGCGAAAATGCCCAGTTTGAGCCGACTCGCCAATATGCACCAAAAACCGGTCATGGCGGTCATGGTCATCACCAGTTGCCAGTGCTTGGTGGTATCTCCGGCAGCAAGAGGGGTAACCACGGCCGCAATCACCGCGCAAGTGGCGGCATCTGGCCCCACAATCAATTGGCGTGAGGTGCCCATCAGCGCATAGACCAACATAGGCAAAACACACGAATAAAGTCCGACAATAGCAGGCACACCGGTTAATTGCGCATAAGCGATCGCTACAGGCAGCGCTACAGCCACAACAGAAAAGGCGGCTCGTACATCATCGGTTAACCATCCCCGTTGATAATCTTTAAATTGGTAGAGTCCGGGGAACCACTGCCTTATCCACACTGCTTTCACTTCATTCACCTGCTGACACTAAATTGTCATTGTTTTGCAACATAGTTTAGCGGAAATTGCCGCAACTCTCTCGCTAATTTACTGATCTAAAATACCCTTATCTACTGATCCAACCACTACAGTTCTAAGTTTTGGGTTCTCTTGCATCAGGAAAGTATGGGACAATAGCCGCTCACTTTTTTGTGGTAGGTGGGAACGATGAACTTTCAAGCCGCTATTTTTGATATGGATGGTCTACTACTCGA';
    //var dnaTest = 'AGACTAACCCCATTAACTCGATGGAAGCGATGGAAGAGACTAACCTACCTAATTGCTTAGGTAGACTAACCCGATGGAAGCGATGGAAGCGATGGAAGCGATGGAAGAGACTAACCTACCTAATCGATGGAAGTACCTAATTGCTTAGGTCGATGGAAGCGATGGAAGCCATTAACTAGACTAACCAGACTAACCAGACTAACCCGATGGAAGTACCTAATCCATTAACTTACCTAATCCATTAACTTGCTTAGGTTACCTAATTGCTTAGGTCGATGGAAGTACCTAATCGATGGAAGCGATGGAAGTACCTAATCCATTAACTTACCTAATAGACTAACCAGACTAACCTGCTTAGGTTACCTAATAGACTAACCCCATTAACTCCATTAACTTACCTAATTGCTTAGGTAGACTAACCCCATTAACTAGACTAACCAGACTAACCCGATGGAAGCCATTAACTCCATTAACTCCATTAACTCCATTAACTCGATGGAAGAGACTAACCAGACTAACCTACCTAATCGATGGAAGTGCTTAGGTTACCTAATCCATTAACTCCATTAACTCCATTAACTTGCTTAGGTTGCTTAGGTTACCTAATAGACTAACCTGCTTAGGTAGACTAACCCGATGGAAGTACCTAATTACCTAATTACCTAATAGACTAACCCGATGGAAGCCATTAACTTACCTAATTACCTAATTGCTTAGGTCCATTAACTTACCTAATCGATGGAAGTGCTTAGGTTACCTAATCGATGGAAGCGATGGAAGCGATGGAAGTGCTTAGGTAGACTAACCAGACTAACCAGACTAACCAGACTAACCTGCTTAGGTTGCTTAGGT';

    initialiseResults();

    var numKmer = document.getElementById('numKmers');
    numKmer.value = 1;

    document.getElementById('stats').style.display = "block";
   // var dnaTest = document.getElementById('dnaInput').value;
   // document.getElementById('dnaLength').innerHTML = dnaTest.length;


    k = parseInt(document.getElementById('kMerLenLT').value);
//timeOutLoop('AGACTAACCCCATTAACTCGATGGAAGCGATGGAAGAGACTAACCTACCTAATTGCTTAGGTAGACTAACCCGATGGAAGCGATGGAAGCGATGGAAGCGATGGAAGAGACTAACCTACCTAATCGATGGAAGTACCTAATTGCTTAGGTCGATGGAAGCGATGGAAGCCATTAACTAGACTAACCAGACTAACCAGACTAACCCGATGGAAGTACCTAATCCATTAACTTACCTAATCCATTAACTTGCTTAGGTTACCTAATTGCTTAGGTCGATGGAAGTACCTAATCGATGGAAGCGATGGAAGTACCTAATCCATTAACTTACCTAATAGACTAACCAGACTAACCTGCTTAGGTTACCTAATAGACTAACCCCATTAACTCCATTAACTTACCTAATTGCTTAGGTAGACTAACCCCATTAACTAGACTAACCAGACTAACCCGATGGAAGCCATTAACTCCATTAACTCCATTAACTCCATTAACTCGATGGAAGAGACTAACCAGACTAACCTACCTAATCGATGGAAGTGCTTAGGTTACCTAATCCATTAACTCCATTAACTCCATTAACTTGCTTAGGTTGCTTAGGTTACCTAATAGACTAACCTGCTTAGGTAGACTAACCCGATGGAAGTACCTAATTACCTAATTACCTAATAGACTAACCCGATGGAAGCCATTAACTTACCTAATTACCTAATTGCTTAGGTCCATTAACTTACCTAATCGATGGAAGTGCTTAGGTTACCTAATCGATGGAAGCGATGGAAGCGATGGAAGTGCTTAGGTAGACTAACCAGACTAACCAGACTAACCAGACTAACCTGCTTAGGTTGCTTAGGT',200,3,0);
    clumps = [];
   // allMFK = [];

    stop = false;

    var clumpSize = parseInt(document.getElementById('ltClumpSize').value); //L
    var ltClumpThresh = parseInt(document.getElementById('ltClumpThresh').value); //T

    var maxMismatch = parseInt(document.getElementById('maxMismatchLT').value);


    w.postMessage(
        {'task' : 'ltClump',
         'useTree' : document.getElementById('treeLT').checked,
         'dna' : dnaMaster,
         'k': k,
         'clumpL' : clumpSize,
         'clumpT' : ltClumpThresh,
         'inclRevCompl': document.getElementById('includeRevComplLT').checked,
         'maxMismatch': maxMismatch
        }); // Start the worker.


   // timeOutLoop(dnaTest, clumpSize, k, 0,document.getElementById('treeLT').checked,document.getElementById('includeRevComplLT').checked);

}

function searchKmer(e) {

    initialiseResults();

    var numKmer = document.getElementById('numKmers');
    numKmer.value = 1;


    var kmer = document.getElementById('kmer').value;
    k = kmer.length;

    var dna = dnaMaster; //document.getElementById('dnaInput').value;

    var inclRevCompl = document.getElementById('includeRevComplKS').checked;

    var maxMismatch = parseInt(document.getElementById('maxMismatchKS').value);

    w.postMessage({'task' : 'searchKmer', 'dna' : dna, 'kmer':kmer, 'inclRevCompl': inclRevCompl,'maxMismatch': maxMismatch}); // Start the worker.

    /*
     setTimeout(function() {
     var indexArray = countKMer(kmer,dna,document.getElementById('includeRevComplKS').checked);
     //alert ('ind: ' + indexArray);
     console.log('normal finished');
     mfk = [indexArray];
     k = kmer.length;
     colourDNA(dna);

     },0);
     */
}


function EdgePrinter(edge,nodePSt,nodePEnd,backNum,forwardNum) {
    this.stX = nodePSt.x;
    this.stY = nodePSt.y - 15;
    this.endX = nodePEnd.x;
    this.endY = nodePEnd.y - 15;
    this.point1 = 0;
    this.point2 = 0;
    this.edge = edge;
    this.nodePSt = nodePSt;
    this.nodePEnd = nodePEnd;

    this.numBackward = backNum;
    this.numForward = forwardNum;

    if (nodePSt.linearNum > nodePEnd.linearNum) {
        this.stY += nodePSt.radius * 2;
        this.endY  += nodePEnd.radius * 2;
        this.point1 = this.stY + (NodePrinter.pointInc * this.numBackward);
        this.point2 = this.endY + (NodePrinter.pointInc * this.numBackward);

    }
    else {
        if( (nodePSt.linearNum == nodePEnd.linearNum - 1) && this.numForward == 0) {//adjacent forward
            this.point1 = this.stY - (NodePrinter.pointInc * this.numForward);
            this.point2 = this.endY - (NodePrinter.pointInc * this.numForward);
            this.stX = this.stX  + nodePSt.radius;
            this.endX = this.endX - nodePEnd.radius;
            this.stY = this.stY + nodePSt.radius;
            this.endY = this.endY + nodePEnd.radius;
            this.point1 = this.stY;
            this.point2 = this.endY;

        }
        else {
            this.point1 = this.stY - (NodePrinter.pointInc * this.numForward);
            this.point2 = this.endY - (NodePrinter.pointInc * this.numForward);

        }
    }

    this.print = function(ctx,walkNum) {


            var col = walkNum % NodePrinter.walkColours.length;

            ctx.strokeStyle = NodePrinter.walkColours[col];

            ctx.moveTo(this.stX,this.stY);
            ctx.bezierCurveTo(this.stX,this.point1,this.endX,this.point2,this.endX,this.endY);



    }



}
function NodePrinter(node,x,y,num) {
    this.x = x;
    this.y = y;
    this.radius = 15;
    this.node = node;

    this.linearNum = num;
    this.numForward = 1;
    this.numBackward = 1;

    this.edgePrinters = [];




    var k = node.dna.length;

    if (k <=3) {
        this.font = '12pt Arial';
    }
    else if (k <= 4) {
        this.font = '10pt Arial';

    }
    else  if (k <= 5) {
        this.font = '8pt Arial';
    }
    else if (k <= 6) {
        this.font = '6pt Arial';
    }
    else {
        this.font = '5pt Arial';
    }




    this.print = function(ctx) {
        ctx.font = this.font;

        ctx.moveTo(this.x + this.radius,this.y);

        ctx.arc(this.x, this.y, this.radius, 0, 2 * Math.PI, false);
        ctx.fillText(this.node.dna,this.x - 15,this.y+5);
        if (this.node.pairedDna) {
            ctx.fillText(this.node.pairedDna,this.x - 15,this.y+15);
        }

    };

    this.join = function(ctx,endPNode,walkNum) {
        var stX = this.x;
        var stY = this.y - 15;
        var endX = endPNode.x;
        var endY = endPNode.y - 15;
        var point1,point2;

        if (this.linearNum > endPNode.linearNum) {
            stY += this.radius * 2;
            endY  += endPNode.radius * 2;
            point1 = stY + (NodePrinter.pointInc * this.numBackward);
            point2 = endY + (NodePrinter.pointInc * this.numBackward);
            ++this.numBackward;
        }
        else {
            if( (this.linearNum == endPNode.linearNum - 1) && this.numForward == 0) {//adjacent forward
                point1 = stY - (NodePrinter.pointInc * this.numForward);
                point2 = endY - (NodePrinter.pointInc * this.numForward);
                stX = stX  + this.radius;
                endX = endX - this.radius;
                stY = stY + this.radius;
                endY = endY + this.radius;
                point1 = stY;
                point2 = endY;
                ++this.numForward;
            }
            else {
                point1 = stY - (NodePrinter.pointInc * this.numForward);
                point2 = endY - (NodePrinter.pointInc * this.numForward);
                ++this.numForward;
            }
        }

        var col = walkNum % NodePrinter.walkColours.length;

        ctx.strokeStyle = NodePrinter.walkColours[col];

        ctx.moveTo(stX,stY);
        ctx.bezierCurveTo(stX,point1,endX,point2,endX,endY);


    };

    this.joinEdges = function(ctx,walkNum) {
        var stX = this.x;
        var stY = this.y - 15;
        var endX = endPNode.x;
        var endY = endPNode.y - 15;
        var point1,point2;

        if (this.linearNum > endPNode.linearNum) {
            stY += this.radius * 2;
            endY  += endPNode.radius * 2;
            point1 = stY + (NodePrinter.pointInc * this.numBackward);
            point2 = endY + (NodePrinter.pointInc * this.numBackward);
            ++this.numBackward;
        }
        else {
            if( (this.linearNum == endPNode.linearNum - 1) && this.numForward == 0) {//adjacent forward
                point1 = stY - (NodePrinter.pointInc * this.numForward);
                point2 = endY - (NodePrinter.pointInc * this.numForward);
                stX = stX  + this.radius;
                endX = endX - this.radius;
                stY = stY + this.radius;
                endY = endY + this.radius;
                point1 = stY;
                point2 = endY;
                ++this.numForward;
            }
            else {
                point1 = stY - (NodePrinter.pointInc * this.numForward);
                point2 = endY - (NodePrinter.pointInc * this.numForward);
                ++this.numForward;
            }
        }

        var col = walkNum % NodePrinter.walkColours.length;

        ctx.strokeStyle = NodePrinter.walkColours[col];

        ctx.moveTo(stX,stY);
        ctx.bezierCurveTo(stX,point1,endX,point2,endX,endY);


    };



    this.path = function(ctx,endPNode,walkNum,edgeNum) {
        var stX = this.x;
        var stY = this.y - 15;
        var endX = endPNode.x;
        var endY = endPNode.y - 15;
        var point1,point2;

        if (this.linearNum > endPNode.linearNum) {
            stY += this.radius * 2;
            endY  += endPNode.radius * 2;
            point1 = stY + (NodePrinter.pointInc * edgeNum);
            point2 = endY + (NodePrinter.pointInc *  edgeNum);
           // ++this.numBackward;
        }
        else {
            if( (this.linearNum == endPNode.linearNum - 1) && edgeNum == 0) {//adjacent forward
                point1 = stY - (NodePrinter.pointInc * edgeNum);
                point2 = endY - (NodePrinter.pointInc * edgeNum);
                stX = stX  + this.radius;
                endX = endX - this.radius;
                stY = stY + this.radius;
                endY = endY + this.radius;
                point1 = stY;
                point2 = endY;
               // ++this.numForward;
            }
            else {
                point1 = stY - (NodePrinter.pointInc * edgeNum);
                point2 = endY - (NodePrinter.pointInc * edgeNum);
              // ++this.numForward;
            }
        }

        var col = walkNum % NodePrinter.walkColours.length;

        ctx.strokeStyle = NodePrinter.walkColours[col];

        ctx.moveTo(stX,stY);
        ctx.bezierCurveTo(stX,point1,endX,point2,endX,endY);




    };


}
NodePrinter.walkColours = ['black','blue','red','green','purple','orange'];
NodePrinter.startY = 80;
NodePrinter.startX = 20;
NodePrinter.xInc = 45;
NodePrinter.pointInc = 15;


function executeMisc(e) {

    initialiseResults();

    var grp = document.getElementsByName('miscRadio');

    var val = '';
    for (var i = 0;i < grp.length; ++i) {
        if (grp[i].checked) {
            val = grp[i].value;
        }
    }


    var resEl = document.getElementById('miscQuickResult');
    var par1El = document.getElementById('miscParam1');
    var par2El = document.getElementById('miscParam2');
    var par3El = document.getElementById('miscParam3');

    switch (val) {
        case '1' :
            var rev = reverseComplement(par1El.value);
            resEl.value = rev;
            break;
        case '2' :
            var ar1 = par1El.value.split('\n');
            var ar2 = ar1.map(function(el) {
                return el.trim();
            });
            var profMat = ar2.map(function(el) {
                var spl = el.split(' ');
                var splFloat = spl.map(function(el2) {
                        return parseFloat(el2);

                    }
                );
                return splFloat;

            });
            var dna = par2El.value;
            var k = parseInt(par3El.value);

            var res = profileMostProbable(profMat,dna,k);
            resEl.value = 'Best kmer: \n' + res[0] + '\nBest prob:\n' + res[1];


            break;
        case '3': //Hamilton path
            var k = parseInt(par2El.value);
            par1El.value = cleanContents(par1El.value);
            var dna = par1El.value;

            resEl.value = 'results\n';
            var hamGraph = new DGraph(dna,DGraph.fromDna,DGraph.hamGraph,k,false);

            resEl.value+= hamGraph.sumInfo();
            resEl.value += '\n' + hamGraph.getAdjList();


            var nodes = createHamNodes(dna,null,k);
            var resArray = [];

           // hamiltonCanvas(nodes);

            for (var i = 0;i < parseInt(par3El.value);++i) {
                resArray.push(hamPath(nodes));

            }

            resEl.value+= '\nresults';
            var eqCount = 0;
            var neButSameLengthCount = 0;

            resArray.forEach(function(el,i) {

                hamiltonCanvas(nodes,el[2]);

                if (dna === el[1]) {
                    resEl.value += '\n' + el[1];
                    resEl.value+= ' Equal ' + i;
                    ++eqCount;
                }
                else if (dna.length == el[1].length) {
                    resEl.value += '\n' + el[1];
                    resEl.value+= ' Nope sl ' + i;
                    ++neButSameLengthCount;
                }
                else if (el[1] === '') {
                    resEl.value+= ' \n' + 'First not found';
                }

            });
            if ((eqCount == 0) && (neButSameLengthCount == 0)) {
                resEl.value+= ' \n' + 'None found';
            }
            //var res = hamPath(nodes);
/*
            resEl.value = 'nodes: ';
            resEl.value += '\n' + res[1];
            if (dna === res[1]) {
                resEl.value+= '\nEqual';
         }
            else if (dna.length == res[1].length) {
                resEl.value+= '\nNot Equal but same len';
            }
            else {
                resEl.value+= '\nNot Equal';
            }
            res[0].forEach(function(el) {
                var pred,succ;
                pred = el.predecessors ? el.predecessors.length : 0;
                succ = el.successors ? el.successors.length : 0;
                resEl.value += '\n' + el.dna + ' p: ' + pred + ' s: ' + succ;

            });
            */
            break;

        case '4': //deb cycle from dna
            var k = parseInt(par2El.value);
            par1El.value = cleanContents(par1El.value);
            var dna = par1El.value;
            //var nodes = createDebNodes(dna,k);
            var debGraph = new DGraph(dna,DGraph.fromDna,DGraph.debGraph,k,true);



            var resArray = [];
            for (var i = 0;i < parseInt(par3El.value);++i) {
                //resArray.push(debPath(nodes));
             //   resArray.push(debGraph.debPath());
                debGraph.debCycle();

                debruijnCanvas(debGraph);

                resArray.push([[],debGraph.edgePathReconstructed()]);
            }
            resEl.value = 'results\n';
            var eqCount = 0;
            var neButSameLengthCount = 0;

           // resEl.value += 'rec: ' + debGraph.edgePathReconstructed();

            resArray.forEach(function(el,i) {


            //resArray.forEach(function(el,i) {
                //resEl.value += el[1];
                if (dna === el[1]) {
                    resEl.value += '\n' + el[1];
                    resEl.value+= ' Equal ' + i;
                    ++eqCount;
                }
                else if (dna.length == el[1].length) {
                    var circ = false;


                    for (var ch = 0;ch < dna.length - k;++ch) {
                        if (dna === el[1].substring(ch) + el[1].substring(0,ch)) {
                            resEl.value += '\n' + el[1];
                            resEl.value+= ' Yep circ ' + i;
                            ++neButSameLengthCount;
                            circ = true;
                            break;

                        }
                    }
                    if (!circ) {
                        resEl.value += '\n' + el[1];
                        resEl.value += ' Nope sl ' + i;
                        ++neButSameLengthCount;
                    }
                }
                else if (el[1] === '') {
                    resEl.value+= ' \n' + 'First not found';
                }
                else {
                    resEl.value+= ' \n' + 'x' + el[1];
                }

            });

            if ((eqCount == 0) && (neButSameLengthCount == 0)) {
                resEl.value+= ' \n' + 'None found';
            }
            break;

        case '5':
            var k = parseInt(par2El.value);
            par1El.value = cleanContents(par1El.value);
            var dna = par1El.value;
            var dist = parseInt(par3El.value);
            var res = kmerPairedComposition(dna,k,dist);

            resEl.value = 'results: ';
            res.forEach(function(el) {
                resEl.value += '\n' + el[0] + ':' +  el[1];
            });
            //resEl.value = 'results: ' + res[0][0];

            break;

        case '6': //deb path from paired dna
            var k = parseInt(par2El.value);
            par1El.value = cleanContents(par1El.value);
            var dna = par1El.value;
            var pairDist = parseInt(par3El.value);
            //var nodes = createDebNodes(dna,k);
            var debGraph = new DGraph(dna,DGraph.fromPairedDna,DGraph.debGraph,k,false,pairDist);


            var resArray = [];
            for (var i = 0;i < parseInt(par3El.value);++i) {
                //resArray.push(debPath(nodes));
                //   resArray.push(debGraph.debPath());
                debGraph.debPath();

                debruijnCanvas(debGraph);

                resArray.push([[],debGraph.edgePathReconstructedPairs()]);
            }
            resEl.value = 'results\n';
            var eqCount = 0;
            var neButSameLengthCount = 0;

            // resEl.value += 'rec: ' + debGraph.edgePathReconstructed();

            resArray.forEach(function(el,i) {


                //resArray.forEach(function(el,i) {
                //resEl.value += el[1];
                if (dna === el[1]) {
                    resEl.value += '\n' + el[1];
                    resEl.value+= ' Equal ' + i;
                    ++eqCount;
                }
                else if (dna.length == el[1].length) {
                    var circ = false;


                    for (var ch = 0;ch < dna.length - k;++ch) {
                        if (dna === el[1].substring(ch) + el[1].substring(0,ch)) {
                            resEl.value += '\n' + el[1];
                            resEl.value+= ' Yep circ ' + i;
                            ++neButSameLengthCount;
                            circ = true;
                            break;

                        }
                    }
                    if (!circ) {
                        resEl.value += '\n' + el[1];
                        resEl.value += ' Nope sl ' + i;
                        ++neButSameLengthCount;
                    }
                }
                else if (el[1] === '') {
                    resEl.value+= ' \n' + 'First not found';
                }
                else {
                    resEl.value+= ' \n' + 'x' + el[1];
                }

            });

            if ((eqCount == 0) && (neButSameLengthCount == 0)) {
                resEl.value+= ' \n' + 'None found';
            }
            break;
        case '7': //debruin path from reads
           // var k = parseInt(par2El.value);
           // par1El.value = cleanContents(par1El.value);
           // var dna = par1El.value;
            var reads = par1El.value.split('\n');

            var k = reads[0].length;

            var debGraph = new DGraph(reads,DGraph.fromReads,DGraph.debGraph,k,false);

            resEl.value = 'results\n';

            var resArray = [];
            for (var i = 0;i < parseInt(par3El.value);++i) {
                //resArray.push(debPath(nodes));
                //   resArray.push(debGraph.debPath());
                debGraph.debPath();

                resStr = '';
                debGraph.nodes.forEach(function(node) {
                    resStr += '\n' + node.dna + ' ->  ';
                    node.successors.forEach(function(suc,i) {
                        resStr += suc.targetNode.dna;
                        if (i < node.successors.length - 1) {
                            resStr += ',';
                        }

                    });

                });
                resEl.value+= resStr;

                debruijnCanvas(debGraph);

                resArray.push([[],debGraph.edgePathReconstructed()]);
            }

            var eqCount = 0;
            var neButSameLengthCount = 0;

            // resEl.value += 'rec: ' + debGraph.edgePathReconstructed();

            resArray.forEach(function(el,i) {


                //resArray.forEach(function(el,i) {
                //resEl.value += el[1];
                resEl.value += '\n' + el[1] + ' ' + i;

            });

            break;


        case '8': //debruijn path from dna
            var k = parseInt(par2El.value);
            par1El.value = cleanContents(par1El.value);
            var dna = par1El.value;
            //var nodes = createDebNodes(dna,k);

            var debGraph = new DGraph(dna,DGraph.fromDna,DGraph.debGraph,k,false);

        resEl.value = 'results\n';

        resEl.value+= debGraph.sumInfo();

          var resArray = [];
            for (var i = 0;i < parseInt(par3El.value);++i) {
                //resArray.push(debPath(nodes));
                //   resArray.push(debGraph.debPath());

                debGraph.debPath();

                var resStr = debGraph.getAdjList();
                /*
                resStr = '';
                debGraph.nodes.forEach(function(node) {
                    resStr += '\n' + node.dna + ' ->  ';
                    node.successors.forEach(function(suc,i) {
                        resStr += suc.targetNode.dna;
                        if (i < node.successors.length - 1) {
                            resStr += ',';
                        }

                    });

                });
                */

                resEl.value+= resStr;


                debruijnCanvas(debGraph);


                resArray.push([[],debGraph.edgePathReconstructed()]);
            }



            var eqCount = 0;
            var neButSameLengthCount = 0;

            // resEl.value += 'rec: ' + debGraph.edgePathReconstructed();

            resArray.forEach(function(el,i) {


                //resArray.forEach(function(el,i) {
                //resEl.value += el[1];
                if (dna === el[1]) {
                    resEl.value += '\n' + el[1];
                    resEl.value+= ' Equal ' + i;
                    ++eqCount;
                }
                else if (dna.length == el[1].length) {
                    var circ = false;


                    for (var ch = 0;ch < dna.length - k;++ch) {
                        if (dna === el[1].substring(ch) + el[1].substring(0,ch)) {
                            resEl.value += '\n' + el[1];
                            resEl.value+= ' Yep circ ' + i;
                            ++neButSameLengthCount;
                            circ = true;
                            break;

                        }
                    }
                    if (!circ) {
                        resEl.value += '\n' + el[1];
                        resEl.value += ' Nope sl ' + i;
                        ++neButSameLengthCount;
                    }
                }
                else if (el[1] === '') {
                    resEl.value+= ' \n' + 'First not found';
                }
                else {
                    resEl.value+= ' \n' + 'x' + el[1];
                }

            });

            if ((eqCount == 0) && (neButSameLengthCount == 0)) {
                resEl.value+= ' \n' + 'None found';
            }
            break;
        case '9': //deb path from paired reads
            var k = parseInt(par2El.value);
            //par1El.value = cleanContents(par1El.value);
           // var dna = par1El.value;
            var reads = par1El.value.split('\n');

            reads = reads.map(function(el) {
                var pair = el.split('|');

                return pair;
            });

            var pairDist = parseInt(par3El.value);
            //var nodes = createDebNodes(dna,k);
            var debGraph = new DGraph(reads,DGraph.fromPairedReads,DGraph.debGraph,k,false,pairDist);


            var resArray = [];
            for (var i = 0;i < 1;++i) {
                //resArray.push(debPath(nodes));
                //   resArray.push(debGraph.debPath());
                debGraph.debPath();

                debruijnCanvas(debGraph);

                resArray.push([[],debGraph.edgePathReconstructedPairs()]);
            }
            resEl.value = 'results\n';
            var eqCount = 0;
            var neButSameLengthCount = 0;

            // resEl.value += 'rec: ' + debGraph.edgePathReconstructed();

            resArray.forEach(function(el,i) {


                //resArray.forEach(function(el,i) {
                //resEl.value += el[1];
                resEl.value += '\n' + el[1] + ' ' + i;
            });

                break;
        case '10': //pattern to number

            par1El.value = cleanContents(par1El.value);
            // var dna = par1El.value;
            resEl.value = 'results\n';
            resEl.value+= kMerToInd(par1El.value);

            break;

        case '11': //number to pattern

            var ind = par1El.value; //string, because can't handle really big ints
            var k = parseInt(par2El.value);
            // var dna =         par1El.value;
            resEl.value = 'results\n';
            resEl.value+= indToKmer(ind,k);

            break;

        case '12': //d-neighbourhood
            par1El.value = cleanContents(par1El.value);

            var d = parseInt(par2El.value);

            // var dna =         par1El.value;
            resEl.value = 'results\n';
            var res = kMersWithMaxDist(par1El.value,d);
            var resStr = '';
            res.forEach(function(el) {
                resStr+=el + '\n';
            });

            resEl.value+= resStr;

            break;

        case '13': //dist pattern strings
            par1El.value = cleanContents(par1El.value);

            var strs = par2El.value.split(' ');

            // var dna =         par1El.value;
            resEl.value = 'results\n';
            var res = patternSequencesDist(par1El.value,strs);

            resEl.value+= res;

            break;




        case '14': //kmer composition
            var k = parseInt(par2El.value);
            par1El.value = cleanContents(par1El.value);
            var dna = par1El.value;
            var dist = parseInt(par3El.value);
            var res = kmerComposition(dna,k);

            resEl.value = 'results: ';
            var resStr = '';
            res.forEach(function(el) {
                resStr += '\n' + el;
            });
            resEl.value += resStr;
            //resEl.value = 'results: ' + res[0][0];

            break;

        case '15': //str from genome path
            var path = par1El.value.split('\n');

            var res = stringFromGenomePath(path);

            resEl.value = 'results: ';

            resEl.value += res;
            //resEl.value = 'results: ' + res[0][0];

            break;

        case '16': //overlap graph

            //par1El.value = cleanContents(par1El.value);
            var reads = par1El.value.split('\n');
            var k = reads[0].length;
            var nodes = createHamNodes(null,reads,k);
            var resArray = [];

            // hamiltonCanvas(nodes);

            var res = hamPath(nodes);

            var resStr = '';

            res[2].forEach(function(node) {
                resStr+='\n'  + node.dna;
                resStr+= ' -> ';
                node.successors.forEach(function(suc) {
                    resStr += suc.dna + ' ';

                });

            });
            resEl.value = 'results';
            resEl.value+= resStr;
            var eqCount = 0;
            var neButSameLengthCount = 0;

            var el = res;

            if (nodes.length < 30) {
                hamiltonCanvas(nodes, el[2]);
            }

            break;
        case '17':
            var k = parseInt(par2El.value);
            //par1El.value = cleanContents(par1El.value);
            // var dna = par1El.value;
            var adjList = par1El.value.split('\n');

            //var nodes = createDebNodes(dna,k);
            var debGraph = new DGraph(adjList,DGraph.fromAdjList,DGraph.debGraph,k,false);
            resEl.value = 'results\n';

            var resArray = [];
            for (var i = 0;i < 1;++i) {
                //resArray.push(debPath(nodes));
                //   resArray.push(debGraph.debPath());
                debGraph.debCycle();

                resEl.value += debGraph.edgePathToText();

                resStr = '';
                /*
                debGraph.nodes.forEach(function(node) {
                    resStr += '\n' + node.dna + ' ->  ';
                    node.successors.forEach(function(suc,i) {
                        resStr += suc.targetNode.dna;
                        if (i < node.successors.length - 1) {
                            resStr += ',';
                        }

                    });

                });
                resEl.value+= resStr;
                */

                debruijnCanvas(debGraph);

                resArray.push([[],debGraph.edgePathReconstructedPairs()]);
            }

            var eqCount = 0;
            var neButSameLengthCount = 0;

            // resEl.value += 'rec: ' + debGraph.edgePathReconstructed();

            resArray.forEach(function(el,i) {


                //resArray.forEach(function(el,i) {
                //resEl.value += el[1];
               // resEl.value += '\n' + el[1] + ' ' + i;
            });

            break;

        case '18':
            var k = parseInt(par2El.value);
            var adjList = par1El.value.split('\n');

            var debGraph = new DGraph(adjList,DGraph.fromAdjList,DGraph.debGraph,k,false);
            resEl.value = 'results\n';

            var resArray = [];
            for (var i = 0;i < 1;++i) {
               debGraph.debPath();
                resEl.value += debGraph.edgePathToText();

                debruijnCanvas(debGraph);

                resArray.push([[],debGraph.edgePathReconstructedPairs()]);
            }


            break;

        case '19':
            var k = parseInt(par1El.value);

            resEl.value = 'results\n';

            var pad = "00000000000000";

            var reads = [];
            for (var i = 0;i < Math.pow(2,k);++i) {
                var binStr = (i >>> 0).toString(2);
                binStr = pad + binStr;
                binStr = binStr.substr(binStr.length - k);
                reads.push(binStr);

               // resEl.value +=  (i >>> 0).toString(2) + ' ';
                resEl.value+= binStr + '\n';
           }


            var debGraph = new DGraph(reads,DGraph.fromReads,DGraph.debGraph,k,false);

            resEl.value+= 'results\n';

            var resArray = [];
            lim = parseInt(par3El.value);
            if (lim == 0) {
                lim = 1;
            }

            for (var i = 0;i < parseInt(par3El.value);++i) {
                //resArray.push(debPath(nodes));
                //   resArray.push(debGraph.debPath());
                debGraph.debPath();

                resStr = '';
                debGraph.nodes.forEach(function(node) {
                    resStr += '\n' + node.dna + ' ->  ';
                    node.successors.forEach(function(suc,i) {
                        resStr += suc.targetNode.dna;
                        if (i < node.successors.length - 1) {
                            resStr += ',';
                        }

                    });

                });
                resEl.value+= resStr;

                debruijnCanvas(debGraph);

                resArray.push([[],debGraph.edgePathReconstructed()]);
            }

            var eqCount = 0;
            var neButSameLengthCount = 0;

            // resEl.value += 'rec: ' + debGraph.edgePathReconstructed();

            resArray.forEach(function(el,i) {


                //resArray.forEach(function(el,i) {
                //resEl.value += el[1];
                resEl.value += '\n' + el[1] + ' ' + i;

            });

            break;

        case '20': //Genome path to string

            var path = par1El.value.split('\n');

            var str = '';
            path.forEach(function(el,i) {
                if (i == 0) {
                    str+=el;

                }
                else {
                    str+=el.substring(el.length - 1);
                }

            });
            resEl.value = 'results\n';

            resEl.value+=str;


            break;

        case '21': //DNA to amino string

            var dna = new DNA(par1El.value);

            var rna = new RNA(dna.rnaTranscript());

            var pep = rna.translate(0,true);

            resEl.value = 'results\n';

            resEl.value = 'len: ' + dna.dna.length + '\n';

            resEl.value +=rna.rna;

            resEl.value +='\n' + pep.toShortString('');
            break;

        case '22': //RNA to amino string

            var rna = new RNA(par1El.value);

            var pep = rna.translate(0,true);

            resEl.value = 'results\n';

            resEl.value = 'len: ' + dna.dna.length + '\n';

            resEl.value +=pep.toShortString('');
            break;

        case '23': //Find Amino Acid Sequence in genome

            var dna = par1El.value;
            var targetSeq = par2El.value;
            var targetDNALen = targetSeq.length * Codon.len;

            var seqs = [];



            for (var i = 0;i < dna.length - targetDNALen; ++i) {

                var forwardFound = false;
                var testDNA = new DNA(dna.substring(i,i + targetDNALen));
                var testRNA = new RNA(testDNA.rnaTranscript());
                var testProt = testRNA.translate(0,true);
                if (testProt.toShortString('') == targetSeq) {
                    forwardFound = true;
                }
                var reverseFound = false;
                var testRevDNA = testDNA.reverseComplement();
                var tst = testRevDNA.rnaTranscript();
                testRevRNA = new RNA(testRevDNA.rnaTranscript());
                testRevProt = testRevRNA.translate(0,true);
                if (testRevProt.toShortString('') == targetSeq) {
                    reverseFound = true;
                }

                if (forwardFound || reverseFound) {
                    seqs.push(testDNA.dna);
                }

            }


            resEl.value = 'results\n';

            resEl.value = 'len: ' + dna.length + '\n';

            seqs.forEach(function(seq) {
                resEl.value +=seq + '\n';
            });

            break;

        case '24': //Spectrum


            var pep = new Peptide(Peptide.AminoArrFromStr(par1El.value));
            var cyc = (par2El.value == 1 ? true : false);

            resEl.value = 'results\n';

            resEl.value += 'Peptide len: ' + pep.peptide.length + '\n';



            var spec = pep.spectrum(cyc);

            resEl.value += 'Spectrum len: ' + spec.length + '\n';

            spec.forEach(function(w) {
                resEl.value += w + ' ';

            });


            break

        case '25': //Num peptides with mass


            var mass = parseInt(par1El.value);


            resEl.value = 'results\n';

            var res = allPossiblePeptidesWithWeight(mass);

            resEl.value +=  res;
            break;

        case '26': //Cyclopeptide sequencing


            var experimentalSpectrum = par1El.value;


            resEl.value = 'results\n';

            var res = cycloPeptideSequencing(experimentalSpectrum);

            resEl.value +=  res;

            break;

        default:
            break;
    }

/*
     mfMotif = [
        [
            [

            ],
            res[0]
        ]
    ];
    dnaMasterStrings.forEach(function(el) {
        mfMotif[0][0].push([[0],[]]);

    });
    colourDNA(dnaMaster,null,false);
*/

}

function initialiseMotifParams() {

    var paramObj = {};

    var dnaStrings = document.getElementById('dnaStrings').value.split('\n');

    dnaStrings = dnaStrings.filter(function(el) {
        if (el.length == 0) {
            return false;
        }
        else {
            return true;
        }
    });
    dnaStrings = dnaStrings.map(function(el) {
        return el.trim().toUpperCase();

    });
    dnaMasterStrings =  dnaStrings;//document.getElementById('dnaStrings').value.split('\n');

    dnaMaster = document.getElementById('dnaStrings').value.replace(/ /g,'').toUpperCase();

    dnaMasterChanged();

    var numKmer = document.getElementById('numKmers');
    numKmer.value = 1;

    paramObj.k = parseInt(document.getElementById('kMerLenMS').value);

    // var dna = dnaMaster; //document.getElementById('dnaInput').value;

    paramObj.laplace = document.getElementById('laplaceMS').checked;

    paramObj.maxMismatch = parseInt(document.getElementById('maxMismatchMS').value);

    paramObj.numIters = parseInt(document.getElementById('numItersMS').value);

    return paramObj;



}


function runMotif(e) {

    var rads = document.getElementsByName('motifMethod');

    var selected;

    for (var i = 0;i < rads.length;++i) {
        if (rads[i].checked) {
            selected = rads[i];
            break;

        }
    }

    switch (selected.id) {

    case 'motifKD':
        motifSearch(e);
        break;
    case 'motifBrute':
        bruteMotifSearch(e);
        break;
    case 'motifMedian':
        medianString(e);
        break;
    case 'motifGreedy':
        greedyMotif(e);
        break;
    case 'motifRandom':
        randomMotif(e);
        break;
    case 'motifGibbs':
        gibbsSampler(e);
         break;

    default:
            break;

    }
}

function motifCompare(e) {
    initialiseResults();

    dnaMasterStrings =  document.getElementById('dnaStrings').value.split('\n');
    dnaMaster = document.getElementById('dnaStrings').value;
    dnaMasterChanged();

    var laplace = document.getElementById('laplaceMS').checked;

    var res = scoreMotifs(dnaMasterStrings,laplace);
    mfMotif = [
                [
                    [

                    ],
                    res[0]
                ]
              ];
    dnaMasterStrings.forEach(function(el) {
        mfMotif[0][0].push([[0],[]]);

    });
    colourDNA(dnaMaster,null,false);

    if (document.getElementById('debugMS').checked) {

        document.getElementById("debugText").value = 'Consensus: ' + res[0];
        document.getElementById("debugText").value += '\nScore: ' +  res[1];
        document.getElementById("debugText").value += '\nEntropy: '  +  res[2];
        var countMatrixStr = '';
       res[3].forEach(function(el) {
            countMatrixStr += '\n' + el;
        });
        var entropyMatrixStr = '';
        res[4].forEach(function(el) {
            entropyMatrixStr += '\n' + el.map(function(e) {
                    return e.toFixed(3);

                });
        });
        var profileMatrixStr = '';
        res[5].forEach(function(el) {
            profileMatrixStr += '\n' + el.map(function(e) {
                    return e.toFixed(3);

                });
        });

        document.getElementById("debugText").value += '\nCount Matrix: '  +  countMatrixStr;
        document.getElementById("debugText").value += '\nEntropy Matrix: '  +  entropyMatrixStr;
        document.getElementById("debugText").value += '\nProfile Matrix: '  +  profileMatrixStr;

        var extraInfoAr = [];

        mfMotif[0][0].forEach(function(el,i) {
            var ind = el[0][0];
            var km = dnaMasterStrings[i].substring(ind,ind+k);

            extraInfoAr.push(km);
        });

        logoCanvas(extraInfoAr);

    }

    document.getElementById("moreDetailText").value = 'Consensus: ' + res[0];
    document.getElementById("moreDetailText").value += '\nScore: ' +  res[1];
    document.getElementById("moreDetailText").value += '\nEntropy: '  +  res[2];
}

function motifSearch(e) {

    initialiseResults();

    dnaMasterStrings =  document.getElementById('dnaStrings').value.split('\n');
    dnaMaster = document.getElementById('dnaStrings').value;
    dnaMasterChanged();

    var numKmer = document.getElementById('numKmers');
    numKmer.value = 1;

    var k = parseInt(document.getElementById('kMerLenMS').value);

   // var dna = dnaMaster; //document.getElementById('dnaInput').value;

    var inclRevCompl = false; //document.getElementById('includeRevComplMS').checked;

    var maxMismatch = parseInt(document.getElementById('maxMismatchMS').value);

    w.postMessage({'task' : 'motifSearch', 'dnaStrings' : dnaMasterStrings, 'k':k, 'inclRevCompl': inclRevCompl,'maxMismatch': maxMismatch}); // Start the worker.

    /*
     setTimeout(function() {
     var indexArray = countKMer(kmer,dna,document.getElementById('includeRevComplKS').checked);
     //alert ('ind: ' + indexArray);
     console.log('normal finished');
     mfk = [indexArray];
     k = kmer.length;
     colourDNA(dna);

     },0);
     */
}


function medianString(e) {

    initialiseResults();

    dnaMasterStrings =  document.getElementById('dnaStrings').value.split('\n');
    dnaMaster = document.getElementById('dnaStrings').value;
    dnaMasterChanged();

    var numKmer = document.getElementById('numKmers');
    numKmer.value = 1;

    var k = parseInt(document.getElementById('kMerLenMS').value);

    // var dna = dnaMaster; //document.getElementById('dnaInput').value;

    var inclRevCompl = false//document.getElementById('includeRevComplMS').checked;

    var maxMismatch = parseInt(document.getElementById('maxMismatchMS').value);

    w.postMessage({'task' : 'medianString', 'dnaStrings' : dnaMasterStrings, 'k':k, 'inclRevCompl': inclRevCompl,'maxMismatch': maxMismatch}); // Start the worker.

    /*
     setTimeout(function() {
     var indexArray = countKMer(kmer,dna,document.getElementById('includeRevComplKS').checked);
     //alert ('ind: ' + indexArray);
     console.log('normal finished');
     mfk = [indexArray];
     k = kmer.length;
     colourDNA(dna);

     },0);
     */
}

function greedyMotif(e) {

    initialiseResults();


    var dnaStrings = document.getElementById('dnaStrings').value.split('\n');

    dnaStrings = dnaStrings.filter(function(el) {
        if (el.length == 0) {
            return false;
        }
        else {
            return true;
        }
    });
    dnaStrings = dnaStrings.map(function(el) {
        return el.trim().toUpperCase();

    });
    dnaMasterStrings =  dnaStrings;//document.getElementById('dnaStrings').value.split('\n');

    dnaMaster = document.getElementById('dnaStrings').value.replace(/ /g,'').toUpperCase();

    dnaMasterChanged();

    var numKmer = document.getElementById('numKmers');
    numKmer.value = 1;

    var k = parseInt(document.getElementById('kMerLenMS').value);

    // var dna = dnaMaster; //document.getElementById('dnaInput').value;

    var laplace = document.getElementById('laplaceMS').checked;

    var maxMismatch = parseInt(document.getElementById('maxMismatchMS').value);

    w.postMessage({'task' : 'greedyMotif', 'dnaStrings' : dnaMasterStrings, 'k':k, 'laplace': laplace,'maxMismatch': maxMismatch}); // Start the worker.

    /*
     setTimeout(function() {
     var indexArray = countKMer(kmer,dna,document.getElementById('includeRevComplKS').checked);
     //alert ('ind: ' + indexArray);
     console.log('normal finished');
     mfk = [indexArray];
     k = kmer.length;
     colourDNA(dna);

     },0);
     */
}


function randomMotif(e) {

    initialiseResults();

    var paramObj = initialiseMotifParams();

    w.postMessage({
        'task': 'randomMotif',
        'dnaStrings': dnaMasterStrings,
        'k': paramObj.k,
        'laplace': paramObj.laplace,
        'numIters': paramObj.numIters,
        'maxMismatch': paramObj.maxMismatch
    }); // Start the worker.

}

function gibbsSampler(e) {

    initialiseResults();

    var paramObj = initialiseMotifParams();

    w.postMessage({
        'task': 'gibbsSampler',
        'dnaStrings': dnaMasterStrings,
        'k': paramObj.k,
        'laplace': paramObj.laplace,
        'numIters': paramObj.numIters,
        'maxMismatch': paramObj.maxMismatch
    }); // Start the worker.

}

function bruteMotifSearch(e) {

    initialiseResults();

    var paramObj = initialiseMotifParams();

    w.postMessage({
        'task': 'bruteMotifSearch',
        'dnaStrings': dnaMasterStrings,
        'k': paramObj.k,
        'laplace': paramObj.laplace,
        'numIters': paramObj.numIters,
        'maxMismatch': paramObj.maxMismatch
    }); // Start the worker.

}

function initialiseSequencingParams() {

    var paramObj = {};

    //clean contents in case apply not pressed yet?
    //....
    //


    //dnaMasterStrings =  dnaStrings;//document.getElementById('dnaStrings').value.split('\n');
    //dnaMaster = document.getElementById('dnaStrings').value.replace(/ /g,'').toUpperCase();
    //dnaMasterChanged();

    var numKmer = document.getElementById('numKmers');
    numKmer.value = 1;

    paramObj.k = parseInt(document.getElementById('kmerLenSA').value);

    paramObj.numIters = parseInt(document.getElementById('numItersSA').value);

    paramObj.debug = document.getElementById('debugSA').checked;

    var rads = document.getElementsByName('seqMethod');

    var methodSelected;

    for (var i = 0;i < rads.length;++i) {
        if (rads[i].checked) {
            methodSelected = rads[i];
            break;

        }
    }

    switch (methodSelected.id) {
        case 'seqComp':
            paramObj.seqMethod = DGraph.seqTypeComp;
            break;
        case 'seqPath':
            paramObj.seqMethod = DGraph.seqTypePath;
            break;
        case 'seqCycle':
            paramObj.seqMethod = DGraph.seqTypeCycle;
            break;
        default:
            break;


    }
    //paramObj.seqMethod = methodSelected.id;

    rads = document.getElementsByName('seqInput');

    var inputSelected;

    for (var i = 0;i < rads.length;++i) {
        if (rads[i].checked) {
            inputSelected = rads[i];
            break;

        }
    }

    switch (inputSelected.id) {
        case 'seqDNA':
            paramObj.seqInput = DGraph.fromDna;
            break;
        case 'seqReads':
            paramObj.seqInput = DGraph.fromReads;
            break;
        case 'seqAdjList':
            paramObj.seqInput = DGraph.fromAdjList;
            break;
        case 'seqPairedDNA':
            paramObj.seqInput = DGraph.fromPairedDna;
            paramObj.pairDist = parseInt(document.getElementById('pairDistSA').value);
            break;
        case 'seqPairedReads':
            paramObj.seqInput = DGraph.fromPairedReads;
            paramObj.pairDist = parseInt(document.getElementById('pairDistSA').value);
            break;
        default:
            break;


    }

    ;


    // paramObj.seqInput = inputSelected.id;


    rads =   document.getElementsByName('seqGraphType');

    var typeSelected;

    for (var i = 0;i < rads.length;++i) {
        if (rads[i].checked) {
            typeSelected = rads[i];
            break;

        }
    }

    switch (typeSelected.id) {
        case 'seqGraphTypeHam':
            paramObj.seqType = DGraph.hamGraph;
            break;
        default:
            paramObj.seqType = DGraph.debGraph;
            break;
    }

   //paramObj.seqType = typeSelected;

    return paramObj;

}

function readsToMfk(grph) {
    var prevKmer = '';
    var currInds = [];

    mfk = [];

    grph.reads.forEach(function(el,i) {

        if (i == 0) {
            prevKmer = el;
        }

        if (el === prevKmer) {
            currInds.push(grph.readPosns[i]);
        }
        else {
            mfk.push([currInds,[],prevKmer]);
            currInds = [grph.readPosns[i]];
            prevKmer = el;
        }

        if (i == grph.reads.length - 1) {
            mfk.push([currInds,[],el]);
        }

    });

    /*
     mfk = res.map(function(el,i) {
     return [[el[1]],[],el[0]];
     });
     */

    numKmer = document.getElementById('numKmers');
    numKmer.max = mfk.length;
    if (mfk.length == 0) {
        numKmer.min = 0;
        numKmer.value = 0;
    }
    else {
        numKmer.min = 1;
        numKmer.value = 1;
    }

    numKmerVal = parseInt(numKmer.value);
}


function runSequencing(e) {

    sequencingInput();

    initialiseResults();

    paramObj = initialiseSequencingParams();

   // var rads = document.getElementsByName('seqMethod');

    /*
    var methodSelected;

    for (var i = 0;i < rads.length;++i) {
        if (rads[i].checked) {
            methodSelected = rads[i];
            break;

        }
    }
    */

    var reads, input;

    var grph;

    var k;

    if (( paramObj.seqInput == DGraph.fromDna) || (paramObj.seqInput == DGraph.fromPairedDna)) {
        k = paramObj.k;
        input = dnaMaster;

    }
    else if ( paramObj.seqInput == DGraph.fromReads)   {

        reads = dnaMasterStrings.split('\n');

        input = reads;

        if (reads[0]) {
            k = reads[0].length;
        }
        else {
            k = 1;
        }
    }
    else if ( paramObj.seqInput == DGraph.fromPairedReads)   {

        reads = dnaMasterStrings.split('\n');


        if (reads[0]) {
            k = reads[0].split('|')[0].length;
        }
        else {
            k = 1;
        }
        reads = reads.map(function(el) {
            var spl = el.split('|');
            return [spl[0],spl[1]];

        });

        input = reads;

    }
    else if (paramObj.seqInput == DGraph.fromAdjList) {
        reads = dnaMasterStrings.split('\n');

        input = reads;

        if (reads[0]) {
            k = reads[0].split(' -> ')[0].length;
        }
        else {
            k = 1;
        }

    }

    grph = new DGraph(input,paramObj.seqInput,paramObj.seqType,k,false,paramObj.pairDist);


    switch (paramObj.seqMethod) {

        case DGraph.seqTypeComp:
            break;

        case DGraph.seqTypePath:

            if (paramObj.seqType == DGraph.hamGraph) {
                grph.hamPath();
            }
            else {
                grph.debPath();
            }
            break;
        case DGraph.seqTypeCycle:
            if (paramObj.seqType == DGraph.hamGraph) {
               // grph.hamCycle();
            }
            else {
                grph.debCycle();
            }
            break;
        default:
            break;



    }



    graphChanged(grph);

    if (paramObj.seqInput == DGraph.fromDna) {

        readsToMfk(grph);
        dna = dnaMaster; // document.getElementById('dnaInput').value;
        colourDNA(dna, null, false);
    }
    else {

    }



    if (document.getElementById('debugSA').checked) {
        var resStr = '';
        grph.reads.forEach(function (el) {
            resStr += '\n' + el;
        });
        document.getElementById('debugText').value = resStr;
    }



/*
    switch (paramObj.seqMethod) {

        case DGraph.seqTypeComp:
            kmerComp(e);
            break;

        case DGraph.seqTypePath:

            findPath(e);
            break;

        case DGraph.seqTypeCycle:
            //pathFromReads(e);
            break;

        default:
            break;

    }

*/

}

function kmerComp(e) {


    var grph = new DGraph(dnaMaster,DGraph.fromDna,DGraph.hamGraph,paramObj.k,false);

    //grph.hamPath();

    graphChanged(grph);

    readsToMfk(grph);

    dna = dnaMaster; // document.getElementById('dnaInput').value;
    colourDNA(dna, null, false);

    if (document.getElementById('debugSA').checked) {
        var resStr = '';
        grph.reads.forEach(function (el) {
            resStr += '\n' + el;
        });
        document.getElementById('debugText').value = resStr;
    }

    /*
    w.postMessage({
        'task': 'randomMotif',
        'dnaStrings': dnaMasterStrings,
        'k': paramObj.k,
        'laplace': paramObj.laplace,
        'numIters': paramObj.numIters,
        'maxMismatch': paramObj.maxMismatch
    }); // Start the worker.
    */

}

function findPath(e) {

    var reads, input;

    var grph;

    var k;

    if (( paramObj.seqInput == DGraph.fromDna) || (paramObj.seqInput == DGraph.fromPairedDna)) {
        k = paramObj.k;
        input = dnaMaster;

    }
    else if ( paramObj.seqInput == DGraph.fromReads)   {

        reads = dnaMasterStrings.split('\n');

        input = reads;

        if (reads[0]) {
            k = reads[0].length;
        }
        else {
            k = 1;
        }
    }
    else if ( paramObj.seqInput == DGraph.fromPairedReads)   {

        reads = dnaMasterStrings.split('\n');

        input = reads;

        if (reads[0]) {
            k = reads[0].split('|')[0].length;
        }
        else {
            k = 1;
        }
    }
    else if (paramObj.seqInput == DGraph.fromAdjList) {
        reads = dnaMasterStrings.split('\n');

        input = reads;

        if (reads[0]) {
            k = reads[0].split(' -> ')[0].length;
        }
        else {
            k = 1;
        }

    }

     grph = new DGraph(input,paramObj.seqInput,paramObj.seqType,k,false,paramObj.pairDist);


    if (paramObj.seqType  == DGraph.hamGraph) {
        grph.hamPath();
    }
    else {
        grph.debPath();
    }


    graphChanged(grph);

    if (paramObj.seqInput == DGraph.fromDna) {

         readsToMfk(grph);
         dna = dnaMaster; // document.getElementById('dnaInput').value;
          colourDNA(dna, null, false);
    }



    if (document.getElementById('debugSA').checked) {
        var resStr = '';
        grph.reads.forEach(function (el) {
            resStr += '\n' + el;
        });
        document.getElementById('debugText').value = resStr;
    }



}

function pathFromReads(e) {

    initialiseResults();

    var paramObj = initialiseSequencingParams();

    var grph;

    var reads = dnaMasterStrings.split('\n');

    var k;
    if (reads[0]) {
        k = reads[0].length;
    }
    else {
        k = 1;
    }

    if (document.getElementById('seqGraphTypeHam').checked) {
        grph = new DGraph(reads,DGraph.fromReads,DGraph.hamGraph,k,false);

        grph.hamPath();

    }
    else {
        grph = new DGraph(reads,DGraph.fromReads,DGraph.debGraph,k,false);

        grph.debPath();
    }


    graphChanged(grph);

    document.getElementById('moreDetailText').value = grph.adjList().join('\n');

    readsToMfk(grph);



    dna = dnaMaster; // document.getElementById('dnaInput').value;
    colourDNA(dna, null, false);

    if (document.getElementById('debugSA').checked) {
        var resStr = '';
        grph.reads.forEach(function (el) {
            resStr += '\n' + el;
        });
        document.getElementById('debugText').value = resStr;
    }



}





//display routines

function collectStats(dna,sk,bc) {

    sk = sk || null;

    bc = bc || null;

    if (bc) {

    }
    else {
        bc = baseCount(dna);
    }

    var stats = 'Stats<br>DNA seq length: ' + dna.length + '<br>';

    var percStats = bc.map(function(el) {
        return [el,el * 1.0 /dna.length * 100.0];

    });

    var gcPerc = gcCount(dna) * 1.0 / dna.length * 100;

    stats+= 'A count: ' + bc[0] +  ' ' + percStats[0][1].toFixed(2) +  '%<br>';
    stats+= 'C count: ' + bc[1] + ' ' + percStats[1][1].toFixed(2) +  '%<br>';
    stats+= 'G count: ' + bc[2] + ' ' + percStats[2][1].toFixed(2) +  '%<br>';
    stats+= 'T count: ' + bc[3] + ' ' + percStats[3][1].toFixed(2) +  '%<br>';
    stats+= 'G-C %: ' + gcPerc.toFixed(2) + '%<br>';

    if (sk) {
        stats+= 'G-C skew min: ' + sk[0] + '[' + sk[3] + ']' +  '<br>';
        stats+= 'G-C skew max: ' + sk[1] +  '[' + sk[4] + ']' +  '<br>';


    }
    else {
        /*
        var basesDone = 0;
        while (basesDone < dna.length) {
            sk = gcSkew(dna, 10, sk);
            basesDone = sk[2].length - 1;
        }
        */
        stats+= 'G-C skew min: ' + '-' + '<br>';
        stats+= 'G-C skew max: ' + '-' + '<br>';



        //var sk = gcSkew(inDNA.value,inDNA.value.length);
    }

    return stats;


}

function initGraphCanvas() {

    var c = document.getElementById("moreCanvas");

    var h = c.height;
    var w = c.width;

    c.setAttribute('width', '660');
    c.setAttribute('height', '150');

    var ctx = c.getContext("2d");

    ctx.clearRect(0,0,w,h);
    ctx.fillStyle = ('#FFFFFF');
    ctx.fillRect(0,0,w,h);
    ctx.fillStyle = ('#000000');

    return ctx;

}

function graphCanvas(grph) {

    var ctx = initGraphCanvas();

    ctx.beginPath();

    var linNum = 1;

    //Build and Print nodes

    var nodePrinters = [];
    grph.nodes.forEach(function(node,i) {
        var nodeAlreadyThere = false;
        for (var j = 0;j < nodePrinters.length;++j) {
            if (nodePrinters[j].node == node) {
                nodeAlreadyThere = true;
                break;
            }
        }
        if (nodeAlreadyThere) {

        }
        else {
            var nodeP = new NodePrinter(node, NodePrinter.startX + ((linNum - 1) * NodePrinter.xInc), NodePrinter.startY,linNum);
            ++linNum;

            nodePrinters.push(nodeP);

        }

    });



    nodePrinters.forEach(function(nodeP) {
        nodeP.print(ctx);

    });

    ctx.stroke();


    //Create  edge printers


    ctx.lineWidth=1;

    ctx.beginPath();
    var stNode,endNode;

    stNode = null;
    endNode = null;
    for (var j = 0;j < nodePrinters.length; ++j) {
        var backCount = 0;
        var forwardCount = 0;

        stNode = nodePrinters[j];



        nodePrinters[j].node.successors.forEach(function (suc) {
            endNode = null;
            for (var i = 0; i < nodePrinters.length; ++i) {
                if (suc.targetNode == nodePrinters[i].node) {
                    endNode = nodePrinters[i];
                    break;
                }
            }

            if (endNode) {

                if (stNode.linearNum > endNode.linearNum) {
                   // stNode.numBackward = backCount;
                    ++backCount;
                }
                else  if ((stNode.linearNum == endNode.linearNum - 1) && (forwardCount == 0)) {
                   // stNode.numForward = forwardCount - 1;
                   // ++forwardCount;

                }
                else {
                   // stNode.numForward = forwardCount;
                    ++forwardCount;

                }

                nodePrinters[j].edgePrinters.push(new EdgePrinter(suc,stNode,endNode,backCount,forwardCount));
                if ((stNode.linearNum == endNode.linearNum - 1) && (forwardCount == 0)) {
                    ++forwardCount;
                }

               // stNode.join(ctx, endNode, 0);

               // ctx.stroke();
            }


        });


    }

    //Print edges

    for (var j = 0;j < nodePrinters.length; ++j) {
            nodePrinters[j].edgePrinters.forEach(function(edgeP) {
                edgeP.print(ctx,0);

            });

    }


    //test
    ctx.stroke();
    ctx.closePath();

    var prevWalkNum = 0;

    ctx.beginPath();
    if (grph.edgePath.length == 0) {

    }
    else {
        grph.edgePath.forEach(function(edge) {
            var srcNode = edge.sourceNode;
            var srcNodeP = null;
            var edgeP = null;
            for (var j = 0; j < nodePrinters.length; ++j) {
                if (nodePrinters[j].node == srcNode) {
                    srcNodeP = nodePrinters[j];
                    break;
                }
            }
            for (var k = 0; k < srcNodeP.edgePrinters.length; ++k) {
                if (srcNodeP.edgePrinters[k].edge == edge) {
                    edgeP = srcNodeP.edgePrinters[k];
                    break;
                }
            }
            if (edgeP) {
                if (edge.walkNum != prevWalkNum) {
                    ctx.stroke();
                     ctx.closePath();
                     ctx.beginPath();
                    prevWalkNum = edge.walkNum;

                }
                edgeP.print(ctx,edge.walkNum);

            }

        });



    }

   // nodePrinters[4].edgePrinters[0].print(ctx,5);



    ctx.stroke();
   // ctx.closePath();
/*

    ctx.beginPath();
    var p = nodePrinters[0];
    p.path(ctx,nodePrinters[2],4,1);

    ctx.stroke();
*/

    ctx.closePath();



}


function hamiltonCanvas(hamNodes,hamPathNodes) {

    var ctx = initGraphCanvas();

    ctx.beginPath();

    var linNum = 1;

    //Build and Print nodes

    var nodePrinters = [];
    hamNodes.forEach(function(hamNode,i) {
        var nodeAlreadyThere = false;
        for (var j = 0;j < nodePrinters.length;++j) {
            if (nodePrinters[j].node == hamNode) {
                nodeAlreadyThere = true;
                break;
            }
        }
        if (nodeAlreadyThere) {

        }
        else {
            var nodeP = new NodePrinter(hamNode, NodePrinter.startX + ((linNum - 1) * NodePrinter.xInc), NodePrinter.startY,linNum);
            ++linNum;

            nodePrinters.push(nodeP);

        }

    });



    nodePrinters.forEach(function(nodeP) {
        nodeP.print(ctx);

    });

    ctx.stroke();


    //Print path

    var backCount = 1;
    var forwardCount = 1;

    for (var i = 0;i < hamPathNodes.length - 1;++i) {
        ctx.beginPath();

        var stNode,endNode;

        stNode = null;
        endNode = null;
        for (var j = 0;j < nodePrinters.length; ++j) {
            if (nodePrinters[j].node == hamPathNodes[i]) {
                stNode = nodePrinters[j];

            }
            if (nodePrinters[j].node == hamPathNodes[i+1]) {
                endNode = nodePrinters[j];

            }

        }
        if (stNode.linearNum > endNode.linearNum) {
            stNode.numBackward = backCount;
            ++backCount;
        }
        else {
            stNode.numForward = forwardCount;
            ++forwardCount;

        }

        stNode.join(ctx,endNode,1);

        ctx.stroke();


    }




    ctx.stroke();
    ctx.closePath();



}

function debruijnCanvas(deb) {

   var ctx = initGraphCanvas();

    ctx.beginPath();

    var linNum = 1;

    //Build and Print nodes

    var nodePrinters = [];
    deb.edgePath.forEach(function(edge,i) {
        var nodeAlreadyThere = false;
        for (var j = 0;j < nodePrinters.length;++j) {
            if (nodePrinters[j].node == edge.sourceNode) {
                nodeAlreadyThere = true;
                break;
            }
        }
        if (nodeAlreadyThere) {

        }
        else {
            var nodeP = new NodePrinter(edge.sourceNode, NodePrinter.startX + ((linNum - 1) * NodePrinter.xInc), NodePrinter.startY,linNum);
            ++linNum;

            nodePrinters.push(nodeP);

        }

    });

    if (deb.edgePath[deb.edgePath.length - 1]) {
        var edge = deb.edgePath[deb.edgePath.length - 1];
        //add target node of last edge
        var nodeAlreadyThere = false;
        for (var j = 0;j < nodePrinters.length;++j) {
            if (nodePrinters[j].node == edge.targetNode) {
                nodeAlreadyThere = true;
                break;
            }
        }
        if (nodeAlreadyThere) {

        }
        else {
            var nodeP = new NodePrinter(edge.targetNode, NodePrinter.startX + ((linNum - 1) * NodePrinter.xInc), NodePrinter.startY,linNum);
            ++linNum;

            nodePrinters.push(nodeP);

        }

    }



    nodePrinters.forEach(function(nodeP) {
        nodeP.print(ctx);

   });

   ctx.stroke();

    //Print path

    for (var i = 0;i < deb.edgePath.length;++i) {
        ctx.beginPath();

        var stNode,endNode;

        stNode = null;
        endNode = null;
        for (var j = 0;j < nodePrinters.length; ++j) {
            if (nodePrinters[j].node == deb.edgePath[i].sourceNode) {
                stNode = nodePrinters[j];

            }
            if (nodePrinters[j].node == deb.edgePath[i].targetNode) {
                endNode = nodePrinters[j];

            }

        }
        stNode.join(ctx,endNode,deb.edgePath[i].walkNum);

        ctx.stroke();


    }



    ctx.stroke();
    ctx.closePath();



}

function logoCanvas(motifs) {


    var charHeightMatrix = calcMotifLogo(motifs,false);


    var res = scoreMotifs(motifs,false);

    var canv = document.getElementById("skewCanvas");
    var ctx = canv.getContext("2d");
    ctx.save();
    var H = canv.height;
    var W = canv.width;
    ctx.clearRect(0,0,W,H);
    ctx.fillStyle = 'white';
    ctx.fillRect(0,0,W,H);

    ctx.font = "30px Arial";
    ctx.scale(1,1.5);


    var consensus = res[0];




    var start = 0;

    var startPixelHeight = 70;

    var measured = [];

    var currPixelHeights = [];


    var maxWidth = 0.0;

    for (var c = 0;c < consensus.length; ++c) {
        var size = charHeightMatrix[0][c][1] * 100;
        ctx.font = size + 'px Arial';
        var wid = ctx.measureText(charHeightMatrix[0][c][0]).width;
        if (wid < 20) {
            wid = 20;
        }
        if (c == 0) {
           measured = [wid];

        }
        else {
           measured.push(measured[c - 1] + wid);
        }

        currPixelHeights.push(startPixelHeight);

        if (wid > maxWidth) {
            maxWidth = wid;
        }
    }

    var evenTot = maxWidth * consensus.length;

    var scale = 1.0;


    if (measured[measured.length - 1] > W) {
        scale = W * 1.0 / measured[measured.length - 1];
        ctx.scale(scale,scale);
    }

    /*
    if (evenTot > W) {
        scale = W * 1.0 / evenTot;
        ctx.scale(scale,scale);
    }
    */



    var charWidth = W  * 1.0 / consensus.length;
    if (charWidth > 40) {
       // charWidth = 40;
    }
    charWidth *=2;

    for (c = 0;c < consensus.length; ++c) {
        for (var r = 0;r < charHeightMatrix.length;++r) {
            size = charHeightMatrix[r][c][1] * 100;
            ctx.font = size + 'px Arial';


            /*
            if (size > 50) {
                size*= 0.3;
            }
            else if (size > 20) {
                size *= 0.7;
            }

            else if (size < 8) {
                size *= 2;
            }
            */


           // ctx.fillText(charHeightMatrix[r][c][0], charWidth * c, 70 + r * 20);
            switch (charHeightMatrix[r][c][0]) {
                case 'A':
                    ctx.fillStyle = 'red';
                    break;
                case 'C':
                    ctx.fillStyle = 'blue';
                    break;
                case 'G':
                   // ctx.fillStyle = 'yellow';
                    ctx.fillStyle = 'rgb(255,215,0)';
                    break;
                case 'T':
                    ctx.fillStyle = 'green';
                    break;
            }
            var thisMeasure;
            var centreOffset = 0;
            //if (r > 0) {
                if (c == 0) {
                    thisMeasure = measured[0];
                }
                else {
                    thisMeasure = measured[c] - measured[c-1];
                }
                var thisWidth = ctx.measureText(charHeightMatrix[r][c][0]).width;
                centreOffset = (thisMeasure - thisWidth) / 2.0;

            //}
            //ctx.fillText(charHeightMatrix[r][c][0], (c == 0) ? 0 + centreOffset : measured[c-1] + centreOffset, 70 + r * 15);
            if (r == 0) {

            }
            else {
                //var h = ctx.measureText(charHeightMatrix[r][c][0]).height;
                var h = ctx.measureText('M').width;
                currPixelHeights[c] += h;


            }
            ctx.fillText(charHeightMatrix[r][c][0], (c == 0) ? 0 + centreOffset : measured[c-1] + centreOffset, currPixelHeights[c]);

           // ctx.fillText(charHeightMatrix[r][c][0], maxWidth * c, 120 + r * 20);  FOR EVEN SPACING
            /*
            if (r == 0) {
                if (c == 0) {
                    measured = [ctx.measureText(charHeightMatrix[r][c][0]).width];

                }
                else {
                    measured.push(measured[c - 1] + ctx.measureText(charHeightMatrix[r][c][0]).width);
                }
            }
            */
        }

    }


    if (scale == 1.0) {

    }
    else {
        ctx.scale(1.0/scale,1.0/scale);
    }

    ctx.scale(1,0.666);

    ctx.restore();



}
function skewCanvas(skewData) {
    var c = document.getElementById("skewCanvas");

    var h = c.height;
    var w = c.width;

    var max = skewData[1];//arrayMax(skewArray);
    var min = skewData[0]; //arrayMin(skewArray);
    var span = max - min;

    var skewArray = skewData[2];

    var unitW;
    var numUnits;

    if ((skewArray.length -1 ) > w) {
        unitW = 1;
        numUnits = w;
    }
    else {
        unitW = Math.floor(w * 1.0 / (skewArray.length - 1));
        numUnits = skewArray.length;
    }


    var unitH;
    var numYUnits;

    if (span > h) {
        unitH = 1;
        numYUnits = h;
    }
    else {
        unitH = Math.floor(h * 1.0 / span);
        numYUnits = span;
    }

    //var unitH = Math.floor(h * 1.0 / (max - min));

    var scale = 1.0;
    if ((skewArray.length - 1) > w) {
        scale = Math.ceil((skewArray.length - 1) / w);
    }

    var scaleY = 1.0;
    if (span > h) {
        scaleY = Math.ceil(span / h);
    }


    var ctx = c.getContext("2d");
    ctx.clearRect(0,0,w,h);

    ctx.font = '10pt Arial';
    ctx.fillStyle = ('#FFFFFF');
    ctx.fillText('G-C Skew',2,15);
    ctx.fillText('Min: ' + min +  '[' + skewData[3] + ']',80    ,15);
    ctx.fillText('Max: ' + max +  '[' + skewData[4] + ']',190,15);




    ctx.beginPath();



    ctx.moveTo(0,0);
    ctx.lineTo(0,h);
    var zeroH = unitH * (max - 0) / scaleY;
    ctx.moveTo(0,zeroH);//h / 2.0);
    ctx.lineTo(w,zeroH);//h / 2.0);
    ctx.moveTo(0,zeroH);//h / 2.0);


    // skewArray.forEach(function(el,i) {

    for (var i = 1;i < numUnits ; ++i) {


        //  ctx.fillRect(unitW * (i+1),unitH * (max - el),4,4);
        var x = unitW * i;
        var y = unitH * ( (max - skewArray[i*scale]) / scaleY);
        ctx.lineTo(x,y);
        ctx.moveTo(x,y);
        //ctx.lineTo(unitW * (i + 1), unitH * (max - skewArray[i*scale]));
        // ctx.moveTo(unitW * (i + 1), unitH * (max - skewArray[i*scale]));

    }
    // });
    ctx.stroke();
    ctx.closePath();

}

function initialiseResults()  {
    var hData = ['k-mers', 'Num found'];
    var tData;
    tData = [['','Searching..']];

    var t = createTable(tData, hData, null, null, null, null, null, 'kmers_');
    t.className += ' codeTable';
    t.id = 'kMerResultsTab';

    var el = document.getElementById('kMerResultsTab');
    el.parentNode.removeChild(el);
    document.getElementById('rightOne').appendChild(t);

    document.getElementById('debugText').value = '';
    document.getElementById('moreDetailText').value = 'Extra Info:';


}

function colourDNATest(dna,mark) {
    var view =  document.getElementById('dnaView');
    view.innerHTML = dna;
}

function colourMotifs(dna,lowThisPage,highThisPage,inclRevCompl,view,bases,mark) {

    var numKmer = document.getElementById('numKmers');

    var totPos = 0;
   // for (var totPos = dnaPage * basesPerPage + dnaPageOffset; totPos <  (dnaPage + 1) * basesPerPage + dnaPageOffset; ++totPos) { //dna.length; ++totPos) {
    for (var strNum = 0;strNum < dna.length;++strNum) {
      //  if (totPos >= dna.length) {
      //      break;
      //  }
     for (var cNum = 0;cNum < dna[strNum].length; ++cNum) {
         var needHighlight = false;
         if (mfMotif.length > 0) {
             //m = mfk[parseInt(numKmer.value) - 1];
             m = mfMotif[parseInt(numKmer.value) - 1][0][strNum].slice(0,2);
            // var k = mfMotif[parseInt(numKmer.value) - 1][0][strNum][2].length;
             var k = mfMotif[parseInt(numKmer.value) - 1][1].length;
             m.forEach(function (mm) {
                 mm.forEach(function (mEntry) {

                     //console.log('totpos: ' + totPos + ' mEntry: ' + mEntry + ' m + k: ' + (mEntry+k));
                     if ((mEntry <= cNum ) && (cNum < (mEntry + k))) {
                         //console.log('aha: ' + totPos);
                         needHighlight = true;
                     }
                 });
             });
         }

         d = dna[strNum][cNum];
         /*
         if (d === '\n') {
             var b = document.createElement('br');
             view.appendChild(b);


         }
         else {
         */

             var s = document.createElement('span');


             //alert('bases dna: ' + bases[d]);
             s.className += bases[d];
             s.className += ' dnaChar';
             if (needHighlight) {
                 s.className += ' dnaCharHighlight';
             }
             if (mark) {
                 if (mark[cNum]) {
                     //  for (var ii = 0;ii < mark.length;++ii) {
                     //      if (mark[ii] == totPos) {
                     s.className += ' dnaClumpHighlight';
                 }

             }
             s.textContent = d;
             view.appendChild(s);
         /*
         }
         */
     }
        var b = document.createElement('br');
        view.appendChild(b);
    }


}



function colourMFK(dna,lowThisPage,highThisPage,inclRevCompl,view,bases,mark) {

    var numKmer = document.getElementById('numKmers');


    for (var totPos = dnaPage * basesPerPage + dnaPageOffset; totPos <  (dnaPage + 1) * basesPerPage + dnaPageOffset; ++totPos) { //dna.length; ++totPos) {
        if (totPos >= dna.length) {
            break;
        }

        var needHighlight = false;
        if (mfk.length > 0) {
            //m = mfk[parseInt(numKmer.value) - 1];
            m = mfk[parseInt(numKmer.value) - 1].slice(0, 2);
            k = mfk[parseInt(numKmer.value) - 1][2].length; //added 10/3/16

            m.forEach(function (mm) {
                mm.forEach(function (mEntry) {

                    //console.log('totpos: ' + totPos + ' mEntry: ' + mEntry + ' m + k: ' + (mEntry+k));
                    if ((mEntry <= totPos ) && (totPos < (mEntry + k))) {
                        //console.log('aha: ' + totPos);
                        needHighlight = true;
                    }
                });
            });
        }

        d = dna[totPos];
        if (d === '\n') {
            var b = document.createElement('br');
            view.appendChild(b);


        }
        else {

            var s = document.createElement('span');


            //alert('bases dna: ' + bases[d]);
            s.className += bases[d];
            s.className += ' dnaChar';
            if (needHighlight) {
                s.className += ' dnaCharHighlight';
            }
            if (mark) {
                if (mark[totPos]) {
                    //  for (var ii = 0;ii < mark.length;++ii) {
                    //      if (mark[ii] == totPos) {
                    s.className += ' dnaClumpHighlight';
                }

            }
            s.textContent = d;
            view.appendChild(s);
        }
    }


}


function colourDNA(dna,mark,inclRevCompl) {

    inclRevCompl = false;
    var motifView = false;

    switch (myParams.tabActive) {
        case 0:
            inclRevCompl = document.getElementById("includeRevComplMF").checked;
            break;
        case 1:
            inclRevCompl = document.getElementById("includeRevComplLT").checked;
            break;
        case 2:
            inclRevCompl = document.getElementById("includeRevComplKS").checked;
            break;
        case 3:
            inclRevCompl = false;// document.getElementById("includeRevComplMS").checked;
            motifView = true;
            break;
        default:
            break;
    }

    console.log('start colour');
    var view = document.getElementById('dnaView');
   /// console.log('start view inner html clear');
    // view.innerHTML = '';//null;

    while (view.firstChild) {
        view.removeChild(view.firstChild);
    }
   /// console.log('start bases defn');

    var lowThisPage = dnaPage * basesPerPage + 1 + dnaPageOffset;
    var highThisPage = (dnaPage + 1) * basesPerPage + dnaPageOffset;
    if (highThisPage > dna.length) {
        highThisPage = dna.length;
    }
    document.getElementById('currRangeDNA').innerHTML = 'Bases ' +  lowThisPage + '-' + highThisPage;
    var bases = {
        'A': 'aChar',
        'C': 'cChar',
        'G': 'gChar',
        'T': 'tChar',
        'a': 'aChar',
        'c': 'cChar',
        'g': 'gChar',
        't': 'tChar'


    };

    var numKmer = document.getElementById('numKmers');

   /// console.log('colour start loop');

    if (motifView) {
        if ((mfMotif.length == 0) || (mfMotif[0][0].length == 0)) {
            colourMotifs(dnaMasterStrings, lowThisPage, highThisPage, inclRevCompl, view, bases, mark);
        }
        else {
            colourMotifs(dnaMasterStrings, lowThisPage, highThisPage, inclRevCompl, view, bases, mark);
        }
    }
    else {
        colourMFK(dna,lowThisPage,highThisPage,inclRevCompl,view,bases,mark);
    }
    /*
    for (var totPos = dnaPage * basesPerPage + dnaPageOffset; totPos <  (dnaPage + 1) * basesPerPage + dnaPageOffset; ++totPos) { //dna.length; ++totPos) {
        if (totPos >= dna.length) {
            break;
        }

        var needHighlight = false;
        if (mfk.length > 0) {
            //m = mfk[parseInt(numKmer.value) - 1];
            m = mfk[parseInt(numKmer.value) - 1].slice(0, 2);
            m.forEach(function (mm) {
                mm.forEach(function (mEntry) {

                    //console.log('totpos: ' + totPos + ' mEntry: ' + mEntry + ' m + k: ' + (mEntry+k));
                    if ((mEntry <= totPos ) && (totPos < (mEntry + k))) {
                        //console.log('aha: ' + totPos);
                        needHighlight = true;
                    }
                });
            });
        }

        d = dna[totPos];
        if (d === '\n') {
            var b = document.createElement('br');
            view.appendChild(b);


        }
        else {

            var s = document.createElement('span');


            //alert('bases dna: ' + bases[d]);
            s.className += bases[d];
            s.className += ' dnaChar';
            if (needHighlight) {
                s.className += ' dnaCharHighlight';
            }
            if (mark) {
                if (mark[totPos]) {
                    //  for (var ii = 0;ii < mark.length;++ii) {
                    //      if (mark[ii] == totPos) {
                    s.className += ' dnaClumpHighlight';
                }

            }
            s.textContent = d;
            view.appendChild(s);
        }
    }
    */


   /// console.log('colour end loop');

    //  if (mfk.length > 0) {
    var hData;
    if (motifView)  {
        hData = ['Consensus', 'Num found'];
    }
    else {
        hData = ['k-mers', 'Num found'];
    }


    var tData;
    if (motifView) {
        if (mfMotif.length == 0) {
            tData = [['', 'None']];
        }
        else {
            tData = resolveMotif(dna, mfMotif, k, inclRevCompl); //[['ACCC','GGGG'],['GG','TT']];

        }

    }
    else {
        if (mfk.length == 0) {
            tData = [['', 'None']];

        }
        else {
            //var mpd = collapseMFK(dna,mfk,k);
            tData = resolveMFK(dna, mfk, k, inclRevCompl); //[['ACCC','GGGG'],['GG','TT']];
        }
    }

   /// console.log('colour start table');
    var t = createTable(tData, hData, null, null, null, null, null, 'kmers_');
    t.className += ' codeTable';
    t.id = 'kMerResultsTab';


   /// console.log('start onclick');
    var rows = t.getElementsByTagName('tbody')[0].getElementsByTagName('tr');
    for (i = 0; i < rows.length; i++) {
        rows[i].onclick = function () {
            //alert(this.rowIndex - 1 );
            var numKmer = document.getElementById('numKmers');
            numKmer.value = this.rowIndex;
            numKmerChanged(event);

        }
    }


    /*
     t.onclick = function(e) {
     e = e || event;
     var eventEl = e.srcElement || e.target,
     parent = eventEl.parentNode,
     isRow = function(el) {
     return el.tagName.match(/tr/i);
     };

     //move up the DOM until tr is reached
     while (parent = parent.parentNode) {
     if (isRow(parent)) {
     //row found, do something with it and return, e.g.
     alert(parent.rowIndex + 1);
     return true;
     }
     }
     return false;
     };
     */

    var el = document.getElementById('kMerResultsTab');
    el.parentNode.removeChild(el);
    document.getElementById('rightOne').appendChild(t);

   /// console.log(resolveMFK(dna, mfk, k));
    //}
    ///console.log('end colour');

}


function colourDNAProgressNew(dna,mark,clumpsDetected,clumpSize) {

    var regionSize = 50;

    var view =  document.getElementById('dnaView');
    view.innerHTML = null;
    var bases = {
        'A': 'aChar',
        'C': 'cChar',
        'G': 'gChar',
        'T': 'tChar',
        'a': 'aChar',
        'c': 'cChar',
        'g': 'gChar',
        't': 'tChar'


    };

    var numKmer = document.getElementById('numKmers');

    /*
     chunks.forEach(function(chunk,chNum) {

     var totPosBegThisLine = chNum*charsPerLine;
     var totPosBegNextLine = (chNum +1)*charsPerLine;
     if (mark[0] < totPosBegThisLine) {
     var s = document.createElement('span');
     //s.className += ' dnaCharHighlight';
     s.className += ' dnaChar';
     s.textContent = chunk;
     view.appendChild(s);

     }
     else
     if (mark[0] >= totPosBegNextLine) {
     var clumpsFoundFirstHalf = false;
     for (var ii = totPosBegThisLine;ii < (totPosBegNextLine - (charsPerLine /2 ));++ii) {
     if (clumpsDetected[ii]) {
     clumpsFoundFirstHalf = true;
     break;
     }

     }

     var clumpsFoundSecondHalf = false;
     for (var ii = totPosBegThisLine + (charsPerLine / 2);ii < totPosBegNextLine;++ii) {
     if (clumpsDetected[ii]) {
     clumpsFoundSecondHalf = true;
     break;
     }

     }


     var s1 = document.createElement('span');
     if (clumpsFoundFirstHalf) {
     s1.className += ' dnaCharHighlight';
     }
     s1.className += ' dnaChar';
     s1.textContent = chunk.substring(0,charsPerLine / 2);
     view.appendChild(s1);

     var s2 = document.createElement('span');
     if (clumpsFoundSecondHalf) {
     s2.className += ' dnaCharHighlight';
     }
     s2.className += ' dnaChar';
     s2.textContent = chunk.substring(charsPerLine / 2);
     view.appendChild(s2);


     }
     */

    var splitPoint = mark[0];

    for (p = 0;p < splitPoint;p+=regionSize) {

        var clumpsFound  = false;
        for (var ii = p;ii < p+regionSize; ++ii) {
            if (clumpsDetected[ii]) {
                clumpsFound = true;
                break;
            }

        }
        var s = document.createElement('span');
        s.className += ' dnaChar';
        var stopHere;
        if (splitPoint < p+regionSize) {
            stopHere = splitPoint;

        }
        else {
            stopHere = p+regionSize;
        }
        s.textContent = dna.substring(p,stopHere);
        if (clumpsFound) {
            s.className += ' dnaCharHighlight';
        }
        view.appendChild(s);


    }


    /*
     var s1 = document.createElement('span');
     s1.className += ' dnaChar';
     s1.textContent = dna.substring(0,splitPoint);
     view.appendChild(s1);
     */

    var s2 = document.createElement('span');
    s2.className += ' dnaChar';
    s2.textContent = dna.substring(splitPoint,splitPoint+1);//clumpSize);
    s2.className += ' dnaCharHighlight';
    view.appendChild(s2);

    var s3 = document.createElement('span');
    s3.className += ' dnaChar';
    s3.textContent = dna.substring(splitPoint+1);//clumpSize);

    view.appendChild(s3);







    /*
     for (var i = 0;i < chunk.length; ++i) {
     var totPos = chNum*charsPerLine + i;
     var needHighlight = false;
     // mfk.forEach(function(m) {
     //console.log('m: ' + m);
     if (mfk.length > 0) {
     m = mfk[parseInt(numKmer.value) - 1];
     m.forEach(function (mEntry) {
     //console.log('totpos: ' + totPos + ' mEntry: ' + mEntry + ' m + k: ' + (mEntry+k));
     if ((mEntry <= totPos ) && (totPos < (mEntry + k))) {
     //console.log('aha: ' + totPos);
     needHighlight = true;
     }
     });
     }


     //  });
     var s = document.createElement('span');
     d = chunk[i];

     //alert('bases dna: ' + bases[d]);
     s.className += bases[d];
     s.className += ' dnaChar';
     if (needHighlight) {
     s.className += ' dnaCharHighlight';
     }
     if (mark) {
     for (var ii = 0;ii < mark.length;++ii) {
     if (mark[ii] == totPos) {
     s.className += ' dnaCharHighlight';
     break;

     }
     }
     }
     s.textContent = d;
     view.appendChild(s);
     // console.log(d);
     }
     */

    //view.appendChild(document.createElement('br'));


    // });


}
/*
function colourDNAProgress(dna,mark,clumpsDetected) {

    var charsPerLine = 100;
    var currPos = 0;
    var chunks = [];
    while (currPos < dna.length) {
        chunks.push(dna.substring(currPos,currPos+charsPerLine));
        currPos += charsPerLine;
    }
    var view =  document.getElementById('dnaView');
    view.innerHTML = null;
    var bases = {
        'A': 'aChar',
        'C': 'cChar',
        'G': 'gChar',
        'T': 'tChar',
        'a': 'aChar',
        'c': 'cChar',
        'g': 'gChar',
        't': 'tChar'


    };

    var numKmer = document.getElementById('numKmers');

    var s1,s2;

    chunks.forEach(function(chunk,chNum) {

        var totPosBegThisLine = chNum*charsPerLine;
        var totPosBegNextLine = (chNum +1)*charsPerLine;
        if (mark[0] < totPosBegThisLine) {
            var s = document.createElement('span');
            //s.className += ' dnaCharHighlight';
            s.className += ' dnaChar';
            s.textContent = chunk;
            view.appendChild(s);

        }
        else
        if (mark[0] >= totPosBegNextLine) {
            var clumpsFoundFirstHalf = false;
            for (var ii = totPosBegThisLine;ii < (totPosBegNextLine - (charsPerLine /2 ));++ii) {
                if (clumpsDetected[ii]) {
                    clumpsFoundFirstHalf = true;
                    break;
                }

            }

            var clumpsFoundSecondHalf = false;
            for (ii = totPosBegThisLine + (charsPerLine / 2);ii < totPosBegNextLine;++ii) {
                if (clumpsDetected[ii]) {
                    clumpsFoundSecondHalf = true;
                    break;
                }

            }


            s1 = document.createElement('span');
            if (clumpsFoundFirstHalf) {
                s1.className += ' dnaCharHighlight';
            }
            s1.className += ' dnaChar';
            s1.textContent = chunk.substring(0,charsPerLine / 2);
            view.appendChild(s1);

            s2 = document.createElement('span');
            if (clumpsFoundSecondHalf) {
                s2.className += ' dnaCharHighlight';
            }
            s2.className += ' dnaChar';
            s2.textContent = chunk.substring(charsPerLine / 2);
            view.appendChild(s2);


        }
        else {
            var splitPoint = mark[0] - totPosBegThisLine;
            s1 = document.createElement('span');
            s1.className += ' dnaChar';
            s1.textContent = chunk.substring(0,splitPoint);
            view.appendChild(s1);

            s2 = document.createElement('span');
            s2.className += ' dnaChar';
            s2.textContent = chunk.substring(splitPoint,splitPoint+1);
            s2.className += ' dnaCharHighlight';
            view.appendChild(s2);

            var s3 = document.createElement('span');
            s3.className += ' dnaChar';
            s3.textContent = chunk.substring(splitPoint+1);

            view.appendChild(s3);


        }


*/

        /*
         for (var i = 0;i < chunk.length; ++i) {
         var totPos = chNum*charsPerLine + i;
         var needHighlight = false;
         // mfk.forEach(function(m) {
         //console.log('m: ' + m);
         if (mfk.length > 0) {
         m = mfk[parseInt(numKmer.value) - 1];
         m.forEach(function (mEntry) {
         //console.log('totpos: ' + totPos + ' mEntry: ' + mEntry + ' m + k: ' + (mEntry+k));
         if ((mEntry <= totPos ) && (totPos < (mEntry + k))) {
         //console.log('aha: ' + totPos);
         needHighlight = true;
         }
         });
         }


         //  });
         var s = document.createElement('span');
         d = chunk[i];

         //alert('bases dna: ' + bases[d]);
         s.className += bases[d];
         s.className += ' dnaChar';
         if (needHighlight) {
         s.className += ' dnaCharHighlight';
         }
         if (mark) {
         for (var ii = 0;ii < mark.length;++ii) {
         if (mark[ii] == totPos) {
         s.className += ' dnaCharHighlight';
         break;

         }
         }
         }
         s.textContent = d;
         view.appendChild(s);
         // console.log(d);
         }
         */
        /*
        view.appendChild(document.createElement('br'));


    })


}
*/


//Events

function dnaMasterChanged() {

    var skewRequired = false;

    switch (myParams.tabActive) {
        case 0:
            skewRequired = true;
            break;
        case 1:
            skewRequired = true;
            break;
        case 2:
            skewRequired = true;
            break;
        case 4:
            skewRequired = true;
            break;

        default:
            break;
    }


    dnaPage = 0;
    dnaPageOffset = 0;
    document.getElementById("jumpToBase").value = 1;


    //var sk = gcSkew(inDNA.value,inDNA.value.length);


    sk = null;


    //skTimeoutLoop(dnaMaster,1000);


    /*
     var sk;
     var basesDone = 0;
     while (basesDone <  inDNA.value.length) {
     sk =  gcSkew(inDNA.value,10,sk);
     basesDone = sk[2].length - 1;
     }
     */

    var stats = collectStats(dnaMaster,sk);
    /*
     var bc = baseCount(inDNA.value);

     var stats = 'Stats<br>DNA seq length: ' + inDNA.value.length + '<br>';

     var percStats = bc.map(function(el) {
     return [el,el * 1.0 /inDNA.value.length * 100.0];

     });

     var gcPerc = gcCount(inDNA.value) * 1.0 / inDNA.value.length * 100;

     stats+= 'A count: ' + bc[0] +  ' ' + percStats[0][1].toFixed(2) +  '%<br>';
     stats+= 'C count: ' + bc[1] + ' ' + percStats[1][1].toFixed(2) +  '%<br>';
     stats+= 'G count: ' + bc[2] + ' ' + percStats[2][1].toFixed(2) +  '%<br>';
     stats+= 'T count: ' + bc[3] + ' ' + percStats[3][1].toFixed(2) +  '%<br>';
     stats+= 'G-C %: ' + gcPerc.toFixed(2) + '%<br>';
     var sk;
     var basesDone = 0;
     while (basesDone <  inDNA.value.length) {
     sk =  gcSkew(inDNA.value,10,sk);
     basesDone = sk[2].length - 1;
     }
     //var sk = gcSkew(inDNA.value,inDNA.value.length);
     stats+= 'G-C skew min: ' + sk[0] + '<br>';
     stats+= 'G-C skew max: ' + sk[1] + '<br>';
     */

    //document.getElementById('dnaView').innerHTML = inDNA.value;

   // document.getElementById('dnaLength').innerHTML = dnaMaster.length;

    mfk = [];
    mfMotif = [];


    colourDNA(dnaMaster);

    document.getElementById('rightTwo').innerHTML = stats;

    //skewCanvas(gcSkew(inDNA.value,Math.min(300,inDNA.value.length - 1)));
    // skewCanvas(gcSkew(inDNA.value,inDNA.value.length));
    if (skewRequired) {
        skTimeoutLoop(dnaMaster, 1000);
        skewCanvas(sk);
    }

    //alert('posn: ' + getCaretCharacterOffsetWithin(document.getElementById('dnaView')));


}


function graphChanged(grph) {

    if (grph.dna) {
        document.getElementById('graphDnaDiv').innerHTML = grph.dna;

    }

    if (grph.reads) {
        if ((grph.sourceType == DGraph.fromPairedDna) || (grph.sourceType == DGraph.fromPairedReads)) {
            var tmp = grph.reads.map(function(el) {
                return '(' + el[0] + '|' + el[1] + ')';

            });
            document.getElementById('graphReadsDiv').innerHTML = tmp;
        }
        else {
            document.getElementById('graphReadsDiv').innerHTML = grph.reads;
        }
    }

    var adj = grph.getAdjList();
    if (adj) {
        var str = '';
        adj.forEach(function(el) {
            str+= el + '<br>';
        });
        document.getElementById('graphAdjListDiv').innerHTML = str;
    }

    document.getElementById('debugText').value = grph.sumInfo();


    graphCanvas(grph);


    if (paramObj.seqMethod == DGraph.seqTypeComp) {
        if (grph.dna) {
            document.getElementById('debugText').value += '\nOriginal: \n' + grph.dna;
        }
    }
    else {
        var recon;
        if (grph.dna) {
            document.getElementById('debugText').value += '\nOriginal: \n' + grph.dna;
            recon =  grph.edgePathReconstructed();

            document.getElementById('debugText').value += '\nReconstructed: \n' + recon;
            if (grph.dna === recon) {
                document.getElementById('debugText').value += '\nMatch!';
            }
            else if (grph.dna.length == recon.length ) {
                document.getElementById('debugText').value += '\nValid but no match!';
            }

        }
        else {
            recon =  grph.edgePathReconstructed();
            document.getElementById('debugText').value += '\nReconstructed: \n' + recon;
        }

        document.getElementById('debugText').value += '\nEdges: ' + grph.edgePathToText();

    }



    var nbps = grph.maximalNonBranchingPaths();

    document.getElementById('debugText').value += '\nMaximal Non Branching:';

    nbps.forEach(function(nbp) {
        document.getElementById('debugText').value += '\n'
        nbp.forEach(function(el,i) {

            document.getElementById('debugText').value += el;
            if (i == nbp.length - 1) {

            }
            else {
                document.getElementById('debugText').value += ' -> ';

            }
        });

    });

    document.getElementById('debugText').value += '\nMaximal Non Branching reconstructed:';

    nbps.forEach(function(nbp) {
        document.getElementById('debugText').value += '\n'
        nbp.forEach(function(el,i) {
            if (i == 0) {
                document.getElementById('debugText').value += el;
            }
            else {
                document.getElementById('debugText').value += el.substring(el.length -1);
            }



        });

    });









}

function cleanContents(contents,contentType) {

    var newContents;

    if (contentType) { //new regime

        switch (contentType) {
            case DGraph.fromDna:case DGraph.fromPairedDna: //dna. Can be on multiple lines - concatenated into one
                if (contents.substring(0, 1) === '>') {

                    var spl = contents.split(/\r\n|\r|\n/);
                    var parts = contents.split(/\r\n|\r|\n/g);
                    parts.shift();
                    var joined = parts.join('\n');
                    parts.shift(); // removes the first item from the array
                    contents = parts.join('_');
                }


                var newContents = contents.replace(/[^ACGT0123456789acgt]/gm, "");
                newContents = newContents.toUpperCase();
                return newContents;

                break;


            case DGraph.fromReads: //reads on multiple lines
                var parts = contents.split(/\r\n|\r|\n/g);
                parts = parts.filter(function(el) {
                    if (el.substring(0, 1) === '>') {
                        return false;
                    }

                    else {
                        return true;
                    }
                });
                parts = parts.map(function(el) {
                    var repl = el.replace(/[^ACGT0123456789acgt]/gm, "");
                    return repl.toUpperCase();

                });
                parts = parts.filter(function(el) {
                    if (el.length == 0) {
                        return false;
                    }
                    else {
                        return true;
                    }

                });
                newContents = parts.join('\n');
                return newContents;

                break;

            case DGraph.fromPairedReads: //reads on multiple lines
                var parts = contents.split(/\r\n|\r|\n/g);
                parts = parts.filter(function(el) {
                    if (el.substring(0, 1) === '>') {
                        return false;
                    }

                    else {
                        return true;
                    }
                });
                parts = parts.map(function(el) {
                    var repl = el.replace(/[^ACGT|0123456789acgt]/gm, "");
                    return repl.toUpperCase();

                });
                parts = parts.filter(function(el) {
                    if (el.length == 0) {
                        return false;
                    }
                    else {
                        return true;
                    }

                });
                newContents = parts.join('\n');
                return newContents;

                break;

            case DGraph.fromAdjList: //reads on multiple lines
                var parts = contents.split(/\r\n|\r|\n/g);
                parts = parts.filter(function(el) {
                    if (el.substring(0, 1) === '>') {
                        return false;
                    }

                    else {
                        return true;
                    }
                });
                parts = parts.map(function(el) {
                    var repl = el.replace(/[^ACGT ->|0123456789acgt]/gm, "");
                    return repl.toUpperCase();

                });
                parts = parts.filter(function(el) {
                    if (el.length == 0) {
                        return false;
                    }
                    else {
                        return true;
                    }

                });
                newContents = parts.join('\n');
                return newContents;

                break;

            default:
                break;
        }

    }
    else {
        if (contents.substring(0, 1) === '>') {

            var spl = contents.split(/\r\n|\r|\n/);
            var parts = contents.split(/\r\n|\r|\n/g);
            parts.shift();
            var joined = parts.join('\n');
            parts.shift(); // removes the first item from the array
            contents = parts.join('_');
        }


        newContents = contents.replace(/[^ACGTacgt]/gm, "");
        // var newContents =  contents.replace(/(\r\n|\n|\r)/gm,"");
        // newContents = newContents.replace(/\s+/g, '');
        // newContents = newContents.toUpperCase();

    }

    return newContents;
}

function dnaInput(e) {
    var inDNA = document.getElementById('dnaInput');


    inDNA.value = cleanContents(inDNA.value);
    inDNA.value = inDNA.value.toUpperCase();
   // inDNA.value = inDNA.value.replace(/(\r\n|\n|\r)/gm,"");
   // inDNA.value = inDNA.value.replace(/\s+/g, '');
   // inDNA.value = inDNA.value.toUpperCase();

    dnaMaster = inDNA.value;
    dnaMasterChanged();


}

function motifsInput(e) {

    var dnaStrings = document.getElementById('dnaStrings').value.split('\n');

    dnaStrings = dnaStrings.filter(function(el) {
        if (el.length == 0) {
            return false;
        }
        else {
            return true;
        }
    });
    dnaStrings = dnaStrings.map(function(el) {
        return el.trim().toUpperCase().replace(/[^ACGTacgt\n]/gm,"");

    });

    dnaMasterStrings =  dnaStrings;//document.getElementById('dnaStrings').value.split('\n');

    dnaMaster = document.getElementById('dnaStrings').value.replace(/ /g,'').toUpperCase();
    dnaMaster = dnaMaster.replace(/[^ACGTacgt\n]/gm,"");

    document.getElementById('dnaStrings').value = dnaMaster;

    dnaMasterChanged();

}

function sequencingInput(e) {

    var seqInp = document.getElementsByName('seqInput');

    var val = '';
    for (var i = 0;i < seqInp.length; ++i) {
        if (seqInp[i].checked) {
            val = seqInp[i].value;
        }
    }
    var cleaned;

        switch (val) {
          case 'seqDNA':
            cleaned =  cleanContents(document.getElementById('sequencingStrings').value,DGraph.fromDna);
            break;
          case 'seqReads':
              cleaned =  cleanContents(document.getElementById('sequencingStrings').value,DGraph.fromReads);
                break;

          case 'seqPairedDNA':
                cleaned =  cleanContents(document.getElementById('sequencingStrings').value,DGraph.fromPairedDna);
                break;
          case 'seqPairedReads':
                cleaned =  cleanContents(document.getElementById('sequencingStrings').value,DGraph.fromPairedReads);
                break;

            case 'seqAdjList':
                cleaned =  cleanContents(document.getElementById('sequencingStrings').value,DGraph.fromAdjList);
                break;


            default:
                break;



    }


    document.getElementById('sequencingStrings').value = cleaned;

    switch (val) {
        case 'seqDNA':case 'seqPairedDNA':
            dnaMaster = cleaned;
            dnaMasterChanged();
            break;
        case 'seqReads':
            dnaMasterStrings = cleaned;
            break;
        case 'seqPairedReads':
            dnaMasterStrings = cleaned;
            break;
        case 'seqAdjList':
            dnaMasterStrings = cleaned;
            break;

        default:
            break;
    }



    /*

    var sequencingStrings = document.getElementById('sequencingStrings').value.split('\n');

    sequencingStrings = sequencingStrings.filter(function(el) {
        if (el.length == 0) {
            return false;
        }
        else {
            return true;
        }
    });
    sequencingStrings = sequencingStrings.map(function(el) {
        return el.trim().toUpperCase().replace(/[^ACGTacgt\n]/gm,"");

    });

    dnaMasterStrings =  sequencingStrings;//document.getElementById('dnaStrings').value.split('\n');

    dnaMaster = document.getElementById('sequencingStrings').value.replace(/ /g,'').toUpperCase();
    dnaMaster = dnaMaster.replace(/[^ACGTacgt\n]/gm,"");

    document.getElementById('dnaStrings').value = dnaMaster;

    dnaMasterChanged();
    */

}


function restrictToACGT(event,allowNewline,otherCharsAllowed)
{
    var keynum;
    var keychar;
    var enttest;
    var validChars;

    allowNewline = allowNewline || false;

    otherCharsAllowed = otherCharsAllowed || '';

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
    validChars = "ACTGactg" + otherCharsAllowed;

    enttest = "\r";

    if (allowNewline) {
        if ((validChars.indexOf(keychar) > -1) || (enttest == keychar)) {
            return true;
        }
        else {
            return false;
        }

    }
    else {
        if (enttest == keychar) {
            event.srcElement.blur();
            return false;
        }
        else if (validChars.indexOf(keychar) > -1) {
            return true;
        }
        else {
            return false;
        }
    }
}

function stopPressed(e) {
    stop = true;

}



function kMerRevComplPressed() {
    var kmer = document.getElementById('kmer').value;
    document.getElementById('revComplKmerLabel').innerHTML = reverseComplement(kmer);

}


function kmerChanged(e) {


}

function kMerLenChanged(e) {
    k = parseInt(e.target.value);
    document.getElementById('numKmers').value = 1;
}
function numKmerChanged(e) {

    switch (myParams.tabActive) {
        case 3: //motif
            colourDNA(dnaMaster,null,inclRevCompl);
            return;
            break;
        default:
            break;
    }
    var numKmer = document.getElementById('numKmers');
    var numKmerVal = parseInt(numKmer.value);
    var dna = dnaMaster; //document.getElementById('dnaInput').value;

    var km = dna.substring(mfk[numKmerVal-1][0][0],mfk[numKmerVal-1][0][0] + k);
    document.getElementById('kMerMostFreq').innerHTML = dna.substring(mfk[numKmerVal-1][0][0],mfk[numKmerVal-1][0][0] + k);
    document.getElementById('kMerRevCompl').innerHTML =  reverseComplement(km);

    document.getElementById('kMerOccurrences').innerHTML = mfk[numKmerVal-1][0].length + mfk[numKmerVal-1][1].length;//mfk[0].length;
    document.getElementById('kMerNumEqualMostFreq').innerHTML = '' + mfk.length;

    var str = '';
    var strRev = '';
    mfk[numKmerVal - 1][0].forEach(function(el) {
        str += el + ' ';

    });
    mfk[numKmerVal - 1][1].forEach(function(el) {
        strRev += el + ' ';

    });

    document.getElementById("moreDetailText").innerHTML
        = km
        + '\n'
        +  str //mfk[numKmerVal-1][0]
       + '\n'
        + strRev;//mfk[numKmerVal-1][1];

    var inclRevCompl = false;
    switch (myParams.tabActive) {
        case 0:
            inclRevCompl = document.getElementById("includeRevComplMF").checked;
            break;
        case 1:
            inclRevCompl = document.getElementById("includeRevComplLT").checked;
            break;
        case 2:
            inclRevCompl = document.getElementById("includeRevComplKS").checked;
            break;
        default:
            break;
    }

    colourDNA(dnaMaster,null,inclRevCompl);

}



function randPressed(event) {
    var randN = document.getElementById('numRand');
    var randNVal = parseInt(randN.value);


   // document.getElementById('dnaInput').value = randomDNA(randNVal);
   // dnaInput();

    dnaMaster =  randomDNA(randNVal);
    dnaMasterChanged();

    document.getElementById('dnaInput').value = '';
}


function randMotifPressed(event) {

    var numSequencesToGen = 10;

    var randN = document.getElementById('numMotifRand');
    var randNVal = parseInt(randN.value);



    var seqLines = '';

    for (var i = 0;i < numSequencesToGen;++i) {
        var dna = randomDNA(randNVal);
        seqLines+=dna + '\n';
    }

    document.getElementById('dnaStrings').value = seqLines;
    motifsInput();


}


function jumpToPressed() {

    var jumpToBase = parseInt(document.getElementById("jumpToBase").value) - 1;

    if (jumpToBase < 0) {
        jumpToBase = 0;
        document.getElementById("jumpToBase").value = jumpToBase + 1;
    }
    if (jumpToBase >= dnaMaster.length) {
        jumpToBase = dnaMaster.length - 1;
        document.getElementById("jumpToBase").value = jumpToBase + 1;
    }


    document.getElementById('currPosInDNA').innerHTML = '';

    var dna = dnaMaster; //document.getElementById('dnaInput').value;

    dnaPage = Math.floor(jumpToBase / basesPerPage);

    dnaPageOffset = jumpToBase % basesPerPage;

    colourDNA(dna);


}
function pagePressed(inc) {

    document.getElementById('currPosInDNA').innerHTML = '';

    var dna = dnaMaster; //document.getElementById('dnaInput').value;

    if (inc == -1000) {
        dnaPage = 0;
        dnaPageOffset = 0;
        colourDNA(dna);
    }
    else if (inc == 1000) {
        dnaPage = Math.floor(dna.length / basesPerPage);
        dnaPageOffset = 0;
        colourDNA(dna);
    }
    else if (inc < 0) {
        if (dnaPage <= 0) {
            dnaPageOffset = 0;
            colourDNA(dna);
        }
        else {
            --dnaPage;
            colourDNA(dna);
        }
    }
    else {
        if ((dnaPage + 1) * basesPerPage > dna.length) {

        }
        else {
            ++dnaPage;
            colourDNA(dna);
        }
    }

}


function motifRadClicked(id) {
    //Motif Radio Button
    switch (id) {
        case 'motifKD':
            document.getElementById('maxMismatchMS').style.display = "block";
            document.getElementById('numItersMS').style.display = "none";
            document.getElementById('laplaceMS').style.display = "none";
            document.getElementById('maxMismatchLabMS').style.display = "block";
            document.getElementById('numItersLabMS').style.display = "none";
            document.getElementById('laplaceLabMS').style.display = "none";

            break;


        case 'motifBrute':
        case 'motifMedian':
            document.getElementById('maxMismatchMS').style.display = "none";
            document.getElementById('numItersMS').style.display = "none";
            document.getElementById('laplaceMS').style.display = "none";
            document.getElementById('maxMismatchLabMS').style.display = "none";
            document.getElementById('numItersLabMS').style.display = "none";
            document.getElementById('laplaceLabMS').style.display = "none";

            break;

        case 'motifGreedy':
            document.getElementById('maxMismatchMS').style.display = "none";
            document.getElementById('numItersMS').style.display = "none";
            document.getElementById('laplaceMS').style.display = "block";
            document.getElementById('maxMismatchLabMS').style.display = "none";
            document.getElementById('numItersLabMS').style.display = "none";
            document.getElementById('laplaceLabMS').style.display = "block";

            break;
        default:
            document.getElementById('maxMismatchMS').style.display = "none";
            document.getElementById('numItersMS').style.display = "block";
            document.getElementById('laplaceMS').style.display = "block";
            document.getElementById('maxMismatchLabMS').style.display = "none";
            document.getElementById('numItersLabMS').style.display = "block";
            document.getElementById('laplaceLabMS').style.display = "block";
            break;

    }

}

function sequencingInputRadClicked(id) {
    //Motif Radio Button
    switch (id) {
        case 'seqDNA':
            document.getElementById('kmerLenSA').style.display = "block";
            document.getElementById('numItersSA').style.display = "none";
            document.getElementById('pairDistSA').style.display = "none";
            document.getElementById('kmerLenLabSA').style.display = "block";
            document.getElementById('numItersLabSA').style.display = "none";
            document.getElementById('pairDistLabSA').style.display = "none";
            document.getElementById('seqGraphTypeHam').style.display = "block";
            document.getElementById('seqGraphTypeHamLab').style.display = "block";
            document.getElementById('seqGraphTypeDeb').style.display = "block";
            document.getElementById('seqGraphTypeDebLab').style.display = "block";



            break;

        case 'seqReads':

            document.getElementById('kmerLenSA').style.display = "none";
            document.getElementById('numItersSA').style.display = "none";
            document.getElementById('pairDistSA').style.display = "none";
            document.getElementById('kmerLenLabSA').style.display = "none";
            document.getElementById('numItersLabSA').style.display = "none";
            document.getElementById('pairDistLabSA').style.display = "none";
            document.getElementById('seqGraphTypeHam').style.display = "block";
            document.getElementById('seqGraphTypeHamLab').style.display = "block";
            document.getElementById('seqGraphTypeDeb').style.display = "block";
            document.getElementById('seqGraphTypeDebLab').style.display = "block";


            break;


        case 'seqPairedReads':

            document.getElementById('kmerLenSA').style.display = "none";
            document.getElementById('numItersSA').style.display = "none";
            document.getElementById('pairDistSA').style.display = "block";
            document.getElementById('kmerLenLabSA').style.display = "none";
            document.getElementById('numItersLabSA').style.display = "none";
            document.getElementById('pairDistLabSA').style.display = "block";
            document.getElementById('seqGraphTypeHam').style.display = "block";
            document.getElementById('seqGraphTypeHamLab').style.display = "block";
            document.getElementById('seqGraphTypeDeb').style.display = "block";
            document.getElementById('seqGraphTypeDebLab').style.display = "block";

            break;

        case 'seqPairedDNA':

            document.getElementById('kmerLenSA').style.display = "block";
            document.getElementById('numItersSA').style.display = "none";
            document.getElementById('pairDistSA').style.display = "block";
            document.getElementById('kmerLenLabSA').style.display = "block";
            document.getElementById('numItersLabSA').style.display = "none";
            document.getElementById('pairDistLabSA').style.display = "block";
            document.getElementById('seqGraphTypeHam').style.display = "block";
            document.getElementById('seqGraphTypeHamLab').style.display = "block";
            document.getElementById('seqGraphTypeDeb').style.display = "block";
            document.getElementById('seqGraphTypeDebLab').style.display = "block";

            break;

        case 'seqAdjList':

            document.getElementById('kmerLenSA').style.display = "none";
            document.getElementById('numItersSA').style.display = "none";
            document.getElementById('pairDistSA').style.display = "none";
            document.getElementById('kmerLenLabSA').style.display = "none";
            document.getElementById('numItersLabSA').style.display = "none";
            document.getElementById('pairDistLabSA').style.display = "none";
            document.getElementById('seqGraphTypeHam').style.display = "none";
            document.getElementById('seqGraphTypeHamLab').style.display = "none";
            document.getElementById('seqGraphTypeDeb').style.display = "block";
            document.getElementById('seqGraphTypeDebLab').style.display = "block";

            document.getElementById('seqGraphTypeDeb').checked = true;

            break;




        default:
            document.getElementById('kmerLenSA').style.display = "block";
            document.getElementById('numItersSA').style.display = "none";
            document.getElementById('kmerLenLabSA').style.display = "block";
            document.getElementById('numItersLabSA').style.display = "none";


            break;

    }

}


function sequencingRadClicked(id) {
    //Motif Radio Button
    switch (id) {
        case 'seqComp':
            document.getElementById('kmerLenSA').style.display = "block";
            document.getElementById('numItersSA').style.display = "none";
            document.getElementById('kmerLenLabSA').style.display = "block";
            document.getElementById('numItersLabSA').style.display = "none";


            break;

        case 'seqPath':

            document.getElementById('kmerLenSA').style.display = "block";
            document.getElementById('numItersSA').style.display = "none";
            document.getElementById('kmerLenLabSA').style.display = "block";
            document.getElementById('numItersLabSA').style.display = "none";


            break;


        default:
            document.getElementById('kmerLenSA').style.display = "block";
            document.getElementById('numItersSA').style.display = "none";
            document.getElementById('kmerLenLabSA').style.display = "block";
            document.getElementById('numItersLabSA').style.display = "none";


            break;

    }

}


function clickInDNA(e) {
    //alert('key down. Posn: ' + getCaretCharacterOffsetWithin(document.getElementById('dnaView')));

    var pos = getCaretCharacterOffsetWithin(document.getElementById('dnaView')) + (dnaPage * basesPerPage) + 1;
    document.getElementById('currPosInDNA').innerHTML =  pos;//getCaretCharacterOffsetWithin(document.getElementById('dnaView'));
}


function tabClickDone(tabNum) {
    //alert('tab clicked: ' + tabNum);
    switch (tabNum) {
        case 3: //Motif tab
             document.getElementById('mfkTools').style.display = "none";
             document.getElementById('motifTools').style.display = "block";
             document.getElementById('sequencingTools').style.display = "none";

            document.getElementById('dnaviewer').style.display = "inline-block";
            document.getElementById('graphviewer').style.display = "none";
            break;
        case 4:  //Sequencing tab
            document.getElementById('mfkTools').style.display = "none";
            document.getElementById('motifTools').style.display = "none";
            document.getElementById('sequencingTools').style.display = "block";

            document.getElementById('dnaviewer').style.display = "none";
            document.getElementById('graphviewer').style.display = "inline-block";
            break;
        default:
            document.getElementById('mfkTools').style.display = "block";
            document.getElementById('motifTools').style.display = "none";
            document.getElementById('sequencingTools').style.display = "none";

            document.getElementById('dnaviewer').style.display = "inline-block";
            document.getElementById('graphviewer').style.display = "none";
            break;

    }
}

//General utilities




function resolveMotif(dna,mfMotif,k,inclRevCompl) {

    var resolved = [];

    var mapped = mfMotif.map(function(el) {
        var rEntry = [];

        var tmp;
        var len = 0;

        tmp = el[1];


        if (inclRevCompl) {
            //tmp+='[' + el[0].length + ']';
            var revTmp = reverseComplement(el[1]);

            tmp += '/' + revTmp; //+ '[' + el[1].length + ']';

        }


        if (el[0].length == 0) {
            len = 0;
        }
        else {
            len = el[0][0].length + el[0][1].length;
        }


        rEntry.push(tmp);
        rEntry.push(len);

        return rEntry;





    });


    /*
    var mapped2 = [];

    var alreadyFound = [];
    for (var i = 0;i < mapped.length;++i) {

        mapped2.push(mapped[i]);

    }

    return mapped2;

    */
    return mapped;

}


function resolveMFK(dna,mfk,k,inclRevCompl) {

    var resolved = [];

    var mapped = mfk.map(function(el) {
        var rEntry = [];

        var tmp;
        var len = 0;

        tmp = el[2];


        if (inclRevCompl) {
            tmp+='[' + el[0].length + ']';
            var revTmp = reverseComplement(el[2]);

            tmp += '/' + revTmp + '[' + el[1].length + ']';

        }


        len = el[0].length + el[1].length;

        /*
        el.forEach(function(indAr,i) {
            if (i == 0) {
                tmp = dna.substring(indAr[0], indAr[0] + k);
                len+= indAr.length;
            }
            else if (i == 1) {
                if (indAr.length == 0) {

                }
                else {
                    var revTmp = reverseComplement(tmp);
                    tmp += '/' + revTmp;//dna.substring(indAr[0], indAr[0] + k);
                }
                len+= indAr.length;
            }


            //if (revComplAlreadyDone.indexOf(tmp) > -1) {
            //    tmp += '/' + reverseComplement(tmp);
            //}
        });
        */
        rEntry.push(tmp);
        rEntry.push(len);

        return rEntry;





    });

    /*
    if (inclRevCompl) {
        mapped.forEach(function (el, i) {
            var test = reverseComplement(el[0]);
            for (var j = i + 1; j < mapped.length; ++j) {
                if (mapped[j][0] === test) {
                    // el[1] += mapped[j][1];
                    el[0] += '/' + mapped[j][0];
                    mapped.splice(j, 1);


                }
            }

        });
    }
    */

    /*
    mfk.forEach(function(el) {
        rEntry = [];
        var tmp= dna.substring(el[0],el[0] + k);
        if (revComplAlreadyDone.indexOf(tmp) > -1) {
            tmp += '/' + reverseComplement(tmp);
        }
        rEntry.push(tmp);
        rEntry.push(el.length);

        resolved.push(rEntry);

    });
    */

    /*
     resolved = mfk.filter(function(el) {
     var tmp= dna.substring(el[0],el[0] + k);
     return dna.substring(el[0],el[0] + k);

     });
     */
/*
    return resolved;
    */

    var mapped2 = [];

    var alreadyFound = [];
    for (var i = 0;i < mapped.length;++i) {
        /*
        var spl = mapped[i][0].split('/');
        if (alreadyFound.indexOf(spl[0]) > -1) {
            continue;
        }
        else {
            spl.forEach(function(el) {
                alreadyFound.push(el);

            });
            mapped2.push(mapped[i]);
        }
        */
        mapped2.push(mapped[i]);

    }

    return mapped2;



}


function collapseMFK(dna,mfk,k) {

    var resolved = [];

    var alreadyHandled = [];
    var mapped = mfk.map(function(el,i) {
        if (alreadyHandled.indexOf(i) > -1) {

        }
        else {

            var mpdEntry = [];
            var tmp = dna.substring(el[0], el[0] + k);
            var revCompl = reverseComplement(tmp);

            var revFound = false;

            for (var j = i + 1; j < mfk.length; ++j) {
                if (revCompl === dna.substring(mfk[j][0], mfk[j][0] + k)) {
                    mpdEntry = [el, mfk[j]];
                    revFound = true;
                    break;
                }
            }
            if (revFound) {

            }
            else {
                mpdEntry = [el];
            }
            return mpdEntry;
        }


    });


    return mapped;



}


function getCaretCharacterOffsetWithin(element) {
    var caretOffset = 0;
    var doc = element.ownerDocument || element.document;
    var win = doc.defaultView || doc.parentWindow;
    var sel;
    if (typeof win.getSelection != "undefined") {
        sel = win.getSelection();
        if (sel.rangeCount > 0) {
            var range = win.getSelection().getRangeAt(0);
            var preCaretRange = range.cloneRange();
            preCaretRange.selectNodeContents(element);
            preCaretRange.setEnd(range.endContainer, range.endOffset);
            caretOffset = preCaretRange.toString().length;
        }
    } else if ( (sel = doc.selection) && sel.type != "Control") {
        var textRange = sel.createRange();
        var preCaretTextRange = doc.body.createTextRange();
        preCaretTextRange.moveToElementText(element);
        preCaretTextRange.setEndPoint("EndToEnd", textRange);
        caretOffset = preCaretTextRange.text.length;
    }
    return caretOffset;
}

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
