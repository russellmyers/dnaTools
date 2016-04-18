
//Global Variables

var k = 2;
var mfk = [];
var mfMotif = []; //used for motif searches

var sk;

var dnaMaster = '';
var dnaPage = 0;
var dnaPageOffset = 0;

var dnaMasterStrings = []; // used for Motif search

var rnaMaster = ''; //used for Transcription/translation
var proteinMaster = '';
var spectrumMaster = '';
var sequencedWeights = '';
var convMaster = '';

var progExtraInfo; // extra info during background processing


var basesPerPage = 20000;

var clumps = [];
var stop = false;
var myParams = Params.getInstance();


//used for expanding/collapsing divs
var expDebugState = true;
var expMultiState = true;

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

/*
    var rands = [];
    for (var i = 0;i < 100; ++i) {
        rands.push(getRandomFromProbDist([0.05,0.05,0.05,0.7,0.05,0.05]));
    }
    */

    /*
    var pep = new Peptide(Peptide.AminoArrFromArr(['A','V','G','134','Q','M','M']));

    var am = Amino.AllUniqueWeightAminos(true);

    //var leaders = 'LAST ALST TLLT TQAS'.split(' ');
    var leaders = 'GFAQHVMEGIGLDVKFTNIISCFFDHEWSTCHCKHHNSINHTMSMVF LIGDDDEADNCMMMVQSIKWKTLLRYGAFFTFPFYSYAILHVFYVLW KPMWWAFIFGFCDMKNCFDAPFWMHNSVQWEQHYRCNDVKMMSQLCW MAPRDIRMYFDKYHETAALDSQWIIQQIYHLMNVRKLNRTNRFTSVG FEKYHQQQILIDAQRVRLVHTVARAGPGWVQTGGWQQTCPRYKPYAW NVNPCERSSPPNFSWFMSFWADNSDYGDVIFCCPSVLRTMEMQSKKG WDTDTFFQKAMLKKDETADQIFNLRPYSLTCHNENILGNDNQEKQAG TLGSGENDKGHTVGAGHKGHPEREFEAPIERHEHPRVMMTKVGCYWI VCGHHHEQTVIMKAFDAWKVGFLGPIVAWVIFPAVYLWGKSLCPWTN YDSPTTYLSTHCHRLTNRMVHENPVICPPQDFAKYLIQSGWEFPLVA KDPINQTGDTNVRNFNVGCFCGCYFQWERHDGTPMHFWFSQKLSLTW HMKKLFWGIMKHHILFDFVNQPAFTNKAKGPTPHKAEELIRNLGQEK FNDRQRLVCHTNQCCAYKNKVVCSGGGSEISTNAHTYHFLALGHQVG MYYSAWTEPYYPPTLQIWWWYWKYGCTACQTGPHTMVFVQPTCKCVH YYGYRQCSWCQRWTVRRMLCWIDVLHKALHWHVCLLFHQALYGFSHE WASIGAIMRSAKDMYESLEFHKTHCTYFVYMVCKEARPGWTFFIEWV'.split(' ');
    leaders = leaders.map(function(el) {
       var pep = new Peptide(Peptide.AminoArrFromStr(el));
       return pep;
    });
  // var spec = '0 71 87 101 113 158 184 188 259 271 372'.split(' ');
    var spec = '0 71 87 97 97 99 101 101 101 101 101 103 103 113 113 113 113 113 113 113 113 114 114 115 115 128 128 129 129 129 129 129 129 131 131 131 137 137 147 147 156 156 156 163 163 163 186 186 198 200 200 202 202 204 214 214 216 216 216 226 228 230 234 242 242 242 242 242 242 244 245 253 257 257 257 259 260 266 266 268 269 271 276 276 276 276 278 283 284 285 287 293 294 299 301 303 313 317 317 317 327 329 331 343 343 347 354 355 355 356 359 363 363 370 370 371 371 372 379 382 384 388 389 397 397 400 405 407 408 408 413 413 415 415 416 418 418 420 422 428 430 430 432 434 439 444 455 456 458 459 473 476 484 484 485 485 487 499 500 501 506 507 510 510 511 513 515 518 522 526 527 528 529 529 531 533 533 537 540 541 544 544 545 547 558 559 562 568 569 571 572 574 585 586 588 597 598 607 610 612 616 619 620 624 624 625 626 631 636 641 644 646 650 657 657 660 662 663 669 669 670 671 672 674 675 678 681 684 685 689 691 692 696 698 699 700 700 701 711 733 733 735 738 738 739 741 753 771 772 772 773 775 775 775 778 779 782 783 783 786 788 789 794 794 797 798 798 800 801 802 804 805 806 808 810 813 815 828 837 840 846 854 862 864 864 866 869 882 884 885 886 888 889 899 901 902 902 903 904 907 908 911 911 914 914 924 925 926 928 935 935 937 937 941 941 945 952 955 961 961 975 975 977 984 987 988 992 995 999 1002 1013 1015 1016 1017 1017 1017 1022 1025 1032 1032 1038 1039 1039 1040 1044 1051 1055 1058 1058 1059 1062 1065 1066 1070 1070 1072 1074 1084 1088 1097 1099 1100 1101 1103 1104 1105 1106 1118 1118 1121 1130 1133 1135 1142 1150 1151 1153 1153 1154 1154 1156 1165 1171 1172 1186 1187 1194 1196 1198 1200 1201 1201 1202 1207 1212 1213 1214 1215 1216 1217 1217 1218 1218 1219 1231 1233 1234 1234 1236 1248 1255 1259 1264 1267 1272 1279 1281 1285 1298 1300 1301 1309 1311 1315 1315 1315 1316 1318 1319 1319 1321 1325 1330 1331 1334 1335 1338 1341 1343 1344 1346 1347 1348 1352 1363 1363 1372 1379 1395 1396 1398 1402 1410 1414 1414 1422 1422 1428 1429 1429 1430 1431 1432 1433 1434 1435 1435 1438 1447 1448 1450 1452 1460 1461 1471 1472 1472 1476 1476 1478 1481 1494 1499 1509 1515 1517 1521 1525 1532 1534 1535 1535 1535 1543 1546 1548 1551 1557 1557 1561 1561 1561 1562 1566 1566 1573 1577 1581 1585 1585 1589 1603 1605 1606 1608 1609 1609 1612 1618 1628 1634 1635 1636 1638 1647 1663 1664 1664 1664 1672 1672 1674 1674 1677 1680 1682 1686 1690 1694 1699 1702 1708 1713 1715 1716 1718 1718 1722 1723 1724 1740 1741 1743 1748 1749 1756 1763 1777 1777 1778 1781 1785 1792 1801 1803 1803 1805 1810 1811 1811 1815 1819 1819 1822 1823 1827 1830 1837 1837 1842 1844 1849 1850 1853 1869 1872 1876 1878 1878 1879 1887 1905 1906 1906 1912 1916 1916 1916 1924 1925 1929 1934 1936 1939 1940 1940 1940 1948 1948 1950 1965 1972 1979 1979 1991 1993 1996 2000 2000 2005 2007 2007 2009 2017 2019 2019 2026 2031 2034 2035 2039 2042 2053 2053 2058 2061 2064 2076 2079 2079 2087 2087 2094 2096 2097 2103 2113 2118 2120 2120 2123 2125 2132 2132 2135 2135 2136 2138 2140 2147 2150 2167 2171 2171 2192 2193 2200 2200 2200 2205 2209 2216 2216 2219 2224 2226 2226 2232 2233 2233 2233 2235 2237 2250 2252 2262 2266 2266 2268 2276 2284 2287 2296 2300 2303 2313 2318 2334 2334 2334 2336 2336 2337 2339 2347 2349 2355 2355 2356 2363 2363 2365 2366 2379 2380 2389 2397 2397 2413 2413 2416 2418 2434 2435 2437 2446 2447 2449 2449 2452 2466 2468 2468 2473 2476 2490 2492 2493 2494 2494 2494 2502 2510 2511 2519 2526 2533 2536 2542 2547 2549 2550 2553 2560 2578 2579 2579 2586 2597 2605 2605 2615 2620 2622 2623 2624 2624 2625 2627 2631 2632 2634 2648 2650 2650 2650 2651 2655 2689 2691 2706 2710 2711 2715 2733 2733 2733 2735 2738 2742 2742 2744 2747 2751 2753 2756 2761 2763 2768 2771 2779 2802 2807 2813 2820 2824 2825 2828 2836 2836 2836 2843 2846 2862 2862 2866 2871 2873 2876 2876 2881 2882 2900 2903 2907 2924 2933 2933 2937 2938 2944 2944 2949 2953 2957 2965 2972 2975 2975 2975 2987 2995 2999 3004 3018 3034 3034 3037 3037 3038 3039 3058 3062 3066 3066 3070 3073 3078 3078 3082 3088 3101 3104 3112 3119 3130 3135 3138 3149 3151 3152 3163 3163 3171 3172 3179 3190 3191 3193 3195 3217 3217 3220 3229 3238 3241 3241 3243 3244 3248 3250 3251 3276 3286 3291 3292 3300 3308 3318 3320 3321 3330 3335 3342 3351 3351 3354 3354 3357 3358 3358 3372 3387 3405 3407 3420 3422 3429 3431 3433 3433 3434 3439 3448 3455 3471 3483 3485 3486 3486 3488 3510 3514 3521 3533 3534 3534 3535 3546 3551 3552 3568 3585 3596 3599 3599 3600 3611 3611 3614 3615 3622 3634 3635 3647 3649 3664 3664 3681 3682 3696 3697 3708 3713 3727 3728 3728 3728 3735 3748 3750 3751 3771 3777 3797 3809 3812 3827 3828 3837 3841 3841 3842 3851 3857 3864 3864 3868 3868 3884 3898 3913 3940 3940 3942 3943 3955 3965 3969 3970 3970 3977 3981 4013 4014 4027 4027 4044 4053 4054 4056 4057 4083 4096 4098 4099 4110 4126 4140 4140 4140 4145 4155 4158 4167 4171 4184 4209 4211 4212 4223 4253 4255 4259 4268 4272 4284 4296 4296 4299 4303 4313 4352 4368 4373 4374 4397 4397 4397 4400 4409 4409 4416 4428 4465 4469 4487 4501 4510 4510 4526 4529 4538 4560 4566 4572 4584 4630 4639 4639 4639 4643 4651 4673 4673 4681 4685 4752 4752 4752 4768 4782 4786 4786 4802 4829 4853 4867 4881 4881 4883 4915 4915 4942 4968 4968 4982 4994 5028 5044 5069 5069 5071 5095 5097 5157 5157 5170 5184 5198 5210 5258 5270 5299 5311 5313 5371 5373 5412 5426 5474 5486 5527 5575 5587 5642 5688 5743 5844'.split(' ');
    spec = spec.map(function(el) {
        return parseInt(el);
    });
   var N = 5;
   var trimmed = trimLeaderboard(leaders,spec,N);

    var str = '';
   trimmed.forEach(function(el) {
      str += el.toShortString('') + ' ';
   });

    console.log('trimmed leader: ' + str);

    var spc = new Spectrum('0 137 186 323');

    var top = spc.topElements(3);

    var conv = spc.convolution();


   // var cyc = cycloPeptideSequencing('0 113 128 186 241 299 314 427');
    var cyc = cyclopeptideSequencing('0 113 114 128 129 242 242 257 370 371 484');


    //var all = allPossiblePeptidesWithWeight(1024);

    */

    /*

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


    var tstA = Amino.AllUniqueWeightAminos();
    */

    /*
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

     */
}

function initialisePage() {

    myParams.tabActive = 5;
    tab_click(7);
    tab_click(6); //tab stuff
    tab_click(5);
    tab_click(4);
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


    document.getElementById('fileSequencingInput')
        .addEventListener('change',readSingleSequencingFile,false);


//  document.getElementById('motifBrute').checked = true;
    motifRadClicked('motifBrute');
    sequencingRadClicked('seqPath');
    sequencingInputRadClicked('seqDNA');

    document.getElementById('trTranscribe').checked = true;
    transRadClicked('trTranscribe');

    document.getElementById('pepIdealSpectrum').checked = true;
    peptideRadClicked('pepIdealSpectrum');


    expDebugState = false;
    expStateChanged('expDebug',false);
    expMultiState = true;
    expStateChanged('expMulti',true);


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

                case 'seqCyclopeptide':
                    /*
                    var resEl = document.getElementById('miscQuickResult');
                    resEl.value = 'results\n';
                    resEl.value += e.data.txtStuff;
                    */
                    processReturnedPeptide(e);
                    break;

                case 'seqLeaderboardCyclopeptide':
                  
                    /*
                    var resEl = document.getElementById('miscQuickResult');
                    resEl.value = 'results\n';
                    resEl.value += e.data.txtStuff;
                    */
                    processReturnedPeptide(e);
                    break;

                case 'seqLeaderboardConvCyclopeptide':

                    /*
                     var resEl = document.getElementById('miscQuickResult');
                     resEl.value = 'results\n';
                     resEl.value += e.data.txtStuff;
                     */
                    processReturnedPeptide(e);
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

    if (e.data.extraInfo) {
        switch (e.data.task) {
            case 'seqLeaderboardCyclopeptide':
                progExtraInfo = e.data.extraInfo;
                colourPep();
                break;
            case 'seqLeaderboardConvCyclopeptide':
                progExtraInfo = e.data.extraInfo;
                colourPep();
                break;
            default:
                break;
        }
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


  //  var numKmerVal = parseInt(numKmer.value);
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


        document.getElementById('debugText').value = e.data.txtStuff;

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

function processReturnedPeptide(e) {

    sequencedWeights = e.data.txtStuff;
    colourDNA(null,null,false);

    if (document.getElementById('debugPS').checked) {
       document.getElementById('debugText').value = e.data.txtStuff;

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

function readSingleSequencingFile(e) {
    //for genome assembly
    var file = e.target.files[0];
    if (!file) {
        return;
    }
    var reader = new FileReader();
    reader.onload = function(e) {
        var contents = e.target.result;
        sequencingInput(e,contents);
    };
    reader.readAsText(file);

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

    if (!dna) {
        return;
    }
   
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

    var dna,k,res,nodes,resArray,eqCount,neButSameLengthCount, splFloat, dist, resStr, debGraph, pairDist;

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
                splFloat = spl.map(function(el2) {
                        return parseFloat(el2);

                    }
                );
                return splFloat;

            });
            dna = par2El.value;
            k = parseInt(par3El.value);

            res = profileMostProbable(profMat,dna,k);
            resEl.value = 'Best kmer: \n' + res[0] + '\nBest prob:\n' + res[1];


            break;
        case '3': //Hamilton path
            k = parseInt(par2El.value);
            par1El.value = cleanContents(par1El.value);
            dna = par1El.value;

            resEl.value = 'results\n';
            var hamGraph = new DGraph(dna,DGraph.fromDna,DGraph.hamGraph,k,false);

            resEl.value+= hamGraph.sumInfo();
            resEl.value += '\n' + hamGraph.getAdjList();


            nodes = createHamNodes(dna,null,k);
            resArray = [];

           // hamiltonCanvas(nodes);

            for (var i = 0;i < parseInt(par3El.value);++i) {
                resArray.push(hamPath(nodes));

            }

            resEl.value+= '\nresults';
            eqCount = 0;
            neButSameLengthCount = 0;

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
            k = parseInt(par2El.value);
            par1El.value = cleanContents(par1El.value);
            dna = par1El.value;
            //var nodes = createDebNodes(dna,k);
             debGraph = new DGraph(dna,DGraph.fromDna,DGraph.debGraph,k,true);



            resArray = [];
            for (var i = 0;i < parseInt(par3El.value);++i) {
                //resArray.push(debPath(nodes));
             //   resArray.push(debGraph.debPath());
                debGraph.debCycle();

                debruijnCanvas(debGraph);

                resArray.push([[],debGraph.edgePathReconstructed()]);
            }
            resEl.value = 'results\n';
            eqCount = 0;
            neButSameLengthCount = 0;

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
            k = parseInt(par2El.value);
            par1El.value = cleanContents(par1El.value);
            dna = par1El.value;
            dist = parseInt(par3El.value);
            res = kmerPairedComposition(dna,k,dist);

            resEl.value = 'results: ';
            res.forEach(function(el) {
                resEl.value += '\n' + el[0] + ':' +  el[1];
            });
            //resEl.value = 'results: ' + res[0][0];

            break;

        case '6': //deb path from paired dna
            k = parseInt(par2El.value);
            par1El.value = cleanContents(par1El.value);
            dna = par1El.value;
            pairDist = parseInt(par3El.value);
            //var nodes = createDebNodes(dna,k);
            debGraph = new DGraph(dna,DGraph.fromPairedDna,DGraph.debGraph,k,false,pairDist);


            resArray = [];
            for (var i = 0;i < parseInt(par3El.value);++i) {
                //resArray.push(debPath(nodes));
                //   resArray.push(debGraph.debPath());
                debGraph.debPath();

                debruijnCanvas(debGraph);

                resArray.push([[],debGraph.edgePathReconstructedPairs()]);
            }
            resEl.value = 'results\n';
            eqCount = 0;
            neButSameLengthCount = 0;

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

            k = reads[0].length;

            debGraph = new DGraph(reads,DGraph.fromReads,DGraph.debGraph,k,false);

            resEl.value = 'results\n';

            resArray = [];
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

            eqCount = 0;
            neButSameLengthCount = 0;

            // resEl.value += 'rec: ' + debGraph.edgePathReconstructed();

            resArray.forEach(function(el,i) {


                //resArray.forEach(function(el,i) {
                //resEl.value += el[1];
                resEl.value += '\n' + el[1] + ' ' + i;

            });

            break;


        case '8': //debruijn path from dna
            k = parseInt(par2El.value);
            par1El.value = cleanContents(par1El.value);
            dna = par1El.value;
            //var nodes = createDebNodes(dna,k);

            debGraph = new DGraph(dna,DGraph.fromDna,DGraph.debGraph,k,false);

        resEl.value = 'results\n';

        resEl.value+= debGraph.sumInfo();

          resArray = [];
            for (var i = 0;i < parseInt(par3El.value);++i) {
                //resArray.push(debPath(nodes));
                //   resArray.push(debGraph.debPath());

                debGraph.debPath();

                resStr = debGraph.getAdjList();
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



            eqCount = 0;
            neButSameLengthCount = 0;

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
            k = parseInt(par2El.value);
            //par1El.value = cleanContents(par1El.value);
           // var dna = par1El.value;
            reads = par1El.value.split('\n');

            reads = reads.map(function(el) {
                var pair = el.split('|');

                return pair;
            });

            pairDist = parseInt(par3El.value);
            //var nodes = createDebNodes(dna,k);
            debGraph = new DGraph(reads,DGraph.fromPairedReads,DGraph.debGraph,k,false,pairDist);


            resArray = [];
            for (var i = 0;i < 1;++i) {
                //resArray.push(debPath(nodes));
                //   resArray.push(debGraph.debPath());
                debGraph.debPath();

                debruijnCanvas(debGraph);

                resArray.push([[],debGraph.edgePathReconstructedPairs()]);
            }
            resEl.value = 'results\n';
            eqCount = 0;
            neButSameLengthCount = 0;

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
            k = parseInt(par2El.value);
            // var dna =         par1El.value;
            resEl.value = 'results\n';
            resEl.value+= indToKmer(ind,k);

            break;

        case '12': //d-neighbourhood
            par1El.value = cleanContents(par1El.value);

            var d = parseInt(par2El.value);

            // var dna =         par1El.value;
            resEl.value = 'results\n';
            res = kMersWithMaxDist(par1El.value,d);
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
            res = patternSequencesDist(par1El.value,strs);

            resEl.value+= res;

            break;




        case '14': //kmer composition
            k = parseInt(par2El.value);
            par1El.value = cleanContents(par1El.value);
            dna = par1El.value;
            var dist = parseInt(par3El.value);
            res = kmerComposition(dna,k);

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

            res = stringFromGenomePath(path);

            resEl.value = 'results: ';

            resEl.value += res;
            //resEl.value = 'results: ' + res[0][0];

            break;

        case '16': //overlap graph

            //par1El.value = cleanContents(par1El.value);
            reads = par1El.value.split('\n');
            k = reads[0].length;
            nodes = createHamNodes(null,reads,k);
            resArray = [];

            // hamiltonCanvas(nodes);

            res = hamPath(nodes);

            resStr = '';

            res[2].forEach(function(node) {
                resStr+='\n'  + node.dna;
                resStr+= ' -> ';
                node.successors.forEach(function(suc) {
                    resStr += suc.dna + ' ';

                });

            });
            resEl.value = 'results';
            resEl.value+= resStr;
            eqCount = 0;
            neButSameLengthCount = 0;

            var el = res;

            if (nodes.length < 30) {
                hamiltonCanvas(nodes, el[2]);
            }

            break;
        case '17':
            k = parseInt(par2El.value);
            //par1El.value = cleanContents(par1El.value);
            // var dna = par1El.value;
            var adjList = par1El.value.split('\n');

            //var nodes = createDebNodes(dna,k);
            debGraph = new DGraph(adjList,DGraph.fromAdjList,DGraph.debGraph,k,false);
            resEl.value = 'results\n';

            resArray = [];
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

            eqCount = 0;
            neButSameLengthCount = 0;

            // resEl.value += 'rec: ' + debGraph.edgePathReconstructed();

            resArray.forEach(function(el,i) {


                //resArray.forEach(function(el,i) {
                //resEl.value += el[1];
               // resEl.value += '\n' + el[1] + ' ' + i;
            });

            break;

        case '18':
            k = parseInt(par2El.value);
            var adjList = par1El.value.split('\n');

            debGraph = new DGraph(adjList,DGraph.fromAdjList,DGraph.debGraph,k,false);
            resEl.value = 'results\n';

            resArray = [];
            for (var i = 0;i < 1;++i) {
               debGraph.debPath();
                resEl.value += debGraph.edgePathToText();

                debruijnCanvas(debGraph);

                resArray.push([[],debGraph.edgePathReconstructedPairs()]);
            }


            break;

        case '19':
            k = parseInt(par1El.value);

            resEl.value = 'results\n';

            var pad = "00000000000000";

            reads = [];
            for (var i = 0;i < Math.pow(2,k);++i) {
                var binStr = (i >>> 0).toString(2);
                binStr = pad + binStr;
                binStr = binStr.substr(binStr.length - k);
                reads.push(binStr);

               // resEl.value +=  (i >>> 0).toString(2) + ' ';
                resEl.value+= binStr + '\n';
           }


            debGraph = new DGraph(reads,DGraph.fromReads,DGraph.debGraph,k,false);

            resEl.value+= 'results\n';

            resArray = [];
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

            eqCount = 0;
            neButSameLengthCount = 0;

            // resEl.value += 'rec: ' + debGraph.edgePathReconstructed();

            resArray.forEach(function(el,i) {


                //resArray.forEach(function(el,i) {
                //resEl.value += el[1];
                resEl.value += '\n' + el[1] + ' ' + i;

            });

            break;

        case '20': //Genome path to string

            var path = par1El.value.split('\n');

            str = '';
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

            dna = new DNA(par1El.value);

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

            resEl.value = 'len: ' + pep.peptide.length + '\n';

            resEl.value +=pep.toShortString('');
            break;

        case '23': //Find Amino Acid Sequence in genome

            dna = par1El.value;
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

            res = allPossiblePeptidesWithWeight(mass);

            resEl.value +=  res;
            break;

        case '26': //Cyclopeptide sequencing

           // var paramObj = initialiseMotifParams();

            var experimentalSpectrum = par1El.value;


            w.postMessage({
            'task': 'seqCyclopeptide',
            'spectrum': experimentalSpectrum
    }); // Start the worker.



            //resEl.value = 'results\n';

           // var res = cycloPeptideSequencing(experimentalSpectrum);

           // resEl.value +=  res;

            break;

        case '27': //Cyclopeptide scoring

            var pep = new Peptide(Peptide.AminoArrFromStr(par1El.value));

            var expSpec = par2El.value.split(' ');
            expSpec = expSpec.map(function(el) {
                return parseInt(el);
            });

            useLinear = parseInt(par3El.value);

            resEl.value = 'results\n';

            res = pep.score(expSpec,useLinear);

            resEl.value +=  res;
            break;

        case '28': //Leaderboard Cyclopeptide sequencing

            var expSpec = par1El.value; //.split(' ');

            var N = parseInt(par2El.value);

            var M = parseInt(par3El.value);

           // expSpec = expSpec.map(function(el) {
           //     return parseInt(el);
           // });

            w.postMessage({
                'task': 'seqLeaderboardCyclopeptide',
                'spectrum': expSpec,
                'M':M,
                'N': N
            }); // Start the worker.

            //var res =leaderboardCyclopeptideSequencing(expSpec,N);

            resEl.value = 'results\n';
            resEl.value+= res.toShortString('') + '\n';

            resEl.value +=  res.toWeightString('');
            break;

        case '29': //Trim leaderboard

            var leaders = par1El.value.split(' ');

            leaders = leaders.map(function(el) {
                return new Peptide(Peptide.AminoArrFromStr(el));

            });
            var spec = par2El.value.split(' ');
            spec = spec.map(function(el) {
                return parseInt(el);
            });

            var N = parseInt(par3El.value);

            res =  trimLeaderboard(leaders,spec,N);



            resEl.value = 'results\n';

            var str = '';

            res.forEach(function(el) {
                str+= el.toShortString('') + ' ';

            });
            resEl.value += str;


            break;

        case '30': //Convolution

            var spec = new Spectrum(par1El.value);
            res =  spec.convolutionStr();

            resEl.value = 'results\n';


            resEl.value += res;


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

    var inclRevCompl = false;//document.getElementById('includeRevComplMS').checked;

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


function initialiseTransParams() {

    var paramObj = {};


    var numKmer = document.getElementById('numKmers');
    numKmer.value = 1;


    paramObj.debug = document.getElementById('debugTT').checked;

    paramObj.useStartCodon = document.getElementById('startCodTT').checked;
    paramObj.useStopCodon = document.getElementById('stopCodTT').checked;
    paramObj.useAllReadingFrames = document.getElementById('readFrameTT').checked;
    paramObj.readingFrameOffset = parseInt(document.getElementById('readFrameOffsetTT').value);
    paramObj.inclRevCompl = document.getElementById('revComplTT').checked;

    var rads = document.getElementsByName('trMethod');

    var methodSelected;

    for (var i = 0;i < rads.length;++i) {
        if (rads[i].checked) {
            methodSelected = rads[i];
            break;

        }
    }

    switch (methodSelected.id) {
        case 'trTranscribe':
            paramObj.trMethod = DNA.TransMethodTranscribe;
            break;
        case 'trTranslate':
            paramObj.trMethod = DNA.TransMethodTranslate;
            break;
        case 'trRetro':
            paramObj.trMethod = DNA.TransMethodRetro;
            break;
        default:
            break;


    }

    rads = document.getElementsByName('trInput');

    var inputSelected;

    for (var i = 0;i < rads.length;++i) {
        if (rads[i].checked) {
            inputSelected = rads[i];
            break;

        }
    }

    switch (inputSelected.id) {
        case 'trDNA':
            paramObj.trInput = DGraph.fromDna;
            break;
        case 'trRNA':
            paramObj.trInput = DGraph.fromRna;
            break;
        case 'trProtein':
            paramObj.trInput = DGraph.fromProtein;
            break;
        default:
            break;


    }


    return paramObj;

}

function initialisePeptideParams() {

    var paramObj = {};


    var numKmer = document.getElementById('numKmers');
    numKmer.value = 1;


    paramObj.debug = document.getElementById('debugPS').checked;


    paramObj.convSize = parseInt(document.getElementById('convPS').value);
    paramObj.leaderboardSize = parseInt(document.getElementById('leaderPS').value);
    
    paramObj.useTheseAminos = document.getElementById('peptideUseTheseAminos').value;

    var rads = document.getElementsByName('pepMethod');

    var methodSelected;

    for (var i = 0;i < rads.length;++i) {
        if (rads[i].checked) {
            methodSelected = rads[i];
            break;

        }
    }

    switch (methodSelected.id) {
        case 'pepIdealSpectrum':
            paramObj.pepMethod = Peptide.PepMethodIdealSpectrum;
            break;
        case 'pepSequenceBrute':
            paramObj.pepMethod = Peptide.PepMethodSequenceBrute;
            break;
        case 'pepSequenceLeaderboard':
            paramObj.pepMethod = Peptide.PepMethodSequenceLeaderboard;
            break;
        case 'pepSequenceLeaderboardConv':
            paramObj.pepMethod = Peptide.PepMethodSequenceLeaderboardConv;
            break;


        default:
            break;


    }

    rads = document.getElementsByName('pepInput');

    var inputSelected;

    for (var i = 0;i < rads.length;++i) {
        if (rads[i].checked) {
            inputSelected = rads[i];
            break;

        }
    }

    switch (inputSelected.id) {
        case 'pepPeptide':
            paramObj.pepInput = DGraph.fromProtein
            break;
        case 'pepSpectrum':
            paramObj.pepInput = DGraph.fromSpectrum;
            break;

        default:
            break;


    }

    rads = document.getElementsByName('pepPeptideType');

    var peptideType;

    for (var i = 0;i < rads.length;++i) {
        if (rads[i].checked) {
            peptideType = rads[i];
            break;

        }
    }

    switch (peptideType.id) {
        case 'pepPeptideTypeLinear':
            paramObj.pepType = Peptide.linear;
            break;
        case 'pepPeptideTypeCircular':
            paramObj.pepType = Peptide.circular;
            break;

        default:
            break;


    }



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

   // sequencingInput();

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
        document.getElementById('debugText').value += '\nReads: ';
        var resStr = '';
        grph.reads.forEach(function (el) {
            resStr += '\n' + el;
        });
        document.getElementById('debugText').value += resStr;
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

function readBreak(e) {
    if (!dnaMasterStrings) {
         return;
    }

    strArr =  dnaMasterStrings.split('\n');
    brokenArr = [];
    brokenDict = {};

    strArr.forEach(function(el,i) {
        var pref = el.substring(0,el.length - 1);
        var suf = el.substring(1,el.length);
        if (pref in brokenDict) {
            brokenDict[pref][0].push(i);
        }
        else {
            var dictEntry = [];
            dictEntry.push([i]);
            dictEntry.push([]);
            brokenDict[pref] = dictEntry;
        }
        if (suf in brokenDict) {
            brokenDict[suf][1].push(i);
        }
        else {
            var dictEntry = [];
            dictEntry.push([]);
            dictEntry.push([i]);
            brokenDict[suf] = dictEntry;
        }


        brokenArr.push(el.substring(0,el.length - 1));
        brokenArr.push(el.substring(1,el.length ));
    });

    brokenArr.sort();

    var key;

    var removeDupArr = [];

    for (key in brokenDict) {
        var dupCount = 0;
        var sameCount = 0; //pref and suf same in single item. need to retain
        var curr = brokenDict[key];

        var count = 0;

       // var matchFound = false;
       // var nonMatchFound = false;
        for (var p = 0;p < curr[0].length;++p) {
            var pref = curr[0][p];
            // curr[0].forEach(function(pref) {
            var ind = curr[1].indexOf(pref);
            if (ind > -1) {
                // matchFound = true;
                curr[0][p] = -1;
                curr[1][ind] = -1;
                ++count;

            }
            else {
                for (var sufInd = 0; sufInd < curr[1].length; ++sufInd) {
                    if (curr[1][sufInd] != pref) {
                        curr[1][sufInd] = -1;
                        curr[0][p] = -1;
                        ++count;
                        break;
                    }

                }
            }

        }

       // curr[1].forEach(function(suf) {

        for (var s = 0;s < curr[1].length;++s) {
            var suf = curr[1][s];
            if (suf == -1) {
                continue;
            }

            var ind = curr[0].indexOf(suf);
            if (ind > -1) {
                // matchFound = true;
                curr[0][ind] = -1;
                curr[1][s] = -1;
                ++count;

            }
            else {
                for (var prefInd = 0; prefInd < curr[0].length; ++prefInd) {
                    if (curr[0][prefInd] != suf) {
                        curr[1][s] = -1;
                        curr[0][prefInd] = -1;
                        ++count;
                        break;
                    }

                }
            }
        }

        curr[0].forEach(function(el) {
            if (el == -1) {

            }
            else {
                ++count;
            }

        });
        curr[1].forEach(function(el) {
            if (el == -1) {

            }
            else {
                ++count;
            }

        });

        for (var i = 0;i < count;++i) {
            removeDupArr.push(key);
        }




    }


    dnaMasterStrings = removeDupArr.join('\n');
    dnaMasterChanged();



}

function pairedToReads(e) {
    if (!dnaMasterStrings) {
        return;
    }

    outArr = [];
    strArr =  dnaMasterStrings.split('\n');
    strArr.forEach(function(el) {
        var split = el.split('|');
        outArr.push(split[0]);
        outArr.push(split[1]);
    });

    dnaMasterStrings = outArr.join('\n');
    dnaMasterChanged();



}




function runTrans(e) {

    transInput();

    initialiseResults();

    paramObj = initialiseTransParams();


    //var input = dnaMaster;


    switch (paramObj.trMethod) {
        case    DNA.TransMethodTranscribe:
            var dnaObj;
            if (paramObj.inclRevCompl) {
                dnaObj =  new DNA(dnaMaster).complement(); //new DNA(dnaMaster);
            }
            else {
                dnaObj =  new DNA(dnaMaster);
            }

            var r = dnaObj.rnaTranscript();
            rnaMaster = new RNA(r);
            var prot = rnaMaster.translate(paramObj.readingFrameOffset,!paramObj.useStartCodon,!paramObj.useStopCodon,paramObj.inclRevCompl);
            proteinMaster = prot;
            break;
        case DNA.TransMethodTranslate:
            //rnaMaster = new RNA(dnaMaster);

            var prot = rnaMaster.translate(paramObj.readingFrameOffset,!paramObj.useStartCodon,!paramObj.useStopCodon,paramObj.inclRevCompl);
            proteinMaster = prot;
            break;
        case DNA.TransMethodRetro:
            break;
        default:
            break;

    }



    // graphChanged(grph);


    colourDNA(null, null, false);


    if (document.getElementById('debugTT').checked) {
        var resStr = '';
        resStr += '\nDebug pressed';
        resStr += '\nDNA: ' + dnaObj.dna;
        resStr += '\nRNA: ' + r;

        document.getElementById('debugText').value = resStr;
    }


}

function runPeptide(e) {



    peptideInput();

    initialiseResults();

    paramObj = initialisePeptideParams();


    //var input = dnaMaster;

    var spec;

    switch (paramObj.pepMethod) {

        case Peptide.PepMethodIdealSpectrum:

            //var prot = rnaMaster.translate(paramObj.readingFrameOffset,!paramObj.useStartCodon,!paramObj.useStopCodon,paramObj.inclRevCompl);
            spectrumMaster = proteinMaster.spectrum(paramObj.pepType == Peptide.circular,true);
            colourDNA(null, null, false);

            if (document.getElementById('debugPS').checked) {
                var resStr = '';
                resStr += '\nDebug pressed';
                resStr += '\nProtein: ' + proteinMaster.toMedString();
                resStr += '\nSpectrum length: ' + spectrumMaster.length;
                resStr += '\nSpectrum: ' + proteinMaster.spectrum(paramObj.pepType == Peptide.circular,false);



                document.getElementById('debugText').value = resStr;
            }

            break;
        case Peptide.PepMethodSequenceBrute:

        w.postMessage({
            'task': 'seqCyclopeptide',
            'spectrum': spectrumMaster,
            'pepType': paramObj.pepType
        }); // Start the worker.

        //var prot = rnaMaster.translate(paramObj.readingFrameOffset,!paramObj.useStartCodon,!paramObj.useStopCodon,paramObj.inclRevCompl);

        break;
        case Peptide.PepMethodSequenceLeaderboard:

            w.postMessage({
                'task': 'seqLeaderboardCyclopeptide',
                'spectrum': spectrumMaster,
                'M':paramObj.convSize,
                'N': paramObj.leaderboardSize,
                'pepType':paramObj.pepType
            }); // Start the worker.

            //var prot = rnaMaster.translate(paramObj.readingFrameOffset,!paramObj.useStartCodon,!paramObj.useStopCodon,paramObj.inclRevCompl);

            break;
        case Peptide.PepMethodSequenceLeaderboardConv:

            w.postMessage({
                'task': 'seqLeaderboardConvCyclopeptide',
                'spectrum': spectrumMaster,
                'useTheseAminos':paramObj.useTheseAminos,
                'M':paramObj.convSize,
                'N': paramObj.leaderboardSize,
                'pepType':paramObj.pepType
            }); // Start the worker.

            //var prot = rnaMaster.translate(paramObj.readingFrameOffset,!paramObj.useStartCodon,!paramObj.useStopCodon,paramObj.inclRevCompl);

            break;

        default:
            break;

    }



    // graphChanged(grph);






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
    if (!dna) {
        return '';
    }


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

    var ctx;
    if (!skewData) {
        ctx = c.getContext("2d");
        ctx.clearRect(0,0,w,h);
        return;
    }

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


    ctx = c.getContext("2d");
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
    if (el) {
        el.parentNode.removeChild(el);
    }
   // document.getElementById('rightOne').appendChild(t);

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

function createSpan(text,spanClass) {
    return '<span class = "' + spanClass + '">' + text + '</span>';
}

function colourTrans() {
    var trans = document.getElementById('transViewer');
    var gap = 10;


    var revCompl = document.getElementById('revComplTT').checked;

    var startPos = -1;
    var endPos = -1;

    var codonNum = 0;

    var dnaTit = 'DNA&nbsp;';
    var revTit = '&nbsp;&nbsp;&nbsp;&nbsp;';
    var rnaTit = 'RNA&nbsp;';
    var proTit = 'PROT';
    var headTit = '&nbsp;&nbsp;&nbsp;&nbsp;';


    var len;
    if (dnaMaster) {
        len = dnaMaster.length;
    }
    else if (rnaMaster) {
        len = rnaMaster.rna.length;
    }
    else {
        len = 0;
    }

    if (proteinMaster) {
        startPos = proteinMaster.startPosInRNA;
        endPos = proteinMaster.peptide.length * Codon.len + startPos -1;
    }

    if (revCompl) {
        var svdSt = startPos;
        startPos =  len - 1 - endPos;
        endPos = len - 1 - svdSt;

    }

    var prefix;
    var protStr = '';
    var shortProtStr = '';

    var headPrimePrefix =  revCompl ? '<<&nbsp;' : '>>&nbsp;';
    var headPrimeSuffix =  revCompl ? '&nbsp;<<' : '&nbsp;>>';
    var primePrefix = '&nbsp;&nbsp;&nbsp;';

    var fivePrime = "5'";
    var threePrime = "3'";
    var sp = "&nbsp;";


    var startPrime = revCompl ? threePrime + sp : fivePrime + sp;
    var endPrime = revCompl ? sp + fivePrime : sp + threePrime;


    var d;
    var rd;
    if (dnaMaster) {
        d = new DNA(dnaMaster);
        rd = d.complement();
    }
    else {
        d = new DNA('');
        rd = d.complement();
    }


    if (proteinMaster) {


        var i;
        prefix = '';
        if (revCompl) {
            for (i = 0; i < startPos; ++i) {
                prefix += '&nbsp;';
            }

        }
        else {
            for (i = 0; i < proteinMaster.startPosInRNA; ++i) {
                prefix += '&nbsp;';
            }

        }

        if (revCompl) {
            protStr = proTit + headPrimePrefix + prefix + proteinMaster.toString(Peptide.Medium,'',revCompl) + headPrimeSuffix;
        }
        else {
            protStr = proTit + headPrimePrefix + prefix + proteinMaster.toMedString('') + headPrimeSuffix;
        }

        var singles;
        if (revCompl) {
            singles = proteinMaster.toString(Peptide.Short,'',revCompl);
        }
        else {
            singles = proteinMaster.toShortString('');
        }

        for (var ch = 0;ch < singles.length;++ch) {

            shortProtStr += singles[ch] + '&nbsp;&nbsp';
        }
        shortProtStr = proTit + primePrefix + prefix + shortProtStr;

        if (singles.length == 0) {
            protStr = proTit + headPrimePrefix + '&nbsp;&nbsp;&nbsp;&nbsp;No protein';
            shortProtStr = '';
        }

    }


   var markStr = '|1';



  // if (rnaMaster) {
       for (var i = 1; i < len; ++i) {


           var pos = i + 1;
           if (pos % gap == 0) {
               var now = pos.toString();
               /*
               if (startPos == i) {
                   markStr += '[' + now;
               }
               else if (i == endPos) {
                   markStr += 'v' + now;
               }
               else {
               */
                   markStr += '|' + now;
              // }
               i += now.length;
           }
           else {
               /*
               if (startPos == i) {
                   markStr += 'v';
               }
               else if (i == endPos) {
                   markStr += 'v';

               }

               else {
               */
                   markStr += '&nbsp;';
              // }
           }


       }
  // }

    if (len > 10) {
        markStr = headTit + primePrefix + markStr;
    }
    else {
        markStr = '';

    }





    var rnaStr = '';




    if (rnaMaster) {
       // var rna = revCompl ?  rnaMaster.rna.split("").reverse().join("") : rnaMaster.rna;
        var rna = rnaMaster.rna;
        var rnaColouredStr = '';
        for (var i = 0;i < rna.length;++i) {


            if ((i < startPos) || i > endPos) {
                rnaColouredStr += rna[i];
            }
            else {
                var codonNum = Math.floor((i - startPos) / Codon.len);
                var txt = createSpan(rna[i], codonNum % 2 == 0 ? 'evenCodon' : 'oddCodon');
                rnaColouredStr += txt;
            }
        }
        //rnaStr = rnaMaster.rna;
        rnaStr = rnaColouredStr;



        rnaStr = rnaTit + primePrefix +  rnaStr;
    }

    var dnaPr = '';
    var dnaPrRev = '';

    if (dnaMaster) {
        dnaPr = dnaMaster;

        if (revCompl) {
           // dnaPrRev = revTit + startPrime + dnaPr.split("").reverse().join("") + endPrime;
            dnaPrRev = revTit + startPrime + rd.dna + endPrime;
            dnaPr = dnaTit + fivePrime + sp  + dnaPr + sp + threePrime;
        }
        else {
            dnaPr = dnaTit + startPrime + dnaPr + endPrime;
        }
    }

    var revDNA = '';
    if (revCompl) {
       revDNA = dnaPrRev + '<BR>';

    }
        trans.innerHTML = markStr + '<BR>' +  dnaPr + '<BR>' + revDNA + rnaStr + '<BR>' + protStr +  '<BR>' + shortProtStr;
}

function colourPep() {
    var trans = document.getElementById('transViewer');

    trans.innerHTML = '';


    if (progExtraInfo) {
        var tab = 34;
        var padStr = '----------------------------------';;
        trans.innerHTML += 'Top 5 sub-candidates:';

        var lenBit = '';
        if (progExtraInfo[0][0]) {
            lenBit = progExtraInfo[0][0].split(' ')[0].length;
            if (lenBit < 10) {
                lenBit = '0' + lenBit;
            }
        }
       // trans.innerHTML += ' (rnd ' + lenBit + ')';

        for (var i = 0;i < 14;++i) {
            trans.innerHTML += '&nbsp';
        }
        trans.innerHTML += 'Best 5 complete matches:' + '<BR>';
        var numLeadersToDisp = 5;
        var leaderStr = '';
        for (var i = 0;i < numLeadersToDisp;++i) {
            if (progExtraInfo[0][i]) {
                leaderStr += progExtraInfo[0][i];
                leaderStr += padStr.substring(0,tab - progExtraInfo[0][i].length);
            }
            else {
                leaderStr += '----------------------------------';
            }
            if (progExtraInfo[1][i]) {
                leaderStr += ' ' + progExtraInfo[1][i];
            }
            else {
                leaderStr += ' ' + '----------------------------------';
            }

            leaderStr += '<BR>';
        }


            /*
             var bestStr = '';
             progExtraInfo[1].forEach(function(el,i) {
             if (i < numLeadersToDisp) {
             bestStr += el + '<BR>';
             }
             });
             */
        trans.innerHTML += leaderStr;
       // trans.innerHTML += 'best: ' + bestStr;
    }

    if (proteinMaster) {
        trans.innerHTML += proteinMaster.toMedString();
    }
    if ((proteinMaster) && (spectrumMaster)) {
        //var spec = proteinMaster.spectrum(document.getElementById('pepPeptideTypeCircular').checked,true);

        var specStr = '';
        spectrumMaster.forEach(function(el) {
            specStr += el[0] + '[' + el[1] + ']' + '&nbsp;';

        });

        trans.innerHTML += '<BR>' + 'Ideal Spectrum: ' + specStr;


    }
    else {
        if (spectrumMaster) {
            trans.innerHTML += '<BR>' + 'Exp Spectrum: ' + spectrumMaster;
        }
        if (convMaster) {
            trans.innerHTML += '<BR>' + 'Top M from Conv: ' + convMaster;
        }
        if (sequencedWeights) {
            var sw = sequencedWeights.split(' ');
            sw = sw.filter(function(el) {
                if (el.length == 0) {
                    return false;

                }
                else {
                    return true;
                }

            });
            var cycNums = [];
            var currCycNum = 1;
            var seqPeps = sw.map(function(el) {
               cycNums.push(0);
               return  new Peptide(Peptide.AminoArrFromWeights(el.split('-')));
            });
            seqPeps.forEach(function(el,i) {
                if (cycNums[i] == 0) {
                    cycNums[i] = currCycNum;
                    var all = el.allCycles();
                    for (var j = i+1;j < seqPeps.length; ++j) {
                        if (all.indexOf(seqPeps[j].toShortString('')) > -1)  {
                            cycNums[j] = currCycNum;
                        }
                    }
                    ++currCycNum;
                }

            });
            var swm = sw.map(function(el,i) {
                var pep = new Peptide(Peptide.AminoArrFromWeights(el.split('-')));
                return el + ' ' + pep.toShortString('-') +  ' ' + pep.getIntegerWeight() + ' ' + cycNums[i] +  '<BR>';
            });
            var swmStr = '';

            swm.forEach(function(el) {
                swmStr+=el;

            });
            if (swm.length == 0) {
                trans.innerHTML += '<BR>' + ' No Peptides match spectrum';
            }
            else {
                trans.innerHTML += '<BR>' + swmStr;
            }
        }
        else {
            if (sequencedWeights == null) {

            }
            else {
                trans.innerHTML += '<BR>' + ' No Peptides match exp spectrum';
            }
        }
    }



}


function colourDNA(dna,mark,inclRevCompl) {

    if (myParams.tabActive == 5) {
        colourTrans();
        return;
    }
    else if (myParams.tabActive == 6) {
        colourPep();
        return;
    }

    inclRevCompl = false;
    var motifView = false;
    var assemblyView = false;

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
        case 4:
            assemblyView = true;
            break;
        case 5:
            inclRevCompl = document.getElementById("revComplTT").checked;
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

    if ((assemblyView) &&(dnaMasterStrings)) {
        var strArray = dnaMasterStrings.split('\n');
        var more = strArray.length - 1;
        if (strArray.length > 20) {
            document.getElementById('graphReadsDiv').innerHTML = strArray[0] + ' + ' + more + ' more';
        }
        else {
            document.getElementById('graphReadsDiv').innerHTML = dnaMasterStrings;
        }
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
    if (el) {
        el.parentNode.removeChild(el);
    }
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
        case 5:
            skewRequired = true;

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

    mfk = [];
    mfMotif = [];

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


    var dispMaxBranchRecon = false;

    if (dispMaxBranchRecon) {
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

        nbps.forEach(function (nbp) {
            document.getElementById('debugText').value += '\n'
            nbp.forEach(function (el, i) {
                if (i == 0) {
                    document.getElementById('debugText').value += el;
                }
                else {
                    document.getElementById('debugText').value += el.substring(el.length - 1);
                }


            });

        });


    }
    else {
        document.getElementById('debugText').value += '\nNum non branch paths: ' + nbps.length;
    }






}

function cleanContents(contents,contentType) {

    var newContents;

    var spl,parts,joined,newContents;

    if (contentType) { //new regime

        switch (contentType) {
            case DGraph.fromDna:case DGraph.fromPairedDna: //dna. Can be on multiple lines - concatenated into one
                if (contents.substring(0, 1) === '>') {

                    spl = contents.split(/\r\n|\r|\n/);
                    parts = contents.split(/\r\n|\r|\n/g);
                    parts.shift();
                    joined = parts.join('\n');
                    parts.shift(); // removes the first item from the array
                    contents = parts.join('_');
                }


                newContents = contents.replace(/[^ACGT0123456789acgt]/gm, "");
                newContents = newContents.toUpperCase();
                return newContents;

                break;


            case DGraph.fromReads: //reads on multiple lines
                parts = contents.split(/\r\n|\r|\n/g);
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
                parts = contents.split(/\r\n|\r|\n/g);
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
                parts = contents.split(/\r\n|\r|\n/g);
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

            case DGraph.fromRna: //rna. Can be on multiple lines - concatenated into one
            if (contents.substring(0, 1) === '>') {

                spl = contents.split(/\r\n|\r|\n/);
                parts = contents.split(/\r\n|\r|\n/g);
                parts.shift();
                joined = parts.join('\n');
                parts.shift(); // removes the first item from the array
                contents = parts.join('_');
            }


            newContents = contents.replace(/[^ACGU0123456789acgu]/gm, "");
            newContents = newContents.toUpperCase();
            return newContents;

            break;

            case DGraph.fromProtein: //protein. Can be on multiple lines - concatenated into one
                if (contents.substring(0, 1) === '>') {

                    spl = contents.split(/\r\n|\r|\n/);
                    parts = contents.split(/\r\n|\r|\n/g);
                    parts.shift();
                    joined = parts.join('\n');
                    parts.shift(); // removes the first item from the array
                    contents = parts.join('_');
                }


                newContents = contents.replace(/[^ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz]/gm, "");
                newContents = newContents.toUpperCase();
                return newContents;

                break;

            case DGraph.fromSpectrum: //spectrum. Can be on multiple lines - concatenated into one
                if (contents.substring(0, 1) === '>') {

                    spl = contents.split(/\r\n|\r|\n/);
                    parts = contents.split(/\r\n|\r|\n/g);
                    parts.shift();
                    joined = parts.join('\n');
                    parts.shift(); // removes the first item from the array
                    contents = parts.join('_');
                }


                newContents = contents.replace(/[,]/gm, " ");
                newContents = newContents.replace(/[^0123456789. ]/gm, "");
                newContents = newContents.trim();
                newContents = newContents.toUpperCase();
                return newContents;

                break;

            default:
                return contents;
                break;
        }

    }
    else {
        if (contents.substring(0, 1) === '>') {

            spl = contents.split(/\r\n|\r|\n/);
            parts = contents.split(/\r\n|\r|\n/g);
            parts.shift();
            joined = parts.join('\n');
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

function sequencingInput(e,input) {

    var seqInp = document.getElementsByName('seqInput');

    var val = '';
    for (var i = 0;i < seqInp.length; ++i) {
        if (seqInp[i].checked) {
            val = seqInp[i].value;
        }
    }

    var source;
    if (input) {
        source = input;

    }
    else {
        source = document.getElementById('sequencingStrings').value;

    }

    var cleaned;

        switch (val) {
          case 'seqDNA':
            cleaned =  cleanContents(source,DGraph.fromDna);
            break;
          case 'seqReads':
              cleaned =  cleanContents(source,DGraph.fromReads);
                break;

          case 'seqPairedDNA':
                cleaned =  cleanContents(source,DGraph.fromPairedDna);
                break;
          case 'seqPairedReads':
                cleaned =  cleanContents(source,DGraph.fromPairedReads);
                break;

            case 'seqAdjList':
                cleaned =  cleanContents(source,DGraph.fromAdjList);
                break;


            default:
                break;



    }


    if (input) {

    }
    else {

        document.getElementById('sequencingStrings').value = cleaned;
    }


    switch (val) {
        case 'seqDNA':case 'seqPairedDNA':
            dnaMaster = cleaned;
            dnaMasterChanged();
            document.getElementById('graphDnaDiv').innerHTML = dnaMaster;
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


    dnaMasterChanged();




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

function transInput(e) {

    var trInp = document.getElementsByName('trInput');

    var val = '';
    for (var i = 0; i < trInp.length; ++i) {
        if (trInp[i].checked) {
            val = trInp[i].value;
        }
    }
    var cleaned;

    switch (val) {
        case 'trDNA':
            cleaned = cleanContents(document.getElementById('transInput').value, DGraph.fromDna);
            break;
        case 'trRNA':
            cleaned = cleanContents(document.getElementById('transInput').value, DGraph.fromRna);
            break;

        case 'trProtein':
            cleaned = cleanContents(document.getElementById('transInput').value, DGraph.fromProtein);
            break;

        default:
            break;


    }


    document.getElementById('transInput').value = cleaned;


    var revCompl = document.getElementById('revComplTT').checked;

    switch (val) {
        case 'trDNA':

            dnaMaster =  cleaned; //revCompl? reverseComplement(cleaned) : cleaned;
            rnaMaster = null;
            proteinMaster = null;
            dnaMasterChanged();
            break;

        default:
            dnaMaster = null;
            rnaMaster = new RNA(cleaned);  // rnaMaster = new RNA(revCompl? reverseComplement(cleaned) : cleaned);
            proteinMaster = null;
            dnaMasterChanged();
            break;
    }
}

function peptideInput(e) {


    var pepInp = document.getElementsByName('pepInput');

    var val = '';
    for (var i = 0; i < pepInp.length; ++i) {
        if (pepInp[i].checked) {
            val = pepInp[i].value;
        }
    }

    var cleaned;

    cleaned = document.getElementById('pepInput').value;

    switch(val) {
        case 'pepPeptide':
            cleaned = cleanContents(cleaned,DGraph.fromProtein);
            break;

        case 'pepSpectrum':
            cleaned = cleanContents(cleaned,DGraph.fromSpectrum);
            break;

    }


    document.getElementById('pepInput').value = cleaned;

    dnaMaster =  null; //revCompl? reverseComplement(cleaned) : cleaned;
    rnaMaster = null;
    proteinMaster = null;
    spectrumMaster = null;
    sequencedWeights = null;
    convMaster = null;

    switch (val) {
        case 'pepPeptide':
            //cleaned = document.getElementById('pepInput').value.toUpperCase();
            proteinMaster = new Peptide(Peptide.AminoArrFromStr(cleaned));
            break;
        case 'pepSpectrum':
            //cleaned = document.getElementById('pepInput').value;
            spectrumMaster = cleaned;

            var spec = new Spectrum(spectrumMaster);
            var topEls = spec.topElements(parseInt(document.getElementById('convPS').value));
            convMaster = topEls;
            break;


    }

    dnaMasterChanged();


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

function transInputRadClicked(id) {
    //Motif Radio Button
    switch (id) {
        case 'trDNA':

            document.getElementById('trTranscribe').checked = true;
            transRadClicked('trTranscribe');


            break;

        case 'trRNA':


            document.getElementById('trTranslate').checked = true;
            transRadClicked('trTranslate');
            break;


        case 'trProtein':

            document.getElementById('trRetro').checked = true;
            transRadClicked('trRetro');

            break;

        default:


            break;

    }

}

function peptideInputRadClicked(id) {
    //Motif Radio Button
    switch (id) {
        case 'pepPeptide':
            document.getElementById('pepIdealSpectrum').checked = true;
            peptideRadClicked('pepIdealSpectrum');
          break;

        case 'pepSpectrum':
            document.getElementById('pepSequenceBrute').checked = true;
            peptideRadClicked('pepSequenceBrute');
            break;


        default:


            break;

    }

}


function transRadClicked(id) {
    //Motif Radio Button
    switch (id) {
        case 'trTranscribe':
            document.getElementById('startCodTT').style.display = "block";
            document.getElementById('stopCodTT').style.display = "block";
            document.getElementById('readFrameTT').style.display = "none";
            document.getElementById('readFrameOffsetTT').style.display = "block";
            document.getElementById('revComplTT').style.display = "block";

            document.getElementById('startCodLabTT').style.display = "block";
            document.getElementById('stopCodLabTT').style.display = "block";
            document.getElementById('readFrameLabTT').style.display = "none";
            document.getElementById('readFrameOffsetLabTT').style.display = "block";
            document.getElementById('revComplLabTT').style.display = "block";

            /*
            document.getElementById('trDNA').style.display = "inline-block";
            document.getElementById('trRNA').style.display = "none";
            document.getElementById('trProtein').style.display = "none";

            document.getElementById('trDNALab').style.display = "inline-block";
            document.getElementById('trRNALab').style.display = "none";
            document.getElementById('trProteinLab').style.display = "none";
            */

            document.getElementById('trDNA').checked = true;


            break;

        case 'trTranslate':

            document.getElementById('startCodTT').style.display = "block";
            document.getElementById('stopCodTT').style.display = "block";
            document.getElementById('readFrameTT').style.display = "none";
            document.getElementById('readFrameOffsetTT').style.display = "block";
            document.getElementById('revComplTT').style.display = "block";

            document.getElementById('startCodLabTT').style.display = "block";
            document.getElementById('stopCodLabTT').style.display = "block";
            document.getElementById('readFrameLabTT').style.display = "none";
            document.getElementById('readFrameOffsetLabTT').style.display = "block";
            document.getElementById('revComplLabTT').style.display = "block";

            /*
            document.getElementById('trDNA').style.display = "none";
            document.getElementById('trRNA').style.display = "inline-block";
            document.getElementById('trProtein').style.display = "none";

            document.getElementById('trDNALab').style.display = "none";
            document.getElementById('trRNALab').style.display = "inline-block";
            document.getElementById('trProteinLab').style.display = "none";
            */

            document.getElementById('trRNA').checked = true;

            break;


        case 'trRetro':

            document.getElementById('startCodTT').style.display = "none";
            document.getElementById('stopCodTT').style.display = "none";
            document.getElementById('readFrameTT').style.display = "none";
            document.getElementById('readFrameOffsetTT').style.display = "none";
            document.getElementById('revComplTT').style.display = "none";

            document.getElementById('startCodLabTT').style.display = "none";
            document.getElementById('stopCodLabTT').style.display = "none";
            document.getElementById('readFrameLabTT').style.display = "none";
            document.getElementById('readFrameOffsetLabTT').style.display = "none";
            document.getElementById('revComplLabTT').style.display = "none";

            /*
            document.getElementById('trDNA').style.display = "none";
            document.getElementById('trRNA').style.display = "none";
            document.getElementById('trProtein').style.display = "inline-block";

            document.getElementById('trDNALab').style.display = "none";
            document.getElementById('trRNALab').style.display = "none";
            document.getElementById('trProteinLab').style.display = "inline-block";
            */

            document.getElementById('trProtein').checked = true;

            alert('Protein reverse translation not yet implemented');

            break;

        default:
            document.getElementById('startCodTT').style.display = "block";
            document.getElementById('stopCodTT').style.display = "block";
            document.getElementById('readFrameTT').style.display = "none";
            document.getElementById('readFrameOffsetTT').style.display = "block";
            document.getElementById('revComplTT').style.display = "block";

            document.getElementById('startCodLabTT').style.display = "block";
            document.getElementById('stopCodLabTT').style.display = "block";
            document.getElementById('readFrameLabTT').style.display = "none";
            document.getElementById('readFrameOffsetLabTT').style.display = "block";
            document.getElementById('revComplLabTT').style.display = "block";

            /*
            document.getElementById('trDNA').style.display = "none";
            document.getElementById('trRNA').style.display = "inline-block";
            document.getElementById('trProtein').style.display = "none";

            document.getElementById('trDNALab').style.display = "none";
            document.getElementById('trRNALab').style.display = "inline-block";
            document.getElementById('trProteinLab').style.display = "none";
            */

            document.getElementById('trRNA').checked = true;

            break;

    }

}

function peptideRadClicked(id) {
    //Peptide Radio Button
    switch (id) {
        case 'pepIdealSpectrum':
            document.getElementById('convPS').style.display = "none";
            document.getElementById('leaderPS').style.display = "none";

            document.getElementById('convLabPS').style.display = "none";
            document.getElementById('leaderLabPS').style.display = "none";


            document.getElementById('pepPeptide').checked = true;


            break;

        case 'pepSequenceBrute':

            document.getElementById('convPS').style.display = "none";
            document.getElementById('leaderPS').style.display = "none";

            document.getElementById('convLabPS').style.display = "none";
            document.getElementById('leaderLabPS').style.display = "none";


            document.getElementById('pepSpectrum').checked = true;

            break;

        case 'pepSequenceLeaderboard':

            document.getElementById('convPS').style.display = "block";
            document.getElementById('leaderPS').style.display = "block";

            document.getElementById('convLabPS').style.display = "block";
            document.getElementById('leaderLabPS').style.display = "block";


            document.getElementById('pepSpectrum').checked = true;

            break;

        case 'pepSequenceLeaderboardConv':

            document.getElementById('convPS').style.display = "block";
            document.getElementById('leaderPS').style.display = "block";

            document.getElementById('convLabPS').style.display = "block";
            document.getElementById('leaderLabPS').style.display = "block";


            document.getElementById('pepSpectrum').checked = true;

            break;
        default:
            document.getElementById('convPS').style.display = "none";
            document.getElementById('leaderPS').style.display = "none";

            document.getElementById('convLabPS').style.display = "none";
            document.getElementById('leaderLabPS').style.display = "none";


            document.getElementById('pepPeptide').checked = true;

            break;

    }

}

function updateExpColl(elId,newState) {
    var el = document.getElementById(elId);

    var imgSrc = newState ? 'img/toggle-collapse-alt_blue.png' : 'img/toggle-expand-alt_blue.png';

    el.src = imgSrc;

    switch (elId) {
        case 'expDebug':

            break;

    }



    var debMax = '25%';
    var debMin = '2%';

    var rOneMax = '27%';
    var rOneMin = '3%';

    var leftMax = '93%';
    var leftMin = '46%';
    var leftWithDeb = '70%';
    var leftWithMulti = '69%';

   // var right = document.getElementById('rightMainDiv');
    var debDiv =  document.getElementById('debugDiv');
    var rOneDiv = document.getElementById('rightOne');
    var left = document.getElementById('leftMainDiv');

    if (expDebugState) {
        document.getElementById('debugLab').innerHTML = 'Debug Info';
    }
    else {
        document.getElementById('debugLab').innerHTML = '';
    }

    if (expMultiState) {
        document.getElementById('rightOneLab').innerHTML = 'Results Selector';
    }
    else {
        document.getElementById('rightOneLab').innerHTML = '';
    }


    if ((expDebugState) && (expMultiState)) {
       // right.style.width = rightMax;
        left.style.width = leftMin;
        debDiv.style.width = debMax;
        rOneDiv.style.width = rOneMax;
    }
    else if (expDebugState)  {
       // right.style.width = rightWithDeb;
        left.style.width = leftWithDeb;
        debDiv.style.width = debMax
        rOneDiv.style.width = rOneMin;
    }
    else if (expMultiState) {
      //  right.style.width = rightWithMulti;
        left.style.width = leftWithMulti;
        debDiv.style.width = debMin;
        rOneDiv.style.width = rOneMax;
    }
    else {
       // right.style.width = rightMin;
        left.style.width = leftMax;
        debDiv.style.width = debMin;
        rOneDiv.style.width = rOneMin;
    }



}
function expCollClicked(elId) {
   // alert('exp coll clicked: ' +  elId);
        var el = document.getElementById(elId);


    var newState;

    switch(elId) {
        case 'expDebug':
            expDebugState = !expDebugState;
            expStateChanged(elId,expDebugState);

            break;
        case 'expMulti':
            expMultiState = !expMultiState;
            expStateChanged(elId,expMultiState);
            break;
        default:
            break;
    }
}

function expStateChanged(elId,newVal) {
    //expCollClicked(elId);
    updateExpColl(elId,newVal);

}



function clickInDNA(e) {
    //alert('key down. Posn: ' + getCaretCharacterOffsetWithin(document.getElementById('dnaView')));

    var pos = getCaretCharacterOffsetWithin(document.getElementById('dnaView')) + (dnaPage * basesPerPage) + 1;
    document.getElementById('currPosInDNA').innerHTML =  pos;//getCaretCharacterOffsetWithin(document.getElementById('dnaView'));
}


function tabClickDone(tabNum) {
    //alert('tab clicked: ' + tabNum);

    initialiseResults();

    document.getElementById('rightTwo').style.display = 'block';

    switch (tabNum) {
        case 3: //Motif tab
             document.getElementById('mfkTools').style.display = "none";
             document.getElementById('motifTools').style.display = "block";
             document.getElementById('sequencingTools').style.display = "none";
             document.getElementById('transTools').style.display = "none";
             document.getElementById('pepTools').style.display = "none";

            document.getElementById('dnaviewer').style.display = "inline-block";
            document.getElementById('graphviewer').style.display = "none";
            document.getElementById('transViewer').style.display = "none";

            expDebugState = false;
            expStateChanged('expDebug',false);
            expMultiState = false;
            expStateChanged('expMulti',false);
            break;
        case 4:  //Sequencing tab
            document.getElementById('mfkTools').style.display = "none";
            document.getElementById('motifTools').style.display = "none";
            document.getElementById('sequencingTools').style.display = "block";
            document.getElementById('transTools').style.display = "none";
            document.getElementById('pepTools').style.display = "none";

            document.getElementById('dnaviewer').style.display = "none";
            document.getElementById('graphviewer').style.display = "inline-block";
            document.getElementById('transViewer').style.display = "none";

            expDebugState = false;
            expStateChanged('expDebug',false);
            expMultiState = false;
            expStateChanged('expMulti',false);
            break;

        case 5:  //Transcription/Translation tab
            document.getElementById('mfkTools').style.display = "none";
            document.getElementById('motifTools').style.display = "none";
            document.getElementById('sequencingTools').style.display = "none";
            document.getElementById('transTools').style.display = "block";
            document.getElementById('pepTools').style.display = "none";

            document.getElementById('dnaviewer').style.display = "none";
            document.getElementById('graphviewer').style.display = "none";
            document.getElementById('transViewer').style.display = "inline-block";

            expDebugState = false;
            expStateChanged('expDebug',false);
            expMultiState = false;
            expStateChanged('expMulti',false);
            break;

        case 6: //Peptide sequencing tab
            document.getElementById('mfkTools').style.display = "none";
            document.getElementById('motifTools').style.display = "none";
            document.getElementById('sequencingTools').style.display = "none";
            document.getElementById('transTools').style.display = "none";
            document.getElementById('pepTools').style.display = "block";

            document.getElementById('dnaviewer').style.display = "none";
            document.getElementById('graphviewer').style.display = "none";
            document.getElementById('transViewer').style.display = "inline-block";

            document.getElementById('rightTwo').style.display = 'none';



            expDebugState = false;
            expStateChanged('expDebug',false);
            expMultiState = false;
            expStateChanged('expMulti',false);
            break;

        case 7: //Misc tab
            document.getElementById('mfkTools').style.display = "block";
            document.getElementById('motifTools').style.display = "none";
            document.getElementById('sequencingTools').style.display = "none";
            document.getElementById('transTools').style.display = "none";
            document.getElementById('pepTools').style.display = "none";

            document.getElementById('dnaviewer').style.display = "inline-block";
            document.getElementById('graphviewer').style.display = "none";
            document.getElementById('transViewer').style.display = "none";

            expDebugState = true;
            expStateChanged('expDebug',true);
            expMultiState = false;
            expStateChanged('expMulti',false);
            break;



        default:  //0, 1, 2 kmer tabs
            document.getElementById('mfkTools').style.display = "block";
            document.getElementById('motifTools').style.display = "none";
            document.getElementById('sequencingTools').style.display = "none";
            document.getElementById('transTools').style.display = "none";
            document.getElementById('pepTools').style.display = "none";

            document.getElementById('dnaviewer').style.display = "inline-block";
            document.getElementById('graphviewer').style.display = "none";
            document.getElementById('transViewer').style.display = "none";

            expDebugState = false;
            expStateChanged('expDebug',false);
            expMultiState = true;
            expStateChanged('expMulti',true);
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
