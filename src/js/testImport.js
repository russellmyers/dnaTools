/**
 * Created by RussellM on 16/12/2015.
 */


var importedVar = 'This is imported';


function scoreMotifs(motifs) {
    //assumes equal length motifs, otherwise returns -1

    var cols = motifsToColumns(motifs);

    var score = 0;
    cols.forEach(function(col,c) {
        var counts = baseCount(col);
        var m = counts.reduce(function(a,b){
            return Math.max(a,b);

        });
        score += (col.length - m);



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
            countMatrix[r][c] = baseCount;

        });


    });

    var profileMatrix = countMatrix.map(function(el,r) {
        return el.map(function(col,c) {
            return col * 1.0  / motifs.length;

        });

    });

    var entropyMatrix  = profileMatrix.map(function(el,r) {
        return el.map(function(col,c) {
            return (col == 0) ? col : col * Math.log2(col);

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
    entropyScore *= -1;

    return [consensus,score,entropyScore,countMatrix,entropyMatrix,profileMatrix];
}
