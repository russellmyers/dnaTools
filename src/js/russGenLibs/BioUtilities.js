/**
 * Created by RussellM on 20/01/2016.
 */


//Bioinformatics routines

var c_Bases = ['A','C','G','T'];
var c_Bases = ['A','C','G','T'];
var c_NumBases = 4;
var c_BaseInds = {'A':0,'C':1,'G':2,'T':3};

var c_Aminos = ['G','A','L','M','F','W','K','Q','E','S','P','V','I','C','Y','H','R','N','D','T','X'];

//  Alignment using Manhattan type grid

function DGEdge(sourceNode,targetNode,edgeType) {

    //edgeType optional

    this.sourceNode = sourceNode;
    this.targetNode = targetNode;
    if (edgeType) {
        this.edgeType = edgeType;
    }
    else {
        this.edgeType = '';
    }

    this.visited = false;
    
    this.edgeWeight = 0;
    
    this.edgeAction = ''; //used for aligning, eg ' A' means insert, 'A ' means del, 'AA' means match

    this.priorityScore = 0; //used for determining which edge to use when equally good

    
}

DGEdge.prototype.edgeLabel = function() {
    return this.sourceNode.label + '->' + this.targetNode.label;
}

DGEdge.prototype.setWeight = function(w) {
    this.edgeWeight = w;
}

function DGFancyEdge(sNode,tNode,fLab) {
    DGEdge.apply(this,[sNode,tNode]);
    this.fancylab = fLab;
    
    this.edgeLabel = function() {
        return this.sourceNode.label + '->' + this.targetNode.label +  ' fancy: ' + this.fancylab;
    }
}

function DGNode(label,repeatNum) {
    this.label = label;

    this.visited = false;
    
    this.cycleNum = 0;

    this.successors = [];
    this.predecessors = [];

    this.longestPathToThisNode = 0;
    this.longestEdgeToThisNode = null;
    this.onSourcePath = true;

    if (repeatNum) {
        this.repeatNum = repeatNum; //for nodes with same label
    }
    else {
        this.repeatNum = 1;
    }
    


}


DGNode.prototype.getTopOrder = function() {
    return parseInt(this.label);
};

DGNode.prototype.setLongestPathToThisNode = function(l,ed) {
    this.longestPathToThisNode = l;
    this.longestEdgeToThisNode = ed;
};

DGNode.prototype.outDegree = function() {
    return this.successors.length;
};
DGNode.prototype.InDegree = function() {
    return this.predecessors.length;
};

DGNode.prototype.degree = function() {
    return this.InDegree() + this.outDegree();
};

DGNode.prototype.resetEdgesVisited = function() {
    this.successors.forEach(function(edge) {
        edge.visited = false;
    });

    this.visited = false; //node visited flag



};

DGNode.prototype.getRandomNonVisitedEdge = function() {

    var nonVisited = [];
    this.successors.forEach(function(el,i) {
        if (!(el.visited)) {
            nonVisited.push(i);
        }

    });
    if (nonVisited.length == 0) {
        return -1;
    }
    else {
        var r = getRandomInt(0, nonVisited.length - 1);
        return nonVisited[r];
    }

};


//Directed Graph routines

function DEdge(sourceNode,targetNode) {
    this.sourceNode = sourceNode;
    this.targetNode = targetNode;
    this.visited = false;

    this.walkNum = 0;

   // this.artificial = artificial; // indicates added to make cycle

    this.edgeLabel = function() {
        if (this.sourceNode.pairedDna) {
            return [this.sourceNode.dna + this.targetNode.dna.substring(this.targetNode.dna.length-1),
                this.sourceNode.pairedDna + this.targetNode.pairedDna.substring(this.targetNode.pairedDna.length-1)];

        }
        return this.sourceNode.dna + this.targetNode.dna.substring(this.targetNode.dna.length-1);
    }
}

function DNode(dna,repeatNum) {
    if (typeof(dna) == 'object') {
        this.dna = dna[0];
        this.pairedDna = dna[1];

    }
    else {
        this.dna = dna;
        this.pairedDna = null;
    }

    this.visited = false;
    
    this.label = '';

    this.successors = [];
    this.predecessors = [];

    if (repeatNum) {
        this.repeatNum = repeatNum; //for nodes with same label
    }
    else {
        this.repeatNum = 1;
    }


    this.overlapSuffixLen = function(compareNode) {
        //determine length of overlap with specified node

        // not implemented for paired reads

        var sLen = this.dna.length;
        var tLen = compareNode.dna.length;

        var done = false;


        var len = sLen;

        while (!done) {
            if (len > tLen) {
                --len;
            }
            else {
                if (this.dna.substring(sLen - len,sLen) == compareNode.dna.substring(0,len)) {
                    done = true;
                }
                else {
                    --len;

                }
            }
            if (len < 1) {
                done == true;
            }

        }

        return len;


    };

    this.overlapPrefixLen = function(compareNode) {
        //determine length of overlap with specified node

        // not implemented for paired reads

        var sLen = this.dna.length;
        var tLen = compareNode.dna.length;

        var done = false;


        var len = sLen;

        while (!done) {
            if (len > tLen) {
                --len;
            }
            else {
                if (this.dna.substring(0,len) == compareNode.dna.substring(sLen - len,sLen)) {
                    done = true;
                }
                else {
                    --len;

                }
            }
            if (len < 1) {
                done == true;
            }

        }

        return len;


    };


    this.prefix = function() {
        if (this.pairedDna) {
            return [this.dna.substring(0, this.dna.length - 1), this.pairedDna.substring(0, this.pairedDna.length - 1)];
        }
        else {
            return this.dna.substring(0, this.dna.length - 1);

        }

    };

    this.suffix  = function() {
        if (this.pairedDna) {
            return [this.dna.substring(1, this.dna.length), this.pairedDna.substring(1, this.pairedDna.length)];
        }
        else {
            return this.dna.substring(1, this.dna.length);

        }

    };




    this.outDegree = function() {
        return this.successors.length;
    };
    this.InDegree = function() {
        return this.predecessors.length;
    };

    this.degree = function() {
        return this.InDegree() + this.outDegree();
    };

    this.resetEdgesVisited = function() {
        this.successors.forEach(function(edge) {
           edge.visited = false;
        });

        this.visited = false; //node visited flag



    };

    this.getRandomNonVisitedEdge = function() {

        var nonVisited = [];
        this.successors.forEach(function(el,i) {
            if (!(el.visited)) {
                nonVisited.push(i);
            }

        });
        if (nonVisited.length == 0) {
            return -1;
        }
        else {
            var r = getRandomInt(0, nonVisited.length - 1);
            return nonVisited[r];
        }

    }


}

function DGAlignSpaceGraph(s,t,alignType,scoreMat,indelPen,mismatchPen,matchScore) {
    var args = [s, t, alignType, scoreMat, indelPen, mismatchPen, matchScore];
    DGAlignGraph.apply(this, args);

    this.lastRow = -1;

    this.middleRow =  (this.rows % 2 == 0) ? Math.floor(this.rows / 2) : Math.floor(this.rows/2) + 1; //actually 1 more than middle

    this.midNode = null;
    this.midEdge = null;
    this.bestBotScore = -99999;
    this.bestTopScore = -99999;
    this.bestScoreThroughMid = -99999;


    this.initGraph = function () {


    };
    
 
    this.quickScore = function(fullPartialFlag) {
        var curRowScores = [];
        var prevRowScores = [];

        var curRowPointers = [];

        var overallBestRowCol = [];
        var overallBestLastRow = [];
        var overallBestLastCol = [];

        var alType = (fullPartialFlag == 'F') ? DGraph.alignTypeGlobal : this.alignType;

        curRowScores.push(0);
        curRowPointers.push('*'); //start

        this.overallBestScore = 0;
        this.overallBestRowCol = [0, 0];
        this.overallBestScoreLastRow = -99999;
        this.overallBestLastRow = [0, 0];
        this.overallBestScoreLastCol = -99999;
        this.overallBestLastCol = [0, 0];



        for (var c = 1;c < this.cols;++c) {
            if (alType == DGraph.alignTypeGlobal) {
                curRowScores.push(this.indelPen * -1 * c);
                curRowPointers.push('>');
            }
            else {
                curRowScores.push(0);
                curRowPointers.push('*');
            }
        }
        
        for (var r = 1;r < this.rows;++r) {
            prevRowScores = curRowScores;
            curRowScores = [];
            curRowPointers = [];
            if (alType == DGraph.alignTypeLocal) {
                curRowScores.push(0);
                curRowPointers.push('*');
            }
            else {
                curRowScores.push(this.indelPen * -1 * r);
                curRowPointers.push('v');
            }

            for (var c = 1; c < this.cols; ++c) {


                if (this.scoreMat) {
                    diag = prevRowScores[c - 1] + this.scoreMat[this.s[c - 1]][this.t[r - 1]];
                }
                else {
                    var diagWeightNoMat;
                    if (this.s[c-1] == this.t[r-1]) {
                        diagWeightNoMat = this.matchScore;
                    }
                    else {
                        diagWeightNoMat = this.mismatchPen * -1;
                    }
                    diag = prevRowScores[c - 1] + diagWeightNoMat;

                }
                var best = diag;
                var pointer = '\\';

                var ins = prevRowScores[c] + indelPen * -1;
                if (ins > best) {
                    best = ins;
                    pointer = 'v';
                }

                var del = curRowScores[c - 1] + indelPen * -1;
                if (del > best) {
                    best = del;
                    pointer = '>';
                }

                if ((alType === DGraph.alignTypeLocal) && (best < 0)) {
                    best = 0;
                    pointer = '*';
                }

                curRowScores.push(best);
                curRowPointers.push(pointer);

                if (best >= this.overallBestScore) {
                    this.overallBestScore = best;
                    this.overallBestRowCol = [r, c];
                }

                if ((this.alignType == DGraph.alignTypeFitting) && (r == this.rows - 1  )) {
                        //last row
                        if (best >= this.overallBestScoreLastRow) {
                            this.overallBestLastRow = [r, c];
                            this.overallBestScoreLastRow = best;
                        }
                }
                if ((this.alignType == DGraph.alignTypeOverlap) && (c == this.cols - 1)) {  //last col

                        if (r == 0) { //ie 1st row - don't count. Need at least something to align
                        }
                        else {
                            if (best >= this.overallBestScoreLastCol) {
                                this.overallBestLastCol = [r, c];
                                this.overallBestScoreLastCol = best;
                            }
                        }
                }



            }
        }

        this.bestNode = new DGNode(this.getTopOrder(this.overallBestRowCol[0],this.overallBestRowCol[1]));
        this.bestNode.longestPathToThisNode = this.overallBestScore;

        if (this.overallBestLastRow) {
            this.bestNodeLastRow = new DGNode(this.getTopOrder(this.overallBestLastRow[0], this.overallBestLastRow[1]));
            this.bestNodeLastRow.longestPathToThisNode = this.overallBestScoreLastRow;
        }

        if (this.overallBestLastCol) {
            this.bestNodeLastCol = new DGNode(this.getTopOrder(this.overallBestLastCol[0], this.overallBestLastCol[1]));
            this.bestNodeLastCol.longestPathToThisNode = this.overallBestScoreLastCol;
        }


        return [curRowScores,curRowPointers,prevRowScores];
    };
 

    this.longestPathsDynamicInitRow = function () {
        // creates nodes and sets weights

        var topOrder;

        var progThreshold = 100000;


        var lim = this.cols * this.rows;

        //  for (var i = 0;i < this.rows; ++ i) {
        var i = this.lastRow;
        for (var j = 0; j < this.cols; ++j) {
            topOrder = (i * this.cols) + j;

            if (this.progressCallback) {
                if ((i  % progThreshold == 0) || (i  == lim - 1)) {
                    this.progressCallback('Init graph', i, lim - 1, 'align');
                }            }


            var node = new DGNode('' + topOrder);
            var matches;
            if (i == 0) {

            }
            else {
                var predUp = (i - 1) * this.cols + j;
                matches = this.findNodes('' + predUp);
                var pred = matches[0];
                var ed = new DGEdge(pred, node);
                // ed.setWeight(this.downWeights[i-1][j]);
                this.setEdgeActionAndWeight(i, j, ed, 'v');
                pred.successors.push(ed);
                node.predecessors.push(ed);
            }

            if (j == 0) {

            }
            else {
                var predLeft = i * this.cols + j - 1;
                matches = this.findNodes('' + predLeft);
                pred = matches[0];
                ed = new DGEdge(pred, node);
                // ed.setWeight(this.rightWeights[i][j-1]);
                this.setEdgeActionAndWeight(i, j, ed, '>');
                pred.successors.push(ed);
                node.predecessors.push(ed);

            }

            if ((i == 0) || (j == 0)) {

            }
            else {
                //diag


                //  if ((this instanceof DGAlignGraph)  || (this.diagWeights)) {}


                // if (this.diagWeights) {
                var predDiag = (i - 1) * this.cols + j - 1;
                matches = this.findNodes('' + predDiag);
                pred = matches[0];
                ed = new DGEdge(pred, node);
                //ed.setWeight(this.diagWeights[i - 1][j - 1]);
                this.setEdgeActionAndWeight(i, j, ed, '\\');

                pred.successors.push(ed);
                node.predecessors.push(ed);
                //   }
            }


            //this.nodes.push(node);
            this.nodes['' + topOrder] = node;

        }
        //}
		
	


    };

    this.longestPathsDynamicLongPathRow = function () {
        // calcs longest paths for a row

        var start = this.lastRow * this.cols;
        var end = (this.lastRow +1) * this.cols;

        this.longestPathsDynamic(null,null,start,end);
    };

    this.delRow = function() {

        var delRow = this.lastRow - 1;
        if (delRow >= 0) {
            var startDel = delRow * this.cols;
            var endDel = (delRow + 1) * this.cols;
            for (var i = startDel;i < endDel;++i) {
                if (this.nodes[i].longestEdgeToThisNode) {
                    this.nodes[i].longestEdgeToThisNode.sourceNode = null;
                    this.nodes[i].longestEdgeToThisNode.targetNode = null;
                    this.nodes[i].longestEdgeToThisNode = null;
                }
                if (this.nodes[i].predecessors) {

                    for (var p = 0;p < this.nodes[i].predecessors.length;++p) {
                        this.nodes[i].predecessors[p].sourceNode = null;
                        this.nodes[i].predecessors[p].targetNode = null;
                        this.nodes[i].predecessors.splice(p,1);
                    }
                }

                this.nodes[i].predecessors = null;
                /*
                 if (this.nodes[i].successors) {

                 for (var s = 0; s < this.nodes[i].successors.length; ++s) {
                 if (this.nodes[i].successors[s] === this.nodes[i]) {
                 this.nodes[i].successors[s] = null;
                 this.nodes[i].successors.splice(s, 1);
                 break;
                 }
                 }
                 }
                 */

                this.nodes[i] = null;
                delete this.nodes[i];
            }

        }


    };

    this.longestPathsDynamicCreateRow = function() {
        this.lastRow += 1;
        this.longestPathsDynamicInitRow();

        this.longestPathsDynamicLongPathRow();

        this.delRow();

    };

    this.longestPathsDynamicCreateRows = function () {


        var progThreshold = 100;

        while ((this.lastRow + 1) < this.rows ) {

                if (this.progressCallback) {
                   if ((this.lastRow   % progThreshold == 0) || (this.lastRow  == this.rows - 1)) {
                    this.progressCallback('Processing rows ', this.lastRow, this.rows, 'align');
                   }
                }


            this.longestPathsDynamicCreateRow();
        }

        

    };

    this.findMiddleEdgeAllInOneQuick = function(fullPartialFlag) {
        var mid = this.middleRow;

        var tTop = this.t.substring(0,mid);
        var tBottom = this.t.substring(mid);
        var tBottomRev =  tBottom.split('').reverse().join('');
        var sRev = this.s.split('').reverse().join('');


        var topG = new DGAlignSpaceGraph(this.s,tTop,fullPartialFlag == 'F' ? DGraph.alignTypeGlobal : this.alignType,this.scoreMat,this.indelPen,this.mismatchPen,this.matchScore);
        topG.progressCallback = this.progressCallback;

        var bottomG = new DGAlignSpaceGraph(sRev,tBottomRev,DGraph.alignTypeGlobal,this.scoreMat,this.indelPen,this.mismatchPen,this.matchScore);
        bottomG.progressCallback = this.progressCallback;

        topG.initGraph();
        bottomG.initGraph();

       // topG.longestPathsDynamicCreateRows();
        //bottomG.longestPathsDynamicCreateRows();
        var topScoresAndPointers = topG.quickScore(fullPartialFlag);
        var bottomScoresAndPointers = bottomG.quickScore(fullPartialFlag);

        var topScores = topScoresAndPointers[0];
        var bottomScores = bottomScoresAndPointers[0];
        var topPointers = topScoresAndPointers[1];
        var topPrevRowScores = topScoresAndPointers[2];

//        var topScores = topG.rowScores(topG.rows - 1);
 //       var bottomScores = bottomG.rowScores(bottomG.rows - 1);

        bottomScores = bottomScores.reverse();

        var scoresMidRowPlusOneNodes = topScores.map(function(el,i) {
            return [el + bottomScores[i],el,bottomScores[i]];
        });

        var bestNodeScore = -99999;
        var bestInd = -1;

        var bestNodeTopScore = -99999; //top contribution
        var bestNodeBotScore = -99999; //bottom contribution

        for (var i = 0;i < this.cols;++i) {
            if (scoresMidRowPlusOneNodes[i][0] > bestNodeScore) {
                bestNodeScore = scoresMidRowPlusOneNodes[i][0];
                bestNodeTopScore = scoresMidRowPlusOneNodes[i][1];
                bestNodeBotScore = scoresMidRowPlusOneNodes[i][2];
                bestInd = i;
            }

        }

        this.bestScoreThroughMid = bestNodeScore;
        this.bestBotScore = bestNodeBotScore;

        this.midPlusOneNode = new DGNode(topG.getTopOrder(mid,bestInd));
        this.midPlusOneNode.longestPathToThisNode = topScores[bestInd];
        var predRow,predCol,edgeAction;
        if (topPointers[bestInd] == '\\') {
            predRow = mid - 1;
            predCol = bestInd - 1;
            edgeAction = topG.s[predCol] + topG.t[predRow];
            this.bestTopScore = topPrevRowScores[predCol];
        }
        else if (topPointers[bestInd] == 'v') {
            predRow = mid - 1;
            predCol = bestInd;
            edgeAction = '-' + topG.t[predRow];
            this.bestTopScore = topPrevRowScores[bestInd];

        }
        else if (topPointers[bestInd] == '>') {
            predRow = mid;
            predCol = bestInd - 1;
            edgeAction = topG.s[predCol] + '-';
            this.bestTopScore = topScores[predCol];
        }
        else if (topPointers[bestInd] == '*') {
            predRow = 0;
            predCol = 0;
            edgeAction = '**';
            this.bestTopScore = 0;
        }
        topG.longestEdgeToThisNode = new DGEdge(this.midNode,this.midPlusOneNode);
        topG.longestEdgeToThisNode.edgeAction = edgeAction;

        this.midNode = new DGNode(topG.getTopOrder(predRow,predCol));
        this.midNode.longestPathToThisNode = this.bestTopScore;
        this.midEdge = topG.longestEdgeToThisNode;


        var topOrderMidPlusOneNode = topG.getTopOrder(mid,bestInd);
        var midNode,midEdgeAction,midNodeRowCol;
       // if (topG.nodes[topOrderMidPlusOneNode].longestEdgeToThisNode) {
        if (topPointers[bestInd] != '*') {
            /*
            this.midNode = topG.nodes[topOrderMidPlusOneNode].longestEdgeToThisNode.sourceNode;
            this.midEdge = topG.nodes[topOrderMidPlusOneNode].longestEdgeToThisNode;
            this.bestTopScore = this.midNode.longestPathToThisNode; // for use in next half division - ie uses
             */                                                       // score up to mid node, not up to end of mid edge
            midEdgeAction =  this.midEdge.edgeAction;
            midNodeRowCol = this.getNodeRowCol(this.midNode.label);
            return [predRow,predCol,mid,bestInd,midEdgeAction[0],midEdgeAction[1],bestNodeScore,this.bestTopScore,bestNodeBotScore];

        }
        else {
            return [-1,-1,mid,bestInd,'-','-',-99999,-99999,bestNodeBotScore];
        }




    };


    this.findMiddleEdgeAllInOne = function(fullPartialFlag) {

        //obsolete - use findMiddleEdgeAllInOneQuick

        var mid = this.middleRow;

        var tTop = this.t.substring(0,mid);
        var tBottom = this.t.substring(mid);
        var tBottomRev =  tBottom.split('').reverse().join('');
        var sRev = this.s.split('').reverse().join('');


        var topG = new DGAlignSpaceGraph(this.s,tTop,fullPartialFlag == 'F' ? DGraph.alignTypeGlobal : this.alignType,this.scoreMat,this.indelPen,this.mismatchPen,this.matchScore);
        topG.progressCallback = this.progressCallback;

        var bottomG = new DGAlignSpaceGraph(sRev,tBottomRev,DGraph.alignTypeGlobal,this.scoreMat,this.indelPen,this.mismatchPen,this.matchScore);
        bottomG.progressCallback = this.progressCallback;

        topG.initGraph();
        bottomG.initGraph();

        topG.longestPathsDynamicCreateRows();
        bottomG.longestPathsDynamicCreateRows();

        //var midNodeData = topG.findMiddleEdge(bottomG);

        var topScores = topG.rowScores(topG.rows - 1);
        var bottomScores = bottomG.rowScores(bottomG.rows - 1);
        bottomScores = bottomScores.reverse();

        var scoresMidRowPlusOneNodes = topScores.map(function(el,i) {
            return [el + bottomScores[i],el,bottomScores[i]];
        });

        var bestNodeScore = -99999;
        var bestInd = -1;

        var bestNodeTopScore = -99999; //top contribution
        var bestNodeBotScore = -99999; //bottom contribution

        for (var i = 0;i < this.cols;++i) {
            if (scoresMidRowPlusOneNodes[i][0] > bestNodeScore) {
                bestNodeScore = scoresMidRowPlusOneNodes[i][0];
                bestNodeTopScore = scoresMidRowPlusOneNodes[i][1];
                bestNodeBotScore = scoresMidRowPlusOneNodes[i][2];
                bestInd = i;
            }

        }

        this.bestScoreThroughMid = bestNodeScore;
        this.bestBotScore = bestNodeBotScore;

        var topOrderMidPlusOneNode = topG.getTopOrder(mid,bestInd);
        var midNode,midEdgeAction,midNodeRowCol;
        if (topG.nodes[topOrderMidPlusOneNode].longestEdgeToThisNode) {
            this.midNode = topG.nodes[topOrderMidPlusOneNode].longestEdgeToThisNode.sourceNode;
            this.midEdge = topG.nodes[topOrderMidPlusOneNode].longestEdgeToThisNode;
            this.bestTopScore = this.midNode.longestPathToThisNode; // for use in next half division - ie uses
                                                                    // score up to mid node, not up to end of mid edge
            midEdgeAction =  this.midEdge.edgeAction;
            midNodeRowCol = this.getNodeRowCol(this.midNode.label);
            return [midNodeRowCol[0],midNodeRowCol[1],mid,bestInd,midEdgeAction[0],midEdgeAction[1],bestNodeScore,this.midNode.longestPathToThisNode,bestNodeBotScore];

        }
        else {
            return [-1,-1,mid,bestInd,'-','-',-99999,-99999,bestNodeBotScore];
        }




    };

    this.findMiddleEdge = function(bottomG) {

        //Obsolete - use findMiddleEdgeAllInOneQuick

        //assumes this graph contains top half. bottomG contains bottom half
        //(actually top half plus 1, so we can get predecessor to middle row as the middle node)

        var topScores = this.rowScores(this.rows - 1);
        var bottomScores = bottomG.rowScores(bottomG.rows - 1);
        bottomScores = bottomScores.reverse();

        var scoresMidRowPlusOneNodes = topScores.map(function(el,i) {
            return [el + bottomScores[i],el,bottomScores[i]];
        });

        var bestNodeScore = -99999;
        var bestInd = -1;

        var bestNodeTopScore = -99999; //top contribution
        var bestNodeBotScore = -99999; //bottom contribution
        
        for (var i = 0;i < this.cols;++i) {
            if (scoresMidRowPlusOneNodes[i][0] > bestNodeScore) {
                bestNodeScore = scoresMidRowPlusOneNodes[i][0];
                bestNodeTopScore = scoresMidRowPlusOneNodes[i][1];
                bestNodeBotScore = scoresMidRowPlusOneNodes[i][2];
                bestInd = i;
            }

        }
        
        this.bestScoreThroughMid = bestNodeScore;
        this.bestBotScore = bestNodeBotScore;
        
        var topOrderMidPlusOneNode = this.getTopOrder(this.rows -1,bestInd);
        var midNode,midEdgeAction,midNodeRowCol;
        if (this.nodes[topOrderMidPlusOneNode].longestEdgeToThisNode) {
            this.midNode = this.nodes[topOrderMidPlusOneNode].longestEdgeToThisNode.sourceNode;
            this.midEdge = this.nodes[topOrderMidPlusOneNode].longestEdgeToThisNode;
            this.bestTopScore = this.midNode.longestPathToThisNode; // for use in next half division - ie uses
                                                                    // score up to mid node, not up to end of mid edge
            midEdgeAction =  this.midEdge.edgeAction;
            midNodeRowCol = this.getNodeRowCol(this.midNode.label);
            return [midNodeRowCol[0],midNodeRowCol[1],this.rows - 1,bestInd,midEdgeAction[0],midEdgeAction[1],bestNodeScore,this.midNode.longestPathToThisNode,bestNodeBotScore];
            
        }
        else {
            return [-1,-1,this.rows - 1,bestInd,'-','-',-99999,-99999,bestNodeBotScore];
        }


        


    };

    /*
    this.findMiddleNode = function(bottomG) {
        //assumes this graph contains top half. bottomG contains bottom half
        var topScores = this.rowScores(this.rows - 1);
        var bottomScores = bottomG.rowScores(bottomG.rows - 1);
        bottomScores = bottomScores.reverse();
  
        var scoresMidRowNodes = topScores.map(function(el,i) {
            return el + bottomScores[i];
        });
        
        var bestNodeScore = -99999;
        var bestInd = -1;
        
        for (var i = 0;i < this.cols;++i) {
            if (scoresMidRowNodes[i] > bestNodeScore) {
                bestNodeScore = scoresMidRowNodes[i];
                bestInd = i;
            }
            
        }
        var topOrderMidNode = this.getTopOrder(this.rows -1,bestInd);
        var midNodePredecessor = this.nodes[topOrderMidNode].longestEdgeToThisNode.sourceNode;
        var topOrderPredecessor = this.getNodeRowCol(midNodePredecessor.label);
                
        return [this.rows - 1,bestInd,topOrderPredecessor[0],topOrderPredecessor[1],bestNodeScore];
    }
    */

}

function DGAffineAlignGraph(s,t,alignType,scoreMat,indelPen,mismatchPen,matchScore,openPen) {
    var args = [s,t,alignType,scoreMat,indelPen,mismatchPen,matchScore];
    DGAlignGraph.apply(this,args);

    this.openPen = openPen; //gap opening penalty

    this.numLevels = 3;
    this.LevelTop = 1;
    this.LevelBottom = 0;
    this.LevelMiddle = 2;

    var sinkNum  =  (this.rows) * (this.cols) * (this.numLevels) - 1;
    this.sinkLab = sinkNum.toString();

    this.parentSetEdgesAndWeights = this.setEdgeActionAndWeight;


    this.getNodeRowColLevel = function(node) {
      return [0,0,0];
    };

    this.getNodeTopOrder = function(r,c,lev) {
      return r * this.cols * this.numLevels + (c * this.numLevels) + lev;
    };


    this.setEdgeActionAndWeight = function(row,col,ed,act) {

        this.parentSetEdgesAndWeights(row,col,ed,act);

        switch (act) {
            case 'C':
                ed.edgeAction = '' + '';
                ed.setWeight(0);
                ed.priorityScore = 5;
                break;
            case 'O>': //open deletion
                ed.edgeAction = this.s[col - 1] + '-';
                ed.setWeight(this.openPen * -1);
                ed.priorityScore = 3;
                break;
            case 'Ov' : //open insertion
                ed.edgeAction =  '-' + this.t[row - 1];
                ed.setWeight(this.openPen * -1);
                ed.priorityScore = 2;
                break;
            case 'v':
                //override for first insertion - actually should have open gap not extension
               // if (row == 1) ed.setWeight(this.openPen);
                break;
            case '>':
                //override for first del
                //if (col == 1) ed.setWeight(this.openPen);
                break;
            default:



                break;

        }

    };

    this.initGraph = function() {
        var topOrder;

        var progThreshold = 100;



        var lim =  (this.cols * this.rows) * this.numLevels;



        for (var i = 0;i < this.rows; ++ i) {
            for (var j = 0; j < this.cols; ++j) {
                for (var lev = 0; lev < this.numLevels; ++lev) {

                    topOrder = (i * this.numLevels * this.cols) + (j * this.numLevels) + lev;

                    if (this.progressCallback) {
                        if ((topOrder % progThreshold == 0) || (topOrder == lim - 1)) {
                            this.progressCallback('Init graph', topOrder, lim - 1, 'align');
                        }
                    }


                    var node = new DGNode('' + topOrder);
                    if (lev == this.LevelBottom) { //bottom = insertions

                        if (i == 0) {

                        }
                        else {
                            var predUp = this.getNodeTopOrder(i - 1, j, lev);
                            matches = this.findNodes('' + predUp);
                            var pred = matches[0];
                            var ed = new DGEdge(pred, node);
                            // ed.setWeight(this.downWeights[i-1][j]);
                            this.setEdgeActionAndWeight(i, j, ed, 'v');
                            //pred.successors.push(ed);
                            node.predecessors.push(ed);

                            var predMiddle = this.getNodeTopOrder(i - 1, j, this.LevelMiddle);
                            matches = this.findNodes('' + predMiddle);
                            var pred = matches[0];
                            var ed = new DGEdge(pred, node);
                            this.setEdgeActionAndWeight(i, j, ed, 'Ov');
                            //pred.successors.push(ed);
                            node.predecessors.push(ed);


                        }
                    }
                    else if (lev == this.LevelTop) { //top = deletions

                        if (j == 0) {

                        }
                        else {
                            var predLeft = this.getNodeTopOrder(i, j - 1, lev);
                            matches = this.findNodes('' + predLeft);
                            var pred = matches[0];
                            var ed = new DGEdge(pred, node);
                            this.setEdgeActionAndWeight(i, j, ed, '>');
                            node.predecessors.push(ed);

                            var predMiddle = this.getNodeTopOrder(i, j - 1, this.LevelMiddle);
                            matches = this.findNodes('' + predMiddle);
                            var pred = matches[0];
                            var ed = new DGEdge(pred, node);
                            this.setEdgeActionAndWeight(i, j, ed, 'O>');
                            //pred.successors.push(ed);
                            node.predecessors.push(ed);


                        }
                    }
                    else {

                        if ((i == 0) || (j == 0)) {
                        }
                        else {
                            //diag

                            var predDiag = this.getNodeTopOrder(i - 1, j - 1, lev);
                            matches = this.findNodes('' + predDiag);
                            pred = matches[0];
                            ed = new DGEdge(pred, node);
                            this.setEdgeActionAndWeight(i, j, ed, '\\');
                            node.predecessors.push(ed);
                        }

                        var predBottom = this.getNodeTopOrder(i, j, this.LevelBottom);
                        matches = this.findNodes('' + predBottom);
                        var pred = matches[0];
                        var ed = new DGEdge(pred, node);
                        this.setEdgeActionAndWeight(i, j, ed, 'C');
                        node.predecessors.push(ed);

                        var predTop = this.getNodeTopOrder(i, j, this.LevelTop);
                        matches = this.findNodes('' + predTop);
                        var pred = matches[0];
                        var ed = new DGEdge(pred, node);
                        this.setEdgeActionAndWeight(i, j, ed, 'C');
                        node.predecessors.push(ed);

                    }

                    this.nodes['' + topOrder] = node;
                }


            }
        }




    }

};

function DG3DAlignGraph(s,t,u,alignType,scoreMat,indelPen,mismatchPen,matchScore) {
    var args = [s,t,alignType,scoreMat,indelPen,mismatchPen,matchScore];
    DGAlignGraph.apply(this,args);

    this.u = u;

    this.numLevels = this.u.length + 1;


    var sinkNum  =  (this.rows) * (this.cols) * (this.numLevels) - 1;
    this.sinkLab = sinkNum.toString();

    this.parentSetEdgesAndWeights = this.setEdgeActionAndWeight;


    this.getNodeRowColLevel = function(node) {
        return [0,0,0];
    };

    this.getNodeTopOrder = function(r,c,lev) {
        return lev * this.cols * this.rows + (r * this.cols) + c;
    };


    this.setEdgeActionAndWeight = function(row,col,lev,ed,act) {

       // this.parentSetEdgesAndWeights(row,col,ed,act);
        ed.priorityScore = 5;

        switch (act) {
            case 'MMM':
                ed.edgeAction = this.s[col - 1] + this.t[row - 1] + this.u[lev - 1];
                if ((this.s[col - 1] == this.t[row - 1]) && (this.s[col - 1] == this.u[lev - 1])) {
                    ed.setWeight(this.matchScore);
                }
                else {
                    ed.setWeight(this.mismatchPen);
                }
                break;

            case '-MM':
                ed.edgeAction = '-' + this.t[row - 1] + this.u[lev - 1];
                ed.setWeight(this.mismatchPen);

                break;

            case 'M-M':
                ed.edgeAction = this.s[col - 1] + '-' + this.u[lev - 1];
                ed.setWeight(this.mismatchPen);

                break;
            case 'MM-':
                ed.edgeAction = this.s[col - 1] + this.t[row - 1] + '-';
                ed.setWeight(this.mismatchPen);

                break;

            case 'M--':
                ed.edgeAction = this.s[col - 1] + '-' + '-';
                ed.setWeight(this.indelPen);

                break;

            case '-M-':
                ed.edgeAction = '-' + this.t[row - 1] + '-';
                ed.setWeight(this.indelPen);

                break;

            case '--M':
                ed.edgeAction = '-' + '-' + this.u[lev - 1];
                ed.setWeight(this.indelPen);

                break;


            default:


                break;


        }

    };

    this.initGraph = function() {
        var topOrder;

        var progThreshold = 100;



        var lim =  (this.cols * this.rows) * this.numLevels;



        for (var lev = 0; lev < this.numLevels; ++lev) {
            for (var i = 0; i < this.rows; ++i) {
                for (var j = 0; j < this.cols; ++j) {


                    topOrder = (lev * this.rows * this.cols) + (i * this.cols) + j;

                    if (this.progressCallback) {
                        if ((topOrder % progThreshold == 0) || (topOrder == lim - 1)) {
                            this.progressCallback('Init graph', topOrder, lim - 1, 'align');
                        }
                    }


                    var node = new DGNode('' + topOrder);

                    if ((i == 0) || (j == 0) || (lev == 0)) {

                    }
                    else {
                        var predSTU = this.getNodeTopOrder(i - 1, j - 1, lev - 1);
                        matches = this.findNodes('' + predSTU);
                        var pred = matches[0];
                        var ed = new DGEdge(pred, node);
                        this.setEdgeActionAndWeight(i, j, lev, ed, 'MMM');
                        node.predecessors.push(ed);


                    }
                    if ((i == 0) || (j == 0)) {

                    }
                    else {
                        var predST = this.getNodeTopOrder(i - 1, j - 1, lev);
                        matches = this.findNodes('' + predST);
                        var pred = matches[0];
                        var ed = new DGEdge(pred, node);
                        this.setEdgeActionAndWeight(i, j, lev, ed, 'MM-');
                        node.predecessors.push(ed);

                    }

                    if ((i == 0) || (lev == 0)) {

                    }
                    else {
                        var predTU = this.getNodeTopOrder(i - 1, j, lev - 1);
                        matches = this.findNodes('' + predTU);
                        var pred = matches[0];
                        var ed = new DGEdge(pred, node);
                        this.setEdgeActionAndWeight(i, j, lev, ed, '-MM');
                        node.predecessors.push(ed);

                    }

                    if ((j == 0) || (lev == 0)) {

                    }
                    else {
                        var predSU = this.getNodeTopOrder(i, j - 1, lev - 1);
                        matches = this.findNodes('' + predSU);
                        var pred = matches[0];
                        var ed = new DGEdge(pred, node);
                        this.setEdgeActionAndWeight(i, j, lev, ed, 'M-M');
                        node.predecessors.push(ed);

                    }

                    if (i == 0) {

                    }
                    else {
                        var predT = this.getNodeTopOrder(i - 1, j, lev);
                        matches = this.findNodes('' + predT);
                        var pred = matches[0];
                        var ed = new DGEdge(pred, node);
                        this.setEdgeActionAndWeight(i, j, lev, ed, '-M-');
                        node.predecessors.push(ed);

                    }

                    if (j == 0) {

                    }
                    else {
                        var predS = this.getNodeTopOrder(i, j - 1, lev);
                        matches = this.findNodes('' + predS);
                        var pred = matches[0];
                        var ed = new DGEdge(pred, node);
                        this.setEdgeActionAndWeight(i, j, lev, ed, 'M--');
                        node.predecessors.push(ed);

                    }

                    if (lev == 0) {

                    }
                    else {
                        var predU = this.getNodeTopOrder(i, j, lev - 1);
                        matches = this.findNodes('' + predU);
                        var pred = matches[0];
                        var ed = new DGEdge(pred, node);
                        this.setEdgeActionAndWeight(i, j, lev, ed, '--M');
                        node.predecessors.push(ed);

                    }

                    this.nodes['' + topOrder] = node;
                }
            }
        }





    }
    
    this.alignStrings = function() {
        /*
         var sourceNum  = 0;
         var source = sourceNum.toString();
         var sinkNum  =  (t.length + 1) * (s.length + 1) - 1;
         var sink = sinkNum.toString();
         */

        // var pathData = this.longestPathBacktrack(sink,source);
        var pathData = this.longestPathBacktrack();
        var path = pathData[1];
        var longest = pathData[0];
        var edges = pathData[2];


        var startPosS = -1;
        var endPosS = -1;
        var startPosT = -1;
        var endPosT = -1;
        var startPosU = -1;
        var endPosU = -1;
        
        if (this.alignType == DGraph.alignTypeGlobal)  {
            //whole strings will be aligned
            startPosS = 0;
            endPosS = this.s.length -1;
            startPosT = 0;
            endPosT = this.t.length - 1;
            startPosU = 0;
            endPosU = this.u.length - 1;


        }
        else {
            var endCoord = this.getNodeRowColLevel(parseInt(path[0]));
            endPosS = endCoord[1] - 1; //col
            endPosT = endCoord[0] - 1; //row
            endPosU = endCoord[2] - 1; //level
            var startCoord = this.getNodeRowColLevel(parseInt(path[path.length -1]));
            startPosS = startCoord[1]; //col
            startPosT = startCoord[0]; //row
            startPosU = startCoord[2]; //level
        }

        var pathStr = '';
        for (var i = path.length - 1;i >=0;--i) {
            pathStr += path[i];
            if (i > 0)  {
                pathStr += '->';
            }

        }

        var sStr = '';
        var tStr = '';
        var uStr = '';
        var lcsStr = '';
        //edges.forEach(function(edge) {
        for (var i = edges.length - 1;i >=0;--i) {
            var edge = edges[i];
            sStr += edge.edgeAction.substring(0,1);
            tStr += edge.edgeAction.substring(1,2);
            uStr += edge.edgeAction.substring(2,3);
            if (( edge.edgeAction.substring(0,1) === edge.edgeAction.substring(1,2) && (edge.edgeAction.substring(0,1) === edge.edgeAction.substring(2,3)))) {
                lcsStr += edge.edgeAction.substring(0,1);
            }

        }



        var line = '';
        
    



        this.sAligned = sStr;
        this.tAligned = tStr;
        this.uAligned = uStr;

        return [longest,sStr,tStr,lcsStr,startPosS,endPosS,startPosT,endPosT,line,uStr,startPosU,endPosU];

    }

};


function DGAlignGraph(s,t,alignType,scoreMat,indelPen,mismatchPen,matchScore) {
    //DGGridGraph.apply(this,null,DGraph.fromGrid);
    var args = [t.length +1, s.length + 1,alignType,null,null,null];
    DGGridGraph.apply(this,args);
    this.s = s;
    this.t = t;
    
    this.scoreMat = scoreMat;

    this.sAligned = '';
    this.tAligned = '';
    
    this.indelPen = indelPen ? indelPen : 0;
    this.mismatchPen = mismatchPen ? mismatchPen : 0;
    this.matchScore = (matchScore == null) ? 1 : matchScore;
    
    
    this.rows = t.length + 1;
    this.cols = s.length + 1;

    /*
    this.initGraph = function() {

        var topOrder;

        for (var i = 0; i < this.rows; ++i) {
            for (var j = 0; j < this.cols; ++j) {
                topOrder = (i * this.cols) + j;
                if ((topOrder % 100) == 0) {
                    console.log('processed init: ' + topOrder);
                }
                var node = new DGNode('' + topOrder);
                var matches;
                if (i == 0) {

                }
                else {
                    var predUp = (i - 1) * this.cols + j;
                    matches = this.findNodes('' + predUp);
                    var pred = matches[0];
                    var ed = new DGEdge(pred, node);

                    ed.edgeAction = '-' + this.t[i-1];
                    ed.setWeight(this.indelPen * -1);
                    ed.priorityScore = 2;
                    pred.successors.push(ed);
                    node.predecessors.push(ed);
                }

                if (j == 0) {

                }
                else {
                    var predLeft = i * this.cols + j - 1;
                    matches = this.findNodes('' + predLeft);
                    pred = matches[0];
                    ed = new DGEdge(pred, node);
                    ed.edgeAction = this.s[j-1] + '-';
                    ed.setWeight(this.indelPen * -1);
                    ed.priorityScore = 1;
                    pred.successors.push(ed);
                    node.predecessors.push(ed);

                }

                if ((i == 0) || (j == 0)) {

                }
                else {
                    //diag


                        var predDiag = (i - 1) * this.cols + j - 1;
                        matches = this.findNodes('' + predDiag);
                        pred = matches[0];
                        ed = new DGEdge(pred, node);
                        ed.edgeAction = this.s[j-1] + this.t[i-1];
                        ed.setWeight(this.s[j-1] == this.t[i-1] ? this.matchScore : this.mismatchPen * -1);
                        if (this.s[j-1] === this.t[i-1]) {
                            ed.priorityScore = 4; //match
                        }
                        else {
                            ed.priorityScore = 3; //mismatch
                        }

                        pred.successors.push(ed);
                        node.predecessors.push(ed);
                    }



                //this.nodes.push(node);
                this.nodes['' + topOrder] = node;

            }
        }
    }
    */



    this.scoreAlignment = function(a,b) {
        
        var totScore = 0;
        
        for (var i = 0;i < a.length;++i) {
            if ((a[i] == '-') || (b[i] == '-')) {
               totScore += this.indelPen * -1;
            }
            else {
                if (this.scoreMat) {
                    totScore += this.scoreMat[a[i]][b[i]];
                }
                else {
                    totScore += a[i] == b[i] ? this.matchScore : this.mismatchPen * -1;
                }
            }
       
        }
            
            return totScore;
    };




this.setEdgeActionAndWeight = function(row,col,ed,act) {


        switch (act) {
            case 'v':
                ed.edgeAction = '-' + this.t[row - 1];
                ed.setWeight(this.indelPen * -1);
                ed.priorityScore = 2;
                break;
            case '>':
                ed.edgeAction = this.s[col - 1] + '-';
                /*
                if ((row == 0) && ((this.alignType == DGraph.alignTypeFitting) || (this.alignType == DGraph.alignTypeOverlap))) {
                    ed.setWeight(0); // free tax ride for deletes of s
                }
                else {
                */
                    ed.setWeight(this.indelPen * -1);
                //}
                ed.priorityScore = 1;
                break;
            case '\\':
                ed.edgeAction = this.s[col - 1] + this.t[row - 1];
                if (this.scoreMat) {
                    ed.setWeight(this.scoreMat[this.t[row - 1]][this.s[col - 1]]);
                }
                else {
                    ed.setWeight(this.s[col - 1] == this.t[row - 1] ? this.matchScore : this.mismatchPen * -1);
                }
                if (this.s[col - 1] === this.t[row - 1]) {
                    ed.priorityScore = 4; //match
                }
                else {
                    ed.priorityScore = 3; //mismatch
                }

                break;
            default:

                ed.edgeAction = '' + ''; //close gap
                ed.setWeight(0);
                ed.priorityScore = 5;

                break;

        }

    };


    
    this.alignStrings = function() {

        /*
        var sourceNum  = 0;
        var source = sourceNum.toString();
        var sinkNum  =  (t.length + 1) * (s.length + 1) - 1;
        var sink = sinkNum.toString();
        */
        
       // var pathData = this.longestPathBacktrack(sink,source);
        var pathData = this.longestPathBacktrack();
        var path = pathData[1];
        var longest = pathData[0];
        var edges = pathData[2];


        var startPosS = -1;
        var endPosS = -1;
        var startPosT = -1;
        var endPosT = -1;
        if (this.alignType == DGraph.alignTypeGlobal)  {
            //whole strings will be aligned
            startPosS = 0;
            endPosS = this.s.length -1;
            startPosT = 0;
            endPosT = this.t.length - 1;


        }
        else {
            var endCoord = this.getNodeRowCol(parseInt(path[0]));
            endPosS = endCoord[1] - 1; //col
            endPosT = endCoord[0] - 1; //row
            var startCoord = this.getNodeRowCol(parseInt(path[path.length -1]));
            startPosS = startCoord[1]; //col
            startPosT = startCoord[0]; //row
        }

        var pathStr = '';
        for (var i = path.length - 1;i >=0;--i) {
            pathStr += path[i];
            if (i > 0)  {
                pathStr += '->';
            }

        }

        var sStr = '';
        var tStr = '';
        var lcsStr = '';
        //edges.forEach(function(edge) {
        for (var i = edges.length - 1;i >=0;--i) {
            var edge = edges[i];
            sStr += edge.edgeAction.substring(0,1);
            tStr += edge.edgeAction.substring(1,2);
            if ( edge.edgeAction.substring(0,1) === edge.edgeAction.substring(1,2)) {
                lcsStr += edge.edgeAction.substring(0,1);
            }

        }



        var sp = '&nbsp';
        var line = '';
        var sExt = '-' + this.s;
        var tExt = '-' + this.t;

        for (var i = 0;i < this.rows + 1;++i) {

            for (var j = 0;j < this.cols + 1; ++j) {
                if ((i == 0) && (j == 0)) {
                    line += sp + sp + sp + sp;
                }
                else if (i == 0) {
                    line += sExt.substring(j - 1, j) + sp + sp + sp;
                }
                else if (j == 0) {
                    line += tExt.substring(i - 1, i) + sp + sp;
                }
                else {
                    var topOrder = (i-1) * this.cols + j -1;
                    var topLab = '' + topOrder;
                    var node =  this.findNodes('' + topOrder)[0];
                    var extraSp =  (node.longestPathToThisNode >= 0) ? '&nbsp;' : '';


                    if (path.indexOf(topLab) > -1) {
                        line += '<b>' + extraSp  +  node.longestPathToThisNode + '</b>' + sp + sp;
                    }
                    else {
                        line += extraSp +  node.longestPathToThisNode + sp + sp;
                    }


                }

            }

            line += '<BR>';


        }




        this.sAligned = sStr;
        this.tAligned = tStr;
        
        return [longest,sStr,tStr,lcsStr,startPosS,endPosS,startPosT,endPosT,line];


    }





}

function DGGridGraph(rows,cols,alignType,downWeights,rightWeights,diagWeights) {

   // DGGraph.apply(this,null,);
    var args = [null,null,alignType];
    DGGraph.apply(this,args);
    
    this.rows = rows;
    this.cols = cols;
    this.downWeights = downWeights;
    this.rightWeights = rightWeights;
    this.diagWeights = diagWeights? diagWeights : null;

    this.alignType = alignType;

    var sourceNum  = 0;
    this.sourceLab = sourceNum.toString();
    var sinkNum  =  (this.rows) * (this.cols) - 1;
    this.sinkLab = sinkNum.toString();


    
    this.initGraph = function() {
        
        var topOrder;

        var progThreshold = 100;

        var lim = this.cols * this.rows;

        for (var i = 0;i < this.rows; ++ i) {
            for (var j = 0;j < this.cols; ++j) {
                topOrder = (i * this.cols) + j;

                if (this.progressCallback) {
                        if ((topOrder % progThreshold == 0) || (topOrder == lim - 1)) {
                            this.progressCallback('Init graph',topOrder,lim - 1,'align');
                        }
                 }


                    var node = new DGNode('' + topOrder);
                var matches;
                if (i == 0) {

                }
                else {
                    var predUp = (i-1) * this.cols + j;
                    matches = this.findNodes('' + predUp);
                    var pred = matches[0];
                    var ed = new DGEdge(pred,node);
                   // ed.setWeight(this.downWeights[i-1][j]);
                    this.setEdgeActionAndWeight(i,j,ed,'v');
                    //pred.successors.push(ed);
                    node.predecessors.push(ed);
                }

                if (j == 0) {

                }
                else {
                    var predLeft = i * this.cols + j-1;
                    matches = this.findNodes('' + predLeft);
                    pred = matches[0];
                    ed = new DGEdge(pred,node);
                   // ed.setWeight(this.rightWeights[i][j-1]);
                    this.setEdgeActionAndWeight(i,j,ed,'>');
                    //pred.successors.push(ed);
                    node.predecessors.push(ed);

                }

                if ((i == 0) || (j == 0)) {

                }
                else {
                    //diag


                  //  if ((this instanceof DGAlignGraph)  || (this.diagWeights)) {}


                   // if (this.diagWeights) {
                        var predDiag = (i - 1) * this.cols + j - 1;
                        matches = this.findNodes('' + predDiag);
                        pred = matches[0];
                        ed = new DGEdge(pred, node);
                        //ed.setWeight(this.diagWeights[i - 1][j - 1]);
                        this.setEdgeActionAndWeight(i,j,ed,'\\');

                        //pred.successors.push(ed);
                        node.predecessors.push(ed);
                 //   }
                }



                //this.nodes.push(node);
                this.nodes['' + topOrder] = node;
                
            }
        }


        
        
    };

    this.setEdgeActionAndWeight = function(row,col,ed,act) {

        ed.edgeAction = act;

        switch (act) {
            case 'v':
                ed.setWeight(this.downWeights[row-1][col]);
                break;
            case '>':
                ed.setWeight(this.rightWeights[row][col-1]);
                break;
            case '\\':
                if (this.diagWeights) {
                    ed.setWeight(this.diagWeights[row - 1][col - 1]);
                }
                else {
                    ed.setWeight(0);
                }
                break;
            default:
                break;

        }



    };

    
    
    this.getNodeRow = function(node) {
        if (node) {
            var topOrder = node.getTopOrder();
            return Math.floor(topOrder / this.cols);
        } else {
            return -1;
        }
    };
    
    this.getNodeCol = function(node) {
        if (node) {
            var topOrder = node.getTopOrder();
            return topOrder % this.cols;
        }
        else {
            return - 1;
        }
    };
    
  

    this.getNodeRowCol = function(topOrder) {
        var row = Math.floor(topOrder / this.cols);
        var col = topOrder % this.cols;
        return [row,col];
    };

    this.getTopOrder = function(row,col) {
        return (row * this.cols) + col;
    };
    
    this.getAbsTopOrderNode = function(node,startRow,startCol,absColSize) {
        var r = this.getNodeRow(node);
        var c = this.getNodeCol(node);
        return ((r + startRow) * absColSize) + c + startCol;
    };
    
    this.rowScores = function(r) {
        var scores = [];
        for (var i = r * this.cols;i < (r+1) * this.cols;++i) {
            scores.push(this.nodes[i].longestPathToThisNode);
        }
        return scores;
    }


}

function DGGraph(source,sourceType,alignType) {
    this.sourceType = sourceType;


    this.sourceLab = '';
    this.sinkLab = '';

    this.nodes = {};
    
    this.bestNode = null;
    this.bestNodeLastRow = null;
    this.bestNodeLastCol = null;
    
  
    this.progressCallback = null;

    this.alignType = alignType;


    this.overallBestScore = -99999; // keep track of best score for use in local align backtrack
    this.overallBestScoreLastRow = -99999; // keep track of best score  in last row for use in fitting align backtrack
    this.overallBestScoreLastCol = -99999; // keep track of best score  in last dol for use in overlap align backtrack


    switch (this.sourceType) {
        case DGraph.fromAdjList:
            this.adjList = source; // array of adjacencies
            break;
        default:
            break;
    }

    this.initGraph = function() {
        switch (this.sourceType) {
            case DGraph.fromAdjList:
                this.buildNodesFromAdjList(this.adjList);
                break;
            default:
                break;
        }
    };

    this.buildNodesFromAdjList = function(adjList) {
        var svdThis = this;

        adjList.forEach(function (adj) {

            var spl = adj.split('->');
            var pref = spl[0];
            var spl2 = spl[1].split(':');
            var suf = spl2[0];
            var weight = spl2[1];

            //var weight = spl2[1];

            var matches = svdThis.findNodes(pref);
            if (matches.length > 0) {
            }
            else {
                var node = new DGNode(pref);
                //svdThis.nodes.push(node);
                svdThis.nodes[pref] = node;
            }
        });

        adjList.forEach(function (adj,i) {

            var spl = adj.split('->');
            var pref = spl[0];
            var spl2 = spl[1].split(':');
            var suf = spl2[0];
            var weight = spl2[1];
 
            var sufNode;
            
            var prefNode = svdThis.findNodes(pref)[0];
    
            var matches = svdThis.findNodes(suf);
            if (matches.length == 0) {
                sufNode = new DGNode(suf);
                //svdThis.nodes.push(sufNode);
                svdThis.nodes[suf] = sufNode;
                matches = svdThis.findNodes(suf);
                
            }
            sufNode = svdThis.findNodes(suf)[0];
  

            var ed = new DGEdge(prefNode, sufNode);
            ed.setWeight(parseInt(weight));
            prefNode.successors.push(ed);
            sufNode.predecessors.push(ed);
            

        });



    };



    this.findNodes = function(label) {
        //return all nodes which match

        var matchingNodes = [];

        /*
        for (var i = 0;i < this.nodes.length;++i ) {
    
            if (this.nodes[i].label === label) {
                    matchingNodes.push(this.nodes[i]);
            }
        }
        */
        if (label in this.nodes) {
            matchingNodes.push(this.nodes[label]);
        }

        return matchingNodes;
    };
    
    this.longestPathsDynamic = function(sourceLab,sinkLab,partStart,partEnd) {
                
        //if source supplied, use that, otherwise use entire graph source
        if (sourceLab) {
            
        }
        else {
            sourceLab = this.sourceLab;
        }
        if (sinkLab) {
            
        }
        else {
            sinkLab = this.sinkLab;
        }
        //create part at a time - used in space efficient alignment
        partStart = partStart ?  partStart : 0;
        partEnd = partEnd ?  partEnd : Object.keys(this.nodes).length + 1;



        var progThreshold = 100;

        
        for (var i = partStart;i < partEnd;++i) {

            if (this.progressCallback) {

                /*
                    if ((i % progThreshold == 0) || (i == partEnd - 1)) {
                        this.progressCallback('Building dynamic paths', i, partEnd - 1, 'align');
                    }
                    */


            }

            //if (i % 100 == 0) {
            //    console.log('dynamic: ' + i);
            //}
            var matches = this.findNodes('' + i);
            if (matches.length == 0) {

            }
            else {

                var node = matches[0];

                var longest = -99999;
                var longestEdge;
                var longestPriority = -1;

                if (node.predecessors.length == 0) {
                    longest = 0;
                    if (node.label === sourceLab) {
                        
                    }
                    else {
                        node.onSourcePath = false;
                    }
                    
                }
                else {

                    var foundSourcePath = false;
                    node.predecessors.forEach(function (predEdge) {

                        if (predEdge.sourceNode.onSourcePath) {
                            var prevNodeLongest = predEdge.sourceNode.longestPathToThisNode;
                            var now = prevNodeLongest + predEdge.edgeWeight;
                            if (now > longest) {
                                longest = now;
                                longestEdge = predEdge;
                                longestPriority = predEdge.priorityScore;
                            }

                           else if (now == longest) {
                                if (predEdge.priorityScore > longestPriority ) {
                                    longest = now;
                                    longestEdge = predEdge;
                                    longestPriority = predEdge.priorityScore;
                                }
                                /*
                                if (predEdge.edgeAction.indexOf(' ') == -1) { //preference for diag, even if eq
                                    longest = now;
                                    longestEdge = predEdge;
                                }
                                */

                            }

                            foundSourcePath = true;
                        }
                    });
                    if ((foundSourcePath) || (node.label === sourceLab)) {
                            if (node.label === sourceLab) {
                                longest = 0;
                            }
                    }
                    else {
                        node.onSourcePath = false;
                    }
                }
                if ((this.alignType == DGraph.alignTypeLocal) && (longest < 0)) {
                    longest = 0;
                      node.setLongestPathToThisNode(longest, null);
                }
                else if (((this.alignType == DGraph.alignTypeFitting) || (this.alignType == DGraph.alignTypeOverlap)) && (longest < 0)) {

                    var coord = this.getNodeRowCol(parseInt(node.label));
                    if (coord[0] == 0) { //1st row free taxi ride
                        longest = 0;
                        node.setLongestPathToThisNode(longest, null);
                    }
                    else {
                        node.setLongestPathToThisNode(longest, longestEdge);
                    }
                }
                else {
                    node.setLongestPathToThisNode(longest, longestEdge);
                }
                if (longest >= this.overallBestScore) {
                    this.overallBestScore = longest;
                    this.bestNode = node;
                }
                if ((this.alignType == DGraph.alignTypeFitting) &&  (i > ((this.rows - 1) * this.cols) )) {
                    //last row
                    if (longest >= this.overallBestScoreLastRow) {
                        this.bestNodeLastRow = node;
                        this.overallBestScoreLastRow = longest;
                    }
                }
                if ((this.alignType == DGraph.alignTypeOverlap) &&  ((i + 1) % this.cols == 0) ) {  //last col

                    if (i == this.cols - 1) { //ie 1st row - don't count. Need at least something to align
                    }
                    else {
                        if (longest >= this.overallBestScoreLastCol) {
                            this.bestNodeLastCol = node;
                            this.overallBestScoreLastCol = longest;
                        }
                    }
                }
            }
        }
        
    }
    
    this.longestPathBacktrack = function(sinkLab,sourceLab) {
        
        sinkLab  = sinkLab ? sinkLab : this.sinkLab;
        sourceLab = sourceLab ? sourceLab : this.sourceLab;
        
        var matches = this.findNodes(sinkLab);
        var sink = null;
        if (matches.length > 0) {
            sink = matches[0];
        }

        var source = null;
        var sourceMatches = this.findNodes(sourceLab);
        if (sourceMatches.length > 0) {
            source = sourceMatches[0];
        }


        if (matches.length > 0) {
            sink = matches[0];
        }

        //Used for fitting/local/overlap:
        var startPosS = -1;
        var startPosT = -1;
        var endPosS = -1;
        var endPosT = -1;

        var startPosU = -1;
        var endPosU = -1;
        

        var curNode;
        if (this.alignType == DGraph.alignTypeLocal) {
            curNode = this.bestNode;
            if (this.u.length > 0) {
                var coord = this.getNodeRowColLevel(parseInt(curNode.label));
                endPosT = coord[0] - 1; //row
                endPosS = coord[1] - 1; //col
                endPosU = coord[2] - 1; //level
            }
            else {
                var coord = this.getNodeRowCol(parseInt(curNode.label));
                endPosT = coord[0] - 1; //row
                endPosS = coord[1] - 1; //col
            }
        }
        else  if (this.alignType == DGraph.alignTypeFitting) {
            curNode = this.bestNodeLastRow;
            var coord = this.getNodeRowCol(parseInt(curNode.label));
            endPosT = coord[0] - 1; //row
            endPosS = coord[1] - 1; //col
        }
        else if (this.alignType == DGraph.alignTypeOverlap) {
            curNode = this.bestNodeLastCol;
            var coord = this.getNodeRowCol(parseInt(curNode.label));
            endPosT = coord[0] - 1; //row
            endPosS = coord[1] -1; //col
        }
        else {
            curNode = matches[0]; //sink;
       }
        var long = curNode.longestPathToThisNode;
        
        var path = []; //nodes in path
        var edgePath = []; //edges in path

        var progThreshold = 1;
        
        path.push(curNode.label);
        while (curNode.label !== source.label) {

            if (this.progressCallback) {
                if (path.length % progThreshold == 0)  {
                    this.progressCallback('Backtrack',path.length,0,'align');
                }
            }

            edgePath.push(curNode.longestEdgeToThisNode);
            curNode = curNode.longestEdgeToThisNode.sourceNode;
            path.push(curNode.label);
            if ((this.alignType == DGraph.alignTypeLocal) && (curNode.longestEdgeToThisNode == null)) {
                break;
            }
            else if ((this.alignType == DGraph.alignTypeFitting) && (curNode.longestPathToThisNode == 0) && curNode.predecessors.length == 1) {
                break;
            }
            else if ((this.alignType == DGraph.alignTypeOverlap) && (curNode.longestPathToThisNode == 0) && curNode.predecessors.length == 1) {
                break;
            }
        }

        if  ((this.alignType == DGraph.alignTypeLocal)
         || (this.alignType == DGraph.alignTypeFitting)
          || (this.alignType == DGraph.alignTypeOverlap) ){
            
            if (this.u.length > 0 ) {
                var coord = this.getNodeRowColLevel(parseInt(curNode.label));
                startPosS = coord[1]; //col
                startPosT = coord[0]; //row
                startPosU = coord[2]; //level
            }
            else {
                var coord = this.getNodeRowCol(parseInt(curNode.label));
                startPosS = coord[1]; //col
                startPosT = coord[0]; //row
            }
        }
        
        return [long,path,edgePath];
    }


}

//
function DGraph(source,sourceType,graphType,k,makeCycle,pairDist,variableOverlap) {

    this.sourceType = sourceType;

    this.graphType = graphType;

    switch (this.sourceType)  {
        case DGraph.fromDna:
        case DGraph.fromPairedDna:
            this.dna = source; //dna
            break;
        case DGraph.fromReads:
        case DGraph.fromPairedReads:
            this.reads = source; //array of reads
            break;
        case DGraph.fromAdjList:
            this.adjList = source; // array of adjacencies
            break;
        default:
            break;
    }


    if (makeCycle) {
        this.makeCycle = makeCycle;
    }
    else {
        this.makeCycle = false;
    }
   // this.dna = dna;

    this.k = k;
    this.pairDist = pairDist ? pairDist : null;

    this.variableOverlap = variableOverlap ? true : false; //if true, constructs overlap graph with all overlaps. Otherwise just uses k-1 overlaps

    /*
    if (this.pairDist) {
        this.dna = this.dna + this.dna.substring(0,this.k + this.pairDist - 1); // extend so pairs work up to end
    }*/

    this.nodes = [];
    this.nodesDict = {};

    this.edgePath = [];

    this.numEdges = function() {
        return this.nodes.map(function(el) {
            return el.outDegree();
        }).reduce(function(a,b) {
            return a + b;
        });
    };

    this.numNodes = function() {
        return this.nodes.length;
    };

    this.isBalanced = function() {
        var balanced = true;

        for (var i = 0;i < this.nodes.length; ++i) {
            if (this.nodes[i].InDegree() == this.nodes[i].outDegree()) {
            }
            else {
                balanced = false;
                break;

            }
        }

        return balanced;



    };

    this.isBalancedLinear = function() {
        //checks if all nodes balanced except one start node (predecessors = successors - 1) and one
        //end node (predecessors = successors + 1)

        var balanced = true;

        var foundStart = false;
        var foundEnd = false;

        for (var i = 0;i < this.nodes.length; ++i) {
            if (this.nodes[i].InDegree() == this.nodes[i].outDegree()) {
            }
            else if (this.nodes[i].InDegree() == this.nodes[i].outDegree() - 1) {
                if (foundStart) {
                    balanced = false; //already found a start node - can't have two
                    break;
                }
                else {
                    foundStart = true;
                }
            }
             else if (this.nodes[i].InDegree() == this.nodes[i].outDegree() + 1) {
                if (foundEnd) {
                    balanced = false; //already found an end node - can't have two
                    break;
                }
                else {
                    foundEnd = true;
                }
            }
            else {
                //completely unbalanced
                balanced = false;
                break;

            }
        }

        return balanced;



    };


    this.getAdjList = function() {

       // var resStr = '';

        var adj = [];


        var adjStr;

        this.nodes.forEach(function(node) {

            adjStr = node.dna + ' -> ';
            node.successors.forEach(function (suc, i) {
                adjStr += suc.targetNode.dna;
                if (i < node.successors.length - 1) {
                    adjStr += ',';
                }

            });
            adj.push(adjStr);
        });
        return adj;
    };

    this.sumInfo = function() {

        var n = this.numNodes();
        var e = this.numEdges();
        var bc =  this.isBalanced() ? 'Cycle Balanced' : 'Cycle Unbalanced';
        var bl =  this.isBalancedLinear() ? 'Linear Balanced' : 'Linear Unbalanced';
        var f =  this.definiteFirst() ? this.definiteFirst().dna : '?';
        var l = this.definiteLast() ? this.definiteLast().dna : '?';

        f = squishString(f,30);
        l = squishString(l,30);
        
        return ('Nodes / Edges: ' + n + ' / ' + e + '\n'
              //  + 'Num Edges: ' + e + '\n'
                + 'First / Last: ' + f   + ' / ' + l + '\n'
               // + 'Last node: ' + l  + '\n'
                + bl + '\n'
                + bc);

    };

    this.resetEdges = function() {

        this.nodes.forEach(function(node) {
            node.resetEdgesVisited();
        });

        this.edgePath = [];



    };



    this.definiteFirst = function() {

        var firstFound = false;
        var first = null;

        var i;

        if (this.graphType == DGraph.hamGraph) {
            for (i = 0; i < this.nodes.length; ++i) {
                if (this.nodes[i].InDegree() == 0) {
                    first = this.nodes[i];
                    firstFound = true;
                    break;
                }
            }
        }
        else {
            for (i = 0; i < this.nodes.length; ++i) {
                if (this.nodes[i].InDegree() == this.nodes[i].outDegree() - 1) {
                    first = this.nodes[i];
                    firstFound = true;
                    break;
                }
            }
        }

        return first;


    };

    this.definiteLast = function() {
        var lastFound = false;
        var last = null;

        var i;

        if (this.graphType == DGraph.hamGraph) {
            for (i = 0; i < this.nodes.length; ++i) {
                if (this.nodes[i].outDegree() == 0) {
                    last = this.nodes[i];
                    lastFound = true;
                    break;
                }
            }

        }
        else {
            for (i = 0; i < this.nodes.length; ++i) {
                if (this.nodes[i].outDegree() == this.nodes[i].InDegree() - 1) {
                    last = this.nodes[i];
                    lastFound = true;
                    break;
                }
            }
        }

        return last;


    };


    this.maximalNonBranchingPaths = function() {


        var svdThis = this;

        var contigs = [];

        var isolatedCycleUsedNodes = [];
        var edgePath = [];

        this.nodes.forEach(function(node) {
            if (  (node.InDegree() == 1) && (node.outDegree() == 1)) {
                edgePath = [];
                if (isolatedCycleUsedNodes.indexOf(node)==-1) {
                    //check for isolated cycle
                    var startNode = node;
                    var cyc = [];
                    cyc.push(startNode);
                    var curNode = startNode.successors[0].targetNode;
                    edgePath.push(new DEdge(startNode,curNode));


                    var isolatedCycleFound = false;
                    while ((curNode.InDegree() == 1) && (curNode.outDegree() == 1)) {
                        if (curNode == startNode) {
                            cyc.push(curNode);
                            isolatedCycleFound = true;
                            break;
                        }
                        else {
                            cyc.push(curNode);
                            edgePath.push(new DEdge(curNode,curNode.successors[0].targetNode))
                            curNode = curNode.successors[0].targetNode;

                        }
                    }
                    if (isolatedCycleFound) {
                        var cycPath = [];
                        cyc.forEach(function (el) {
                            cycPath.push(el.dna);
                            isolatedCycleUsedNodes.push(el);
                        });
                       // contigs.push(cycPath);
                        contigs.push(edgePath);
                    }
                }

            }
            else {

                node.successors.forEach(function(suf) {
                    var curNode;
                    var contig = [];
                    edgePath = [];
                    edgePath.push(new DEdge(suf.sourceNode,suf.targetNode));
                    contig.push(suf.sourceNode.dna);
                    contig.push(suf.targetNode.dna);
                    curNode = suf.targetNode;
                    while ((curNode.outDegree() == 1) && (curNode.InDegree() == 1)) {
                        contig.push(curNode.successors[0].targetNode.dna);
                        edgePath.push(new DEdge(curNode,curNode.successors[0].targetNode));
                        curNode = curNode.successors[0].targetNode;
                    }
                    //contigs.push(contig);
                    contigs.push(edgePath);

                });
            }

        });


        contigs.sort(function(a, b){
            // ASC  -> a.length - b.length
            // DESC -> b.length - a.length
            return b.length - a.length;
        });

        return contigs;



    };


    this.edgePathToText = function() {
        var str = '';
        for (var i=0;i < this.edgePath.length; ++i) {
            if (i == 0) {
                str += this.edgePath[i].sourceNode.dna + '->' + this.edgePath[i].targetNode.dna;
            }
            else {
                str += '->' + this.edgePath[i].targetNode.dna;
            }

        }
        return str;


    };

    this.edgePathReconstructedPairs = function(edgeP) {

        //if edgeP supplied, reconstuct that, else reconstruct graph's own edgepath

        var currEdgePath;
        if (edgeP) {
            currEdgePath = edgeP;
        }
        else {
            currEdgePath = this.edgePath;
        }

        //var svdThis = this;

        var origStr = '';
        var pairStr = '';

        currEdgePath.forEach(function(ed,i) {
            var lab = ed.edgeLabel();



            if (i == 0) {
                origStr += lab[0];
                pairStr += lab[1];
            }
            else {
                origStr += lab[0].substring(lab[0].length - 1);
                pairStr += lab[1].substring(lab[1].length - 1);
            }



        });

       // var totLen = origStr.length + this.k + this.pairDist;

        var extra =  this.k + this.pairDist;

        var padding = '';
        if (pairStr.length < extra) {
            for (var i = pairStr.length;i < extra;++i) {
                padding += 'x';
            }

        }

       // var overlapLen = origStr.length - extra;
        var mismatch = false;
        if (origStr.substring(extra) == pairStr.substring(0,origStr.length - extra)) {

        }
        else {
            mismatch = true;
        }


        return origStr + padding + pairStr.substring(pairStr.length - extra) + (mismatch ? ' Mismatch' : '');



        /*
        var len = this.edgePath.length + this.k - 1  + this.k + this.pairDist;
        var reconAr = [];
        for (var i = 0;i < len; ++i) {
            reconAr[i] = 'x';
        }

        var svdThis = this;

        this.edgePath.forEach(function(ed,i) {
            var lab = ed.edgeLabel();
            for (var j = 0; j < lab[0].length;++j) {
                reconAr[i + j] = lab[0].substring(j,j+1);
                reconAr[i + j + lab[0].length + svdThis.pairDist] = lab[1].substring(j,j+1);
            }
        });

        var str = '';
        reconAr.forEach(function(el) {
            str+=el;

        });
        return str;
        */


    };

    this.edgePathReconstructed = function(edgeP) {
        //if edgeP supplied, reconstuct that, else reconstruct graph's own edgepath

        var currEdgePath;
        if (edgeP) {
            currEdgePath = edgeP;
        }
        else {
            currEdgePath = this.edgePath;
        }

        if ((this.sourceType == DGraph.fromPairedReads) || (this.sourceType == DGraph.fromPairedDna)) {
            return this.edgePathReconstructedPairs(edgeP);
        }

        var str = '';
        for (var i=0;i < currEdgePath.length; ++i) {
            if (i == 0) {
                str += currEdgePath[i].edgeLabel();
            }
           // else if ((this.makeCycle) && (i >= this.edgePath.length - this.k + 1)) {
            else if ((this.isBalanced()) && (i >= currEdgePath.length - this.k + 1)) {


            }
            else   {
                //str += this.edgePath[i].edgeLabel().substring(this.edgePath[i].edgeLabel().length - 1);
                str+= currEdgePath[i].targetNode.dna.substring(currEdgePath[i].targetNode.dna.length - 1);
            }


        }

        return str;


    };

    this.nodesVisited = function() {
        return this.nodes.filter(function(node) {
            return node.visited;
        });
    };

    this.visitedNodesWithUnvisitedEdges = function() {

        var vis = this.nodesVisited();

        return vis.filter(function(node) {
                var r = node.getRandomNonVisitedEdge();
                return r != -1;
        });


    };


    this.makeCyclical = function() {
       // var f = this.definiteFirst();
       // var l = this.definiteLast();

       // if (f || l) {

        /*
        if (this.pairedDist) {
            return; //already added to dna to make enough
        }
        */

        var kmers;
        if (this.pairDist)  {
            kmers = kmerPairedComposition(this.dna.substring(this.dna.length - (this.k - 1)) + this.dna.substring(0,this.k-1),this.k,this.pairDist);
        }
        else {
            kmers = kmerComposition(this.dna.substring(this.dna.length - (this.k - 1)) + this.dna.substring(0,this.k-1),this.k);
        }
        var svdThis = this;
        kmers.forEach(function(kmer) {
            var pref,suf;
            if (pairDist) {
                pref = [kmer[0].substring(0,svdThis.k - 1),kmer[1].substring(0,svdThis.k - 1)];
                suf = [kmer[0].substring(1,svdThis.k),kmer[1].substring(1,svdThis.k)];
            }
            else {
                pref = kmer.substring(0,svdThis.k - 1);
                suf = kmer.substring(1,svdThis.k);

            }

            var prefNode = svdThis.findNode(pref);
            var sufNode = svdThis.findNode(suf);
            if (!prefNode) {
                prefNode = new DNode(pref);
                svdThis.nodes.push(prefNode);
            }

            if (!sufNode) {
                sufNode = new DNode(suf);
                svdThis.nodes.push(sufNode);
            }

            prefNode.successors.push(new DEdge(prefNode,sufNode));
            sufNode.predecessors.push(prefNode);




            });

     //   }

    };


    this.findOverlapNodesWithPrefix = function(node)   {
        var matches = [];

        for (var i = 0;i < this.nodes.length;++i ) {
            if (this.nodes[i] === node) {

            }
            else {
                if ((this.sourceType == DGraph.fromPairedDna) || (this.sourceType == DGraph.fromPairedReads)) {
                    if ((this.nodes[i].suffix()[0] == node.prefix()[0]) && (this.nodes[i].suffix()[1] == node.prefix()[1])) {
                        matches.push(this.nodes[i]);
                    }
                }
                else {
                    if (this.nodes[i].suffix() == node.prefix()) {
                        matches.push(this.nodes[i]);
                    }
                }

            }
        }
        return matches;


    };

    this.findOverlapNodesWithSuffix = function(node)   {
        var matches = [];

        for (var i = 0;i < this.nodes.length;++i ) {
            if (this.nodes[i] === node) {

            }
            else {
                if ((this.sourceType == DGraph.fromPairedDna) || (this.sourceType == DGraph.fromPairedReads)) {
                    if ((this.nodes[i].prefix()[0] == node.suffix()[0]) && (this.nodes[i].prefix()[1] == node.suffix()[1])) {
                        matches.push(this.nodes[i]);
                    }

                }
                else {
                    if (this.nodes[i].prefix() == node.suffix()) {
                        matches.push(this.nodes[i]);
                    }
                }
            }


        }
        return matches;


    };

    this.findNodes = function(dna) {
        //return all nodes which match

        var matchingNodes = [];

        for (var i = 0;i < this.nodes.length;++i ) {
            if (typeof dna == 'object') {
                if (this.nodes[i].dna === dna[0] && this.nodes[i].pairedDna === dna[1]) {
                    matchingNodes.push(this.nodes[i]);
                }
            }
            else {
                if (this.nodes[i].dna === dna) {
                    matchingNodes.push(this.nodes[i]);


                }
            }


        }
        return matchingNodes;

    };

    this.findNode = function(dna,repNum) {

        if ((this.graphType == DGraph.hamGraph) || repNum) {
            //still uses array of nodes, not dict of nodes


            for (var i = 0; i < this.nodes.length; ++i) {
                if (typeof dna == 'object') {
                    if (this.nodes[i].dna === dna[0] && this.nodes[i].pairedDna === dna[1]) {
                        if (repNum) {
                            if (repNum == this.nodes[i].repeatNum) {
                                return this.nodes[i];
                            }
                        }
                        else {
                            return this.nodes[i];
                        }
                    }


                }
                else {
                    if (this.nodes[i].dna === dna) {
                        if (repNum) {
                            if (repNum == this.nodes[i].repeatNum) {
                                return this.nodes[i];
                            }
                        }
                        else {
                            return this.nodes[i];
                        }
                    }
                }


            }
            return null;
        }
        else {
            //new: use dict of nodes for speed
            var dictKey;
            if (typeof dna == 'object') {
                dictKey = dna[0] + dna[1];
            }
            else {
                dictKey = dna;
            }

            if (dictKey in this.nodesDict) {
                return this.nodesDict[dictKey];

            }
            else {
                return null;
            }
        }





    };

    this.chooseStartCycle = function(linear) {

        //linear flag: used if looking for linear path rather than cycle

        //var nodesAvail = [];

      //  if (this.nodesVisited().length == 0) {
        if (this.edgePath.length == 0) {
             //starting - get first linear node

            var r;

            if (linear) {
                var f = this.definiteFirst();
                if (f) {
                    return [f, null];
                }
                else {
                    return [null, null]; // no definite first. Won't work
                }
            }
            else {
                r = getRandomInt(0,this.nodes.length - 1);
                return [this.nodes[r],null];
                return [this.nodes[r],null];
            }


        }
        else {
            // choose an unexplored path along nodes already visited
            var currNode = this.edgePath[0].sourceNode;
            for (var i = 0;i < this.edgePath.length;++i) {
                r = currNode.getRandomNonVisitedEdge();
                if (r == -1) {
                    currNode = this.edgePath[i].targetNode;
                }
                else {
                    return [currNode,i];
                }

            }
            return [null,null];
           //nodesAvail = this.visitedNodesWithUnvisitedEdges();
        }



    };

    this.debCycle = function() {
        this.resetEdges();

        var allFinished = false;

        var walkNum = 0;

        while (!allFinished) {
            ++walkNum;
            var res = this.chooseStartCycle();
            var startNode = res[0];
            var edgePathInsert = res[1];

            var currNode;

            if (startNode) {
                currNode = startNode;
                currNode.visited = true;
            }
            else {
                allFinished = true;
                break;
            }

            var done = false;
            var newEdgePath = [];
            while (!done) {
                //var nonVisited = [];
                var r = currNode.getRandomNonVisitedEdge();
                if (r == -1) {
                    done = true;
                }
                else {
                    //this.edgePath.push(currNode.successors[r]);
                    newEdgePath.push(currNode.successors[r]);

                    currNode.successors[r].visited = true;
                    currNode.successors[r].walkNum = walkNum;
                    currNode = currNode.successors[r].targetNode;
                    currNode.visited = true;
                }
            }
            if (this.edgePath.length == 0) {
                this.edgePath = newEdgePath.map(function(el) {
                    return el;
                });
            }
            else {
                var svdThis = this;
                newEdgePath.forEach(function(el,i) {
                    svdThis.edgePath.splice(i + edgePathInsert,0,el);
                });
            }



        }





    };


    this.hamPath = function() {
        this.resetEdges();

       // var allFinished = false;

        var walkNum = 0;

      //  while (!allFinished) {
            ++walkNum;
            var res = this.chooseStartCycle(true);
            var startNode = res[0];
            var edgePathInsert = res[1];

            var currNode;

            if (startNode) {
                currNode = startNode;
                currNode.visited = true;
            }
            else {
               // allFinished = true;
              //  break;
            }

            var done = false;
            var newEdgePath = [];

            while (!done) {
                //var potentialNext = currNode.successors;
                var nonVisited = [];
                currNode.successors.forEach(function (el, i) {
                    if (!el.targetNode.visited) {
                        nonVisited.push(i);
                    }

                });
                if (nonVisited.length == 0) {
                    done = true;
                }
                else {
                    var r = getRandomInt(0, nonVisited.length - 1);
                    currNode.successors[nonVisited[r]].visited = true;
                    currNode.successors[nonVisited[r]].walkNum = walkNum;

                    newEdgePath.push(currNode.successors[nonVisited[r]]);
                    currNode = currNode.successors[nonVisited[r]].targetNode;
                    currNode.visited = true;

                }

            }

            if (this.edgePath.length == 0) {
                this.edgePath = newEdgePath.map(function (el) {
                    return el;
                });
            }
            else {
                var svdThis = this;
                newEdgePath.forEach(function (el, i) {
                    svdThis.edgePath.splice(i + edgePathInsert, 0, el);
                });
            }


        //}

    };



    this.debPath = function() {
        this.resetEdges();

        var allFinished = false;

        var walkNum = 0;

        while (!allFinished) {
            ++walkNum;
            var res = this.chooseStartCycle(true);
            var startNode = res[0];
            var edgePathInsert = res[1];

            var currNode;

            if (startNode) {
                currNode = startNode;
                currNode.visited = true;
            }
            else {
                allFinished = true;
                break;
            }

            var done = false;
            var newEdgePath = [];
            while (!done) {
                //var nonVisited = [];
                var r = currNode.getRandomNonVisitedEdge();
                if (r == -1) {
                    done = true;
                }
                else {
                    //this.edgePath.push(currNode.successors[r]);
                    newEdgePath.push(currNode.successors[r]);

                    currNode.successors[r].visited = true;
                    currNode.successors[r].walkNum = walkNum;
                    currNode = currNode.successors[r].targetNode;
                    currNode.visited = true;
                }
            }
            if (this.edgePath.length == 0) {
                this.edgePath = newEdgePath.map(function(el) {
                    return el;
                });
            }
            else {
                var svdThis = this;
                newEdgePath.forEach(function(el,i) {
                    svdThis.edgePath.splice(i + edgePathInsert,0,el);
                });
            }



        }





    };

    this.debPathOrig = function() {
        var first;
        var path = [];

        var r;

        //var definiteFirstFound = false;
       // var definiteLastFound = false;

        this.resetEdges();

        first = this.definiteFirst();
        if (first) {
            //definiteFirstFound = true;
        }
        else {

            r = getRandomInt(0,this.nodes.length - 1);
            first = this.nodes[r];
            //return [[],''];
        }

        var last = this.definiteLast();
        if (last) {
           // definiteLastFound = true;
        }



        var done = false;
        //path.push(first.dna);

        var curr = first;
        curr.visited = true;

        while (!done) {
           // var potentialNext = curr.successors;
           // var nonVisited = [];
            r = curr.getRandomNonVisitedEdge();
            if (r == -1) {
                    done = true;
            }
            else {
                path.push(curr.successors[r].edgeLabel());
                curr.successors[r].visited = true;
                curr = curr.successors[r].targetNode;
                curr.visited = true;
            }
        }

      //  var splitPoint = -1;

        var reconstructed = '';
        path.forEach(function(el,i) {
            if (i == 0) {
                reconstructed+=el;
            }
            else {
                reconstructed+=el.substring(el.length-1);
            }


        });

        if (reconstructed.length > this.dna.length) {
            reconstructed = reconstructed.substring(0, reconstructed.length - this.k + 1);
        }


        return [path,reconstructed];

    };

    this.buildNodesFromKmers = function() {
        var svdThis = this;

        if (this.graphType == DGraph.hamGraph) {
            svdThis = this;
            this.reads.forEach(function(kmer) {
                var node = new DNode(kmer);
                svdThis.nodes.push(node);
            });

            this.nodes.forEach(function(node) {
                var suffs = svdThis.findOverlapNodesWithSuffix(node);
                var prefs = svdThis.findOverlapNodesWithPrefix(node);

                suffs.forEach(function(sufNode) {
                    node.successors.push(new DEdge(node,sufNode));

                });
                prefs.forEach(function(prefNode) {
                    node.predecessors.push(prefNode);

                });


            });


        }
        else {

            this.reads.forEach(function (kmer) {
                var pref, suf;

                var dictKeyPref,dictKeySuf;
                if (svdThis.pairDist) {
                    pref = [kmer[0].substring(0, svdThis.k - 1), kmer[1].substring(0, svdThis.k - 1)];
                    suf = [kmer[0].substring(1, svdThis.k), kmer[1].substring(1, svdThis.k)];
                    dictKeyPref = pref[0] + pref[1];
                    dictKeySuf = suf[0] + suf[1];
                }
                else {
                    pref = kmer.substring(0, svdThis.k - 1);
                    suf = kmer.substring(1, svdThis.k);
                    dictKeyPref = pref;
                    dictKeySuf = suf;

                }

                var prefNode = svdThis.findNode(pref);
                if (!prefNode) {
                    prefNode = new DNode(pref);
                    //svdThis.nodes.push(prefNode);

                    svdThis.nodesDict[dictKeyPref] = prefNode;
                }

                var sufNode = svdThis.findNode(suf);


                if (!sufNode) {
                    sufNode = new DNode(suf);
                    //svdThis.nodes.push(sufNode);
                    svdThis.nodesDict[dictKeySuf] = sufNode;
                }

                prefNode.successors.push(new DEdge(prefNode, sufNode));
                sufNode.predecessors.push(prefNode);
            });

        }

        if (this.graphType == DGraph.hamGraph) {

        }
        else {

            for (var key in this.nodesDict) {
                this.nodes.push(this.nodesDict[key]);
            }
        }




        };

    this.buildNodesFromAdjList = function(adjList) {
        var svdThis = this;

        adjList.forEach(function (adj) {

            var spl = adj.split(' -> ');
            var pref = spl[0];

            var node = new DNode(pref);

            var matches = svdThis.findNodes(pref);
            if (matches.length > 0) {

                node.repeatNum = matches.length + 1;
            }
            svdThis.nodes.push(node);



        });

        adjList.forEach(function (adj,i) {

            var spl = adj.split(' -> ');
            //var pref = spl[0];
            var sufs = spl[1].split(',');
            var sufNode;



            var prevSuf;
            var rep = 0;
            sufs.forEach(function(suf,ii) {
                if (ii == 0) {
                    prevSuf = suf;
                    rep = 0;
                }
                if (suf == prevSuf) {
                    ++rep;
                }
                else {
                    prevSuf = suf;
                    rep = 1;
                }

                if (svdThis.graphType == DGraph.debGraph) {
                    //no repeat nodes (glue together)
                    rep = 1;
                }

                sufNode = svdThis.findNode(suf,rep);
                if (!sufNode) {
                    sufNode = new DNode(suf);
                    sufNode.repeatNum = rep;
                    svdThis.nodes.push(sufNode);
                }

                svdThis.nodes[i].successors.push(new DEdge(svdThis.nodes[i], sufNode));
                sufNode.predecessors.push(svdThis.nodes[i]);
            });

        });



            /*
            adjList.forEach(function (adj) {
                var spl = adj.split(' -> ');
                var pref = spl[0];
                var sufs = spl[1].split(',');
                var suf;

                var prefNode = svdThis.findNode(pref);

                if (!prefNode) {
                    prefNode = new DNode(pref);
                    svdThis.nodes.push(prefNode);
                }

                sufs.forEach(function(suf) {
                    var sufNode = svdThis.findNode(suf);
                    if (!sufNode) {
                        sufNode = new DNode(suf);
                        svdThis.nodes.push(sufNode);
                    }
                    prefNode.successors.push(new DEdge(prefNode, sufNode));
                    sufNode.predecessors.push(prefNode);

                });


            });
            */



    };


    this.initGraph = function() {
        var kmers;

        switch (this.sourceType)  {
            case DGraph.fromDna:
                kmers = kmerComposition(this.dna, this.k,true);
                this.reads = [];
                this.readPosns = [];
                var svdThis = this;
                kmers.forEach(function(el) {
                    svdThis.reads.push(el[0]);
                    svdThis.readPosns.push(el[1]);
                });
                //this.reads = kmers;
                this.buildNodesFromKmers();
                break;
            case DGraph.fromReads:
                //kmers = this.reads;

                this.buildNodesFromKmers();
                break;
            case DGraph.fromPairedDna:
                kmers = kmerPairedComposition(this.dna, this.k, this.pairDist);
                this.reads = kmers;
                this.buildNodesFromKmers();
                break;

            case DGraph.fromPairedReads:
                //kmers = this.reads;
                this.buildNodesFromKmers();
                break;

            case DGraph.fromAdjList:
                this.buildNodesFromAdjList(this.adjList);
                break;
            default:
                break;
        }




        /*
        if (reads) {

        }
        else {
            if (this.pairDist) {
                kmers = kmerPairedComposition(this.dna, this.k, this.pairDist);
            }
            else {
                kmers = kmerComposition(this.dna, this.k);
            }
        }
        */

        /*

        var svdThis = this;
        kmers.forEach(function (kmer) {
            var pref, suf;
            if (pairDist) {
                pref = [kmer[0].substring(0, svdThis.k - 1), kmer[1].substring(0, svdThis.k - 1)];
                suf = [kmer[0].substring(1, svdThis.k), kmer[1].substring(1, svdThis.k)];
            }
            else {
                pref = kmer.substring(0, svdThis.k - 1);
                suf = kmer.substring(1, svdThis.k);

            }

            var prefNode = svdThis.findNode(pref);
            var sufNode = svdThis.findNode(suf);
            if (!prefNode) {
                prefNode = new DNode(pref);
                svdThis.nodes.push(prefNode);
            }

            if (!sufNode) {
                sufNode = new DNode(suf);
                svdThis.nodes.push(sufNode);
            }

            prefNode.successors.push(new DEdge(prefNode, sufNode));
            sufNode.predecessors.push(prefNode);


        });

        */

        /*
        var balBefore = this.isBalanced();
        var balLinearBefore = this.isBalancedLinear();
        */

        if (this.makeCycle) {
            this.makeCyclical(); // make graph into cycle if not already
        }

        //var balAfter = this.isBalanced();
        //var balLinearAfter = this.isBalancedLinear();


    };


    this.initGraph();


}

DGraph.fromDna = 1;
DGraph.fromReads = 2;
DGraph.fromAdjList = 3;
DGraph.fromPairedDna = 4;
DGraph.fromPairedReads = 5;
DGraph.fromRna = 6; //not used for Graph, but used in transcrip/translate
DGraph.fromProtein = 7; //not used for Graph, but used in transcrip/translate
DGraph.fromSpectrum = 8; //not used for Graph, but used in peptide seq


DGraph.fromGrid = 9; // used in manhattan grid

DGraph.fromBreakpoint = 10; // used in chromosome distances

DGraph.fromDistanceMatrix = 11; //used in phylogeny
DGraph.fromFastaDna = 12;//used in phylogeny

DGraph.hamGraph = 1;
DGraph.debGraph = 2;

DGraph.seqTypeComp = 1;
DGraph.seqTypePath = 2;
DGraph.seqTypeCycle = 3;

DGraph.alignTypeGlobal = 1;
DGraph.alignTypeLocal = 2;
DGraph.alignTypeFitting = 3;
DGraph.alignTypeOverlap = 4;


//used for tree graph view styles
DGraph.styleUnrooted = 1;
DGraph.styleRooted = 2;



//probably should be in gen utilities

DGraph.infinity = 1000000;

function DBRect(x,y,w,h) {
    //origin is top left, not bottom left
    this.w = w;
    this.h = h;
    this.x = x;
    this.x = x;
    this.y = y;

    this.area = function() {
        return this.w * this.h;
    }

    this.containsPoint = function(pos) {
        //assuming y goes down not up
        return  ((pos.x >= this.x && (pos.x <= (this.x + this.w))
             && (pos.y >= this.y && (pos.y <= (this.y  + this.h)))));

    };
}

function DBPos(x,y) {
    this.x = x;
    this.y = y;

}

function DBEdgeView(startPos,angle,w,ctx,absEnd,lenFactor) {

    this.text = "blah";
    this.ctx = ctx;
    this.angle = angle;
    this.w = w;
    if (lenFactor == null) {
        this.len = 100;
    }
    else {
        this.len = lenFactor * w;
    }

    this.absEnd = absEnd;

    this.lenFactor = lenFactor;

    this.highlight = false;

    this.textOffsetX = 3;


    this.startPos = startPos;

    if (this.absEnd) {
         this.xLen = Math.abs(absEnd.x - startPos.x);
         this.yLen = Math.abs(absEnd.y - startPos.y);
         this.len = Math.sqrt((Math.pow(this.xLen,2) + Math.pow(this.yLen,2)));
         this.endPos = this.absEnd;
    }
    else {

        this.xLen = this.len * Math.cos(this.angle);
        this.yLen = this.len * Math.sin(this.angle);

        this.endPos = new DBPos(this.startPos.x + this.xLen, this.startPos.y - this.yLen);
    }


    this.display = function() {


        this.ctx.beginPath();
        var svdStrokeStyle = this.ctx.strokeStyle;
        if (this.highlight) {
            this.ctx.strokeStyle = "#FF0000";
        }
        this.ctx.moveTo(this.startPos.x,this.startPos.y);
        this.ctx.lineTo(this.endPos.x, this.endPos.y);
        this.ctx.stroke();
        this.ctx.strokeStyle = svdStrokeStyle;
        var txtPosXSt = Math.min(this.startPos.x, this.endPos.x);
        var txtPosXEnd = Math.max(this.startPos.x, this.endPos.x);
        var txtPosX = txtPosXSt +  ((txtPosXEnd - txtPosXSt) / 2);
        var txtPosYSt = Math.min(this.startPos.y, this.endPos.y);
        var txtPosYEnd = Math.max(this.startPos.y, this.endPos.y);
        var txtPosY = txtPosYSt +  ((txtPosYEnd - txtPosYSt) / 2);
        this.ctx.fillText(parseFloat(this.w).toFixed(3),txtPosX, txtPosY);


    }

}
/*
function DBEdgeView(ctx,startInd,endInd) {

    this.text = "blah";
    this.ctx = ctx;
    if (startInd) {
        this.startInd = startInd;
    }
    else {
        this.startInd = null;
    }

    if (endInd) {
        this.endInd = endInd;
    }
    else {
        this.endInd = null;
    }
    this.textOffsetX = 3;


    this.startPos = null;
    this.endPos = null;

    this.display = function() {


        this.ctx.beginPath();
        this.ctx.moveTo(this.startPos.x,this.startPos.y);
        this.ctx.lineTo(this.endPos.x, this.endPos.y);
        this.ctx.stroke();
        var txtPosXSt = Math.min(this.startPos.x, this.endPos.x);
        var txtPosXEnd = Math.max(this.startPos.x, this.endPos.x);
        var txtPosX = txtPosXSt +  ((txtPosXEnd - txtPosXSt) / 2);
        var txtPosYSt = Math.min(this.startPos.y, this.endPos.y);
        var txtPosYEnd = Math.max(this.startPos.y, this.endPos.y);
        var txtPosY = txtPosYSt +  ((txtPosYEnd - txtPosYSt) / 2);
        this.ctx.fillText(this.text,txtPosX, txtPosY);


    }

}
*/

function DBNodeView(ctx,r,style) {

    this.text = "blah";
    this.tag = -1;
    this.ctx = ctx;
    this.style = this.style ? this.style : DGraph.styleUnrooted;
    this.startAngleRad = 0;
    this.incRad = 90 / 360 * 2 * Math.PI;
    this.isLeaf = false;

    this.highlight = false;
    this.edgeHighlight = false;

    if (r) {
        this.r = r;
    }
    else {
        this.r = new DBRect(0,0,10,10);
    }

    this.absPos = false; //absolute positioning

    this.textOffsetX = 3;

    this.edgeViews = [];

    this.centre = function() {
       return new DBPos(this.r.x + (this.r.w/2),this.r.y + (this.r.h/2));
    };

    this.addEdge = function(nodeTo,w,lenFactor) {
        var ang = 0;
        var leafAdj = 0;


        if ((this.style == DGraph.styleRooted) && nodeTo.isLeaf) {
            leafAdj = 20 / 360 * 2 * Math.PI;
        }

        if (this.edgeViews.length == 0) {
            ang = this.startAngleRad - leafAdj;
        }
        else if (this.edgeViews.length == 1) {
            ang =  this.style == DGraph.styleRooted ? this.startAngleRad - this.incRad + leafAdj : this.startAngleRad + this.incRad;
        }
        else {
            ang = this.startAngleRad + (Math.PI * 1.25);
        }

        var edge;
        if (nodeTo.absPos) {
            var absEnd = new DBPos(nodeTo.r.x + (10/2),nodeTo.r.y + (10/2));
            edge = new DBEdgeView(this.centre(), ang, w, this.ctx,absEnd);
        }
        else {
            edge = new DBEdgeView(this.centre(), ang, w, this.ctx,null,lenFactor);//(this.x,this.y,this.len,ang,this.edges.length,w);
            nodeTo.r.x = edge.endPos.x - (10 / 2);
            nodeTo.r.y = edge.endPos.y - (10 / 2);
        }
        nodeTo.startAngleRad = this.style == DGraph.styleRooted ? -45 / 360 * 2 * Math.PI  :  ang - (this.incRad/2);

        if ((this.edgeHighlight) && (nodeTo.edgeHighlight)) {
            edge.highlight = true;
        }

        this.edgeViews.push(edge);
        //nodeTo.edges.push(edge);
    }


    this.display = function() {
        this.ctx.beginPath();
        var svdStrokeStyle = this.ctx.strokeStyle;
        if (this.highlight) {
            this.ctx.strokeStyle = "#FF0000";
        }
        this.ctx.arc(this.centre().x,this.centre().y,2,0,Math.PI*2,true);
        this.ctx.stroke();
        this.ctx.strokeStyle = svdStrokeStyle;

        var svdFont = this.ctx.font;

        this.ctx.font = '8px Arial';
        if (this.edgeViews.length == 0) {
            //no edges going out


            this.ctx.font = 'bold ' + this.ctx.font;
        }

        this.ctx.fillText(this.text,this.r.x + this.textOffsetX,this.r.y);
        this.ctx.font = svdFont;

        this.displayEdges();




    };

    this.displayEdges = function() {
        for (var i = 0;i < this.edgeViews.length; ++i) {
            this.edgeViews[i].display();
        }

    };

    this.nodeCentre = function() {
        return new DBPos(this.r.x,this.r.y);
    };

}

function DBTestSimpleController() {
    this.nodeViews = function() {
        var nvs = [];
        for (var i = 0;i < 6;++i) {
            var nv = new DBNodeView();
            nv.text = 'n' + i;
            nvs.push(nv);
        }

        return nvs;
        
    }

    this.numNodes = function() {
        return 6;
    }
    
    this.nodeViewOrder = function() {
        var pairs = [];
        pairs.push([2,4,3]);
        pairs.push([2,5,2]);
        pairs.push([2,3,4]);
        pairs.push([3,1,3]);
        pairs.push([3,0,7]);
        return pairs;




    }
    
}

function DBGraphViewController(g) {

    //"implements" gvDataSource
    //-numNodes: returns number of nodes in graph
    //-nodeViews: returns a view node for every node in graph
    //-nodeViewOrder: returns an array of node to successor node index pairs, plus weights
    //-longestPathToLeaf: returns longest path from first view node, ie central or root,  to all leaves
    //-preferredStyle: returns style
    //-getComments: returns comments about graph

    this.g = g;
    this.numNodes = function() {
        return g.numNodes();
    }
    

    this.nodeViews = function() {
        var vNodes = [];
        var svdThis = this;
        this.g.nodes.forEach(function(el,i) {
            var vNode = new DBNodeView();
            vNode.text = el.label;
            vNode.highlight = el.freshNode ? true : false;

            if (svdThis.g.builder.nodeBelongsToFreshEdge(el)) {
                vNode.edgeHighlight = true;
            }

            if ((el.sequence.length > 0) && (el.sequence != el.label))  {
                var seq = el.sequence;
                var thresh = (svdThis.preferredStyle() == DGraph.styleRooted) ? 8 : 20;
                if (el.sequence.length > thresh) {
                    seq = el.sequence.substring(0, thresh) + '..';
                }

                vNode.text += ':';
                vNode.text += seq;
            }
            vNode.tag = i;
            vNode.isLeaf = el.isLeaf(svdThis.g.checkDirected());
            vNodes.push(vNode);
        });
        return vNodes;

    }

    /*
    this.nodeViewAtIndex = function(i) {
        var nv = new DBNodeView();
        nv.text = g.nodes[i].label;
        for (var j = 0;j < g.nodes[i].successors.length;++j) {
            var edV = new DBEdgeView();
            edV.text = g.nodes[i].successors[j].edgeWeight;
            edV.startInd = g.getNodeIndexForNode(g.nodes[i].successors[j].sourceNode);
            edV.endInd = g.getNodeIndexForNode(g.nodes[i].successors[j].targetNode);
            nv.edgeViews.push(edV);
        }
        return nv;
    }
    */


    this.longestPathToLeaf = function() {

        var lvs = g.leaves();

        var startNode = this.findStartNode();

        var longest = -1;
        lvs.forEach(function(leafNode) {
            var pathLen = parseFloat(g.findPathBetweenNodes(startNode,leafNode,startNode,null).pathLen);
            if (pathLen > longest) {
                longest = pathLen;
            }
        });

        return longest;


    };

    //private
    this.findStartNode = function() {
        var startNodeLab;
        var startNode;

        if (g.isRooted()) {
            startNodeLab = g.findRoot().label;
        }
        else {
            startNodeLab = g.centralNode();
        }
        // var ind = g.getNodeIndexForNode(centralNode);

        var done = false;
        var startNode = g.getNodeFromLabel(startNodeLab);

        return startNode;


    };

    this.nodeViewOrder = function() {

        var startNode = this.findStartNode();

        var nodesToProcess = [startNode];

        var nodePairs = [];


        var nodesProcessed = [];

        var count = 0;

        while (nodesToProcess.length > 0) {

            var newNodesToProcess = [];
           // console.log('nodes to proc: ' + nodesToProcess);
            nodesToProcess.forEach(function(node) {
                //var nexts = findInAdjList(nodeNum,nodesProcessed);
                //console.log('nexts: ' + nexts);
                //var node = g.nodes[nodeNum];
                node.getSuccessors().forEach(function(nx) {
                    var targNode = nx.getTargetNode(node);
                   if (nodesProcessed.indexOf(targNode) > -1) {

                   }
                    else {
                       newNodesToProcess.push(targNode);
                       nodePairs.push([g.nodes.indexOf(node),g.nodes.indexOf(targNode),nx.edgeWeight]);
                   }
                });

            });
            nodesProcessed = nodesProcessed.concat(nodesToProcess);
            console.log('nodes proc: ' + nodesProcessed);
            console.log('new nodes to proc: ' + newNodesToProcess);
            nodesToProcess = newNodesToProcess;

            ++count;
            //if (count > 5) break;


        }

        return nodePairs;

    }

    this.preferredStyle = function() {
        if (this.g.isRooted()) {
            return DGraph.styleRooted;
        }
        else {
            return DGraph.styleUnrooted;
        }

    };

    this.getComments = function() {
        return g.comments;
    }



}

function DBGraphView(canv,viewPort,gvDataSource,ctx) {
    //canv = canvas element

    // graphDataSource "implements" DBGraphSourceIF

    var svdThis = this;

    this.canv = canv;
    if (ctx) {
        this.ctx = ctx;
    }
    else {
        this.ctx = canv.getContext('2d');
    }
    this.gvDataSource = gvDataSource;
    this.margin = 3;

    this.viewPort = viewPort;
    

    this.nvBeingDragged = null;

    this.nodesOffset = 20;

    this.nvs = [];

    //this.ctx.translate(0,-150); perhaps work out topmost point and translate to it


    this.r = new DBRect(this.margin,this.margin,this.ctx.canvas.clientWidth - (this.margin * 2),this.ctx.canvas.clientHeight - (this.margin*2));
   // this.r = new DBRect(this.margin,this.margin,this.viewPort.w - (this.margin * 2),this.viewPort.h  - (this.margin*2));





    this.arrange = function(order) {
        var svdThis = this;

        var longest = this.gvDataSource.longestPathToLeaf();
        var lenFactor;
        if (gvDataSource.preferredStyle() == DGraph.styleRooted) {
            //lenFactor = 500 / longest; // make longest 300
            lenFactor = null; //uncomment above line if reimplementing this
        }
        else {
            lenFactor = null;
        }


        this.nvs.forEach(function(nv) {
            nv.edgeViews = [];
        });

        order.forEach(function(pair) {
            svdThis.nvs[pair[0]].addEdge(svdThis.nvs[pair[1]],pair[2],lenFactor);
        });

        var highest = this.findHighestNodeView();
        var leftest = this.findLeftmostNodeView();
        var lowest =  this.findLowestNodeView();

        var amtToOffset = 0;
        ///if (highest.r.y < this.margin) {

        //}
        //else {

        amtToOffset = highest.r.y - this.margin;
        this.nvs[order[0][0]].r.y -= amtToOffset;
        this.nvs[order[0][0]].r.y += 10;

        this.nvs.forEach(function(nv) {
            nv.edgeViews = [];
        });

        order.forEach(function(pair) {
            svdThis.nvs[pair[0]].addEdge(svdThis.nvs[pair[1]],pair[2],lenFactor);
        });

        //this.arrange(order);
        //}

        if (leftest.r.x < this.margin) {
            amtToOffset = leftest.r.x - this.margin;

            this.nvs[order[0][0]].r.x -= amtToOffset;
            this.nvs[order[0][0]].r.x += 10;
        }

        //}

        this.nvs.forEach(function(nv) {
            nv.edgeViews = [];
        });

        order.forEach(function(pair) {
            svdThis.nvs[pair[0]].addEdge(svdThis.nvs[pair[1]],pair[2],lenFactor);
        });
        //this.arrange(order);


        //Ultrametric - all leaf nodes aligned at bottom
        // only attempt for small trees

        if ((gvDataSource.preferredStyle() == DGraph.styleRooted) && this.nvs.length < 20) {
            this.nvs.forEach(function (nv) {
                if (nv.edgeViews.length == 0) {
                    if (nv.r.y < lowest.r.y) {
                        nv.r.y = lowest.r.y;
                        nv.absPos = true;
                    }
                }
            });

            this.nvs.forEach(function (nv) {
                nv.edgeViews = [];
            });

            order.forEach(function (pair) {
                svdThis.nvs[pair[0]].addEdge(svdThis.nvs[pair[1]], pair[2],lenFactor);
            });
        }


    };

    this.findLeftmostNodeView = function() {

        var leftestNV = null;
        var leftestX = DGraph.infinity;

        this.nvs.forEach(function(nv) {
            if (nv.r.x < leftestX) {
                leftestX = nv.r.x;
                leftestNV = nv;
            }
        });

        return leftestNV;
    };

    this.findHighestNodeView = function() {
      var highestNV = null;
      var highestY = DGraph.infinity;

      this.nvs.forEach(function(nv) {
          if (nv.r.y < highestY) {
              highestY = nv.r.y;
              highestNV = nv;
          }
      });

      return highestNV;

    };

    this.findLowestNodeView = function() {
        var lowestNV = null;
        var lowestY = DGraph.infinity * -1;

        this.nvs.forEach(function(nv) {
            if (nv.r.y  > lowestY) {
                lowestY = nv.r.y;
                lowestNV = nv;
            }
        });

        return lowestNV;

    };


    this.refreshDisplay = function() {
        var order = this.gvDataSource.nodeViewOrder();
        this.arrange(order);

        this.ctx.clearRect(0, 0, this.canv.width, this.canv.height);
        this.ctx.strokeRect(this.r.x,this.r.y,this.r.w,this.r.h);

        this.nvs.forEach(function(nv) {
            nv.display();
        });


    };

    this.display = function() {
       //this.ctx.translate(0,-200);
       this.ctx.clearRect(0, 0, this.canv.width, this.canv.height);
       this.ctx.strokeRect(this.r.x,this.r.y,this.r.w,this.r.h);

       if (this.gvDataSource.getComments().length > 0) {
           this.ctx.fillText(this.gvDataSource.getComments(),15,15);
       }
       var cols = 2;
       var rows = this.gvDataSource.numNodes() / cols;
       var nodeSpacingX = (this.r.w - this.nodesOffset) / cols;
       var nodeSpacingY = (this.r.h - this.nodesOffset) / rows;// this.gvDataSource.numNodes();
       //nodeSpacing = 30;

       var order = this.gvDataSource.nodeViewOrder();
       var row = 0;var col = 0;

       this.nvs = this.gvDataSource.nodeViews();
       var svdThis = this;



       this.nvs.forEach(function(nv) {
          nv.ctx = svdThis.ctx;
          nv.style = svdThis.gvDataSource.preferredStyle();
       });

       //this.nvs[order[0][0]].r = new DBRect(this.r.x + (this.r.w / 2),this.r.y  + (this.r.h/2),10,10); //central node
        this.nvs[order[0][0]].r = new DBRect(this.viewPort.x + (this.viewPort.w / 2),this.viewPort.y  + (this.viewPort.h/2),10,10); //central n

       if (this.gvDataSource.preferredStyle() == DGraph.styleRooted) {
          this.nvs[order[0][0]].startAngleRad = -45 / 360 * 2 * Math.PI;
          this.nvs[order[0][0]].r.y = 20;
       }

       this.arrange(order);
        /*
       var highest = this.findHighestNodeView();
       var leftest = this.findLeftmostNodeView();

        var amtToOffset = 0;
       ///if (highest.r.y < this.margin) {

        //}
        //else {
        amtToOffset = highest.r.y - this.margin;
        this.nvs[order[0][0]].r.y -= amtToOffset;
        this.nvs[order[0][0]].r.y +=10;
        this.nvs.forEach(function(nv) {
            nv.edgeViews = [];
        });
        this.arrange(order);
        //}

        if (leftest.r.x < this.margin) {
            amtToOffset = leftest.r.x - this.margin;

            this.nvs[order[0][0]].r.x -= amtToOffset;
            this.nvs[order[0][0]].r.x += 10;
        }

        //}

        this.nvs.forEach(function(nv) {
            nv.edgeViews = [];
        });
        this.arrange(order);

       */

       this.nvs.forEach(function(nv) {
           nv.display();
       });
        /*
       for (var i = 0;i < this.gvDataSource.numNodes();++i) {
           var nv = this.gvDataSource.nodeViewAtIndex(i);
           nv.ctx = this.ctx;
           for (var j = 0;j < nv.edgeViews.length;++j) {
               nv.edgeViews[j].ctx = this.ctx;
           }
           nv.r = new DBRect(this.r.x + (col * nodeSpacingX) + this.nodesOffset,this.r.y + (row * nodeSpacingY) + this.nodesOffset,10,10);


           nv.display();
           ++col;
           if (col > 1) {
             col = 0;
             ++row;
           }
           nvs.push(nv);
       }

        for (var i = 0;i < nvs.length;++i) {
            nv = nvs[i];
            for (var j = 0;j < nv.edgeViews.length;++j) {
                var ev = nv.edgeViews[j];
                ev.startPos = nvs[ev.startInd].nodeCentre();
                ev.endPos = nvs[ev.endInd].nodeCentre();
                ev.display();
            }

        }
        */



    };

    this.determineNodeViewClicked = function(clickEvent) {
        //determine which node has been clicked on depending on click event
        var pos = new DBPos(clickEvent.offsetX,clickEvent.offsetY);

        var nvClicked = null;
        
        this.nvs.forEach(function(nv) {
            if (nv.r.containsPoint(pos)) {
                console.log('node clicked: ' + nv.text);
                nvClicked = nv;
            }
        });
        
        return nvClicked;


    };


    this.mouseDown = function(e) {
        console.log('mouse down');
        //this.nvBeingDragged = this.determineNodeViewClicked.call(this,e);
        this.nvBeingDragged = this.determineNodeViewClicked(e);

    };

    this.mouseUp = function(e) {
        console.log('mouse up');
        if (this.nvBeingDragged) {
            this.nvBeingDragged.r.x = e.offsetX;
            this.nvBeingDragged.r.y = e.offsetY;
            this.nvBeingDragged.absPos = true;
            this.refreshDisplay();

            this.nvBeingDragged = null;
        }

    };


    this.mouseMove = function(e) {
        //console.log('mouse move');
        if (this.nvBeingDragged) {
            this.nvBeingDragged.r.x = e.offsetX;
            this.nvBeingDragged.r.y = e.offsetY;
            this.nvBeingDragged.absPos = true;
            this.refreshDisplay();


        }


    };

    this.mouseClicked = function(e) {

        var x = 1;
    };


   // this.canv.graphView = this;

   // this.canv.addEventListener("click", svdThis.mouseClicked);

    //this.canv.addEventListener("click", svdThis.mouseClicked.bind(null, e,test));
    this.canv.addEventListener("click", this.mouseClicked.bind(this), false);
    this.canv.addEventListener("mousedown",this.mouseDown.bind(this),false);
    this.canv.addEventListener("mouseup", this.mouseUp.bind(this),false);
    this.canv.addEventListener("mousemove", this.mouseMove.bind(this),false);



}

function DBEdge(sourceNode,targetNode,w,dirFlag) {


    this.sourceNode = sourceNode;
    this.targetNode = targetNode;

    this.visited = false;
    this.realigned = false; // used when making edge directed temporarily, ie adding root to tree
                            // this var used to remember original undirected source/target node
    this.edgeWeight = w ? w : 0;

    this.directedEdge = (dirFlag == null) ? true : dirFlag;

    this.freshEdge = false;

    this.edgeLabel = function () {
        var joiner = this.directedEdge ? '->' : '<->';
        return this.sourceNode.label + joiner + this.targetNode.label;
    }

    this.setWeight = function (w) {
        this.edgeWeight = w;
    }

    this.setDirected = function (dirFlag,realignedTarg) {
        this.directedEdge = dirFlag;

        if (realignedTarg) { //swap to orient in correct direction towards targ
            if (realignedTarg === this.targetNode) {

            }
            else {
                this.sourceNode = this.targetNode;
                this.targetNode = realignedTarg;
                this.realigned = true;
            }
        }
    }


    this.isDirected = function () {
        /*
         var succ = this.targetNode.successors;

         if (succ.length == 0) {
         return true;
         }

         var dirFlag = true;
         var svdThis = this;
         succ.forEach(function(el) {
         if (el.targetNode === svdThis.sourceNode) {
         dirFlag = false;
         }
         });
         return dirFlag;
         }
         */
        return this.directedEdge;

    }

    this.getTargetNode = function(nodeStart) {
        //get target from viewpoint of start node - for undirected edges, could be "backward"

        if (this.isDirected()) {
            return this.targetNode;
        }
        else {
            return (nodeStart === this.sourceNode) ? this.targetNode : this.sourceNode;
        }

    }

    this.getSourceNode = function(nodeEnd) {
        //get source from viewpoint of end node - for undirected edges, could be "backward"

        if (this.isDirected()) {
            return this.sourceNode;
        }
        else {
            return (nodeEnd === this.targetNode) ? this.sourceNode : this.targetNode;
        }

    }



}

function DBTreeEdge(sourceNode,targetNode,w,dirFlag) {
    DBEdge.apply(this, [sourceNode, targetNode, w,dirFlag]);
    
    this.isLimb = function() {
        if (sourceNode.isLeaf() || targetNode.isLeaf())  {
            return true;
        }
        
        return false;
    }
    
    this.isInternalEdge = function() {
        if  (sourceNode.isLeaf() || targetNode.isLeaf())  {
            return false;
            
        }
        return true;
    }

    this.neighbouringEdges = function() {
        //return neighbouring edges

        if (this.isDirected()) {
            return null;
            //only implemented for undirected
        }



        var nFrom = this.sourceNode;
        var nTo = this.targetNode;
        var edNodes = [nFrom,nTo];

        var fromNeighbours = nFrom.getSuccessors(nTo);
        var toNeighbours = nTo.getSuccessors(nFrom);

        return [fromNeighbours,toNeighbours];



    };

    this.removeAllConnectingEdges = function() {

        var svdThis = this;

        //remove all edges  from source node except this one
        var targs = [];
        for (var i = this.sourceNode.edges.length -1;i >= 0;--i) {
            var targ = this.sourceNode.edges[i].getTargetNode(this.sourceNode);
            if (targ == this.targetNode) {

            }
            else {
                this.sourceNode.edges.splice(i,1);
                targs.push(targ);
            }
        }
        targs.forEach(function(t) {
            for (var i = t.edges.length -1;i >= 0;--i) {
                var targ = t.edges[i].getTargetNode(t);
                if (targ == svdThis.sourceNode) {
                    t.edges.splice(i,1);

                }
            }

        });


        //remove all edges  from target node except this one
        targs = [];
        for (var i = this.targetNode.edges.length -1;i >= 0;--i) {
            var targ = this.targetNode.edges[i].getTargetNode(this.targetNode);
            if (targ == this.sourceNode) {

            }
            else {
                this.targetNode.edges.splice(i,1);
                targs.push(targ);
            }
        }
        targs.forEach(function(t) {
            for (var i = t.edges.length -1;i >= 0;--i) {
                var targ = t.edges[i].getTargetNode(t);
                if (targ == svdThis.targetNode) {
                    t.edges.splice(i,1);

                }
            }

        });



    };

}



function DBNode(label) { //basic node
    this.label = label;

    this.visited = false;

    this.cycleNum = 0;

    this.freshNode = false; //used to highlight newly added nodes

    //this.successors = [];
    //this.predecessors = [];

    this.edges = [];


    this.getSuccessorNodes = function(onlyNonVisited) {
         var succEdges = this.getSuccessors();
        var succNodes = [];
        var svdThis = this;
        succEdges.forEach(function(edge) {
            if (edge.sourceNode === svdThis) {
                succNodes.push(edge.targetNode);
            }
            else {
                succNodes.push(edge.sourceNode);
            }
        });

        return succNodes;


    };

    this.getSuccessors = function (awayFromNode) {
        //gets successor edges of a node. If awayFromNode is supplied for undirected graph, only returns
        //edges in other direction of awayFromNode
        
        
        var sucs = [];
        for (var i = 0; i < this.edges.length; ++i) {
            var edge = this.edges[i];
            if (edge.isDirected()) {
                if (edge.sourceNode === this) {
                    sucs.push(edge); //only outgoing edges
                }
            }
            else {
                if (awayFromNode && ((edge.sourceNode == awayFromNode) || (edge.targetNode == awayFromNode))) {
                    
                }
                else {
                    sucs.push(edge); //all edges are outgoing (and incoming)
                }
            }
        }

        return sucs;

    };

    this.getSuccessors = function (awayFromNode) {
        //gets successor edges of a node. If awayFromNode is supplied for undirected graph, only returns
        //edges in other direction of awayFromNode


        var sucs = [];
        for (var i = 0; i < this.edges.length; ++i) {
            var edge = this.edges[i];
            if (edge.isDirected()) {
                if (edge.sourceNode === this) {
                    sucs.push(edge); //only outgoing edges
                }
            }
            else {
                if (awayFromNode && ((edge.sourceNode == awayFromNode) || (edge.targetNode == awayFromNode))) {

                }
                else {
                    sucs.push(edge); //all edges are outgoing (and incoming)
                }
            }
        }

        return sucs;

    };


    this.getPredecessors = function () {
        var preds = [];
        for (var i = 0; i < this.edges.length; ++i) {
            var edge = this.edges[i];

            if (edge.isDirected()) {
                if (edge.targetNode === this) {
                    preds.push(edge.sourceNode); //only outgoing edges
                }
            }
            else {
                var predNode = (edge.sourceNode == this) ? edge.targetNode : edge.sourceNode;
                preds.push(predNode); //all edges are incoming (and outgoing)
            }
        }
        return preds;

    };

    this.outDegree = function () {

        //return this.successors.length;

        //new regime:
        return this.getSuccessors().length;

    };

    this.inDegree = function () {

       // return this.predecessors.length;


        //new regime:
        return this.getPredecessors().length;
    };

    this.degree = function () {


        /*
        if (this.successors.length > 0) {
            if (this.successors[0].isDirected()) {
                return this.outDegree(); // directed, so don't count twice
            }
        }

        return this.inDegree() + this.outDegree();
        */


        //new regime:
        return this.edges.length;

    };

    this.resetEdgesVisited = function () {
        this.successors.forEach(function (edge) {
            edge.visited = false;
        });

        this.visited = false; //node visited flag

    };

    this.getRandomNonVisitedEdge = function () {

        var nonVisited = [];
        this.getSuccessors().forEach(function (el, i) {
            if (!(el.visited)) {
                nonVisited.push(i);
            }

        });
        if (nonVisited.length == 0) {
            return -1;
        }
        else {
            var r = getRandomInt(0, nonVisited.length - 1);
            return nonVisited[r];
        }
    };

    this.connectsTo = function (node2) {
        //checks if node2 is a successor neighbour of this node

        /*
        for (var i = 0; i < this.successors.length; ++i) {
            if (this.successors[i].targetNode == node2) {
                return true;

            }
        }
        */
        //new regime:
        for (var i = 0;i < this.edges.length;++i) {
            if (this.edges[i].directedEdge) {
                if (this.edges[i].targetNode === node2) {
                    return true;
                }
            }
            else {
                if ((this.edges[i].targetNode === node2) || (this.edges[i].sourceNode === node2)) {
                    return true;
                }

            }

        }
        return false;

    }

    this.connectsBothWaysTo = function (node2) {
        if (this.connectsTo(node2) && node2.connectsTo(this)) {
            return true;
        }
        return false;
    }
    
    
    this.removeEdge = function(ed) {
      
        for (var i = this.edges.length-1;i >=0;--i) {
           // if (this.edges[i] === ed) {  // not sure why this isn't working
            if ((this.edges[i].sourceNode == ed.sourceNode) && (this.edges[i].targetNode == ed.targetNode)) {
                this.edges.splice(i,1);
                break;
            }
        }
    };


}


function DBTreeNode(label) {
    var labSplit = label.split(';');
    
    DBNode.apply(this,[labSplit[0]]);
    
    this.sequence = labSplit.length > 1 ? labSplit[1] : '';

    this.parsimonyScore = null;
   
    this.isLeaf = function(isDirected) {
        if (isDirected) {
            return (this.outDegree() == 0);
        }
        else {
            return (this.outDegree() == 1);
        }
    };

    this.nodeOrderToLeaves = function() {
        //returns list of nodes on way from curr node to leaves
        var currNodes = [this];
        var allNodes = [this];
        var edgesProcessed = [];
        while (currNodes.length > 0) {
            var newNodes = [];

            currNodes.forEach(function(node) {
               var nextEdges = node.getSuccessors();
                nextEdges.forEach(function(ed) {
                    if (edgesProcessed.indexOf(ed) > -1) {

                    }
                    else {
                        edgesProcessed.push(ed);
                        var nextNode = ed.getTargetNode(node);
                        newNodes.push(nextNode);
                        allNodes.push(nextNode);
                    }
                });
            });
            currNodes = newNodes;
        }

        return allNodes;


    };
    
    this.bestParsimonyBasedOnChar = function(ch,childScore) {
        //if ch is empty, return absolute best. Used for root
        var bestScore = DGraph.infinity;
        var bestCh = '';

        if (childScore == null) {
            console.log(' childscore null in bestparbased. ch is: ' + ch);
        }

       // console.log('childscore len: ' + childScore.length + ' ' +  childScore);

        var keys = Object.keys(childScore).sort();
        keys.forEach(function(key) {
        //Object.keys(childScore).forEach(function(key) {
            var currScore;
            if (ch === '') {
                currScore = childScore[key];
            }
            else {
                currScore =  key == ch ? childScore[key] : childScore[key] + 1;
            }

            if (currScore < bestScore) {
                bestScore = currScore;
                bestCh = key;
            }
        });

        return [bestScore,bestCh];
        
    };

    this.calcParsimonyScore = function(child1Score,child2Score,useAminos) {
        var score = [];
        var len = child1Score.length;

        console.log('node: ' + this.label + ' child1 sc len: ' + child1Score.length + ' chid2 sc len: ' + child2Score.length);


        if (this.label == '20') {
            console.log('node 20 children');
        }
        //console.log('child1 sc: ' + child1Score + ' child2sc: ' + child2Score);

        //If useAminos flag is set, uses Amino acid alphabet instead of Bases alphabet
        if (useAminos == null) {
            useAminos = false;
        }

        var alphabet = c_Bases;
        if (useAminos) {
           alphabet = c_Aminos;
        }


        var svdThis = this;

        if (child1Score == null) {
            console.log('whole child1score null');
        }
        if (child2Score == null) {
            console.log('whole child2score null');
        }


        for (var i = 0;i < len;++i) {
            var scoreDict = {};
            if (child1Score[i] == null) {
                console.log('ch1, i: ' + i);
            }
            if (child2Score[i] == null) {
                console.log('ch2, i: ' + i);
            }

            alphabet.forEach(function (ch) {
               //check for chil2score being null. seems to be a problem here
                scoreDict[ch] = svdThis.bestParsimonyBasedOnChar(ch,child1Score[i])[0] + svdThis.bestParsimonyBasedOnChar(ch,child2Score[i])[0];
            });
            score.push(scoreDict);
        }
        if (this.label == 'I1') {
            console.log(' score len I1: ' + score.length);

        }
        return score;
    };

}
    



function DNodeBuilder() {

    this.createNode = function(lab) {
       return new DBNode(lab);
    }
}

function DTreeNodeBuilder() {
    DNodeBuilder.apply(this,[]);
    
    this.createNode = function(lab) {
        return new DBTreeNode(lab);
    }

}

function DEdgeBuilder() {

    this.createEdge = function(sourceNode,targetNode,w,dirFlag) {
        return new DBEdge(sourceNode,targetNode,w,dirFlag);
    }
}

function DTreeEdgeBuilder() {
    DEdgeBuilder.apply(this,[]);

    this.createEdge = function(sourceNode,targetNode,w,dirFlag) {
        return new DBTreeEdge(sourceNode,targetNode,w,dirFlag);
    }

}


function DGraphBuilder(source,nodeBuilder, edgeBuilder) {
    
    this.nodes = [];
    this.source = source;

    this.progress = []; //used to store state of graphs as being built
    

    this.internalNodePrefix = 'I';
    this.internalNodeNextNum = 1;

    if (!nodeBuilder) {
        this.nodeBuilder  = new DNodeBuilder();
    }
    else {
        this.nodeBuilder = nodeBuilder;
    }

    if (!edgeBuilder) {
        this.edgeBuilder  = new DEdgeBuilder();
    }
    else {
        this.edgeBuilder = edgeBuilder;
    }
   
    this.buildNodes = function() {
        //subclass this

        this.buildNodesFromAdjList(this.source.split('\n'));
        return this.nodes;
    };

    this.addNodeBetween = function(node1,node2,dist1,dist2) {

        //adds new node between two neighbours and joins up

        var directedFlag;
        var direction;
        if (node1.connectsBothWaysTo(node2)) {
            directedFlag = false;
        }
        else {
            directedFlag = true;
            direction = node1.connectsTo(node2) ? 'R' : 'L';
        }



       // var newNode = this.addNode(this.internalNodePrefix + this.internalNodeNextNum);
       // ++this.internalNodeNextNum;
        var newNode = this.addInternalNode();

        //remove existing edge(s) connecting node1, node2

        //new regime:
        for (var i = node1.edges.length - 1;i >= 0;--i) {
            var ed = node1.edges[i];
            if (ed.getTargetNode(node1) === node2){
                node1.edges.splice(i,1);
            }
        }

        for (var i = node2.edges.length - 1;i >= 0;--i) {
            var ed = node2.edges[i];
            if ((ed.getTargetNode(node2) === node1) ){
                node2.edges.splice(i,1);
            }
        }


        //old regime
        /*
        for (var i = 0;i < node1.successors.length;++i) {
            if (node1.successors[i].targetNode == node2) {
                node1.successors.splice(i,1);
                break;
            }
        }
        for (i = 0;i < node2.successors.length;++i) {
            if (node2.successors[i].targetNode == node1) {
                node2.successors.splice(i,1);
                break;
            }
        }

        for (i = 0;i < node2.predecessors.length;++i) {
            if (node2.predecessors[i] == node1) {
                node2.predecessors.splice(i,1);
                break;
            }
        }
        for (i = 0;i < node1.predecessors.length;++i) {
            if (node1.predecessors[i] == node2) {
                node1.predecessors.splice(i,1);
                break;
            }
        }
        */

        if (directedFlag) {
            if (direction == 'R') {
                this.connectNodes(node1, newNode, directedFlag, dist1);
                this.connectNodes(newNode, node2, directedFlag, dist2);
            }
            else {
                this.connectNodes(node2, newNode, directedFlag, dist2);
                this.connectNodes(newNode, node1, directedFlag, dist1);
            }
        }
        else {
            this.connectNodes(node1, newNode, directedFlag, dist1);
            this.connectNodes(newNode, node2, directedFlag, dist2);
        }

        return newNode;


    };


    this.connectNodes = function(node1,node2,dirFlag,w) {
        if (!w) {
            w = 0;
        }

       // node1.successors.push(new DBEdge(node1,node2,w));

        //If directed, first check if a directed edge already exits the other way. If so, join and make undir

        var processed = false;
        node2.getSuccessors().forEach(function(ed) {
           if ((ed.isDirected()  && (ed.getTargetNode(node2) === node1) )){
               ed.setDirected(false);
               processed = true;

           }
        });

        if (processed) {
            return;
        }

        var newEdge = this.edgeBuilder.createEdge(node1,node2,w,dirFlag);

        /*
        node1.successors.push(newEdge);
        node2.predecessors.push(node1);
        */

        //new regime
        node1.edges.push(newEdge);
        node2.edges.push(newEdge);

        if (dirFlag) {

        }
        else {
            //node2.successors.push(new DBEdge(node2,node1,w));
            /*
            node2.successors.push(this.edgeBuilder.createEdge(node2,node1,w,dirFlag));
            node1.predecessors.push(node2)
            */

        }

        return newEdge;


    };
    
    
    this.swapEdges = function(ed1,ed2,fromNode1,fromNode2) {
      //  assumes fromNode1 is attached to ed1, fromNode2 is attached to ed2.
      // attaches ed1 to fromNode2 and ed2 to fromNode1
        
      var targNode1 = ed1.getTargetNode(fromNode1);
      var targNode2 = ed2.getTargetNode(fromNode2);
      
      fromNode1.removeEdge(ed1);
      targNode1.removeEdge(ed1);  
      fromNode2.removeEdge(ed2);
      targNode2.removeEdge(ed2);
      this.connectNodes(fromNode1,targNode2,ed1.isDirected(),17);
      this.connectNodes(fromNode2,targNode1,ed2.isDirected(),17);
        
        
    };



    this.addInternalNode = function(nodeNum) {

        if (nodeNum == null) {
            nodeNum = this.internalNodeNextNum;
            ++this.internalNodeNextNum;
        }

        var newNode = this.addNode(this.internalNodePrefix + nodeNum);

        return newNode;
    }
    
    this.addNode = function(lab) {
      //  var node = new DBNode(lab);
        var node  = this.nodeBuilder.createNode(lab);
        this.nodes.push(node);
        return node;
        
    };
    
    
    
    this.findNodes = function(lab) {

        var matchingNodes = [];

        for (var i = 0; i < this.nodes.length; ++i) {

            if (this.nodes[i].label === lab) {
                matchingNodes.push(this.nodes[i]);
            }
        }
   
        return matchingNodes;
        
        
    };

    this.findEdges = function(edge) {
        //return all edges which match start label and end label

        var matches = [];

        this.nodes.forEach(function(grNode) {

            grNode.edges.forEach(function (grEdge) {
                if (edge.directedEdge == grEdge.directedEdge) {
                    if (grEdge.directedEdge) {
                        if ((grEdge.sourceNode.label == edge.sourceNode.label) && (grEdge.targetNode.label == edge.targetNode.label)) {
                            matches.push(grEdge);
                        }
                    }
                    else {

                        if (((grEdge.sourceNode.label == edge.sourceNode.label) && (grEdge.targetNode.label == edge.targetNode.label)) ||
                            (grEdge.sourceNode.label == edge.targetNode.label) && (grEdge.targetNode.label == edge.sourceNode.label)) {
                            matches.push(grEdge);
                        }
                    }
                }

            });
        });

        return matches;



    };



    this.buildNodesFromAdjList = function(adjList,saveProgFlag) {
        //built in default

        if (saveProgFlag == null) {
            saveProgFlag = false;
        }

        var svdThis = this;
        
        var isDirected;

        adjList.forEach(function (adj) {

            isDirected = false;

            var spl = adj.split(' <-> ');
            if (spl.length < 2) {
                spl = adj.split('<->');
            }
            if (spl.length < 2) {
                //isDirected = true;
                spl = adj.split(' -> ');
                if (spl.length < 2) {
                    spl = adj.split('->');
                    if (spl.length < 2) {
                        spl = adj.split(' ');
                    }
                    else {
                        isDirected = true;
                    }
                }
                else {
                    isDirected = true;
                }
            }

            var pref = spl[0];

            var remainder = spl[1].split(':');
            var sufs = remainder[0].split(',');
            var weight = 0;

            if (remainder.length > 1) {

                weight = parseFloat(remainder[1]);

            }
            var prefNode = null;
            var matches = svdThis.findNodes(pref.split(';')[0]);
            if (matches.length > 0) {
                prefNode = matches[0];
            }
            else {
                prefNode = svdThis.addNode(pref);
            }
            sufs.forEach(function (suf) {

                sufNode = null;
                matches = svdThis.findNodes(suf.split(';')[0]);
                var sufNode;
                if (matches.length > 0) {
                    sufNode = matches[0];

                }
                else {
                    sufNode = svdThis.addNode(suf);
                }


                svdThis.connectNodes(prefNode,sufNode,isDirected,weight);

                if (saveProgFlag) {
                    svdThis.saveProgress([prefNode, sufNode]);
                }

            });


        });
        



    }

    this.saveProgress = function(freshNodes,freshEdges) {

        var g = new DBGraph(new DGraphFromNodesBuilder(this.nodes,this.nodeBuilder,this.edgeBuilder,true));



        if (freshNodes) {
            freshNodes.forEach(function(node) {
                var f = g.builder.findNodes(node.label)[0];
                f.freshNode = true;

            });
        }

        if (freshEdges) {
            freshEdges.forEach(function(edge) {
               var f = g.builder.findEdges(edge)[0];
                f.freshEdge = true;
            });
        }


        this.progress.push(g);

    };

    this.nodeBelongsToFreshEdge = function(node) {
      var found = false;
      node.edges.forEach(function(edge) {
          if (edge.freshEdge) {

                  found = true;

          }
      });
      return found;
    };
   // this.buildNodes();
    
}

function DGraphFromNodesBuilder(source,nodeBuilder,edgeBuilder,copyFlag) {
    DGraphBuilder.apply(this, [source, nodeBuilder, edgeBuilder]);

    this.copyFlag = copyFlag;


    this.buildNodes = function () {

        if (this.copyFlag) {


            var newNodes = [];

            var svdThis = this;

            var oldEdgeList = [];
            var newEdgeList = [];

            this.source.forEach(function (node) {
                var newNode = svdThis.addNode(node.label);

            });

            this.source.forEach(function (node, i) {
                node.edges.forEach(function (edge) {
                    if (oldEdgeList.indexOf(edge) > -1) {
                        // svdThis.nodes[i].edges.push(newEdgeList[oldEdgeList.indexOf(edge)]);
                    }
                    else {
                        var newFrom = null;
                        var newTo = null;
                        for (var j = 0; j < svdThis.nodes.length; ++j) {
                            if (svdThis.nodes[j].label == edge.sourceNode.label) {
                                newFrom = svdThis.nodes[j];
                            }
                            if (svdThis.nodes[j].label == edge.targetNode.label) {
                                newTo = svdThis.nodes[j];
                            }

                        }

                        var newEdge = svdThis.connectNodes(newFrom, newTo, edge.directedEdge, edge.edgeWeight);
                        oldEdgeList.push(edge);
                        newEdgeList.push(newEdge);


                    }


                });

                //return source;


            });

            return this.nodes;
        }
        else {
            return source;
        }


    };
}

function DGraphFromSpectrumBuilder(source) {
	
	 DGraphBuilder.apply(this, [source, new DNodeBuilder(), new DEdgeBuilder()]);
	
	 this.buildNodes = function () {
 
            var svdThis = this;
			
			var srcAr = this.source.split(' ');


            srcAr.forEach(function (spec) {
                var newNode = svdThis.addNode(spec);
                
			});
			
			srcAr = srcAr.map(function(el) {
				return parseInt(el);
			});
			
			var amW = Amino.AminoWeights();
			
			var foundWs = [];
			
			for (var i = 0;i < this.nodes.length;++i) {
				for (var j = i+1;j < this.nodes.length;++j) {
					var iW = parseInt(this.nodes[i].label);
					var jW = parseInt(this.nodes[j].label);
					if  (amW.indexOf(jW - iW) > -1) {
						foundWs.push(jW - iW);
						
						var newEdge = this.connectNodes(this.nodes[i],this.nodes[j],true,jW - iW);
					}
				}
			}
			
			
			
			return this.nodes;
	 }
	
}


function DGraphFromSpectralVectorBuilder(source) {
	
	 DGraphBuilder.apply(this, [source, new DNodeBuilder(), new DEdgeBuilder()]);
	
	 this.buildNodes = function () {
 
            var svdThis = this;
			
			var srcAr = this.source.split(' ');

 			
			svdThis.addNode('0:0');
            srcAr.forEach(function (spec,i) {
                var newNode = svdThis.addNode('' + (i+1) + ':' + spec);
                
			});
			
			/*
			srcAr = srcAr.map(function(el) {
				return parseInt(el);
			});
			*/
			
			var amW = Amino.AminoWeights();
			
			var foundWs = [];
			
			for (var i = 0;i < this.nodes.length;++i) {
				for (var j = i+1;j < this.nodes.length;++j) {
					var iW = parseInt(this.nodes[i].label.split(':')[0]);
					var jW = parseInt(this.nodes[j].label.split(':')[0]);
					if  (amW.indexOf(jW - iW) > -1) {
						foundWs.push(jW - iW);
						
						var newEdge = this.connectNodes(this.nodes[i],this.nodes[j],true,jW - iW);
					}
				}
			}
			
			
			
			return this.nodes;
	 }
	
}

function DGraphTreeBuilder(source) {

    DGraphBuilder.apply(this, [source, new DTreeNodeBuilder(), new DTreeEdgeBuilder()]);

    this.saveProgress = function(freshNodes,freshEdges) {
        var  g = new DBTreeGraph(new DGraphTreeFromNodesBuilder(this.nodes,true));

        if (freshNodes) {
            freshNodes.forEach(function(node) {
                var f = g.builder.findNodes(node.label)[0];
                f.freshNode = true;

            });
        }

        if (freshEdges) {
            freshEdges.forEach(function(edge) {
                var f = g.builder.findEdges(edge)[0];
                f.freshEdge = true;
            });
        }

        this.progress.push(g);
    };



}

function DGraphTreeFromNodesBuilder(source,copyFlag) {
    DGraphTreeBuilder.apply(this, [source, new DTreeNodeBuilder(), new DTreeEdgeBuilder()]);

    this.copyFlag = copyFlag;

    this.buildNodes = function () {


        if (this.copyFlag) {


            var newNodes = [];

            var svdThis = this;

            var oldEdgeList = [];
            var newEdgeList = [];

            this.source.forEach(function (node) {
                var newNode = svdThis.addNode(node.label);
                newNode.sequence = node.sequence;

            });

            this.source.forEach(function (node, i) {
                node.edges.forEach(function (edge) {
                    if (oldEdgeList.indexOf(edge) > -1) {
                        // svdThis.nodes[i].edges.push(newEdgeList[oldEdgeList.indexOf(edge)]);
                    }
                    else {
                        var newFrom = null;
                        var newTo = null;
                        for (var j = 0; j < svdThis.nodes.length; ++j) {
                            if (svdThis.nodes[j].label == edge.sourceNode.label) {
                                newFrom = svdThis.nodes[j];
                            }
                            if (svdThis.nodes[j].label == edge.targetNode.label) {
                                newTo = svdThis.nodes[j];
                            }

                        }

                        var newEdge = svdThis.connectNodes(newFrom, newTo, edge.directedEdge, edge.edgeWeight);
                        oldEdgeList.push(edge);
                        newEdgeList.push(newEdge);


                    }


                });

                //return source;


            });

            return this.nodes;
        }
        else {
            return source;
        }


        //return source;
    }
}

function DGraphTreeFromDistMatrixBuilder(source) {
    DGraphTreeBuilder.apply(this, [source, new DTreeNodeBuilder(), new DTreeEdgeBuilder()]);


    this.matrix = [];
    this.leafLabels = [];

    this.isAdditive = false;

    this.checkAdditive = function() {

        for (var i = 0;i < this.matrix.length; ++i) {
            for (var j = 0;j < this.matrix.length;++j) {
                if (i == j) continue;
                for (var m = 0;m < this.matrix.length;++m) {
                    if ((m == i) || (m == j)) continue;
                    for (var n = 0;n < this.matrix.length;++n) {
                        if ((n == i) || (n == j) || (n == m)) continue;
                        var left = this.matrix[i][j] + this.matrix[m][n];
                        var right = Math.max((this.matrix[i][m] + this.matrix[j][n]),(this.matrix[j][m] + this.matrix[i][n]));
                        if (left <= right) {

                        }
                        else {
                            return false;
                        }
                    }
                }
            }
        }
        return true;

    };



    this.parseSource = function () {
        var dataArray = source.split('\n');
        this.leafLabelString = dataArray.shift();
        this.leafLabelString = this.leafLabelString.replace(/\s\s+/g, ' ');
        this.leafLabelString = this.leafLabelString.trim(' ');
        this.leafLabels = this.leafLabelString.split(' ');






        this.matrix = [];

        var svdThis = this;
        dataArray.forEach(function (r) {
            r = r.replace(/\s\s+/g, ' ');
            var rowData = r.split(' ');
            var rowInts = rowData.map(function (rowEntry) {
                return parseInt(rowEntry);
            });
            svdThis.matrix.push(rowInts);
        });

         this.isAdditive = this.checkAdditive();

        if (this.leafLabels.length == 1) {
            var numLeaves = this.matrix.length;
            this.leafLabels = [];
            for (var i = 0; i < numLeaves; ++i) {
                this.leafLabels.push('' + i);

            }
        }

        this.internalNodeNextNum = this.leafLabels.length;
        this.internalNodePrefix = '';


    };

    this.trimMatrix = function(m,leafLabs,ind) {
        // var ind = m.length - 1;

        //var lab = '' + ind; //pick last row
        var lab = leafLabs[ind];

        //var limLen = this.limbLength(lab,m,leafLabs);

        var dTrimmed = [];
        m.forEach(function(row,i) {
            if (i == ind) {

            }
            else {
                var newRow = [];
                row.forEach(function(el,j) {
                    if (j == ind) {

                    }
                    else {
                        newRow.push(el);
                    }

                });
                dTrimmed.push(newRow);
            }
        });
        var newLabs = [];
        leafLabs.forEach(function(el,i) {
            if (i == ind) {

            }
            else {
                newLabs.push(el);
            }
        });

        return [dTrimmed,newLabs];

    };
    
    this.minDistInMatrix = function(m) {
        var minDist = 1000000;
        var lowI = 0;
        var lowJ = 0;
        for (var i = 0;i < m.length;++i) {
            for (var j = 0;j < m.length;++j) {
                if (i == j) {
                    
                }
                else {
                    if (m[i][j] < minDist) {
                        lowI = i;
                        lowJ = j;
                        minDist = m[i][j];
                    }
                }
            }
        }
        return [minDist,lowI,lowJ];
        
    };
    
    this.parseSource();
}

function DGraphUltraTreeFromDistBuilder(source) {
    //ultrametric tree from dist matrix
    DGraphTreeFromDistMatrixBuilder.apply(this,[source]);

    this.clusterDistance = function(clust1,clust2) {
        var totDist = 0;
        var svdThis = this;

        clust1.forEach(function(iVal) {
            clust2.forEach(function(jVal) {
               var pairDist = svdThis.matrix[iVal][jVal];
                totDist += pairDist;

            });
        });
        return totDist / (clust1.length * clust2.length);
    };

    this.buildNodes = function() {

        var svdThis = this;
        this.leafLabels.forEach(function(lab) {
           var newNode = svdThis.addNode(lab);
            newNode.age = 0;
        });

        var done = false;
        var mat = this.matrix.map(function(el) {
            var row = el.map(function(el2) {
                return el2;
            });
            return row;
        });
        var clusters = {};
        var svdThis = this;
        mat.forEach(function(el,i) {
           clusters[svdThis.leafLabels[i].split(';')[0]] = [i];
        });
        var leafLabs = this.leafLabels.map(function(el) {
           return el;
        });


        while (!done) {
            var minD = this.minDistInMatrix(mat);
            var lab1 = leafLabs[minD[1]].split(';')[0];
            var lab2 = leafLabs[minD[2]].split(';')[0];
            var newNode = this.addInternalNode();
            var node1 = this.findNodes(lab1)[0];
            var node2 = this.findNodes(lab2)[0];
            var av = minD[0]/2.0;

            newNode.age = av;

            this.connectNodes(newNode,node1,true,av - node1.age);
            this.connectNodes(newNode,node2,true,av - node2.age);

            this.saveProgress([newNode,node1,node2]);

            //mat = this.trimMatrix(mat,minD[1],minD[2],av);
            clusters[newNode.label] = clusters[lab1].concat(clusters[lab2]);
            //var tmp = this.clusterDistance(clusters['i'],clusters['4']);
            delete clusters[lab1];
            delete clusters[lab2];
            var newMat = [];

            var newFigures = [];

            mat.forEach(function(row,i) {
                if ((i == minD[1] || (i == minD[2]))) {

                }
                else {
                    var newRow = [];
                    row.forEach(function (col, j) {
                        if ((j == minD[1] || (j == minD[2]))) {

                        }
                        else newRow.push(col);

                    });
                    var tmp = svdThis.clusterDistance(clusters[leafLabs[i].split(';')[0]],clusters[newNode.label]);
                    newRow.push(tmp);
                    newFigures.push(tmp);
                    newMat.push(newRow);
                }

            });
            newFigures.push(0);
            newMat.push(newFigures);
            mat = newMat;

            leafLabs = leafLabs.filter(function(el,i) {
                if ((i == minD[1]) || (i == minD[2])) {
                   return false;

                }
                return true;
            });
            leafLabs.push(newNode.label);

            if (leafLabs.length <= 1) {
                done = true;
            }


        }

        return this.nodes;
        
    }
    
}

 function DGraphTreeFromDistBuilder(source) {
     //tree from additive dist matrix using additive phylogeny
     DGraphTreeFromDistMatrixBuilder.apply(this, [source]);


     // this.matrix = [];
     // this.leafLabels = [];

     // this.buildingNodes = []; // for recursive building

     /*
      this.parseSource = function() {
      var dataArray = source.split('\n');
      this.leafLabelString = dataArray.shift();
      this.leafLabelString = this.leafLabelString.replace( /\s\s+/g, ' ' );
      this.leafLabelString = this.leafLabelString.trim(' ');
      this.leafLabels = this.leafLabelString.split(' ');

      this.matrix = [];

      var svdThis = this;
      dataArray.forEach(function(r) {
      r = r.replace( /\s\s+/g, ' ' );
      var rowData = r.split(' ');
      var rowInts = rowData.map(function(rowEntry) {
      return parseInt(rowEntry);
      });
      svdThis.matrix.push(rowInts);
      });

      if (this.leafLabels.length == 1) {
      var numLeaves = this.matrix.length;
      this.leafLabels = [];
      for (var i = 0 ; i < numLeaves; ++i) {
      this.leafLabels.push('' + i);
      }
      }

      this.internalNodeNextNum = this.leafLabels.length;
      this.internalNodePrefix = '';





      };
      */


     this.findAttachmentPointAndAttach = function (mat, leafLabs, ind, limLen) {
         //look for ik such that ij = jk

         var j = ind;

         var found = false;

         for (var i = 0; i < mat.length; ++i) {
             if (i == ind) continue;
             for (var k = 0; k < mat.length; ++k) {
                 if ((k == i) || (k == ind)) continue;

                 if (mat[j][i] + mat[j][k] == mat[i][k]) {
                     found = true;
                     break;

                 }

             }
             if (found) break;
         }

         if (!found) {
             //error

         }
         else {

             var treeSoFar = new DBTreeGraph(new DGraphTreeFromNodesBuilder(this.nodes));
             var nodeI = treeSoFar.getNodeFromLabel(leafLabs[i].split(';')[0]);
             var nodeK = treeSoFar.getNodeFromLabel(leafLabs[k].split(';')[0]);
             var path = treeSoFar.findPathBetweenNodes(nodeI, nodeK, nodeI, null);

             var travelled = 0;
             var attachEdge = null;
             var attachPointOnEdge = 0;
             var attachSource = null;
             var attachTarget = null;
             var newNodeRequired = false;
             var attachNode = null;

             var source = nodeI;

             for (var n = 0; n < path.edgeList.length; ++n) {
                 var targ = path.edgeList[n].getTargetNode(source);

                 if (mat[j][i] < (travelled + path.edgeList[n].edgeWeight)) {
                     attachEdge = path.edgeList[n];
                     attachPointOnEdge = mat[j][i] - travelled;
                     newNodeRequired = true;
                     attachSource = source;
                     attachTarget = targ;
                     break;
                 }
                 else if (mat[j][i] == (travelled + path.edgeList[n].edgeWeight)) {

                     //attachment point already exists
                     //attachNode =   targ;    may need to reinstate this?

                     attachEdge = path.edgeList[n];
                     attachPointOnEdge = mat[j][i] - travelled;
                     newNodeRequired = true;
                     attachSource = source;
                     attachTarget = targ;

                     break;
                 }
                 travelled += path.edgeList[n].edgeWeight;

                 source = targ;

             }

             if (newNodeRequired) {
                 var dist1 = attachPointOnEdge;
                 var dist2 = attachEdge.edgeWeight - attachPointOnEdge;
                 attachNode = this.addNodeBetween(attachSource, attachTarget, dist1, dist2);

             }

             var newLeaf = this.addNode(leafLabs[ind]);
             this.connectNodes(attachNode, newLeaf, false, limLen);

             this.saveProgress([attachNode,newLeaf]);

         }

     };

     this.additivePhylogeny = function (mat, leafLabs) {

         if (mat.length <= 2) {
             var n1 = this.addNode(leafLabs[0]);
             var n2 = this.addNode(leafLabs[1]);
             this.connectNodes(n1, n2, false, mat[0][1]);
             this.saveProgress([n1,n2]);

             return [leafLabs[0] + leafLabs[1]];
         }

         var ind = mat.length - 1;
         var ret = this.trimMatrix(mat, leafLabs, ind);
         //ret = this.trimMatrix(ret[1],ret[2]);
         //ret = this.trimMatrix(ret[1],ret[2]);
         var addRet = this.additivePhylogeny(ret[1], ret[2]);
         this.findAttachmentPointAndAttach(ret[0], leafLabs, ind, ret[3]);




         addRet.push(leafLabs[ind] + '-' + ret[3]);
         return addRet;
     };

     this.trimMatrix = function (m, leafLabs, ind) {
         // var ind = m.length - 1;

         //var lab = '' + ind; //pick last row
         var lab = leafLabs[ind];

         var limLen = this.limbLength(lab, m, leafLabs);
         var dBald = m.map(function (row, i) {
             var newRow = row.map(function (el, j) {
                 if (i == j) {
                     return el;
                 }
                 if ((i == ind) || j == ind) {
                     return el - limLen;
                 }
                 return el;

             });
             return newRow;
         });
         var dTrimmed = [];
         dBald.forEach(function (row, i) {
             if (i == ind) {

             }
             else {
                 var newRow = [];
                 row.forEach(function (el, j) {
                     if (j == ind) {

                     }
                     else {
                         newRow.push(el);
                     }

                 });
                 dTrimmed.push(newRow);
             }
         });
         var newLabs = [];
         leafLabs.forEach(function (el, i) {
             if (i == ind) {

             }
             else {
                 newLabs.push(el);
             }
         });

         return [dBald, dTrimmed, newLabs, limLen];

     };

     this.buildNodes = function () {
         //overriding parent
         if (this.isAdditive) {
             var ret = this.additivePhylogeny(this.matrix, this.leafLabels);
         }
         return this.nodes;
     };

     this.limbLength = function (limbLabel, mat, leafLabs) {
         if (!mat) {
             mat = this.matrix;
         }
         if (!leafLabs) {
             leafLabs = this.leafLabels;
         }
         var j = leafLabs.indexOf(limbLabel);
         var min = 99999999;
         var ln = -1;
         for (var i = 0; i < leafLabs.length; ++i) {
             if (i == j) {
                 continue;
             }
             for (var k = i + 1; k < leafLabs.length; ++k) {
                 if (k == j) {
                     continue;
                 }
                 var ij = mat[i][j];
                 var jk = mat[j][k];
                 var ik = mat[i][k];
                 ln = (ij + jk - ik) / 2.0;
                 if (ln < min) {
                     min = ln;
                 }
             }
         }
         return min;

     }

     //this.parseSource();


 }


function DGraphTreeNJFromDistBuilder(source) {
    //tree from additive or non-additive dist matrix using Neighbour Joining algorithm

    DGraphTreeFromDistMatrixBuilder.apply(this,[source]);

    this.nextInternalNode = this.matrix.length;


    this.neighbourJoiningMatrix = function(mat,leafLabs) {
      var n = mat.length;  
        
        var totDists = mat.map(function(el) {
            return el.reduce(function(el1,el2) {
                return el1 + el2;
            });
            
        });
        
        var njMat = mat.map(function(line,i) {
            var newLine = line.map(function(el,j) {
                if (i == j) {
                    return 0;
                }
                return (n - 2) * mat[i][j] - totDists[i] - totDists[j];
            });
            return newLine;
        });
        
        return [njMat,totDists];
    };

    this.addParentToMatrix = function(mat,leafLabs,iLeaf,jLeaf) {
        var newMat = [];
        var newLine;
        mat.forEach(function(row,r) {
            newLine = [];
            row.forEach(function(col,c) {
                newLine.push(col);
            });
            if ((r == iLeaf) || (r == jLeaf)) {
                newLine.push(0); //will be removed anyway in trim
            }
            else {
                newLine.push( (mat[r][iLeaf] + mat[r][jLeaf] - mat[iLeaf][jLeaf]) / 2);
            }

            newMat.push(newLine);

        });
        newLine = [];
        var mInd = newMat[0].length - 1;
        for (var k = 0;k < mat.length;++k) {
            newLine.push(newMat[k][mInd]);
        }
        newLine.push(0);
        newMat.push(newLine);

        var newLabs = leafLabs.map(function(el) {
            return el;
        });
        newLabs.push('' + this.nextInternalNode);
        ++this.nextInternalNode;

        return [newMat,newLabs];
    };

    this.neighbourJoiningPhylogeny = function(mat,leafLabs,origLeafLabs) {
        //recursive

        if (mat.length <= 2) {
            var n1,n2;
            if (origLeafLabs.indexOf(leafLabs[0]) > -1) {
                //leaf
                n1 = this.addNode(leafLabs[0]);
            }
            else {
                n1 = this.addInternalNode(leafLabs[0]);
            }
            if (origLeafLabs.indexOf(leafLabs[1]) > -1) {
                //leaf
                n2 = this.addNode(leafLabs[1]);
            }
            else {
                n2 = this.addInternalNode(leafLabs[1]);
            }

           //var n1 = this.addNode(leafLabs[0]);
           // var n2 = this.addNode(leafLabs[1]);
            this.connectNodes(n1, n2, false, mat[0][1]);
            this.saveProgress([n1,n2]);
            return [leafLabs[0] + leafLabs[1]];
        }

        var njm = this.neighbourJoiningMatrix(mat, leafLabs);

        var minDist = this.minDistInMatrix(njm[0]);
        var i = minDist[1];
        var j = minDist[2];

        var totDists = njm[1];
        var deltaIJ = (totDists[i] - totDists[j]) / (mat.length - 2);
        var limbLenI = (mat[i][j] + deltaIJ) / 2;
        var limbLenJ = (mat[i][j] - deltaIJ) / 2;

        var addedMat = this.addParentToMatrix(mat, leafLabs, i, j);


        var high = Math.max(i, j);
        var low = Math.min(i, j);
        var trimmedMat = this.trimMatrix(addedMat[0], addedMat[1], high);
        trimmedMat = this.trimMatrix(trimmedMat[0], trimmedMat[1], low);


        this.neighbourJoiningPhylogeny(trimmedMat[0], trimmedMat[1],origLeafLabs);

        var parentNode = this.findNodes(addedMat[1][addedMat[1].length - 1])[0];

        var newNode1,newNode2;
        if (origLeafLabs.indexOf(addedMat[1][i]) > -1) {
            // leaf
            newNode1 = this.addNode(addedMat[1][i]);
        }
        else {
            //var newLeaf = this.addNode(addedMat[1][i]);
            newNode1 = this.addInternalNode(addedMat[1][i]);
        }
        this.connectNodes(parentNode, newNode1, false, limbLenI);

        if (origLeafLabs.indexOf(addedMat[1][j]) > -1) {
            //leaf
            newNode2 = this.addNode(addedMat[1][j]);
        }
        else {
            newNode2 = this.addInternalNode(addedMat[1][j]);
        }

        this.connectNodes(parentNode, newNode2, false, limbLenJ);
        this.saveProgress([parentNode,newNode1,newNode2]);
    };


    this.buildNodes = function() {
        //overriding parent
        var origLeafLabels = this.leafLabels.map(function(el) {
            return el;
        });

        var ret = this.neighbourJoiningPhylogeny(this.matrix,this.leafLabels,origLeafLabels);
        return this.nodes;
    };



}



DGraphTreeFromDistBuilder.AlignmentsToDistMatrix = function(alignments) {
    var alArray = alignments.split('\n');
    var alDict = {};
    
    var labs = '';
    var seqs = [];
    alArray.forEach(function(el,i) {
        labs += el.split(' ')[0] + ';' + el.split(' ')[1];
        if  (i < alArray.length - 1) {
            labs += ' ';
        }
    });
    
    alArray.forEach(function(el,i) {
        seqs.push(el.split(' ')[1]);
        
    });


    var mat = [];
    mat.push(labs);
    seqs.forEach(function(rowSeq,i) {
        var matRow = '';
        seqs.forEach(function(colSeq,j) {
            if (i == j) {
                matRow += '0';
                
            }
            else {
                matRow += hamDist(rowSeq,colSeq);
            }
            if (j < seqs.length - 1) {
                matRow += ' ';
            }
        });
        mat.push(matRow);
        
    });


    var matString = '';
    mat.forEach(function(el,i) {
        matString += el;
        if (i < mat.length - 1) {
            matString += '\n';

        }
    });
    
    return matString;
    
};

function DBGraph(builder,comments) {
    this.nodes = builder.buildNodes();
    this.builder = builder;


    if (comments) {
        this.comments = comments;
    }
    else {
        this.comments = '';
    }
    


    this.checkDirected = function() {
        var undirected = false;
        
        for (var i = 0;i < this.nodes.length;++i) {
            var node = this.nodes[i];

            //for (var j = 0;j < node.successors.length;++j) {
            //new regime:
            for (var j = 0;j < node.edges.length;++j) {
                if (node.edges[j].isDirected()) {

                }
                else {
                    undirected = true;
                    break;
                }

                /* old regime
                var suc = node.successors[j].targetNode;
                for (var k = 0;k < suc.successors.length;++k ) {
                    if (suc.successors[k].targetNode == node) {
                        undirected = true;
                        break;
                        
                    }
                }
                if (undirected) {
                    break;
                }
                */
                
            }
            if (undirected) {
                break;
            }
            
        }
        return (!undirected);
    };


    this.numEdges = function() {
        return this.nodes.map(function(el) {
            return el.outDegree();
        }).reduce(function(a,b) {
            return a + b;
        });
    };

    this.numNodes = function() {
        return this.nodes.length;
    };

   
    //end of "Interface" implementation

    this.resetNodesVisited = function() {
       this.nodes.forEach(function(node) {
           node.visited = false;
       });  
    };
    
    this.getNodeFromLabel = function(lab) {

        //get first node which matches label

        for (var i = 0;i < this.nodes.length; ++i) {
            if (this.nodes[i].label == lab) {
                return this.nodes[i];
            }
        }
        return null;
    };

    this.getNodeIndexForNode = function(node) {

        var ind = this.nodes.indexOf(node);
        return ind;

    };
    
    this.getConnected = function(node) {
        if (node.visited) {
            return [];
        }

        node.visited = true;
        var succs = node.getSuccessorNodes();
        var svdThis = this;
        var allConnected = [node];
        succs.forEach(function(n) {
            allConnected = allConnected.concat(svdThis.getConnected(n));
        });
        return allConnected;

        
    };

    this.isBalanced = function() {
        var balanced = true;

        for (var i = 0;i < this.nodes.length; ++i) {
            if (this.nodes[i].inDegree() == this.nodes[i].outDegree()) {
            }
            else {
                balanced = false;
                break;

            }
        }

        return balanced;



    };

    this.edgeList = function() {
        var edges = [];

        
        
        this.nodes.forEach(function(node) {
           node.edges.forEach(function(ed) {
               if (edges.indexOf(ed) > -1) {
                   
               }
               else {
                   edges.push(ed);
               }
               
           })
        });
        
        return edges;

    };
    
    this.toText = function() {
      
        var txt = '';
        
        this.nodes.forEach(function(node) {
            txt += node.label;
            txt +='-';
            var succs = node.getSuccessors();
            succs.forEach(function(suc,i) {
                txt += suc.getTargetNode(node).label;
                if (i == succs.length - 1) {
                }
                else {
                    txt+=',';
                }

                
            });
            txt+= '\n';
        });
        return txt;
        
    };
    
    this.toAdjacencyList = function(numDec,showSeqs,showDirAsUndir,suppressWeights) {

        if (numDec == null) {
            numDec = 3;
        }

        if (showDirAsUndir == null) {
            showDirAsUndir = false;
        }

        if (showSeqs == null) {
            showSeqs = 1;  //1 - label, 2- seq, 3 -both
        }
        if (suppressWeights == null) {
            suppressWeights = false;
        }

        var adjList = '';
        var svdThis = this;

        /*
        this.nodes.forEach(function(node,i) {
            var sucs = node.getSuccessors();
            sucs.forEach(function(edge,j) {
                adjList+= edge.sourceNode.label;
                adjList += '->' + edge.targetNode.label + ':' + edge.edgeWeight;
                if ((i == svdThis.nodes.length - 1) && (j == sucs.length -1)) {

                }
                else {
                    adjList +='\n';
                }
            });
            
        });
        */

        var edgeList = this.edgeList();

        edgeList.forEach(function(ed,i) {
            var sourceLabPart = '';
            var targetLabPart = '';
            if (showSeqs == 1) {
                sourceLabPart = ed.sourceNode.label;
                targetLabPart = ed.targetNode.label;
            }
            else if (showSeqs == 2) {
                sourceLabPart = ed.sourceNode.sequence;
                targetLabPart = ed.targetNode.sequence;
            }
            else {

                sourceLabPart = ed.sourceNode.label;
                if (ed.sourceNode.sequence.length > 0) {
                    sourceLabPart += ';' + ed.sourceNode.sequence;
                }
                targetLabPart = ed.targetNode.label;
                if (ed.targetNode.sequence.length > 0) {
                    targetLabPart += ';' + ed.targetNode.sequence;
                }
            }

            adjList+= sourceLabPart;
            var weightPart = '';

            if (suppressWeights) {

            }
            else {
                weightPart =':' + ed.edgeWeight.toFixed(numDec);
            }

            adjList += '->' + targetLabPart + weightPart;

            if (ed.isDirected() && (!showDirAsUndir)) {
                if(i == edgeList.length - 1 ) {

                }
                else {
                    adjList +='\n';
                }
            }
            else {
                adjList += '\n';
                adjList += targetLabPart;
                adjList += '->' + sourceLabPart +  weightPart;

                if (i == edgeList.length - 1) {

                }
                else {
                    adjList += '\n';
                }
            }
        });
        return adjList;
    };
    
    
    this.copyGraph = function() {
        var adjList = this.toAdjacencyList();
        var  build = new DGraphBuilder(adjList,this.builder.nodeBuilder,this.builder.edgeBuilder);
        var copiedGraph = new DBGraph(this.builder);
        return copiedGraph;
    };

    
   
    
    this.deletePredecessor = function(node,pred) {

        /* old regime
        for (var i = node.predecessors.length-1;i >=0;--i)  {
            if (node.predecessors[i] == pred) {
                node.predecessors.splice(i,1);
            }
        }
        */

        //new regime:

        //var preds = node.getPredecessors();
        for (var i = node.edges.length -1;i >=0;--i) {
            var ed = node.edges[i];
            if (ed.isDirected()) {
                if (ed.sourceNode === pred) {
                    node.edges.splice(i,1);
                }
            }
            else {
                if ((ed.sourceNode == pred) || (ed.targetNode == pred)) {
                    node.edges.splice(i,1);
                }
            }
        }

    };
    
    this.deleteSuccessor = function(node,succ) {

        /*
        for (var i = node.successors.length-1;i >=0;--i) {
            if (node.successors[i].targetNode == succ) {
                node.successors.splice(i,1);
            }
        }
        */
        for (var i = node.edges.length -1;i >=0;--i) {
            var ed = node.edges[i];
            if (ed.isDirected()) {
                if (ed.targetNode === succ) {
                    node.edges.splice(i,1);
                }
            }
            else {
                if ((ed.targetNode == succ) || (ed.sourceNode == succ)) {
                    node.edges.splice(i,i);
                }
            }
        }


    }

    this.deleteNode = function(node) {
        var svdThis = this;
        node.getSuccessors().forEach(function(succ) {
            var targ = succ.getTargetNode(node);
            svdThis.deletePredecessor(targ,node);

        });
        node.successors = null;


        node.getPredecessors().forEach(function(pred) {
            svdThis.deleteSuccessor(pred,node);
        });
        /*
        node.predecessors.forEach(function(pred) {
            svdThis.deleteSuccessor(pred,node);
            
        });
        */

        node.edges = null;
        
        
        for (var i = this.nodes.length-1;i>=0;--i) {
            if (this.nodes[i] == node) {
                this.nodes.splice(i,1);
            }
        }
        



    };
	
	this.getSourceNode = function() {
		//assumes DAG
		for (var i = 0;i < this.nodes.length;++i) {
			if (this.nodes[i].inDegree() == 0) {
				return this.nodes[i];
			}
		}
		return null;
	};
	
	this.getSinkNode = function() {
		//assumes DAG
		for (var i = this.nodes.length - 1;i >= 0;--i) {
			if (this.nodes[i].outDegree() == 0) {
				return this.nodes[i];
			}
		}
		return null;
	};	
	
	this.longestPathNodeWeighted = function(sourceNode,sinkNode) {
                //longest path in Graph based on node weights not edge weights)
				//assumes DAG
				//assumes nodes are in order, ie later nodes are all successors of prev ones
        
		if (sourceNode) {
			
     
		}
		else {
			sourceNode = this.getSourceNode();
		}
		
		if (sinkNode) {
			
		}
		else {
			sinkNode = this.getSinkNode();
		}
		
		var done = false;
		
		var bests = [[0,0]];
		
		var currNode = sourceNode;
		var currNodeNum = parseInt(currNode.label.split(':')[0]);
		
		var sinkNodeNum = parseInt(sinkNode.label.split(':')[0]);
		
		for (var i = currNodeNum+1;i <= sinkNodeNum;++i) {
			var currNode = this.nodes[i];
			var predNodes = currNode.getPredecessors();
			var bestPred = DGraph.infinity * -1;
			var bestNum = -1;
			
			predNodes.forEach(function(pred) {
				var predNum  = parseInt(pred.label.split(':')[0]);
				var predScore = bests[predNum][0];
				if (predScore  > bestPred) {
					bestPred = predScore;
					bestNum = predNum;
					
				}
					
				
			});
			
			var thisNodeWeight = parseInt(currNode.label.split(':')[1]);
			if (bestPred == DGraph.infinity * -1) {
				thisNodeWeight = 0;
				bestNum = -1;
			}
			bests.push([bestPred + thisNodeWeight,bestNum]);
			
			
		}
		
		
		//backtrack
		var currNum = bests.length - 1;
		
		var diffAr = [];
		
		while (currNum > 0) {
			var diff = currNum - bests[currNum][1];
			diffAr.unshift(diff);
			currNum = bests[currNum][1];
		}
		
		var pep = new Peptide(Peptide.AminoArrFromWeights(diffAr));
		
		return [pep.toShortString(''),bests];
	 
	};
	


      
	
	this.allPathsSourceToSink = function() {
		
         //assumes DAG
		 
		 //iterative queue based
		 
		 var done = false;
		 
		 var sourceNode = this.getSourceNode();
		 var sinkNode = this.getSinkNode();
		 
		 var paths = [];
		 
		 sourceNode.getSuccessors().forEach(function(el) {
			var path = [];
			path.push(el);
            paths.push(path);				
		 });
		 
		 var svdThis = this;
		 
		 while (!done) {
			 var newPaths = [];
			 var done = true;
			 
			 paths.forEach(function(path) {
				if (path[path.length-1].targetNode == svdThis.getSinkNode()) {
					newPaths.push(path);
				}
				else {
					done = false;
			    	var succs = path[path.length-1].targetNode.getSuccessors();
					succs.forEach(function(suc) {
						var newPath = path.map(function(p) {
							return p;
						});
						newPath.push(suc);
						newPaths.push(newPath);
					});
				}
				

				
			 });
			 
			 paths = newPaths;
			 
		 }
		 
		 return paths;
		 
		 
	};



    this.isDirected = this.checkDirected();
}

function DBTreeGraph(builder,comments) {
    DBGraph.apply(this,[builder,comments]);


    this.leaves = function() {
        var lvs = [];
        var svdThis = this;
        lvs = this.nodes.filter(function(node) {
            return (node.isLeaf(svdThis.checkDirected()));
        });

        return lvs;
        
    };
    
    this.internalNodes = function() {
        var iNodes = [];
        var svdThis = this;
        iNodes = this.nodes.filter(function(node) {
           return (!node.isLeaf(svdThis.checkDirected())); 
        });
        
        return iNodes;
    };
    
    this.internalEdges = function() {
        var iEdges = [];
        var svdThis = this;
        var iEdges = this.edgeList().filter(function(ed) {
            if ( (ed.sourceNode.isLeaf(svdThis.checkDirected())
                    || (ed.targetNode.isLeaf(svdThis.checkDirected())))) {
                return false;
            }
            else {
                return true;
            }
        });
        
        return iEdges;
        
    };

    this.stripLeaves = function() {
        var lvs = this.leaves();
        for (var i = this.nodes.length-1;i >=0;--i) {
            if (lvs.indexOf(this.nodes[i]) > -1) {
                this.deleteNode(this.nodes[i]);
            }
        }
       

    };
    
    this.findRoot = function() {
        var rootsFound = 0; //should only find 1. If more, not rooted
        var rootNode = null;
        this.nodes.forEach(function(el) {
            if ((el.outDegree() == 2) && (el.inDegree() == 0)) {
                ++rootsFound;
                rootNode = el;
            }
        });
        
        if (rootsFound == 1) {
            return rootNode;
        }
        else {
            return null;
        }
    };
    
    this.isRooted = function() {
        if (this.findRoot() == null) {
            return false;
        }
        else {
            return true;
        }
        
        
    };

    this.centralNode = function() {

        if (this.nodes.length <= 2) {
            return this.nodes[0].label;
        }

        var strippingTree = this.copyGraph();
        var done = false;

        while (!done) {
            strippingTree.stripLeaves();
       
            if (strippingTree.nodes.length <= 2) {
                done = true;
                return strippingTree.nodes[0].label;
            }

        }


    };
    
    this.unRoot = function() {
        
        var root = this.findRoot();
        var rootEdges = root.getSuccessors();
        var child1 = rootEdges[0].getTargetNode(root);
        var child2 = rootEdges[1].getTargetNode(root);
        this.builder.connectNodes(child1,child2,true,rootEdges[0].edgeWeight + rootEdges[1].edgeWeight);

        this.deleteNode(root);

        var svdThis = this;
        this.nodes.forEach(function(node) {
            node.edges.forEach(function(ed) {
                ed.setDirected(false);
                if (ed.realigned) {
                    var src = ed.sourceNode;
                    ed.sourceNode = ed.targetNode;
                    ed.targetNode = src;
                    ed.realigned = false;
                }
            });
        });
        
    };
    
    
    
    this.addRoot = function(edge) {
        
        //add root to undirected, unrooted tree

      var svdThis = this;
        
      if (this.checkDirected()) {
          return null;
      }
        
        if (this.isRooted()) {
            return null;
      }
        
      var rootNode = this.builder.addNodeBetween(edge.sourceNode,edge.targetNode,edge.edgeWeight /2 , edge.edgeWeight / 2);

      var done = false;
      //rootNode.predecessors = [];
      //rootNode.


      var nodesProcessed = [];
      var nodesToProcess = [rootNode];
      while (nodesToProcess.length > 0) {
          var nextNodesToProcess = [];
          svdThis = this;
          for (var i = 0;i < nodesToProcess.length;++i) {
              var node = nodesToProcess[i];
              //nodesToProcess.forEach(function(node) {
              //node.getSuccessors().forEach(function(suc) {
              for (var j = 0; j < node.getSuccessors().length; ++j) {
                  var suc = node.getSuccessors()[j];
                  var targ = suc.getTargetNode(node);
                  if (nodesProcessed.indexOf(targ) > -1) {

                  }
                  else {
                      nextNodesToProcess.push(targ);
                      suc.setDirected(true, targ);
                  }
              }
          }
          nodesProcessed = nodesProcessed.concat(nodesToProcess);
          nodesToProcess = nextNodesToProcess;

      }

    };

    this.swapNeighbours = function(ed) {
        var neighb = ed.neighbouringEdges();
        var fromNode1 = ed.sourceNode;
        var fromNode2 = ed.targetNode;
        this.builder.swapEdges(neighb[0][1],neighb[1][0],fromNode1,fromNode2);
        
        
      
        
    };



    this.findRipeNodes = function() {

        var ripes = [];

        this.nodes.forEach(function(node) {

            var ripe = true;

            var ripeDets = [node];

            var childEdges = node.getSuccessors();
            var children = [];
            childEdges.forEach(function(ed) {
                children.push(ed.getTargetNode(node));
            });

            if (node.parsimonyScore) {
                ripe = false; //already done
            }

            if (children.length == 0) {
                ripe = false;
            }

            children.forEach(function(child) {

                if(child.parsimonyScore) {
                    ripeDets.push(child.parsimonyScore);

                }

                else {

                    ripe = false;

                }

            });

            if (ripe) {

                ripes.push(ripeDets);

            }

        });

        return ripes;

    };
    
    this.largeParsimony = function(useAminos) {
      
        var svdThis = this;

        var bestEdge = null;
        var  switches = [[0,2,1,3],[0,3,1,2],[0,1,2,3]]; //last one reverts to orig
        var bestSwitch = -1;
        var bestScore = DGraph.infinity;
        var bestNeighbNodes = [];

        var prevBestScore = DGraph.infinity;

        var done = false;
        while (!done) {

            prevBestScore = bestScore;

            //var swTree = this.copyGraph();

            this.internalEdges().forEach(function (ed) {

                var neighb = ed.neighbouringEdges();
                var neighbNodes = [neighb[0][0].getTargetNode(ed.sourceNode),
                    neighb[0][1].getTargetNode(ed.sourceNode),
                    neighb[1][0].getTargetNode(ed.targetNode),
                    neighb[1][1].getTargetNode(ed.targetNode)];


                var fromNode1 = ed.sourceNode;
                var fromNode2 = ed.targetNode;

                var score = DGraph.infinity;


                switches.forEach(function (sw, i) {
                    ed.removeAllConnectingEdges();
                    svdThis.builder.connectNodes(ed.sourceNode, neighbNodes[sw[0]], svdThis.checkDirected(), 17);
                    svdThis.builder.connectNodes(ed.sourceNode, neighbNodes[sw[1]], svdThis.checkDirected(), 17);
                    svdThis.builder.connectNodes(ed.targetNode, neighbNodes[sw[2]], svdThis.checkDirected(), 17);
                    svdThis.builder.connectNodes(ed.targetNode, neighbNodes[sw[3]], svdThis.checkDirected(), 17);
                    score = svdThis.smallParsimony(null,useAminos);
                    if (score < bestScore) {
                        bestEdge = ed;
                        bestScore = score;
                        bestSwitch = i;
                        bestNeighbNodes = neighbNodes;
                    }

                });


            });

            if (bestScore < prevBestScore) {
                var thisBestEdge = null;
                var edgeL = this.edgeList();
                for (var i = 0;i < edgeL.length;++i) {
                    if (((edgeL[i].sourceNode.label == bestEdge.sourceNode.label) && (edgeL[i].targetNode.label == bestEdge.targetNode.label))
                        || ((edgeL[i].sourceNode.label == bestEdge.targetNode.label) && (edgeL[i].targetNode.label == bestEdge.sourceNode.label)))
                    {
                       thisBestEdge = edgeL[i];
                    }
                }

                var sw = switches[bestSwitch];
                //var neighb = thisBestEdge.neighbouringEdges();
                //var neighbNodes = [neighb[0][0].getTargetNode(thisBestEdge.sourceNode),
                //    neighb[0][1].getTargetNode(thisBestEdge.sourceNode),
                //    neighb[1][0].getTargetNode(thisBestEdge.targetNode),
                //    neighb[1][1].getTargetNode(thisBestEdge.targetNode)];

                var neighbNodes = bestNeighbNodes;
                thisBestEdge.removeAllConnectingEdges();
                svdThis.builder.connectNodes(thisBestEdge.sourceNode, neighbNodes[sw[0]], svdThis.checkDirected(), 17);
                svdThis.builder.connectNodes(thisBestEdge.sourceNode, neighbNodes[sw[1]], svdThis.checkDirected(), 17);
                svdThis.builder.connectNodes(thisBestEdge.targetNode, neighbNodes[sw[2]], svdThis.checkDirected(), 17);
                svdThis.builder.connectNodes(thisBestEdge.targetNode, neighbNodes[sw[3]], svdThis.checkDirected(), 17);
                svdThis.smallParsimony(null,useAminos); //added 18/3/17
                svdThis.builder.saveProgress([thisBestEdge.sourceNode,thisBestEdge.targetNode,neighbNodes[sw[0]],neighbNodes[sw[1]],neighbNodes[sw[2]],neighbNodes[sw[3]]],[thisBestEdge]);

            }
            else {
                done = true;
            }

        }
        return bestScore;

    };

    this.prepareSmallParsimony = function(rooted,useAminos) {

        var alphabet = c_Bases;
        if (useAminos) {
            alphabet = c_Aminos;
        }

        rooted.nodes.forEach(function(node) {
            node.parsimonyScore = null;
        });

        var lvs = rooted.leaves();
        lvs.forEach(function(leaf) {
            leaf.parsimonyScore = [];
            if (leaf.sequence.length == 0) {
                leaf.sequence = leaf.label;
            }
            seqLength = leaf.sequence.length;
            for (var i = 0;i < seqLength;++i) {
                var scoreDict = {};
                alphabet.forEach(function (ch) {
                    var score = (ch == leaf.sequence[i]) ? 0 : DGraph.infinity;
                    scoreDict[ch] = score;
                });
                leaf.parsimonyScore.push(scoreDict);

            }


        });

        if (seqLength < 346) {
            console.log(' seq len: ' + seqLength);
        }
        return seqLength;


    };
    
    this.smallParsimony = function(stepNum,useAminos) {
      //find parsimonious sequences for internal nodes given sequences for leaves (current species)
      //if useAminos is set, assumes sequences are amino acides, not bases

        var seqLength = 0;

       


        var tempRoot = false;
        var rooted = this.copyGraph();

        if (!this.isRooted()) {
            tempRoot = true; //temporarily make tree rooted
            rooted.addRoot(rooted.nodes[0].edges[0]);
        }



        
        //do parsimony here

        var svdThis = this;


        seqLength = this.prepareSmallParsimony(rooted,useAminos);

        var ripes  = rooted.findRipeNodes();




        while (ripes.length > 0) {
            ripes.forEach(function(ripeDets) {
                var ripeNode = ripeDets[0];
                console.log('calcing pars score: ' + ripeNode.label);
                ripeNode.parsimonyScore = ripeNode.calcParsimonyScore(ripeDets[1],ripeDets[2],useAminos);

            });

            ripes = rooted.findRipeNodes();

        }


        var nodeList = rooted.findRoot().nodeOrderToLeaves();

        var seq = '';

        var totParsimonyScore = 0;

        var edgesProcessed = [];

        nodeList.forEach(function(node,n) {
            var seq = '';
           for (var i = 0;i <seqLength;++i )  {
               if (n == 0) {
                   charToUse = '';
               }
               else {
                   charToUse = node.getPredecessors()[0].sequence[i];
               }

               if (node.parsimonyScore[i] == null) {
                   console.log('i: ' + i + ' node par null. char: ' + charToUse + ' node: ' + node.label + ' nodepars: ' + node.parsimonyScore.length + ' seqlen: '  +seqLength);
               }
               var best = node.bestParsimonyBasedOnChar(charToUse,node.parsimonyScore[i]);
               seq += best[1];
               if (n == 0) {
                   totParsimonyScore += best[0];
               }
           }
            node.sequence = seq;
        });



        if (tempRoot) {
            rooted.unRoot();
        }

        var svdThis = this;
        this.nodes.forEach(function(node) {
           var rootedNode =  rooted.getNodeFromLabel(node.label);
            node.sequence = rootedNode.sequence;
            node.parsimonyScore = rootedNode.parsimonyScore;
        });

        //calc weights

        this.nodes.forEach(function(node) {
            node.edges.forEach(function(edge) {
                edge.edgeWeight = hamDist(edge.sourceNode.sequence,edge.targetNode.sequence);
            });

        });

        this.nodes.forEach(function(node) {



            var i = 0;

            var ch = '%';
            var refCh = '';

            var ref = 'IPKTQKAMVFYKNGGPLKYEDIPVPKPKPSEILINVRYSGVCHTDLHAWKGDWPLPTKLPLVGGHEGAGVVVACGSEVKNFKVGDYAGIKWLNGSCMGCEYCMQGAEPNCPKADLSGYTHDGSFQQYATADAVQAAHIPQGTDLAAAAPILCAGVTVYKALKTADLRPGQWVAISGAGGGLGSLAVQYAKAMGLRVVGIDGGSEKKELATKLGAEEFIDFTQVSDVVKXMQNVTNGGPHGVINVSVSPRAMSQSVEYVRTLGKVVLVGLPADAVVQTKVFDHVIKSIQIRGSYVGNREDTAEALDFFERGLVHSPIKVVGLSDLPKVFSLMEKGKIAGRYVLDTSK';

            while (ch != '') {

                ch = node.sequence.substring(i, i+1);
                refCh = ref.substring(i,i+1);
                if (node.label.substring(0,8) == 'B.bruxel') {

                    console.log('ch: ' + ch + ' ref: ' + refCh);
                }

                ++i;



            }
            var manNum = i+1;
           console.log('sm pars node: ' + node.label + ' seq: ' + node.sequence.length + 'man: ' +  manNum + ' ' + node.sequence.substring(0,10));
        });
        
        return totParsimonyScore;
    };

    this.copyGraph = function() {
        var adjList = this.toAdjacencyList(3,3);
        var  build = new DGraphBuilder(adjList,this.builder.nodeBuilder,this.builder.edgeBuilder);
        var copiedGraph = new DBTreeGraph(build);
        return copiedGraph;
    };

    this.findPathBetweenNodes = function(node1,node2,startNode,prevNode) {
        if (node1 == node2) {
            return {'pathLen':0,'edgeList':[]}; //lenSoFar; //found it
        }
        else if (node1.isLeaf() && (node1 != startNode)) {
            //dead end
            return {'pathLen':-1,'edgeList':[]};
        }
        else {
            var succs = node1.getSuccessors();
            for (var i = 0;i < succs.length;++i) {
                var suc = succs[i];
                if (suc.getTargetNode(node1) == prevNode) {

                }
                else {
                    var subPath = this.findPathBetweenNodes(suc.getTargetNode(node1),node2,startNode,node1);
                    if (subPath.pathLen  >= 0) {
                        subPath.edgeList.unshift(suc);
                        subPath.pathLen += suc.edgeWeight;
                        return subPath;
                        //return subPath + suc.edgeWeight;
                        
                    }
                    
                }
                
            }
            return {'pathLen':-1,'edgeList':[]}; //did't find
            

        }

    }
    
    this.distanceMatrixFromTree = function() {
        
        //assumes all edges are already weighted
        
        var leaves = this.leaves();
        var rows = [];
        for (var i = 0;i < leaves.length;++i) {
            var row = [];
            for (var j = 0;j < leaves.length;++j) {
                if (i == j) {
                    row.push(0);
                }
                else {
                    row.push(this.findPathBetweenNodes(leaves[i],leaves[j],leaves[i],null).pathLen);
                }
                
            }
            rows.push(row);
        }
        return rows;
    }
}


function DBasicGraph(source,sourceType) {

    if (source) {
        this.source = source.split('\n');
    }
    else {
        this.source = null;
    }
    this.sourceType = sourceType;

    this.nodes = {};

    this.numNodes = 0;


    this.stringToChrom = function (sourceStr) {
        sourceStr = sourceStr.replace(/\(/g, '');
        var chromesAr = sourceStr.split(')');

        chromesAr.pop();

        chromesAr = chromesAr.map(function (el) {
            return el.split(' ');
        });

        return chromesAr;


    };

    this.findNonVisitedNeighbour = function (node) {
        for (var i = 0; i < node.successors.length; ++i) {
            if (node.successors[i].targetNode.visited) {
            }
            else {
                return node.successors[i].targetNode;
            }
        }

        return null;


    };

  

    this.traverseCycle = function (st) {

        var finished = false;

       // var numNodesInCycle = 1;
        var nodesInCycle = [st];

        var node = st;

        node.visited = true;

        while (!(finished)) {

            var nextUnvisited = this.findNonVisitedNeighbour(node);
            if (nextUnvisited) {
                node = nextUnvisited;
                node.visited = true;
               // ++numNodesInCycle;
                nodesInCycle.push(node);
            }
            else {
                finished = true;
            }

        }

        return nodesInCycle;

    };
    
    this.graphToGenome = function(edgeType) {

        this.resetNodesVisited();

        var genStr = '';



        var ind = 0;
        if ((edgeType) && (edgeType == 'blue') ) {
            ind = 1;
        }

        for (var i = 1;i < this.numNodes;++i) {

            if (this.nodes['' + i].visited) {
                continue;
            }

            genStr += '(';
            var node = this.nodes['' + i];
            var startNode = node;
            var startNodeLabNum = Math.ceil(parseInt(startNode.label) / 2);
            node.visited = true;
            var nodeSign = (parseInt(node.label) % 2 == 0) ? '-' : '+';
            var nodeOtherEnd = nodeSign == '-' ? parseInt(node.label) - 1 : parseInt(node.label) + 1;
            this.nodes['' + nodeOtherEnd].visited = true;

            var nodeLabNum = Math.ceil(parseInt(node.label) / 2);
            genStr += nodeSign + nodeLabNum;
            node = this.nodes['' + nodeOtherEnd].successors[ind].targetNode;
            var nodeLabNum = Math.ceil(parseInt(node.label) / 2);
            while (nodeLabNum != startNodeLabNum) {
                var nodeSign = (parseInt(node.label) % 2 == 0) ? '-' : '+';
                genStr += ' ' + nodeSign + nodeLabNum;
                var nodeOtherEnd = nodeSign == '-' ? parseInt(node.label) - 1 : parseInt(node.label) + 1;
                node.visited = true;
                this.nodes['' + nodeOtherEnd].visited = true;

                node = this.nodes['' + nodeOtherEnd].successors[ind].targetNode;
                var nodeLabNum = Math.ceil(parseInt(node.label) / 2);
            }
            genStr+=')';

        }
        return genStr;
    };

    this.findEdgeInNonTrivialCycle = function () {
        
        var genStages = [];
        
        var edgeInd = 1; //blue edge


        var genStr = this.graphToGenome('red');
        genStages.push(genStr);

        for (var i = 1; i <= this.numNodes; ++i) {
            var node = this.nodes['' + i];
            var blueEdge = node.successors[1];
            var redEdge = node.successors[0];
            if (redEdge.targetNode == blueEdge.targetNode) {
                //trivial
            }
            else {
                var otherNode = blueEdge.targetNode;
                var otherRedEdge = otherNode.successors[0];
                var pair1 = [node,otherNode]; // makes trivial edge
                var pair2 = [redEdge.targetNode,otherRedEdge.targetNode];
                //create new red edges:
                var newRedEdge1 = new DGEdge(node,otherNode,'red');
                var newRedEdge2 = new DGEdge(redEdge.targetNode,otherRedEdge.targetNode);
                //replace existing red edges:
                this.connectNodesUndirected(node,otherNode,'red',true);
                this.connectNodesUndirected(redEdge.targetNode,otherRedEdge.targetNode,'red',true);
                genStr = this.graphToGenome('red');
                genStages.push(genStr);

            }

        }
        return genStages;

    };

    this.resetNodesVisited = function () {
        var keys = Object.keys(this.nodes);

        var svdThis = this;

        keys.forEach(function (el) {
            svdThis.nodes[el].visited = false;
            svdThis.nodes[el].cycleNum = 0;
        });
    };

            this.twoBreakDistance = function() {
                return this.numBlocks() - this.numCycles();
            };

            this.numBlocks = function() {
                return this.numNodes / 2;
            };

            this.numCycles = function(minLength) {

                if (minLength == null) { //optional
                   minLength = 0;  // min number of nodes in a cycle, otherwise doesn't count as a cycle
                }

                var visitedNodes = [];

                this.resetNodesVisited();

                var cycles = 0;

              //  for (var i = 1;i <= this.numNodes;++i) {

               // var keys = Object.keys(this.nodes);
               // keys.sort();
                
                for (var node in this.nodes) {
                //keys.forEach(function(node) {
                    if (this.nodes[node].visited) {

                    }
                    else {
                        var nodesInCycle = this.traverseCycle(this.nodes[node]);
                        if (nodesInCycle.length > minLength) {
                            ++cycles;
                            nodesInCycle.forEach(function(el) {
                                el.cycleNum = cycles;
                            });
                        }
                        else {
                            nodesInCycle.forEach(function(el) {
                                el.cycleNum = 0;
                            });
                            
                        }

                    }
                }


                return cycles;




            };


            this.connectNodesUndirected = function(node1,node2,edgeType,replaceFlag) {
                var ed1 = new DGEdge(node1,node2,edgeType);
                var ed2 = new DGEdge(node2,node1,edgeType);

                var ind = 0;
                if ((edgeType) && (edgeType == 'blue')) {
                    ind = 1;
                }

                if (replaceFlag) {
                   node1.successors[ind] = ed1;
                   node1.predecessors[ind] = ed2;
                   node2.successors[ind] = ed2;
                   node2.predecessors[ind] = ed1;
                }
                else {
                    node1.successors.push(ed1);
                    node1.predecessors.push(ed2);
                    node2.successors.push(ed2);
                    node2.predecessors.push(ed1);
                }

            };

            this.initFromBreakpoint = function() {
                var pSource = this.source[0];
                var qSource = this.source[1];

                //construct nodes using p
                var p = this.stringToChrom(pSource);
                var q = this.stringToChrom(qSource);

                var prev,first;

                var svdThis = this;

                //      p[0].forEach(function(el,i) { //assumes only 1 p chromosome
                p.forEach(function(pChr,a) { //assumes could be multiple q chromosomes
                    pChr.forEach(function (el, i) {

                        var sign = el[0];
                        var num = parseInt(el.substring(1));
                        var stNum, endNum;
                        if (sign == '+') {
                            stNum = num * 2 - 1;
                            endNum = num * 2;
                        }
                        else {
                            stNum = num * 2;
                            endNum = num * 2 - 1;
                        }

                        var stNode = new DGNode('' + stNum);
                        var endNode = new DGNode('' + endNum);


                        svdThis.nodes['' + stNum] = stNode;
                        svdThis.nodes['' + endNum] = endNode;

                        if (i == 0) {
                            first = stNode;

                        }
                        else {
                            svdThis.connectNodesUndirected(prev, stNode,'red');
                        }
                        if (i == pChr.length - 1) {
                            svdThis.connectNodesUndirected(endNode, first,'red');
                        }


                        prev = endNode;
                    });

                });

                this.numNodes = Object.keys(this.nodes).length;

                q.forEach(function(qChr,b) { //assumes could be multiple q chromosomes
                    qChr.forEach(function(el,i) {
                        var sign = el[0];
                        var num = parseInt(el.substring(1));
                        var stNum,endNum;
                        if (sign == '+') {
                            stNum = num * 2 - 1;
                            endNum = num * 2;
                        }
                        else {
                            stNum = num * 2;
                            endNum = num * 2 - 1;
                        }

                        var stNode = svdThis.nodes['' + stNum];
                        var endNode = svdThis.nodes['' + endNum];

                        if (i == 0) {
                            first = stNode;

                        }
                        else {
                            svdThis.connectNodesUndirected(prev,stNode,'blue');
                        }
                        if (i == qChr.length - 1)  {
                            svdThis.connectNodesUndirected(endNode,first,'blue');
                        }


                        prev = endNode;


                    });
                });

            };

            this.initGraph = function() {

                switch (this.sourceType) {
                    case DGraph.fromBreakpoint:
                        this.initFromBreakpoint();

                        break;


                    default:
                        break;
                }


            };





    }



    /*
    //Hamilton path routines

    function HamNode(dna,k) {
        this.dna = dna;
        this.k = k;

        this.successors = [];
        this.predecessors = [];

        this.suffix = function() {
            return this.dna.substring(1);
        };
        this.prefix = function() {
            return this.dna.substring(0,k-1);
        };

        this.visited = false;

    }

    function hamPath(nodes) {
        var first;
        var path = []; //path containing kmers

        var nodePath = []; //path containing nodes themselves


        var definiteFirstFound = false;
        var definiteLastFound = false;

        nodes.forEach(function(node) {
            node.visited = false;
        });


        for (var i = 0;i < nodes.length; ++i) {
            if (nodes[i].predecessors.length == 0) {
                first = nodes[i];
                definiteFirstFound = true;
                break;
            }
        }
        if (!first) {
            var r = getRandomInt(0,nodes.length - 1);
            first = nodes[r];
            //return [[],''];
        }

        var done = false;
        path.push(first.dna);
        nodePath.push(first);

        first.visited = true;
        var curr = first;

        while (!done) {
            var potentialNext = curr.successors;
            var nonVisited = [];
            curr.successors.forEach(function(el,i) {
                if (!el.visited) {
                    nonVisited.push(i);
                }

            });
            if (nonVisited.length == 0) {
                done = true;
            }
            else {
                var r = getRandomInt(0,nonVisited.length-1);
                curr  = curr.successors[nonVisited[r]];
                curr.visited = true;
                path.push(curr.dna);
                nodePath.push(curr);
            }

        }

        var reconstructed = '';
        path.forEach(function(el,i) {
            if (i == 0) {
                reconstructed+=el;
            }
            else {
                reconstructed+=el.substring(el.length-1);
            }


        });

        return [path,reconstructed,nodePath];
    }

    function createHamNodes(dna,reads,k) {
        var kmers;
        if (dna) {
            kmers = kmerComposition(dna, k);
        }
        else {
            kmers = reads;
        }

        var nodes = kmers.map(function(el) {
            var node = new HamNode(el,k);
            return node;
        });

        nodes.forEach(function(node,i) {
            var succ = [];
            var pred = [];
            nodes.forEach(function(searchNode,j) {
                if (i == j) {

                }
                else
                    if (node.suffix() ===  searchNode.prefix()) {
                        succ.push(searchNode);
                    }
                    if (node.prefix() ===  searchNode.suffix()) {
                       pred.push(searchNode);
                    }


            });
            node.successors = succ;
            node.predecessors = pred;

        });


        return nodes;


       // var nodes = kmers.map(function(el) {


        //});

    }
    */

function stringFromGenomePath(path) {
    // not very useful - assumes as input array of kmers which all overlap previous by 1

    var str = '';

    path.forEach(function(el,i) {

        if (i == 0) {
            str = el;
        }
        else {
            str+=el[el.length-1];
        }

    });

    return str;

}



function Amino(amino) {

    this.short = amino;

    this.init = function() {

    };

    this.medText = function() {
        if (this.short in Amino.transTable) {
            return Amino.transTable[this.short][0];
        }
        else {
            return this.short;
        }
    };

    this.longText = function() {
        if (this.short in Amino.transTable)  {
            return Amino.transTable[this.short][1];
        }
        else {
            return this.short;
        }
    };
    this.getIntegerWeight = function() {
        if (this.short in Amino.transTable) {
            return Amino.transTable[this.short][2];
        }
        else {
            //assumes the amino short title contains the weight
            if (isNaN(this.short)) {
                return 0;
            }
            else {
                return parseInt(this.short);
            }
            
        }

    };

    this.init();
}


Amino.transTable =  {
    'K' : ['Lys','Lycine',128],
    'V' : ['Val','Valine',99],
    'Q' : ['Gln','Glutamine',128],
    'N' : ['Asn','Asparagine',114],
    'F' : ['Phe','Phenylalanine',147],
    'L' : ['Leu','Leucine',113],
    'I' : ['Ile','Isoleucine',113],
    'M' : ['Met','Methionine',131],
    'S' : ['Ser','Serine',87],
    'P' : ['Pro','Proline',97],
    'T' : ['Thr','Threonine',101],
    'A' : ['Ala','Alanine',71],
    'Y' : ['Tyr','Tyrosine',163],
    '*' : ['STP','STP',0],
    'H' : ['His','Histidine',137],
    'D' : ['Asp','Aspartic Acid',115],
    'E' : ['Glu','Glutamic Acid',129],
    'C' : ['Cys','Cysteine',103],
    'W' : ['Trp','Tryptophan',186],
    'R' : ['Arg','Arginine',156],
    'G' : ['Gly','Glysine',57]

};

// Toy testing:

/*
Amino.transTable =  {
    'X' : ['Lys','Lycine',4],
    'Z' : ['Val','Valine',5]
}

*/

Amino.AminoWeights = function() {
    var aminoWeights = [];

    for (var amEntry in Amino.transTable) {
        if (Amino.transTable.hasOwnProperty(amEntry)) {
            if (amEntry == '*') {

            }
            else {
                if (aminoWeights.indexOf(Amino.transTable[amEntry][2]) == -1) {
                    aminoWeights.push(Amino.transTable[amEntry][2]);
                }
            }
        }

    }
    aminoWeights.sort(function(a,b) {
        return a - b;

    });
    return aminoWeights;
};


Amino.CodeForWeight = function(w) {
	  
	for (var amEntry in Amino.transTable) {
        if (Amino.transTable.hasOwnProperty(amEntry)) {
           if (Amino.transTable[amEntry][2] == w) {
			   return amEntry;
		   }
            
        }
	}
	return '-';
		
}

Amino.AllAminos = function() {

    var amList = [];

    for (var amEntry in Amino.transTable) {
        if (Amino.transTable.hasOwnProperty(amEntry)) {
            if (amEntry == '*') {

            }
            else {
                amList.push(amEntry);


            }
        }

    }
    amList.sort();

    return amList;

};

Amino.AllUniqueWeightAminos = function(include200) {
    
    //include200 - flag to include all other potential Amino Acids between weights 57 - 200
    
    var amList = Amino.AllAminos();

    var amUniqueList = [];

    var foundWeights = [];
    
    if (include200) {
        for (var w = Amino.LowestWeight;w <= Amino.HighestWeight;++w) {
            if (foundWeights.indexOf(w) > -1) {

            }
            else {
                var foundAm = false;
                for (var amEntry in Amino.transTable) {
                    if (Amino.transTable.hasOwnProperty(amEntry)) {
                        if (Amino.transTable[amEntry][2] == w) {
                            if (foundAm) {

                            }
                            else {
                                amUniqueList.push(amEntry);
                                foundAm = true;
                            }

                        }
                    }
                }
                if (foundAm) {

                }
                else {
                    amUniqueList.push(w.toString());
                }
                foundWeights.push(w);
            }

        }
        return amUniqueList;
    }
    else {

        amList.forEach(function (am) {
            var w = Amino.transTable[am][2];
            if (foundWeights.indexOf(w) > -1) {

            }
            else {
                foundWeights.push(w);
                amUniqueList.push(am);
            }
        });

        amUniqueList.sort();

        return amUniqueList;
    }

};

Amino.LowestWeight = 57;
Amino.HighestWeight = 200;


function Codon(codon) {
    this.codon = codon;

    this.init = function() {

    };

    this.isStartCodon = function() {
        return (Codon.StartCodons.indexOf(this.codon) > -1);
    };

    this.isStopCodon = function() {
        return (Codon.StopCodons.indexOf(this.codon) > -1);
    };

    this.translate = function() {
        if (this.codon in Codon.transTable) {
            return new Amino(Codon.transTable[this.codon]);
        }
        else {
            return null;
        }


    };

    this.init();

}

Codon.len = 3;
Codon.StartCodons = ['AUG'];
Codon.StopCodons = ['UAA','UAG','UGA'];

Codon.transTable = {
    'AAA' : 'K',
    'AAC' : 'N',
    'AAG' : 'K',
    'AAU' : 'N',
    'ACA' : 'T',
    'ACC' : 'T',
    'ACG' : 'T',
    'ACU' : 'T',
    'AGA' : 'R',
    'AGC' : 'S',
    'AGG' : 'R',
    'AGU' : 'S',
    'AUA' : 'I',
    'AUC' : 'I',
    'AUG' : 'M',
    'AUU' : 'I',
    'CAA' : 'Q',
    'CAC' : 'H',
    'CAG' : 'Q',
    'CAU' : 'H',
    'CCA' : 'P',
    'CCC' : 'P',
    'CCG' : 'P',
    'CCU' : 'P',
    'CGA' : 'R',
    'CGC' : 'R',
    'CGG' : 'R',
    'CGU' : 'R',
    'CUA' : 'L',
    'CUC' : 'L',
    'CUG' : 'L',
    'CUU' : 'L',
    'GAA' : 'E',
    'GAC' : 'D',
    'GAG' : 'E',
    'GAU' : 'D',
    'GCA' : 'A',
    'GCC' : 'A',
    'GCG' : 'A',
    'GCU' : 'A',
    'GGA' : 'G',
    'GGC' : 'G',
    'GGG' : 'G',
    'GGU' : 'G',
    'GUA' : 'V',
    'GUC' : 'V',
    'GUG' : 'V',
    'GUU' : 'V',
    'UAA' : '*',
    'UAC' : 'Y',
    'UAG' : '*',
    'UAU' : 'Y',
    'UCA' : 'S',
    'UCC' : 'S',
    'UCG' : 'S',
    'UCU' : 'S',
    'UGA' : '*',
    'UGC' : 'C',
    'UGG' : 'W',
    'UGU' : 'C',
    'UUA' : 'L',
    'UUC' : 'F',
    'UUG' : 'L',
    'UUU' : 'F'

};

function Peptide(aminoArr,startPosInRNA) {

    this.peptide = aminoArr;

    this.startPosInRNA = startPosInRNA ? startPosInRNA : 0;

    this.init = function () {

    };

    this.toString = function (detailLevel, separator,reverse,cycleStart) {

        //cycleStart: offset of peptide string to make a different cyclical vers
        cycleStart = cycleStart ? cycleStart : 0;

        var sep;
        if (separator == null) {
            sep = '-'
        }
        else {
            sep = separator;
        }

        if (!detailLevel) {
            detailLevel = Peptide.Short;
        }

        var str = '';
        
        var pepStrAr = [];

        var pepAr = this.peptide.slice(cycleStart).concat(this.peptide.slice(0,cycleStart));
        
        //this.peptide.forEach(function (am, i) {
        pepAr.forEach(function (am,i) {
            if (i > 0) {
                str += sep;
                pepStrAr.push(sep);
            }

            if (am) {
                switch (detailLevel) {
                    case Peptide.Short:
                        str += am.short;
                        pepStrAr.push(am.short);
                        break;

                    case Peptide.Medium:
                        str += am.medText();
                        pepStrAr.push(am.medText());
                        break;
                    case Peptide.Long:
                        str += am.longText();
                        pepStrAr.push(am.longText());
                        break;
                    default:
                        str += am.short;
                        pepStrAr.push(am.short);
                        break;

                }


            }
            else {
                str += 'unk';
                pepStrAr.push('unk');
            }
        });
        
        if (reverse) {
            str = '';
            for (var i = pepStrAr.length - 1;i >=0;--i) {
                str += pepStrAr[i];
            }
       }
        return str;

    };


    this.toMedString = function (separator) {
        return this.toString(Peptide.Medium, separator);
    };

    this.toLongString = function (separator) {
        return this.toString(Peptide.Long, separator);
    };

    this.toShortString = function (separator) {
        return this.toString(Peptide.Short, separator);
    };

    this.toWeightString = function (separator) {

        var sep = '-';

        if (separator) {
            if (separator.length == 0) {
                sep = '';
            }
            else {
                sep = separator;
            }
        }
        var mapped = this.peptide.map(function (am) {
            return am ? am.getIntegerWeight() : 0;
        });

        var str = '';
        mapped.forEach(function (el, i) {
            str += el;
            if (i < mapped.length - 1) {
                str += sep;
            }

        });
        return str;

    };

    this.toStringArray = function () {

        var arr;
        arr = this.peptide.map(function (am) {
            return am.short;
        });

        return arr;

    };


    this.getIntegerWeight = function () {
        var mapped = this.peptide.map(function (am) {
            return am ? am.getIntegerWeight() : 0;
        });


        if ((mapped) && mapped.length > 0) {
            return mapped.reduce(function (a, b) {
                return a + b;
            });
        }
        else {
            return 0;
        }
    };
	
	this.prefixMasses = function() {
		
		var prefMasses = [];
		
		var totMass = 0;
		
		
		
		this.peptide.forEach(function(am) {
			var w = am.getIntegerWeight();
			totMass += w;
			prefMasses.push(totMass);
		});
		
		return prefMasses;
		
	};
	
	
	this.toPeptideVector = function() {
		//prefix masses
		
		var prefMasses = this.prefixMasses();
		
		var vect = [];
		
		for (var i = 1;i <=		 prefMasses[prefMasses.length -1];++i) {
		   	if (prefMasses.indexOf(i) > -1) {
				vect.push(1);
			}
			else {
				vect.push(0);
			}
		}
		
		return vect;
		
		
		
	};


    
    
    this.spectrum = function (cyclic,retSubPeptides,prefSufOnly) {
        
        //retSubPeptides: flag to indicate whether to also return the sub peptide strings
        //prefSufOnly - only include prefixes and suffices, not "interior" sub-peptides

        var pepStr = this.toShortString('');
        
        var prefixMasses = [0];
        this.peptide.forEach(function (am, i) {
            var prefMass = am.getIntegerWeight() + prefixMasses[i];
            prefixMasses.push(prefMass);

        });

        var peptideMass = prefixMasses[prefixMasses.length - 1];

        var w;
        var spec = [[0,'']];

        for (var st = 0; st < this.peptide.length + 1; ++st) {
            for (var end = st + 1; end < this.peptide.length + 1; ++end) {
                w = prefixMasses[end] - prefixMasses[st];
                if (prefSufOnly) {
                    if ((st == 0) || (end == this.peptide.length)) {
                        spec.push([w,pepStr.substring(st,end)]);
                    }
                }
                else {
                    spec.push([w, pepStr.substring(st, end)]);
                }
                if (st > 0 && end < this.peptide.length) {
                    if (cyclic) {
                        spec.push([peptideMass - w,pepStr.substring(end) + pepStr.substring(0, st)]); //wrap around bit
                    }
                }
            }
        }

        spec.sort(function(a,b) {
            return a[0] - b[0];
        });

        if (retSubPeptides) {
            return spec;
        }
        else {
            return spec.map(function(el) {
                return el[0];
                
            });
        }
       

    };

    this.linearSpectrum = function (prefSufOnly) {
        return this.spectrum(false,false,prefSufOnly);
    };

    this.cyclicSpectrum = function () {
        return this.spectrum(true);
    };


    this.linearScore = function (expSpectrum,prefSufOnly) {

        return this.score(expSpectrum, true,prefSufOnly);

    };

    this.score = function (expSpectrum, linear,prefSufOnly) {

        var isFloatSpectrum;
        var floatErr = 0;
        if (expSpectrum[0] % 1 === 0) {
            isFloatSpectrum = false;
        }
        else {
            isFloatSpectrum = true;
            floatErr = 0.4;
        }

        //linear is optional. Defaults to cyclic

        var sc = 0;

        var theor;
        if (linear) {
            theor = this.linearSpectrum(prefSufOnly);
        }
        else {
            theor = this.cyclicSpectrum();
        }


        var exp = expSpectrum.map(function (el) {
            return el;
        });

        if (isFloatSpectrum) {
            for (var tm = 0;tm < theor.length; ++tm) {
                var theorMass = theor[tm];

                for (var i = 0; i < exp.length; ++i) {
                    if ((theorMass + 1 >= exp[i] - floatErr) && (theorMass + 1 <= exp[i] + floatErr)) {
                 //   if (theorMass + 1 == Math.round(exp[i])) {
                        ++sc;
                        exp.splice(i, 1);
                        break;
                    }
                    else if (theorMass + 1 < exp[i] - floatErr) {
                   // else if (theorMass + 1 < Math.round(exp[i])) {
                        //assumes exp is sorted in ascending ordeer
                        break;
                    }
                }

            }


        }
        else {
            theor.forEach(function (theorMass) {
                var ind = exp.indexOf(theorMass);
                if (ind > -1) {
                    ++sc;
                    exp.splice(ind, 1); //remove, already used
                }

            });

        }


        return sc;

    };
	
	this.scoreAgainstSpectralVector = function(specVector) {
		
		var pepVector = this.toPeptideVector();
		
		if (pepVector.length == specVector.length) {
			
			var score = 0;
			
			specVector.forEach(function(el,i) {
				score += (el * pepVector[i]);
			});
			
			return score;
			
		}
		else {
			return DGraph.infinity * -1;
		
		}
	};

    this.allCycles = function() {
      var cyclesAr = [];  
      for (var i = 0;i < this.peptide.length;++i) {
          cyclesAr.push(this.toString(Peptide.short,'',false,i));
          cyclesAr.push(this.toString(Peptide.short,'',true,i)); //reverse
      }
      return cyclesAr;
    };

    this.init();


}


Peptide.Short = 1;
Peptide.Medium = 2;
Peptide.Long = 3;

Peptide.AminoArrFromStr = function(str) {
    var amAr = [];

    for (var i = 0;i < str.length; ++i) {
        var am = new Amino(str.substring(i,i+1));
        amAr.push(am);
    }

    return amAr;

};

Peptide.AminoArrFromArr = function(arr) {

    var amAr;

    amAr = arr.map(function(el) {
        return new Amino(el);
        
    });
  
    return amAr;


};

Peptide.AminoArrFromWeights = function(weights) {

    var amArr;
    
    amArr = weights.map(function(w) {
        for (var amEntry in Amino.transTable) {
            if (Amino.transTable.hasOwnProperty(amEntry)) {
                if (Amino.transTable[amEntry][2] == w) {
                    return new Amino(amEntry);
                }
            }
        }
        return new Amino(w);

    });
    return amArr;

};

Peptide.AminoArrFromVector = function(vector) {
	
	var amWeightArr = [];
	
	var vecArr = vector.split(' ');
	
	var prevMass = 0;
	
	vecArr.forEach(function(el,i) {
		
		if (el == 1) {
			var thisMass = i + 1 - prevMass;
			amWeightArr.push(thisMass);
			
			prevMass = i + 1;
		}
			
			
			
	});
	
	return Peptide.AminoArrFromWeights(amWeightArr);
	
	
};


Peptide.linear = 1;
Peptide.circular = 2;

Peptide.PepMethodIdealSpectrum = 1;
Peptide.PepMethodSequenceBrute = 2;
Peptide.PepMethodSequenceLeaderboard = 3;
Peptide.PepMethodSequenceLeaderboardConv = 4;
Peptide.PepMethodSequenceGraphBrute = 5;
Peptide.PepMethodSequenceGraphVector = 6;


function RNA(rna) {

    this.rna = rna;

    this.init = function() {

    };

    this.codons = function(frameOffset,rev) {
        if (!frameOffset) {
            frameOffset = 0;
        }
        var arr = [];
        var cod,i;

        if (rev) {
            var strRev = this.rna.split('').reverse().join('');
            for (i = frameOffset; i < this.rna.length - Codon.len + 1; i = i + Codon.len) {
                cod = new Codon(strRev.substring(i, i + Codon.len));
                arr.push(cod);

            }

        }
        else {
            for (i = frameOffset; i < this.rna.length - Codon.len + 1; i = i + Codon.len) {
                cod = new Codon(this.rna.substring(i, i + Codon.len));
                arr.push(cod);

            }

        }
        return arr;

    };

    this.translate = function(frameOffset,forceTranslateStart,forceTranslateStop,rev) {
        var cods = this.codons(frameOffset,rev);
        
        var startPos = frameOffset;
        
        var started = false;
        var stopped = false;

        if (forceTranslateStart) {
            started = true; //ie don't wait for start codon
            
        }

        var tranCods = cods.filter(function(cod,codNum) {

            if (stopped) {
                return false;
            }

            if (!started) {
                if (cod.isStartCodon()) {
                    started = true;
                    startPos = frameOffset + (codNum * Codon.len);
                }
            }

            if (started) {
                if (cod.isStopCodon()) {
                    if (forceTranslateStop) {
                        return true;
                    }
                    else {
                        stopped = true;
                        return false;
                    }
                }
                else {
                    return true;
                }
            }
            else {
                return false;
            }

        });

        var aminos = tranCods.map(function(cod) {
            return cod.translate();

        });
        return new Peptide(aminos,startPos);
    };

    this.init();
}

RNA.Uracil = 'U';

function DNA(dna) {

    this.dna = dna;

    this.init = function() {

    };

    this.rnaTranscript = function() {
        return this.dna.split(DNA.Thymine).join(RNA.Uracil);
    };

    this.complement = function() {

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

        var compl = '';

        for (var i = 0;i < this.dna.length;++i) {
            compl += revCompl[this.dna[i]];
        }

        return new DNA(compl);

    };

    this.reverseComplement = function() {

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

        for (var i = this.dna.length-1;i >= 0;--i) {
            rev += revCompl[this.dna[i]];
        }

        return new DNA(rev);
    };

    this.init();
}

DNA.Thymine = 'T';

DNA.TransMethodTranscribe = 1;
DNA.TransMethodTranslate = 2;
DNA.TransMethodRetro = 3;



function Spectrum(spectrumString) {
    //floatSpectrum flag = spectrum is float not int

    if (spectrumString.indexOf('.') > -1) {
        this.isFloatSpectrum = true;
        this.floatErr = 0.3;
    }
    else {
        this.isFloatSpectrum = false;
        this.floatErr = 0;
    }

    var svdThis = this;
    this.specAr = spectrumString.split(' ');
    this.specAr = this.specAr.map(function(el) {
           return svdThis.isFloatSpectrum ?  parseFloat(el) : parseInt(el);
          });
    
    this.init = function() {
        
    };
    
    this.convolutionStr = function() {
        var convAr = this.convolution();
        
        var str = '';
        convAr.forEach(function(el) {
            str+= el + ' ';
         });
        return str;
    };

    this.topElements = function(N) {
        
        var topEls = [];
        
        var conv = this.convolution();
        
        var prev = conv[0];
        
        var currCount = 0;
        var prevCount = -1;
        var NCount = -1;
        
        for (var i =0;i < conv.length;++i) {
            if (conv[i] == prev) {
                ++currCount;
                
            }
            else {
  
                
                if ((prev >= Amino.LowestWeight) && (prev <= Amino.HighestWeight)) {
                    if (topEls.length > N) {
                        if (currCount == NCount) {
                            topEls.push(prev); //tie
                        }
                        else {
                            break;
                        }
                    }
                    else {
                        if (topEls.length == N) {
                            NCount = currCount;
                        }
                        topEls.push(prev);
                    }
                }

                prev = conv[i];
                prevCount = currCount;
                currCount = 1;
            }
            
        }

        return topEls;
        

    };

    this.intsBetween = function(a,b) {

        if ((a == b) && (Math.ceil(a) == a)) {
            //if a and b are both the same integer, then return that integer
            return [a];
        }

         var nums = [];
         for (var i = Math.ceil(a);i < b;++i) {
            nums.push(i);
        }
        return nums;


    };
    
    this.convolution = function() {
        
        //var svdThis = this;
        
        var convo = [];
        
        //this.specAr.forEach(function(el1) {
        for (var i = 0;i < this.specAr.length; ++i) {
            for (var j = i+1;j < this.specAr.length;++j) {
                var diff =  Math.abs(this.specAr[i] - this.specAr[j]);
                if (this.isFloatSpectrum) {
                    var possNums = this.intsBetween(diff - this.floatErr*2,diff + this.floatErr*2);
                    possNums.forEach(function(el) {
                        if (el != 0) {
                            convo.push(el);
                        }

                    });

                }
                else {
                    if (diff != 0 ) {
                        convo.push(diff);
                    }

                }

            }
        }

       // convo.sort(function(a,b) {
       //     return a - b;
       // });

         var convCounts = {};

        convo.forEach(function(el) {
            if (el in convCounts) {
                ++convCounts[el];

            }
            else {
                convCounts[el] = 1;
            }

        });
        var convWithCounts = convo.map(function(el) {
            return [el,convCounts[el]];

        });



        convWithCounts.sort(function(a,b) {
            var diff = b[1] - a[1];
            if (diff == 0) {
                return a[0] - b[0];
            }
            else {
                return b[1] - a[1];
            }

        });

        convo = convWithCounts.map(function(el) {
            return el[0];
        });


        return convo;
    };
    
    this.init();
    
}

Spectrum.PepMethodIdealSpectrum = 1;
Spectrum.PepMethodSequenceBrute = 2;


function Aligner(p,t,m) {
    //p = pattern, t = whole text, m = num mismatches

    this.p = p;
    this.t = t;
    this.pLen = this.p.length;
    this.tLen = this.t.length;

    this.m = m ? m : 0;
    
    this.debugComparisons = 0;
    this.debugTime = 0;

    this.initBadChars = function(p) {
        //this.badChars = [];
        var badChars = [];
        for (var i = 0; i < p.length;++i) {
            var thisChar = p[i];
            var shiftDict = {};
            c_Bases.forEach(function(el) {
                shiftDict[el] = i+1;
            });
            for (var j = i;j >= 0;--j) {
                if (shiftDict[p[j]] < i+1) {
                    //already found
                }
                else {
                    shiftDict[p[j]] = i - j;
                }
            }
            badChars.push(shiftDict);

        }

        return badChars;


    };

    this.initGoodSuffs = function(p) {
        //this.goodSuffs = [];

        //goodsuffs has an extra cell at the end for empty substring (no match)
        goodSuffs = [];

        for (var i = 0;i < p.length+1;++i) {
            goodSuffs[i] = -1;
        }

        goodSuffs[p.length] = 0; //empty

        for (var i = p.length-1;i >= 0;--i) {
            var suff = p.substring(i,p.length);
            var latestMatch = -1;
            for (var j = 0;j < p.length - suff.length;++j) {
                if (p.substring(j,j+suff.length) == suff) {
                    latestMatch = j;

                }
            }
            if (latestMatch == -1) {
                //check for prefixes
                var foundPrefix = false;
                for (var k = suff.length - 1;k >= 1;--k) {
                    if (p.substring(0,k) == suff.substring(suff.length - k,suff.length)) {
                        // this.goodSuffs[suff] = i + (suff.length - k);
                        goodSuffs[i] = i + (suff.length - k);
                        foundPrefix = true;
                        break;
                    }

                }
                if (!foundPrefix) {
                    //this.goodSuffs[suff] = this.p.length;
                    goodSuffs[i] = p.length;
                }
                //this.goodSuffs[suff] = -1;
            }
            else {
                //this.goodSuffs[suff] = i - latestMatch;
                goodSuffs[i] = i - latestMatch;
            }
        }

        return goodSuffs;


    };

    this.initPenaltyDictEdit = function() {

        var penDict = {};

        var bases = [];

        c_Bases.forEach(function(base) {
            bases.push(base);
        });

        bases.push('_');

        bases.forEach(function(rowBase) {
            var rowDict = {};
            bases.forEach(function(colBase) {
                if (rowBase == colBase) {
                    rowDict[colBase] = 0;
                }
                else {
                    rowDict[colBase] = -1;
                }
            });
            penDict[rowBase] = rowDict;

        });

        return penDict;




    };

    this.initPenaltyDictLCS = function() {

        var penDict = {};

        var bases = [];

        c_Bases.forEach(function(base) {
            bases.push(base);
        });

        bases.push('_');

        bases.forEach(function(rowBase) {
            var rowDict = {};
            bases.forEach(function(colBase) {
                if (rowBase == colBase) {
                    rowDict[colBase] = 1;
                }
                else {
                    rowDict[colBase] = 0;
                }
            });
            penDict[rowBase] = rowDict;

        });

        return penDict;




    };

    this.initPenaltyDictDNA = function(indelPen) {
        var dnaInputColHeadsStr = "A C G T _";
        var dnaInputColHeads = dnaInputColHeadsStr.split(' ');
        var dnaInputStr =
                "A 0 -4 -2 -4"
            + "\nC -4 0 -4 -2"
            + "\nG -2 -4 0 -4"
            + "\nT -4 -2 -4 0";


        var penDict = {};

        var dnaInput = dnaInputStr.split('\n');
        dnaInput = dnaInput.map(function(el) {
            var els = el.substring(2,el.length).split(' ');

            els = els.map(function(e) {
                return parseInt(e);
            });
            els.push(indelPen);
            return [el[0],els];

        });
        var inRow = {};
        dnaInputColHeads.forEach(function(colHead) {
            inRow[colHead] = indelPen;
        });


        dnaInput.forEach(function(el) {
            var colDict = {};
            dnaInputColHeads.forEach(function(colHead,i) {
                colDict[colHead] = el[1][i];
            });
            penDict[el[0]] = colDict;

        });
        penDict['_'] = inRow;

        return penDict;

    };

    this.initPenaltyDictDNALoc = function(indelPen) {
        var dnaInputColHeadsStr = "A C G T _";
        var dnaInputColHeads = dnaInputColHeadsStr.split(' ');
        var dnaInputStr =
            "A 2 -4 -4 -4"
            + "\nC -4 2 -4 -4"
            + "\nG -4 -4 2 -4"
            + "\nT -4 -4 -4 2";


        var penDict = {};

        var dnaInput = dnaInputStr.split('\n');
        dnaInput = dnaInput.map(function(el) {
            var els = el.substring(2,el.length).split(' ');

            els = els.map(function(e) {
                return parseInt(e);
            });
            els.push(indelPen);
            return [el[0],els];

        });
        var inRow = {};
        dnaInputColHeads.forEach(function(colHead) {
            inRow[colHead] = indelPen;
        });


        dnaInput.forEach(function(el) {
            var colDict = {};
            dnaInputColHeads.forEach(function(colHead,i) {
                colDict[colHead] = el[1][i];
            });
            penDict[el[0]] = colDict;

        });
        penDict['_'] = inRow;

        return penDict;

    };

    this.initPenaltyDictBlosum = function(indelPen) {

        var bloInputColHeadsStr = "A C D E F G H I K L M N P Q R S T V W Y _";
        var bloInputColHeads = bloInputColHeadsStr.split(' ');
        var blosumInputStr =
            "A 4 0 -2 -1 -2 0 -2 -1 -1 -1 -1 -2 -1 -1 -1 1 0 0 -3 -2"
        + "\nC 0 9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2"
        + "\nD -2 -3 6 2 -3 -1 -1 -3 -1 -4 -3 1 -1 0 -2 0 -1 -3 -4 -3"
        + "\nE -1 -4 2 5 -3 -2 0 -3 1 -3 -2 0 -1 2 0 0 -1 -2 -3 -2"
        + "\nF -2 -2 -3 -3 6 -3 -1 0 -3 0 0 -3 -4 -3 -3 -2 -2 -1 1 3"
        + "\nG 0 -3 -1 -2 -3 6 -2 -4 -2 -4 -3 0 -2 -2 -2 0 -2 -3 -2 -3"
        + "\nH -2 -3 -1 0 -1 -2 8 -3 -1 -3 -2 1 -2 0 0 -1 -2 -3 -2 2"
        + "\nI -1 -1 -3 -3 0 -4 -3 4 -3 2 1 -3 -3 -3 -3 -2 -1 3 -3 -1"
        + "\nK -1 -3 -1 1 -3 -2 -1 -3 5 -2 -1 0 -1 1 2 0 -1 -2 -3 -2"
        + "\nL -1 -1 -4 -3 0 -4 -3 2 -2 4 2 -3 -3 -2 -2 -2 -1 1 -2 -1"
        + "\nM -1 -1 -3 -2 0 -3 -2 1 -1 2 5 -2 -2 0 -1 -1 -1 1 -1 -1"
        + "\nN -2 -3 1 0 -3 0 1 -3 0 -3 -2 6 -2 0 0 1 0 -3 -4 -2"
        + "\nP -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2 7 -1 -2 -1 -1 -2 -4 -3"
        + "\nQ -1 -3 0 2 -3 -2 0 -3 1 -2 0 0 -1 5 1 0 -1 -2 -2 -1"
        + "\nR -1 -3 -2 0 -3 -2 0 -3 2 -2 -1 0 -2 1 5 -1 -1 -3 -3 -2"
        + "\nS 1 -1 0 0 -2 0 -1 -2 0 -2 -1 1 -1 0 -1 4 1 -2 -3 -2"
        + "\nT 0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1 0 -1 -1 -1 1 5 0 -2 -2"
        + "\nV 0 -1 -3 -2 -1 -3 -3 3 -2 1 1 -3 -2 -2 -3 -2 0 4 -3 -1"
        + "\nW -3 -2 -4 -3 1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11 2"
        + "\nY -2 -2 -3 -2 3 -3 2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1 2 7";

        var penDict = {};

        var blosumInput = blosumInputStr.split('\n');
        blosumInput = blosumInput.map(function(el) {
            var els = el.substring(2,el.length).split(' ');

            els = els.map(function(e) {
                return parseInt(e);
            });
            els.push(indelPen);
            return [el[0],els];

        });
        var inRow = {};
        bloInputColHeads.forEach(function(colHead) {
            inRow[colHead] = indelPen;
        });


        blosumInput.forEach(function(el) {
            var colDict = {};
           bloInputColHeads.forEach(function(colHead,i) {
               colDict[colHead] = el[1][i];
           });
           penDict[el[0]] = colDict;

        });
        penDict['_'] = inRow;

        return penDict;
    };

    this.pen = function(ch1,ch2,type,localFlag) {

        if (type == 'blo') {
            return this.penaltyDictBlosum[ch1][ch2];
        }
        
        if (type == 'LCS') {
            return ch1 == ch2 ? 1 : 0;
        }
        
        if (type == 'dna') {
            if (localFlag) {
                return this.penaltyDictDNALoc[ch1][ch2];
            }
            else {
                return this.penaltyDictDNA[ch1][ch2];
            }
        }

        if ( ((c_Bases.indexOf(ch1) > -1)  || (ch1 == '_')) &&
            ((c_Bases.indexOf(ch2)  > -1) || (ch2 == '_'))) {


        }
        else {
            //contains non DNA. Use simple
            if (ch1 == ch2) {
                return 0;
            }
            else {
                return -1;
            }
        }

        return this.penaltyDict[ch1][ch2];

    };

    this.alignInit = function(localFlag,alignType,indelPen) {

        //Alignment pre-process
        this.scoreMatrix = [];
        
        for (var i = 0;i < this.p.length + 1;++i) {
            var rowAr = [];
            if (i == 0) {

                for (var j = 0;j < this.t.length+1;++j) {
                    rowAr.push(localFlag ? 0 : j * indelPen);
                }
            }
            else {
                if (alignType ==Aligner.EditDistType) {
                    rowAr.push(i * indelPen);
                }
                else {
                    if (localFlag) {
                        rowAr.push(0);
                    }
                    else {
                        rowAr.push(i * indelPen);
                    }
                }
            }
            this.scoreMatrix.push(rowAr);

        }


    };

    this.editDistInit = function(localFlag) {
        for (var i = 0;i < this.p.length + 1;++i) {
            var rowAr = [];
            if (i == 0) {

                for (var j = 0;j < this.t.length+1;++j) {
                    rowAr.push(localFlag ? 0 : j * -1);
                }
            }
            else {
                rowAr.push(i * -1);
            }
            this.scoreMatrix.push(rowAr);

        }

    };

    this.score = function(ch1,ch2) {


    };

    this.init = function() {

 

   
       // this.editDistInit(true);
        /*
         for (var i = 0;i < this.p.length + 1;++i) {
         var rowAr = [];
         for (var j = 0;j < this.t.length)
         }
         */
        
        //Boyer Moore pre-process
        
        this.badChars = this.initBadChars(this.p);
        this.goodSuffs = this.initGoodSuffs(this.p);

        this.badCharsPart = [];
        this.goodSuffsPart = [];
        this.parts = [];
        this.partOffsets = [];

        if (this.m > 0) {

            var numParts = this.m + 1;
            var partSize = Math.floor(this.p.length / numParts);
            var mod = this.p.length % numParts;
            for (var i = 0;i < numParts;++i) {
                var len;
                if (i == numParts - 1) {
                    len = partSize + mod;
                }
                else {
                    len = partSize;
                }
                var part = this.p.substring(i*partSize,i*partSize + len);
                this.badCharsPart.push(this.initBadChars(part));
                this.goodSuffsPart.push(this.initGoodSuffs(part));
                this.parts.push(part);
                this.partOffsets.push(i*partSize);

            }
        }

        else {
            //whole only, no mismatches
            this.badCharsPart.push(this.badChars);
            this.goodSuffsPart.push(this.goodSuffs);
            this.parts.push(this.p);
            this.partOffsets.push(0);

        }
        
    };
    
    this.toText = function() {
        return 'len p: ' + this.p.length + ' len t: ' + this.t.length;
    };
    
    this.naive = function() {

        return this.naiveWithMismatch(0);
        
    };
    
    this.matchString = function(posns) {
        var svdThis = this;
        return posns.map(function(el) {
            return svdThis.t.substring(el,el+svdThis.p.length);
        });
        
        
    };
    
    this.naiveWithMismatch = function(maxMismatch) {
        var start = performance.now();
   
        var posns = [];
        this.debugComparisons = 0;

        for (var i = 0;i < this.t.length - this.p.length +1;++i) {
            var mismatches = 0;
            for (var j = 0;j < this.p.length;++j) {
                ++ this.debugComparisons;
                if (this.p[j] != this.t[i+j]) {
                    ++mismatches;
                    if (mismatches > maxMismatch) {
                        break;
                    }
                }
            }
            if (mismatches <= maxMismatch) {
                posns.push(i);
            }
        }

        var end = performance.now();
        this.debugTime = end - start;
        

        return posns;       
        
    };

    this.bmFind = function(part,partNum) {

        var posns = [];
        var pLen = part.length;

        var currPos = pLen - 1;

        while (currPos <= this.tLen - 1) {
           // console.log('now: ' + this.t.substring(currPos-part.length+1,currPos+1));
            var pOffset = pLen - 1;
            var match = true;
            for (var tOffset = currPos; tOffset >= currPos - pLen + 1;--tOffset) {
                ++this.debugComparisons;
                if (this.t[tOffset] != part[pOffset]) {
                    var skipBad = this.badCharsPart[partNum][pOffset][this.t[tOffset]];
                    var skipGoodSuff = this.goodSuffsPart[partNum][pOffset+1];
                    //console.log('skip bad: ' + skipBad + ' skip good: ' + skipGoodSuff);
                    /*
                     if (pOffset == this.pLen - 1) {
                     //nothing matched
                     skipGoodSuff = 0;
                     }
                     else {
                     //skipGoodSuff = this.goodSuffs[this.p.substring(pOffset+1,this.p.length)];
                     skipGoodSuff = this.goodSuffs[pOffset+1];
                     }
                     */
                    currPos += Math.max(skipBad,skipGoodSuff);
                    match = false;
                    break;
                }
                --pOffset;
            }
            if (match) {
                posns.push(currPos - pLen + 1);
                ++currPos;
            }

        }

        return posns;


    };

    this.boyerMoore = function() {
        //var start = new Date().getTime();

        var start = performance.now();
        this.debugComparisons = 0;

        var posns = [];
        var allPosns = new Set();
        this.debugComparisons = 0;

        var svdThis = this;

        this.parts.forEach(function(prt,i) {
            var posns = svdThis.bmFind(prt,i);
            posns.forEach(function(pos) {
               var tStart = pos - svdThis.partOffsets[i];
                var tEnd = tStart + svdThis.p.length;
                if (tEnd > svdThis.t.length) {

                }
                else {
                    if (hamDist(svdThis.t.substring(tStart, tStart + svdThis.p.length), svdThis.p) <= svdThis.m) {
                        //if (svdThis.t.substring(tStart,svdThis.p.length) == svdThis.p ) {
                        //allPosns.push(tStart);
                        allPosns.add(tStart);
                    }
                }


            });
            //allPosns = allPosns.concat(posns);
        });

       
        var end = performance.now();
        this.debugTime = end - start;

        var allPosnsList = Array.from(allPosns);
        return allPosnsList;


    };

    this.editDelta = function(a,b) {
        return a == b ? 0 : 1;

    };

    this.align = function(localFlag,alignType,indelPen) {
 
        //var editDistType = 1;
        //var alignBlo = 2;
        //var alignDNA = 3;
        
        var type;
        switch (alignType) {
            case Aligner.EditDistType:
                type = '';
                break;
            case Aligner.AlignBlo:
                type = 'blo';
                break;
            case Aligner.AlignDNA:
                type = 'dna';
                break;
            case Aligner.LCS:
                type = 'LCS';
                break;
            default:
                break;
  
        }
   
        
        this.penaltyDictBlosum = this.initPenaltyDictBlosum(indelPen);
        this.penaltyDictDNA = this.initPenaltyDictDNA(indelPen);
        this.penaltyDictDNALoc = this.initPenaltyDictDNALoc(indelPen);
        this.penaltyDictLCS = this.initPenaltyDictLCS();
        this.penaltyDict = this.initPenaltyDictEdit();

        
        this.alignInit(localFlag,alignType,alignType == Aligner.EditDistType ? -1 : indelPen);

        for (var i = 1;i < this.p.length + 1;++i) {
            for (var j = 1;j < this.t.length + 1;++j) {
                var curr = Math.max(
                    this.scoreMatrix[i][j-1]+ this.pen('_',this.t[j-1],type,localFlag),
                    this.scoreMatrix[i-1][j]+ this.pen(this.p[i-1],'_',type,localFlag),
                    this.scoreMatrix[i-1][j-1] + this.pen(this.p[i-1],this.t[j-1],type,localFlag)

                );
                if (localFlag) {
                    if (alignType == Aligner.EditDistType) {

                    }
                    else {
                        curr = Math.max(0,curr);
                    }
                }

                this.scoreMatrix[i].push(curr);


            }
        }

        if ((localFlag) && (alignType == Aligner.EditDistType)) {
            var highestScore = -99999;
            var highestCol = -99999;
            var highestCols = [];
            for (var c = 0; c < this.t.length + 1; ++c) {
                if (this.scoreMatrix[this.p.length][c] > highestScore) {
                    highestScore = this.scoreMatrix[this.p.length][c];
                    highestCols = [c];
                }
                else if (this.scoreMatrix[this.p.length][c] == highestScore) {
                    highestCols.push(c);
                }
            }
            var res = [];
            var svdThis = this;
            highestCols.forEach(function(highestCol) {
                res.push(svdThis.backtrack(svdThis.p.length, highestCol,localFlag,type));
            });

            //return [lowestScore, res[0], res[1],res[2]];
            return [highestScore,res];
        }
        else if (localFlag) {
            var highestScore = -99999;
            var highestCol = -99999;
            var highestRow = -99999;
            var highestPosns = [];
            for (var r = 0;r < this.p.length + 1;++r) {
                for (var c = 0; c < this.t.length + 1; ++c) {
                    if (this.scoreMatrix[r][c] > highestScore) {
                        highestScore = this.scoreMatrix[r][c];
                        highestPosns = [[r,c]];
                    }
                    else if (this.scoreMatrix[r][c] == highestScore) {
                        highestPosns.push([r,c]);
                    }
                }
            }
            var res = [];
            var svdThis = this;
            highestPosns.forEach(function(highestPos) {
                res.push(svdThis.backtrack(highestPos[0], highestPos[1],localFlag,type));
            });

            //return [lowestScore, res[0], res[1],res[2]];
            return [highestScore,res];


        }
        else {
            var res = this.backtrack(this.p.length, this.t.length,localFlag,type);
            return [this.scoreMatrix[this.p.length][this.t.length], res[0], res[1],res[2]];

        }



    };

    this.editDist = function(localFlag) {
      this.editDistInit(localFlag);

      for (var i = 1;i < this.p.length + 1;++i) {
          for (var j = 1;j < this.t.length + 1;++j) {
              /*
              var curr = Math.min(
                  this.scoreMatrix[i][j-1]+1,
                  this.scoreMatrix[i-1][j]+1,
                  this.scoreMatrix[i-1][j-1] + this.editDelta(this.p[i-1],this.t[j-1])
              );
              */
              var curr = Math.max(
                  this.scoreMatrix[i][j-1]+ this.pen('_',this.p[j-1]),
                  this.scoreMatrix[i-1][j]+ this.pen(this.p[i-1],'_'),
                  this.scoreMatrix[i-1][j-1] + this.pen(this.p[i-1],this.t[j-1])
                  // this.scoreMatrix[i-1][j-1] + this.penaltyDict[this.p[i-1]][this.t[j-1]]
                 // this.scoreMatrix[i][j-1]+ this.penaltyDict['_'][this.t[j-1]],
                 // this.scoreMatrix[i-1][j]+ this.penaltyDict[this.p[i-1]]['_'],
                 // this.scoreMatrix[i-1][j-1] + this.penaltyDict[this.p[i-1]][this.t[j-1]]
              );
              this.scoreMatrix[i].push(curr);

          }
      }

      if (localFlag) {
          var highestScore = -99999;
          var highestCol = -99999;
          var highestCols = [];
          for (var c = 0; c < this.t.length + 1; ++c) {
              if (this.scoreMatrix[this.p.length][c] > highestScore) {
                  highestScore = this.scoreMatrix[this.p.length][c];
                  highestCols = [c];
              }
              else if (this.scoreMatrix[this.p.length][c] == highestScore) {
                  highestCols.push(c);
              }
          }
          var res = [];
          var svdThis = this;
          highestCols.forEach(function(highestCol) {
              res.push(svdThis.backtrack(svdThis.p.length, highestCol,localFlag));
          });

          //return [lowestScore, res[0], res[1],res[2]];
          return [highestScore,res];
      }
      else {
          var res = this.backtrack(this.p.length, this.t.length,localFlag);
          return [this.scoreMatrix[this.p.length][this.t.length], res[0], res[1],res[2]];

      }

    };

    this.backtrack = function(r,c,localFlag,type) {
        //backtrack from r, c
        var i = r;
        var j = c;
        var pAligned = '';
        var tAligned = '';
        
        var tStart = 0; // used for local edit dist

        while ((i >= 1) || (j >=1)) {

            if (localFlag) {
                if (type == '') {

                }
                else {
                    if (this.scoreMatrix[i][j] == 0) {
                        //stop when 0 found if doing local alignment (except for edit dist)
                        tStart = j;
                        break;
                    }
                }
            }
            var pGap,tGap,match;
            if (j == 0) {
                pGap = -999999;
                tGap =  this.scoreMatrix[i-1][j]; //vertical
                match = -999999;
            }
            else if (i == 0) {
                pGap = this.scoreMatrix[i][j - 1]; //horizontal
                tGap = -999999;
                match = -999999;
            }
            else {
                pGap = this.scoreMatrix[i][j - 1]; //horizontal
                tGap = this.scoreMatrix[i - 1][j]; //vertical
                match = this.scoreMatrix[i - 1][j - 1];
            }
            if (j == 0) {
                pAligned = this.p[i - 1] + pAligned;
                tAligned = '-' + tAligned;
                --i;
                continue;
            }
            if (i == 0) {
                if (localFlag) {
                    //No need for P gaps at start
                    tStart = j;
                    break;
                }
                else {
                    pAligned = '-' + pAligned;
                    tAligned = this.t[j - 1] + tAligned;

                }
             --j;
             continue;

            }


            var pGapScore = pGap + this.pen('_',this.t[j-1],type,localFlag);
            var tGapScore = tGap + this.pen(this.p[i-1],'_',type,localFlag);
            var matchScore = match + this.pen(this.p[i-1],this.t[j-1],type,localFlag);

           // var max = Math.max(pGap,tGap,match);
           // if (max == match) { //first priority is match
           if (matchScore == this.scoreMatrix[i][j]) { //first priority is a match
                pAligned = this.p[i-1] + pAligned;
                tAligned = this.t[j-1] + tAligned;
                --i;
                --j;
            }
           // else if (max == pGap) {
           else if (pGapScore == this.scoreMatrix[i][j]) {

                if ((localFlag) && (i == 0)) {
                    //No need for P gaps at start
                    tStart = j;
                    break;
                }
                else {
                    pAligned = '-' + pAligned;
                    tAligned = this.t[j - 1] + tAligned;

                }



                --j;
            }
            else {
                pAligned = this.p[i-1] + pAligned;
                tAligned = '-' + tAligned;

                --i;
            }
        }

        return [pAligned,tAligned,tStart];

    };
    
    

    this.init();
}

Aligner.EditDistType = 1;
Aligner.AlignBlo = 2;
Aligner.AlignDNA = 3;
Aligner.LCS = 4;

function allPossiblePeptidesWithWeight(w) {

    var i;

    var wMatrix = [];
    for (i = 0;i <= w;++i) {
        wMatrix[i] = 0;
    }

    var aminoWeights = Amino.AminoWeights();


    //var possWeightsPrevLen = [];
   // var totWeightsMatched = 0;

    var done = false;

    //for (var len = 1;len < 25; ++len) {
    var round = 0;



    while (!done) {
        ++round;
        console.log('round: ' + round);
       // var wMatrixThisLen  = [];
      //  var weightsMatchedThisLen = 0;

        var wNewMatrix = [];
        for (i = 0;i <= w;++i) {
            wNewMatrix[i] = 0;
        }

        if (round == 1) {
            aminoWeights.forEach(function(el) {
                if (el > w) {
                  wNewMatrix[el] = 0;
                }
                else {
                    ++wNewMatrix[el];
                }

            });

        }
        else {

            wNewMatrix = wMatrix.map(function(el) {
                return el;
            });

            wMatrix.forEach(function(matCount,matInd) {
                if (matCount == 0) {
                    //wNewMatrix[matInd] = 0;

                }
                else if (matInd == w) {
                    /*
                    if (wNewMatrix[matInd] == 0) {
                        wNewMatrix[matInd] = wMatrix[matInd];
                    }
                    */
                }
                else {
                    var subtr = 0;
                    aminoWeights.forEach(function (amW,amInd) {
                        var thisW = matInd + amW;
                        if ((thisW > w) && (amInd == 0)) {
                            //assumes aminoweights are sorted, thus as soon as one blows, all blow

                                wNewMatrix[matInd] = 0;

                        }
                        else if (thisW > w) {

                        }
                        else {
                            if (thisW == w) {
                                console.log('thisw = w: ' + w +  ' matind: ' + matInd + ' amw: ' + amW + ' matCount: ' + matCount);
                            }
                         //   if (wNewMatrix[thisW] == 0) {
                                wNewMatrix[thisW]  = wNewMatrix[thisW] + matCount;
                                subtr = matCount;

                          //  }
                           // else {
                           //     wNewMatrix[thisW]  = wNewMatrix[thisW] + matCount;
                           // }

                        }
                    });
                    wNewMatrix[matInd] -= subtr;
                }



            });

        }

        var totLeft = 0;

        wMatrix = wNewMatrix.map(function(el,j) {
            if (j == wNewMatrix.length - 1) {

            }
            else {
                totLeft += el;
            }
           return el;
        });


        console.log('count at end of round ' + round + ' ' + wMatrix[w]);

        if (totLeft == 0) {
            done = true;
        }
    }

    return wMatrix[w];
}

function cyclopeptideSequencing(expSpectrum,cyclicFlag,progCallback,prefSufOnly) {
    //actually caters for linear peptide  (cyclic flag false) or cyclic peptide (cyclic flag true)
    
    var expSpectrumAr = expSpectrum.split(' ');
    expSpectrumAr = expSpectrumAr.map(function (el) {
        return parseInt(el);
    });



    var done = false;

    var goodOnes = [];

    var parentWeight = expSpectrumAr[expSpectrumAr.length - 1];

    var expMatrix = [];
    for (var i = 0;i < parentWeight + 1;++i) {
        var yep = (expSpectrumAr.indexOf(i) > -1);
        if (yep) {
            expMatrix.push(1);
        }
        else {
            expMatrix.push(0);
        }
    }


    //var candidates = Amino.AllAminos();

    var allAminoWeights = Amino.AminoWeights();


    var candidates = Amino.AminoWeights();

    candidates = candidates.filter(function (el) {
       // if (expSpectrumAr.indexOf(el) > -1) {
        return (expMatrix[el] == 1);

    });

    candidates = candidates.map(function(el) {
        return [el];
    });




    var round = 0;
    var progThreshold = 100;


    while (!done) {
        ++round;



        //bound

        /*
        candidates = candidates.filter(function (cand,candInd) {
            var candW = cand.reduce(function (a, b) {
                return a + b;
            });
            // var pep = new Peptide(Peptide.AminoArrFromStr(cand));

            if (progCallback) {
                if ((round % progThreshold == 0) ) {
                    progCallback('Process bo [' + candInd + '/' + candidates.length + '] round ',round,'unknown');
                }
            }

            if (expSpectrumAr.indexOf(candW) > -1) {

                if (candW == parentWeight) {
                    console.log('aha');
                    var candPep = new Peptide(Peptide.AminoArrFromWeights(cand));
                    var candSpec = candPep.cyclicSpectrum();
                    var match = true;
                    if(candSpec.length !== expSpectrumAr.length) {
                        match = false;

                    }
                    else {
                        for(var ii = candSpec.length; ii--;) {
                            if(candSpec[ii] !== expSpectrumAr[ii]) {
                                match = false;
                                break;
                            }
                        }

                    }
                    if (match) goodOnes.push(cand);

                    return false;
                }
                else {
                    return true;
                }
            }
            else {
                return false;

            }
        });

        if (candidates.length == 0) {
            done = true;
            break;

        }
        */

        //branch


        var newCandidates = [];
       // candidates.forEach(function (cand,candInd) {
        for (var c = 0;c< candidates.length; ++c) {
            var cand = candidates[c];
            var candInd = c;

            if (progCallback) {
                if ((candInd % progThreshold == 0)) {
                    progCallback('Process [' + candInd + '/' + candidates.length + '] round ', round, '');
                }
            }

            //allAminoWeights.forEach(function (amW) {
            for (var am = 0; am < allAminoWeights.length; ++am) {
                var amW = allAminoWeights[am];
                var newCandW = 0;

                var newCand = new Array(cand.length);
                var jj = cand.length;
                while (jj--) {
                    newCand[jj] = cand[jj];
                    newCandW += cand[jj];
                }
                /*
                 var newCand = cand.map(function (el) {
                 newCandW += el;
                 return el;
                 });
                 */
                newCand.push(amW);
                newCandW += amW;

                /*
                 var newCandW = newCand.reduce(function (a, b) {
                 return a + b;
                 });
                 */

                if (newCandW > parentWeight) {
                    break; //assumes amino weights are in ascending order, so if one blows, the rest blow
                }

                // if (expSpectrumAr.indexOf(newCandW) > -1) {
                if (expMatrix[newCandW] == 1) {

                    if (newCandW == parentWeight) {
                       // console.log('aha');
                        var candPep = new Peptide(Peptide.AminoArrFromWeights(newCand));
                        var candSpec =  candPep.spectrum(cyclicFlag,false,prefSufOnly); ///candPep.cyclicSpectrum();
                        var match = true;
                        if (candSpec.length !== expSpectrumAr.length) {
                            match = false;

                        }
                        else {
                            for (var ii = candSpec.length; ii--;) {
                                if (candSpec[ii] !== expSpectrumAr[ii]) {
                                    match = false;
                                    break;
                                }
                            }

                        }
                        if (match) goodOnes.push(newCand);


                    }
                    else {
                        newCandidates.push(newCand);
                    }
                }
                else {
                    //do nothing
                }
                //});
            }


            //});
        }

        if (newCandidates.length == 0) {
            done = true;
        }
        else {
            candidates = newCandidates.map(function (el) {
                return el;
            });
        }

    }


    var goodStrings = '';

    var str;
    goodOnes.forEach(function(goodOne) {
        str = '';
        goodOne.forEach(function(weight,weightInd) {
            str+=weight;

            if (weightInd == goodOne.length - 1) {

            }
            else {
                str+='-';

            }


        });
        goodStrings+= str + ' ';

    });


    return goodStrings;
        //done = true;


        //aminoWeights.forEach(function(am) {

        //});




}

function leaderboardCyclopeptideSequencing(expSpectrum,M,N,includeAll200,cyclicFlag,useTheseAminos,progCallback,prefSufOnly) {

    var linearFlag = !cyclicFlag;

    var isFloatFlag;
    var floatErr = 0.0;

    if (expSpectrum.indexOf('.') > -1) {
        isFloatFlag = true;
        floatErr = 0.3;
    }
    else {
        isFloatFlag = false;
    }

    var expSpectrumAr = expSpectrum.split(' ');
    expSpectrumAr = expSpectrumAr.map(function (el) {
        return isFloatFlag ? parseFloat(el) : parseInt(el);
    });
    expSpectrumAr.sort(function(a,b) {
        return a - b;
    });

    
    var  allAms = Amino.AllUniqueWeightAminos(includeAll200);

    var spec = new Spectrum(expSpectrum);
    var topEls = spec.topElements(M);

    if ((useTheseAminos) && useTheseAminos.length > 0) {
        var usePep = new Peptide(Peptide.AminoArrFromStr(useTheseAminos));
        var useWeights = usePep.toWeightString('-').split('-');
        useWeights = useWeights.map(function(el) {
            return parseInt(el);
        });
        topEls = useWeights.slice(0,M);
    }


    if ((includeAll200) && M <144){
        //reduce to top M
        allAms = allAms.filter(function (am) {

            var w = new Amino(am).getIntegerWeight();
            return (topEls.indexOf(w) > -1);

        });
    }

    var done = false;

    //var goodOnes = [];

    var parentWeight = expSpectrumAr[expSpectrumAr.length - 1];

    /*
    var expMatrix = [];
    for (var i = 0;i < parentWeight + 1;++i) {
        var yep = (expSpectrumAr.indexOf(i) > -1);
        if (yep) {
            expMatrix.push(1);
        }
        else {
            expMatrix.push(0);
        }
    }
    */

    var candidates = [new Peptide(Peptide.AminoArrFromStr(''))];


/*
    //var candidates = Amino.AllAminos();

    var allAminoWeights = Amino.AminoWeights();


    var candidates = Amino.AminoWeights();

    candidates = candidates.filter(function (el) {
        // if (expSpectrumAr.indexOf(el) > -1) {
        if (expMatrix[el] == 1) {
            return true;
        }
        else {
            return false;
        }

    });

    candidates = candidates.map(function(el) {
        return [el];
    });

*/


    var round = 0;
    var progThreshold = 10;

    var highScore = 0;
    var bestCands = [];


    while (!done) {
        ++round;


        //branch

        var candProgAr = candidates.map(function(el) {
            return el.toShortString('') + ' ' + el.getIntegerWeight() + ' ' +  el.linearScore(expSpectrumAr,prefSufOnly);
        });
        var bestProgAr = bestCands.map(function(el) {
            return el.toShortString('') + ' '  + el.getIntegerWeight() + ' '+ el.score(expSpectrumAr,false,prefSufOnly);
        });


        var newCandidates = [];
        // candidates.forEach(function (cand,candInd) {
        for (var c = 0; c < candidates.length; ++c) {




            var cand = candidates[c];
            var candInd = c;

            if (progCallback) {
                if (((candInd % progThreshold == 0)) || candInd == candidates.length - 1) {
                    /*
                    var candProgAr = candidates.map(function(el) {
                        return el.toShortString('') + ' ' + el.getIntegerWeight() + ' ' +  el.linearScore(expSpectrumAr);
                    });
                    var bestProgAr = bestCands.map(function(el) {
                        return el.toShortString('') + ' '  + el.getIntegerWeight() + ' '+ el.score(expSpectrumAr,false);
                    });
                    */
                    progCallback('Process branch [' + candInd + '/' + candidates.length + '] round ', round, '','seqLeaderboardCyclopeptide',[candProgAr,bestProgAr]);
                }
            }

            for (var am = 0; am < allAms.length; ++am) {
                var amStr =  allAms[am];

                var newCand;
                if (includeAll200) {
                    var arr = cand.toStringArray();
                    arr.push(amStr);
                    newCand = new Peptide(Peptide.AminoArrFromArr(arr));
                }
                else {
                    newCand = new Peptide(Peptide.AminoArrFromStr(cand.toShortString('') + amStr));
                }

                var parWeight;

                var newCandW = newCand.getIntegerWeight();
                var parTolerance = 0; //used for float spectrum where parent weight not known
                var parMinTolerance = 0;
                if (isFloatFlag) {
                    parWeight = Math.round(parentWeight - 1);
                    parTolerance = 1;
                    parMinTolerance = 0;
                }
                else {
                    parWeight = parentWeight;
                }

               // if (newCandW > (parentWeight + floatErr)) {
                //if (newCandW > parWeight) {
                if (newCandW > parWeight + parTolerance) {
                   // break; // did assume amino weights are in ascending order, so if one blows, the rest blow. Removed this assumpgtion
                }
              //  else if ((newCandW >= parentWeight - floatErr) && (newCandW <= parentWeight + floatErr)) {
               // else if (newCandW == parWeight) {
                else if ((newCandW >= parWeight + parMinTolerance) && (newCandW <= parWeight + parTolerance)) {
                    if (newCand.score(expSpectrumAr,linearFlag,prefSufOnly) >  highScore) {
                        bestCands = [newCand];
                        highScore = newCand.score(expSpectrumAr,linearFlag,prefSufOnly);
                    }
                    else if (newCand.score(expSpectrumAr,linearFlag,prefSufOnly) == highScore) {
                        bestCands.push(newCand);
                    }
                    newCandidates.push(newCand); //12/4/16 - keep full peptides  on leaderboard!
                }
                else {
                    newCandidates.push(newCand);

                }
            }


        }

        //bound

        newCandidates = trimLeaderboard(newCandidates,expSpectrumAr,N,round,progCallback);


        if (newCandidates.length == 0) {
            done = true;
        }
        else {
            candidates = newCandidates.map(function (el) {
                return el;
            });
        }


    }

/*
    var goodStrings = '';

    var str;
    goodOnes.forEach(function(goodOne) {
        str = '';
        goodOne.forEach(function(weight,weightInd) {
            str+=weight;

            if (weightInd == goodOne.length - 1) {

            }
            else {
                str+='-';

            }


        });
        goodStrings+= str + ' ';

    });


    return goodStrings;
    //done = true;


    //aminoWeights.forEach(function(am) {

    //});
*/
    return bestCands;



}


function trimLeaderboard(leaderBoard,spectrum,N,round,progCallback,prefSufOnly) {
    //trim leaderboard in cyclopeptide sequencing

    //leaderboard is an array of [[peptide,score],[peptide,score]..]

    var progThreshold = 100;


    leaderBoard = leaderBoard.map(function(el,i) {
        if (progCallback) {
            if (((i % progThreshold == 0)) || i  == leaderBoard.length - 1) {
                progCallback('Process bound [' + i + '/' + leaderBoard.length + '] round ', round, '');
            }
        }
       return [el,el.linearScore(spectrum,prefSufOnly)]; 
       // return [el,el.score(spectrum,false)]; //uses cyclic score on Rosalind tedxtbook track problem
    });

    leaderBoard.sort(function(a,b) {
        //return b.linearScore(spectrum) - a.linearScore(spectrum);
        return b[1] - a[1];
    });

    if (leaderBoard.length  <= N) {
        return leaderBoard.map(function(el) {
            return el[0];
        });
    }

    //var limScore = leaderBoard[N - 1].linearScore(spectrum);
    var limScore = leaderBoard[N-1][1];
    var i = N;
   // while ((i < leaderBoard.length) && (leaderBoard[i].linearScore(spectrum) == limScore)) {
    while ((i < leaderBoard.length) && (leaderBoard[i][1] == limScore)) {
        ++i; //add ties
    }

    return leaderBoard.slice(0,i).map(function(el) {
        return el[0];
    });
}

function kmerComposition(dna,k,includePositions) {

    //includePositions = flag to determine whether the positions of each kmer are also returned

    var kmers = [];


    for (var i = 0; i < dna.length - k + 1; ++i) {
        if (includePositions) {
            kmers.push([dna.substring(i, i + k),i]);
        }
        else {
            kmers.push(dna.substring(i, i + k));
        }

    }
    return kmers.sort();

}

function kmerPairedComposition(dna,k,dist) {
    var kmers = [];
    for (var i = 0;i < dna.length - (k*2)  - dist + 1;++i) {
        kmers.push([dna.substring(i,i+k),dna.substring(i+ k + dist,i+ k + dist+k)]);
    }
    return kmers.sort();
}

function hamDist(dna1,dna2,limit) {
    //assumes equal length, otherwise returns -1
    //limit is max num mismatches before not bothering to continue

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

                if (dist >= limit) {
                    return dist;
                }
                        

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
    
    if (!dna) {
        return   [-1,-1,[],[],[]];
    }

    //skewAlready contains skews already processed for part of DNA

    var skewArray = [0];
    var min = 0;
    var max = 0;
    //var minPos = -1;
    //var maxPos = -1;
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

    var i,sum,rem;

    if (typeof(ind) === 'string') {
        //assuming big int

        var curr = BigInteger(ind);
        var pow;

        var inds = '';

        for (i = 0;i < k;++i) {
            pow = k - i - 1;
            var currPow = BigInteger(c_NumBases).pow(BigInteger(pow));
            sum = curr.divide(currPow);
           // var sum = Math.floor(curr / Math.pow(c_NumBases,pow));
            rem = curr.remainder(currPow);
            inds+= c_Bases[sum.valueOf()];
            curr = rem;

        }


    }

    else {
        curr = ind;


        inds = '';

        for (i = 0;i < k;++i) {
            pow = k - i - 1;
            sum = Math.floor(curr / Math.pow(c_NumBases,pow));
            rem = curr % Math.pow(c_NumBases,pow);
            inds+= c_Bases[sum];
            curr = rem;

        }


    }

    return inds;


}

function kMerToInd(kmer) {

    var pow = 0;
    //var sum = 0;
    var vals = {'A':0,'C':1,'G':2,'T':3};


    var big = BigInteger(0);
    for (var i = kmer.length-1;i >= 0; --i) {
        var placeVal = BigInteger(vals[kmer[i]]);
        var bigPow = BigInteger(c_NumBases);
        bigPow = bigPow.pow(BigInteger(pow));
        var curr = placeVal.multiply(bigPow);
        big = big.add(curr);
        //sum+= vals[kmer[i]] * (Math.pow(c_NumBases,pow));
        ++pow;

    }

   // return sum;
    return big.toString();


}

function kMerToIndSmall(kmer) {

    var pow = 0;
    //var sum = 0;
    var vals = {'A':0,'C':1,'G':2,'T':3};


    var summ = 0;
    for (var i = kmer.length-1;i >= 0; --i) {
        summ+= vals[kmer[i]] * (Math.pow(c_NumBases,pow));
        ++pow;
    }

    // return sum;
    return summ;


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
    cols.forEach(function(col) {
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

    c_Bases.forEach(function() {
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

    var profileMatrix = countMatrix.map(function(el) {
        return el.map(function(col) {
            if (laplace) {
                return col * 1.0 / (motifs.length + c_NumBases);
            }
            else {
                return col * 1.0 / motifs.length;
            }

        });

    });

    var entropyMatrix  = profileMatrix.map(function(el) {
        return el.map(function(col) {
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

    var c;
    var r;

    for(c = 0;c < consensus.length;++c) {
        var colEnt = 0;
        for (r = 0;r < entropyMatrix.length; ++r) {
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
    for (c = 0;c < consensus.length;++c) {
        var colRows = [];
        for (r = 0; r < profileMatrix.length;++r) {
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

    return transpose(colMatrix);

   

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

    var r;

    for (var i = 0;i < dna.length - k +1 ; ++i) {
        var prob = 1.0;
        var kmer = dna.substring(i,i+k);
        for (var c = 0;c < kmer.length;++c) {
             r = c_BaseInds[kmer[c]];
            prob*= profMatrix[r][c];
        }
        probDist.push(prob);
   }
   r = getRandomFromProbDist(probDist);

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

function numBreakpointsGenome(g) {

    var n = g.length;

    var newG = g.map(function(el) {
        return el;
    });
    var top = newG.length + 1
    newG.push("+" + top);
    newG.unshift("+0");

    var prev;
    var numAdj = 0;

    newG.forEach(function(el,i) {
        if (i == 0) {

        }
        else {
            if (parseInt(el) == prev + 1) {
                ++numAdj;
            }
        }
        prev = parseInt(el);

    });

    var numBP = n + 1 - numAdj;
    return numBP;

}

function reverseGenome(g,st,end) {
    //inclusive
    gRev = [];
    iNew = 0;

    for (var i = 0;i < st;++i) {
        gRev.push(g[i]);

    }
    for (var i = end;i >=st;--i) {
        var signRev = g[i][0] == '+' ? '-' : '+';
        gRev.push(signRev + g[i].substring(1));
    }
    for (var i = end+1;i < g.length;++i) {
        gRev.push(g[i]);
    }
    return gRev;
}

function findInGenome(g,num,st) {
    for (var i = st;i < g.length;++i) {
        if (parseInt(g[i].substring(1)) == num) {
            return i;
        }
    }
    return -1;

}

function greedyReversal(genome) {

    var bp = numBreakpointsGenome(genome);

    var steps = [];
    for (pos = 0;pos < genome.length - 1;++pos) {
        var sign = genome[pos][0];
        var num = parseInt(genome[pos].substring(1));
        if (pos+1 == num) {
            //correct pos already
            if (sign == '+') {
                //correct sign already
            }
            else {
                var genomeNew = genome.map(function(el) {
                    return el;
                });
                genomeNew[pos] = genomeNew[pos].replace('-','+');
                steps.push(genomeNew);
                genome = genomeNew.map(function(el) {
                    return el;
                });
            }

        }
        else {
            var actualPos = findInGenome(genome,pos + 1,pos+1);
            var genomeNew = reverseGenome(genome,pos,actualPos);
            steps.push(genomeNew);
            genome = genomeNew.map(function(el) {
                return el;
            });

            if (genome[pos][0] == '+') {
                //correct sign already
            }
            else {
                var genomeNew = genome.map(function(el) {
                    return el;
                });
                genomeNew[pos] = genomeNew[pos].replace('-','+');
                steps.push(genomeNew);
                genome = genomeNew.map(function(el) {
                    return el;
                });

            }

        }
    }

    if (genome[genome.length - 1][0] == '+') {
        //correct sign already
    }
    else {
        var genomeNew = genome.map(function(el) {
            return el;
        });
        genomeNew[genomeNew.length - 1] = genomeNew[genomeNew.length - 1].replace('-','+');
        steps.push(genomeNew);
        genome = genomeNew.map(function(el) {
            return el;
        });
    }
    return steps;
    
}

function partSharedKmers(k,tOffset,tNumToProcess,s,t,progressCallback) {

    var sharedAr = [];

    var tDict = {};

    var progThreshold = 1000;

    for (var i = tOffset;(i <  tOffset + tNumToProcess) && (i < t.length - k + 1);++i) {

        if (progressCallback) {
            if ((i  % progThreshold == 0) || (i  == t.length - k)) {
                progressCallback('T Processing', i, t.length - k);
            }
        }


        var kmer = t.substring(i,i+k);
        if (kmer in tDict) {
            tDict[kmer].push(i);
        }
        else {
            tDict[kmer] = [i];
        }
    }
    var soFar = i;

    var numShared = 0;

    for (i = 0;i < s.length - k + 1;++i) {

        if (progressCallback) {
            if ((i  % progThreshold == 0) || (i  == s.length - k)) {
                progressCallback('S Processing [T = ' + soFar + ']', i, s.length - k);
            }
        }

        kmer = s.substring(i, i + k);
        var kmerRev = reverseComplement(kmer);

        if (kmer in tDict) {
            for (var j = 0;j < tDict[kmer].length; ++j) {
                var el = tDict[kmer][j];
                // tDict[kmer].forEach(function (el) {
               // sharedAr.push('(' + i + ', ' + el + ')');
                sharedAr.push([i,el,'C:R']);
                ++numShared;
            }
        }


        if (kmer == kmerRev) {
        }
        else {
            if (kmerRev in tDict) {
                for (var j = 0;j < tDict[kmerRev].length;++j) {
                    var el = tDict[kmerRev][j];
                    //  tDict[kmerRev].forEach(function (el) {
                    //sharedAr.push('(' + i + ', ' + el + ')');
                    sharedAr.push([i,el,'C:B']);
                    ++numShared;
                }
            }
        }

    }

    return sharedAr;



}

function initPlot(c,points) {
    var h = c.height;
    var w = c.width;

    var ctx;

    var params = {};


    var maxX = -999999999;
    var maxY = -999999999;
    var minX = 999999999;
    var minY = 999999999;
    points.forEach(function(el) {
        if (el[0] < minX) {
            minX = el[0];
        }
        if (el[0] > maxX) {
            maxX = el[0];
        }
        if (el[1] < minY) {
            minY = el[1];
        }
        if (el[1] > maxY) {
            maxY= el[1];
        }
    });


    var spanY = maxY - minY;

    var spanX = maxX - minX;

    var unitW;
    var numUnits;

    if (spanX > w) {
        unitW = 1;
        numUnits = w;
    }
    else {
        unitW = Math.floor(w * 1.0 / ( spanX + 2));
        numUnits = spanX;
    }


    var unitH;
    // var numYUnits;

    if (spanY > h) {
        unitH = 1;
        // numYUnits = h;
    }
    else {
        unitH = Math.floor(h * 1.0 / spanY);
        // numYUnits = span;
    }

    //var unitH = Math.floor(h * 1.0 / (max - min));

    var scaleX = 1.0;
    if (spanX  > w) {
        scaleX = Math.ceil(spanX / w);
    }

    var scaleY = 1.0;
    if (spanY > h) {
        scaleY = Math.ceil(spanY / h);
    }

    params.minX = minX;
    params.maxX = maxX;
    params.minY = minY;
    params.maxY = maxY;
    params.spanX = spanX;
    params.spanY = spanY;
    params.unitH = unitH;
    params.unitW = unitW;
    params.scaleX = scaleX;
    params.scaleY = scaleY;

    ctx = c.getContext("2d");
    ctx.clearRect(0,0,w,h);

    return params;


}

function plot(c,points,params) {
    //plot on supplied canvas
   
    var h = c.height;
    var w = c.width;

    var ctx;
    if (points) {
    }
    else {
        ctx = c.getContext("2d");
        ctx.clearRect(0, 0, w, h);
        return;
    }

    ctx = c.getContext("2d");

    if (params) {

    }
    else {
        params = initPlot(c,points);
        ctx.clearRect(0,0,w,h);

    }



    ctx.font = '10pt Arial';
    ctx.fillStyle = ('#FFFFFF');
   // ctx.fillText('G-C Skew',2,15);
   // ctx.fillText('Min: ' + min +  '[' + skewData[3] + ']',80    ,15);
   // ctx.fillText('Max: ' + max +  '[' + skewData[4] + ']',190,15);




    ctx.beginPath();



    ctx.moveTo(0,0);
    ctx.lineTo(0,h);
    var zeroH = params.unitH * (params.maxY - 0) / params.scaleY;
    ctx.moveTo(0,zeroH);//h / 2.0);
    ctx.lineTo(w,zeroH);//h / 2.0);
    ctx.moveTo(0,zeroH);//h / 2.0);


    // skewArray.forEach(function(el,i) {

   // for (var i = 1;i < numUnits ; ++i) {
    for (var n = 0;n < points.length;++n) {


        //  ctx.fillRect(unitW * (i+1),unitH * (max - el),4,4);
        var x = params.unitW * ((points[n][0]  - params.minX) / params.scaleX);
        var y = params.unitH * ( (params.maxY - points[n][1] )/ params.scaleY);

        //ctx.lineTo(x,y);
        //ctx.moveTo(x,y);
        if (points[n][2] == 'C:U') {
            ctx.fillStyle = ('#0000FF');

        }
        else if (points[n][2] == 'C:W') {
            ctx.fillStyle = ('#FFFFFF');
        }
        else if (points[n][2] == 'C:R') {
            ctx.fillStyle = ('#FF0000');
        }
        else if (points[n][2] == 'C:G') {
            ctx.fillStyle = ('#00FF00');
        }
        else if (points[n][2] == 'C:B') {
            ctx.fillStyle = ('#000000');
        }
        else if (points[n][2] == 'C:Y') {
            ctx.fillStyle = ('#FFFF00');
        }
        else if (points[n][2] == 'C:O') {
            ctx.fillStyle = ('#FF4500');
        }
        else {
            ctx.fillStyle = ('#FF0000');
        }

        ctx.fillRect(x,y,2,2); //2,2

    }
    // });
    ctx.stroke();
    ctx.closePath();


}

function lowestCoins(amt,coins) {

    var lowest = coins.reduce(function (a, b) {
        return a < b ? a : b;
    });

    lowArr = [];
    lowArr.push(0);

    for (var i = 1; i < amt + 1; ++i) {
        if (i < lowest) {
            lowArr.push(-1);
        }
        else if (i == lowest) {
            lowArr.push(1);
        }
        else {
            var lowCands = [];
            coins.forEach(function (el) {
                if ((i - el) < 0) {
                    //lowCands.push(-1);
                }
                else {
                    var tmp = i - el;

                    lowCands.push(lowArr[i - el] + 1);

                }


            });
            if (lowCands.length < 2) {
                lowArr.push(lowCands[0])
            }
            else {
                var bestCand = lowCands.reduce(function (a, b) {
                    return a < b ? a : b;

                });
                lowArr.push(bestCand);
            }
            ;

        }


    }

    return lowArr[lowArr.length - 1];

    var lowest = coins.reduce(function (a, b) {
        return a < b ? a : b;
    });

    if (amt < lowest) {
        return 999999;
    }
    if (coins.indexOf(amt) > -1) {
        return 1;
    }

    var lowests = coins.map(function (el) {
        return [el, lowestCoins(amt - el, coins)];
    });

    var best = lowests.reduce(function (a, b) {
            if (a[1] < b[1]) {
                return a;
            }
            else {
                return b;
            }
        }
    );

    return lowestCoins(amt - best[0], coins) + 1;
}

function closestCentre(point,centres) {
      var minDist = 999999999999;
      var minInd = -1;
   centres.forEach(function(c,i) {
      var d = 0;
      if (point[0] == 18.1) {
         console.log(' aha point found');
         
      }
      for (var ind = 0;ind < c.length;++ind) {
         var thisD = (c[ind] - point[ind]) * (c[ind] - point[ind]);
         
         d+= thisD;
      }
      //console.log('d before sqrt: ' + d);
      
      d = Math.sqrt(d);
      
     // console.log('centre ind: ' + i  + ' point: ' + point +  ' dist to point: ' + d );
     
      if (d < minDist) {
         minDist = d;
         minInd = i;
      }
      
      
   });
   
  if (point[0] == 18.1) {
     console.log('mindist: ' + minDist + ' minInd: ' + minInd);
  }
 
   
   return [minDist,minInd];

}

function peptideIdentification(spectralVectorString,proteome) {
	
	      var vecAr = spectralVectorString.split(' ');
		  
	      vecAr = vecAr.map(function(el) {
			 return parseInt(el); 
		  });
		  
     	  var pepWeight = vecAr.length;
		  
		  var candidates = [];
		  
		  		  
		  for (var stPos = 0;stPos <= proteome.length - 1;++stPos) {
			  //find candidate peptide starting at stPos with correct weight
			  var candFound = false;
			  for (var endPos = stPos;endPos <= proteome.length - 1;++endPos) {
				  var pepStr = proteome.substring(stPos,endPos+1);
				  var pep = new Peptide(Peptide.AminoArrFromStr(pepStr));
				  var w = pep.getIntegerWeight();
				  if (w == pepWeight) {
					  candFound = true;
					  candidates.push(pep);
				  }
				  if (w >= pepWeight) {
					  break;
				  }
			  }
		  }
		  
		  
		  var bestCand = '';
		  var bestScore = DGraph.infinity * -1;
		
		  candidates.forEach(function(cand,i) {
			  var score = cand.scoreAgainstSpectralVector(vecAr);
			  if (score > bestScore) {
				  bestCand = cand.toShortString('');
				  bestScore = score;
			  }
			  
		  });
		  
	 
		  
		  return [bestCand,bestScore];
	

	
}

function spectralDictSize(vecAr,t,grid,useToy) {
	
	if (t < 0) {
		return 0;
	}
	var i = vecAr.length;
	
	if ((i == 0) && (t == 0)) {
		return 1;
	}
	
	
	var aminoSet;
	if (useToy) {
	  aminoSet = ['X','Z'];
	}
	else {
	  
	  aminoSet = c_Aminos.slice(0,c_Aminos.length -1);  //c_aminos contains an extra 'X' at end
	}
	
	var sumPrevious = 0;
	
	var currPeak = vecAr[vecAr.length - 1];
	
	aminoSet.forEach(function(el) {
		
		var w;
		if (useToy) {
			w = 0;
			if (el == 'X') {
				w = 4;
			}
			else if (el == 'Z') {
				w = 5;
			}
		}
		else {
			var am = new Amino(el);
		    w = am.getIntegerWeight();
		}
		
		var prevI = i - w;
		if (prevI < 0) {
			
		}
		else {
			var prevT = t - currPeak;
			var key = prevI + ':' + prevT;
			if (key in grid) {
				sumPrevious += grid[key];
			}
			else {
				prev = spectralDictSize(vecAr.slice(0,prevI),prevT,grid,useToy);
				grid[key] = prev;
				sumPrevious += prev;
			}
		}
		
	});
	
	return sumPrevious;
	
	
}