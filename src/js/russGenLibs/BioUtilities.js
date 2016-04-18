/**
 * Created by RussellM on 20/01/2016.
 */


//Bioinformatics routines

var c_Bases = ['A','C','G','T'];
var c_NumBases = 4;
var c_BaseInds = {'A':0,'C':1,'G':2,'T':3};



//Directed Graph routines

function DEdge(sourceNode,targetNode,artificial) {
    this.sourceNode = sourceNode;
    this.targetNode = targetNode;
    this.visited = false;

    this.walkNum = 0;

    this.artificial = artificial; // indicates added to make cycle

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

    this.successors = [];
    this.predecessors = [];

    if (repeatNum) {
        this.repeatNum = repeatNum; //for nodes with same label
    }
    else {
        this.repeatNum = 1;
    }


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

//
function DGraph(source,sourceType,graphType,k,makeCycle,pairDist) {

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

    /*
    if (this.pairDist) {
        this.dna = this.dna + this.dna.substring(0,this.k + this.pairDist - 1); // extend so pairs work up to end
    }*/

    this.nodes = [];

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


        //var svdThis = this;

        var contigs = [];

        var isolatedCycleUsedNodes = [];

        this.nodes.forEach(function(node) {
            if (  (node.InDegree() == 1) && (node.outDegree() == 1)) {

                if (isolatedCycleUsedNodes.indexOf(node)==-1) {
                    //check for isolated cycle
                    var startNode = node;
                    var cyc = [];
                    cyc.push(startNode);
                    var curNode = startNode.successors[0].targetNode;


                    var isolatedCycleFound = false;
                    while ((curNode.InDegree() == 1) && (curNode.outDegree() == 1)) {
                        if (curNode == startNode) {
                            cyc.push(curNode);
                            isolatedCycleFound = true;
                            break;
                        }
                        else {
                            cyc.push(curNode);
                            curNode = curNode.successors[0].targetNode;
                        }
                    }
                    if (isolatedCycleFound) {
                        var cycPath = [];
                        cyc.forEach(function (el) {
                            cycPath.push(el.dna);
                            isolatedCycleUsedNodes.push(el);
                        });
                        contigs.push(cycPath);
                    }
                }

            }
            else {

                node.successors.forEach(function(suf) {
                    var curNode;
                    var contig = [];
                    contig.push(suf.sourceNode.dna);
                    contig.push(suf.targetNode.dna);
                    curNode = suf.targetNode;
                    while ((curNode.outDegree() == 1) && (curNode.InDegree() == 1)) {
                        contig.push(curNode.successors[0].targetNode.dna);
                        curNode = curNode.successors[0].targetNode;
                    }
                    contigs.push(contig);

                });
            }

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

    this.edgePathReconstructedPairs = function() {

        //var svdThis = this;

        var origStr = '';
        var pairStr = '';

        this.edgePath.forEach(function(ed,i) {
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

        var totLen = origStr.length + this.k + this.pairDist;

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

    this.edgePathReconstructed = function() {

        if ((this.sourceType == DGraph.fromPairedReads) || (this.sourceType == DGraph.fromPairedDna)) {
            return this.edgePathReconstructedPairs();
        }

        var str = '';
        for (var i=0;i < this.edgePath.length; ++i) {
            if (i == 0) {
                str += this.edgePath[i].edgeLabel();
            }
           // else if ((this.makeCycle) && (i >= this.edgePath.length - this.k + 1)) {
            else if ((this.isBalanced()) && (i >= this.edgePath.length - this.k + 1)) {


            }
            else   {
                //str += this.edgePath[i].edgeLabel().substring(this.edgePath[i].edgeLabel().length - 1);
                str+= this.edgePath[i].targetNode.dna.substring(this.edgePath[i].targetNode.dna.length - 1);
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
                if (r == -1) {
                    return false;
                }
                else {
                    return true;
                }
        });


    };


    this.makeCyclical = function() {
       // var f = this.definiteFirst();
       // var l = this.definiteLast();

       // if (f || l) {

        if (this.pairedDist) {
            return; //already added to dna to make enough
        }

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

        for (var i = 0;i < this.nodes.length;++i ) {
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
                svdThis = this;
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

                if (svdThis.pairDist) {
                    pref = [kmer[0].substring(0, svdThis.k - 1), kmer[1].substring(0, svdThis.k - 1)];
                    suf = [kmer[0].substring(1, svdThis.k), kmer[1].substring(1, svdThis.k)];
                }
                else {
                    pref = kmer.substring(0, svdThis.k - 1);
                    suf = kmer.substring(1, svdThis.k);

                }

                var prefNode = svdThis.findNode(pref);
                if (!prefNode) {
                    prefNode = new DNode(pref);
                    svdThis.nodes.push(prefNode);
                }

                var sufNode = svdThis.findNode(suf);


                if (!sufNode) {
                    sufNode = new DNode(suf);
                    svdThis.nodes.push(sufNode);
                }

                prefNode.successors.push(new DEdge(prefNode, sufNode));
                sufNode.predecessors.push(prefNode);
            });

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

        var balBefore = this.isBalanced();
        var balLinearBefore = this.isBalancedLinear();

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

DGraph.hamGraph = 1;
DGraph.debGraph = 2;

DGraph.seqTypeComp = 1;
DGraph.seqTypePath = 2;
DGraph.seqTypeCycle = 3;


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
            for (i = pepStrAr.length - 1;i >=0;--i) {
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

        var arr = this.peptide.map(function (am) {
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


    
    
    this.spectrum = function (cyclic,retSubPeptides) {
        
        //retSubPeptides: flag to indicate whether to also return the sub peptide strings

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
                spec.push([w,pepStr.substring(st,end)]);
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

    this.linearSpectrum = function () {
        return this.spectrum(false);
    };

    this.cyclicSpectrum = function () {
        return this.spectrum(true);
    };


    this.linearScore = function (expSpectrum) {

        return this.score(expSpectrum, true);

    };

    this.score = function (expSpectrum, linear) {

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
            theor = this.linearSpectrum();
        }
        else {
            theor = this.cyclicSpectrum();
        }


        var exp = expSpectrum.map(function (el) {
            return el;
        });

        if (isFloatSpectrum) {
            for (var tm = 0;tm < theor.length; ++tm) {
                theorMass = theor[tm];

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

    var amAr = [];

    amAr = arr.map(function(el) {
        var am = new Amino(el);
        return am;
    });
  
    return amAr;


};

Peptide.AminoArrFromWeights = function(weights) {

    var amArr = weights.map(function(w) {
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

Peptide.linear = 1;
Peptide.circular = 2;

Peptide.PepMethodIdealSpectrum = 1;
Peptide.PepMethodSequenceBrute = 2;
Peptide.PepMethodSequenceLeaderboard = 3;
Peptide.PepMethodSequenceLeaderboardConv = 4;


function RNA(rna) {

    this.rna = rna;

    this.init = function() {

    };

    this.codons = function(frameOffset,rev) {
        if (!frameOffset) {
            frameOffset = 0;
        }
        var arr = [];

        if (rev) {
            var strRev = this.rna.split('').reverse().join('');
            for (var i = frameOffset; i < this.rna.length - Codon.len + 1; i = i + Codon.len) {
                var cod = new Codon(strRev.substring(i, i + Codon.len));
                arr.push(cod);

            }

        }
        else {
            for (var i = frameOffset; i < this.rna.length - Codon.len + 1; i = i + Codon.len) {
                var cod = new Codon(this.rna.substring(i, i + Codon.len));
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
        var weightsMatchedThisLen = 0;

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

function cyclopeptideSequencing(expSpectrum,cyclicFlag,progCallback) {
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
                amW = allAminoWeights[am];
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
                        var candSpec =  candPep.spectrum(cyclicFlag); ///candPep.cyclicSpectrum();
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

function leaderboardCyclopeptideSequencing(expSpectrum,M,N,includeAll200,cyclicFlag,useTheseAminos,progCallback) {

    var linearFlag = !cyclicFlag;

    var isFloatFlag;
    var floatErr = 0;

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
            if (topEls.indexOf(w) > -1) {
                return true;
            }
            else {
                return false;
            }

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
            return el.toShortString('') + ' ' + el.getIntegerWeight() + ' ' +  el.linearScore(expSpectrumAr);
        });
        var bestProgAr = bestCands.map(function(el) {
            return el.toShortString('') + ' '  + el.getIntegerWeight() + ' '+ el.score(expSpectrumAr,false);
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
                    if (newCand.score(expSpectrumAr,linearFlag) >  highScore) {
                        bestCands = [newCand];
                        highScore = newCand.score(expSpectrumAr,linearFlag);
                    }
                    else if (newCand.score(expSpectrumAr,linearFlag) == highScore) {
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


function trimLeaderboard(leaderBoard,spectrum,N,round,progCallback) {
    //trim leaderboard in cyclopeptide sequencing

    //leaderboard is an array of [[peptide,score],[peptide,score]..]

    var progThreshold = 100;


    leaderBoard = leaderBoard.map(function(el,i) {
        if (progCallback) {
            if (((i % progThreshold == 0)) || i  == leaderBoard.length - 1) {
                progCallback('Process bound [' + i + '/' + leaderBoard.length + '] round ', round, '');
            }
        }
       return [el,el.linearScore(spectrum)]; 
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

    kmers = [];


    for (i = 0; i < dna.length - k + 1; ++i) {
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
    for (i = 0;i < dna.length - (k*2)  - dist + 1;++i) {
        kmers.push([dna.substring(i,i+k),dna.substring(i+ k + dist,i+ k + dist+k)]);
    }
    return kmers.sort();
}

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

    for(c = 0;c < consensus.length;++c) {
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
    for (c = 0;c < consensus.length;++c) {
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