'''
PROGRAM TO MAP CORE REGULATORY CIRCUITRY
VERSION 1.0, December 2015
SOFTWARE AUTHORS: Violaine Saint-Andre, Alexander J. Federation, Charles Y. Lin
REFERENCE: Models of Human Core Transcriptional Regulatory Circuitries.
Violaine Saint-Andre, Alexander J. Federation, Charles Y. Lin,  Brian J. Abraham, Jessica Reddy, Tong Ihn Lee, James E. Bradner, Richard A. Young
CONTACT: violaine.saint-andre@curie.fr
Developed using Python 2.7.3
'''

#==================================================================
#=========================DEPENDENCIES=============================
#==================================================================


import os
import sys
import utils
import string
import numpy
import scipy
import scipy.stats
from string import upper
from subprocess import call
from random import randrange
import networkx as nx
from networkx.algorithms.clique import find_cliques_recursive


#==================================================================
#=========================FUNCTIONS================================
#==================================================================


def calculatePromoterActivity(annotationFile, bamFile, projectName, projectFolder, refseqToNameDict):
    '''
    calculates the level of H3K27ac at each promoter from a H3K27ac bam file
    '''

    print 'IDENTIFY EXPRESSED GENES'

    annotTable = utils.parseTable(annotationFile, '\t')
    output = []
    counter = 0

    bam = utils.Bam(bamFile)

    startDict = utils.makeStartDict(annotationFile)

    tssLoci = []
    for gene in startDict:
        tssLoci.append(utils.makeTSSLocus(gene,startDict,1000,1000))
    tssCollection = utils.LocusCollection(tssLoci,50)

    gff = utils.locusCollectionToGFF(tssCollection)


    outputname = projectFolder + projectName + '_TSS.gff'
    utils.unParseTable(gff, outputname, '\t')

    # run bamToGFF.py to quantify signal at each TSS +/- 1kb

    mappingCmd = 'python ./bamToGFF.py'
    mappingCmd += ' -r '
    mappingCmd += ' -d '
    mappingCmd += ' -o ' + projectFolder + 'matrix.gff'
    mappingCmd += ' -m 1 -f 0 -e 200 '
    mappingCmd += ' -i ' + projectFolder + projectName + '_TSS.gff'
    mappingCmd += ' -b ' + bamFile

    call(mappingCmd, shell=True)

    print  mappingCmd

def createSuperLoci(superTable, Enumber='super'):
    '''
    takes as input a ROSE SuperEnhancer table 
    output a table of loci for SuperEnhancers
    '''

    print 'CREATING SUPER-ENHANCER LOCUS COLLECTION'

    output = []

    if Enumber == 'super':
        for line in superTable[6:]:
            if line[-1] == '1':
                locus = utils.Locus(line[1], line[2], line[3], '.', line[0], (float(line[6])-float(line[7])))
                output.append(locus)
    else:
        end = 6+int(Enumber)
        for line in superTable[6:end]:
            locus = utils.Locus(line[1], line[2], line[3], '.', line[0], (float(line[6])-float(line[7])))
            output.append(locus)

    return output

def createExpressionDict(annotationFile, projectFolder, projectName, refseqToNameDict,expressionTable):
    '''
    takes as input an activity table with refseq NMID in first column and expression or promoter
    acetylation level in a second column
    output a dictionary keyed by refseq containing activity
    '''

    print 'CREATING EXPRESSION DICTIONARY'

    annotTable = utils.parseTable(annotationFile, '\t')
    for line in annotTable:
        gid = line[1]
        genename = upper(line[12])
        refseqToNameDict[gid] = genename

    expresionFilename = projectFolder + 'matrix.gff'
    expressionTable = utils.parseTable(expresionFilename, '\t')

    expressionDictNM = {}
    expressionDictGene = {}

    for line in expressionTable[1:]:
        trid = line[0]
        geneName = refseqToNameDict[trid]
        if len(expressionTable[1]) == 3: #when expressionTable is an output from bamToGFF.py
                exp = float(line[2])
        else: #when expressionTable is passed as an option (2 columns)
                exp = float(line[1])

        # Store the expression value for each NMid in a dict, keep higher value if multiple identical NMIDs
        if trid in expressionDictNM and exp > expressionDictNM[trid]:
            expressionDictNM[trid] = exp
        elif trid not in expressionDictNM:
            expressionDictNM[trid] = exp

        # Store the highest value of transcript expression for each gene
        if geneName in expressionDictGene and exp > expressionDictGene[geneName]:
            expressionDictGene[geneName] = exp
        elif geneName not in expressionDictGene:
            expressionDictGene[geneName] = exp

    # Calculate the cutoff H3K27ac signal value to consider top 2/3 of genes expressed 
    # or the percentile of genes considered expressed passed in option
    cutoff = numpy.percentile(expressionDictGene.values(), 33)
    print 'Expression cutoff: ' + str(cutoff)

    # Select all NMids that are above the computed cutoff
    expressedGenes = []
    expressedNM = []
    for trid in expressionDictNM:
        if float(expressionDictNM[trid]) >= cutoff:
            expressedGenes.append(refseqToNameDict[trid])
            expressedNM.append(trid)
    expressedGenes = utils.uniquify(expressedGenes)

   # Output the list of transcripts considered expressed
    NMfilename = projectFolder + projectName + '_EXPRESSED_TRANSCRIPTS.txt'

   # Output the list of genes considered expressed
    Genefilename = projectFolder + projectName + '_EXPRESSED_GENES.txt'

    utils.unParseTable(expressedNM, NMfilename, '')
    utils.unParseTable(expressedGenes, Genefilename, '')

    return expressedNM

def findCanidateTFs(annotationFile, superLoci, expressedNM, TFlist, refseqToNameDict, projectFolder, projectName):
    '''
    find all TFs within 1Mb of the super-enhancer center that are considered expressed 
    return a dictionary keyed by TF that points to a list of super-enhancer loci
    '''

    print 'FINDING CANIDATE TFs'

    startDict = utils.makeStartDict(annotationFile)

    # Find the location of the TSS of all transcripts (NMid) considered expressed
    tssLoci = []
    for geneID in expressedNM:
        tssLoci.append(utils.makeTSSLocus(geneID,startDict,0,0))
    tssCollection = utils.LocusCollection(tssLoci,50)

    # Assign all transcripts (NMid) that are TFs to a super-enhancer if it is the closest gene
    seAssignment = []
    seAssignmentGene = []
    TFandSuperDict = {}

    for superEnh in superLoci:

        seCenter = (superEnh.start() + superEnh.end()) / 2 

        # Find all transcripts whose TSS occur within 1Mb of the SE center
        searchLocus = utils.Locus(superEnh.chr(), superEnh.start()-1000000, superEnh.end()+1000000, '.')
        allEnhancerLoci = tssCollection.getOverlap(searchLocus)
        allEnhancerGenes = [locus.ID() for locus in allEnhancerLoci]

        # Find the transcript that is closest to the center
        if allEnhancerGenes:
            distList = [abs(seCenter - startDict[geneID]['start'][0]) for geneID in allEnhancerGenes]
            closestGene = allEnhancerGenes[distList.index(min(distList))]
        else:
            closestGene = ''

        seAssignment.append([superEnh.chr(), superEnh.start(), superEnh.end(), closestGene])

        # Select the transcript if it is a TF, and allow for a TF to have multiple SEs
        if closestGene in TFlist and closestGene not in TFandSuperDict.keys():
            TFandSuperDict[closestGene] = [superEnh]
        elif closestGene in TFlist and closestGene in TFandSuperDict.keys():
            TFandSuperDict[closestGene].append(superEnh)

        # Convert the selected TF NMids to gene names
        if closestGene != '':
            geneName = refseqToNameDict[closestGene]
            seAssignmentGene.append([superEnh.chr(), superEnh.start(), superEnh.end(), geneName])

    # Output the list of SE-assigned transcripts (NMids)
    seAssignmentFile = projectFolder + projectName + '_SE_ASSIGNMENT_TRANSCRIPT.txt'
    utils.unParseTable(seAssignment, seAssignmentFile, '\t')

    # Output the list of SE-assigned genes
    seAssignmentGeneFile = projectFolder + projectName + '_SE_ASSIGNMENT_GENE.txt'
    utils.unParseTable(seAssignmentGene, seAssignmentGeneFile, '\t')

    print 'Number of canidate TFs:', len(TFandSuperDict)

    return TFandSuperDict

def formatOutput(TFandSuperDict, refseqToNameDict, projectName, projectFolder):

    '''
    takes as input the dictionary mapping TFs to all proximal super-enhancers
    returns a file that lists each candidate TFs
    and gives the coordinates of the super-enhancers around them
    '''

    print 'CREATE CANDIDATE TFs AND SE TABLE'

    output = [['TF_refseq', 'TF_name', 'chr', 'start', 'stop', 'SuperID', 'Super_Load' ]]

    used = []
 
    for gene in TFandSuperDict.keys():
        for superEnh in TFandSuperDict[gene]:

            check = (refseqToNameDict[gene], superEnh.chr(), superEnh.start(), superEnh.end())

            if check not in used:
                newline = [gene, refseqToNameDict[gene]]
                newline.append(superEnh.chr())
                newline.append(superEnh.start())
                newline.append(superEnh.end())
                newline.append(superEnh.ID())
                newline.append(superEnh.score())
                output.append(newline)

                used.append(check)

    # Output the list of SE-assigned TFs and the associated super-enhancer loci
    outputname = projectFolder + projectName + '_CANIDATE_TF_AND_SUPER_TABLE.txt'

    utils.unParseTable(output, outputname, '\t')

    return 1

def generateSubpeakFASTA(TFandSuperDict, subpeaks, genomeDirectory, projectName, projectFolder, motifExtension):
    '''
    takes as input a BED file of constituents
    outputs a FASTA  file of merged extended super-enhancer consituents and associated formated name
    '''

    print 'MAKE FASTA'

    subpeakDict = {}
    subpeakBED = [['track name=' + projectName + ' color=204,0,204']]
    subpeakTable = utils.parseTable(subpeaks, '\t')

    subpeakLoci = [utils.Locus(l[0], int(l[1]), int(l[2]), '.') for l in subpeakTable]
    subpeakCollection = utils.LocusCollection(subpeakLoci, 50)

    for gene in TFandSuperDict.keys():
        subpeakDict[gene] = []
        for region in TFandSuperDict[gene]:
            overlaps = subpeakCollection.getOverlap(region)
            extendedOverlaps = [utils.makeSearchLocus(x, motifExtension, motifExtension) for x in overlaps]

            overlapCollectionTemp = utils.LocusCollection(extendedOverlaps, 50)
            overlapCollection = overlapCollectionTemp.stitchCollection()
            for overlap in overlapCollection.getLoci():
                subpeakBED.append([overlap.chr(), overlap.start(), overlap.end()])
                subpeakDict[gene].append(overlap)

    bedfilename = projectFolder + projectName + '_subpeaks.bed'
    utils.unParseTable(subpeakBED, bedfilename, '\t')

    fasta = []

    for gene in subpeakDict:
        for subpeak in subpeakDict[gene]:

            fastaTitle = gene + '|'  + subpeak.chr() + '|' + str(subpeak.start()) + '|' + str(subpeak.end())
            fastaLine = utils.fetchSeq(genomeDirectory, subpeak.chr(), int(subpeak.start()+1), int(subpeak.end()+1))

            fasta.append('>' + fastaTitle)
            fasta.append(upper(fastaLine))

    # Output the fasta file of extended SE constituents
    outname = projectFolder + projectName + '_SUBPEAKS.fa'

    utils.unParseTable(fasta, outname, '')

def findMotifs(candidateGenes, projectFolder, projectName, motifConvertFile, motifDatabaseFile):
    '''Run the motif search on the extended SE constituents with FIMO
    '''

    print 'MOTIF SEARCH'

    # Create a dictionary of motif keyed on each TF
    motifDatabase = utils.parseTable(motifConvertFile, '\t')
    motifDatabaseDict = {}
    motifNames = [line[1] for line in motifDatabase]
    for line in motifDatabase:
        motifDatabaseDict[line[1]] = []
    for line in motifDatabase:
        motifDatabaseDict[line[1]].append(line[0])

    canidateMotifs = []
    for gene in candidateGenes:
        if gene in motifNames:
            canidateMotifs.append(gene)

    print 'Number of annotated candidate TFs that have motifs: ' + str(len(canidateMotifs))
    canidateMotifs = sorted(canidateMotifs)

    # Create a backgroud sequence file to use with FIMO
    bgCmd = 'fasta-get-markov -m 1 < ' + projectFolder + projectName + '_SUBPEAKS.fa > ' + projectFolder + projectName + '_bg.meme'
    call(bgCmd, shell=True)

    # Run the motif search with FIMO
    fimoCmd = 'fimo'
    for motif in canidateMotifs:
        for x in motifDatabaseDict[motif]:
            fimoCmd += ' --motif ' + "'%s'" % (str(x))
    fimoCmd += ' -verbosity 1'
    fimoCmd += ' -text'
    fimoCmd += ' -oc ' + projectFolder
    fimoCmd += ' --bgfile ' + projectFolder + projectName + '_bg.meme'
    fimoCmd += ' ' + motifDatabaseFile + ' '
    fimoCmd += projectFolder + projectName + '_SUBPEAKS.fa'
    fimoCmd += ' > '+ projectFolder + 'fimo.txt'
    print fimoCmd

    fimoOutput = call(fimoCmd, shell=True)

    return fimoCmd

def buildNetwork(projectFolder, projectName, candidateGenes, refseqToNameDict, motifConvertFile):
    '''takes as input the FIMO output file
    identify TF-TF interactions, define candidate TFs as nodes and draw all edges
    '''

    print 'IDENTIFY TF-TF INTERACTIONS'

    motifDatabase = utils.parseTable(motifConvertFile, '\t')
    motifDatabaseDict = {}
    motifNames = [line[1] for line in motifDatabase]
    for line in motifDatabase:
        motifDatabaseDict[line[0]] = line[1]

    fimoFile =  projectFolder + 'fimo.txt'
    fimoTable = utils.parseTable(fimoFile, '\t')

    graph = nx.DiGraph(name=projectName)
    graph.add_nodes_from(candidateGenes)

    motifDictSE = {}

    for gene in candidateGenes:
        motifDictSE[gene] = []

    edgeCountDictSE = {}

    for line in fimoTable[1:]:

        source = motifDatabaseDict[line[0]]
        # line[1] changed to line[2] to adapt to the output of the new version of fimo
        region = line[2].split('|')      
        target = refseqToNameDict[region[0]]
        location = (region[1], int(region[2]), int(region[3]))


        # Count the number of motifs in SEs

        # Initialize the dictionary
        if (source, target) not in edgeCountDictSE.keys():
            edgeCountDictSE[(source,target)] = 0

        # Count unique motifs
        # line[2] changed to line[3] and line[3] changed to line[4] to adapt to the output of the new version of fimo
        if (region[1], int(region[2]) + int(line[3]), int(region[2]) + int(line[4])) not in motifDictSE[source]:  
            edgeCountDictSE[(source, target)] += 1
            motifDictSE[source].append((region[1], int(region[2]) + int(line[3]), int(region[2]) + int(line[4])))

    # Draw an edge if there are at least 3 motif instances in the sum of the merged extended SE constituents
    for connection in edgeCountDictSE.keys():
        if edgeCountDictSE[connection] > 2:
            graph.add_edge(connection[0], connection[1])

    # Output a bedfile of motif locations for each candidate TF
    for gene in motifDictSE.keys():
        if motifDictSE[gene]:
            bed = []
            for loc in motifDictSE[gene]:
                bed.append([loc[0], loc[1], loc[2]])
            filename = projectFolder + gene + '_' + projectName + '_motifs.bed'
            utils.unParseTable(bed, filename, '\t')

    return graph

def formatNetworkOutput(graph, projectFolder, projectName, candidateGenes):
    '''
    takes as input the TF-TF interactions
    Outputs all possible CRCs
    '''

    print 'IDENTIFYING CRCs'

    # Output the list of autoregulated TFs
    autoreg = graph.selfloop_edges()
    selfLoops = [x for x,y in autoreg]
    selfLoopFile = projectFolder + projectName + '_AUTOREG.txt'
    utils.unParseTable(selfLoops, selfLoopFile, '')

    # Recover all bidirectional edges and create a file of TF-TF interactions
    pairs = []
    for n in selfLoops:
        for m in selfLoops:
            if n != m:
                if graph.has_edge(n,m) and graph.has_edge(m,n):
                    pairs.append([n,m])

    #fill up the graph
    G=nx.Graph()
    G.add_nodes_from(selfLoops)
    G.add_edges_from(pairs)
    cliques = find_cliques_recursive(G)
    cliqueList = list(cliques)

    print 'Number of possible CRCs:'
    print len(cliqueList)

    #Score the CRCs

    #count the occurences of the TFs accross the loops
    dicoTFinloopsCounts={}
    for clique in cliqueList:
        for TF in clique:

            if dicoTFinloopsCounts.has_key(TF):
                dicoTFinloopsCounts[TF]+=1

            else:
                dicoTFinloopsCounts[TF]=1

    #calculate a score by CRC
    cliqueRanking = []

    for clique in cliqueList:
        cliqueScore=0

        for TF in clique:
            cliqueScore = (float(cliqueScore) + (float(dicoTFinloopsCounts[TF])))
        cliqueRanking.append((clique, cliqueScore/len(clique), len(clique)))

    # Output a file containing all possible ranked CRCs
    sortCliqueRanking = sorted(cliqueRanking, reverse=True, key=lambda x:x[1])
    cliqueFile = projectFolder + projectName + '_CRC_SCORES.txt'
    utils.unParseTable(sortCliqueRanking, cliqueFile, '\t')

   # Print the top CRC to the standard output
    print 'Top CRC:'
    print sortCliqueRanking[0]


#==================================================================
#=========================MAIN=====================================
#==================================================================

def main():

    from optparse import OptionParser

    usage = "usage: %prog [options] -e [ENHANCER_FILE] -b [BAM_FILE] -g [GENOME] -o [OUTPUTFOLDER] -n [NAME] -s [SUBPEAKS] -x [EXP_CUTOFF] -l [EXTENSION_LENGTH]"
    parser = OptionParser(usage = usage)

    # Required flags                                                                                                                                        
    parser.add_option("-e","--enhancer_file", dest="enhancers",nargs = 1, default=None,
                      help = "Provide a ROSE generated enhancer table (_AllEnhancers.table.txt)")
    parser.add_option("-b","--bam_file",dest="bam",nargs =1, default = None,
                      help = "Provide a sorted indexed bam file for H3K27ac sequencing reads")
    parser.add_option("-g","--genome",dest="genome",nargs =1, default = None,
                      help = "Provide the build of the genome to be used for the analysis. Currently supports HG19, HG18 and MM9")
    parser.add_option("-f","--fasta",dest="fasta",nargs =1, default = None,
                      help = "Enter location of the fasta files for the genome version used")
    parser.add_option("-s","--subpeaks", dest="subpeaks",nargs=1,default=None,
                      help = "Enter a bedfile of peaks output from MACS used to identify SE constituents")
    parser.add_option("-x","--exp_Cutoff", dest="expCutoff",nargs=1,default=33,
                      help = "Enter the percentage of transcripts that are not considered expressed, default=33")
    parser.add_option("-l","--extension_length", dest="extension",nargs = 1, default = 500,
                      help = "Enter the length (in bp) to extend constituents for motif search, default=500")
    parser.add_option("-n","--name",dest="name",nargs =1, default = None,
                      help = "Enter the sample name")
    parser.add_option("-o","--output",dest="output",nargs =1, default = None,
                      help = "Enter directory to be used for storing output")

    # Options                                                                                                               
    parser.add_option("-a","--activity", dest="activity",nargs = 1, default=None,
                      help = "Enter a two column table with refseq in the first column and the associated activity (expression or promoter acetylation level) in the second column")
    parser.add_option("-E","--enhancer_number", dest="Enumber",nargs = 1, default='supers',
                      help = "Enter the number of top ranked enhancers to include in the anlaysis, default = supers")

    (options,args) = parser.parse_args()

    print(options)


    if options.enhancers and options.bam and options.genome and options.fasta and options.subpeaks and options.expCutoff and options.extension and options.name and options.output:

        # Set parameters

        genomeDirectory = options.fasta

        genome = options.genome
        genome = upper(genome)

        if genome == 'HG19':
            annotationFile = './annotation/hg19_refseq.ucsc'
            TFfile = './TFlist_NMid_hg.txt'

        if genome == 'HG18':
            annotationFile = './annotation/hg18_refseq.ucsc'
            TFfile = './TFlist_NMid_hg.txt'

        if genome == 'MM9':
            annotationFile = './annotation/mm9_refseq.ucsc'
            TFfile = './TFlist_NMid_ms.txt'

        motifConvertFile = './MotifDictionary.txt'
        motifDatabaseFile = './VertebratePWMs.txt'

        TFtable = utils.parseTable(TFfile, '\t')
        TFlist = [line[0] for line in TFtable]
        TFlistGene = [line[1] for line in TFtable]

        superFile = options.enhancers
        superTable = utils.parseTable(superFile, '\t')

        bamFile = options.bam
        bam = utils.Bam(bamFile)

        subpeaks = options.subpeaks

        expCutoff = int(options.expCutoff)

        motifExtension = int(options.extension)

        projectName = options.name

        projectFolder = options.output

        refseqToNameDict = {}
        expressionFile = options.activity
        if expressionFile:
            expressionTable = utils.parseTable(expressionFile, '\t')
        else:
            calculatePromoterActivity(annotationFile, bamFile, projectName, projectFolder, refseqToNameDict)
            expresionFilename = projectFolder + 'matrix.gff'
            expressionTable = utils.parseTable(expresionFilename, '\t')
        if options.Enumber != 'super':
            enhancerNumber = options.Enumber
        else:
            enhancerNumber = 'super'

        # Run the program

        superLoci = createSuperLoci(superTable)

        expressedNM = createExpressionDict(annotationFile, projectFolder, projectName, refseqToNameDict, expressionTable)

        TFandSuperDict = findCanidateTFs(annotationFile, superLoci, expressedNM, TFlist, refseqToNameDict, projectFolder, projectName)

        formatOutput(TFandSuperDict, refseqToNameDict, projectName, projectFolder)

        candidateGenes = [upper(refseqToNameDict[x]) for x in TFandSuperDict.keys()]

        candidateGenes = utils.uniquify(candidateGenes)

        generateSubpeakFASTA(TFandSuperDict, subpeaks, genomeDirectory, projectName, projectFolder, motifExtension)

        findMotifs(candidateGenes, projectFolder, projectName, motifConvertFile, motifDatabaseFile)

        graph = buildNetwork(projectFolder, projectName, candidateGenes, refseqToNameDict, motifConvertFile)

        formatNetworkOutput(graph, projectFolder, projectName, candidateGenes)

    # Return help

    else:
        parser.print_help()
        sys.exit()


if __name__ == '__main__':
    main()

