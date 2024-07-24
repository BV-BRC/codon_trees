import sys
import re
from math import sqrt
import os
import os.path
import glob
import argparse
import subprocess
import json
try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3
from time import time, localtime, strftime                        
import phylocode
import patric_api
from Bio import SeqIO

def genomeIdFromFigId(figId):
    m = re.match("fig\|(\d+\.\d+)", figId)
    if m:
        return m.group(1)
    return None

parser = argparse.ArgumentParser(description="Codon-oriented aligment and tree analysis of PATRIC protein families", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--parametersJson", metavar="file.json", type=str, help="parameters in json format (command line overrides)")
parser.add_argument("--outputBase", metavar="filebase", type=str, help="base name for output files, def=codontree")
parser.add_argument("--outputDirectory", type=str, metavar="out_dir", help="for output, create if needed")
parser.add_argument("--genomeIdsFile", metavar="file", type=str, nargs="*", help="file with PATRIC genome IDs, one per line (or first column of TSV)")
parser.add_argument("--genomeGroup", metavar="name", type=str, nargs="*", help="name of user's genome group at PATRIC")
parser.add_argument("--genomeObjectFile", metavar="file", type=str, help="genome object (json file)")
parser.add_argument("--homologIdsFile", metavar="file", type=str, help="use these homologs (eg PGFams), instead of searching for single-copy ones")
parser.add_argument("--genomePgfamGeneFile", metavar="file", type=str, help="read geneIDs per PGFam per genome from this file")
parser.add_argument("--optionalGenomeIdsFile", metavar="file", type=str, help="optional genome ids, one per line (or first column of TSV)")
parser.add_argument("--homologScope", metavar="global/local", choices=['global', 'local'], default='global', help="use PGFams (global) or PLFams (local}")
parser.add_argument("--maxGenes", metavar="#", type=int, default=50, help="number of genes in concatenated alignment")
parser.add_argument("--excessGenesProp", metavar="prop", type=float, default=0.5, help="multiplier of maxGenes to add to filter out low-scoring alignments")
parser.add_argument("--excessGenesFixed", metavar="#", type=int, default=20, help="fixed excess genes to add to filter out low-scoring alignments")
parser.add_argument("--bootstrapReps", metavar="#", type=int, help="number of raxml 'fast boostrap' replicates")
parser.add_argument("--maxGenomesMissing", metavar="#", type=int, default=0, help="genomes allowed to lack a member of any homolog group")
parser.add_argument("--maxAllowedDups", metavar="#", type=int, default=0, help="duplicated gene occurrences allowed within homolog group")

parser.add_argument("--aligner", metavar="program", type=str, choices=('muscle', 'mafft'), default="mafft", help="program to align protein sequences")
parser.add_argument("--endGapTrimThreshold", metavar="maxPropGaps", type=float, default=0.5, help="stringency of end-gap trimming, 0-1.0, lower for less trimming")
parser.add_argument("--raxmlExecutable", metavar="program_name", type=str, default="raxml", help="program to call, possibly with path")
parser.add_argument("--rateModel", metavar="model", type=str, choices = ['CAT', 'GAMMA'], default="CAT", help="variable rate category model CAT|GAMMA")
parser.add_argument("--proteinModel", metavar="model", type=str, default="AUTO", help="raxml protein substitution model")
parser.add_argument("--protein_sample_prop", metavar="0.XX", type=float, default=0.5, help="proportion of positions to sample for discovery of optimal model")
parser.add_argument("--analyzeCodons", action='store_true', help="analyze only codon nucleotides")
parser.add_argument("--analyzeProteins", action='store_true', help="analyze only amino acids")
parser.add_argument("--threads", "-t", metavar="T", type=int, default=2, help="threads for raxml")
parser.add_argument("--deferRaxml", action='store_true', help="does not run raxml")
parser.add_argument("--exitBeforeAlignment", action='store_true', help="Don't align, useful for getting PGFam matrix")
parser.add_argument("--writePgfamAlignments", action='store_true', help="write fasta alignment per homolog used for tree")
parser.add_argument("--writePgfamAlignmentsDNA", action='store_true', help="write fasta alignment per homolog used for tree")
parser.add_argument("--writePgfamMatrix", type=float, help="write table of gene_id per homolog per genome (present in > proportion of genomes)")
parser.add_argument("--writePgfamCountMatrix", type=float, help="write table of counts per homolog per genome (present in > proportion of genomes)")
parser.add_argument("--writePhyloxml", action='store_true', help="write tree in phyloxml format")
parser.add_argument("--phyloxmlFields", type=str, metavar='data fields', help="comma-sparated genome fields for phyloxml")
parser.add_argument("--pathToFigtreeJar", type=str, metavar="path", help="not needed if figtree executable on path")
parser.add_argument("--universalRolesFile", type=str, metavar="path", help="path to file with universal roles to select conserved genes")
parser.add_argument("--focusGenome", metavar="id", type=str, help="to be highlighted in color in Figtree")
parser.add_argument("--debugMode", action='store_true', help="more output to log file")
parser.add_argument("--authToken", metavar="STRING", type=str, help="patric authentication token")
parser.add_argument("--ignoreAuthEnv", action='store_true', help="turn off authorization by environmental variable")
parser.add_argument("--ignoreAuthRC", action='store_true', help="turn off authorization by file")
#parser.add_argument("--enableGenomeGenePgfamFileReuse", action='store_true', help="read genes and homologs from stored file matching genomeIdsFile if it exists")
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
starttime = time()
genomeIds = set() # list of genome IDs for tree building (enforcing maxGenomesMissing and maxAllowedDupes)
optionalGenomeIds = set()

if args.parametersJson:
    # read parameters in json file, avoiding overwriting any specified on command line
    params = json.load(open(args.parametersJson))
    if "output_path" in params and not args.outputDirectory:
        args.outputDirectory = os.path.abspath(params["output_path"])
    if "output_file" in params and not args.outputBase:
        args.outputBase = params["output_file"]
    if "genome_group" in params: #not necessary if UI flattens these out to ids
        if not args.genomeGroup:
            args.genomeGroup = []
        args.genomeGroup.extend(params["genome_group"])
    if "genome_ids" in params: 
        genomeIds |= set(params["genome_ids"])
    if "optional_genome_ids" in params:
        optionalGenomeIds |= set(params["optional_genome_ids"])
    if "number_of_genes" in params:
        args.maxGenes = params["number_of_genes"]
    if "bootstraps" in params:
        if params["bootstraps"]:
            args.bootstrapReps = 100
    if "max_genomes_missing" in params:
        args.maxGenomesMissing = params["max_genomes_missing"]
    if "max_allowed_dups" in params:
        args.maxAllowedDups = params["max_allowed_dups"]
    
if not args.outputDirectory:
    args.outputDirectory="./"
if not args.outputDirectory.endswith("/"):
    args.outputDirectory += "/"
if not os.path.exists(args.outputDirectory):
    os.makedirs(args.outputDirectory)
if not args.outputBase:
    args.outputBase = "codontree"

logfileName = os.path.basename(sys.argv[0])
logfileName = re.sub(r"\..*", "", logfileName)
logfileName += ".log"
logfileName = os.path.join(args.outputDirectory, logfileName)

global LOG 
#LOG = open(logfileName, 'w')
LOG = sys.stderr
LOG.write("starting %s\n"%sys.argv[0])
LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(starttime))+"\n")
LOG.write("args= "+str(args)+"\n\n")
LOG.flush()
phylocode.LOG = LOG
patric_api.LOG = LOG
patric_api.PatricUser = None
if args.debugMode:
    patric_api.setDebug(True)
    phylocode.Debug = True

def authenticate():
    if args.authToken:
        if args.debugMode:
            LOG.write("authorizing by passed token: %s\n"%args.authToken)
        patric_api.authenticateByString(args.authToken)
    if not patric_api.PatricUser and not args.ignoreAuthEnv:
        if "KB_AUTH_TOKEN" in os.environ:
            if args.debugMode:
                LOG.write("reading auth key from environment\n")
            patric_api.authenticateByString(os.environ.get('KB_AUTH_TOKEN'))
        elif args.debugMode:
            LOG.write("environment variable for authorization not found\n")
    if not patric_api.PatricUser and not args.ignoreAuthRC:
        tokenFile = os.path.join(os.environ.get('HOME'), ".patric_token")
        if os.path.exists(tokenFile):
            LOG.write("reading auth key from file %s\n"%tokenFile)
            with open(tokenFile) as T:
                tokenString = T.read().rstrip()
                patric_api.authenticateByString(tokenString)
        elif args.debugMode:
            LOG.write("authorization file %s not found\n"%tokenFile)
    LOG.write("Patric User = %s\n"%patric_api.PatricUser)

environmentTests = []
environmentTests.append(("phylocode.which(args.raxmlExecutable) != None", phylocode.which(args.raxmlExecutable) != None))
environmentTests.append(("phylocode.which({}) != None".format(args.aligner), phylocode.which(args.aligner) != None))

authenticate()

# update genomeIds set
if args.genomeIdsFile:
    for filename in args.genomeIdsFile:
        with open(filename) as F:
            for line in F:
                m = re.match(r"(\d+\.\d+)", line)
                if m:
                    genomeIds.add(m.group(1))
    LOG.write("from %s got %d genome IDs\n%s\n"%(",".join(args.genomeIdsFile), len(genomeIds), "\t".join(list(genomeIds))))
if args.genomeObjectFile:
    LOG.write("genome object file: %s\n"%args.genomeObjectFile)
LOG.flush()

groupsPerGenome = None
if args.genomeGroup:
    groupsPerGenome = {}
    if not patric_api.PatricUser:
        raise Exception("No patric user is defined, required for using genome group.\n")
    for groupPath in args.genomeGroup:
        LOG.write("requesting genome IDs for user group %s\n"%groupPath)
        groupIds = patric_api.getGenomeGroupIds(groupPath)
        groupName = groupPath.split("/")[-1] # file name with path removed (basename)
        for genomeId in groupIds:
            if genomeId not in groupsPerGenome:
                groupsPerGenome[genomeId] = []
            groupsPerGenome[genomeId].append(groupName)
        genomeIds.update(set(groupIds))
        LOG.write("got %d ids for group %s, total IDs = %d\n"%(len(groupIds), groupName, len(genomeIds)))

if args.optionalGenomeIdsFile:
    with open(args.optionalGenomeIdsFile) as F:
        for line in F:
            m = re.match(r"(\d+\.\d+)", line)
            if m:
                optionalGenomeIds.add(m.group(1))
    LOG.write("got %d optionalGenomeIds\n%s\n"%(len(optionalGenomeIds), "\t".join(optionalGenomeIds)))
    LOG.flush()

# if either codons or proteins is specified, analyze just that, otherwise analyze both
if (args.analyzeCodons or args.analyzeProteins):
    if args.analyzeCodons:
        LOG.write("analyzing just codon nucleotides, not proteins\n")
    else:
        LOG.write("analyzing just proteins, not nucleotides\n")
else:
    args.analyzeCodons = args.analyzeProteins = True
    LOG.write("analyzing both codons and proteins\n")
LOG.flush()


# this is where we gather the list of Pgfam genes for each genome ID
homologMatrix = {}
if args.genomePgfamGeneFile: # reserve for future use (convenient for debugging)
    if not os.path.exists(args.genomePgfamGeneFile):
        raise Exception("genomePgfamGeneFile specified as %s does not exist"%args.genomePgfamGeneFile)
    with open(args.genomePgfamGeneFile) as F:
        homologMatrix = patric_api.read_homolog_gene_matrix(F)

genomeObject = None
genomeObject_genomeId = None
genomeObject_name = None
genomeObjectProteins = None
genomeObjectGeneDna = None
if args.genomeObjectFile:
    genomeObject = json.load(open(args.genomeObjectFile))
    genomeObjectProteins = patric_api.getGenomeObjectProteins(genomeObject)
    genomeObjectGeneDna = patric_api.getGenomeObjectGeneDna(genomeObject)
    genomeObjectGeneHomologs = patric_api.get_homologs_from_genome_object(genomeObject)
    for row in genomeObjectGeneHomologs:
        genomeId, gene, homolog = row
        if homolog not in homologMatrix:
            homologMatrix[homolog] = {}
        if genomeId not in homologMatrix[homolog]:
            homologMatrix[homolog][genomeId] = []
        homologMatrix[homolog][genomeId].append(gene)
    
    genomeObject_genomeId = genomeObject['id']
    genomeObject_name = genomeObject['scientific_name']
    genomeIds.add(genomeObject_genomeId)
    args.focusGenome = genomeObject_genomeId
    LOG.write("parsed json file %s"%args.genomeObjectFile)
    LOG.flush()

allGenomeIds = genomeIds
allGenomeIds.update(optionalGenomeIds)
genomesWithData = set()
for homolog in homologMatrix:
    genomesWithData.update(set(homologMatrix[homolog].keys()))
genomesWithoutData = allGenomeIds - genomesWithData
if len(genomesWithoutData):
    if args.homologIdsFile:
        homologList = []
        with open(args.homologIdsFile) as F:
            for line in F:
                m = re.match("\s*(\S+).*", line)
                if m:
                    homologList.append(m.group(1))
        LOG.write("use homologs listed in {}\n".format(args.homologIdsFile))
        homologMatrix = patric_api.getHomologGenomeMatrix(genomesWithoutData, homologList, homologMatrix)
    else:
        LOG.write("getPgfamGenomeMatrix()\n")
        homologMatrix = patric_api.get_homolog_gene_matrix(genomesWithoutData, homologMatrix, scope=args.homologScope)

if args.writePgfamMatrix:
    with open(os.path.join(args.outputDirectory, args.outputBase+".homologMatrix.txt"), 'w') as F:
        patric_api.write_homolog_gene_matrix(homologMatrix, F, args.writePgfamMatrix)

if args.writePgfamCountMatrix:
    with open(os.path.join(args.outputDirectory, args.outputBase+".homologCountMatrix.txt"), 'w') as F:
        patric_api.write_homolog_count_matrix(homologMatrix, F, args.writePgfamCountMatrix)

if args.exitBeforeAlignment:
    sys.exit(0)

genesPerGenome = {}
for genome in genomeIds:
    genesPerGenome[genome] = 0
for homolog in homologMatrix:
    for genome in homologMatrix[homolog]:
        genesPerGenome[genome] += 1

deletedGenomes = set()
for genome in genomeIds:
    if genesPerGenome[genome] == 0:
        deletedGenomes.add(genome)
        LOG.write("Deleting genome %s for lack of data.\n"%genome)
num_requested_genomes = len(genomeIds)
genomeIds -= deletedGenomes

LOG.write("got homologs for genomes, len=%d\n"%len(homologMatrix))
LOG.flush()

environmentTests.append(("Need at least 4 genomes to build a tree", len(genomeIds) + len(optionalGenomeIds) >= 4))

if args.maxGenomesMissing:
    LOG.write("allowing %d genomes missing per PGfam (out of %d total)\n"%(args.maxGenomesMissing, len(genomeIds)))
environmentTests.append(("Num genomes - maxGenomesMissing >= 4", len(genomeIds) - args.maxGenomesMissing >= 4))

# call to getSingleCopyHomologs uses main genome, optional genomes are not involved in selecting single copy homologs
singleCopyHomologs = phylocode.selectSingleCopyHomologs(homologMatrix, genomeIds, requiredGenome=args.focusGenome, maxGenomesMissing=args.maxGenomesMissing, maxAllowedDups=args.maxAllowedDups)

LOG.write("got single copy homologs, num=%d\n"%len(singleCopyHomologs))
LOG.flush()
singleCopyGenesPerGenome = {}
for genome in genesPerGenome:
    singleCopyGenesPerGenome[genome] = 0
for homolog in singleCopyHomologs:
    for genome in homologMatrix[homolog]:
        singleCopyGenesPerGenome[genome] += 1

filesToMoveToDetailsFolder =  []
filesToDelete = []

environmentTests.append(("single copy homologs > 0", len(singleCopyHomologs) > 0))
## perform preflight test

maxGenesWithExcess = int(args.maxGenes * (1+args.excessGenesProp) + args.excessGenesFixed)
if len(singleCopyHomologs) > maxGenesWithExcess:
    singleCopyHomologs=singleCopyHomologs[0:maxGenesWithExcess]
    LOG.write("\tevaluating alignments of %d single-family genes, excess of %d\n"%(len(singleCopyHomologs), maxGenesWithExcess - args.maxGenes))

proteinAlignments = {}
proteinAlignmentStats = {}
alignmentScore = {}
perGenomeLogFreq = {}
genomeAlignmentCount = {}
for genomeId in genomeIds:
    perGenomeLogFreq[genomeId] = 0
    genomeAlignmentCount[genomeId] = 0  # for denominator of perSeqMeanLogFreq
alignedTaxa=set()
protein_alignment_time = time()
for homologId in singleCopyHomologs:
        LOG.write("aligning {}\n".format(homologId))
        geneIdSet = set()
        for genome in homologMatrix[homologId]:
            for proteinId in homologMatrix[homologId][genome]:
                if not "undefined" in proteinId:
                    geneIdSet.add(proteinId)
            #geneIdSet.update(set(homologMatrix[homologId][genome]))

        proteinFasta = patric_api.getSequenceOfFeatures(geneIdSet, 'protein')
        # replace bad characters for good while it is still text (e.g., 'J' to 'X')
        lines = proteinFasta.split("\n")
        for i in range(len(lines)):
            if not lines[i].startswith(">"):
                lines[i] = lines[i].replace("J", "X") #shouldn't occur, but avoid crashing downstream analysis programs
        proteinFasta = "\n".join(lines)

        seqRecords = SeqIO.parse(StringIO(proteinFasta), "fasta") #, alphabet=) #IUPAC.extended_protein)
        proteinSeqDict = SeqIO.to_dict(seqRecords)

        for genomeId in homologMatrix[homologId]:
            for geneId in homologMatrix[homologId][genomeId]:
                if genomeId == genomeObject_genomeId:
                    proteinSeqDict[geneId] = genomeObjectProteins[geneId]
                if geneId in proteinSeqDict:
                    proteinSeqDict[geneId].annotations["genome_id"] = genomeId
        all_ok = True
        for feature_id in proteinSeqDict:
            if "undefined" in feature_id:
                problem_seq = proteinSeqDict[feature_id] # seqrecord with "undefined" in id
                comment = "homology group {} has undefined identfier {}, skipping it\nrecord={}\n".format(homologId, feature_id, problem_seq)
                LOG.write(comment)
                #sys.stderr.write(comment)
                all_ok = False
        if not all_ok:
            continue # skip this homology group
        if args.debugMode:
            LOG.write("protein set for %s has %d seqs\n"%(homologId, len(proteinSeqDict)))
        if args.aligner == 'muscle':
            proteinAlignment = phylocode.alignSeqRecordsMuscle(proteinSeqDict.values())
        else:
            proteinAlignment = phylocode.alignSeqRecordsMafft(proteinSeqDict.values())
        if not proteinAlignment:
            continue
        proteinAlignment = phylocode.resolveDuplicatesPerPatricGenome(proteinAlignment)
        proteinAlignment.sort()
        proteinAlignments[homologId] = proteinAlignment
        alignmentStats = phylocode.calcAlignmentStats(proteinAlignment)
        # dividing by sqrt(alignment length) yields result intermediate between sum_squared_freq and mean_squared_freq (where squared_freq is the within-column sum of squared state (letter, amino acid) frequency
        if alignmentStats['num_pos']:
            alignmentScore[homologId] = alignmentStats['sum_squared_freq'] / sqrt(alignmentStats['num_pos'])
        else:
            alignmentScore[homologId] = 0
        for figId in alignmentStats['perseq_meanlogfreq']:
            genomeId = genomeIdFromFigId(figId)
            perGenomeLogFreq[genomeId] += alignmentStats['perseq_meanlogfreq'][figId]
            genomeAlignmentCount[genomeId] += 1
        proteinAlignmentStats[homologId] = alignmentStats

protein_alignment_time = time() - protein_alignment_time
LOG.write("Protein alignments completed. Number = %d, time = %.1f\n"%(len(proteinAlignments), protein_alignment_time))

# normalize perGenomeLogFreq by how many alignments there were per genome
for genomeId in perGenomeLogFreq:
    if genomeAlignmentCount[genomeId]:
        perGenomeLogFreq[genomeId] /= genomeAlignmentCount[genomeId]

F = open("alignment_perseq_meanlogfreq.txt", "w")
F.write("PGFam\tgenome\tfigtail\tgenome_mean\tgene_lf\tscore\n")
# consider perGenomeLogFreq as a filter to disqualify alignments with one or more really bad sequence
for homologId in proteinAlignmentStats:
    alignmentStats = proteinAlignmentStats[homologId]
    worstScore = 0
    for figId in alignmentStats['perseq_meanlogfreq']:
        genomeId = genomeIdFromFigId(figId)
        seqLogFreq =  alignmentStats['perseq_meanlogfreq'][figId]
        badSeqScore = seqLogFreq / (perGenomeLogFreq[genomeId] - 1e-5)
        if badSeqScore > worstScore:
            worstScore = badSeqScore
        F.write("{}\t{}\t{}\t{:.5f}\t{:.5f}\t{:.5f}\n".format(homologId, genomeId, figId[-7:], perGenomeLogFreq[genomeId], seqLogFreq, badSeqScore))
        if badSeqScore >  8 : # scores are positive, positive excursion is bad
            if 'bad_seqs' not in proteinAlignmentStats[homologId]:
                alignmentStats['bad_seqs'] = {}
            alignmentStats['bad_seqs'][genomeId] = badSeqScore
            LOG.write("homolog {} has bad sequence {} scoring {} relative to mean {}\n".format(homologId, figId, badSeqScore, perGenomeLogFreq[genomeId]))
    alignmentStats['worstSeqScore'] = worstScore
F.close()

# select top alignments by score
singleCopyHomologs = []
homologsWithOneBadSeq = []
for homologId in sorted(alignmentScore, key=alignmentScore.get, reverse=True):
    if 'bad_seqs' in proteinAlignmentStats[homologId]:
        proteinAlignmentStats[homologId]['disqualification'] = 'bad seqs'
        if len(proteinAlignmentStats[homologId]['bad_seqs']) == 1:
            homologsWithOneBadSeq.append(homologId)
    else:
        if len(singleCopyHomologs) < args.maxGenes:
            singleCopyHomologs.append(homologId)
        else:
            proteinAlignmentStats[homologId]['disqualification'] = 'low alignment score'
LOG.write("\tselected top %d single-family genes based on alignment score\n"%len(singleCopyHomologs))
LOG.write("\talignments with one bad seq: {}\n".format(len(homologsWithOneBadSeq)))

filteredSingleCopyGenesPerGenome = {}
for genome in genesPerGenome:
    filteredSingleCopyGenesPerGenome[genome] = 0
for homolog in singleCopyHomologs:
    for genome in homologMatrix[homolog]:
        filteredSingleCopyGenesPerGenome[genome] += 1

genesPerGenomeFile = os.path.join(args.outputDirectory, args.outputBase+".genesPerGenome.txt")
filesToMoveToDetailsFolder.append(args.outputBase+".genesPerGenome.txt")
with open(genesPerGenomeFile, 'w') as F:
    F.write("Genome\tAllGenes\tSingleCopy\tFiltered_SingleCopy\n")
    for genome in sorted(genesPerGenome, key=genesPerGenome.get):
        F.write("%s\t%d\t%d\t%d\n"%(genome, genesPerGenome[genome], singleCopyGenesPerGenome[genome], filteredSingleCopyGenesPerGenome[genome]))

passed = True
LOG.write("Initial checks of environment:\n")
for test in environmentTests:
    passed &= test[1]
    LOG.write(test[0] + "\t" + str(test[1]) + "\n")
LOG.write("All tests passed\t"+str(passed)+"\n")

LOG.flush()
codonAlignments = {}
sum_protein_positions = 0
for homologId in singleCopyHomologs:
    proteinAlignment = proteinAlignments[homologId]
    if args.debugMode:
        LOG.write("alignment for %s has %d seqs\n"%(homologId, len(proteinAlignment)))
        LOG.write("   ids: ")
        for seqrecord in proteinAlignment:
            LOG.write(", "+seqrecord.id)
        LOG.write("\n")
    codonAlignment = phylocode.gapCdsToProteins(proteinAlignment, genomeObjectGeneDna)
    phylocode.relabelSequencesByGenomeId(codonAlignment)
    phylocode.relabelSequencesByGenomeId(proteinAlignment)
    if codonAlignment.get_alignment_length() % 3:
        raise Exception("codon alignment length not multiple of 3 for %s\n"%homologId)
    if args.endGapTrimThreshold:
        codonAlignment = phylocode.trimEndGaps(codonAlignment, args.endGapTrimThreshold)
    codonAlignments[homologId] = codonAlignment
    if args.debugMode:
        LOG.write("dna alignment for %s has %d seqs\n"%(homologId, len(codonAlignment)))
    for seqRecord in proteinAlignment:
        alignedTaxa.add(seqRecord.id)
    # trim protein alignment after codon alignment if trimming enabled
    if args.endGapTrimThreshold:
        proteinAlignment = phylocode.trimEndGaps(proteinAlignment, args.endGapTrimThreshold)
        proteinAlignments[homologId] = proteinAlignment 

numTaxa=len(alignedTaxa)

LOG.write("Codon alignments completed. Number = %d\n"%(len(codonAlignments)))
LOG.flush()

phyloFileBase = args.outputBase 

# change to output directory to simplify file naming
os.chdir(args.outputDirectory)

proteinPositions=0
codonPositions = 0
protein_substitution_model = args.proteinModel
rate_model = args.rateModel
if numTaxa > 50:
    rate_model = 'CAT'  # faster with this many taxa, per RAxML manual

raxmlCommand = []
raxml_command_lines = []
raxml_analysis_goal = []
raxml_process_time = []
partitionFile = ''
if len(proteinAlignments):
# write the genes included in each homology group (and those paralogs excluded)
    homologsAndGenesIncludedInAlignmentFile = phyloFileBase+".homologsAndGenesIncludedInAlignment.txt"
    filesToMoveToDetailsFolder.append(homologsAndGenesIncludedInAlignmentFile)
    with open(homologsAndGenesIncludedInAlignmentFile, 'w') as F:
        for homologId in singleCopyHomologs:
            if homologId in proteinAlignments:
                proteinPositions += proteinAlignments[homologId].get_alignment_length()
                genesNotIncluded = set()
                for genome in homologMatrix[homologId]:
                    genesNotIncluded.update(set(homologMatrix[homologId][genome]))
                genesIncluded = set()
            F.write(homologId+"\tProteins\t")
            for seqRecord in proteinAlignments[homologId]:
                originalId = seqRecord.annotations['original_id']
                F.write("\t"+originalId)
                if originalId not in genesNotIncluded:
                    LOG.write("Problem: originalId %s not in genesForHomologs for %s\n"%(originalId, homologId))
                else:
                    genesNotIncluded.remove(originalId)
                genesIncluded.add(originalId)
            if len(genesNotIncluded):
                F.write("\tdeletedParalogs: \t"+"\t".join(genesNotIncluded))
            F.write("\n")
            if homologId in codonAlignments:
                F.write(homologId+"\tCodons\t")
                codonPositions += codonAlignments[homologId].get_alignment_length()
                for seqRecord in codonAlignments[homologId]:
                    originalId = seqRecord.annotations['original_id']
                    F.write("\t"+originalId)
                    genesIncluded.remove(originalId)
                if len(genesIncluded):
                    F.write("\tlackingCodonAlignment: \t"+"\t".join(genesIncluded))
                F.write("\n")
            else:
                F.write(homologId+"\nNo codon alignment\n")
        F.write("\n")
    F.close()

    alignmentStatsFile = phyloFileBase+".homologAlignmentStats.txt"
    filesToMoveToDetailsFolder.append(alignmentStatsFile)
    with open(alignmentStatsFile, "w") as F:
        #first = True
        keystats = ['gaps', 'sum_squared_freq', 'mean_squared_freq', 'num_pos', 'num_seqs', 'worstSeqScore']
        F.write("PGFam\t"+"\t".join(keystats)+"used_in_tree\n")
        for homolog in sorted(alignmentScore, key=alignmentScore.get, reverse=True): #proteinAlignments:
            stats = proteinAlignmentStats[homolog]
            #if first:
            #    F.write("PGFam\t"+"\t".join(sorted(stats.keys()))+"\tUsedInAnalysis\n")
            #    first = False
            F.write(homolog)
            for key in keystats:
                F.write("\t{}".format(stats[key]))
            F.write("\t"+str(homolog in singleCopyHomologs))
            F.write("\n")
            if 'bad_seqs' in stats:
                F.write('bad seqs: {}\n'.format(stats['bad_seqs']))
            if 'disqualification' in stats:
                F.write('disqualification: {}\n'.format(stats['disqualification']))

# change proteinAlignments to only include selected ones
    selectedAlignments = {}
    for homolog in singleCopyHomologs:
        selectedAlignments[homolog] = proteinAlignments[homolog]
    proteinAlignments = selectedAlignments

    if args.writePgfamAlignments:
        for homolog in proteinAlignments:
            SeqIO.write(proteinAlignments[homolog], homolog+".afa", "fasta")
    if args.writePgfamAlignmentsDNA:
        for homolog in proteinAlignments:
            if homolog in codonAlignments:
                SeqIO.write(codonAlignments[homolog], homolog+".afn", "fasta")

# it is possible all codon alignments failed, remove from intent to analyze
    if len(codonAlignments) == 0:
        args.analyzeCodons = False

# finally, output concatenated protein and/or DNA alignment and partitions and raxml command to appropriate files
    startTree = None
    if args.analyzeProteins and args.proteinModel == "AUTO" and not args.deferRaxml:
        # analyze just protein data with model AUTO, parse output to find best protein model
        # use output tree as starting tree to save time
        LOG.write("running raxml on proteins with model = 'AUTO' to find best model for proteins.\n")
        positions_to_sample = proteinPositions
        if args.protein_sample_prop > 1.0 or args.protein_sample_prop < 0:
            args.protein_sample_prop = 1.0
        sampled_prots = phylocode.sample_and_concatenate_alignments(proteinAlignments, args.protein_sample_prop)
        proteinFileBase = phyloFileBase+"_proteins"
        proteinAlignmentFile = proteinFileBase+".afa"
        F = open(proteinAlignmentFile, 'w')
        for seq_id in sorted(sampled_prots):
            F.write(">"+seq_id+"\n"+sampled_prots[seq_id]+"\n")
        F.close()
        raxmlCommand = [args.raxmlExecutable, "-s", proteinAlignmentFile, "-n", proteinFileBase, "-m", "PROTCATAUTO", "-p", "12345", "-T", str(args.threads), '-e', '10']
        LOG.write("command = "+" ".join(raxmlCommand)+"\n")
        raxml_command_lines.append(" ".join(raxmlCommand))
        raxml_analysis_goal.append("Analyze proteins with model 'AUTO' to find best substitution model.")
        proc_start_time = time()
        return_code = subprocess.call(raxmlCommand, shell=False) #, stdout=LOG, stderr=LOG)
        raxml_process_time.append(time() - proc_start_time)
        LOG.write("return code from 'AUTO' run was "+str(return_code)+"\n")
        if return_code != 0:
            LOG.write("raxml run to find best protein model failed.\n")

        filesToMoveToDetailsFolder.append(proteinAlignmentFile)
        if os.path.exists(proteinAlignmentFile+".reduced"):
            filesToDelete.append(proteinAlignmentFile+".reduced")
        protein_substitution_model = "LG" # default if we don't find answer
        if os.path.exists("RAxML_info."+proteinFileBase):
            filesToMoveToDetailsFolder.append("RAxML_info."+proteinFileBase)
            F = open("RAxML_info."+proteinFileBase)
            bestModel = None
            for line in F:
                m = re.match(r"\s+Partition: 0 best-scoring AA model:\s+(\S+)", line)
                if m:
                    bestModel = m.group(1)
            if bestModel:
                protein_substitution_model = bestModel
                LOG.write("Analysis of proteins found best model: "+bestModel+"\n")
            startTree = "RAxML_bestTree."+proteinFileBase
        else:
            LOG.write("Could not find output from running raxml on proteins alone to find best protein model. Defaulting to '{}'\n".format(protein_substitution_model))
    raxmlCommand=[]
    if args.analyzeProteins and args.analyzeCodons:
        alignmentFile = phyloFileBase+".afa"
        filesToMoveToDetailsFolder.append(alignmentFile)
        F = open(alignmentFile, "w")
        seqDict = phylocode.concatenate_codons_proteins(codonAlignments, proteinAlignments)
        for seqId in seqDict:
            F.write(">"+seqId+"\n"+seqDict[seqId]+"\n")
        F.close()
        
        partitionFile = phyloFileBase+".partitions"
        filesToMoveToDetailsFolder.append(partitionFile)
        with open(partitionFile, 'w') as PART:
            for i in range(1,4):
                PART.write("DNA, codon%d = %d-%d\\3\n"%(i, i, codonPositions))
            PART.write("%s, proteins = %d-%d\n"%(protein_substitution_model, codonPositions+1, codonPositions+proteinPositions))
        raxmlCommand = [args.raxmlExecutable, "-s", alignmentFile, "-n", phyloFileBase, "-m",  "GTR"+rate_model, "-q",  partitionFile,  "-p", "12345", "-T", str(args.threads)]

    elif args.analyzeCodons:
        alignmentFile = phyloFileBase+".phy"
        filesToMoveToDetailsFolder.append(alignmentFile)
        phylocode.writeConcatenatedAlignmentsPhylip(codonAlignments, alignmentFile)
        with open(phyloFileBase+".partitions", 'w') as PartitionFile:
            for i in range(1,4):
                PartitionFile.write("DNA, codon%d = %d-%d\\3\n"%(i, i, codonPositions))
        raxmlCommand = [args.raxmlExecutable, "-s", alignmentFile, "-n", phyloFileBase, "-m",  "GTR%s"%rate_model, "-q",  phyloFileBase+".partitions",  "-p", "12345", "-T", str(args.threads)]
        filesToMoveToDetailsFolder.append(phyloFileBase+".partitions")

    elif args.analyzeProteins:
        alignmentFile = phyloFileBase+".phy"
        filesToMoveToDetailsFolder.append(alignmentFile)
        phylocode.writeConcatenatedAlignmentsPhylip(proteinAlignments, alignmentFile)
        raxmlCommand = [args.raxmlExecutable, "-s", alignmentFile, "-n", phyloFileBase, "-m",  "PROT%s%s"%(rate_model, protein_substitution_model), "-p", "12345", "-T", str(args.threads)]
    #raxmlCommand.extend(["-e", "0.1"]) #limit on precision, faster than default 0.1

    if args.bootstrapReps > 0:
        raxmlCommand.extend(["-f", "a", "-x", "12345", "-N", str(args.bootstrapReps)]) 
    else:
        raxmlCommand.extend(["-f", "D"]) 

    raxmlCommandFile = phyloFileBase+".raxmlCommand.sh"
    with open(raxmlCommandFile, 'w') as F:
        F.write(" ".join(raxmlCommand)+"\n")
    filesToMoveToDetailsFolder.append(raxmlCommandFile)
    raxml_command_lines.append(" ".join(raxmlCommand))
    raxml_analysis_goal.append("Find best tree.")

genomeIdToName = patric_api.getNamesForGenomeIdsByN(allGenomeIds)
svgTreeImage = None
program_return_value = 0 
if len(proteinAlignments) and not args.deferRaxml:
    goal = "Find best tree."
    if args.bootstrapReps > 0:
        goal = "Find best tree and perform 'fast bootstrapping'."
    raxml_analysis_goal.append(goal)
    LOG.write("prepare to run raxml.\n")
    #remove RAxML files that clash in name, their existence blocks raxml from running
    for fl in glob.glob("RAxML_*"+phyloFileBase):
        os.remove(fl)
    raxml_start_time = time()
    result_code = subprocess.call(raxmlCommand, shell=False) #, stdout=LOG, stderr=LOG)
    raxml_process_time.append(time() - raxml_start_time)
    
    LOG.write("raxml completed: elapsed seconds = %f\n"%(time() - raxml_start_time))
    LOG.flush()
    if genomeObject:
        genomeIdToName[genomeObject_genomeId] = genomeObject_name+" "+genomeObject_genomeId
    originalNewick = ""
    raxmlNewickFileName = "RAxML_bestTree."+phyloFileBase
    if args.bootstrapReps > 0:
        raxmlNewickFileName = "RAxML_bipartitions."+phyloFileBase
    elif os.path.exists("RAxML_rellBootstrap."+phyloFileBase):
        fileBase = "tree_with_rell_support"
        raxmlCommand = [args.raxmlExecutable, "-s", alignmentFile, "-n", fileBase, "-m",  "GTRCAT", "-f", "b", "-t", "RAxML_bestTree."+phyloFileBase, "-z", "RAxML_rellBootstrap."+phyloFileBase]
        raxml_command_lines.append(" ".join(raxmlCommand))
        raxml_analysis_goal.append("Map RELL support values onto best tree.")
        proc_start_time = time()
        result_code = subprocess.call(raxmlCommand, shell=False) #/, stdout=LOG, stderr=LOG)
        raxml_process_time.append(time() - proc_start_time)
        LOG.write("Mapped ultrafast bootstrap support values onto tree.\n")
        if os.path.exists("RAxML_bipartitions."+fileBase):
            raxmlNewickFileName = "RAxML_bipartitions."+fileBase
            filesToMoveToDetailsFolder.append("RAxML_info."+fileBase)
        if os.path.exists(alignmentFile+".reduced"):
            filesToDelete.append(alignmentFile+".reduced")
        
    treeWithGenomeIdsFile = None
    if not os.path.exists(raxmlNewickFileName):
        LOG.write("Expected tree file %s not found\n" % (raxmlNewickFileName) );
        program_return_value = 1 # return this to signal no tree was generated == failure
    else:
        program_return_value = 0 # means success
        # generate rooted or balanced tree
        LOG.write("rootTree({})\n".format(raxmlNewickFileName))
        output_suffix = phyloFileBase + "_rooted"
        raxmlCommand = [args.raxmlExecutable, '-t', raxmlNewickFileName, '-f', 'I', '-m', 'GTRCAT', '-n', output_suffix] 
        LOG.write("command:\n{}\n".format(" ".join(raxmlCommand)))
        raxml_command_lines.append(" ".join(raxmlCommand))
        raxml_analysis_goal.append("Root the tree by RAxML 'balance' metric (like mid-point rooting).")
        proc_start_time = time()
        result_code = subprocess.call(raxmlCommand, shell=False) #/, stdout=LOG, stderr=LOG)
        raxml_process_time.append(time() - proc_start_time)
        if os.path.exists("RAxML_rootedTree."+output_suffix):
            os.rename("RAxML_rootedTree."+output_suffix, phyloFileBase + "_tree_rooted.nwk")
            filesToMoveToDetailsFolder.append(phyloFileBase + "_tree_rooted.nwk")
            LOG.write("rename {} to {}\n".format("RAxML_rootedTree."+output_suffix, phyloFileBase + "_rooted_tree.nwk"))
        else:
            LOG.write("generating balanced tree failed\n")

        F = open(raxmlNewickFileName)
        originalNewick = F.read()
        F.close()
        treeWithGenomeIdsFile = phyloFileBase + "_tree.nwk"
        #filesToMoveToDetailsFolder.append(treeWithGenomeIdsFile)
        F = open(treeWithGenomeIdsFile, 'w')
        F.write(originalNewick)
        F.close()
        renamedNewick = phylocode.relabelNewickTree(originalNewick, genomeIdToName)
        renamedNewickFile = phyloFileBase+"_treeWithGenomeNames.nwk"
        filesToMoveToDetailsFolder.append(renamedNewickFile)
        F = open(renamedNewickFile, 'w')
        F.write(renamedNewick)
        F.close()
        LOG.write("codonTree newick relabeled with genome names written to "+renamedNewickFile+"\n")
        LOG.flush()


        if args.writePhyloxml and treeWithGenomeIdsFile:
            LOG.write("writePhyloxml\n")
            command = ["p3x-newick-to-phyloxml","-l","genome_id"]
            if args.phyloxmlFields and len(args.phyloxmlFields):
                command.extend(("-g", args.phyloxmlFields))
            if groupsPerGenome:
                with open("groupsPerGenome.tsv", "w") as F:
                    F.write("genomeId\tGenome_Group\n")
                    for genomeId in groupsPerGenome:
                        F.write("{}\t{}\n".format(genomeId, ",".join(groupsPerGenome[genomeId])))
                command.extend(["--annotationtsv", "groupsPerGenome.tsv"])
                filesToMoveToDetailsFolder.append("groupsPerGenome.tsv")
            command.append(treeWithGenomeIdsFile)
            if args.debugMode:
                sys.stderr.write("calling p3x-newick-to-phyloxml, command line:\n"+" ".join(command)+"\n")
            result_code = subprocess.call(command, shell=False) #/, stdout=LOG, stderr=LOG)
            phyloxml_file = treeWithGenomeIdsFile.replace(".nwk", ".phyloxml")


        figtree_found = phylocode.checkCommandline("figtree")
        if figtree_found or (args.pathToFigtreeJar and os.path.exists(args.pathToFigtreeJar)):
            nexusFilesWritten = phylocode.generateNexusFile(originalNewick, phyloFileBase, nexus_template = None, align_tips = "no", focus_genome = args.focusGenome, genomeIdToName=genomeIdToName)
            LOG.write("nexus file written to %s\n"%(", ".join(nexusFilesWritten)))
            filesToMoveToDetailsFolder.append(nexusFilesWritten[0])
            if len(nexusFilesWritten) > 1:
                filesToDelete.append(nexusFilesWritten[1]) # get rid of codontree_tipsAligned.nex
            for nexusFile in nexusFilesWritten:
                #for imageFormat in ("SVG"): #("PNG", "SVG", "PDF"):
                imageFormat = "SVG"
                imageFile = phylocode.generateFigtreeImage(nexusFile, numTaxa=len(allGenomeIds), figtreeJar=args.pathToFigtreeJar, imageFormat=imageFormat)
                LOG.write("created figtree figure: {}\n".format(imageFile))
                if os.path.exists(imageFile):
                    if "_tipsAligned" in nexusFile:
                        filesToMoveToDetailsFolder.append(imageFile)
                    if imageFormat == "SVG" and not svgTreeImage:
                        svgTreeImage = imageFile
                        LOG.write("svg image for report: {}\n".format(imageFile))
                        filesToMoveToDetailsFolder.append(svgTreeImage)
                else:
                    LOG.write("image file {} does not exist.\n".format(imageFile))

LOG.write(strftime("%a, %d %b %Y %H:%M:%S", localtime(time()))+"\n")
LOG.write("Total job duration %d seconds\n"%(time()-starttime))
        
analysisStatsFile = phyloFileBase+".analysisStats"
OUT = open(analysisStatsFile, 'w')
OUT.write("Statistics for CodonTree Analysis\n")
OUT.write("Num_genomes\t%s\n"%numTaxa)
OUT.write("Num_protein_alignments\t%s\n"%len(proteinAlignments))
OUT.write("Num_aligned_amino_acids\t%s\n"%proteinPositions)
OUT.write("Num_CDS_alignments\t%s\n"%len(proteinAlignments))
OUT.write("Num_aligned_nucleotides\t%s\n"%codonPositions)
OUT.write("PGFams\t%s\n"%",".join(sorted(proteinAlignments)))
OUT.write("command_line\t%s\n"%" ".join(raxmlCommand))
OUT.write("Total job duration %d seconds\n"%(time()-starttime))
OUT.close()
filesToMoveToDetailsFolder.append(analysisStatsFile)


htmlFile = os.path.abspath(args.outputBase+"_report.html")
HTML = open(htmlFile, 'w')
HTML.write("<h1>Phylogenetic Tree Report</h1>\n")
if svgTreeImage is not None and os.path.exists(svgTreeImage):
    HTML.write("<h2>Rendered Tree</h2>\n")
    HTML.write(open(svgTreeImage).read()+"\n\n")
alternateSvgFile = "%s_tipsAligned.svg"%phyloFileBase
if os.path.exists(alternateSvgFile):
    HTML.write("<p><a target='_parent' href='detail_files/%s'>Alternate View</a></p>"%alternateSvgFile)

HTML.write("<p><a target='_parent' href='detail_files/'>Files with more details</a></p>")


HTML.write("<h2>Tree Analysis Statistics</h2>\n<table border='1'>\n")
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Requested genomes", num_requested_genomes))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Genomes with data", len(genomeIds)))
if len(deletedGenomes):
    HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Genomes lacking data", len(deletedGenomes)))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Max allowed deletions", args.maxGenomesMissing))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Max allowed duplications", args.maxAllowedDups))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Single-copy genes requested", args.maxGenes))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Single-copy genes found", len(singleCopyHomologs)))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Num protein alignments", len(proteinAlignments)))
HTML.write("<tr><td><b>%s</b></td><td>%s</td></tr>\n"%("Alignment program", args.aligner))
HTML.write("<tr><td><b>%s</b></td><td>%.1f seconds</td></tr>\n"%("Protein alignment time", protein_alignment_time))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Num aligned amino acids", proteinPositions))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Num CDS alignments", len(codonAlignments)))
HTML.write("<tr><td><b>%s</b></td><td>%d</td></tr>\n"%("Num aligned nucleotides", codonPositions))
if args.analyzeProteins:
    title = "Protein substitution model"
    if protein_substitution_model != args.proteinModel:
        title = "Best protein model found by RAxML"
    HTML.write("<tr><td><b>%s</b></td><td>%s</td></tr>\n"%(title, protein_substitution_model))
m = re.search("/awe/work/../../([^_]+)_/", os.getcwd())
if m:
    patricJobId = m.group(1)
    HTML.write("<tr><td><b>PATRIC Job Idx=</b></td><td>"+patricJobId+"</td></tr>\n")
# tree analysis may or may not have completed, report only if appropriate
branch_support_method = "RELL (Resampling Estimated Log Likelihoods)"
if args.bootstrapReps > 0:
    branch_support_method = "RAxML Fast Bootstrapping"
HTML.write("<tr><td><b>%s</b></td><td>%s</td></tr>\n"%("Branch support method", branch_support_method))

raxmlInfoFile = "RAxML_info."+phyloFileBase
raxmlWarnings = []
if os.path.exists(raxmlInfoFile):
    filesToMoveToDetailsFolder.append(raxmlInfoFile)
    try:
        raxmlInfo = open(raxmlInfoFile).read()
        m = re.search(r"Final.*core of best tree ([-\d\.]+)", raxmlInfo)
        if not m:
            m = re.search(r"Final ML Optimization Likelihood: ([-\d\.]+)", raxmlInfo)
        if m:
            raxmlLikelihood = m.group(1) 
            raxmlLikelihood = float(raxmlLikelihood)
            HTML.write("<tr><td><b>%s</b></td><td>%.4f</td></tr>\n"%("RAxML likelihood", raxmlLikelihood))
        m = re.search(r"RAxML version ([\d\.]+)", raxmlInfo)
        if m:  
            raxmlVersion = m.group(1)
            HTML.write("<tr><td><b>%s</b></td><td>%s</td></tr>\n"%("RAxML version", raxmlVersion))
        for m in re.finditer("IMPORTANT WARNING: (Sequences.*?identical)", raxmlInfo):
            raxmlWarnings.append(m.group(1))
        m = re.search(r"Overall execution time.*: ([\d\.]*) secs", raxmlInfo) 
    except Exception as e:
        LOG.write("Exception parsing %s: %s\n"%(raxmlInfoFile, str(e))) 
raxmlDuration = 0 
for seconds in raxml_process_time:
    raxmlDuration += seconds
HTML.write("<tr><td><b>%s</b></td><td>%.1f seconds</td></tr>\n"%("RAxML time", raxmlDuration))
HTML.write("<tr><td><b>%s</b></td><td>%.1f seconds</td></tr>\n"%("Total time", time()-starttime))
HTML.write("</table>\n\n")
if raxml_command_lines and len(raxml_command_lines):
    HTML.write("<h2>RAxML Command Line</h2>")
    for i, line in enumerate(raxml_command_lines):
        if i < len(raxml_analysis_goal):
            HTML.write("<p>Goal: "+raxml_analysis_goal[i]+"<br>\n")
        HTML.write(line+"<br>\n")
        if i < len(raxml_process_time):
            HTML.write("Process time: {:.3f} seconds<br>\n".format(raxml_process_time[i]))
else:
    HTML.write("<h2>RAxML Not Run</h2>\n")

if os.path.exists(partitionFile):
    HTML.write("<h2>RAxML Codon and Amino Acid Partitions</h2>\n")
    HTML.write("<pre>\n"+open(partitionFile).read()+"</pre>\n")

if len(raxmlWarnings):
    HTML.write("<h2>RAxML Warnings</h2>\n")
    for line in raxmlWarnings:
        HTML.write(line+"<br>\n")

HTML.write("<h2>Genome Statistics</h2>\n<table border='1'>\n")
HTML.write("<tr><th>GenomeId</th><th>Total Genes</th><th>Single Copy</th><th>Used</th><th>Name</th></tr>\n")
for genome in sorted(genesPerGenome, key=genesPerGenome.get):
    if genome not in genomeIds:
        LOG.write("%s not in genomeIds\n"%genome)
    elif genome not in singleCopyGenesPerGenome:
        LOG.write("%s not in singleCopyGenesPerGenome\n"%genome)
    elif genome not in genomeIdToName:
        LOG.write("%s not in genomeIdToName\n"%genome)
    else:
        HTML.write("<tr><td>%s</td><td>%d</td><td>%d</td><td>%d</td>"%(genome, genesPerGenome[genome], singleCopyGenesPerGenome[genome], filteredSingleCopyGenesPerGenome[genome]))
        if genome not in genomeIdToName:
            genomeIdToName[genome] = genome
        HTML.write("<td><a htarget='_blank' ref='https://patricbrc.org/view/Genome/%s'>%s</a></td></tr>\n"%(genome, genomeIdToName[genome]))
for genome in sorted(deletedGenomes):
    HTML.write("<tr><td>%s</td><td>0</td><td>0</td><td>0</td><td>Deleted</td></tr>\n"%(genome)) 

HTML.write("</table>\n\n")

if len(alignmentScore):
    HTML.write("<h2>Gene Family Statistics</h2>\n<table border='1'>\n")
    HTML.write("<p>Gene families are ranked by alignment score combining mean per-position variability, alignment length, and gappiness.</p>\n")
    homologProduct = patric_api.getProductsForPgfamsByN(alignmentScore.keys())
    statsToShow = ['num_pos', 'num_seqs', 'mean_squared_freq', 'prop_gaps']
    altHeadings = ['Align.<br>Length', 'Num<br>Seqs', 'Mean<br>Sqr Freq', 'Prop<br>Gaps']
    HTML.write("<tr><th>PGFam</th><th>Align.<br>Score</th><th>"+"</th><th>".join(altHeadings)+"</th><th>Used In<br>Analysis</th><th>Product</th></tr>\n")
    for homolog in sorted(alignmentScore, key=alignmentScore.get, reverse=True): #proteinAlignments:
        stats = proteinAlignmentStats[homolog]
        HTML.write("<tr><td>%s</td><td>%.2f</td>"%(homolog, alignmentScore[homolog]))
        for key in statsToShow:
            val = stats[key]
            if isinstance(val,int):
                HTML.write("<td>%d</td>"%stats[key])
            else:
                HTML.write("<td>%.3f</td>"%stats[key])
        HTML.write("<td>%s</td>\n"%str(homolog in singleCopyHomologs))
        HTML.write("<td>%s</td></tr>\n"%homologProduct[homolog])
    HTML.write("</table>\n")

if args.debugMode or len(singleCopyHomologs) < args.maxGenes:
    HTML.write("<h2>Strategies to Increase Single-Copy Gene Number</h2>\n")
    if len(singleCopyHomologs) < args.maxGenes:
        HTML.write("Number of single-copy genes (%d) was less than requested (%d).<br>\n"%(len(singleCopyHomologs), args.maxGenes))
    HTML.write("Examining the number of genes per genome in the 'Genome Statistics' table above may indicate incomplete or plasmid entries with few genes which may be removed.<br>\n")
    HTML.write("Criteria for calling single copy genes can be made more lenient by increasing Max Allowed Deletions and/or Max Allowed Duplications.<br>\n")
    #compute single copy counts for subsets of taxa, helps user identify optimal subsets for further analysis
    genomeSubsetSingleCopy = phylocode.countSingleCopyForGenomeSubsets(homologMatrix, genomeIds, maxAllowedDups = args.maxAllowedDups)
    if args.debugMode:
        LOG.write("len(genomeSubsetSingleCopy) = %d\n"%len(genomeSubsetSingleCopy))
    if len(genomeSubsetSingleCopy) > 0:
        HTML.write("<p>Omitting one of the following sets of genomes will provide the approximate boost in single-copy genes:</p>\n")
        HTML.write("<table border='1'>\n")
        HTML.write("<tr><th>scGene Boost</th><th>Genomes to Omit</th></tr>\n")
        maxToWrite = 5
        for subsetTuple in sorted(genomeSubsetSingleCopy, key=genomeSubsetSingleCopy.get, reverse=True):
            omissionSet = genomeIds - set(subsetTuple)
            HTML.write("<tr><td>%d</td><td>%s</td></tr>\n"%(genomeSubsetSingleCopy[subsetTuple], " ".join(sorted(omissionSet))))
            maxToWrite -= 1
            if not maxToWrite:
                break
        HTML.write("</table>\n")
HTML.write("</html>\n")
HTML.close()


LOG.write("output written to directory %s\n"%args.outputDirectory)
sys.stdout.write("\n")
detailsDirectory = "detail_files"
LOG.write("details files moved to subdirectory %s\n"%detailsDirectory)
if not os.path.isdir(detailsDirectory):
    if os.path.exists(detailsDirectory):
        os.remove(detailsDirectory)
    os.mkdir(detailsDirectory)
numMoved=0
for fn in filesToMoveToDetailsFolder:
    LOG.write("\t"+fn+"\t"+os.path.join(detailsDirectory, fn)+"\n")
    if os.path.exists(fn):
        os.rename(fn, os.path.join(detailsDirectory, fn))
        numMoved += 1
    else:
        LOG.write("tried to move {} to detailsDirectory, but not found.\n".format(fn))
LOG.write("files moved: %d\n"%numMoved)
filesToDelete.extend(glob.glob("RAxML*"))
if args.debugMode:
    numDeleted = 0
    for fn in filesToDelete:
        os.remove(fn)
        numDeleted += 1
        LOG.write("\t"+fn+"\n")
    LOG.write("files deleted: %d\n"%numDeleted)
logfileName = os.path.basename(logfileName)
# finally, move the log file into the detailsDirectory
if os.path.exists(logfileName):
    LOG.write("finally, will move this file, %s, to %s\n"%(logfileName, detailsDirectory))
    LOG.close()
    os.rename(logfileName, os.path.join(detailsDirectory, logfileName))

# return value indicates whether tree was constructed or not
sys.exit(program_return_value)
