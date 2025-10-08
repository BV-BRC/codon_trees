import sys
import re
import subprocess
import os
import os.path
import patric_api
from math import log
from random import sample
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonExperimentalWarning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import codonalign
try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3

Debug = False #shared across functions defined here
LOG = sys.stderr

genetic_code_table11 = {
        "TTT": "F",      "TCT": "S",      "TAT": "Y",      "TGT": "C",  
        "TTC": "F",      "TCC": "S",      "TAC": "Y",      "TGC": "C",  
        "TTA": "L",      "TCA": "S",      "TAA": "*",      "TGA": "*",  
        "TTG": "L",      "TCG": "S",      "TAG": "*",      "TGG": "W",  

        "CTT": "L",      "CCT": "P",      "CAT": "H",      "CGT": "R",  
        "CTC": "L",      "CCC": "P",      "CAC": "H",      "CGC": "R",  
        "CTA": "L",      "CCA": "P",      "CAA": "Q",      "CGA": "R",  
        "CTG": "L",      "CCG": "P",      "CAG": "Q",      "CGG": "R",  

        "ATT": "I",      "ACT": "T",      "AAT": "N",      "AGT": "S",  
        "ATC": "I",      "ACC": "T",      "AAC": "N",      "AGC": "S",  
        "ATA": "I",      "ACA": "T",      "AAA": "K",      "AGA": "R",  
        "ATG": "M",      "ACG": "T",      "AAG": "K",      "AGG": "R",  

        "GTT": "V",      "GCT": "A",      "GAT": "D",      "GGT": "G",  
        "GTC": "V",      "GCC": "A",      "GAC": "D",      "GGC": "G",  
        "GTA": "V",      "GCA": "A",      "GAA": "E",      "GGA": "G",  
        "GTG": "V",      "GCG": "A",      "GAG": "E",      "GGG": "G"
        }

def which(program):
    """ Can use to test for existence of needed files on path
    based on code from https://stackoverflow.com/users/20840/jay """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def getPgfamDistribution(genomeGenePgfamList):
    """ Given list of genome-gene-pgfam tuples, 
        tabulate counts per genome per pgfam, 
        num duplications per pgfam, 
        num single-copy.
        This information can be used to identify PGFams suitable for phylogenetic inference for these genomes.
        Alternatively, can use to identify genomes with sparse single-copy coverage to remove to create denser phylogenetic matrix. 
    """
    ggpMat = {} # genome-gene-pgfam matrix (really just a dictionary)
    for row in genomeGenePgfamList:
        genome, gene, pgfam = row
        if pgfam not in ggpMat:
            ggpMat[pgfam] = {}
        if genome not in ggpMat[pgfam]:
            ggpMat[pgfam][genome] = []
        ggpMat[pgfam][genome].append(gene)
    return ggpMat

def selectSingleCopyHomologs(pgfamMatrix, genomeIdList, requiredGenome=None, maxGenomesMissing=0, maxAllowedDups=0):
    # given a genome-gene-pgfam matrix
    # find the set of pgfam_ids which satisfy the maxGenomesMissing & maxAllowedDups criteria for specified genomes
    pgfamScore = {}
    if not len(genomeIdList) > 3:
        raise Exception("getSingleCopyPgfams: number of genome IDs too low: %d"%len(genomeIdList))
    minGenomesPresent = len(genomeIdList) - maxGenomesMissing
    if minGenomesPresent <= 0:
        raise Exception("getSingleCopyPgfams: maxGenomesMissing too large: %d"%maxGenomesMissing)
    for pgfam in pgfamMatrix:
        if requiredGenome and requiredGenome not in pgfamMatrix[pgfam]:
            continue # skip pgfam if it does not include required genome
        pgfamGenomeCount = 0
        numDups = 0
        for genome in genomeIdList:
            if genome in pgfamMatrix[pgfam]:
                pgfamGenomeCount += 1
                numDups += len(pgfamMatrix[pgfam][genome]) - 1
        if pgfamGenomeCount >= minGenomesPresent and numDups <= maxAllowedDups:
            pgfamScore[pgfam] = pgfamGenomeCount - numDups # combined score of both factors, for sorting
    suitablePgfamList = sorted(pgfamScore, key=pgfamScore.get, reverse=True)
    return suitablePgfamList

def countSingleCopyForGenomeSubsets(pgfamMatrix, genomeIds, maxAllowedDups=0):
    """ Return dict of how many single copy genes would be found if 
    certain taxa were omitted from the analysis. 
    Takes into consideration maxAllowedDups, but not maxGenomesMissing.
    Return dict of genome tuples to count of single copy genes IF that set is omitted.
    Record subsets down to half the original set size.
    """
    scForSubset = {}
    minSetSize = int(len(genomeIds) * 0.75)
    for pgfam in pgfamMatrix:
        numDups = 0
        present = []
        missing = 0
        for genome in sorted(pgfamMatrix[pgfam]):
            x = len(pgfamMatrix[pgfam][genome])
            if x == 0:
                missing += 1
            else:
                present.append(genome)
                if x > 1:
                    numDups += (x - 1)
        if False and Debug:
            LOG.write("countSingleCopyForGenomeSubsets, pgfam=%s, missing=%d, present=%d, nd=%d\n"%(pgfam, missing, len(present), numDups))

        if numDups <= maxAllowedDups and (len(present) >= minSetSize) and (len(present) < len(genomeIds)):
            presentTuple = tuple(present)
            if presentTuple not in scForSubset:
                scForSubset[presentTuple] = 0
            scForSubset[presentTuple] += 1
    return scForSubset

def getGenesForPgfams(genomeGenePgfam, genomeIdList, singleCopyPgfams):
    genesForPgfams={}
    for pgfam in singleCopyPgfams:
        genesForPgfams[pgfam] = []
    for row in genomeGenePgfam:
        genome, gene, pgfam = row
        if genome in genomeIdList and pgfam in genesForPgfams:
            genesForPgfams[pgfam].append(gene)
    return genesForPgfams

def relabelNewickTree(newick, labelDict):
# newick file has tip labels eg genomIDs
# return value has replaced these according to dictionary passed as labelDict
    retval = newick
    for oldLabel in labelDict:
        retval = re.sub("\\b"+oldLabel+"\\b", '"'+labelDict[oldLabel]+'"', retval)
    return retval

def writeTranslatedNexusTree(outFile, newick, labelDict, figtreeParameters=None, highlightGenome=None):
    # intended to use with FigTree to generate figures, but can be used without a figtreeParameters object
#write taxa block
    taxonIds = re.findall(r"[(,]([^:(),]+)[,:)]", newick)
    outFile.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%d;\n\ttaxlabels\n"%len(taxonIds))
    for tax in taxonIds:
        if tax not in labelDict:
            labelDict[tax] = tax
        outFile.write("\t\"%s %s\""%(labelDict[tax], tax))
        if tax == highlightGenome:
            outFile.write("[&!color=#ff0000]")
        outFile.write("\n")
    outFile.write(";\nend;\n\n")
# prefix support values with "[&label="
    newick = re.sub(r"\)(\d+):", r")[&support=\1]:", newick)
# write trees block
    outFile.write("begin trees;\n")
    if labelDict:
        outFile.write("\ttranslate\n")
        for tax in taxonIds:
            translate = tax
            if tax in labelDict:
                translate = labelDict[tax]
            outFile.write("\t\t%s \"%s %s\",\n"%(tax, translate, tax))
        outFile.write("\t;\n")
    outFile.write("\ttree one = [&U] %s\n"%newick)
    outFile.write("\nend;\n\n")
# write figtree block
    if figtreeParameters:
        figtreeParameters["nodeLabels.isShown"]="true"
        figtreeParameters["nodeLabels.displayAttribute"]="\"support\""
        outFile.write("begin figtree;\n")
        for param in sorted(figtreeParameters):
            outFile.write("\tset %s=%s;\n"%(param, str(figtreeParameters[param])))
        outFile.write("\nend;\n\n")
    return

def readFigtreeParameters(filename):
    infile = open(filename)
    retval = {}
    inFigtreeBlock = False
    for line in infile:
        if re.search("begin figtree", line, re.IGNORECASE):
            inFigtreeBlock = True
        if re.search(r"^end;", line):
            inFigtreeBlock = False
        if inFigtreeBlock:
            m = re.search(r"set\s+(\S.*)=(\S.*);", line, re.IGNORECASE)
            if m:
                retval[m.group(1)] = m.group(2)
    return retval

def generateNexusFile(newick, outfileBase, nexus_template = None, align_tips = "both", focus_genome = None, genomeIdToName=None):
    figtreeParams={}
    if not nexus_template:
        LOG.write("Look for figtree.nex template in sys.path directories\n")
        for dirname in sys.path: # should be in ../lib
            if os.path.isfile(os.path.join(dirname, "figtree.nex")):
                nexus_template = os.path.join(dirname, "figtree.nex")
                LOG.write("Found figtree template file: %s\n"%nexus_template)
    # read a model figtree nexus file
    if nexus_template and os.path.exists(nexus_template):
        LOG.write("Read figtree template file: %s\n"%nexus_template)
        figtreeParams = readFigtreeParameters(nexus_template)
    genomeIds = re.findall("[(,]([^(,):]+)[,:)]", newick)
    if not genomeIdToName:
        genomeIdToName = patric_api.getNamesForGenomeIds(genomeIds)
    figtreeParams["trees.rooting"]="true"
    figtreeParams["trees.rootingType"]="Midpoint"
    figtreeParams["trees.order"]="true"
    figtreeParams["trees.orderType"]="increasing"
    figtreeParams["tipLabels.fontName"]="sanserif"
    figtreeParams["tipLabels.fontSize"]=14
    figtreeParams["tipLabels.fontStyle"]=0
    figtreeParams["tipLabels.isShown"]="true"
    filesWritten=[]
    if align_tips in ("no", "both"):
        figtreeParams['rectilinearLayout.alignTipLabels'] = 'false'
        nexusOut = open(outfileBase+".nex", "w")
        writeTranslatedNexusTree(nexusOut, newick, genomeIdToName, figtreeParameters=figtreeParams, highlightGenome=focus_genome)
        nexusOut.close()
        filesWritten.append(outfileBase+".nex")
    if align_tips in ("yes", "both"):
        figtreeParams['rectilinearLayout.alignTipLabels'] = 'true'
        nexusOut = open(outfileBase+"_tipsAligned.nex", "w")
        writeTranslatedNexusTree(nexusOut, newick, genomeIdToName, figtreeParameters=figtreeParams, highlightGenome=focus_genome)
        nexusOut.close()
        filesWritten.append(outfileBase+"_tipsAligned.nex")
    return filesWritten

def generateFigtreeImage(nexusFile, numTaxa=0, figtreeJar=None, imageFormat="PDF"):
    if Debug:
        LOG.write("generateTreeFigure(%s, %d, figtreeJarFile=%s, imageFormat=%s)\n"%(nexusFile, numTaxa, figtreeJar, imageFormat))
    if imageFormat not in ('PDF', 'SVG', 'PNG', 'JPEG'):
        raise Exception("imageFormat %s not in ('PDF', 'SVG', 'PNG', 'JPEG')"%imageFormat)
    imageFileName = nexusFile
    imageFileName = nexusFile.replace(".nex", ".")
    imageFileName += imageFormat.lower()
    if figtreeJar and os.path.exists(figtreeJar):
        figtreeCommand = ['java',  '-jar', figtreeJar, '-graphic', imageFormat]
    else:
        # assume a figtree executable is on the path
        figtreeCommand = ['figtree', '-graphic', imageFormat]
    if numTaxa > 40:
        height = 600 + 15 * (numTaxa - 40) # this is an empirical correction factor to avoid taxon name overlap
        figtreeCommand.extend(['-height', str(int(height))])
    figtreeCommand.extend([nexusFile, imageFileName])
    if Debug:
        LOG.write("running this command:\n%s\n"%" ".join(figtreeCommand))
    subprocess.call(figtreeCommand, stdout=LOG)
    return imageFileName

def checkCommandline(program):
    found = False
    LOG.write("checking for {} on commandline:\n".format(program))
    try:
        subprocess.check_call(['which', program], stdout=LOG)
        found = True
    except Exception:
        pass
    LOG.write("found = {}\n".format(found))
    return found

def alignSeqRecordsMuscle(seqRecords):
    try:
        #python3
        muscleProcess = subprocess.Popen(['muscle', '-quiet'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    except TypeError:
        # python2
        muscleProcess = subprocess.Popen(['muscle', '-quiet'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    SeqIO.write(seqRecords, muscleProcess.stdin, 'fasta')
    muscleProcess.stdin.close()
    #alphabet=None
    #for s in seqRecords:
    #    alphabet=s.seq.alphabet
    #    break
    try:
        alignment = AlignIO.read(muscleProcess.stdout, "fasta") #, alphabet=) #alphabet)
        alignment.sort()
    except ValueError as ve:
        alignment = None
        comment = "Error in muscle alignment: {}\n".format(ve)
        sys.stderr.write(comment)
        LOG.write(comment)
    return(alignment)

def alignSeqRecordsMafft(seqRecords):
    mafftProcess = subprocess.Popen(['mafft', '--quiet', '--auto', '-'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    SeqIO.write(seqRecords, mafftProcess.stdin, 'fasta')
    mafftProcess.stdin.close()
    #alphabet=None
    #for s in seqRecords:
    #    alphabet=s.seq.alphabet
    #    break # only need first one
    alignment = None
    try:
        alignment = AlignIO.read(mafftProcess.stdout, "fasta") #, alphabet=alphabet)
        alignment.sort()
    except ValueError as ve:
        alignment = None
        comment = "Error in mafft alignment: {}\n".format(ve)
        sys.stderr.write(comment)
        LOG.write(comment)
    return(alignment)

def calcAlignmentStats(alignment):
    # analyze a BioPython MultipleSeqAlignment object to describe conservation levels
    stats = {'num_pos':0, 'num_seqs':0, 'sum_squared_freq':0, 'mean_squared_freq':0, 'gaps':0, 'prop_gaps':0}
    num_pos = alignment.get_alignment_length()
    num_seqs = len(alignment)
    if num_pos == 0 or num_seqs == 0: 
        return stats
    numGaps = 0
    numNonGaps = 0
    sumSquaredFreq = 0
    perSeqLogP = [0.0] * num_seqs
    avgEntropy = 0
    for pos in range(0, num_pos):
        stateCount = {}
        states = alignment[: , pos]
        for residue in states:
            if residue == '-':
                numGaps += 1
            else:
                numNonGaps += 1
            if residue not in stateCount:
                stateCount[residue] = 0
            stateCount[residue] += 1
        stateFreq = {}
        for residue in stateCount:
            freq = stateCount[residue]/float(num_seqs)
            stateFreq[residue] = freq
            if not residue == '-': # don't let gaps contribute to squared frequency 
                squaredFreq = freq*freq
                sumSquaredFreq += squaredFreq
        if len(stateCount) > 1:    
            #print("calc per seq mean log(freq)")
            for i, residue in enumerate(states):
                if residue in stateFreq:
                    perSeqLogP[i] += log(stateFreq[residue], 2) #log base 2
    stats['perseq_meanlogfreq'] = {}
    for i, record in enumerate(alignment):
        stats['perseq_meanlogfreq'][record.id] = perSeqLogP[i]/num_pos
            
    stats['num_pos'] = num_pos
    stats['num_seqs'] = num_seqs
    stats['sum_squared_freq'] = sumSquaredFreq
    stats['mean_squared_freq'] = sumSquaredFreq/num_pos
    stats['gaps'] = numGaps
    stats['prop_gaps'] = numGaps/float(numGaps+numNonGaps)
    return stats

def suggestAlignmentDeletions(alignment):
    """ 
    analyze a BioPython MultipleSeqAlignment object to calculate how many gaps could be avoided by omitting sequences
    return dict of tuple of seqIds to improvement score of deleting that set of seqs
    """
    delsets = {}
    for pos in range(0, alignment.get_alignment_length()):
        states = alignment[: , pos]
        gaps = 0.0
        nongap_seqs = []
        for i, residue in enumerate(states):
            if residue == '-':
                gaps += 1
            else:
                nongap_seqs.append(alignment[i].id)
        if gaps and gaps/len(alignment) >= 0.75: # a gap-heavy position
            key = tuple(nongap_seqs)
            if key not in delsets:
                delsets[key] = 0
            delsets[key] += 1
    return delsets

def calcSumAlignmentDistance(alignment, querySeq):
    sumDist = 0
    for record in alignment:
        for posPair in zip(querySeq, record.seq):
            if not posPair[0] == posPair[1]:
                sumDist += 1
        #sumDist += Levenshtein.distance(str(querySeq), str(record.seq))
    return(sumDist)

def trimEndGaps(alignment, trimThreshold=0.5):
    # trim leading and trailing gaps from alignment
    # trim columns where proportion of gaps is > trimThreshold, stop at first column below threshold
    # higher threshold means more trimming
    # return tuple with number of columns trimmed front and back
    if trimThreshold < 0 or trimThreshold >= 1.0:
        raise Exception("trimEndGaps: trimThreshold (%.2f) out of range 0.0-1.0"%trimThreshold)
    trimmableNumberOfGaps = int(trimThreshold * len(alignment)) # number of sequences with gaps allowed
    if trimmableNumberOfGaps == 0:
        if Debug:
            LOG.write("trimEndGaps: trimThreshold (%.2f) so lenient no trimming possible with %d sequences\n"%(trimThreshold, len(alignment)))
        #alignment.annotation["endgaps_trimmed"] = (0,0)
        return(alignment)
    leadGaps={}
    endGaps={}
    for record in alignment:
        leadGapLen = 0
        m = re.match("-+", str(record.seq))
        if m:
            leadGapLen=len(m.group(0))
        if leadGapLen not in leadGaps:
            leadGaps[leadGapLen] = 0
        leadGaps[leadGapLen] += 1
        endGapLen = 0
        m = re.search("-+$", str(record.seq))
        if m:
            endGapLen=len(m.group(0))
        if endGapLen not in endGaps:
            endGaps[endGapLen] = 0
        endGaps[endGapLen] += 1
    leadTrimLen = 0
    gappedSeqs = len(alignment) # start considering all sequence as gaps (at pos -1)
    for l in sorted(leadGaps):
        gappedSeqs -= leadGaps[l]
        if gappedSeqs <= trimmableNumberOfGaps:
            leadTrimLen = l
            break
    gappedSeqs = len(alignment) # start considering all sequence as gaps (at pos len+1)
    endTrimLen=0
    for l in sorted(endGaps):
        gappedSeqs -= endGaps[l]
        if gappedSeqs <= trimmableNumberOfGaps:
            endTrimLen = l
            break
    if leadTrimLen > 0 or endTrimLen > 0:
        endPos = alignment.get_alignment_length() - endTrimLen
        trimmedAlignment = alignment[:,leadTrimLen:endPos]
        for i in range(0, len(alignment)):
            trimmedAlignment[i].annotations = alignment[i].annotations
        alignment = trimmedAlignment
    #alignment.annotation["endgaps_trimmed"] = (leadTrimLen, endTrimLen)
    return(alignment)

def resolveDuplicatesPerPatricGenome(alignment):
# accepts Bio.MultipleSeqAlignment, returns list of seqIds to keep
# calculate average similarity/distance of each seqToResolve to entire alignment
# identify best seq per genome (highest avg similarity to rest of alignment) and save that one, remove others from seqIds list
# return list of seqIds to keep
    if not alignment:
        return None
    seqIds=list()
    genomesToResolve=set()
    seqsPerGenome = {}
    initialNumRecords=len(alignment)
    for record in alignment:
        seqIds.append(record.id)
        # assume PATRIC gene identifier like this: fig|1399771.3.peg.1094
        genomeId = ".".join(record.id.split(".")[:2]).split("|")[1]
        #LOG.write("%s\t%s\n"%(record.id, genomeId))
        if genomeId not in seqsPerGenome:
            seqsPerGenome[genomeId] = []
        else:
            genomesToResolve.add(genomeId)
        seqsPerGenome[genomeId].append(record.id)

    for genomeId in genomesToResolve:
        setToResolve = seqsPerGenome[genomeId]
        bestAlignDist = alignment.get_alignment_length()
        seqToKeep = None
        for seqId in setToResolve:
            dist = calcSumAlignmentDistance(alignment, record.seq)
            if seqToKeep == None:
                bestAlignDist=dist
                seqToKeep=seqId
            else:
                if dist < bestAlignDist:
                    bestAlignDist = dist
                    seqToKeep = seqId
        for seqId in setToResolve:
            if not seqId == seqToKeep:
                seqIds.remove(seqId)
    if len(seqIds) < initialNumRecords:
        if Debug:
            LOG.write("after resolveDups num seqIds is %d, versus prev %d\n"%(len(seqIds), initialNumRecords))
        reducedSeqs = []
        for record in alignment:
            if record.id in seqIds:
                record.seq = Seq(str(record.seq).replace('-', '')) #, record.seq.alphabet) # remove gaps
                reducedSeqs.append(record)
        alignment = alignSeqRecordsMuscle(reducedSeqs)
    return alignment

def gapCdsToProteins(proteinAlignment, extraDnaSeqs=None):
    """ to replace proteinToCodonAlignment() """
    protSeqDict = {}
    for seqRecord in proteinAlignment:
        protSeqDict[seqRecord.id] = seqRecord
    dnaFasta = patric_api.getSequenceOfFeatures(protSeqDict.keys(), 'dna')
    #if Debug:
    #     LOG.write("dnaFasta sample: %s\n"%dnaFasta[:100])

    dnaSeqDict = SeqIO.to_dict(SeqIO.parse(StringIO(dnaFasta), "fasta")) #, alphabet=IUPAC.IUPACAmbiguousDNA()))
    for seqId in protSeqDict:
        if extraDnaSeqs and seqId in extraDnaSeqs:
            dnaSeqDict[seqId] = extraDnaSeqs[seqId]
            if Debug:
                LOG.write("appending extra DNA seq %s\n"%seqId)
    if set(dnaSeqDict.keys()) != set(protSeqDict.keys()):
        raise Exception("Protein and DNA sets differ:\nProteins: %s\nDNA: %s\n"%(", ".join(sorted(protSeqDict)), ", ".join(sorted(dnaSeqDict))))
    dnaAlignFasta = StringIO()
    prot_align_len = proteinAlignment.get_alignment_length()
    for seqId in dnaSeqDict:
        dnaSeq = dnaSeqDict[seqId].seq
        if len(dnaSeq) < 3 * prot_align_len:
            # this is to handle cases where protein exists but DNA does not
            dnaSeq += '---' * (prot_align_len - len(dnaSeq))
        protSeq = protSeqDict[seqId].seq
        dnaAlignFasta.write(">"+seqId+"\n")
        dnaSeqPos = 0
        codon_to_aa_mismatch = 0
        for protPos in range(0, len(protSeq)):
            if protSeq[protPos] == '-':
                codon = '---'
            else:
                codon = str(dnaSeq[dnaSeqPos:dnaSeqPos+3]).upper()
                if len(codon) != 3:
                    LOG.write(" failed to get full codon for {} at aa-pos {} of {}\n".format(seqId, protPos, prot_align_len))
                    codon = '---'
                dnaSeqPos += 3
                # check whether codon codes for amino acid
                if Debug and (codon in genetic_code_table11):
                    if (genetic_code_table11[codon] != protSeq[protPos]):
                        codon_to_aa_mismatch += 1
                        LOG.write(" aa to codon mismatch at protpos {} in {}, codon is {} => {}\n".format(protPos, seqId, codon, genetic_code_table11[codon]))
            dnaAlignFasta.write(codon)
        protPos += 1 # should now be equal to prot_align_len
        if Debug:
            LOG.write(seqId+"\tprotPos={}\torig_prot={}\tdnaPos={}\torig_DNA={}\tgcmisses={}\n".format(protPos, len(protSeq), dnaSeqPos, len(dnaSeq), codon_to_aa_mismatch))
        if protPos < prot_align_len:
            dnaAlignFasta.write(''.join("---"*(prot_align_len - protPos)))
            LOG.write("padding short seq {0}, of {1} pos out to {2}, orig_DNA_len={3}, orig_prot_len={4}\n".format(seqId, protPos, prot_align_len, len(dnaSeq), len(protSeq)))
        
        dnaAlignFasta.write("\n")
    dnaAlignFasta_text = dnaAlignFasta.getvalue()
    retval = AlignIO.read(StringIO(dnaAlignFasta_text), 'fasta')
    return retval

def proteinToCodonAlignment(proteinAlignment, extraDnaSeqs = None):
    protSeqDict = {}
    for seqRecord in proteinAlignment:
        protSeqDict[seqRecord.id] = seqRecord
    dnaFasta = patric_api.getSequenceOfFeatures(protSeqDict.keys(), 'dna')
    #if Debug:
    #     LOG.write("dnaFasta sample: %s\n"%dnaFasta[:100])

    dnaSeqDict = SeqIO.to_dict(SeqIO.parse(StringIO(dnaFasta), "fasta")) #, alphabet=IUPAC.IUPACAmbiguousDNA()))
    for seqId in protSeqDict:
        if extraDnaSeqs and seqId in extraDnaSeqs:
            dnaSeqDict[seqId] = extraDnaSeqs[seqId]
            if Debug:
                LOG.write("appending extra DNA seq %s\n"%seqId)
    if set(dnaSeqDict.keys()) != set(protSeqDict.keys()):
        raise Exception("Protein and DNA sets differ:\nProteins: %s\nDNA: %s\n"%(", ".join(sorted(protSeqDict)), ", ".join(sorted(dnaSeqDict))))
    for seqId in dnaSeqDict:
        if not len(dnaSeqDict[seqId].seq):
            #del(dnaSeqDict[seqId])
            LOG.write("warning: seqId %s length of dna was zero\n"%seqId)
    dnaSeqRecords=[]
    for proteinSeq in proteinAlignment:
        dnaSeqRecords.append(dnaSeqDict[proteinSeq.id])

    if Debug:
        LOG.write("dna seqs has %d seqs\n"%(len(dnaSeqRecords)))
        #LOG.write("DNA seq ids: %s\n"%(", ".join(sorted(dnaSeqDict))))
        #LOG.write("pro seq ids: %s\n"%(", ".join(sorted(protSeqDict))))
        #LOG.write("first two aligned DNA seqs:\n")
        #SeqIO.write(dnaSeqRecords[:2], LOG, "fasta")
        #LOG.flush()
   
    """
    # now check length of protein vs dna sequences, extend dna if needed to make match in numbers of codons
    for i, protRec in enumerate(proteinAlignment):
        protSeq = str(protRec.seq)
        protSeq.replace('-','')
        protLen = len(protSeq)
        if len(dnaSeqs[i].seq) < protLen*3:
            shortfall = (protLen*3) - len(dnaSeqs[i].seq)
            if Debug:
                LOG.write("DNA seq for %s is too short for protein, shortfall = %d\n"%(protRec.id, shortfall))
            # extend on both ends to be safe
            dnaSeqs[i].seq = "N"*shortfall + dnaSeqs[i].seq + "N"*shortfall
    """
    returnValue = None
    #with warnings.catch_warnings():
        #warnings.simplefilter('ignore', BiopythonWarning)
        #try:
    #ambiguous_nucleotide_values = {'K': 'GT', 'M': 'AC', 'N': 'ACGT', 'S': 'CG', 'R': 'AG', 'W': 'AT', 'Y': 'CT'}
    #ambiguous_protein_values = {'X': 'ACDEFGHIKLMNOPQRSTVWY', 'J': 'IL', 'B': 'DN', 'Z': 'EQ'}
    #ambiguous_codon_table = CodonTable.AmbiguousCodonTable(CodonTable.ambiguous_dna_by_name["Standard"], IUPAC.IUPACAmbiguousDNA(), ambiguous_nucleotide_values, IUPAC.protein, ambiguous_protein_values)
    #returnValue = codonalign.build(pro_align=proteinAlignment, nucl_seqs=dnaSeqRecords, codon_table=ambiguous_codon_table, max_score=1000)
    returnValue = codonalign.build(pro_align=proteinAlignment, nucl_seqs=dnaSeqRecords, max_score=1000)
    for dnaSeq in returnValue:
        proteinRecord = protSeqDict[dnaSeq.id]
        if proteinRecord.annotations:
            dnaSeq.annotations = proteinRecord.annotations.copy()

        #except Exception as e:
        #    LOG.write("problem in codonalign, skipping\n%s\n"%str(e))
        #    raise(e)
    return returnValue
    
def relabelSequencesByGenomeId(seqRecordSet):
    #rename sequences by genome instead of sequence Id
    for seqRecord in seqRecordSet:
        originalId = seqRecord.id
        if originalId.startswith(">"):
            orignalId = originalId[1:]
        genomeId = ".".join(seqRecord.id.split(".")[:2]).split("|")[1]
        seqRecord.id = genomeId
        #stash original ID in the annotations dictionary of the seqRecord
        if not seqRecord.annotations:
            seqRecord.annotations={}
        seqRecord.annotations['original_id'] = originalId

def trimAlignments(alignmentDict, endGapTrimThreshold=0.5):
    if endGapTrimThreshold == 0:
        return
    for alignmentId in alignmentDict:
        alignmentDict[alignmentId] = trimEndGaps(alignmentDict[alignmentId], endGapTrimThreshold)
    return

def writeOneAlignmentPhylip(alignment, destination, idList, outputIds=True, codonPos=0):
    nameFieldWidth = 10
    if outputIds:
        for seqid in idList:
            nameFieldWidth = max(nameFieldWidth, len(seqid))
    theseIds = set(rec.id for rec in alignment)
    if len(theseIds) < len(idList):
        #nullSeq = UnknownSeq(alignment.get_alignment_length(), alphabet=alignment._alphabet)
        nullSeq = '-'*alignment.get_alignment_length()
        if codonPos:
            nullSeq = nullSeq[:len(nullSeq)/3]
    seqDict = {}
    for record in alignment:
        seqDict[record.id] = record.seq
    if not outputIds:
        destination.write("\n")
    for seqid in idList:
        if outputIds:
            destination.write("{:{width}}  ".format(seqid, width = nameFieldWidth))
        if seqid in theseIds:
            seqStr = str(seqDict[seqid])
            if codonPos:
                codonPosStr = ''
                for i in range(len(seqStr)):
                    if ((i+1) % 3) == codonPos:
                        codonPosStr += seqStr[i]
                seqStr = codonPosStr
            destination.write(seqStr+"\n")
        else:
            destination.write(str(nullSeq)+"\n")

def writeConcatenatedAlignmentsPhylip(alignments, destination):
# alignments is dictionary of Bio.multipleSeqAlignments
    if type(destination) == str:
        destination = open(destination, "w")
    taxonIdSet = set(seq.id for aln in alignments for seq in alignments[aln])
    totalLength = 0
    for alignmentId in alignments:
        totalLength += alignments[alignmentId].get_alignment_length()
    destination.write("%d\t%d\n"%(len(taxonIdSet), totalLength))
    taxonIdList=sorted(taxonIdSet)
    alignmentIdList = sorted(alignments)
    writeOneAlignmentPhylip(alignments[alignmentIdList[0]], destination, taxonIdList, outputIds=True)
    for alignmentId in alignmentIdList[1:]:
        writeOneAlignmentPhylip(alignments[alignmentId], destination, taxonIdList, outputIds=False)

def sample_and_concatenate_alignments(alignments, sample_prop=1.0):
# alignments is dictionary of Bio.multipleSeqAlignments
    id_set = set(seq.id for aln in alignments for seq in alignments[aln])
    totalLength = 0
    taxonIdList=sorted(id_set)
    alignmentIdList = sorted(alignments)
    seqDict = {}
    for seqid in id_set:
        seqDict[seqid] = ""
    for alignmentId in alignmentIdList:
        alignment = alignments[alignmentId]
        columns_to_write = None
        sample_size = alignment.get_alignment_length()
        if sample_prop < 1.0:
            sample_size = int(alignment.get_alignment_length()*sample_prop)
            columns_to_write = sample(range(alignment.get_alignment_length()), sample_size)
        temp_id_set = set(id_set)
        for record in alignment:
            if columns_to_write:
                for i in columns_to_write:
                    seqDict[record.id] += str(record.seq[i]) 
            else:
                seqDict[record.id] += str(record.seq)
            temp_id_set.remove(record.id)
        for missing_id in temp_id_set:
            seqDict[missing_id] += '-' * sample_size
    return seqDict

def outputCodonsProteinsPhylip(codonAlignments, proteinAlignments, destination):
    if Debug:
        LOG.write("outputCodonsProteinsPhylip, ncodonAln=%d, nprotAln=%d, destType=%s\n"%(len(codonAlignments), len(proteinAlignments), type(destination)))
    if len(codonAlignments) == 0:
        LOG.write("outputCodonsProteinsPhylip() called with zero codonAlignments()\n")
        return
    if type(destination) == str: # or type(destination) == unicode:
        if Debug:
            LOG.write("outputCodonsProteinsPhylip opening file %s\n"%destination)
        destination = open(destination, "w")
    codonPositions = 0
    taxonSet=set()
    for alignmentId in codonAlignments:
        codonPositions += codonAlignments[alignmentId].get_alignment_length()
        for seqRecord in codonAlignments[alignmentId]:
            taxonSet.add(seqRecord.id)
    proteinPositions = 0
    for alignmentId in proteinAlignments:
        proteinPositions += proteinAlignments[alignmentId].get_alignment_length()
        for seqRecord in proteinAlignments[alignmentId]:
            taxonSet.add(seqRecord.id)
    taxonIdList = sorted(taxonSet)
    #destination = open(directory+fileBase+"+codonsAndProteins.phy", 'w')
    destination.write("%d\t%d\n"%(len(taxonIdList), codonPositions+proteinPositions))
    alignmentIdList = sorted(codonAlignments)
    writeOneAlignmentPhylip(codonAlignments[alignmentIdList[0]], destination, taxonIdList, outputIds=True)
    for alignmentId in alignmentIdList[1:]:
        writeOneAlignmentPhylip(codonAlignments[alignmentId], destination, taxonIdList, outputIds=False)
    for alignmentId in sorted(proteinAlignments):
        writeOneAlignmentPhylip(proteinAlignments[alignmentId], destination, taxonIdList, outputIds=False)
    destination.write("\n")
    destination.close()

    return()

def concatenate_codons_proteins(codonAlignments, proteinAlignments, codonPos=0):
    if Debug:
        LOG.write("output_codons_proteins_fasta, ncodonAln=%d, nprotAln=%d\n"%(len(codonAlignments), len(proteinAlignments)))
    if len(codonAlignments) + len(proteinAlignments) == 0:
        LOG.write("output_codons_proteins_fasta() called with no alignments()\n")
        return
    taxonCodons={}
    taxonSet = set()
    for alignmentId in codonAlignments:
        for seqRecord in codonAlignments[alignmentId]:
            if not seqRecord.id in taxonCodons:
                taxonCodons[seqRecord.id] = {}
            fullseq = str(seqRecord.seq)
            if codonPos:
                codonPosSeq = ''
                for i in range(0, len(fullseq), 3):
                    codonPosSeq +=fullseq[i]
                taxonCodons[seqRecord.id][alignmentId] = codonPosSeq
            else:
                taxonCodons[seqRecord.id][alignmentId] = fullseq
            taxonSet.add(seqRecord.id)
    taxonProteins={}
    for alignmentId in proteinAlignments:
        for seqRecord in proteinAlignments[alignmentId]:
            if not seqRecord.id in taxonProteins:
                taxonProteins[seqRecord.id] = {}
            taxonProteins[seqRecord.id][alignmentId] = str(seqRecord.seq)
            taxonSet.add(seqRecord.id)
    taxonIdList = sorted(taxonSet)
    #destination = open(directory+fileBase+"+codonsAndProteins.phy", 'w')
    #destination.write("%d\t%d\n"%(len(taxonIdList), codonPositions+proteinPositions))
    alignmentIdList = sorted(codonAlignments)
    proteinIdList = sorted(proteinAlignments)
    seqDict = {}
    for taxon in taxonSet:
        seqDict[taxon] = ""
        for alignmentId in alignmentIdList:
            if taxon in taxonCodons and alignmentId in taxonCodons[taxon]:
                seqDict[taxon] += taxonCodons[taxon][alignmentId]
            else:
                al_len = codonAlignments[alignmentId].get_alignment_length()
                if codonPos:
                    al_len /= 3
                seqDict[taxon] += "-"*al_len+"\n"
        for alignmentId in sorted(proteinAlignments):
            al_len = proteinAlignments[alignmentId].get_alignment_length()
            if taxon in taxonProteins and alignmentId in taxonProteins[taxon]:
                seqDict[taxon] += taxonProteins[taxon][alignmentId]
            else:
                seqDict[taxon] += "-"*al_len+"\n"
    return seqDict

