import os
import sys
import requests
try: # first try python3 method
    import urllib.parse as urllib_parse
except: # fall back to python2
    import urllib as urllib_parse
import copy
import time
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from requests.packages.urllib3.exceptions import InsecureRequestWarning
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
import subprocess

Debug = False #shared across functions defined here
LOG = sys.stderr
#Base_url="https://www.patricbrc.org/api/"
Base_url="https://p3.theseed.org/services/data_api/";

Session = requests.Session()
Session.headers.update({ 'accept': "text/tsv" })
Session.headers.update({ "Content-Type": "application/rqlquery+x-www-form-urlencoded" })

PatricUser = None

def setDebug(state):
    Debug = state
    LOG.write("setting Debug to {}\n".format(Debug))
    Session.verify= not Debug # turn off authentiction if in debug mode

def authenticateByFile(tokenFile=None):
    if not tokenFile:
        tokenFile = os.path.join(os.environ.get('HOME'), ".patric_token")
    if os.path.exists(tokenFile):
        LOG.write("reading auth key from file %s\n"%tokenFile)
        with open(tokenFile) as F:
            tokenString = F.read().rstrip()
            authenticateByString(tokenString)

def authenticateByEnv():
    if "KB_AUTH_TOKEN" in os.environ:
        LOG.write("reading auth key from environment\n")
        authenticateByString(os.environ['KB_AUTH_TOKEN'])

def authenticateByString(tokenString):
    Session.headers.update({ 'Authorization' : tokenString })
    if "Authorization" in Session.headers:
        global PatricUser
        PatricUser = Session.headers["Authorization"].split(r"|")[3].split("=")[1]
        LOG.write("Patric user = %s\n"%PatricUser)

def getGenomeIdsNamesByName(name, limit='10'):
    query = "eq(genome_name,%s)"%name
    query += "&select(genome_id,genome_name)"
    query += "&limit(%s)"%limit
    ret = Session.get(Base_url+"genome/", params=query)
    if Debug:
        LOG.write(ret.url+"\n")
    return(ret.text.replace('"', ''))

def getGenomeGroupIds(genomeGroupName):
    LOG.write("getGenomeGroupIds({}), PatricUser={}, debug={}\n".format(genomeGroupName, PatricUser, Debug))
    genomeGroupSpecifier = urllib_parse.quote_plus(genomeGroupName)   #.replace("/", "%2f")
    genomeGroupSpecifier = genomeGroupSpecifier.replace('+', '%20')
    query = "in(genome_id,GenomeGroup("+genomeGroupSpecifier+"))"
    query += "&select(genome_id)"
    query += "&limit(10000)"
    if 1 or Debug:
        LOG.write("debug = {}\n".format(Debug))
        LOG.write("query =  {}\n".format(query))
    ret = Session.get(Base_url+"genome/", params=query)
    if Debug:
        LOG.write(ret.url+"\n")
    return(ret.text.replace('"', '').split("\n"))[1:-1]

def getGenomeGroupIdsCLI(genomeGroupName):
    # use Patric/BVBRC CLI to get genome group members
    LOG.write("getGenomeGroupIdsCLI({}), PatricUser={}, debug={}\n".format(genomeGroupName, PatricUser, Debug))
    command = ['p3-get-genome-group', genomeGroupName]
    proc = subprocess.run(command, shell=False, capture_output = True, text = True)
    retval = []
    for line in proc.stdout.split("\n")[1:]:
        retval.append(line)
    return retval

def getNamesForGenomeIds(genomeIds):
#    return getDataForGenomes(genomeIdSet, ["genome_id", "genome_name"])
    retval = {}
    for genome in genomeIds:
        retval[genome] = ""
    query="in(genome_id,("+",".join(genomeIds)+"))&select(genome_id,genome_name)"
    response = Session.get(Base_url+"genome/", params=query) #, 
    if Debug:
        LOG.write("    response URL: %s\n"%response.url)
        LOG.write("    len(response.text)= %d\n"%len(response.text))
    if not response.ok:
        LOG.write("Error code %d returned by %s in getNamesForGenomeIds\n"%(response.status_code, response.url))
    for line in response.text.split("\n"):
        line = line.replace('"','')
        row = line.split("\t", 1)
        if len(row) >= 2:
            genome, name = row
            retval[genome] = name
    return retval

def getNamesForGenomeIdsByN(genomeIds, n=5):
    """ For some reason, grabbing them in bulk misses some, so grab N at a time.
    """
    retval = {}
    i = 0
    genomeIds = list(genomeIds)
    while i < len(genomeIds):
        subset = genomeIds[i:i+n]
        retval.update(getNamesForGenomeIds(subset))
        i += n
    return retval


def getGenomeIdByFieldValue(queryField, queryValue):
    query = "eq(%s,%s)"%(queryField, queryValue)
    query += "&select(genome_id)"
    req = Session.get(Base_url+"genome/", params=query) 
    if Debug:
        LOG.write("getGenomeIdsByQuery: "+req.url+"\n")
        LOG.write(req.text+"\n")
    data = req.text.split("\n")
    genomeId = ""
    if len(data) > 1:
        genomeId = data[1]
        genomeId = genomeId.replace('\"', '')
    return genomeId

def getDataForGenomes(genomeIdSet, fieldNames):
    query = "in(genome_id,(%s))"%",".join(genomeIdSet)
    if fieldNames:
        query += "&select(%s)"%",".join(fieldNames)
    query += "&limit(%s)"%len(genomeIdSet)

    response = Session.get(Base_url+"genome/", params=query)
    if Debug:
        LOG.write("getDataForGenomes:\nurl="+response.url+"\nquery="+query+"\n")
    if not response.ok:
        errorMessage = "Error code %d returned by %s in getDataForGenomes\n"%(response.status_code, response.url)
        LOG.write(errorMessage)
        LOG.write("length of query was %d\n"%len(query))
        LOG.write("url="+response.url+"\nquery="+query+"\n")
        raise Exception(errorMessage)
    data = response.text.replace('"','') #get rid of quotes
    rows = data.split("\n")[:-1] # leave off empty last element
    retval = []
    for row in rows:
        fields = row.split("\t")
        #if len(fields) != len(fieldNames):
         #   continue
        retval.append(fields)
    return(retval)

def getSequenceOfFeatures(feature_ids, seq_type='dna'):
    # can specify dna or protein, replaces getProteinFastaForPatricIds and getDNAFastaForPatricIds
    max_per_query=100 # avoid 414 error if too many ids are passed on one query
    start = 0
    retval = ""
    feature_ids = list(feature_ids) # so we can index over them
    num_tries = 10
    output_ids = set() # to exclude duplicate IDs in output (can happen with 'undefined')
    while start < len(feature_ids):
        end = start + max_per_query
        end = min(end, len(feature_ids))
        #query="in(patric_id,("+",".join(map(urllib.quote, feature_ids[start:end]))+"))"
        query="in(patric_id,("+",".join(feature_ids[start:end])+"))"
        query += "&limit(%d)"%len(feature_ids)
        sys.stderr.write(f"getSequenceOfFeatures(), type={seq_type}\n")
        if Debug:
            sys.stderr.write("query =: {}\n".format(query))
        
        response=Session.get(Base_url+"genome_feature/", params=query, headers={'Accept': 'application/%s+fasta'%seq_type}, verify= not Debug)
        if response.ok:
            seqId = None
            for line in response.text.split("\n"):
                if line.startswith(">"):
                    seqId = line.split(" ", 1)[0] # space delimited
                    parts = line.split("|")
                    if len(parts) > 2:
                        seqId = "|".join(parts[:2])
                    if seqId in output_ids:
                        LOG.write("duplicate ID in getSequenceOfFeatures: {}\n".format(seqId))
                        seqId = None # invalidate a duplicated ID
                    else:
                        output_ids.add(seqId)
                        line = seqId
                if seqId:
                    retval += line+"\n"
            start += max_per_query
        else: # response not ok
            num_tries -= 1
            errorMessage= "Error code %d returned by %s in getSequenceOfFeatures\nnumber of tries left = %d\n"%(response.status_code, Base_url, num_tries)
            LOG.write(errorMessage)
            LOG.flush()
            if num_tries <= 0:
                raise Exception(errorMessage)
            time.sleep(100 - (num_tries*10)) # grow sleep time each time
    return retval

def getProteinsFastaForGenomeId(genomeId):
    query="in(genome_id,("+genomeId+"))"
    query += "&limit(25000)"
    response=Session.get(Base_url+"genome_feature/", params=query, headers={'Accept': 'application/protein+fasta'})
    if Debug:
        LOG.write("getProteinsFastaForGenomeId:\nurl="+response.url+"\nquery="+query+"\n")
    if not response.ok:
        LOG.write("Error code %d returned by %s in getProteinsFastaForGenomeId\n"%(response.status_code, Base_url))
        errorMessage= "Error code %d returned by %s in getProteinsFastaForGenomeId\nlength of query was %d\n"%(response.status_code, Base_url, len(query))
        LOG.write(errorMessage)
        LOG.flush()
        raise Exception(errorMessage)
    idsFixedFasta=""
    for line in response.text.split("\n"):
        if line.startswith(">"):
            parts = line.split("|")
            if len(parts) > 2:
                line = "|".join(parts[:2])+"\n"
        idsFixedFasta += line
    return idsFixedFasta

def getProductsForPgfams(pgfams):
    retval = {}
    for pgfam in pgfams:
        retval[pgfam] = ""
    query="in(family_id,("+",".join(pgfams)+"))&select(family_id,family_product)"
    response = Session.get(Base_url+"protein_family_ref/", params=query) #, 
    if Debug:
        LOG.write("    response URL: %s\n"%response.url)
        LOG.write("    len(response.text)= %d\n"%len(response.text))
    if not response.ok:
        LOG.write("Error code %d returned by %s in getProductsForPgfams\n"%(response.status_code, response.url))
    for line in response.text.split("\n"):
        line = line.replace('"','')
        row = line.split("\t", 1)
        if len(row) >= 2:
            pgfam, product = row
            retval[pgfam] = product
    return retval

def getProductsForPgfamsByN(pgfams, n=5):
    """ For some reason, grabbing them in bulk misses some, so grab N at a time.
    """
    retval = {}
    i = 0
    pgfams = list(pgfams)
    while i < len(pgfams):
        subset = pgfams[i:i+n]
        retval.update(getProductsForPgfams(subset))
        i += n
    return retval

def getPatricGenePosForGenome(genomeId):
    if Debug:
        LOG.write("getPatricGenesPosForGenome() called for %s\n"%genomeId)
    retval = []
    query="and(%s,%s,%s)"%("eq(genome_id,(%s))"%genomeId, "eq(feature_type,CDS)", "eq(pgfam_id,PGF*)")
    query += "&select(genome_id,patric_id,pgfam_id,accession,start,end,strand)"
    query += "&limit(25000)"
    response = Session.get(Base_url+"genome_feature/", params=query) 
    if Debug:
        LOG.write("    response URL: %s\n"%response.url)
        LOG.write("    len(response.text)= %d\n"%len(response.text))
    for line in response.text.split("\n"):
        line = line.replace('"','')
        row = line.split("\t")
        if len(row) != 7:
            continue
        retval.append(row)
    return retval

def get_homologs_for_genomes(genomeIdSet, scope='global'):
    if Debug:
        LOG.write("get_homologs_for_genomes() called for %d genomes\n"%len(genomeIdSet))
        LOG.write("    Session headers=\n"+str(Session.headers)+"\n")
    retval = []
    # one genome at a time, so using 'get' should be fine
    target_family_type = ('plfam_id', 'pgfam_id')[scope == 'global'] #test is index into alternatives
    for genomeId in genomeIdSet:
        query="and(%s,%s,%s)"%("eq(genome_id,(%s))"%genomeId, "eq(feature_type,CDS)", "eq("+target_family_type+",P*)") 
        # select genes from specified genome, be a CDS, and have a pgfam or plfam id starting with 'P'
        query += "&select(genome_id,patric_id,"+target_family_type+")"
        query += "&limit(25000)"
        response = Session.get(Base_url+"genome_feature/", params=query, verify= not Debug) #, 
        """
        req = requests.Request('POST', Base_url+"genome_feature/", data=query)
        prepared = Session.prepare_request(req) #req.prepare()
        response=Session.send(prepared, verify=False)
        """
        if Debug:
            LOG.write("    response URL: %s\n"%response.url)
            LOG.write("    len(response.text)= %d\n"%len(response.text))
        curLen = len(retval)
        for line in response.text.split("\n"):
            line = line.replace('"','')
            row = line.split("\t")
            if len(row) != 3:
                continue
            if not row[2].startswith("P"):
                continue
            retval.append(row)
        if Debug:
            LOG.write("    got %d %ss for that genome\n"%((len(retval)-curLen), target_family_type))
    return(retval)

def get_homolog_gene_matrix(genomeIdSet, ggpMat = None, scope='global'):
    """ Given list of genome ids: 
        tabulate genes per genome per homolog (pgfam or plfam) 
        (formats data from get_homologs_for_genomes as table)
    """
    genome_gene_homolog_list = get_homologs_for_genomes(genomeIdSet, scope=scope)
    if not ggpMat: # if a real value was passed, extend it
        ggpMat = {} # genome-gene-homolog matrix (really just a dictionary)
    for row in genome_gene_homolog_list:
        genome, gene, homolog = row
        if homolog not in ggpMat:
            ggpMat[homolog] = {}
        if genome not in ggpMat[homolog]:
            ggpMat[homolog][genome] = set()
        ggpMat[homolog][genome].add(gene)
    return ggpMat

def get_homolog_count_matrix(genomeIdSet, ggpMat = None, scope='global'):
    """ Given list of genome ids: 
        tabulate counts per genome per homology class 
        (formats data from get_homologs_for_genomes as table)
    """
    genomeGenePgfamList = get_homologs_for_genomes(genomeIdSet, scope=scope)
    prevMat = None
    if ggpMat: # if an existaing matrix was passed, extend it
        prevMat = copy.deepcopy(ggpMat)
    else:
        ggpMat = {} # genome-gene-homolog matrix (really just a dictionary)
    for row in genomeGenePgfamList:
        genome, gene, homolog = row
        if prevMat and homolog in prevMat and genome in prevMat[homolog]:
            continue # don't double-count
        if homolog not in ggpMat:
            ggpMat[homolog] = {}
        if genome not in ggpMat[homolog]:
            ggpMat[homolog][genome] = 0
        ggpMat[homolog][genome] += 1
    return ggpMat

def getHomologGenomeMatrix(genomes, homologs, ggpMat=None, scope='global'):
    """ Get homolog (eg PGFam) matrix specifying both homologs and genomes. """
    if Debug:
        LOG.write("patric_api.getHomologGenomeMatrix() called with {} genomes and {} homologs\n".format(len(genomes), len(homologs)))
    familyType = 'pgfam_id'
    if scope == 'local':
        familyType = 'plfam_id'
    if not ggpMat: # if an existaing matrix was passed, extend it
        ggpMat = {} # genome-gene-pgfam matrix (really just a dictionary)
    # consider breaking query up by slices of genomes or homologs or both, if needed
    query = "in(genome_id,({}))".format(",".join(genomes)) #, "in(product,(%s))"%",".join(roles))
    query += "&in({},({}))".format(familyType, ",".join(homologs))
    query += "&select(genome_id,patric_id,{})".format(familyType)
    query += "&limit(25000)"
    response = Session.get(Base_url+"genome_feature/", params=query) #, 
    if Debug:
        LOG.write("query= %s\n"%response.url)
    for line in response.text.split("\n")[1:]: # skip first line header
        line = line.replace('"','')
        if Debug:
            LOG.write("row= "+line+"\n")
        row = line.split("\t")
        if len(row) != 3:
            continue
        genome, gene, pgfam = row
        if pgfam not in ggpMat:
            ggpMat[pgfam] = {}
        if genome not in ggpMat[pgfam]:
            ggpMat[pgfam][genome] = set()
        ggpMat[pgfam][genome].add(gene)
    return ggpMat

def write_homolog_gene_matrix(ggpMat, fileHandle, minGenomes=0):
    """ write out homologGeneMatrix to file handle 
    data is list of genes per homolog per genome
    rows are homologs
    cols are genomes
    column headers identify genomes
    genes are comma-separated
    write only rows where > minProp*num_genomes genomes have genes
    if minProp is > 1, write rows where > minProp genomes have genes
    """
    # first collect set of all genomes
    genomeSet = set()
    for homolog in ggpMat:
        genomeSet.update(set(ggpMat[homolog].keys()))
    genomes = sorted(genomeSet)
    min_genomes_required =  minGenomes
    if minGenomes < 1:
        min_genomes_required = minGenomes * len(genomes)
    fileHandle.write("PGFam\t"+"\t".join(genomes)+"\n")
    for homolog in ggpMat:
        if len(ggpMat[homolog]) >= min_genomes_required:
            fileHandle.write(homolog)
            for genome in genomes:
                gene = ""
                if genome in ggpMat[homolog]:
                    gene = ",".join(ggpMat[homolog][genome])
                fileHandle.write("\t"+gene)
            fileHandle.write("\n")

def write_homolog_count_matrix(ggpMat, fileHandle, minGenomes=0):
    """ write out matrix of counts per homolog per genome to file handle 
    data is count of genes per homolog per genome (integers)
    rows are homologs
    cols are genomes
    column headers identify genomes
    write only rows where > minProp*num_genomes genomes have genes
    if minProp is > 1, write rows where > minProp genomes have genes
    """
    # first collect set of all genomes
    genomeSet = set()
    for homolog in ggpMat:
        genomeSet.update(set(ggpMat[homolog].keys()))
    genomes = sorted(genomeSet)
    min_genomes_required =  minGenomes
    if minGenomes < 1:
        min_genomes_required = minGenomes * len(genomes)
    fileHandle.write("PGFam\t"+"\t".join(genomes)+"\n")
    for homolog in ggpMat:
        if len(ggpMat[homolog]) >= min_genomes_required:
            fileHandle.write(homolog)
            for genome in genomes:
                count = 0
                if genome in ggpMat[homolog]:
                    count = len(ggpMat[homolog][genome])
                fileHandle.write("\t%d"%count)
            fileHandle.write("\n")

def read_homolog_gene_matrix(fileHandle):
    """ read homologGeneMatrix from file handle
    Data are list of genes (comma-delimited) per genome per homolog
    rows are homologs, cols are genomes, column headers identify genomes
    """
    # genome ids are headers in first line
    header = fileHandle.readline().rstrip()
    genomes = header.split("\t")[1:] # first entry is placeholder for homolog rownames
    pgMat = {} # genome-gene-homolog matrix (really just a dictionary)
    for row in fileHandle:
        fields = row.rstrip().split("\t")
        homolog = fields[0]
        pgMat[homolog] = {}
        data = fields[1:]
        for i, genome in enumerate(genomes):
            if len(data[i]):
                pgMat[homolog][genome] = data[i]
    return pgMat

def read_homolog_count_matrix(fileHandle):
    """ read homologCountMatrix from file handle
    rows are homologs
    cols are genomes
    data are integer counts of that homolog in that genome
    column headers identify genomes
    """
    # genome ids are headers in first line
    header = fileHandle.readline().rstrip()
    genomes = header.split("\t")[1:] # first entry is placeholder for homolog rownames
    pcMat = {} # homolog count matrix (really just a dictionary)
    for row in fileHandle:
        fields = row.rstrip().split("\t")
        homolog = fields[0]
        pcMat[homolog] = {}
        data = fields[1:]
        for i, genome in enumerate(genomes):
            pcMat[homolog][genome] = int(float(data[i]))
    return pcMat

def get_homologs_from_genome_object(genomeObject, scope='global'):
# parse a PATRIC genome object (read from json format) for PGFams
# scope can be 'global' or 'local' for pgfams or plfams
    retval = [] # a list of tupples of (genomeId, Pxfam, geneId)
    genomeId = genomeObject['id']
    target_family_type = ('PLFAM', 'PGFAM')[scope == 'global'] #test is index into alternatives
    for feature in genomeObject['features']:
        if 'family_assignments' in feature:
            for family_assignment in feature['family_assignments']:
                if family_assignment[0] == target_family_type:
                    retval.append((genomeId, feature['id'], family_assignment[1]))
    return retval

def getGenomeObjectProteins(genomeObject):
# return dictionary of patricId -> BioPython.SeqRecord
    genomeId = genomeObject['id']
    retval = {}
    for feature in genomeObject['features']:
        patricId, product, genomeId, aa_sequence = '', '', '', ''
        patricId = feature['id']
        if "protein_translation" in feature:
            aa_sequence = feature["protein_translation"]
            aa_sequence = aa_sequence.replace("J", "X") # because RAxML considers 'J' illegal
        if 'function' in feature:
            product = feature['function']
        simpleSeq = Seq(aa_sequence) #, IUPAC.extended_protein)
        seqRecord = SeqRecord(simpleSeq, id=patricId, description=product)
        seqRecord.annotations["genome_id"] = genomeId
        retval[patricId] = seqRecord
    return retval

def getGenomeObjectGeneDna(genomeObject):
# return dictionary of patricId -> BioPython.SeqRecord
    genomeId = genomeObject['id']
    contigSeq = {}
    for contig in genomeObject['contigs']:
        contigSeq[contig['id']] = contig['dna']
    retval = {} # dict of SeqRecords
    for feature in genomeObject['features']:
        geneId = feature['id']
	# ?? patricIds not defined here. Guessing this was a copy and paste issue. RDO.
        #if geneId not in patricIds:
        #    continue
        product = ''
        if 'product' in feature:
            product = feature['function']
        if not 'location' in feature:
            continue
        contig, start, ori, length = feature['location'][0] # this should be an array of (contig, start, orientation, length)
        start = int(float(start))
        length = int(float(length))
        if ori == '+':
            start -= 1
            simpleSeq = Seq(contigSeq[contig][start:start+length]) #, IUPAC.ambiguous_dna)
        if ori == '-':
            simpleSeq = Seq(contigSeq[contig][start-length:start]) #, IUPAC.ambiguous_dna)
            simpleSeq = simpleSeq.reverse_complement()

        seqRecord = SeqRecord(simpleSeq, id=geneId, description=product)
        seqRecord.annotations["genome_id"] = genomeId
        retval[geneId] = seqRecord
    return retval

