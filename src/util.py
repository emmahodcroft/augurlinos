import numpy as np

def generic_argparse(desc):
    import argparse
    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--path', default='',
                        help="Prefix for output files (e.g., 'zika' for output like 'auspice/zika_meta.json')")
    return parser


def write_fasta(seqs, fname, ungap=False):
    from Bio import SeqIO
    with open(fname, 'w') as ofile:
        for seq in seqs:
            seq.description=""
            seq.name=seq.id
            if ungap:
                seq.seq.ungap('-')
            seq.seq.upper()
            SeqIO.write(seq, ofile, format='fasta')


def read_sequence_meta_data(path):
    from filenames import meta_file_name
    import pandas as pd
    df = pd.read_csv(meta_file_name(path), sep='\t').fillna('')
    #reads blanks in TSV as blanks, which plays nice with TreeTime
    return {m[0]:m.to_dict() for mi, m in df.iterrows()}


def write_sequence_meta_data(path, df):
    from filenames import meta_file_name
    df.to_csv(meta_file_name(path), sep='\t', index=False)


def read_tree_meta_data(path):
    from filenames import tree_meta_file_name
    # import pandas as pd
    # df = pd.read_csv(tree_meta_file_name(path), sep='\t')
    # return {m[0]:m.to_dict() for mi, m in df.iterrows()}
    import json
    with open(tree_meta_file_name(path)) as ifile:
        meta = json.load(ifile)
    return meta


def write_tree_meta_data(path, meta_dic, indent=1):
    from filenames import tree_meta_file_name
    # import pandas as pd
    # header = set()
    # for m in meta_dic.values():
    #     header.update(m.keys())
    # header = sorted(header, key=lambda x:x!='name')
    # df = pd.DataFrame(meta_dic.values(), columns = header)
    # df.to_csv(tree_meta_file_name(path), sep='\t', index=False)
    import json
    with open(tree_meta_file_name(path), 'w') as ofile:
        json.dump(meta_dic, ofile, indent=indent)



def collect_tree_meta_data(T, fields, isvcf=False, meta=None):
    def mutation_format(muts):
        if isvcf:
            #converts from Python numbering [0] to standard [1] for output
            return ",".join(['%s%d%s'%(x[0], x[1]+1, x[2]) for x in muts])
        else:
            return ",".join(['%s%d%s'%(x[0], x[1], x[2]) for x in muts])

    if meta is None:
        meta = {}
    for n in T.find_clades():
        meta_dic = {'name':n.name}
        for field in fields:
            if hasattr(n,field):
                if 'mutations' in field:
                    meta_dic[field]=mutation_format(n.__getattribute__(field))
                else:
                    meta_dic[field]=n.__getattribute__(field)

        if n.name in meta:
            meta[n.name].update(meta_dic)
        else:
            meta[n.name] = meta_dic

    return meta

def read_in_vcf(vcf_file, ref_file):
    """
    Reads in a vcf/vcf.gz file and associated
    reference sequence fasta (to which the VCF file is mapped)

    Parses mutations, insertions, and deletions and stores them in a nested dict
    with the format:
    {'reference':'AGCTCGA..A',
     'sequences': { 'seq1':{4:'A', 7:'-'}, 'seq2':{100:'C'} },
     'insertions': { 'seq1':{4:'ATT'}, 'seq3':{1:'TT', 10:'CAG'} },
     'positions': [1,4,7,10,100...] }

    Calls with heterozygous values 0/1, 0/2, etc and no-calls (./.) are
    replaced with Ns at the associated sites. 
    
    Positions are stored to correspond the location in the reference sequence
    in Python (numbering is transformed to start at 0)

    Args
    ----
    vcf_file : string
        Path to the vcf or vcf.gz file to be read in
    ref_file : string
        Path to the fasta reference file to be read in

    Returns
    --------
    compress_seq : nested dict
        Contains the following keys:

        references : string
            String of the reference sequence read from the Fasta
        sequences : nested dict
            Dict containing sequence names as keys which map to dicts
            that have position as key and the single-base mutation (or deletion)
            as values
        insertions : nested dict
            Dict in the same format as the above, which stores insertions and their
            locations. The first base of the insertion is the same as whatever is
            currently in that position (Ref if no mutation, mutation in 'sequences'
            otherwise), so the current base can be replaced by the bases held here
            without losing that base.
        positions : list
            Python list of all positions with a mutation, insertion, or deletion.

    Note on VCF Format
    -------------------
    'Insertion where there are also deletions' (special handling)
        Ex:
          REF     ALT         Seq1    Seq2
          GC      GCC,G       1/1     2/2
        Insertions formatted differently - don't know how many bp match
        the Ref (unlike simple insert below). Could be mutations, also.
        
    'Deletion'
        Ex:
          REF     ALT 
          GC      G   
        Alt does not have to be 1 bp - any length shorter than Ref.
        
    'Insertion'
        Ex:
          REF     ALT 
          A       ATT
        First base always matches Ref.

    'No indel'
        Ex:
          REF     ALT 
          A       G
          
        
    EBH 4 Dec 2017
    """
    #define here, so that all sub-functions can access them
    sequences = {}
    insertions = {} #Currently not used, but kept in case of future use.
    positions = []
    
    #In future, if TreeTime handles 2-3 base ambig codes, this will allow that.
    def getAmbigCode(bp1, bp2, bp3=""):
        bps = [bp1,bp2,bp3]
        bps.sort()
        key = "".join(bps)

        return {
            'CT': 'Y',
            'AG': 'R',
            'AT': 'W',
            'CG': 'S',
            'GT': 'K',
            'AC': 'M',
            'AGT': 'D',
            'ACG': 'V',
            'ACT': 'H',
            'CGT': 'B'
        }[key]

    #Parses a 'normal' (not hetero or no-call) call depending if insertion+deletion, insertion,
    #deletion, or single bp subsitution
    def parseCall(pos, ref, alt):
        
        #Insertion where there are also deletions (special handling)
        if len(ref) > 1 and len(alt)>len(ref):
            if seq not in insertions.keys():
                insertions[seq] = {}
            for i in xrange(len(ref)):
                #if the pos doesn't match, store in sequences
                if ref[i] != alt[i]:
                    sequences[seq][pos+i] = alt[i] if alt[i] != '.' else 'N' #'.' = no-call
                #if about to run out of ref, store rest:
                if (i+1) >= len(ref):
                    insertions[seq][pos+i] = alt[i:]
        #Deletion
        elif len(ref) > 1:
            for i in xrange(len(ref)):
                #if ref is longer than alt, these are deletion positions
                if i+1 > len(alt):
                    sequences[seq][pos+i] = '-'
                #if not, there may be mutations
                else:
                    if ref[i] != alt[i]:
                        sequences[seq][pos+i] = alt[i] if alt[i] != '.' else 'N' #'.' = no-call
        #Insertion
        elif len(alt) > 1:
            if seq not in insertions.keys():
                insertions[seq] = {}
            insertions[seq][pos] = alt
        #No indel
        else:
            sequences[seq][pos] = alt
    
    
    #Parses a 'bad' (hetero or no-call) call depending on what it is
    def parseBadCall(pos, ref, ALT):
        #Deletion
        #   REF     ALT     Seq1    Seq2    Seq3
        #   GCC     G       1/1     0/1     ./.
        # Seq1 (processed by parseCall, above) will become 'G--'
        # Seq2 will become 'GNN'
        # Seq3 will become 'GNN'
        if len(ref) > 1:
            #Deleted part becomes Ns
            if seq not in sequences.keys():
                sequences[seq] = {}
                
            if gen[0] == '0' or gen[0] == '.':
                if gen[0] == '0':   #if het, get first bp
                    alt = str(ALT[int(gen[2])-1])
                else: #if no-call, there is no alt, so just put Ns after 1st ref base
                    alt = ref[0]
                for i in xrange(len(ref)):
                    #if ref is longer than alt, these are deletion positions
                    if i+1 > len(alt):
                        sequences[seq][pos+i] = 'N'
                    #if not, there may be mutations
                    else:
                        if ref[i] != alt[i]:
                            sequences[seq][pos+i] = alt[i] if alt[i] != '.' else 'N' #'.' = no-call

        #If not deletion, need to know call type
        #if het, see if proposed alt is 1bp mutation
        elif gen[0] == '0':
            alt = str(ALT[int(gen[2])-1])
            if len(alt)==1:
                #alt = getAmbigCode(ref,alt) #if want to allow ambig
                alt = 'N' #if you want to disregard ambig
                if seq not in sequences.keys():
                    sequences[seq] = {}
                sequences[seq][pos] = alt
            #else a het-call insertion, so ignore.

        #else it's a no-call; see if all alts have a length of 1
        #(meaning a simple 1bp mutation)
        elif len(ALT)==len("".join(ALT)):
            alt = 'N'
            if seq not in sequences.keys():
                sequences[seq] = {}
            sequences[seq][pos] = alt
        #else a no-call insertion, so ignore.
 
        
    #House code is *much* faster than pyvcf because we don't care about all info
    #about coverage, quality, counts, etc, which pyvcf goes to effort to parse
    #(and it's not easy as there's no standard ordering). Custom code can completely
    #ignore all of this.
    import gzip
    from Bio import SeqIO
    import numpy as np

    nsamp = 0
    posLoc = 0
    refLoc = 0
    altLoc = 0
    sampLoc = 9

    #Use different openers depending on whether compressed
    opn = gzip.open if vcf_file.endswith(('.gz', '.GZ')) else open 

    with opn(vcf_file, mode='rt') as f:
        for line in f:
            if line[0] != '#':
                #actual data - most common so first in 'if-list'!
                line = line.strip()
                dat = line.split('\t')
                POS = int(dat[posLoc])
                REF = dat[refLoc]
                ALT = dat[altLoc].split(',')
                calls = np.array(dat[sampLoc:])

                #get samples that differ from Ref at this site
                recCalls = {}
                for sname, sa in zip(samps, calls):
                    if ':' in sa: #if proper VCF file (followed by quality/coverage info)
                        gt = sa.split(':')[0]
                    else: #if 'pseudo' VCF file (nextstrain output, or otherwise stripped)
                        gt = sa
                    if '/' in gt and gt != '0/0':  #ignore if ref call: '.' or '0/0', depending on VCF
                        recCalls[sname] = gt

                #store the position and the alt
                for seq, gen in recCalls.iteritems():
                    ref = REF
                    pos = POS-1     #VCF numbering starts from 1, but Reference seq numbering
                                    #will be from 0 because it's python!
                                    
                    #Accepts only calls that are 1/1, 2/2 etc. Rejects hets and no-calls
                    if gen[0] != '0' and gen[2] != '0' and gen[0] != '.' and gen[2] != '.':
                        alt = str(ALT[int(gen[0])-1])   #get the index of the alternate
                        if seq not in sequences.keys():
                            sequences[seq] = {}
                            
                        parseCall(pos, ref, alt)
                        
                    #If is heterozygote call (0/1) or no call (./.)
                    else:
                        #alt will differ here depending on het or no-call, must pass original
                        parseBadCall(pos, ref, ALT)

            elif line[0] == '#' and line[1] == 'C':
                #header line, get all the information
                header = line.strip().split('\t')
                posLoc = header.index("POS")
                refLoc = header.index('REF')
                altLoc = header.index('ALT')
                sampLoc = header.index('FORMAT')+1
                samps = header[sampLoc:]
                nsamp = len(samps)

            #else you are a comment line, ignore.
            
    #Get list of all mutation positions       
    for seq, muts in sequences.iteritems():
        positions = positions + list( set(muts.keys()) - set(positions) )

    #One or more seqs are same as ref! (No non-ref calls) So haven't been 'seen' yet   
    if nsamp > len(sequences): 
        missings = set(samps).difference(sequences.keys())
        for s in missings:
            sequences[s] = {}

    refSeq = SeqIO.parse(ref_file, format='fasta').next()
    refSeqStr = str(refSeq.seq)

    compress_seq = {'reference':refSeqStr,
                    'sequences': sequences,
                    'insertions': insertions,
                    'positions': positions}

    return compress_seq

def read_in_translate_vcf(vcf_file, ref_file, compressed=True):
    """
    Reads in a vcf file where TRANSLATIONS have been stored and associated
    reference sequence fasta (to which the VCF file is mapped)
    This is the file output by "write_VCF_translation" below
    
    Very simple compared to the above as will never be insertion or deletion

    Returns a nested dict in the same format as is *input* in "write_VCF_translation" below,
    with a nested dict for each gene, which contains 'sequences', 'positions', and 'reference'
    """
    from Bio import SeqIO
    import numpy as np

    prots = {}

    posLoc = 0
    refLoc = 0
    altLoc = 0
    sampLoc = 9

    with open(vcf_file) as f:
        for line in f:
            if line[0] != '#':
                #actual data
                line = line.strip()
                dat = line.split('\t')
                POS = int(dat[posLoc])
                REF = dat[refLoc]
                ALT = dat[altLoc].split(',')
                GEN = dat[0] #'CHROM' or the gene name here
                calls = np.array(dat[sampLoc:])

                #get samples that differ from Ref at this site
                recCalls = {}
                for sname, sa in zip(samps, calls):
                    if sa != '.':
                        recCalls[sname] = sa

                #store position and the altLoc
                for seq, gen in recCalls.iteritems():
                    alt = str(ALT[int(gen[0])-1])   #get the index of the alternate
                    ref = REF
                    pos = POS-1     #VCF numbering starts from 1, but Reference seq numbering
                                    #will be from 0 because it's python!
                    gen = GEN       #from CHROM, gene name

                    if gen not in prots.keys():
                        prots[gen] = {}
                        prots[gen]['sequences'] = {}
                        prots[gen]['positions'] = []
                        prots[gen]['reference'] = ''
                    if seq not in prots[gen]['sequences'].keys():
                        prots[gen]['sequences'][seq] = {}

                    #will never be insertion or deletion! because translation.
                    prots[gen]['sequences'][seq][pos] = alt
                    prots[gen]['positions'].append(pos)

            elif line[0] == '#' and line[1] == 'C':
                #header line, get all the information
                header = line.strip().split('\t')
                posLoc = header.index("POS")
                refLoc = header.index("REF")
                altLoc = header.index("ALT")
                sampLoc = header.index("FORMAT")+1
                samps = header[sampLoc:]
                nsamp = len(samps)

    for refSeq in SeqIO.parse(translation_ref_file(path), format='fasta'):
        prots[refSeq.name]['reference'] = str(refSeq.seq)
        posN = np.array(prots[refSeq.name]['positions'])
        posN = np.unique(posN)
        prots[refSeq.name]['positions'] = np.sort(posN)

    return prots


def read_in_DRMs(drm_file):
    """
    Reads in a tab-delim file giving information on Drug Resistance Mutations.
    Format at the moment must include columns titled:
    GENOMIC_POSITION (bp of the mutation)
    ALT_BASE    (the alt-base conferring resistance)
    SUBSTITUTION    (the aminoacid-codonloc-aminoacid change)
    DRUG    (the drug of resistance)
    GENE    (the gene the mutation is in)

    Other columns can be present but will be ignored

    Returns a nested dict containing:
        drmPositions: an array of all the positions of DRMs

        DRMs: another nested dict with information about the mutations
            'position' is the key, and maps to a nested dict that contains:
            base, drug, AA, gene
            base and AA are arrays as one bp position can have multiple 'alt' bases

    """
    import pandas as pd
    import numpy as np

    DRMs = {}
    drmPositions = []

    df = pd.read_csv(drm_file, sep='\t')
    for mi, m in df.iterrows():
        pos = m.GENOMIC_POSITION-1 #put in python numbering
        drmPositions.append(pos)

        if pos in DRMs:
            DRMs[pos]['base'].append(m.ALT_BASE)
            DRMs[pos]['AA'].append(m.SUBSTITUTION)
        else:
            DRMs[pos] = {}
            DRMs[pos]['base'] = [m.ALT_BASE]
            DRMs[pos]['drug'] = m.DRUG
            DRMs[pos]['AA'] = [m.SUBSTITUTION]
            DRMs[pos]['gene'] = m.GENE

    drmPositions = np.array(drmPositions)
    drmPositions = np.unique(drmPositions)
    drmPositions = np.sort(drmPositions)

    DRM_info = {'DRMs': DRMs,
            'drmPositions': drmPositions}

    return DRM_info


#####################################################
# date parsing and conversions
#####################################################
def parse_date(datein, fmt, max_min_year=None):
    from datetime import datetime
    import numpy as np
    try:
        if 'XX' in datein:
            min_date, max_date = ambiguous_date_to_date_range(datein, fmt, max_min_year)
            n_date = np.array((numerical_date(min_date), numerical_date(max_date)))
        else:
            tmp = datetime.strptime(datein, fmt).date()
            n_date = numerical_date(tmp)
    except:
        print("Can't parse ",datein)
        n_date=None

    return n_date


def numerical_date(date):
    from datetime import datetime
    days_in_year = date.toordinal()- datetime(year=date.year, month=1, day=1).date().toordinal()
    return date.year + days_in_year/365.25


def ambiguous_date_to_date_range(mydate, fmt, max_min_year=None):
    #Can now be passed a 1- or 2-element list containing the minimum or minimum and maximum
    #dates to be used in the range.
    from datetime import datetime
    sep = fmt.split('%')[1][-1]
    min_date, max_date = {}, {}
    today = datetime.today().date()

    for val, field  in zip(mydate.split(sep), fmt.split(sep+'%')):
        f = 'year' if 'y' in field.lower() else ('day' if 'd' in field.lower() else 'month')
        if 'XX' in val:
            if f=='year':
                if max_min_year:
                    min_date[f]=max_min_year[0]
                    if len(max_min_year)>1:
                        max_date[f]=max_min_year[1]
                    elif len(max_min_year)==1:
                        max_date[f]=4000 #will be replaced by 'today' below.
                else:
                    return None, None
            elif f=='month':
                min_date[f]=1
                max_date[f]=12
            elif f=='day':
                min_date[f]=1
                max_date[f]=31
        else:
            min_date[f]=int(val)
            max_date[f]=int(val)
    max_date['day'] = min(max_date['day'], 31 if max_date['month'] in [1,3,5,7,8,10,12]
                                           else 28 if max_date['month']==2 else 30)
    lower_bound = datetime(year=min_date['year'], month=min_date['month'], day=min_date['day']).date()
    upper_bound = datetime(year=max_date['year'], month=max_date['month'], day=max_date['day']).date()
    return (lower_bound, upper_bound if upper_bound<today else today)

##########################################
# IO
##########################################

def write_json(data, file_name, indent=1):
    import json
    import os

    #in case auspice folder does not exist yet
    if not os.path.exists(os.path.dirname(file_name)):
        try:
            os.makedirs(os.path.dirname(file_name))
        except OSError: #Guard against race condition
            if not os.path.isdir(os.path.dirname(file_name)):
                raise
    try:
        handle = open(file_name, 'w')
    except IOError:
        raise
    else:
        json.dump(data, handle, indent=indent)
        handle.close()

def write_VCF_style_alignment(tree_dict, file_name, compress=False):
    """
    Writes out a VCF-style file (which seems to be minimally handleable
    by vcftools and pyvcf) of the alignment from the input of a dict
    in a similar format to what's created from the read_in_vcf function above.
    
    For a sequence like:
    Pos     1 2 3 4 5 6
    Ref     A C T T A C
    Seq1    A C - - - G
    
    In a dict it is stored:
    Seq1:{3:'-', 4:'-', 5:'-', 6:'G'}  (Numbering from 1 for simplicity)
    
    In a VCF it needs to be:
    POS REF     ALT     Seq1
    2   CTTA    C       1/1
    6   C       G       1/1
    
    If a position is deleted (pos 3), need to get invariable position preceeding it
    
    However, in alternative case, the base before a deletion is mutant, so need to check
        that next position isn't a deletion (as otherwise won't be found until after the
        current single bp mutation is written out)
        
    When deleted position found, need to gather up all adjacent mutant positions with deletions,
        but not include adjacent mutant positions that aren't deletions (pos 6)
        
    Don't run off the 'end' of the position list if deletion is the last thing to be included
        in the VCF file

    EBH 7 Dec 2017
    """
    sequences = tree_dict['sequences']
    ref = tree_dict['reference']
    positions = tree_dict['positions']
    
    def handleDeletions(i, pi, pos, ref, delete, pattern):
        refb = ref[pi]
        if delete: #Need to get the position before
            i-=1    #As we'll next go to this position again
            pi-=1
            pos = pi+1
            refb = ref[pi]
            #re-get pattern
            pattern = []
            for k,v in sequences.iteritems():
                try:
                    pattern.append(sequences[k][pi])
                except KeyError, e:
                    pattern.append(ref[pi])
            pattern = np.array(pattern)

        sites = []
        sites.append(pattern)

        #Gather all positions affected by deletion - but don't run off end of position list
        while (i+1) < len(positions) and positions[i+1] == pi+1:
            i+=1
            pi = positions[i]
            pattern = []
            for k,v in sequences.iteritems():
                try:
                    pattern.append(sequences[k][pi])
                except KeyError, e:
                    pattern.append(ref[pi])
            pattern = np.array(pattern)

            #Stops 'greedy' behaviour from adding mutations adjacent to deletions 
            if any(pattern == '-'): #if part of deletion, append
                sites.append(pattern)
                refb = refb+ref[pi]
            else: #this is another mutation next to the deletion!
                i-=1    #don't append, break this loop

        #Rotate them into 'calls'
        sites = np.asarray(sites)
        align = np.rot90(sites)
        align = np.flipud(align)

        #Get rid of '-', and put '.' for calls that match ref
        #Only removes trailing '-'. This breaks VCF convension, but the standard
        #VCF way of handling this* is really complicated, and the situation is rare.
        #(*deletions and mutations at the same locations)
        fullpat = []
        for pt in align:
            gp = len(pt)-1
            while pt[gp] == '-':
                pt[gp] = ''
                gp-=1
            pat = "".join(pt)
            if pat == refb:
                fullpat.append('.')
            else:
                fullpat.append(pat)

        pattern = np.array(fullpat)

        return i, pi, pos, refb, pattern

        
    #prepare the header of the VCF & write out
    header=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+sequences.keys()
    with open(file_name, 'w') as the_file:
        the_file.write( "##fileformat=VCFv4.2\n"+
                        "##source=NextStrain\n"+
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        the_file.write("\t".join(header)+"\n")

    vcfWrite = []
    errorPositions = []
    explainedErrors = 0

    #Why so basic? Because we sometimes have to back up a position!
    i=0
    while i < len(positions):
        #Get the 'pattern' of all calls at this position.
        #Look out specifically for current (this pos) or upcoming (next pos) deletions
        #But also distinguish these two, as handled differently.
        
        pi = positions[i]
        pos = pi+1 #change numbering to match VCF, not python, for output
        refb = ref[pi] #reference base at this position
        delete = False #deletion at this position - need to grab previous base (invariable)
        deleteGroup = False #deletion at next position (mutation at this pos) - do not need to get prev base

        #try/except is much more efficient than 'if' statements for constructing patterns,
        #as on average a 'variable' location will not be variable for any given sequence
        pattern = []
        #pattern2 gets the pattern at next position to check for upcoming deletions
        #it's more efficient to get both here rather than loop through sequences twice!
        pattern2 = [] 
        for k,v in sequences.iteritems():
            try:
                pattern.append(sequences[k][pi])
            except KeyError, e:
                pattern.append(ref[pi])

            try:
                pattern2.append(sequences[k][pi+1])
            except KeyError, e:
                pattern2.append(ref[pi+1])

        pattern = np.array(pattern)
        pattern2 = np.array(pattern2)

        #If a deletion here, need to gather up all bases, and position before
        if any(pattern == '-'):
            if pos != 1:
                deleteGroup = True
                delete = True
            else:
                #If theres a deletion in 1st pos, VCF files do not handle this well.
                #Proceed keeping it as '-' for alt (violates VCF), but warn user to check output.
                #(This is rare)
                print "WARNING: You have a deletion in the first position of your alignment. VCF format does not handle this well. Please check the output to ensure it is correct."
        else:
            #If a deletion in next pos, need to gather up all bases
            if any(pattern2 == '-'):
                deleteGroup = True

        #If deletion, treat affected bases as 1 'call':
        if delete or deleteGroup:
            i, pi, pos, refb, pattern = handleDeletions(i, pi, pos, ref, delete, pattern)
        #If no deletion, replace ref with '.', as in VCF format
        else: 
            pattern[pattern==refb] = '.'

        #Get the list of ALTs - minus any '.'!
        uniques = np.unique(pattern)
        uniques = uniques[np.where(uniques!='.')]

        #Convert bases to the number that matches the ALT
        j=1
        for u in uniques:
            pattern[np.where(pattern==u)[0]] = str(j)
            j+=1
        #Now convert these calls to #/# (VCF format)
        calls = [ j+"/"+j if j!='.' else '.' for j in pattern ]

        #What if there's no variation at a variable site??
        #This can happen when sites are modified by TreeTime - see below.
        printPos = True
        if len(uniques)==0:
            #If we expect it (it was made constant by TreeTime), it's fine.
            if 'inferred_const_sites' in tree_dict and pi in tree_dict['inferred_const_sites']:
                explainedErrors += 1
                printPos = False #and don't output position to the VCF
            else:
                #If we don't expect, raise an error
                errorPositions.append(str(pi))

        #Write it out - Increment positions by 1 so it's in VCF numbering
        #If no longer variable, and explained, don't write it out
        if printPos:  
            output = ["MTB_anc", str(pos), ".", refb, ",".join(uniques), ".", "PASS", ".", "GT"] + calls
            vcfWrite.append("\t".join(output))

        i+=1

    #Note: The number of 'inferred_const_sites' passed back by TreeTime will often be longer
    #than the number of 'site that were made constant' that prints below. This is because given the site:
    # Ref   Alt     Seq
    # G     A       AANAA
    #This will be converted to 'AAAAA' and listed as an 'inferred_const_sites'. However, for VCF
    #purposes, because the site is 'variant' against the ref, it is variant, as expected, and so
    #won't be counted in the below list, which is only sites removed from the VCF.

    if 'inferred_const_sites' in tree_dict and explainedErrors != 0:
        print "Sites that were constant except for ambiguous bases were made constant by TreeTime. This happened {} times. These sites are now excluded from the VCF.".format(explainedErrors)

    if len(errorPositions) != 0:
        print ("\n***UNEXPECTED ERROR!! util.py: write_VCF_style_alignment"
            "\n{} sites were found that had no alternative bases. If this data has been "
            "run through TreeTime and contains ambiguous bases, try calling get_tree_dict with "
            "var_ambigs=True to see if this clears the error."
            "\nIf it does not, or if the data hasn't come from TreeTime, something unexpected "
            "is happening and you should debug. In TreeTime, this can be caused by overwriting "
            "variants in tips with small branch lengths."
            "\nThese are the positions affected (numbering starts at 0):").format(str(len(errorPositions)))
        print ",".join(errorPositions)

    with open(file_name, 'a') as the_file:
        the_file.write("\n".join(vcfWrite))

    if compress:
        import os
        call = ["gzip", file_name]
        os.system(" ".join(call))



########################################
# translation
#######################################3
nuc_alpha = 'ACGT-N'
aa_alpha = 'ACDEFGHIKLMNPQRSTVWY*-X'
TINY = 1e-12

def safe_translate(sequence, report_exceptions=False):
    """Returns an amino acid translation of the given nucleotide sequence accounting
    for gaps in the given sequence.

    Optionally, returns a tuple of the translated sequence and whether an
    exception was raised during initial translation.

    >>> safe_translate("ATG")
    'M'
    >>> safe_translate("ATGGT-")
    'MX'
    >>> safe_translate("ATG---")
    'M-'
    >>> safe_translate("ATGTAG")
    'M*'
    >>> safe_translate("")
    ''
    >>> safe_translate("ATGT")
    'M'
    >>> safe_translate("ATG", report_exceptions=True)
    ('M', False)
    >>> safe_translate("ATGA-G", report_exceptions=True)
    ('MX', True)
    """
    from Bio.Seq import Seq, CodonTable
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import generic_dna
    from Bio.Data.CodonTable import TranslationError
    translation_exception = False

    try:
        # Attempt translation by extracting the sequence according to the
        # BioPhython SeqFeature in frame gaps of three will translate as '-'
        translated_sequence = str(Seq(sequence).translate(gap='-'))
    except TranslationError:
        translation_exception = True

        # Any other codon like '-AA' or 'NNT' etc will fail. Translate codons
        # one by one.
        codon_table  = CodonTable.ambiguous_dna_by_name['Standard'].forward_table
        str_seq = str(sequence)
        codons = np.fromstring(str_seq[:len(str_seq) - len(str_seq) % 3], dtype='S3')
        assert len(codons) > 0
        aas = []

        for c in codons:
            # Parse result of single codon translation, add amino acids as
            # appropriate.
            try:
                aa = codon_table.get(c)
                if aa is None:
                    if c == '---':
                        aas.append('-')
                    else:
                        aas.append('X')
                else:
                    aas.append(aa)
            except (TranslationError, ValueError):
                aas.append('X')

        translated_sequence = "".join(aas)

    if report_exceptions:
        return translated_sequence, translation_exception
    else:
        return translated_sequence


### get genes and translations
def get_genes_and_alignments(path, tree=True):
    import os.path
    from filenames import translation_ref_file, translation_vcf_file

    genes = []

    if os.path.isfile(translation_ref_file(path)): #vcf file, use different method!
        from Bio import SeqIO
        for seq in SeqIO.parse(translation_ref_file(path), format='fasta'):
            genes.append((seq.name, translation_vcf_file(path)))

    else:
        import glob
        if tree:
            from filenames import tree_sequence_alignment as func
        else:
            from filenames import ref_alignment as func

        mask = func(path, prot='*')
        aln_files = glob.glob(mask)
        for aln_fname in aln_files:
            gene = aln_fname.rstrip(mask.split("*")[-1]).lstrip(mask.split('*')[0])
            genes.append((gene, aln_fname))

    return genes


def calc_af(aln_array, alpha):
    af = np.zeros((len(alpha), aln_array.shape[1]))
    for ai, state in enumerate(alpha):
        af[ai] += (aln_array==state).mean(axis=0)
    af[-1] = 1.0 - af[:-1].sum(axis=0)
    return af

def diversity_statistics(fname, nuc=True):
    from Bio import AlignIO
    aln_array = np.array(AlignIO.read(fname, 'fasta'))
    af = calc_af(aln_array, nuc_alpha if nuc else aa_alpha)
    tmp_af = af[:-2]/(af[:-2].sum(axis=0)+TINY)
    entropy = -(tmp_af*np.log(tmp_af+TINY)).sum(axis=0)

    return entropy


def load_features(reference, feature_names=None):
    #read in appropriately whether GFF or Genbank
    #checks explicitly for GFF otherwise assumes Genbank
    features = {}

    if '.gff' in reference.lower():
        #looks for 'gene' and 'gene' as best for TB
        from BCBio import GFF
        limit_info = dict( gff_type = ['gene'] )

        in_handle = open(reference)
        for rec in GFF.parse(in_handle, limit_info=limit_info):
            for feat in rec.features:
                if "gene" in feat.qualifiers:
                    fname = feat.qualifiers["gene"][0]
                else:
                    fname = feat.qualifiers["locus_tag"][0]
                if feature_names is None or fname in feature_names:
                    features[fname] = feat

        if feature_names is not None:
            for fe in feature_names:
                if fe not in features:
                    print "Couldn't find gene {} in GFF or GenBank file".format(fe)

        in_handle.close()

    else:
        from Bio import SeqIO
        for feat in SeqIO.read(reference, 'genbank').features:
            if feat.type=='CDS':
                if "locus_tag" in feat.qualifiers:
                    fname = feat.qualifiers["locus_tag"][0]
                    if feature_names is None or fname in feature_names:
                        features[fname] = feat

    return features


def write_VCF_translation(prot_dict, vcf_file_name, ref_file_name, compress=False):
    """
    Writes out a VCF-style file (which seems to be minimally handleable
    by vcftools and pyvcf) of the AA differences between sequences and the reference.
    This is a similar format created/used by read_in_vcf except that there is one
    of these dicts (with sequences, reference, positions) for EACH gene.

    Also writes out a fasta of the reference alignment.

    EBH 12 Dec 2017
    """

    #for the header
    seqNames = prot_dict[prot_dict.keys()[0]]['sequences'].keys()

    #prepare the header of the VCF & write out
    header=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+seqNames
    with open(vcf_file_name, 'w') as the_file:
        the_file.write( "##fileformat=VCFv4.2\n"+
                        "##source=NextStrain_Protein_Translation\n"+
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        the_file.write("\t".join(header)+"\n")

    refWrite = []
    vcfWrite = []

    #go through for every gene/protein
    for fname, prot in prot_dict.iteritems():
        sequences = prot['sequences']
        ref = prot['reference']
        positions = prot['positions']

        #write out the reference fasta
        refWrite.append(">"+fname)
        refWrite.append(ref)

        #go through every variable position
        #There are no deletions here, so it's simpler than for VCF nuc sequenes!
        i=0
        while i < len(positions):
            pi = positions[i]
            pos = pi+1 #change numbering to match VCF not python
            refb = ref[pi] #reference base at this position

            pattern = np.array([ sequences[k][pi] if pi in sequences[k].keys() else '.' for k,v in sequences.iteritems() ])

            #get the list of ALTs - minus any '.'!
            uniques = np.unique(pattern)
            uniques = uniques[np.where(uniques!='.')]

            #Convert bases to the number that matches the ALT
            j=1
            for u in uniques:
                pattern[np.where(pattern==u)[0]] = str(j)
                j+=1
            #Now convert these calls to #/# (VCF format)
            calls = [ j+"/"+j if j!='.' else '.' for j in pattern ]
            if len(uniques)==0:
                print "UNEXPECTED ERROR WHILE CONVERTING TO VCF AT POSITION {}".format(str(pi))
                break

            #put it all together and write it out!
            #increment positions by 1 so it's in VCF numbering not python numbering
            output = [fname, str(pos), ".", refb, ",".join(uniques), ".", "PASS", ".", "GT"] + calls

            vcfWrite.append("\t".join(output))

            i+=1

    #write it all out
    with open(ref_file_name, 'w') as the_file:
        the_file.write("\n".join(refWrite))

    with open(vcf_file_name, 'a') as the_file:
        the_file.write("\n".join(vcfWrite))

    if compress:
        import os
        call = ["gzip", vcf_file_name]
        os.system(" ".join(call))

def load_lat_long_defs():
    places = {}
    with open('../fauna/source-data/geo_lat_long.tsv', 'r') as latlongfile:
        header = latlongfile.readline().strip().split('\t')
        for line in latlongfile:
            try:
                place, country, latitude, longitude = line.strip().split('\t')
                places[place] = {'country_code':country,
                                 'latitude':float(latitude),
                                 'longitude':float(longitude)}
            except:
                print("trouble parsing", line)
    return places