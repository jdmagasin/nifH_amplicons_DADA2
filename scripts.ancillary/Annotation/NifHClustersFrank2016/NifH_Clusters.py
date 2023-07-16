
## Copyright (C) 2016 Ildiko E. Frank

#-----------------------------------------------------------------------------------------------------------
# Python script to classify NifH protein sequences into main clusters and
# subclusters based on CART models [Frank et al. 2016]. It requires the
# Biopython Module.
#
# Command:  python NifH_Clusters.py FILE PULL
#    e.g.   python NifH_Clusters.py YogevPythonIN 3E
#
# Input parameters:
# FILE:     name of input fasta file of aligned NifH protein sequences 
#            first sequence must be NifH from Azotobacter vinelandii,
#            specifically it must be identical to WP_012698955.1 (but
#            need not be the entire sequence).
#
# Dropped!  POSITION: starting position in aligned NifH sequences according to
#           Azotobacter vinelandii e.g. STRLIL...  sequence starts at position
#           45
#
# PULL: optional input to specify which predicted main cluster or subcluster
#           records should be written on the output file if PULL is omitted then
#           all records are written on the output file
#
# Output:
# Fasta file that has the same name as the input FILE with _Clusters extension;
# e.g. YogevPythonIN_Clusters.fasta.  Each record contains the input fasta
# description augmented with main cluster and subcluster prediction.  If both
# primary and surrogate positions are missing then the cluster prediction is NA
#
# Changes by jmagasin:
#  - Small changes for swtich from Python 2 --> 3
#
#  - No longer require the input multiple alignment (A2M format) to be in the
#    script directory.
#
#  - FILE must be named as the operating system sees it. Don't omit the ".fasta"
#  - (if any).
#
#  - Check that the first sequence appears to be NifH from A.vinelandii,
#    specifically WP_012698955.1. That enables the script to identify the
#    STRLI residues so we no longer need the POSITION parameter (which used
#    to be after FILE).
#
#  - If the A.vinelandii sequence is predicted to be other than cluster 1 and
#    subcluster 1G, give a big warning.
#
#-----------------------------------------------------------------------------------------

from Bio import SeqIO
import sys

# jmagasin: Use the A.vinelandii NifH sequence WP_012698955.1 for identifying
# the column in the multiple alignment that corresponds to residue 45 at STRLI.
refAzoNifH = "MALRQCAIYGKGGIGKSTTTQNLVAALAEAGKKVMIVGCDPKADSTRLILHSKAQGTVMEMAASAGSVEDLELEDVLQIGFGGVKCVESGGPEPGVGCAGRGVITAINFLEEEGAYSDDLDFVFYDVLGDVVCGGFAMPIRENKAQEIYIVCSGEMMAMYAANNIAKGIVKYAHSGSVRLGGLICNSRKTDREDELIMALAAKIGTQMIHFVPRDNVVQHAEIRRMTVIEYDPKAGQADEYRALARKIVDNKLLVIPNPASMEELEELLMEFGIMEVEDESVVGKAAAEG"
assert(refAzoNifH.find("STRLILHSKA")+1 == 45), "Bug: Reference A.vinelandii NifH lacks STRLI at 45."

# jmagasin: First sequence is required to be a NifH from A.vinelandii so that
# the script can determine how to map multiple alignment positions to NifH
# residues in the CART model.  This function does a simple check for
# A.vinelandii NifH by searching for STRLI... in refAzoNifH.  Then it returns
# the offset that decider() must use to map Azo residues to string positions.
def GetAzoResidueToStringMap(seq_in):
    strliPos = seq_in.find("STRLILHSKA")
    emsg = "The first sequence does not appear to be NifH from A.vinelandii (WP_012698955.1) since STRLI could not be found."
    assert (strliPos != -1), emsg
    return 45-strliPos

# List non-relevant positions based on Azotobacter sequence
def censor_maker(seq_in):
    L_out = []
    for i in range(len(seq_in)):
        if seq_in[i] == "-":
            L_out.append(i)
    return L_out

# Compress sequence to be classified, i.e. get rid of non-relevant positions 
def ref_edit(seq_in, pos_list): 
    seq_out = ""
    for i in range(len(seq_in)):
        if i not in pos_list:
            seq_out += seq_in[i]
    return seq_out

# Evaluate primary and surrogate positions in specified decision node of specified tree
# azoStart is determined near "Main routine"
def decider(seq_in, Node, icluster, inode):
    p1 = Node[0] - residuePosToStringPos # jmagasin: Subtract the num Azo residues that
    p2 = Node[2] - residuePosToStringPos # are not in the multiple alignment.
    n = len(seq_in) - 1  
    if (p1 <= n):  
        if seq_in[p1] in Node[1]:  
# primary exists, go left
            return True
        elif seq_in[p1] in ['-']:  
# primary missing 
            if p2 <= n:
               if seq_in[p2] in Node[3]:  
# surrogate exists, go left
                  return True
               elif seq_in[p2] in ['-']:
                   print("Primary missing and Surrogate missing in sequence " + record.description, repr(icluster), repr(inode))
                   return "NA"
               else: 
# surrogate exists, go right
                   return False
            else:
                print("Primary missing and Surrogate out of range in sequence " + record.description,  repr(icluster), repr(inode))
                return "NA"
        else:  
# primary exists, go right
            return False
    elif p1 > n and p2 <= n:  
# primary out of range
        if seq_in[p2] in Node[3]: 
# surrogat exists, go left
            return True
        elif seq_in[p2] in ['-']:
            print("Primary out of range and Surrogate missing in sequence " + record.description, repr(icluster),  repr(inode))
            return "NA"
        else:  
# surrogate exists, go right
            return False
    else:
        print("Primary and Surrogate out of range in sequence " + record.description, repr(icluster), repr(inode))
        return "NA"
        
# Classify sequences into main clusters: 1 or 2 or 3 or 4
def primary_tree(seq_in):
    Node1 =  [109, ['F', 'W', 'Y', 'S'], 52, ['A', 'C', 'F', 'L', 'M' 'P', 'S', 'T', 'Y']]
    Node2 =  [49, ['A', 'D', 'I'], 110, ['C', 'M', 'T', 'V']]
    Node3 =  [53, ['L', 'M', 'W'], 106, ['L', 'P', 'S']]

    #jmagasin: Next line helpful for debugging.
    #print("primary_tree is looking at " + seq_in )
    node1_val = decider(seq_in, Node1, 0, 1)
    if node1_val == "NA":
        return "NA"
    node2_val = decider(seq_in, Node2, 0, 2)
    if node2_val == "NA":
        return "NA"
    node3_val = decider(seq_in, Node3, 0, 3)
    if node3_val == "NA":
        return "NA"
         
    if node1_val:
        return "1"
    elif node2_val:
        return "2"
    elif node3_val:
        return "3"
    else:
        return "4"

# Classify sequences within cluster 1: 12 categories, 1A, 1B, etc.
def fork_1(seq_in):
    Node1 = [63, ['H', 'L', 'R'], 60, ['D', 'V']]
    Node2 = [62, ['F', 'I', 'M', 'V'], 78, ['C', 'K', 'N', 'R', 'W']]
    Node3 = [124, ['T'], 59, ['I']]
    Node4 = [65, ['F', 'I', 'K', 'L', 'R', 'W'], 78, ['A', 'C', 'D', 'F', 'G', 'I', 'L', 'P', 'Q']]
    Node5 = [79, ['M', 'V'], 122, ['Y']]
    Node6 = [56, ['A', 'P', 'T', 'V'], 52, ['A', 'C', 'F', 'G', 'N', 'R', 'S', 'Y']]
    Node7 = [124, ['A', 'C', 'G', 'P', 'S', 'V'], 56, ['A', 'C', 'D', 'G', 'P', 'Q', 'S']]
    Node8 = [56, ['C', 'E', 'N', 'S', 'T'], 59, ['A', 'D', 'G', 'W']]
    Node9 = [58, ['A', 'E', 'M', 'V'], 77, ['L', 'P', 'S', 'V']]
    Node10 = [61, ['A', 'K', 'M', 'T'], 60, ['E', 'Q', 'V']]
    Node11 = [58, ['I', 'T'], 78, ['K', 'R']]
    
    node1_val = decider(seq_in, Node1, 1, 1)
    if node1_val == "NA":
        return "NA"
    node2_val = decider(seq_in, Node2, 1, 2)
    if node2_val == "NA":
        return "NA"
    node3_val = decider(seq_in, Node3, 1, 3)
    if node3_val == "NA":
        return "NA"
    node4_val = decider(seq_in, Node4, 1, 4)
    if node4_val == "NA":
        return "NA"
    node5_val = decider(seq_in, Node5, 1, 5)
    if node5_val == "NA":
        return "NA"
    node6_val = decider(seq_in, Node6, 1, 6)
    if node6_val == "NA":
        return "NA"
    node7_val = decider(seq_in, Node7, 1, 7)
    if node7_val == "NA":
        return "NA"
    node8_val = decider(seq_in, Node8, 1, 8)
    if node8_val == "NA":
        return "NA"
    node9_val = decider(seq_in, Node9, 1, 9)
    if node9_val == "NA":
        return "NA"
    node10_val = decider(seq_in, Node10, 1, 10)
    if node10_val == "NA":
        return "NA"
    node11_val = decider(seq_in, Node11, 1, 11)
    if node11_val == "NA":
        return "NA"

    if node1_val:
        if node2_val:
            return "1A"
        else:
            return "1C"
    elif node3_val:
        return "1D"
    elif node4_val:
        if node5_val:
            return "1"
        elif node6_val:
            return "1B"
        else:
            return "1E"
    elif node7_val:
        if node8_val:
            return "1F"
        elif node9_val:
            return "1J"
        else: 
            return "1K"
    elif node10_val:
        return "1G"
    elif node11_val:
        return "1O"
    else: 
        return "1P"

# Classify sequences within cluster 2: 5 categories, 2A, 2B, etc.
def fork_2(seq_in):
    Node1 = [54, ['H', 'N', 'S'], 56, ['K', 'Q']]
    Node2 = [67, ['D', 'E', 'K', 'R'], 113, ['H', 'L', 'M', 'Q', 'R', 'Y']]
    Node3 = [115, ['A'], 51, ['G']]
    Node4 = [117, ['G', 'S'], 73, ['A', 'E', 'K', 'Q']]
    
    node1_val = decider(seq_in, Node1, 2, 1)
    if node1_val == "NA":
        return "NA"
    node2_val = decider(seq_in, Node2, 2, 2)
    if node2_val == "NA":
        return "NA"
    node3_val = decider(seq_in, Node3, 2, 3)
    if node3_val == "NA":
        return "NA"
    node4_val = decider(seq_in, Node4, 2, 4)
    if node4_val == "NA":
        return "NA"
    
    if node1_val:
        return "2A"
    elif node2_val:
        if node3_val:
            if node4_val:
                return "2"
            else:
                return "2D"
        else:
            return "2B"
    else:
        return "2C"

# Classify sequences within cluster 3: 18 categories, 3A, 3B, etc.
def fork_3(seq_in):
    Node1 = [50, ['T'], 77, ['Q']]
    Node2 = [85, ['A', 'C', 'E', 'F', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W'], 64, ['A', 'D', 'E', 'G', 'I', 'K', 'L', 'M', 'N', 'Q', 'T', 'V', 'Y']]
    Node3 = [84, ['I'], 49, ['G', 'I', 'V']]
    Node4 = [76, ['A', 'M', 'V'], 75, ['A', 'E', 'L', 'M', 'T', 'V']]
    Node5 = [57, ['L', 'S'], 112, ['A', 'H', 'P', 'S']]
    Node6 = [77, ['A', 'C', 'L', 'T'], 85, ['I', 'K', 'Q', 'R', 'S', 'T', 'V']]
    Node7 = [76, ['I', 'L', 'N', 'S', 'T'], 109, ['I', 'M', 'R', 'T', 'V']]
    Node8 = [61, ['A', 'M', 'N', 'P', 'T'], 82, ['A', 'C', 'D', 'G', 'H', 'K', 'L', 'N', 'Q', 'S', 'T']]
    Node9 = [117, ['K'], 78, ['I', 'S']]
    Node10 = [87, ['A', 'T'], 85, ['C', 'F', 'L', 'M', 'Q', 'W']]
    Node11 = [84, ['V'], 51, ['D', 'E', 'N']]
    Node12 = [76, ['E', 'V'], 117, ['A', 'H', 'S', 'T']]
    Node13 = [77, ['C', 'M', 'R', 'V'], 85, ['F', 'K', 'L', 'S', 'V', 'W']]
    Node14 = [85, ['A', 'F', 'H', 'I', 'K', 'L', 'N', 'P', 'Q', 'R', 'S', 'W'], 74, ['D', 'E', 'F', 'G', 'M', 'N', 'Q', 'S']]
    Node15 = [142, ['E', 'Q'], 85, ['A', 'F', 'H', 'I', 'K', 'L', 'N', 'P', 'Q', 'R', 'W']]
    Node16 = [51, ['H', 'N'], 84, ['C']]
    Node17 = [77, ['C', 'E', 'L', 'M', 'Q', 'R', 'S'], 85, ['A', 'F', 'H', 'I', 'K', 'L', 'N', 'P', 'Q', 'S', 'W']]
    Node18 = [79, ['A', 'C', 'E', 'F', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V'], 84, ['C', 'T', 'V']]

    node1_val = decider(seq_in, Node1, 3, 1)
    if node1_val == "NA":
        return "NA"
    node2_val = decider(seq_in, Node2, 3, 2)
    if node2_val == "NA":
        return "NA"
    node3_val = decider(seq_in, Node3, 3, 3)
    if node3_val == "NA":
        return "NA"
    node4_val = decider(seq_in, Node4, 3, 4)
    if node4_val == "NA":
        return "NA"
    node5_val = decider(seq_in, Node5, 3, 5)
    if node5_val == "NA":
        return "NA"
    node6_val = decider(seq_in, Node6, 3, 6)
    if node6_val == "NA":
        return "NA"
    node7_val = decider(seq_in, Node7, 3, 7)
    if node7_val == "NA":
        return "NA"
    node8_val = decider(seq_in, Node8, 3, 8)
    if node8_val == "NA":
        return "NA"
    node9_val = decider(seq_in, Node9, 3, 9)
    if node9_val == "NA":
        return "NA"
    node10_val = decider(seq_in, Node10, 3, 10)
    if node10_val == "NA":
        return "NA"
    node11_val = decider(seq_in, Node11, 3, 11)
    if node11_val == "NA":
        return "NA"
    node12_val = decider(seq_in, Node12, 3, 12)
    if node12_val == "NA":
        return "NA"
    node13_val = decider(seq_in, Node13, 3, 13)
    if node13_val == "NA":
        return "NA"
    node14_val = decider(seq_in, Node14, 3, 14)
    if node14_val == "NA":
        return "NA"
    node15_val = decider(seq_in, Node15, 3, 15)
    if node15_val == "NA":
        return "NA"
    node16_val = decider(seq_in, Node16, 3, 16)
    if node16_val == "NA":
        return "NA"
    node17_val = decider(seq_in, Node17, 3, 17)
    if node17_val == "NA":
        return "NA"
    node18_val = decider(seq_in, Node18, 3, 18)
    if node18_val == "NA":
        return "NA"

    if node1_val:
        return "3G"
    elif node2_val:
        if node3_val:
            if node4_val:
                return "3A"
            else:
                return "3C"
        elif node5_val:
            if node6_val:
                return "3K"
            elif node7_val:
                return "3E"
            else:
                return "3R"
        elif node8_val:
            if node9_val:
                return "3"
            elif node10_val:
                if node11_val:
                    return "3B"
                elif node12_val:
                    return "3H"
                elif node13_val:
                    return "3Q"
                else:
                    return "3P"
            elif node14_val:
                if node15_val:
                    if node16_val:
                        return "3J"
                    elif node17_val:
                        if node18_val:
                            return "3I"
                        else:
                            return "3Q"
                    else:
                        return "3N"
                else:
                    return "3L"
            else:
                return "3T"
        else:
            return "3S"
    else:
        return "3M"
    
# Classify sequences within cluster 4: 8 categories, 4A, 4B, etc.
def fork_4(seq_in):
    Node1 = [150, ['T'], 93, ['A', 'Q', 'S']]
    Node2 = [144, ['C', 'F', 'H', 'I', 'L', 'M', 'R', 'Y'], 55, ['A', 'I', 'L', 'M', 'N', 'P', 'R', 'S', 'T', 'V']]
    Node3 = [93, ['D', 'E', 'K', 'P', 'R'], 85, ['A', 'D', 'E', 'G', 'H', 'K', 'M', 'R', 'S', 'T', 'W', 'Y']]
    Node4 = [57, ['A', 'S'], 50, ['A', 'T']]
    Node5 = [147, ['A', 'G', 'K', 'M', 'S', 'T', 'V'], 61, ['E', 'H', 'Y']]
    Node6 = [154, ['C', 'E', 'G'], 93, ['E']]
    Node7 = [104, ['N', 'S'], 48, ['I', 'T', 'V']]
 
    node1_val = decider(seq_in, Node1, 4, 1)
    if node1_val == "NA":
        return "NA"
    node2_val = decider(seq_in, Node2, 4, 2)
    if node2_val == "NA":
        return "NA"
    node3_val = decider(seq_in, Node3, 4, 3)
    if node3_val == "NA":
        return "NA"
    node4_val = decider(seq_in, Node4, 4, 4)
    if node4_val == "NA":
        return "NA"
    node5_val = decider(seq_in, Node5, 4, 5)
    if node5_val == "NA":
        return "NA"
    node6_val = decider(seq_in, Node6, 4, 6)
    if node6_val == "NA":
        return "NA"
    node7_val = decider(seq_in, Node7, 4, 7)
    if node7_val == "NA":
        return "NA"
    
    if node1_val:
        return "4B"
    elif node2_val:
        if node3_val:
            if node4_val:
                return "4"
            elif node5_val:
                return "4A"
            elif node6_val:
                return "4D"
            else:
                return "4I"
        elif node7_val:
            return "4C"
        else:
            return "4F"
    else:
        return "4G"

# Main routine
assert (len(sys.argv) >= 2), "Must specify a multiple alignment file (A2M format)."
inname = str(sys.argv[1])
pull = " "        # Default is to not restrict output to a cluster or subcluster
if len(sys.argv) == 3:  pull = str(sys.argv[2])

from os import getcwd, path
outname = path.splitext(path.basename(inname))[0]           # Get file name, then drop extension if any
outname = path.join(getcwd(), outname + "_Clusters.fasta")  # Out file in current dir and renamed

#                                      get only first record for Azoto screening
infile = open(inname, "rt")
first_record = next(SeqIO.parse(infile, "fasta"))
infile.close()
refcensor = censor_maker(first_record.seq)
## jmagasin: Get the mapping after removing gaps.
residuePosToStringPos = GetAzoResidueToStringMap(ref_edit(first_record.seq, refcensor))
print("Good, the first sequence appears to be NifH from A.vinelandii.")
infile = open(inname, "rt")
outfile = open(outname, "w")
count = 0
for record in SeqIO.parse(infile, "fasta"):
    main_cluster = "ERROR"
    sub_cluster = "ERROR"
    workingseq = ref_edit(record.seq, refcensor)       
    main_cluster = primary_tree(workingseq)

    if main_cluster == "1":
        sub_cluster = fork_1(workingseq)
        
    if main_cluster == "2":
        sub_cluster = fork_2(workingseq)
            
    if main_cluster == "3":
        sub_cluster = fork_3(workingseq)

    if main_cluster == "4":
        sub_cluster = fork_4(workingseq)

    if count == 0 and not (main_cluster == "1" and sub_cluster == "1G"):
        # This must be the A.vinelandii NifH. It better be cluster 1, subcluster 1G.
        print("WARNING: The first sequence in the multiple alignment should be A.vinelandii\n",
              "NifH (e.g. WP_012698955.1) and thus in cluster 1, subcluster 1G.  It was not\n",
              "classified that way.  Probably you need to change the first sequence or make\n",
              "sure it has all the positions used by the CART model.\n")
        
    record.description = record.description + " main cluster = " + main_cluster + " subcluster = " + sub_cluster
    if (pull == main_cluster) or (pull == sub_cluster) or (pull == " "):
        SeqIO.write(record, outfile, "fasta")
        count = count + 1
print(count, "records on the output file")
infile.close()
outfile.close()
