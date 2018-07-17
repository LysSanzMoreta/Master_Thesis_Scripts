#!/usr/bin/env python
#####IMPORT ALL NECESSARY MODULES#####



#### ONLY CREATES THE EQUIVALENT POSITIONS DATAFRAME. NO PyMOL, no GUI

#EXAMPLE TO RUN
#python Link_Your_Sites_python_3_Updated.py --PDB "C:/Users/Lys Sanz Moreta/Documents/New_genes_10_23_2017/PINK1.pdb" --Gene "C:/Users/Lys Sanz Moreta/Dropbox/THESIS/GENEaln/New_IDS/TGUT00000017792.map2finch.fasta" --M8 "C:/Users/Lys Sanz Moreta/Documents/New_genes_10_23_2017/PINK1_TGUT00000017792.map2finch.fasta/M8/out" --Full_PDB_sequence "C:/Users/Lys Sanz Moreta/Documents/New_genes_10_23_2017/PINK1.fasta" --domains "list(range(1,229))" --format fasta --sequence_number 1
#python Link_Your_Sites_python_3_Updated.py --PDB "C:/Users/Lys Sanz Moreta/Documents/New_genes_10_23_2017/PINK1.pdb" --Gene "C:/Users/Lys Sanz Moreta/Dropbox/THESIS/GENEaln/New_IDS/TGUT00000017792.map2finch.fasta" --M8 "C:/Users/Lys Sanz Moreta/Documents/New_genes_10_23_2017/PINK1_TGUT00000017792.map2finch.fasta/M8/out" --Full_PDB_sequence "C:/Users/Lys Sanz Moreta/Documents/New_genes_10_23_2017/PINK1.fasta" --domains "list(range(1,229))" --sequence_number 1

import sys, ast, json
import os
import argparse
import numpy as np
import pandas as pd
import itertools
import ntpath
import numbers
import decimal

#Biopython
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from Bio.Alphabet import generic_protein
from Bio import SeqRecord,Alphabet,SeqIO
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import Bio.PDB as PDB
from Bio.Seq import MutableSeq
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.pairwise2 import format_alignment



#Command Line Arguments
parser = argparse.ArgumentParser()

#COMPULSORY files: Give error if not present
parser.add_argument("--PDB",required = True,help ='Path to Protein DataBank file (use "" if errors)')
parser.add_argument('--Gene',required = True,help = 'Path to Single nucleotide gene sequence or multiple alignment nucleotide sequences')
parser.add_argument('--M8',required =True,help ='Path to Out file from the Codeml model M8')
parser.add_argument('--Full_PDB_sequence',required =True,help = 'Path to the Complete Fasta file of the PDB sequence (download from PDB as well)')
#OPTIONAL arguments:
parser.add_argument('--domains',help = 'List of funcional domain positions retrieved from the PDB protein webpage starting in index 1 eg "[1,4,6]" or "list(range(1,30))+list(range(70,80))" (quotes are essential!), default=[]', default=[])
parser.add_argument('--file_list_domains',help = 'Path to file with functional domain positions in format 1\n1\n , default = []', default=[] )
parser.add_argument('--file_domains',help = 'Path to file with protein domains sequences (all in same format as alignment), default = []', default=[] )
parser.add_argument('--format', help = 'Sequence or Multiple Alignment File Format (default fasta), do not write it with quotes',default = 'fasta')
parser.add_argument("--prob", help = 'Choice of level of posterior probability on the sites, 95% or 99% from M8 Out file', default= int(99)) 	# prob : Posterior probability percentage of the positive sites , 99 or 95, default= 99
parser.add_argument("--missing_data", help = 'Decide if the missing data("N") should be removed from the nucleotide sequence, default = yes, recommended for the alignment', default='yes') 	# missing_data : remove or not missing data from the sequence alignment: yes or no, default=yes
parser.add_argument("--sequence_number", help = 'Number of the sequence in the multiple alignment file', default=int(0)) 	# sequence_number : When using a multiple alignment file indicate with sequence number to use, if single sequence use 0, default = 0
parser.add_argument("--print_alignment",help= 'Choose to visualize the PDB file sequence aligned with the gene, default = no', default='no')
parser.add_argument("--PDB_chains",help= 'Chain(s) to extract from the PDB sequence(s), write separated by spaces: A B C, default = all ',nargs = '*',default='all')
args = parser.parse_args()
PDB_file,Full_PDB_sequence, Gene, M8,domains,List_domains,File_domains,Gene_file_format,prob,missing_data,sequence_number,print_alignment,chains = [args.PDB, args.Full_PDB_sequence,args.Gene, args.M8,args.domains,args.file_list_domains,args.file_domains,args.format,args.prob,args.missing_data,args.sequence_number,args.print_alignment,args.PDB_chains]
try:
    sequence_number = eval(sequence_number)
except:
    pass
try: #LINUX USERS
    prob =eval(prob)
except:
    pass
#Check all the necessary files are present and the right arguments are given:
File_formats = ['fasta','phylip_sequential','clustal','embl','genebank','gb','abi','ace','fastq-sanger','fastq','fastq-solexa','fastq-illumina','ig','imgt','pdb-seqres','pdb-atom','phd','phylip','pir','seqxml','sff','stockholm','swiss','tab','qual','uniprot-xml']
if not PDB_file or not Gene or not M8 or not Full_PDB_sequence:
    parser.error('Missing one of the 4 required arguments: PDB file / Gene / Out file / Full PDB Sequence')
elif prob not in [95,99]:
    parser.error('Probability values can only be 95  or 99')
elif missing_data not in ['no','yes']:
    parser.error('Choose to keep missing data: yes or to remove it: no')
elif not isinstance(sequence_number, numbers.Number):
    parser.error('The number of the sequence in the alignment needs to be an integer')
elif Gene_file_format not in File_formats:
    parser.error('Invalid Gene File Format.Check available SeqIO biopython file formats ')
elif List_domains and File_domains:
    parser.error('Duplicated information, choose either --domains or --file_domains')
else:
    if not List_domains and not File_domains:
        print('Not List of Domains or Path to Text file with domains specified, using empty list')
        print('Building dataframe...')
    else:
        print('Building dataframe...')



def Global_alignment(chain_A,chain_B):
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
    alignment_global = pairwise2.align.globalms(chain_A, chain_B, 2, -1, -5, -1)
    #print(format_alignment(*alignment_global[0]))
    alignment_info_global = alignment_global[0]
    Aligned_A, Aligned_B, score, begin, end = alignment_info_global
    return alignment_info_global,alignment_global

def Local_alignment(chain_A,chain_B):
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
    alignment_local = pairwise2.align.localms(chain_A, chain_B, 2, -1, -5, -1)
    #print(format_alignment(*alignment_global[0]))
    alignment_info_local = alignment_local[0]
    Aligned_A, Aligned_B, score, begin, end = alignment_info_local
    return alignment_info_local,alignment_local

def fasta_to_sequence(Fasta_file,Format):
    from Bio import SeqRecord, SeqIO
    fasta_sequences = SeqIO.parse(open(Fasta_file), Format)
    for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            return sequence,name

def fasta_to_sequences(Fasta_file,Format):
    from Bio import SeqRecord, SeqIO
    List_of_sequences = []
    List_of_names = []
    fasta_sequences = SeqIO.parse(open(Fasta_file),Format)
    for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            List_of_sequences.append(sequence)
            List_of_names.append(name)
    return List_of_sequences,List_of_names

def Translate_sequence(Fasta_file,Format,sequence_number):
    if sequence_number != 0:
        sequence = fasta_to_sequences(Fasta_file,Format)[0][sequence_number]
        aa_seq = Seq(sequence).translate(stop_symbol="X")

    else:
        sequence= fasta_to_sequence(Fasta_file,Format)[0]
        aa_seq = Seq(sequence).translate(stop_symbol="X")
    return aa_seq

def Translate_and_Remove_missing_data(Fasta_file,Format,sequence_number):
    '''Remove the missing data ('X') from the sequences after being translated, otherwise the codons are affected'''
    if sequence_number != 0:
        sequence =fasta_to_sequences(Fasta_file,Format)[0][sequence_number]
        clean_sequence = Seq(sequence).translate(stop_symbol="X")
        clean_sequence = clean_sequence.ungap('X')
    else:
        sequence=fasta_to_sequence(Fasta_file,Format)[0]
        clean_sequence = Seq(sequence).translate(stop_symbol="X")
        clean_sequence = clean_sequence.ungap('X')
    return clean_sequence

def Extract_sequence_from_PDB(PDB_file,chains):
    ''' Returns both the sequence contained in the PDB file and the residues coordinates'''
    parser = PDB.PDBParser()
    Name = ntpath.basename(PDB_file).split('.')[0]
    #Name = PDB_file.split('/')[-1].split('.')[0]
    structure = parser.get_structure('%s' %(Name),PDB_file)
    ############## Iterating over residues to extract all of them even if there is more than 1 chain
    sequence = []
    Residues_ID = []
    if chains == 'all':
        for chain in structure.get_chains():
            for residue in chain:
                if is_aa(residue.get_resname(), standard=True):
                    sequence.append(residue.get_resname())
                    Residues_ID.append(residue.get_id()[1])
    else :
        accumulated_residues = []
        accumulated_ids = []
        for letter in chains:
           try: #in case the chain requested does not exits
                for residue in structure[0][letter]:
                    #print(letter)
                    if is_aa(residue.get_resname(), standard=True):
                        accumulated_residues.append(residue.get_resname())
                        accumulated_ids.append(residue.get_id()[1])
           except:
               pass
        sequence.append(''.join(accumulated_residues))
        Residues_ID.append(accumulated_ids)

    joined_sequence = ''.join(sequence)
    PDB_sequence = seq1(joined_sequence) #3 letter code into 1 letter code
    try:
        Residues_ID = list(itertools.chain.from_iterable(Residues_ID))
    except:
        pass
    #print(PDB_sequence)
    #print(Residues_ID)

    return Residues_ID,PDB_sequence


def equivalent_positions(chain_A, chain_B,Aligned_A, Aligned_B,Residues_ID = None, Domain_Positions =None): #chainA = PDB #chainB=Transcript #Positions = List of Positive Selected sites (optional argument because we will use this function in 2 manners)
    ''' This function returns the corresponding coordinates of the gene residues in the PDB sequence, the domain positions in the trimmed PDB sequence
    or among the gene with missing data and without'''
    import numpy as np
    ##Equivalent positions
    List = [] #Will store the Positions where there are no gaps in the aligned sequences, these are still not the PDB positions !!
    for index,residue in enumerate(Aligned_A):
        if Aligned_A[index] != '-' and Aligned_B[index] != '-':
                #OriginalA[:index] += 1
                List.append(index +1) #In index 1
        else:
            pass

    Final_Positions = ['nan']*len(chain_B) #should have the lenght of the transcript, will replace some nan with the equivalent PDB positions
    #####Finding the equivalent PDB Position of the residues in the gene ignoring the GAPS generated in the alignment#######
    Position_in_PDB =[]
    gaps_first_segment = ''.join(Aligned_A[0:List[0]]).count('-')
    Position_in_PDB.append(List[0] - gaps_first_segment)
    Position_in_PDB.append(List[1] - gaps_first_segment - ''.join(Aligned_A[List[0]:List[1]]).count('-'))

    accumulated_number_gaps_in_this_segment = gaps_first_segment
    for i in range(0,len(List)): #we skip the first interval, already performed #I changed 1 back to 0 again,i think I did it when residues ID was giving problems
        try:
            accumulated_number_gaps_in_this_segment += ''.join(Aligned_A[List[i]:List[i+1]]).count('-')
            Position_in_PDB.append(List[i+1] - accumulated_number_gaps_in_this_segment)
        except:
            pass
    #####Finding out the equivalent positions in the transcript of the PDB residues######
    Position_in_Transcript = []
    gaps_first_segment_Transcript = ''.join(Aligned_B[0:List[0]]).count('-')
    Position_in_Transcript.append(List[0] - gaps_first_segment_Transcript)
    Position_in_Transcript.append(List[1] - gaps_first_segment_Transcript - ''.join(Aligned_B[List[0]:List[1]]).count('-'))

    accumulated_number_gaps_in_this_segment_transcript = gaps_first_segment_Transcript
    for i in range(0, len(List)):  # we skip the first interval
        try:
            accumulated_number_gaps_in_this_segment_transcript += ''.join(Aligned_B[List[i]:List[i + 1]]).count('-')
            Position_in_Transcript.append(List[i+1] - accumulated_number_gaps_in_this_segment_transcript) # plus on otherwise negative numbers
        except:
            pass
    Equivalent_Domain_positions = []
    if not Domain_Positions:
        for position_transcript,position_PDB in zip(Position_in_Transcript,Position_in_PDB):
                 if Residues_ID:

                    Final_Positions[position_transcript - 1] = Residues_ID[position_PDB - 1]

                 else:
                    Final_Positions[position_transcript-1]  = position_PDB

        return Final_Positions
    else:
        for position_transcript, position_PDB in itertools.zip_longest(Position_in_Transcript, Position_in_PDB):
            Final_Positions[position_transcript - 1] = position_PDB

        for residue_id, domain_position in itertools.zip_longest( Residues_ID,Domain_Positions):
            try:
                specific_corresponding_domain_position = Final_Positions[domain_position - 1]
                try:
                    Equivalent_Domain_positions.append(Residues_ID[specific_corresponding_domain_position - 1])
                except:
                    pass
            except:
                pass

        return Equivalent_Domain_positions


def List_of_positions_of_Positive_Sites(file,prob):  # file = Out file from M8 #prob = 95 or 99
        '''Will return the Positions of the selected sites with regards to the ALIGNMENT, later use other fucntion to find the equivalent in the clean from missing data sequence'''
        length_alignment = []
        with open(file, 'r') as f:
            data = f.read().split('\n')
            positions = []
            line_start = []
            line_stop = []  # Line number reference to stop collecting the info related to positive selected sites
            for number, line in enumerate(data, 1):
                if 'Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)' in line:
                    line_start = number
                if 'The grid' in line:
                    line_stop = number
                #if number == 1:
                    #length_alignment = line.split()[1]
            if 'Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)' in data:
                diff = line_stop - line_start - 3
                for i in range(6, int(diff)):  # Start at line 6 after match
                    position = data[data.index(
                        'Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)') + i]
                    if prob == 99 :
                        if str(position.split()[2]).endswith('**'):  # * > 95% confidence ** > 99 confidence
                            # print(position)
                            position = position.split()[0]
                            positions.append(int(position))  # Needs to be an int for later comparison
                    else:
                        if str(position.split()[2]).endswith('*'):  # * > 95% confidence ** > 99 confidence
                            # print(position)
                            position = position.split()[0]
                            positions.append(int(position))  # Needs to be an int for later comparison
            return positions

def Corresponding_positions_missing_notmissing_data(Missing_data,Clean): #Missing_data = WITH missing data ; Clean = NO missing data
    ''' Returns list of the equivalent positions among 2 sequences, in this case the same sequence with and without missing data'''

    alignment_info_global, alignment_global = Global_alignment(Missing_data,Clean)
    Aligned_A, Aligned_B, score, begin, end = alignment_info_global #A = Gapped; B = Ungapped
    #print(format_alignment(*alignment_global[0]))
    List = equivalent_positions(Missing_data,Clean, Aligned_A, Aligned_B) #All corresponding positions of the gapped to the ungapped
    return List

def Corresponding_functional_positions(PDB,Full_PDB_sequence,Residues_ID,functional_list):
    alignment_info_global, alignment_global = Global_alignment(PDB,Full_PDB_sequence) #CHECK THIS is in the right order and it should be the other way round
    Aligned_A, Aligned_B, score, begin, end = alignment_info_global
    # print(format_alignment(*alignment_global[0]))
    List = equivalent_positions(PDB,Full_PDB_sequence,Aligned_A,Aligned_B,Residues_ID=Residues_ID,Domain_Positions=functional_list)  # All corresponding positions of the gapped to the ungapped
    return List


def Corresponding_Coordinates_and_labels_PDB_Gene(PDB,Gene,Full_PDB_sequence,List_of_Positive_Positions,functional_list,Residues_ID,basepath,print_alignment,Clean_positions = None):

    '''Performs the alignment, retrieves the equivalent coordinates and labels among the PDB sequence and the gene for each of the PDB positions
    as 'Not', which refers to not being positively selected or a functional domain, 'Selected', which stands for positive selected, 'Domain',
    belongs to the functional domain or 'Selected_and_Domain' which stands for both'''

    alignment_info_global, alignment_global = Global_alignment(PDB, Gene)
    Aligned_A, Aligned_B, score, begin, end = alignment_info_global #A = PDB; B = Gene(with or without missing data)
    if print_alignment == 'yes':
        print(format_alignment(*alignment_global[0]))
    else:
        pass
    List_positions = list(range(1, len(Gene) + 1))
    #Extract the corresponding positions of the positive selected sites in the clean of missing data gene sequence
    List_positive_positions = []
    if Clean_positions: #If the gene sequence has been cleaned and we have the corresponding positions
        for element in List_of_Positive_Positions:
            try:
                List_positive_positions.append(Clean_positions.index(element)) #Append the index of the residue that has the same id as the one in positive selected sites
            except:
                pass
    else: #we don't do anything
        List_positive_positions = List_of_Positive_Positions
    #print(List_of_Positive_Positions,List_positive_positions)
    #For the dataframe we can label the positions from 1 to length of gene sequence (with or without missing data)
    positions_dataframe = pd.DataFrame(pd.Series(List_positions))
    positions_dataframe.rename(columns={positions_dataframe.columns[0]: "Gene_Position"}, inplace=True)

    List = equivalent_positions(PDB, Gene, Aligned_A, Aligned_B,Residues_ID = Residues_ID) #List of equivalent positions of each of the gene residues into the PDB
    ###The functional domains require also processing to find their equivalent residues in the PDB
    alignment_info_global, alignment_global = Global_alignment(PDB, Full_PDB_sequence)
    Aligned_A, Aligned_B, score, begin, end = alignment_info_global  # A = PDB; B = Gene(with or without missing data)
    #functional_list = equivalent_positions(PDB,Full_PDB_sequence,Aligned_A,Aligned_B,Residues_ID=Residues_ID,Domain_Positions=functional_list)
    functional_list = equivalent_positions(PDB,Full_PDB_sequence, Aligned_A, Aligned_B, Residues_ID=Residues_ID,Domain_Positions=functional_list)
    #print(functional_list)
    positions_dataframe['PDB_Position'] = List
    Label_1 = positions_dataframe['Gene_Position'].isin(List_positive_positions) #Check if the position in the gene sequence is + selected
    Label_2 = positions_dataframe['PDB_Position'].isin(functional_list) #Check if the PDB position is a functional domain

    positions_dataframe['Label_1'] = Label_1  # Create a Column were the positive positions have a label True
    positions_dataframe['Label_2'] = Label_2 #Create a column were the functional domain positions recieve the label True
    positions_dataframe['Label_1'] = positions_dataframe['Label_1'].replace(True, 'Selected')
    positions_dataframe['Label_1'] = positions_dataframe['Label_1'].replace(False, 'Not')
    positions_dataframe['Label_2'] = positions_dataframe['Label_2'].replace(True, 'Domain')
    positions_dataframe['Label_2'] = positions_dataframe['Label_2'].replace(False, 'Not')
    positions_dataframe['Label'] = pd.Series(List_positions)
    for index, row in positions_dataframe.iterrows():

        if positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label_1')] == 'Selected' and positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label_2')] == 'Domain':
            positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label')] = 'Selected_and_Domain'
        elif positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label_1')] == 'Selected':
            positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label')] = 'Selected'
        elif positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label_2')] == 'Domain':
            positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label')] = 'Domain'
        else:
            positions_dataframe.iloc[index, positions_dataframe.columns.get_loc('Label')] = 'Not'
    positions_dataframe.drop(['Label_1', 'Label_2'], 1, inplace=True)
    # print(positions_dataframe.columns)
    # positions_dataframe.drop(positions_dataframe.columns[0],1,inplace=True)
    Directory = os.path.dirname(basepath)  # Obtain absolute path from the gene file
    #Gene_name = ntpath.basename(Gene)
    PDB_name = ntpath.basename(basepath)
    positions_dataframe.to_csv(Directory + "/%s_Positions" % (PDB_name), index=False, sep='\t')
    print('Dataframe Ready at %s!' % (Directory))


def Read_List_of_domains(List_domains,domains,File_domains):
    """Handles 3 possible types of input of the residues that conform the protein domain and returns their positions with regards to the Full sequence of the protein provided in the PDB website"""
    if not File_domains:
        if not domains:
            if not List_domains:
                List_of_domains = []

            else:
                List_of_domains = [line.rstrip('\n') for line in open(List_domains)]
                List_of_domains = list(map(int, List_of_domains))
        else:
            # Convert the input list in string format to a string by 'executing' it
            List_of_domains = eval(domains)
        #List_domains = Read_Domains(List_domains,domains)  # domains is the manually introduced numbers, list domains is the fle with the number line by line
    else:#If there is a file with the domain sequences we get the equivalent positions with regards to the complete PDB sequence, not the one in the PDB file yet
        Full_PDB = ''.join(fasta_to_sequences(Full_PDB_sequence, 'fasta')[0])
        Sequences_domains = ''.join(fasta_to_sequences(File_domains,Gene_file_format)[0])
        #alignment_info_global, alignment_global = Global_alignment(Full_PDB,Sequences_domains)
        #Aligned_Full, Aligned_Domain, score, begin, end = alignment_info_global  # A = PDB; B = Gene(with or without missing data)
        ########################LOCAL ALIGNMENT FOR THE DOMAIN SEQUENCES , THEY ARE TOO SHORT FOR A GLOBAL
        alignment_info_local, alignment_local = Local_alignment(Full_PDB,Sequences_domains)
        Aligned_Full, Aligned_Domain, score, begin, end = alignment_info_local  # A = PDB; B = Gene(with or without missing data)
        List_of_domains = equivalent_positions(Full_PDB_sequence,Sequences_domains,Aligned_Full,Aligned_Domain)
    return List_of_domains
List_domains = Read_List_of_domains(List_domains,domains,File_domains) #Choose upon the 3 options to get the numbers

def Wrapper_of_all_functions(PDB_file,Gene,Full_PDB_sequence,M8,List_Domains,Format,prob,Sequence_number,missing_data,print_alignment,chains):
    '''Calling all the functions in the corresponding order/combination according to the optional arguments'''
    Residues_ID,PDB_sequence = Extract_sequence_from_PDB(PDB_file,chains)
    Full_PDB_sequence = ''.join(fasta_to_sequences(Full_PDB_sequence,'fasta')[0]) #We use all the chains in the file (fasta file provided by PDB)
    List_of_Positive_Positions =List_of_positions_of_Positive_Sites(M8,prob)

    #Extract sequence no matter if single sequence or multiple alignment

    Gene_missing_data = fasta_to_sequences(Gene,Format)[0][Sequence_number]
    #Checking if the user wants to perform the alignment with or without missing data in the gene
    if missing_data == 'no':
        Clean_protein_sequence = Translate_and_Remove_missing_data(Gene,Format,Sequence_number)
        Protein_missing_data = Translate_sequence(Gene, Format, Sequence_number)  # Translate
        Clean_positions = Corresponding_positions_missing_notmissing_data(Protein_missing_data, Clean_protein_sequence)
        Corresponding_Coordinates_and_labels_PDB_Gene(PDB_sequence,Clean_protein_sequence,Full_PDB_sequence, List_of_Positive_Positions,List_domains,Residues_ID,PDB_file,print_alignment, Clean_positions)

    else: #Gene_missing_data is our sequence
        Protein_missing_data = Translate_sequence(Gene,Format,Sequence_number) #Translate
        Corresponding_Coordinates_and_labels_PDB_Gene(PDB_sequence,Protein_missing_data,Full_PDB_sequence,List_of_Positive_Positions,List_domains,Residues_ID,PDB_file,print_alignment)


if __name__ == "__main__":
    Wrapper_of_all_functions(PDB_file,Gene,Full_PDB_sequence,M8,List_domains,Gene_file_format,prob,sequence_number,missing_data,print_alignment,chains)
