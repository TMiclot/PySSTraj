#!/usr/bin/env python3
##################################################
##################################################
# NAME      PySSTraj
# AUTHOR    Tom MICLOT
# DATE      02/2003
# VERSION   1.0
# PYTHON    3.10.8
# LICENCE   MIT Licence
#
# DEPENDANCES   mdtraj, pandas
#
# DESCRIPTION
# DSSP code:
#    alpha_helix    H
#    beta_bridge    B
#    strand         E
#    helix_3        G
#    helix_5        I
#    turn           T
#    bend           S
#    loop          ' ' 
##################################################
##################################################





##################################################
##################################################
## Import Packages
import os
import argparse
from datetime import datetime
import mdtraj
import pandas
##################################################
##################################################





##################################################
##################################################
TIME_start = datetime.now()
##################################################
##################################################





##################################################
##################################################
## definition of functions
def sequence_compare(SEQUENCE_A, SEQUENCE_B):    # because analisis is performed on the same chain SEQUENCE_A and SEQUENCE_B have exactly the same lenght
    COUNT_similar = 0
    SEQUENCE_A_lenght = len(SEQUENCE_A)
    SEQUENCE_B_lenght = len(SEQUENCE_B)
    #
    ## for each position in sequence check if values are the same
    for POSITION in range (0,min(SEQUENCE_A_lenght,SEQUENCE_B_lenght)):
        if SEQUENCE_A[POSITION] == SEQUENCE_B[POSITION]:      # if values are the same
            COUNT_similar += 1    # is increased by 1
    #
    ## return similitude counted, proportion os similarity (0-1)
    return COUNT_similar, (COUNT_similar/SEQUENCE_A_lenght)

##################################################
##################################################





##################################################
##################################################
## Software informations (display in hel)
PARSER = argparse.ArgumentParser(
                    prog = 'PySSTraj.py',
                    description = 'What the program does',
                    epilog = "MIT License, Copyright (c) 2023 Tom Miclot")
#
#
## Parse command line arguments
PARSER.add_argument("-p", "--top", help="Topology file. If you want a result per chain, use a pdb file.")
PARSER.add_argument("-t", "--traj", help="Length of time series to fit the ARIMA model")
PARSER.add_argument("-s", "--str", default=1, type=int, help="Stride. Value must be an integer.")
PARSER.add_argument("-d", "--dir", default=".", help="Working directory.")
PARSER.add_argument("-v", "--ver", action="version", version='%(prog)s version 1.0 02/2023')
#
#
## Take arguments
ARGUMENTS = vars(PARSER.parse_args())
##################################################
##################################################





##################################################
##################################################
## Set up working directory
os.chdir(ARGUMENTS["dir"])
#
#
## set argument to variable
TOP_file = ARGUMENTS["top"]
TRAJ_file = ARGUMENTS["traj"]
STRIDE = ARGUMENTS["str"]
#
#
## Create majour output files
FILE_LOGGING = open('pysstraj.log', 'w')                                # Create logging file
FILE_all_chains_statistics = open('all_chains_statistics.txt', 'w')     # Create report file for all chains statistics
##################################################
##################################################






##################################################
##################################################
## Write head of the log file
FILE_LOGGING.write( '            >>>> Welcome to PySSTraj: Python DSSP on MD Trajectory <<<<' + '\n'*2)
FILE_LOGGING.write( '\n' + 'LOGGING FILE')
FILE_LOGGING.write( '\n' + '    Creation date    ' + str(datetime.now()) )     # writes file creation date/time to the log file
FILE_LOGGING.write( '\n'*3 )
FILE_LOGGING.write( '\n' + 'PARAMETERS')
FILE_LOGGING.write( '\n' + '    Topology            ' + TOP_file )             # writes topology file to the log file
FILE_LOGGING.write( '\n' + '    Trajectory          ' + TRAJ_file )            # writes trajectory file to the log file
FILE_LOGGING.write( '\n' + '    Stride              ' + str(STRIDE) )          # writes stride used information to the log file
FILE_LOGGING.write( '\n' + '    Wording directory   ' + ARGUMENTS["dir"] )     # writes work dir to the log file
##################################################
##################################################






##################################################
##################################################
## Write head of the all_chains_statistics file
FILE_all_chains_statistics.write('INFORMATIONS')
FILE_all_chains_statistics.write( '\n' + '    Count = Number of analyzed frames')
FILE_all_chains_statistics.write( '\n' + '    The values of the statistics are given in number of residues in the chain sequence.')
##################################################
##################################################





##################################################
##################################################
## Load trajectory and write system info in the log file
FILE_LOGGING.write( '\n'*4 + 'SYSTEM')
FILE_LOGGING.write( '\n' + '    INFO: Loading trajectory    ' + str(datetime.now()))    # writes advance to the log file
#
TRAJECTORY = mdtraj.load(TRAJ_file, top=TOP_file, stride=STRIDE)                          # load topology and trajectory
#
FILE_LOGGING.write( '\n' + '    INFO: Trajectory loaded     ' + str(datetime.now()))      # writes advance to the log file
FILE_LOGGING.write( '\n' + '    Chains      ' + str(TRAJECTORY.topology.n_chains) )       # writes system information to the log file
FILE_LOGGING.write( '\n' + '    Residues    ' + str(TRAJECTORY.topology.n_residues) )     # writes system information to the log file
FILE_LOGGING.write( '\n' + '    Atoms       ' + str(TRAJECTORY.topology.n_atoms) )        # writes system information to the log file
FILE_LOGGING.write( '\n' + '    Bonds       ' + str(TRAJECTORY.topology.n_bonds) )        # writes system information to the log file
FILE_LOGGING.write( '\n' + '    Frames      ' + str(TRAJECTORY.n_frames) )                # writes system information to the log file
##################################################
##################################################





##################################################
##################################################
# Set the variable to "chainid" 0 
INDEX=0
#
## Perform analyses for each chain of the system
while INDEX < TRAJECTORY.topology.n_chains:    # wile INDEX < number of chains
    ## Set up variables | they are reset for each chain
    DSSP_sequence = []
    FRAME_chain_data = []
    DATA_frame_residue = []
    DATA_matrix_convergence = []
    #
    #
    ## select protein in the correspond chainID
    SELECTION = "protein and chainid " + str(INDEX)                 # define selection as protein residue af the chain INDEX
    PROTEIN_selection = TRAJECTORY.topology.select(SELECTION)       # apply selection in the topology
    PROTEIN_trajectory = TRAJECTORY.atom_slice(PROTEIN_selection)   # apply selection in the traj
    #
    RESIDUES_SEQUENCE = TRAJECTORY.topology.to_fasta(INDEX)         # get the residue sequence in fasta format
    #
    #
    ## write curent work in the log file
    FILE_LOGGING.write( '\n'*4 + 'CHAIN_ID ' + str(INDEX) + '    (' + str(len(RESIDUES_SEQUENCE)) + ' residues)' )
    #
    #
    ## perform DSSP on the selection for each frame
    FILE_LOGGING.write( '\n' + '    INFO: Compute DSSP'  + ' '*38 + str(datetime.now()))    # writes advance to the log file
    DSSP = mdtraj.compute_dssp(PROTEIN_trajectory,simplified=False)                             # perfom DSSP 
    FILE_LOGGING.write( '\n' + '    INFO: End DSSP'  + ' '*42 + str(datetime.now()))    # writes advance to the log file
    #
    #
    ## convert DSSP (type np.ndarray of np.ndarray) into a list of list
    for ITEM in range(len(DSSP)):
        DSSP_sequence.append(DSSP[ITEM].tolist())
    FILE_LOGGING.write( '\n' + '    INFO: Convert DSSP datta'  + ' '*32 + str(datetime.now()))    # writes advance to the log file
    #
    #
    ## count number of each SS for each frame
    FILE_LOGGING.write( '\n' + '    INFO: Count SS/Frame'  + ' '*36 + str(datetime.now()))    # writes advance to the log file
    # 
    for FRAME in range(len(DSSP_sequence)):
        SS_SEQUENCE_chain = DSSP_sequence[FRAME]     # select the frame in DSSP_sequence list
        FRAME_chain_data.append([FRAME+1, SS_SEQUENCE_chain.count('H'), SS_SEQUENCE_chain.count('B'), SS_SEQUENCE_chain.count('E'), SS_SEQUENCE_chain.count('G'), SS_SEQUENCE_chain.count('I'), SS_SEQUENCE_chain.count('T'), SS_SEQUENCE_chain.count('S'), SS_SEQUENCE_chain.count(' ')]) # append data_frame_chain list with the new counted values
    #
    FILE_LOGGING.write( '\n' + '    INFO: End counting SS/Frame'  + ' '*29 + str(datetime.now()) )    # writes advance to the log file  
    FILE_LOGGING.write( '\n' + '    INFO: Analysis done with SS sequence length = ' + str(len(DSSP_sequence[0])) )
    if len(TRAJECTORY.topology.to_fasta(INDEX)) != len(DSSP_sequence[0]):
        FILE_LOGGING.write( '\n' + '    WARNING: SS sequence length (' + str(len(DDSSP_sequence[0])) + ' is diffrent from Chain sequence length (' + str(len(TRAJECTORY.topology.to_fasta(INDEX))) + ')'  + ' '*6 + str(datetime.now()) )
    #
    #
    ## Create database of SS/Frame value and perform statistics
    FILE_LOGGING.write( '\n' + '    INFO: Create SS/Frame database'  + ' '*26 + str(datetime.now()) )    # writes advance to the log file
    DATABASE_chain = pandas.DataFrame(FRAME_chain_data, columns=['Frame','Alpha_Helix','Beta_Bridge','Strand','Helix_3','Helix_5','Turn','Bend','None'])    # create the database with all SS count for each frame
    FILE_LOGGING.write( '\n' + '    INFO: SS/Frame database created'  + ' '*25 + str(datetime.now()) )    # writes advance to the log file
    #
    #
    ## write the DATABASE_chain of SS/frame in a file for the chainID
    FILE_LOGGING.write( '\n' + '    INFO: Export database into file' + ' '*25 + str(datetime.now()) )    # writes advance to the log file
    DATABASE_chain.to_csv('chain_' + str(INDEX) + '_SSframe.csv', index=None) # writes DATABASE_chain to the log file
    FILE_LOGGING.write( '\n' + '    INFO: Exportation done' + ' '*34 + str(datetime.now()) )    # writes advance to the log file
    #
    #
    ## write SS/frame statistics of the chain in the all_chains_statistics file
    FILE_LOGGING.write( '\n' + '    INFO: Export database statistics into file' + ' '*14 + str(datetime.now()) )    # writes advance to the log file
    FILE_all_chains_statistics.write('\n'*4 + '> CHAIN_ID ' + str(INDEX) + '    (' + str(len(RESIDUES_SEQUENCE)) + ' residues)' + '\n')
    FILE_all_chains_statistics.write(DATABASE_chain.drop('Frame', axis=1).describe().to_string() ) # remove the 'Frame' column and write performed stitistic, then convert the table to string
    FILE_LOGGING.write( '\n' + '    INFO: Exportation done' + ' '*34 + str(datetime.now()) )    # writes advance to the log file
    # 
    #
    ## Search SS structure of each residue in each frame
    FILE_LOGGING.write( '\n' + '    INFO: Search each residue SS/Frame'  + ' '*22 + str(datetime.now()) )    # writes advance to the log file
    for RESID in range(len(RESIDUES_SEQUENCE)):
        SS_SEQUENCE_residues = []    # setup the variable and reset it for each RESID
        for FRAME in range(len(DSSP_sequence)):
            SS_SEQUENCE_residues += DSSP_sequence[FRAME][RESID]
        # count number of each residues SS in the traj and divide by the number of frame to have a proportion
        DATA_frame_residue.append([RESID+1,RESIDUES_SEQUENCE[RESID],SS_SEQUENCE_residues.count('H')/len(DSSP_sequence), SS_SEQUENCE_residues.count('B')/len(DSSP_sequence), SS_SEQUENCE_residues.count('E')/len(DSSP_sequence), SS_SEQUENCE_residues.count('G')/len(DSSP_sequence), SS_SEQUENCE_residues.count('I')/len(DSSP_sequence), SS_SEQUENCE_residues.count('T')/len(DSSP_sequence), SS_SEQUENCE_residues.count('S')/len(DSSP_sequence), SS_SEQUENCE_residues.count(' ')/len(DSSP_sequence)])
    FILE_LOGGING.write( '\n' + '    INFO: Search done'  + ' '*39 + str(datetime.now()) )    # writes advance to the log file
    #
    #
    ## Create DATABASE_residues of SS count for each residue in all traj frames
    FILE_LOGGING.write( '\n' + '    INFO: Create residue SS database'  + ' '*24 + str(datetime.now()) )    # writes advance to the log file
    DATABASE_residues = pandas.DataFrame(DATA_frame_residue, columns=['ResID','Residue','Alpha_Helix','Beta_Bridge','Strand','Helix_3','Helix_5','Turn','Bend','None'])    # create the database with column names
    DATABASE_residues['SS_Count'] = DATABASE_residues.drop(['ResID','Residue'], axis=1).astype(bool).sum(axis=1)     # count number of SS type for each residue
    FILE_LOGGING.write( '\n' + '    INFO: Residue SS database created'  + ' '*23 + str(datetime.now()) )    # writes advance to the log file
    #
    #
    ## write the DATABASE_residues of residue SS in a file, for the chainID
    FILE_LOGGING.write( '\n' + '    INFO: Export residues database into file' + ' '*16 + str(datetime.now()) )    # writes advance to the log file
    DATABASE_residues.to_csv('chain_' + str(INDEX) + '_residues_SScount.csv', index=None)
    FILE_LOGGING.write( '\n' + '    INFO: Exportation done' + ' '*34 + str(datetime.now()) )    # writes advance to the log file
    #
    #
    ## Compute SS convergence matrix
    FILE_LOGGING.write( '\n' + '    INFO: Compute SS convergence matrix' + ' '*21 + str(datetime.now()) )    # writes advance to the log file
    for FRAME_A in range(len(DSSP_sequence)):
        for FRAME_B in range(len(DSSP_sequence)):
            COUNT, SIMILARITY = sequence_compare(DSSP_sequence[FRAME_A], DSSP_sequence[FRAME_B])
            DATA_matrix_convergence.append([FRAME_A, FRAME_B, COUNT, SIMILARITY])
    #
    #
    ## Create DATABASE_matrix of SS sequence convergence
    FILE_LOGGING.write( '\n' + '    INFO: Create SS convergence matrix database'  + ' '*13 + str(datetime.now()) )    # writes advance to the log file
    DATABASE_matrix = pandas.DataFrame(DATA_matrix_convergence, columns=['Frame_A','Frame_B','Count','Similarity'])
    FILE_LOGGING.write( '\n' + '    INFO: SS convergence matrix database created'  + ' '*12 + str(datetime.now()) )    # writes advance to the log file
    #
    #
    ## write the DATABASE_matrix in a file, for the chainID
    FILE_LOGGING.write( '\n' + '    INFO: Export SS convergence matrix into file' + ' '*12 + str(datetime.now()) )    # writes advance to the log file
    DATABASE_matrix.to_csv('chain_' + str(INDEX) + '_SSconv_matrix.csv', index=None)
    FILE_LOGGING.write( '\n' + '    INFO: Exportation done' + ' '*34 + str(datetime.now()) )    # writes advance to the log file
    # 
    #
    ## END of work for the chain
    FILE_LOGGING.write( '\n' + '    *END: Analysis was completed normally' + ' '*19 + str(datetime.now()) )
    INDEX += 1    # increate INDEX by 1 to explore the next chain ID
##################################################
##################################################







##################################################
##################################################
## End text
TIME_end = datetime.now()
FILE_LOGGING.write( '\n'*4 + 'END OF FILE')
FILE_LOGGING.write( '\n' + "    PySSTraj.py was completed normally !" )
FILE_LOGGING.write( '\n' + "    Execution finished    " + str(TIME_end) )
FILE_LOGGING.write( '\n' + "    Execution time        " + str(TIME_end - TIME_start) )
#
#
## Close files
FILE_all_chains_statistics.close()     # close all all_chains_statistics file
FILE_LOGGING.close()                   # close logging file
##################################################
##################################################
