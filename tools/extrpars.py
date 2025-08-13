#!/usr/bin/env python3
import re

"""
 For each file, extract the "jobs", delineated at the start by either
   "RESTRICTED RIGHTS LEGEND:"  or   "Link1:  Proceeding to internal job step number"
 and at the end by
   " termination"
"""
def extract_jobs(input_file):
  jobs = []  # List to store extracted jobs
  current_job = []  # Temporary storage for the current job
  inside_job = False  # Flag to track if we are inside a job

  with open(input_file, 'r') as file:
    for line in file:
      if "RESTRICTED RIGHTS LEGEND" in line:  # Start of a new job
        inside_job = True
        current_job = []  # Reset the current job
      elif "Link1:  Proceeding to internal job step number" in line:   # Start of linked job; treat as new job
        inside_job = True
        current_job = []  # Reset the current job
      elif " termination" in line:  # End of the current job
        if inside_job:
          jobs.append(current_job)  # Save the completed job
          inside_job = False
      elif inside_job:
        current_job.append(line.strip())  # Add lines to the current job
  return jobs

def parse_job_method(job):
  method = ""
  jobtype = "sp"  # default jobtype in Gaussian is single-point energy
# First pass through job to get jobtype and number of atoms (determine geometry format),
#  method (determines energy format)
  num_atoms = 0
  flag_cp_ver = False
  for line in job:
# try to find method; once this is set, stop checking
    if method == "":
      if re.search(r'^#.*CCSD', line, flags=re.IGNORECASE):  # this must appear before CCSDT so CCSDT can overwrite
        method = "CCSD"
      if re.search(r'^#.*CCSD\(T\)', line, flags=re.IGNORECASE):
        method = "CCSDT"
      elif re.search(r'^#.*QCISD', line, flags=re.IGNORECASE):
        method = "QCISD"
      elif re.search(r'^#.*MP2', line, flags=re.IGNORECASE):
        method = "MP2"
      elif re.search(r'^#.*B3LYP', line, flags=re.IGNORECASE):
        method = "B3LYP"
      elif re.search(r'^#.*MN15', line, flags=re.IGNORECASE):
        method = "MN15"
      elif re.search(r'^#.*HF', line, flags=re.IGNORECASE):
        method = "HF"

# try to find jobtype; once this is set, stop checking
    if jobtype == "sp":
      if re.search(r'#.*opt*', line, flags=re.IGNORECASE):
        jobtype = "opt"  # first in list so it can be overwritten by more specific opt jobtypes
#        if re.search(r'#.*opt.*freq*', line, flags=re.IGNORECASE):
#          jobtype = "opt_freq"
        if re.search(r'#.*opt.*modredundant', line, flags=re.IGNORECASE):
          jobtype = "opt_mod"
        elif re.search(r'#.*opt.*z-mat*', line, flags=re.IGNORECASE):
          jobtype = "opt_zmat"
      elif re.search(r'Scan the potential surface.', line, flags=re.IGNORECASE):
        jobtype = "scan_zmat"
      elif re.search(r'#.*freq*', line, flags=re.IGNORECASE) and not re.search(r'#.*opt*', line, flags=re.IGNORECASE):
        jobtype = "freq"
# check if counterpoise job
    if re.search(r'#.*counterpoise*', line, flags=re.IGNORECASE):
      flag_cp_ver = True
# get number of atoms
    if re.search("NAtoms.*NActive.*NUniq", line, flags=0):
      num_atoms = line.split()[1]
  return method, jobtype, flag_cp_ver, num_atoms


"""
 For each job in the current file, extract the "blocks" as delneated at the start by 
   the first line of the job   or   the first line after the last line of the previous block
 and at the end by
   "Optimized Parameters" followed by two instances of "---------------------- [] -----"
     for jobtypes opt, opt=z-matrix (with and without scan), freq
   " Variable Step   Value"  followed by two instances of "---------------------- [] -----"
       or
   " Scan completed."         for jobtype scan 
   "GINC"                     for all others
"""
def extract_blocks(job,jobtype):
  blocks = []  # List to store extracted blocks
  current_block = []  # Temporary storage for the current block
  inside_block = False  # Flag to track if we are inside a block
#  print("lines in job = ",len(job))
#  print("jobtype = ",jobtype)
  for idx, line in enumerate(job):
    if inside_block == False:  # Start of a new block
      inside_block = True
      sig_stop = 0
#      print("sigstop = ",sig_stop," line in job = ",idx)
      current_block = [line.strip()]  # Reset the current block
    if "opt" in jobtype or jobtype == "scan_zmat":
      if sig_stop == 2:
        if "------------------------------------------------------------------------" in line:
#          print("lines in block = ",len(current_block))
          blocks.append(current_block)  # Save the completed block
          inside_block = False
#          print("sigstop = ",sig_stop," line in job = ",idx)
      elif sig_stop == 1:
        if "------------------------------------------------------------------------" in line:
#          print("sigstop = ",sig_stop," line in job = ",idx)
          sig_stop = 2
      elif sig_stop == 0:
        if "opt" in jobtype:
          if "Optimized Parameters" in line:  # Cue for approaching end of the current block
#            print("sigstop = ",sig_stop," line in job = ",idx)
            sig_stop = 1
        elif jobtype == "scan_zmat":
          if "Variable Step   Value" in line:  # Cue for approaching end of the current block
#            print("sigstop = ",sig_stop," line in job = ",idx)
            sig_stop = 1
          elif  "Scan completed." in line:
#            print("lines in block = ",len(current_block))
            blocks.append(current_block)  # Save the completed block
            inside_block = False
    elif "GINC" in line:    # block termination for all other job types; block should not read archive text in logfile (which can mislead keystring search) so terminate block no lower in file than this
      current_block.append(line.strip())  # Add line to the current block
#      print("lines in block = ",len(current_block))
      blocks.append(current_block)  # Save the completed block
#      inside_block = False
    if inside_block:
      current_block.append(line.strip())  # Add line to the current block
  return blocks


def parse_block(block):
  flag_stat = 0
  for line in block:
    if "Stationary" in line:
      flag_stat = 1
  return flag_stat


def extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines):
  i_lines = 0
  num_lines = int(num_skiplines) + int(num_readlines)  # number of lines to read into read_lines
  lines = []  # clear lines
  read_lines = []  # clear read_lines
  current_lines = []  # clear current_lines
  flag_anch1 = False   # flag_anch1 is set if anchor1 is found, should be False after loop if anchor1=noanchor
  flag_read = False   # flag_read is set when anchor2 is encountered and unset num_lines later
  flag_break = False   # flag_break ends the loop once the string is read for a "getany" case
#  print("anchor1 = ",anchor1)
#  print("anchor2 = ",anchor2)
  for line in block:
    if anchor1 != "getall":
      if re.search(anchor1, line, flags=0):
#        print("D")
        flag_anch1 = True
        current_lines = read_lines   # store most recent read_lines in holding list
      if re.search(anchor2, line, flags=0):
#        print("A")
        flag_read = True 
        read_lines = []  # clear read_lines
        i_lines = 0
      if flag_read:
#        print("B")
        read_lines.append(line)   # save line in read_lines
        i_lines += 1
        if i_lines == num_lines:
#          print("C")
          flag_read = False
          if anchor1 == "getany":   # read the first lines triggered by anchor2 and then stop
            flag_break = True
      if flag_break:
        break
    elif anchor1 == "getall":   # read all lines that match anchor2
      if re.search(anchor2, line, flags=0):
        read_lines.append(line)
  if not flag_anch1:
    current_lines = read_lines
  if anchor1 != "getall" and anchor1 != "getany":
    lines = current_lines[int(num_skiplines) : num_lines ]
  elif anchor1 == "getall" or anchor1 == "getany":
    lines = read_lines
#  print("num_lines = ",num_lines)
#  print(read_lines)
#  print(current_lines)
#  print(lines)
  return lines

