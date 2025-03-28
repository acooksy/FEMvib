#!/bin/python3
# This script parses Gaussian log files.  
# Converted and upgraded from Cooksy's es_parse.csh and femvib_parse.csh.
# 
# The strategy is to take each file and divide it into blocks for each distinct geometry
#  and/or method of interest.  As examples:
#   Each relaxed geometry optimization in a scan job is a separate block.
#   Each job in a series of Link1 steps from a single input file is a separate block, including opt+freq.
#   Each step in a rigid z-matrix-based scan is a separate block.
# For scans or groups of files with the same method, the option exists to create the
#  FEMvib geometry and coordinates files.
#
# Usage:  es_parse.py  [input]  {options}
#  input can be a single filename or a directory containing a set of files to be processed in one run
#
# Known failure modes:
#  -F: PES coordinates are interpreted as geometric parameters that appear as either "Scan" or "Frozen"
#      in log file.  That means that any other frozen parameters may be interpreted as PES coordiates.
#      These have to be removed after running es_parse.
#  -F -o: if optfile has the same extension and is in the same directory as the PES file(s), then it
#      will be analyzed twice (once as the optfile and once as one of the logfiles) and will appear
#      twice in the output files.

import argparse, extrpars, glob, os, re, sys

# initialize parser and read arguments
parser = argparse.ArgumentParser()
parser.add_argument("inputfile", help="location (directory or filename) with PES energies and geometries")
parser.add_argument("-b", "--basis_size", action='store_true', help="retrieve number of basis functions")
parser.add_argument("-c", "--counterpoise", action='store_true', help="retrieve counterpoise-corrected energy and BSSE")
parser.add_argument("-e", "--extension", help="change log file extension")
parser.add_argument("-F", "--FEMvib", action='store_true', help="create FEMvib coordinate and energyP files")
parser.add_argument("-f", "--frequencies", action='store_true', help="create a semicolon-delimited list of vibrational constants and IR intensities")
parser.add_argument("-g", "--geometries", action='store_true', help="create a semicolon-delimited list of atomic XYZ coordinates")
parser.add_argument("-H", "--HartreeFock", action='store_true', help="get HF energy")
parser.add_argument("-o", "--optfile", help="file with optimized geometry (for FEMvib input set-up)")
parser.add_argument("-s", "--stoichiometry", action='store_true', help="retrieve the molecular formula, charge, and multiplicity")
parser.add_argument("-t", "--thermochemistry", action='store_true', help="retrieve the ZPE, enthalpy, and Gibbs free energy")
parser.add_argument("-z", "--zmatrix", action='store_true', help="calculate nearest-neighbor distances, angles, and dihedrals")

args = parser.parse_args()

# set the input file locations
inputfiles = []
inputfile = args.inputfile
if args.extension is not None:   # change the input file extension from "log"
  ext = args.extension
else:
  ext = "log"
file_root = os.path.splitext(inputfile)[0]
if os.path.isdir(inputfile):
  if "/" in file_root:
    file_root = file_root.strip("/")   # if inputfile is a directory, save to [dir].csv, not [dir]/.csv
  inputdir = './' + inputfile
  inputfiles = glob.glob(inputdir+'/*.'+ext)  
else:
  if os.path.isfile(inputfile):
    inputfiles.append(inputfile)
  else:
    print("Input file location ", inputfile," not found.  Exiting.")
    sys.exit(1)
num_files = len(inputfiles) 
if num_files == 0:
  print("No input files found.  Exiting.")
  sys.exit(1)
print("inputfiles = ", inputfiles)

# set other filename variables, option flags
flag_b = args.basis_size
flag_c = args.counterpoise
flag_f = args.frequencies
flag_F = args.FEMvib
flag_g = args.geometries
flag_H = args.HartreeFock
flag_o = False
if args.optfile is not None:
  flag_o = True
flag_s = args.stoichiometry
flag_t = args.thermochemistry
flag_z = args.zmatrix
if flag_g:
  geomsfile = file_root + "_geoms.csv"
if flag_f:
  freqsfile = file_root + "_freqs.csv"
if flag_F:
  if flag_o:
    optfile = args.optfile
    if os.path.isfile(optfile):
      print("Optimized geometry will be taken from ", optfile)
    elif os.path.isdir(inputfile):
      optfile = inputdir + "/" + optfile 
      if os.path.isfile(optfile):
        print("Optimized geometry will be taken from ", optfile)
    else:
      print("Optimization job file ", optfile, " not found.  Exiting.")
      sys.exit(1)
  energyfile = file_root + ".energies"
  coordsfile = file_root + ".coords"
# set output filenames
outfile = file_root + ".csv"

# Initialize keystrings (the hope is these can be changed for different programs)
kstrng_jobstart = "Entering Gaussian System"
kstrng_jobstop = "Normal termination"
kstrng_atoms = "NAtoms="
kstrng_basis = "basis functions\,"
kstrng_bsse = "BSSE energy"
kstrng_cntrps = "Counterpoise corrected energy"
kstrng_cp_ver = "doing calculation for fragment   1 using the basis set of the full-system"
kstrng_stoich = "Stoichiometry"
kstrng_stat = "Stationary"
kstrng_geom = "Standard orientation"
kstrng_freq = "Frequencies --"
kstrng_IRint = "IR Inten    --"
kstrng_ZPE = "Zero-point correction="
kstrng_thermo = "Sum of electronic and thermal Energies="
noanchor = "text that doesnt appear anywhere in any of the logfiles" # this is to read to end of the block, since we don't know what strings may terminate the block
num_field = 0

# Open output files as needed
file_out = open(outfile,"w")
if flag_g:
  file_geoms = open(geomsfile,"w")  
if flag_f:
  file_freqs = open(freqsfile,"w")  
# Write CSV file header based on flags
file_out.write("filename; jobtype; method; ")
if flag_s:
  file_out.write("stoichiometry; ") 
if flag_b:
  file_out.write("contracteds; primitives; ") 
if flag_H:
  file_out.write("HF energy; ") 
file_out.write("energy; ") 
if flag_c:
  file_out.write("counterpoise; BSSE; ") 
if flag_t:
  file_out.write("zero-point; thermal E; enthalpy; free energy; ") 
file_out.write("status; \n")

# Initialize lists for FEMvib coords and energies
energies = []   # list of energies for all blocks
geometries = []   # list of geometry tuples for all blocks
geometry_opt = []
params = []   # list of scan coordinates 
coord_values = []   # list of coordinate tuples for all blocks
if flag_o:
  inputfiles.append(optfile)  
num_null_files = 0

# Loop over inputfiles and break into blocks
num_tot_files = 0
num_tot_blocks = 0
for file in inputfiles:
  print("processing file ",file)

# Separate the current file into blocks as described above.
# First level: distinct Link1 jobs, distinct steps in a scan job, the 2 halves of opt+freq 
#   appear in the logfile as separate jobs.
  with open(file) as f:
    jobs = extrpars.extract_jobs(file)
    for i, job in enumerate(jobs):
      print(f"processing job {i+1}:")

# First pass through job to get jobtype and number of atoms (determine geometry format), 
#  method (determines energy format)
      flag_cp_ver = False
      method, jobtype, flag_cp_ver, num_atoms_str = extrpars.parse_job_method(job)
      num_atoms = int(num_atoms_str)
      print("Method = ",method)
      match jobtype:
        case "sp":
          print("Apparent job type: single-point energy.")
        case "opt":
          print("Apparent job type: vanilla geometry optimization.")
#        case "opt_freq":
#          print("Apparent job type: geometry optimization + frequency calculation.")
        case "opt_mod":
          print("Apparent job type: relaxed scan using redundant internal coordinates.")
        case "opt_zmat":
          print("Apparent job type: relaxed scan using stacked Z-matrix optimization jobs.")
        case "scan_zmat":
          print("Apparent job type: rigid scan using Z-matrix.")
        case "freq":
          print("Apparent job type: vibrational frequency calculations.")

# Second level: divide each job (if necessary) into distinct blocks such that each block
#  gives results for one geometry with one method and basis set
      blocks = extrpars.extract_blocks(job,jobtype)
      print("Found ", len(blocks), " block(s) in job.")
      if len(blocks) == 0:
        num_null_files += 1
      ilines = 0
      for j, block in enumerate(blocks):
        print(f"Processing block {j+1}:")
        num_lines_block = len(block)
#        print("# lines in block = ",num_lines_block)

# write filename, jobtype, method to output
        file_out.write("%s; %s; %s; " % (file, jobtype, method))

# get stoichiometry
        if flag_s:
          anchor1 = "getany"
          anchor2 = kstrng_stoich
          num_skiplines = 0
          num_readlines = 1
          lines = extrpars.extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines)
          file_out.write("%s; " % lines[0].split()[1])

# get numbers of basis functions
        if flag_b:
          anchor1 = "getany"
          if flag_cp_ver:
            anchor1 = kstrng_cp_ver
          anchor2 = kstrng_basis
          num_skiplines = 0
          num_readlines = 1
          lines = extrpars.extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines)
          file_out.write("%s; %s; " % (lines[0].split()[0], lines[0].split()[3] )) # contracteds, primitives

# get energy, including HF energy if requested
        if "opt" in jobtype:
          anchor1 = "Stationary"
        elif "opt" not in jobtype:
          anchor1 = noanchor
# For counterpoise, the last E is for a fragment so we use a different anchor string to get cluster energy
        if flag_cp_ver:
          anchor1 = kstrng_cp_ver
        if flag_H:
          anchor2 = r'E\(.*HF\)'
          num_field = 4
          num_skiplines = 0
          num_readlines = 1
          lines = extrpars.extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines)
          if len(lines) == 1:  # DFT jobs may not report the HF energy, so check
            file_out.write("%s; " % lines[0].split()[num_field]) 
          else:
            file_out.write("; ") 
        if method == "HF":   # selecting flag_H when method=HF will just print the HF energy twice
          anchor2 = r'E\(.*HF\)'
          num_field = 4
        if method == "MN15":
          anchor2 = r'E\(.*MN15\)'
          num_field = 4
        if method == "B3LYP":
          anchor2 = r'E\(.*B3LYP\)'
          num_field = 4
        if method == "MP2":
          anchor2 = r'E.*MP2 ='
          num_field = 5
        if method == "QCISD":
          anchor2 = r'E\(Corr\)='
          num_field = 4
        if method == "CCSD":
          anchor2 = r'E\(Corr\)='
          num_field = 4
        if method == "CCSDT":
          anchor2 = r'CCSD\(T\)= '  # the trailing space is to avoid finding archor2 in the archive block
          num_field = 1
        num_skiplines = 0
        num_readlines = 1
        lines = extrpars.extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines)
        # Sometimes the energy is missing; at least one case with TD-DFT jobs
        print("lines = ",lines)
        if lines == []:
          print("Missing energy for this block")
          energy = ""
        else:
          energy_str = lines[0].split()[num_field]
        # if energy is in D exponential form (e.g., MP2 energies), change to E form for Python
          if re.search("D", energy_str, flags=re.IGNORECASE):
            t = str.maketrans('D','E')
            energy_str = energy_str.translate(t)
          energy = float(energy_str)
        file_out.write("%s; " % energy)
        if flag_F and jobtype != "freq" and (not flag_c or not flag_cp_ver):  # for FEMvib, use the counterpoise energy if available
          energies.append(float(lines[0].split()[num_field]))
          if flag_o and file == optfile:
            energy_opt = float(lines[0].split()[num_field])

# get PES coordinates for FEMvib files
        if flag_F and jobtype != "freq":
          current_values = []  # the list of PES coord values for the current block
          if jobtype == "scan_zmat":
            anchor1 = "Scan the potential surface"
            anchor2 = "Variables:"
            num_skiplines = 1
            num_readlines = 3 * num_atoms - 5  # not strictly right but should be adequate
            lines = extrpars.extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines)
            for line in lines:
              if "Scan" in line:
                num_field = 0
                param_i = line.split()[num_field]
                params.append(param_i)
                num_field = 1
                value = float(line.split()[num_field]) # scan_zmat gives first set of values in unique format so we read it here
                if param_i[0] == "D" or param_i[0] == "A": # convert any angles in initial values to radians
                  value = value / 57.2957795
                current_values.append(value) 
            if current_values != []: # This should be triggered only on the first block in the scan_zmat file
              coord_values.append(current_values)
          if jobtype == "opt_zmat" or jobtype == "opt_mod":
            if params == []:   # This should be triggered only on the first block in the opt file
              anchor1 = "Number of steps in this run="
              anchor2 = "!    Initial Parameters    !"
              num_skiplines = 5
              num_readlines = 3 * num_atoms - 5  # not strictly right but should be adequate
              num_field = 1
              lines = extrpars.extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines)
              for line in lines:
                if "Scan" in line:
                  params.append(line.split()[num_field])
                elif "Frozen" in line:
                  params.append(line.split()[num_field])
          num_params = len(params)
  # check if any of the coordinates are angular, since these need to be converted to radians
          flag_ang = []
          for i,param_i in enumerate(params):
            flag_ang.append(False)
            if param_i[0] == "D" or param_i[0] == "A":
              flag_ang[i] = True
            
  # extract the values of the PES coordinates for this block and add to a running list
          current_values = []
          if jobtype == "scan_zmat":
            anchor1 = "Z-MATRIX \(ANGSTROMS AND DEGREES\)"
            anchor2 = "Variable"
            num_skiplines = 2
            num_readlines = num_params 
            num_field = 2
            lines = extrpars.extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines)
            for i, line in enumerate(lines):
              if flag_ang[i] == False:
                value = float(line.split()[num_field]) * 0.52918
              elif flag_ang[i] == True:
                value = float(line.split()[num_field]) 
              current_values.append(value)
          if jobtype != "scan_zmat":
#          if jobtype == "opt_zmat" or jobtype == "opt_mod":
            for line in block:
              for i,param in enumerate(params):
                search_string = r"!.*" + param + ".*-DE/DX" 
                if re.search(search_string, line, flags=re.IGNORECASE):
                  num_field = 3
                  if jobtype == "opt_zmat":
                    num_field = 2
                  if flag_ang[i] == False:
                    value = float(line.split()[num_field])
                  elif flag_ang[i] == True:
                    value = float(line.split()[num_field]) / 57.2957795
                  current_values.append(value)
          if current_values != []:
            coord_values.append(current_values)
    
#          print("params = ",params," num = ",num_params,"  ang = ",flag_ang)
#          print("current_values = ",current_values)
#          print("coord_values = ",coord_values)

# get counterpoise-corrected energy and BSSE
        if flag_c:
          anchor1 = "complexation energy ="
          anchor2 = kstrng_cntrps   # for counterpoise energy
          num_skiplines = 0
          num_readlines = 1
          num_field = -1
          lines = extrpars.extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines)
          if len(lines) == 1:  # flag_c may be set where some jobs are not counterpoise 
            file_out.write("%s; " % lines[0].split()[num_field]) 
          else:
            file_out.write("; ") 
          if flag_F and flag_cp_ver and jobtype != "freq":
            energies.append(float(lines[0].split()[num_field]))
            if flag_o and file == optfile:
              energy_opt = float(lines[0].split()[num_field])
          anchor2 = kstrng_bsse   # for counterpoise energy
          lines = extrpars.extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines)
          if len(lines) == 1:  # flag_c may be set where some jobs are not counterpoise 
            file_out.write("%s; " % lines[0].split()[num_field]) 
          else:
            file_out.write("; ") 

# get thermochem
        if flag_t: 
          if jobtype == "freq":
            anchor1 = "getany"
            anchor2 = kstrng_ZPE
            num_skiplines = 0
            num_readlines = 1
            num_field = -2
            lines = extrpars.extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines)
            for line in lines:
              file_out.write("%s; " % (line.split()[num_field]))   # zero-point energy; use next-to-last last field 
            anchor2 = kstrng_thermo
            num_skiplines = 0
            num_readlines = 3
            num_field = -1
            lines = extrpars.extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines)
            for line in lines:   # thermal energy; enthalpy; Gibbs free energy; use last field
              file_out.write("%s; " % (line.split()[num_field]))
          elif flag_t and jobtype != "freq":
            file_out.write("; ; ; ;") 

# get geometry
        if flag_g or flag_F:
          if flag_g:
            file_geoms.write("%s; %s; %s; " % (file, jobtype, method))
          anchor1 = noanchor
          anchor2 = kstrng_geom
          num_skiplines = 5
          num_readlines = num_atoms
          lines = extrpars.extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines)
          for line in lines:
            if flag_g:
              file_geoms.write("\n%s; %s; %s; %s" % (line.split()[1],line.split()[3],line.split()[4],line.split()[5]))   
  # also save the geometries for FEMvib coordinates file
            if flag_F and jobtype != "freq":
              geometries.append((line.split()[3],line.split()[4],line.split()[5]))   
    # if this is the optimized geometry, extract that geometry 
              if flag_o and file == optfile:
                geometry_opt.append((line.split()[1],line.split()[3],line.split()[4],line.split()[5]))   
          if flag_g:
            file_geoms.write("\n")
    # if this is the optimized geometry, extract the coordinate values
        if flag_F and flag_o and file == optfile and jobtype != "freq":  
          anchor1 = noanchor
          anchor2 = kstrng_stat
          num_skiplines = 6
          num_readlines = 3 * num_atoms - 5
          lines = extrpars.extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines)
          for line in lines:
            for i,param in enumerate(params):
              search_string = r"!.*" + param + ".*-DE/DX" 
              num_field = 3
              if jobtype == "opt_zmat":
                num_field = 2
              if re.search(search_string, line, flags=re.IGNORECASE):
                current_value = float(line.split()[num_field])
                if flag_ang[i]:
                  current_value = current_value / 57.2957795
#                current_values.append(current_value)
#                print("current_values = ",current_values)
          coord_values.append(current_values)

# get frequencies
        if flag_f and jobtype == "freq":
          file_freqs.write("%s; %s; %s; " % (file, jobtype, method))
          anchor1 = "getall"
          anchor2 = kstrng_freq
          num_skiplines = 0
          num_readlines = 1
          lines = extrpars.extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines)
          file_freqs.write("\n")
          for line in lines:  # There are 1, 2, or 3 freqs per line
            if len(line.split()) == 3:
              file_freqs.write("%s; " % (line.split()[2]))
            if len(line.split()) == 4:
              file_freqs.write("%s; %s ;" % (line.split()[2],line.split()[3]))
            if len(line.split()) == 5:
              file_freqs.write("%s; %s; %s; " % (line.split()[2],line.split()[3],line.split()[4]))   
          file_freqs.write("\n")
# repeat to get intensities
          anchor2 = kstrng_IRint
          lines = extrpars.extract_lines(block,anchor1,anchor2,num_skiplines,num_readlines)
          for line in lines:
            if len(line.split()) == 4:
              file_freqs.write("%s; " % (line.split()[3]))
            if len(line.split()) == 5:
              file_freqs.write("%s; %s; " % (line.split()[3],line.split()[4]))
            if len(line.split()) == 6:
              file_freqs.write("%s; %s; %s; " % (line.split()[3],line.split()[4],line.split()[5]))   
          file_freqs.write("\n")
 
# last task for current block: check if any errors detected
        status = "completed"
        if energy == []:
          status = "incomplete"
        for line in block:
          if "Error termination" in line:
            status = "crashed"
        file_out.write("%s \n" % status)

        num_tot_blocks += 1
      # end loop over blocks

    # end loop over jobs

  num_tot_files += 1
# end loop over files

# FEMvib file set-up
if flag_F:
  # Write files output from stored lists
  file_energy = open(energyfile,"w")  
  file_coords = open(coordsfile,"w")  
  file_energy.write("%i\n" % num_tot_blocks)
  file_coords.write("%i %i %i\n" % (num_params, num_atoms, num_tot_blocks))
   
  if flag_o:
    energy_min = energy_opt
  else:
    energy_min = min(energies)
#  print("energy_min = ",energy_min)
#  print("energies= ",energies)
#  print("coord_values= ",coord_values)
  if flag_o:
    for j in range(num_atoms):
      file_coords.write(' '.join(str(s) for s in geometry_opt[j]) + "\n")
#  print("geometries = ",geometries)
  for i, energy in enumerate(energies):
#    print("i coord_values energy =",i,coord_values[i])
    file_energy.write(' '.join(str(s) for s in coord_values[i]) + " ")
    file_coords.write(' '.join(str(s) for s in coord_values[i]) + "\n")
    energy_cm = (energy - energy_min) * 219475
    file_energy.write("%f\n" % energy_cm)
    j = 0
    for j in range(num_atoms):
#      file_coords.write("%i %i %i \n" % (i, j, num_atoms))
      file_coords.write(' '.join(str(s) for s in geometries[num_atoms*i+j]) + "\n")
  print("Coordinates saved to ",coordsfile)
  print("Energies saved to ",energyfile)
  if flag_o:
    print("  NOTE: you must change atomic numbers in ",coordsfile," to atomic masses.")
  file_coords.close()
  file_energy.close()

# close out
file_out.close()
print("Output data saved to ",outfile)
print(num_tot_files," file(s) processed in total.")
if num_null_files > 0:
  print(num_null_files, " files were NOT parsed; check for errors.")
print(num_tot_blocks," block(s) processed in total.")

sys.exit(0)

"""
Notes:
Based on the cases below, to etract a relevant string of data we need to call 
the extraction routine with 
 anchor1 = keystring whose last occurrence signals completion of current calculation, 
   often report of optimized geometry; if anchor1 is empty, use bottom of block as termination of search;
   also set anchor1 = "getall" to simply collect *all* lines that contain string anchor2
 anchor2 = keystring whose last occurrence *before* anchor1 signals start of section with
   appropriate data (needed for geometries, energies)
 num_skiplines = number of lines to skip after anchor2 (needed for geometries)
 num_readlines = number of lines to read
Ideally the extraction routine returns all and only the appropriate lines from the file,
 and we subsequently pull the fields from those lines (so that the extraction routine can be as
 general as possible).

Cases: 
independent jobs end with "Normal termination" including opt in opt+freq job
              begin with " Entering Gaussian System" 
  Beware: opt job for opt+freq terminates with "Normal termination" but 
          freq job opens with " Link1:  Proceeding to internal job step number"
          Also, freq job in opt+freq does *not* have its own keyword line " # freq ...".

job types:
  opt:
  opt=z-mat with or without scanned params, no relaxed params:
    geom in last "Standard orientation" prior to "Stationary", skip 5 lines, read next NAtoms lines; fields 2 4 5 6
  scan (rigid scan using Z-matrix):
    no "Stationary" statements
    geom in each "Standard orientation", skip 5 lines; fields 2 4 5 6
  freq
    only job type with freq and thermochem data
    not necessarily any geometry data -- use other jobs for that
    to read freq data, all lines with " Frequencies --"; fields 3 4 5
                       all lines with " IR Inten    --"; fields 4 5 6 
  thermochem
    only in freq job type
    read each line with:
      " Zero-point correction="; field 3
      " Sum of electronic and thermal Energies="; field 7   and next two lines:
      " Sum of electronic and thermal Enthalpies="; field 7
      " Sum of electronic and thermal Free Energies="; field 8

methods:
  HF:
    opt:
    opt=z-mat:
      energy in last "E(.*HF" prior to "Stationary"; field 5
    scan
      energy in each "E(.*HF"; field 5
  DFT:
    energy in last "E(.*B3LYP"; field 5
  MP2:
    energy in last "EUMP2 = "; field 6
  QCISD:
    energy in last "E(Corr)="; field 5
  CCSD
    energy in last "E(Corr)="; field 5
  CCSD(T)
    energy in last "CCSD(T)="; field 2

FEMvib coords and energies
  scan_zmat:
    Gaussian jobtype=scan is a mess for parsing.  Using nearly the same format, it prints the 
      first value of the variables in A and degrees, and thereafter in a0 and radians.  I assume 
      that the expectation is that one would read the geomatry and energy data from the summary 
      at the end of the scan job, but that would mean (a) no way to process crashed jobs; (b) no 
      way to extract other properties such as dipole moments at each PES point; (c) treating scan 
      jobs with a completely different algorithm than all other jobtypes.

misc:
  NAtoms
    any " NAtoms="; field 2
  Stoichiometry
    any " Stoichiometry"; field 2
  Scanned parameters (PES coords)
    opt=modredundant
      any " ! .* Scan"; field 2 (param name) 4 (initial value)
      any " ! [param name]"; field 4 (new value)
"""
 
