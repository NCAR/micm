#shim for python 2
from __future__ import print_function

import json
import sys
import shutil

relative_source_directory_for_includes = "./generated/user-defined/"

forcing_inc_file = relative_source_directory_for_includes + "forcing.inc"
jacobian_inc_file = relative_source_directory_for_includes + "jacobian.inc"
k_rateConst_inc_file = relative_source_directory_for_includes + "k_rateConst.inc"
model_name_inc_file = relative_source_directory_for_includes + "model_name.inc"
molec_info = relative_source_directory_for_includes + "molec_info.json"
p_rates_inc_file = relative_source_directory_for_includes + "prates.inc"
rates_inc_file = relative_source_directory_for_includes + "rates.inc"

# Did user provide a file name?
if(len(sys.argv) > 1):
  source_data_file = sys.argv[1]
  print("Copying data from "+sys.argv[1])
else:
  print("Usage:  python "+sys.argv[0]+" /full/path/to/file/file.json")
  sys.exit(1)

# Open source data file
try:
  json_data = open(source_data_file)
except:
  print("souce data file: "+source_data_file+" could not be opened.");
  print("Usage:  python "+sys.argv[0]+" /full/path/to/source/data/file/file.json")
  sys.exit(1)


# Load data from source file as json (and close file)
try:
  data_complete = json.load(json_data)
  data = data_complete["mechanism"]
  json_data.close()
except:
  print("Unable to load json data from "+source_data_file)

#
# Extract forcing from JSON and put it in chosen directory
#


# Forcing
try:
  forcing_inc_file_handle = open(forcing_inc_file,"w")
except:
  print("Unable to open "+forcing_inc_file+" in which the fortran snippet is to be placed. ")

# Write a string in the file
try:
  #print(data['forcingInc'])
  forcing_inc_file_handle.write(data['forcingInc'])
except:
  print("Unable to place forcing in forcing.inc")
  sys.exit(1)



# Jacobian
try:
  jacobian_inc_file_handle = open(jacobian_inc_file,"w")
except:
  print("Unable to open "+jacobian_inc_file+" in which the fortran snippet is to be placed. ")
  sys.exit(1)

# Write a string in the file
try:
  #print(data['jacobianInc'])
  jacobian_inc_file_handle.write(data['jacobianInc'])
except:
  print("Unable to place jacobian in jacobian.inc")
  sys.exit(1)



# K Rate Constants (Rateconstants for non-photolysis reactions
try:
  k_rateConst_inc_file_handle = open(k_rateConst_inc_file,"w")
except:
  print("Unable to open "+k_rateConst_inc_file+" in which the fortran snippet is to be placed. ")
  sys.exit(1)

# Write a string in the file
try:
  #print(data['k_rateConstInc'])
  k_rateConst_inc_file_handle.write(data['k_rateConstInc'])
except:
  print("Unable to place k rate constants in k_rateConst.inc")
  sys.exit(1)



# molecule_info
try:
  shutil.copy2(sys.argv[1],molec_info)
except:
  print("Unable to copy json file to ./generated/user-defined/ ")
  sys.exit(1)




# TUV rates to prates mapping
try:
  p_rates_inc_file_handle = open(p_rates_inc_file,"w")
except:
  print("Unable to open "+p_rates_inc_file+" in which the fortran snippet is to be placed. ")
  sys.exit(1)

# Write a string in the file
try:
  #print(data['modelNameInc'])
  p_rates_inc_file_handle.write(data['pratesInc'])
except:
  print("Unable to place p_rate_mapping in "+p_rates_inc_file);
  sys.exit(1)


# Rates
try:
  rates_inc_file_handle = open(rates_inc_file,"w")
except:
  print("Unable to open "+rates_inc_file+" in which the fortran snippet is to be placed. ")
  sys.exit(1)

# Write a string in the file
try:
  #print(data['ratesInc'])
  rates_inc_file_handle.write(data['ratesInc'])
except:
  print("Unable to place rates in rates.inc")
  sys.exit(1)






