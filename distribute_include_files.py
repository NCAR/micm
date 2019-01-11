#shim for python 2
from __future__ import print_function

import json
import sys
import shutil
import os

troe_shim = """! Included shim as a number of WRF functions depend on this
!-------------------------------------------
! Troe equilibrium reactions (as in Stockwell et al, 1997)

    real(kind=r8) FUNCTION TROEE(A, B, k0_300K,n,kinf_300K,m,temp,cair)

    INTRINSIC LOG10

    real(kind=r8), INTENT(IN) :: temp      ! temperature [K]
    real(kind=r8), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    real(kind=r8),     INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    real(kind=r8),     INTENT(IN) :: n         ! exponent for low pressure limit
    real(kind=r8),     INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    real(kind=r8),     INTENT(IN) :: m         ! exponent for high pressure limit
    real(kind=r8),     INTENT(IN) :: A, B
    real(kind=r8)             :: zt_help, k0_T, kinf_T, k_ratio, troe


    zt_help = 300._r8/temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    troe   = k0_T/(1._r8+k_ratio)*0.6_r8**(1._r8/(1._r8+LOG10(k_ratio)**2))

    TROEE = A * EXP( - B / temp) * troe



  END FUNCTION TROEE
"""


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
  sys.exit(1)


# get name of mechanism
try:
  mechanism_name = data["modelNameInc"].replace('"', '');
  relative_source_directory_for_includes = "./generated/"+mechanism_name+"/"
except:
  print("Cannot find name of mechanism");
try:
  os.mkdir(relative_source_directory_for_includes)
except:
  print("Unable to make directory; Perhaps this mechanism has already been created?")
  sys.exit(1)


# names of all files to which data will be distributed
model_name_inc_file = relative_source_directory_for_includes + "model_name.inc"
forcing_inc_file = relative_source_directory_for_includes + "forcing.inc"
jacobian_inc_file = relative_source_directory_for_includes + "jacobian.inc"
k_rateConst_inc_file = relative_source_directory_for_includes + "k_rateConst.inc"
model_name_inc_file = relative_source_directory_for_includes + "model_name.inc"
molec_info = relative_source_directory_for_includes + "molec_info.json"
p_rates_inc_file = relative_source_directory_for_includes + "prates.inc"
rates_inc_file = relative_source_directory_for_includes + "rates.inc"
rates_functions_inc_file = relative_source_directory_for_includes + "rate_functions.inc"


# Model name file
try:
  model_name_inc_file_handle = open(model_name_inc_file,"w")
except:
  print("Unable to open "+model_name_inc_file+" in which the fortran snippet is to be placed. ")

# Write a string in the file
try:
  #print(data['forcingInc'])
  model_string = "model_name = '"+mechanism_name+"'";
  model_name_inc_file_handle.write(model_string)
except:
  print("Unable to place model name in "+model_name_inc_file)
  sys.exit(1)



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




# Rate Functions
try:
  rates_functions_inc_file = relative_source_directory_for_includes + "rate_functions.inc"
  rates_functions_inc_file_handle = open(rates_functions_inc_file, "w")
except:
  print("Unable to open "+rates_functions_inc_file+" in which the fortran snippet is to be placed. ")
  sys.exit(1)

# Write a string in the file
try:
  #print(data['ratesInc'])
  function_array = data['custom_rates']
  l = "! number of Functions: " + str(len(function_array)) + "\n\n"
  rates_functions_inc_file_handle.write(l)
except:
  print("Unable to write summary string to: "+  rates_functions_inc_file) 
  sys.exit(1)



try:
  for function in function_array :
    rates_functions_inc_file_handle.write( "\n\n")
    rates_functions_inc_file_handle.write( function['code'] )
  rates_functions_inc_file_handle.write( "\n\n")
  rates_functions_inc_file_handle.write( troe_shim )
except:
  print("Unable to write a particular rate_function to "+  rates_functions_inc_file)
  sys.exit(1)




# Suite Definition File
relative_source_directory_for_chemistry_suite = "../MusicBox_host/suites/";
chemistry_suite_file = relative_source_directory_for_chemistry_suite + mechanism_name+ ".xml";
suite_name="MICM_"+"mechanism_name"

try:
  chemistry_suite_file_handle = open(chemistry_suite_file,"w")
except:
  print("Unable to open chemistry suite definition file: "+ chemistry_suite_file)
  sys.exit(1)


header = """<?xml version=\"1.0\" encoding=\"UTF-8\"?>
"""

suite_name = "<suite name=\""+mechanism_name+"\" lib=\"micmchem\" ver=\"1.0.0\">"

footer= """
  <group name=\"time_vary\">
    <subcycle loop=\"1\">
      <scheme>mass_quantities_util</scheme>
      <scheme>k_rateConst</scheme>
      <scheme>tuv_photolysis</scheme>
      <scheme>photolysis_interstitial</scheme>
      <scheme>chemistry_driver_moz</scheme>
    </subcycle>
  </group>
  <finalize>mass_quantities_util</finalize>
  <finalize>k_rateConst</finalize>
  <finalize>tuv_photolysis</finalize>
  <finalize>photolysis_interstitial</finalize>
  <finalize>chemistry_driver_moz</finalize>
</suite>
"""
chemistry_suite_default= header+suite_name+footer;

print(chemistry_suite_default);

try:
  chemistry_suite_file_handle.write(chemistry_suite_default);
except:
  print("Unable to write chemistry suite definition to "+ chemistry_suite_file)
  sys.exit(1)



print("python script reached the end")
print("Look in "+relative_source_directory_for_includes+" for the resulting include files")
