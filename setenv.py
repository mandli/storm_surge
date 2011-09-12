#
# setenv.py 
#
# Creates a shell script that sets up the environment to use the multi-layer
# shallow water code.  It is assumed that the Clawpack environment has already
# been appropriately set.
#
# Environment variables set:
#   DATA_PATH = Path where the output is put, defaults to the current
#     directory.  If a variable with this name already exists it is not set.
#   ML_PATH = Path to this code repository, defaults to current directory from
#     which this script has been run
#   ML_SRC = Path to the source files
#   ML_PYTHON = Path to the python files, this will be appended to the 
#     environment variable PYTHONPATH

import os
import sys

def write_csh_line(file_handle,var,value):
    file_handle.write('setenv %s "%s"\n' % (var.upper(),value))

def write_bash_line(file_handle,var,value):
    file_handle.write('export %s="%s"\n' % (var.upper(),value))


if __name__ == "__main__":
    print " "
    print "------------------------------------------------------------"
    
    # Dictionary of environment variables to be written out
    var_dict = {}
    
    # Set base directory path
    base_path = os.path.abspath(os.curdir)
    if len(sys.argv) > 1:
        base_path = os.path.abspath(sys.argv[1])
    print "Full path to multi-layer code directory should be:"
    print "      $ML_PATH = ",base_path
    var_dict['ML_PATH'] = base_path
    
    # Set DATA_PATH variable if not present
    if len(sys.argv) > 2:
        data_path = os.path.abspath(sys.arv[2])
    else:
        data_path = "./"

    # PYTHONPATH modification
    ml_python_path = os.path.join(base_path,"src","python")
    python_path = ":".join((ml_python_path,"${PYTHONPATH}"))
    var_dict['PYTHONPATH'] = python_path
    var_dict['ML_PYTHON'] = ml_python_path
    
    # Source directory
    var_dict['ML_SRC'] = os.path.join(base_path,"src")

    # Write out files
    csh_file = open(os.path.join(base_path,"setenv.csh"),'w')
    bash_file = open(os.path.join(base_path,"setenv.bash"),'w')
    for (var,value) in var_dict.iteritems():
        write_bash_line(bash_file,var,value)
        write_csh_line(csh_file,var,value)
    
    # Add DATA_PATH check
    csh_file.write('if ($?DATA_PATH == 0) then\n')
    csh_file.write("    setenv 'DATA_PATH' '%s'\n"% data_path)
    csh_file.write('endif\n')
    bash_file.write('if [ -z "${DATA_PATH}" ]; then\n')
    bash_file.write('export DATA_PATH="%s"\n'% data_path)
    bash_file.write('fi\n')
    
    csh_file.close()
    bash_file.close()

    print "------------------------------------------------------------"
    print "The files setenv.csh and setenv.bash contain the appropriate"
    print "commands to set environment variables for csh or bash shells"
    print "  and also some aliases you may find convenient             "
    print "------------------------------------------------------------"
    print " "

    sys.exit(0)
