Input files here are from Chemistry Cafe. The tar file included contains a directory called `cesm`, which is where these files are from.
All of the files output from the tar file were renamed to TS1.<whatever-original-extension-was>

The only modification was the inclusion of the BEGSIM section, which will look like this in `input/TS1.in.cmake`:

```
BEGSIM
output_unit_number = 7
output_file        = musica-perf.out
procout_path       = ../musica/
src_path           = ${CAMCHEM_SRC_PATH}

... rest of the file ...

ENDSIM
```

This is a cmake file so that cmake can correctly overwrite `CAMCHEM_SRC_PATH` to point to the preprocessor source files...

Also, the input files MUST be in the input directory. the preprocessor tries to put output
into "../output" which apparently can be overriden by setting `sim_dat_path`, but that 
doesn't seem to work.