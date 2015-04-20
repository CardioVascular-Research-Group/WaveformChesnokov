# WaveformChesnokov
This version has been stripped of Microsoft specific code and will compile on Linux using gcc.
It still has some HelloWorld scraps left over which should be removed soon.

usage: ecg.exe <filters dir> physionetfile.dat <output filepath> leadnum
filters dir - directory containing the filter files
physionetfile.dat - data file of the record,
   must be in the local directory, 
   or the directory in the environment variable WFDB,
   or the variable DEFWFDB (defined in /opt/wfdb/include/wfdb/wfdblib.h)
leadnum - which lead to analyze, zero to analyze all available.
