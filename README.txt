# WaveformChesnokov
Chesnokov's AutoQRS algorithm, Linux version 1.01, May 8, 2015. created by Michael Shipway
This version has been stripped of Microsoft specific code and will compile on Linux using gcc.
It still has some HelloWorld scraps left over which should be removed soon.

usage: ecg.exe <filters dir> physionetfile.dat <output filepath> leadnum "verbose"
filters dir - directory containing the filter files
physionetfile.dat - full path/filename of the .hea file of the record,
                    (May also be full path/recordName)
leadnum - which lead to analyze, zero to analyze all available.
"verbose" - optional, if 5th parameter is the word "verbose", progress messages will be printed to console.
