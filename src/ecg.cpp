// ecg.cpp : Defines the entry point for the console application.
// Modified by Stephen Granite (sgranite@jhu.edu) to extend functionality
// version 1.0:  Computes corrected QT along with QT and RR statistical functions for all leads in the input set
// version 1.1:  Stores computed values in XML format (autoqrs.xsd schema) to allow for query and cross-comparison
// version 1.2:  Added filters directory path to command line parameters and added computation of QTVI_log to output
// version 1.3:  Revised code to compute HR and RR statistics properly
// version 1.4:  Revised code to compute minimums and maximums for QT, HR and RR, and QT Dispersion
// version 1.5:  Revised code to add bioportal URLs to tag the values
// version 1.6:  Revised code to replace all sr values with 1000, as it is that value for Chesnokov computations, rather than the sampling rate
// version 1.65: Revised code for clarity, converting sample-rate to msec or Heart Rate for all values before mean, variance and standard deviation are calculated.
//               Added documentation and regrouped some variable declarations.

#include "lib/stdafx.h"
#include "lib/lib.h"


//#include </opt/wfdb/include/wfdb/wfdb.h>

void help();
//int runChesnokov(char *fltdir, char *physionetdatfile, char *outputfilepath, int leadNumber);
int runChesnokov(char *fltdir, char *physionetdatfile, char *outputfilepath);
class Signal signal;
bool verbose;

vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

int main(int argc, char* argv[])
{
	verbose = false;
	if(argc>4){ //any value for argv[4] works
		verbose = true;
	}
	if(verbose) {
		printf("\n|---------------------------------------------------------------------------------------------|\n"
			"| Chesnokov's AutoQRS algorithm, Linux version 1.3, May 13, 2015. created by Michael Shipway |\n"
			"| Chesnokov's AutoQRS algorithm, Linux version 1.65x, December 2, 2015. modified by Stephen Granite |\n"
			"|---------------------------------------------------------------------------------------------|\n"
			"ecg main()    \n "
			"argv[0], <command>: %s  \n "
			"argv[1], <filters dir>: %s  \n "
			"argv[2], <physionetfile.dat>: %s  \n "
			"argv[3], <output filepath>: %s  \n"
			"argv[4], <verbose>: %s  \n",
			argv[0], argv[1], argv[2], argv[3], argv[4]);
	}
//***********************************************************
    if( argc < 4 ) {
        help();
    }
    else
    {
        int leadNumber = 0;
        if( argc == 4+1 ) {
            leadNumber =  *argv[4]  - 1;
            if( leadNumber < 0 ) leadNumber = 0;
        }

        vector<string> datFileInfo = split(argv[2],'/');
		string path = "";
		for (unsigned int i = 1; i < datFileInfo.size()-1; ++i) {
			path.append("/" + datFileInfo[i]);
		}
		char *charPath = new char[path.length() + 1];
		strcpy(charPath, path.c_str());
		setwfdb(charPath);
		char *datFileName = new char[datFileInfo[datFileInfo.size()-1].length()+ 1];
		strcpy(datFileName, datFileInfo[datFileInfo.size()-1].c_str());
		signal.PopulateSignal(datFileName, verbose);
        runChesnokov(argv[1],argv[2],argv[3]); //<filters dir> , physionetfile.dat, <output filepath>, leadnum
    }
}

int runChesnokov(char *fltdir, char *physionetdatfile, char *outputfilepath){
	{
        if( signal.GetLeadsNum() != 0){
        	if(signal.GetLength()*(1000/(int)signal.GetSR()) > 1000){ // must be greater than 1 second duration
				int leadNumber=0;
				int annTotalCount = 0;
				FILE *fp_output = fopen( outputfilepath, "wt" );
				int size = signal.GetLength();

				double sr = signal.GetSR();
				int h,m,s,ms;
				int msec = int( ((long double)size/sr) * 1000.0 );
				signal.mSecToTime(msec,h,m,s,ms);
				wchar_t leads[18][6] =  {L"I",L"II",L"III",L"aVR",L"aVL",L"aVF",L"v1",L"v2",L"v3",L"v4",L"v5",L"v6",L"MLI",L"MLII",L"MLIII",L"vX",L"vY",L"vZ"};

				fwprintf( fp_output,L"<?xml version = \"1.0\" encoding = \"UTF-8\"?>\n");
				fwprintf( fp_output,L"<autoQRSResults xmlns=\"http://www.cvrg.org/1/AutoQRSService\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.cvrg.org/1/AutoQRSService\">\n");
				fwprintf( fp_output,L"\t<jobIdentifier/>\n");
				fwprintf( fp_output,L"\t<FileAnalyzed>%s</FileAnalyzed>\n", physionetdatfile);
				fwprintf( fp_output,L"\t<LeadCount>%d</LeadCount>\n", signal.GetLeadsNum() );
				fwprintf( fp_output,L"\t<SR>%.2f</SR>\n", sr );
				fwprintf( fp_output,L"\t<UmV>%d</UmV>\n", signal.GetUmV() );
				fwprintf( fp_output,L"\t<Length>%02d:%02d:%02d.%03d</Length>\n", h,m,s,ms );

				for( leadNumber = 0; leadNumber < signal.GetLeadsNum(); leadNumber++) {
					signal.GetData( leadNumber );

					if(verbose) {
						printf("\tlead: %ls\n", leads[leadNumber] );
						printf("\tsamplingRate: %.2f\n", sr );
						printf("\tbits: %d\n", signal.GetBits() );
						printf("\tumv: %d\n", signal.GetUmV() );
						printf("\tsize: %d\n", size );
						printf("\tduration: %02d:%02d:%02d.%03d\n", h,m,s,ms );
					}

					long double* data = signal.GetData();

					//annotation
					class Annotation ann;  //default annotation params

					//or add your custom ECG params to annotation class from lib.h
					// ANNHDR hdr;
					//  hdr.minbpm = 30;
					//  etc...
					// class Annotation ann( &hdr );

					fwprintf( fp_output,L"\t<LeadResults id=\"%d\">\n", leadNumber + 1 );

					if(verbose) { printf( "\tgetting QRS complexes... \n"); }
					int** qrsAnn = ann.GetQRS( data, size, sr, fltdir );       //get QRS complexes

					if( qrsAnn ) {
						fwprintf( fp_output,L"\t\t<Lead bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000000150\">%ls</Lead>\n", leads[leadNumber] );
						fwprintf( fp_output,L"\t\t<TotalBeatCount bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001083\">%d</TotalBeatCount>\n", ann.GetQRSnum() );
						ann.GetECT( qrsAnn, ann.GetQRSnum(), sr );      //label Ectopic beats
						fwprintf( fp_output,L"\t\t<EctopicBeatCount bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001084\">%d</EctopicBeatCount>\n", ann.GetECTnum() );
						if(verbose) { printf("\tgetting P, T waves... \n" ); }
						int annNum = 0;
						int** ANN = ann.GetPTU( data, size, sr, fltdir, qrsAnn, ann.GetQRSnum() );   //find P,T waves
						if(ANN) {
							annNum = ann.GetANNnum();
							if(verbose) { printf("\t%d annotations detected.\n\n", annNum ); }
						}
						else {
							ANN = qrsAnn;
							annNum = 2 * ann.GetQRSnum();
							printf("\tP, T waves detection failed.\n" );
						}
						annTotalCount += annNum;

						char anncodes [51][10] =  {"notQRS", "N",       "LBBB",    "RBBB",     "ABERR", "PVC",
									"FUSION", "NPC",     "APC",     "SVPB",     "VESC",  "NESC",
									"PACE",   "UNKNOWN", "NOISE",   "q",        "ARFCT", "Q",
									"STCH",   "TCH",     "SYSTOLE", "DIASTOLE", "NOTE",  "MEASURE",
									"P",      "BBB",     "PACESP",  "T",        "RTM",   "U",
									"LEARN",  "FLWAV",   "VFON",    "VFOFF",    "AESC",  "SVESC",
									"LINK",   "NAPC",    "PFUSE",   "(",        ")",     "RONT",

			//user defined beats//
									"(p",     "p)",      "(t",      "t)",       "ECT",
									"r",      "R",       "s",       "S"};

						//printing out annotation
						/*
						for( int i = 0; i < annNum; i++ ) {
							int smpl = ANN[i][0];
							int type = ANN[i][1];

							msec = int( ((long double)smpl/sr) * 1000.0 );
							signal.mSecToTime(msec,h,m,s,ms);
							//fwprintf( fp_output,L"  %02d:%02d:%02d.%03d   %s\n", h,m,s,ms, anncodes[type] );
						}
						*/
						//saving RR seq
						vector<long double> rrs;
						vector<int> rrsPos;
						vector<long double> qts;
						vector<int> qtsPos;
						double meanRR = 0.0, meanQT = 0.0, QTc = 0.0;
						double minQT = 10000000000.0, maxQT = 0.0, minRR = 10000000000.0, maxRR = 0.0;
						double minHR = 10000000000.0, maxHR = 0.0;
						double varRR = 0.0, sdRR = 0.0, varQT = 0.0, sdQT = 0.0;
						double meanHR = 0.0, sdHR = 0.0, varHR = 0.0;
						//double dNaN = std::numeric_limits<double>::quiet_NaN();

						if( ann.GetRRseq( ANN, annNum, sr, &rrs, &rrsPos) ) {
							for( int i = 0; i < (int)rrs.size(); i++ ) {
								meanHR = meanHR + rrs[i];
								if( rrs[i] > maxHR) maxHR = rrs[i];
								if( rrs[i] < minHR) minHR = rrs[i];
								meanRR = meanRR + 60000/rrs[i];
							}
							maxRR = 60000/minHR;
							minRR = 60000/maxHR;
							//fclose( fp );
						}
						meanRR = meanRR / rrs.size();

						for( int i = 0; i < (int)rrs.size(); i++ ) {
							varHR = varHR + (rrs[i] - meanHR) * (rrs[i] - meanHR);
							varRR = varRR + (((sr * (60/rrs[i])) - meanRR) * ((sr * (60/rrs[i])) - meanRR));
						}
						varRR = varRR / (rrs.size() - 1);
						varHR = varHR / (rrs.size() - 1);
						sdRR = sqrt(varRR);
						sdHR = sqrt(varHR);

						if( ann.GetQTseq( ANN, annNum, sr, &qts, &qtsPos) ) {
							//FILE *fp = _wfopen( L"qts.txt", L"wt" );
							for( int i = 0; i < (int)qts.size(); i++ ) {
								//fwprintf( fp, L"%Lf\n", qts[i] );
								meanQT = meanQT + qts[i];
								if( qts[i] > maxQT) maxQT = qts[i];
								if( qts[i] < minQT) minQT = qts[i];

							}
							//fclose( fp );
						}
						ann.GetQTseq( ANN, annNum, sr, &qts, &qtsPos);
						meanQT = meanQT / qts.size();
						for( int i = 0; i < (int)qts.size(); i++ )
							varQT = varQT + (qts[i] - meanQT) * (qts[i] - meanQT);
						varQT = varQT / (qts.size() - 1);
						sdQT = sqrt(varQT);

						if( (qts.size() < 1) || (rrs.size() < 1) ) {
							fwprintf( fp_output,L"\t\t<QTCorrected_Bazett bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000000701\">NaN</QTCorrected_Bazett>\n");
							fwprintf( fp_output,L"\t\t<QTVI_log bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001070\">NaN</QTVI_log>\n");
							fwprintf( fp_output,L"\t\t<QT_Dispersion bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001078\">NaN</QT_Dispersion>\n");
						} else {
							fwprintf( fp_output,L"\t\t<QTCorrected_Bazett bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000000701\">%.5lf</QTCorrected_Bazett>\n", (double)meanQT/sqrt(meanRR/1000));
							fwprintf( fp_output,L"\t\t<QTVI_log bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001070\">%.5lf</QTVI_log>\n", log((varQT/(meanQT*meanQT))/(varHR/(meanHR*meanHR))));
							fwprintf( fp_output,L"\t\t<QT_Dispersion bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001078\">%.5lf</QT_Dispersion>\n", maxQT-minQT);
						}
						fwprintf( fp_output,L"\t\t<RRIntervalResults>\n");
						fwprintf( fp_output,L"\t\t\t<RRIntervalCount bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001086\">%d</RRIntervalCount>\n", (int)rrs.size() );
						if( rrs.size() < 1 ) {
							fwprintf( fp_output,L"\t\t\t<RRMean bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001086\">NaN</RRMean>\n");
							fwprintf( fp_output,L"\t\t\t<RRMin bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001091\">NaN</RRMin>\n");
							fwprintf( fp_output,L"\t\t\t<RRMax bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001090\">NaN</RRMax>\n");
							fwprintf( fp_output,L"\t\t\t<RRVariance bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001094\">NaN</RRVariance>\n");
							fwprintf( fp_output,L"\t\t\t<RRStandardDeviation bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001096\">NaN</RRStandardDeviation>\n");
						} else {
							fwprintf( fp_output,L"\t\t\t<RRMean bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001086\">%.5f</RRMean>\n", meanRR);
							fwprintf( fp_output,L"\t\t\t<RRMin bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001091\">%.5f</RRMin>\n", minRR );
							fwprintf( fp_output,L"\t\t\t<RRMax bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001090\">%.5f</RRMax>\n", maxRR );
							fwprintf( fp_output,L"\t\t\t<RRVariance bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001094\">%.5f</RRVariance>\n", varRR );
							fwprintf( fp_output,L"\t\t\t<RRStandardDeviation bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001096\">%.5f</RRStandardDeviation>\n", sdRR );
						}
						fwprintf( fp_output,L"\t\t</RRIntervalResults>\n");
						fwprintf( fp_output,L"\t\t<QTIntervalResults>\n");
						fwprintf( fp_output,L"\t\t\t<QTIntervalCount bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001085\">%d</QTIntervalCount>\n", (int)qts.size() );
						if( qts.size() < 1 ) {
							fwprintf( fp_output,L"\t\t\t<QTMean bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001089\">NaN</QTMean>\n");
							fwprintf( fp_output,L"\t\t\t<QTMin bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001093\">NaN</QTMin>\n");
							fwprintf( fp_output,L"\t\t\t<QTMax bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001092\">NaN</QTMax>\n");
							fwprintf( fp_output,L"\t\t\t<QTVariance bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001095\">NaN</QTVariance>\n");
							fwprintf( fp_output,L"\t\t\t<QTStandardDeviation bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001097\">NaN</QTStandardDeviation>\n");
						} else {
							fwprintf( fp_output,L"\t\t\t<QTMean bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001089\">%.5f</QTMean>\n", meanQT*1000 );
							fwprintf( fp_output,L"\t\t\t<QTMin bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001093\">%.5f</QTMin>\n", minQT*1000 );
							fwprintf( fp_output,L"\t\t\t<QTMax bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001092\">%.5f</QTMax>\n", maxQT*1000 );
							fwprintf( fp_output,L"\t\t\t<QTVariance bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001095\">%.5f</QTVariance>\n", varQT*1000 );
							fwprintf( fp_output,L"\t\t\t<QTStandardDeviation bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001097\">%.5f</QTStandardDeviation>\n", sdQT*1000 );
						}
						fwprintf( fp_output,L"\t\t</QTIntervalResults>\n");
						fwprintf( fp_output,L"\t\t<HRIntervalResults>\n");
						fwprintf( fp_output,L"\t\t\t<HRIntervalCount bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001103\">%d</HRIntervalCount>\n", (int)rrs.size() );
						if( rrs.size() < 1 ) {
							fwprintf( fp_output,L"\t\t\t<HRMean bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001099\">NaN</HRMean>\n");
							fwprintf( fp_output,L"\t\t\t<HRMin bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001100\">NaN</HRMin>\n");
							fwprintf( fp_output,L"\t\t\t<HRMax bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001098\">NaN</HRMax>\n");
							fwprintf( fp_output,L"\t\t\t<HRVariance bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001101\">NaN</HRVariance>\n");
							fwprintf( fp_output,L"\t\t\t<HRStandardDeviation bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001102\">NaN</HRStandardDeviation>\n");
						} else {
							fwprintf( fp_output,L"\t\t\t<HRMean bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001099\">%.5f</HRMean>\n", meanRR );
							fwprintf( fp_output,L"\t\t\t<HRMin bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001100\">%.5f</HRMin>\n", minRR );
							fwprintf( fp_output,L"\t\t\t<HRMax bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001098\">%.5f</HRMax>\n", maxRR );
							fwprintf( fp_output,L"\t\t\t<HRVariance bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001101\">%.5f</HRVariance>\n", varRR );
							fwprintf( fp_output,L"\t\t\t<HRStandardDeviation bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001102\">%.5f</HRStandardDeviation>\n", sdRR );
						}
						fwprintf( fp_output,L"\t\t</HRIntervalResults>\n");
						fwprintf( fp_output,L"\t</LeadResults>\n");
					} // if( qrsAnn )
					else {
						wprintf( L"could not get QRS complexes. make sure you have got \"filters\" directory in the ecg application dir." );
						exit(1);
					}
				} // for( leadNumber...
				fwprintf( fp_output,L"</autoQRSResults>\n");
				fclose( fp_output );
				if(verbose) {
					printf("chesnokov output finished, %d total annotations detected, %d leads.\n\n", annTotalCount, signal.GetLeadsNum()  );
				}
				printf("success");
			}else{  // must be greater than 1 second duration
				printf("Could not analyze %s, duration is less than 1 second (%d:%d:%d).", physionetdatfile, signal.GetH(),signal.GetM(), signal.GetS());
			}
        } else {
            wprintf( L"failed to read %s file.", physionetdatfile );
            exit(1);
        }

    }

	return 0;
}

void help()
{
    printf( "\n|---------------------------------------------------------------------------------------------|\n"
    		"| Chesnokov's AutoQRS algorithm, Linux version 1.04, May 20, 2015. created by Michael Shipway |\n"
    		"|---------------------------------------------------------------------------------------------|\n"
    		"usage: \n"
//    		"chesnokov <filters dir> <physionetfile.dat> <output filepath> <leadnum> <verbose>\n"
    		"chesnokov <filters dir> <physionetfile.dat> <output filepath> [verbose]\n"
    		"1) filters dir - directory containing the filter files\n"
    		"2) physionetfile.dat - full path/filename of the .hea file of the record,\n"
    		"                    (May also be full path/recordName) \n"
//    		"   or the directory in the environment variable WFDB,\n"
//    		"   or the variable DEFWFDB (defined in /opt/wfdb/include/wfdb/wfdblib.h)\n"
//    		"3) leadnum - which lead to analyze, zero to analyze all available.\n"
    		"3) output filepath - full path/filename of the XML output file.\n"
    		"4) [verbose] - optional, if there is any 4th parameter, progress messages will be printed to console.\n"
    		"Do not forget the filters directory must be present.\n"
    		);
}

