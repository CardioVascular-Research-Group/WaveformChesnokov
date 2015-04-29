// ecg.cpp : Defines the entry point for the console application.
// Modified by Stephen Granite (sgranite@jhu.edu) to extend functionality
// version 1.0:  Computes corrected QT along with QT and RR statistical functions for all leads in the input set
// version 1.1:  Stores computed values in XML format (autoqrs.xsd schema) to allow for query and cross-comparison
// version 1.2:  Added filters directory path to command line parameters and added computation of QTVI_log to output
// version 1.3:  Revised code to compute HR and RR statistics properly
// version 1.4:  Revised code to compute minimums and maximums for QT, HR and RR, and QT Dispersion

#include "lib/stdafx.h"
#include "lib/lib.h"


//#include </opt/wfdb/include/wfdb/wfdb.h>


void help();
int runChesnokov(char *fltdir, char *physionetdatfile, char *outputfilepath, int leadNumber);
class Signal signal;

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
	printf("ecg main()\n argv[0]: %s  \n argv[1]: %s  \n argv[2]: %s  \n argv[3]: %s\n", argv[0], argv[1], argv[2], argv[3]);
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
		signal.PopulateSignal(datFileName);
//        signal.PopulateSignal(argv[2]);
        runChesnokov(argv[1],argv[2],argv[3],leadNumber); //<filters dir> , physionetfile.dat, <output filepath>, leadnum
    }
}

int runChesnokov(char *fltdir, char *physionetdatfile, char *outputfilepath, int leadNumber){
	{
        if( signal.GetLeadsNum() != 0){ // pretending that the WFDB input file was read... // !signal.ReadFile( physionetdatfile )

        	cout << "** Open outputfile: " << outputfilepath <<" **\n";
			FILE *fp_output = fopen( outputfilepath, "wt" );
            int size = signal.GetLength();
		    cout << "** GetSR() **\n";

            double sr = signal.GetSR();
            int h,m,s,ms;
            int msec = int( ((long double)size/sr) * 1000.0 );
            signal.mSecToTime(msec,h,m,s,ms);
//		    char leads[18][6] =  {"I","II","III","aVR","aVL","aVF","v1","v2","v3","v4","v5","v6","MLI","MLII","MLIII","vX","vY","vZ"};
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
				printf("** Get signal data for lead:%d **\n", leadNumber );
				signal.GetData( leadNumber );
				fwprintf( fp_output,L"\t<LeadResults id=\"%d\">\n", leadNumber + 1 );

				printf("\tlead: %ls\n", leads[leadNumber] );
	            printf("\tsamplingRate: %.2f\n", sr );
	            printf("\tbits: %d\n", signal.GetBits() );
				printf("\tumv: %d\n", signal.GetUmV() );
				printf("\tsize: %d\n", size );
				printf("\tduration: %02d:%02d:%02d.%03d\n", h,m,s,ms );
				long double* data = signal.GetData();
	            
				//annotation
				class Annotation ann;  //default annotation params

				//or add your custom ECG params to annotation class from lib.h
				// ANNHDR hdr;
				//  hdr.minbpm = 30;
				//  etc...
				// class Annotation ann( &hdr );


				wprintf( L"\tgetting QRS complexes... \n" );
				int** qrsAnn = ann.GetQRS( data, size, sr, fltdir );       //get QRS complexes

				if( qrsAnn ) {
					fwprintf( fp_output,L"\t\t<Lead>%ls</Lead>\n", leads[leadNumber] );
					fwprintf( fp_output,L"\t\t<TotalBeatCount>%d</TotalBeatCount>\n", ann.GetQRSnum() );
					ann.GetECT( qrsAnn, ann.GetQRSnum(), sr );      //label Ectopic beats
					fwprintf( fp_output,L"\t\t<EctopicBeatCount>%d</EctopicBeatCount>\n", ann.GetECTnum() );
//					cout << "\tgetting P, T waves... \n";
					printf("\tgetting P, T waves... \n" );
					int annNum = 0;
					int** ANN = ann.GetPTU( data, size, sr, fltdir, qrsAnn, ann.GetQRSnum() );   //find P,T waves
					if(ANN) {
						annNum = ann.GetANNnum();
						printf("\t%d total annotations detected.\n\n", annNum );
					}
					else {
						ANN = qrsAnn;
						annNum = 2 * ann.GetQRSnum();
						printf("\tfailed.\n" );
					}

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
				//	ann.GetRRseq( ANN, annNum, sr, &rrs, &rrsPos);
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
						fwprintf( fp_output,L"\t\t<QTCorrected_Bazett>NaN</QTCorrected_Bazett>\n");
						fwprintf( fp_output,L"\t\t<QTVI_log>NaN</QTVI_log>\n");
						fwprintf( fp_output,L"\t\t<QT_Dispersion>NaN</QT_Dispersion>\n");
					} else {
						fwprintf( fp_output,L"\t\t<QTCorrected_Bazett>%.5f</QTCorrected_Bazett>\n", (double)meanQT/sqrt(meanRR) );
						fwprintf( fp_output,L"\t\t<QTVI_log>%.5f</QTVI_log>\n", log((varQT/(meanQT*meanQT))/(varHR/(meanHR*meanHR))) );
						fwprintf( fp_output,L"\t\t<QT_Dispersion>%.5f</QT_Dispersion>\n", maxQT-minQT );
					}
					fwprintf( fp_output,L"\t\t<RRIntervalResults>\n");
					fwprintf( fp_output,L"\t\t\t<RRIntervalCount>%d</RRIntervalCount>\n", (int)rrs.size() );
					if( rrs.size() < 1 ) {
						fwprintf( fp_output,L"\t\t\t<RRMean>NaN</RRMean>\n");
						fwprintf( fp_output,L"\t\t\t<RRMin>NaN</RRMin>\n");
						fwprintf( fp_output,L"\t\t\t<RRMax>NaN</RRMax>\n");
						fwprintf( fp_output,L"\t\t\t<RRVariance>NaN</RRVariance>\n");
						fwprintf( fp_output,L"\t\t\t<RRStandardDeviation>NaN</RRStandardDeviation>\n");
					} else {
						fwprintf( fp_output,L"\t\t\t<RRMean>%.5f</RRMean>\n", meanRR );
						fwprintf( fp_output,L"\t\t\t<RRMin>%.5f</RRMin>\n", minRR );
						fwprintf( fp_output,L"\t\t\t<RRMax>%.5f</RRMax>\n", maxRR );
						fwprintf( fp_output,L"\t\t\t<RRVariance>%.5f</RRVariance>\n", varRR );
						fwprintf( fp_output,L"\t\t\t<RRStandardDeviation>%.5f</RRStandardDeviation>\n", sdRR );
					}
					fwprintf( fp_output,L"\t\t</RRIntervalResults>\n");
					fwprintf( fp_output,L"\t\t<QTIntervalResults>\n");
					fwprintf( fp_output,L"\t\t\t<QTIntervalCount>%d</QTIntervalCount>\n", (int)qts.size() );
					if( qts.size() < 1 ) {
						fwprintf( fp_output,L"\t\t\t<QTMean>NaN</QTMean>\n");
						fwprintf( fp_output,L"\t\t\t<QTMin>NaN</QTMin>\n");
						fwprintf( fp_output,L"\t\t\t<QTMax>NaN</QTMax>\n");
						fwprintf( fp_output,L"\t\t\t<QTVariance>NaN</QTVariance>\n");
						fwprintf( fp_output,L"\t\t\t<QTStandardDeviation>NaN</QTStandardDeviation>\n");
					} else {
						fwprintf( fp_output,L"\t\t\t<QTMean>%.5f</QTMean>\n", meanQT );
						fwprintf( fp_output,L"\t\t\t<QTMin>%.5f</QTMin>\n", minQT );
						fwprintf( fp_output,L"\t\t\t<QTMax>%.5f</QTMax>\n", maxQT );
						fwprintf( fp_output,L"\t\t\t<QTVariance>%.5f</QTVariance>\n", varQT );
						fwprintf( fp_output,L"\t\t\t<QTStandardDeviation>%.5f</QTStandardDeviation>\n", sdQT );
					}
					fwprintf( fp_output,L"\t\t</QTIntervalResults>\n");
					fwprintf( fp_output,L"\t</LeadResults>\n");

				}
				else {
					wprintf( L"could not get QRS complexes. make sure you have got \"filters\" directory in the ecg application dir." );
					exit(1);
				}
			}
			fwprintf( fp_output,L"</autoQRSResults>\n");
			fclose( fp_output );
			cout << "output finished.\n";
        } else {
            wprintf( L"failed to read %s file.", physionetdatfile );
            exit(1);
        }

    }

	return 0;
}

void help()
{
    printf( "\n-------------------------\n"
    		"usage: ecg.exe <filters dir> physionetfile.dat <output filepath> leadnum\n"
    		"filters dir - directory containing the filter files\n"
    		"physionetfile.dat - full path/filename of the .hea file of the record,\n"
    		"   (May also be full path/recordName, \n"
//    		"   or the directory in the environment variable WFDB,\n"
//    		"   or the variable DEFWFDB (defined in /opt/wfdb/include/wfdb/wfdblib.h)\n"
    		"leadnum - which lead to analyze, zero to analyze all available.\n"
    		"Do not forget the filters directory must be present.\n" );
}

