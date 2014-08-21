// ecg.cpp : Defines the entry point for the console application.
// Modified by Stephen Granite (sgranite@jhu.edu) to extend functionality
// version 1.0:  Computes corrected QT along with QT and RR statistical functions for all leads in the input set
// version 1.1:  Stores computed values in XML format (autoqrs.xsd schema) to allow for query and cross-comparison
// version 1.2:  Added filters directory path to command line parameters and added computation of QTVI_log to output
// version 1.3:  Revised code to compute HR and RR statistics properly, as those were reversed in the getRRseq function
// version 1.4:  Revised code to compute minimums and maximums for QT, HR and RR, and QT Dispersion
// version 1.5:  Revised code to add bioportal URLs to tag the values
// version 1.6:  Revised code to replace all sr values with 1000, as it is that value for Chesnokov computations, rather than the sampling rate

#include "stdafx.h"
#include "lib\lib.h"

void help();

int _tmain(int argc, _TCHAR* argv[])
{
    //add your code here
    if( argc < 3 ) {
        help();
    }
    else
    {
        int leadNumber = 0;
        if( argc == 4+1 ) {
            leadNumber = _wtoi( argv[4] ) - 1;
            if( leadNumber < 0 ) leadNumber = 0;
        }
		//_wsetlocale(LC_ALL, L"English");
        class Signal signal;
        if( signal.ReadFile( argv[2] ) ) {

			wchar_t *fltdir = argv[1];
			FILE *fp_output = _wfopen( argv[3], L"wt" );
			//FILE *fp = _wfopen( L"annotations.txt", L"wt" );
            int size = signal.GetLength();
            long double sr = signal.GetSR();
            int h,m,s,ms;
            int msec = int( ((long double)size/sr) * 1000.0 );
            signal.mSecToTime(msec,h,m,s,ms);
//		    char leads[18][6] =  {"I","II","III","aVR","aVL","aVF","v1","v2","v3","v4","v5","v6","MLI","MLII","MLIII","vX","vY","vZ"};
		    wchar_t leads[18][6] =  {L"I",L"II",L"III",L"aVR",L"aVL",L"aVF",L"v1",L"v2",L"v3",L"v4",L"v5",L"v6",L"MLI",L"MLII",L"MLIII",L"vX",L"vY",L"vZ"};
			//wchar_t lead[6] = L" ";

			fwprintf( fp_output,L"<?xml version = \"1.0\" encoding = \"UTF-8\"?>\n");
			fwprintf( fp_output,L"<autoQRSResults xmlns=\"http://www.cvrg.org/1/AutoQRSService\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.cvrg.org/1/AutoQRSService\">\n");
			fwprintf( fp_output,L"\t<jobIdentifier/>\n");
			fwprintf( fp_output,L"\t<FileAnalyzed>%s</FileAnalyzed>\n", argv[2]);
            fwprintf( fp_output,L"\t<LeadCount>%d</LeadCount>\n", signal.GetLeadsNum() );
            fwprintf( fp_output,L"\t<SR>%.2lf</SR>\n", sr );
            fwprintf( fp_output,L"\t<UmV>%d</UmV>\n", signal.GetUmV() );
            fwprintf( fp_output,L"\t<Length>%02d:%02d:%02d.%03d</Length>\n", h,m,s,ms );

			for( leadNumber = 0; leadNumber < signal.GetLeadsNum(); leadNumber++) {
				fwprintf( fp_output,L"\t<LeadResults id=\"%d\">\n", leadNumber + 1 );
				signal.GetData( leadNumber );
				wprintf( L"\n\tlead: %ls\n", leads[leadNumber] );
	            wprintf( L"\tsamplingRate: %.2lf\n", sr );
	            wprintf( L"\tbits: %d\n", signal.GetBits() );
				wprintf( L"\tumv: %d\n", signal.GetUmV() );
				wprintf( L"\tsize: %d\n", size );
				long double* data = signal.GetData();
	            
				//annotation
				class Annotation ann;  //default annotation params

				//or add your custom ECG params to annotation class from lib.h
				// ANNHDR hdr;
				//  hdr.minbpm = 30;
				//  etc...
				// class Annotation ann( &hdr );


				wprintf( L"\tgetting QRS complexes... \n" );
//				int** qrsAnn = ann.GetQRS( data, size, sr, L"filters" );       //get QRS complexes
				int** qrsAnn = ann.GetQRS( data, size, sr, fltdir );       //get QRS complexes
				//int** qrsAnn = ann.GetQRS(data,size,sr,fltdir,qNOISE);    //get QRS complexes if signal is quite noisy

				if( qrsAnn ) {
					fwprintf( fp_output,L"\t\t<Lead>%s</Lead>\n", leads[leadNumber] );
					fwprintf( fp_output,L"\t\t<TotalBeatCount>%d</TotalBeatCount>\n", ann.GetQRSnum() );
					ann.GetECT( qrsAnn, ann.GetQRSnum(), sr );      //label Ectopic beats
					fwprintf( fp_output,L"\t\t<EctopicBeatCount>%d</EctopicBeatCount>\n", ann.GetECTnum() );
					wprintf( L"\tgetting P, T waves... \n" );
					int annNum = 0;
					int** ANN = ann.GetPTU( data, size, sr, fltdir, qrsAnn, ann.GetQRSnum() );   //find P,T waves
					if(ANN) {
						annNum = ann.GetANNnum();
						wprintf( L"\t%d total annotations detected.\n", annNum );
					}
					else {
						ANN = qrsAnn;
						annNum = 2 * ann.GetQRSnum();
						wprintf( L"\tfailed.\n" );
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
					for( int i = 0; i < annNum; i++ ) {
						int smpl = ANN[i][0];
						int type = ANN[i][1];

						msec = int( ((long double)smpl/sr) * 1000.0 );
						signal.mSecToTime(msec,h,m,s,ms);
						//fwprintf( fp,L"  %02d:%02d:%02d.%03d   %s\n", h,m,s,ms, anncodes[type] );
					}

					//saving RR seq
					vector<long double> rrs;
					vector<int> rrsPos;
					vector<long double> qts;
					vector<int> qtsPos;
					long double meanRR = 0.0, meanQT = 0.0, QTc = 0.0;
					long double minQT = 10000000000.0, maxQT = 0.0, minRR = 10000000000.0, maxRR = 0.0;
					long double minHR = 10000000000.0, maxHR = 0.0;
					long double varRR = 0.0, sdRR = 0.0, varQT = 0.0, sdQT = 0.0;
					long double meanHR = 0.0, sdHR = 0.0, varHR = 0.0;
					//double dNaN = std::numeric_limits<double>::quiet_NaN();

					if( ann.GetRRseq( ANN, annNum, sr, &rrs, &rrsPos) ) {
						//FILE *fp = _wfopen( L"rrs.txt", L"wt" );
						for( int i = 0; i < (int)rrs.size(); i++ ) {
							//fwprintf( fp, L"%lf\n", rrs[i] );       
							meanHR = meanHR + rrs[i];
							if( rrs[i] > maxHR) maxHR = rrs[i];
							if( rrs[i] < minHR) minHR = rrs[i];
							//rrs[i] = 60/((r2-r1)/sr) so RR = sr * (60/rrs[i])
							//meanRR = meanRR + (sr * (60/rrs[i]));
							meanRR = meanRR + 60000/rrs[i];
						}
//						maxRR = sr * (60/minHR);
//						minRR = sr * (60/maxHR);
						maxRR = 60000/minHR;
						minRR = 60000/maxHR;
						//fclose( fp );
					}
					ann.GetRRseq( ANN, annNum, sr, &rrs, &rrsPos);
					meanRR = meanRR / rrs.size();
					meanHR = meanHR / rrs.size();
					for( int i = 0; i < (int)rrs.size(); i++ ) {
						varHR = varHR + (rrs[i] - meanHR) * (rrs[i] - meanHR);
//						varRR = varRR + (((sr * (60/rrs[i])) - meanRR) * ((sr * (60/rrs[i])) - meanRR));
						varRR = varRR + ((60000/rrs[i] - meanRR) * (60000/rrs[i] - meanRR));
					}
					varRR = varRR / (rrs.size() - 1);
					varHR = varHR / (rrs.size() - 1);
					sdRR = sqrt(varRR);
					sdHR = sqrt(varHR);

					if( ann.GetQTseq( ANN, annNum, sr, &qts, &qtsPos) ) {
						//FILE *fp = _wfopen( L"qts.txt", L"wt" );
						for( int i = 0; i < (int)qts.size(); i++ ) {
							//fwprintf( fp, L"%lf\n", qts[i] );                    
							meanQT = meanQT + qts[i];
							if( qts[i] > maxQT) maxQT = qts[i];
							if( qts[i] < minQT) minQT = qts[i];

						}
//						maxQT = maxQT * sr;
//						minQT = minQT * sr;
						maxQT = maxQT * 1000;
						minQT = minQT * 1000;
						//fclose( fp );
					}
					ann.GetQTseq( ANN, annNum, sr, &qts, &qtsPos);
//					meanQT = meanQT * sr / qts.size();
					meanQT = meanQT * 1000 / qts.size();
					for( int i = 0; i < (int)qts.size(); i++ ) 
//						varQT = varQT + (qts[i] * sr - meanQT) * (qts[i] * sr - meanQT);
						varQT = varQT + (qts[i] * 1000 - meanQT) * (qts[i] * 1000 - meanQT);
					varQT = varQT / (qts.size() - 1);
					sdQT = sqrt(varQT);
					std::string bioportalURLStart = "http://purl.bioontology.org/ontology/ECGT/";

					if( (qts.size() < 1) || (rrs.size() < 1) ) {
						fwprintf( fp_output,L"\t\t<QTCorrected_Bazett bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000000701\">NaN</QTCorrected_Bazett>\n");
						fwprintf( fp_output,L"\t\t<QTVI_log bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001070\">NaN</QTVI_log>\n");
						fwprintf( fp_output,L"\t\t<QT_Dispersion>NaN</QT_Dispersion>\n");
					} else {
//						fwprintf( fp_output,L"\t\t<QTCorrected_Bazett bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000000701\">%.5lf</QTCorrected_Bazett>\n", meanQT/sqrt(meanRR/sr) );
						fwprintf( fp_output,L"\t\t<QTCorrected_Bazett bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000000701\">%.5lf</QTCorrected_Bazett>\n", meanQT/sqrt(meanRR/1000) );
						fwprintf( fp_output,L"\t\t<QTVI_log bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001070\">%.5lf</QTVI_log>\n", log((varQT/(meanQT*meanQT))/(varHR/(meanHR*meanHR))) );
						fwprintf( fp_output,L"\t\t<QT_Dispersion>%.5lf</QT_Dispersion>\n", maxQT-minQT );
					}
					fwprintf( fp_output,L"\t\t<RRIntervalResults>\n");
					fwprintf( fp_output,L"\t\t\t<RRIntervalCount>%d</RRIntervalCount>\n", (int)rrs.size() );
					if( rrs.size() < 1 ) {
						fwprintf( fp_output,L"\t\t\t<RRMean bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000000042\">NaN</RRMean>\n");
						fwprintf( fp_output,L"\t\t\t<RRMin>NaN</RRMin>\n");
						fwprintf( fp_output,L"\t\t\t<RRMax>NaN</RRMax>\n");
						fwprintf( fp_output,L"\t\t\t<RRVariance>NaN</RRVariance>\n");
						fwprintf( fp_output,L"\t\t\t<RRStandardDeviation>NaN</RRStandardDeviation>\n");
					} else {
						fwprintf( fp_output,L"\t\t\t<RRMean bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000000042\">%.5lf</RRMean>\n", meanRR );
						fwprintf( fp_output,L"\t\t\t<RRMin>%.5lf</RRMin>\n", minRR );
						fwprintf( fp_output,L"\t\t\t<RRMax>%.5lf</RRMax>\n", maxRR );
						fwprintf( fp_output,L"\t\t\t<RRVariance>%.5lf</RRVariance>\n", varRR );
						fwprintf( fp_output,L"\t\t\t<RRStandardDeviation>%.5lf</RRStandardDeviation>\n", sdRR );
					}
					fwprintf( fp_output,L"\t\t</RRIntervalResults>\n");
					fwprintf( fp_output,L"\t\t<HRIntervalResults>\n");
					fwprintf( fp_output,L"\t\t\t<HRIntervalCount>%d</HRIntervalCount>\n", (int)rrs.size() );
					if( rrs.size() < 1 ) {
						fwprintf( fp_output,L"\t\t\t<HRMean bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000000833\">NaN</HRMean>\n");
						fwprintf( fp_output,L"\t\t\t<HRMin>NaN</HRMin>\n");
						fwprintf( fp_output,L"\t\t\t<HRMax>NaN</HRMax>\n");
						fwprintf( fp_output,L"\t\t\t<HRVariance>NaN</HRVariance>\n");
						fwprintf( fp_output,L"\t\t\t<HRStandardDeviation>NaN</HRStandardDeviation>\n");
					} else {
						fwprintf( fp_output,L"\t\t\t<HRMean bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000000833\">%.5lf</HRMean>\n", meanHR );
						fwprintf( fp_output,L"\t\t\t<HRMin>%.5lf</HRMin>\n", minHR );
						fwprintf( fp_output,L"\t\t\t<HRMax>%.5lf</HRMax>\n", maxHR );
						fwprintf( fp_output,L"\t\t\t<HRVariance>%.5lf</HRVariance>\n", varHR );
						fwprintf( fp_output,L"\t\t\t<HRStandardDeviation>%.5lf</HRStandardDeviation>\n", sdHR );
					}
					fwprintf( fp_output,L"\t\t</HRIntervalResults>\n");
					fwprintf( fp_output,L"\t\t<QTIntervalResults>\n");
					fwprintf( fp_output,L"\t\t\t<QTIntervalCount>%d</QTIntervalCount>\n", (int)qts.size() );
					if( qts.size() < 1 ) {
						fwprintf( fp_output,L"\t\t\t<QTMean bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000000682\">NaN</QTMean>\n");
						fwprintf( fp_output,L"\t\t\t<QTMin>NaN</QTMin>\n");
						fwprintf( fp_output,L"\t\t\t<QTMax>NaN</QTMax>\n");
						fwprintf( fp_output,L"\t\t\t<QTVariance>NaN</QTVariance>\n");
						fwprintf( fp_output,L"\t\t\t<QTStandardDeviation>NaN</QTStandardDeviation>\n");
					} else {
						fwprintf( fp_output,L"\t\t\t<QTMean bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000000682\">%.5lf</QTMean>\n", meanQT );
						fwprintf( fp_output,L"\t\t\t<QTMin>%.5lf</QTMin>\n", minQT );
						fwprintf( fp_output,L"\t\t\t<QTMax>%.5lf</QTMax>\n", maxQT );
						fwprintf( fp_output,L"\t\t\t<QTVariance>%.5lf</QTVariance>\n", varQT );
						fwprintf( fp_output,L"\t\t\t<QTStandardDeviation>%.5lf</QTStandardDeviation>\n", sdQT );
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
			//fclose( fp );
        } else {
            wprintf( L"failed to read %s file", argv[2] );
            exit(1);
        }

    }

	return 0;
}

void help()
{
    wprintf( L"usage: ecg.exe <filters dir> physionetfile.dat <output filepath> leadnum\n\
       do not forget about \\filters dir to be present." );
}

