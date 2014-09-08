/*

ecg.cpp - ECG annotation console app based lib.cpp ECG library.
Copyright (C) 2007 YURIY V. CHESNOKOV

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


You may contact the author by e-mail (chesnokov.yuriy@gmail.com or chesnokov_yuriy@mail.ru)
or postal mail (Unilever Centre for Molecular Sciences Informatics,
University Chemical Laboratory, Cambridge University, Lensfield Road, Cambridge, CB2 1EW, UK)

*/
// ecg.cpp : Defines the entry point for the console application.
// Modified by Stephen Granite (sgranite@jhu.edu) to extend functionality
// version 1.0:  Computes corrected QT along with QT and RR statistical functions for all leads in the input set
// version 1.1:  Stores computed values in XML format (autoqrs.xsd schema) to allow for query and cross-comparison
// version 1.2:  Added filters directory path to command line parameters and added computation of QTVI_log to output
// version 1.3:  Revised code to compute HR and RR statistics properly, as those were reversed in the getRRseq function
// version 1.4:  Revised code to compute minimums and maximums for QT, HR and RR, and QT Dispersion
// version 1.5:  Revised code to add bioportal URLs to tag the values
// version 1.6:  Revised code to replace all sr values with 1000, as it is that value for Chesnokov computations, rather than the sampling rate
// version 1.65: Revised code for clarity, converting sample-rate to msec or Heart Rate for all values before mean, variance and standard deviation are calculated.
//               Added documentation and regrouped some variable declarations.
#include "stdafx.h"
#include "lib\lib.h"

wchar_t params[_MAX_PATH] = L"params";

void tic();
void toc();
void saveXMLheader(FILE *fp_output, _TCHAR* dataFile, int LeadsNum, double sr, int UmV, int h, int m, int s, int ms);
void saveXMLleadResultsHeader(FILE *fp_output, int leadID, wchar_t *lead, int totalBeatCount, int ECTnum, 
						double qtsSize, double rrsSize, 
						double QTCorrected_Bazett, double QTVI_log, double QT_Dispersion );
void saveXMLleadResultsRRinterval(FILE *fp_output,  double rrsSize,
			double meanRR,  double minRR, double maxRR, double varRR, double sdRR);
void saveXMLleadResultsHRinterval(FILE *fp_output,  double rrsSize,
			double meanHR,  double minHR, double maxHR, double varHR, double sdHR);
void saveXMLleadResultsQTinterval(FILE *fp_output, double qtsSize, 
			double meanQT, double varQT, double sdQT, double maxQT,  double minQT);
//void saveXMLWaveMeanResults(FILE *fp_output, EcgAnnotation &ann);
void saveXMLleadResultsFooter(FILE *fp_output);
//void saveXMLleadResults(FILE *fp_output, int leadID, wchar_t *lead, int totalBeatCount, double ECTnum, double qtsSize, double rrsSize, 
//						double meanHR, double varHR, 
//						double meanRR, double varRR, double sdRR, double maxRR, double minRR,
//						double meanQT, double varQT, double sdQT, double maxQT,  double minQT, 
//						int countQ, double qAmp, double qDur,
//						int countR, double rAmp, double rDur, 
//						int countS, double sAmp, double sDur);
void saveXMLfooter(FILE *fp_output);
void help();
int parse_params(class EcgAnnotation &ann);
void change_extension(wchar_t* path, const wchar_t* ext);
float codeVersion = 1.02;

int _tmain(int argc, _TCHAR* argv[])
{
    //add your code here
    if( argc < 3 ) {
        help();
    }
    else
    {
        int leadNumber = 0;
		_TCHAR* dataFile = argv[2];
		_TCHAR* outFile = argv[3];
        if( argc == 4+1 ) {
            leadNumber = _wtoi( argv[4] ) - 1;
            if( leadNumber < 0 ) leadNumber = 0;
        }
		//_wsetlocale(LC_ALL, L"English");
        class Signal signal;
        if( signal.ReadFile( dataFile ) ) {

			wchar_t *fltdir = argv[1];
			FILE *fp_output = _wfopen( outFile, L"wt" );
			//FILE *fp = _wfopen( L"annotations.txt", L"wt" );
            int size = signal.GetLength();
            long double sr = signal.GetSR();
			int leadsNum = signal.GetLeadsNum();
			int umV = signal.GetUmV();
            int h,m,s,ms;
            int msec = int( ((long double)size/sr) * 1000.0 );
//			int ECTnum=0, QRSnum=0;
            signal.mSecToTime(msec,h,m,s,ms);
//		    char leads[18][6] =  {"I","II","III","aVR","aVL","aVF","v1","v2","v3","v4","v5","v6","MLI","MLII","MLIII","vX","vY","vZ"};
		    wchar_t leads[18][6] =  {L"I",L"II",L"III",L"aVR",L"aVL",L"aVF",L"v1",L"v2",L"v3",L"v4",L"v5",L"v6",L"MLI",L"MLII",L"MLIII",L"vX",L"vY",L"vZ"};
			//wchar_t lead[6] = L" ";
			saveXMLheader(fp_output, dataFile, leadsNum, sr, umV, h, m, s, ms);
			/*
			fwprintf( fp_output,L"<?xml version = \"1.0\" encoding = \"UTF-8\"?>\n");
			fwprintf( fp_output,L"<autoQRSResults xmlns=\"http://www.cvrg.org/1/AutoQRSService\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.cvrg.org/1/AutoQRSService\">\n");
			fwprintf( fp_output,L"\t<jobIdentifier/>\n");
			fwprintf( fp_output,L"\t<FileAnalyzed>%s</FileAnalyzed>\n", argv[2]);
            fwprintf( fp_output,L"\t<LeadCount>%d</LeadCount>\n", signal.GetLeadsNum() );
            fwprintf( fp_output,L"\t<SR>%.2lf</SR>\n", sr );
            fwprintf( fp_output,L"\t<UmV>%d</UmV>\n", signal.GetUmV() );
            fwprintf( fp_output,L"\t<Length>%02d:%02d:%02d.%03d</Length>\n", h,m,s,ms );
			*/
			wprintf( L"\nWaveformChesnokov Version  %.2lf\n", codeVersion);

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
	//					fwprintf( fp_output,L"\t\t<Lead>%s</Lead>\n", leads[leadNumber] );
	//					fwprintf( fp_output,L"\t\t<TotalBeatCount>%d</TotalBeatCount>\n", ann.GetQRSnum() );
						int QRSnum = ann.GetQRSnum();
	//					fwprintf( fp_output,L"\t\t<TotalBeatCount>%d</TotalBeatCount>\n", QRSnum );
						ann.GetECT( qrsAnn, ann.GetQRSnum(), sr );      //label Ectopic beats
	//					fwprintf( fp_output,L"\t\t<EctopicBeatCount>%d</EctopicBeatCount>\n", ann.GetECTnum() );
	//					ann.GetECT( qrsAnn, QRSnum, sr );      //label Ectopic beats
						int ECTnum = ann.GetECTnum();

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
						/*
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
						*/
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
						//long double QTc = 0.0;
						long double meanQT = 0.0; // QT interval mean in milliseconds.
						long double meanRR = 0.0; // RR interval mean in milliseconds.
						long double meanHR = 0.0; // Heart Rate  mean interval in milliseconds.
						long double minQT = 10000000000.0, maxQT = 0.0; // QT interval min/max in milliseconds.
						long double minRR = 10000000000.0, maxRR = 0.0; // RR interval min/max in milliseconds.
						long double minHR = 10000000000.0, maxHR = 0.0; // Heart Rate  min/max in beats per minute.
						long double varQT = 0.0, sdQT = 0.0; // QT interval variance and standard deviation.
						long double varRR = 0.0, sdRR = 0.0; // RR interval variance and standard deviation.
						long double sdHR = 0.0, varHR = 0.0; // Heart Rate  variance and standard deviation.

						long double HRsamp, HRbpm; // intermediate conversion value
						long double RRmsec, QTmsec; // intermediate conversion value
						//double dNaN = std::numeric_limits<double>::quiet_NaN();

						if( ann.GetRRseq( ANN, annNum, sr, &rrs, &rrsPos) ) { // gets array of RR intervals in instantanious Beats per Minute?
							//FILE *fp = _wfopen( L"rrs.txt", L"wt" );
							for( int i = 0; i < (int)rrs.size(); i++ ) {
								//fwprintf( fp, L"%lf\n", rrs[i] );    

	//							meanHR = meanHR + rrs[i];
	//							if( rrs[i] > maxHR) maxHR = rrs[i];
	//							if( rrs[i] < minHR) minHR = rrs[i];

								HRsamp = rrs[i];		// rrs[] is measured in "count of samples",
//								HRbpm = (60*sr)/HRsamp;	//   this converts it to Beats per Minute.
								HRbpm = HRsamp;	//   this is already Beats per Minute??
								meanHR = meanHR + HRbpm;
								if( HRbpm > maxHR) maxHR = HRbpm;
								if( HRbpm < minHR) minHR = HRbpm;

								//rrs[i] = 60/((r2-r1)/sr) so RR = sr * (60/rrs[i])
								//meanRR = meanRR + (sr * (60/rrs[i]));
//								RRmsec = (rrs[i]/sr)*1000; // rrs[] is measured in "count of samples", this converts it to milliseconds.
								RRmsec = 60000/rrs[i]; // rrs[] is already Beats per Minute??, this converts it to milliseconds.
								meanRR = meanRR + RRmsec;
								if( RRmsec > maxRR) maxRR = RRmsec;
								if( RRmsec < minRR) minRR = RRmsec;

								if(i<2){
									cout << "rrs[" << i << "] = " << rrs[i] << " : " << RRmsec << " msec, HR=" << HRbpm << "bpm, minRR=" << minRR <<  ", maxRR=" << maxRR << "\n";
								}
								//meanRR = meanRR + (rrs[i]/sr)/1000; // rrs[] is measured in "count of samples", this converts it to milliseconds.
							}
	//						maxRR = sr * (60/minHR);
	//						minRR = sr * (60/maxHR);
	//						maxRR = 60000/minHR;
	//						minRR = 60000/maxHR;
							//fclose( fp );
						}
						ann.GetRRseq( ANN, annNum, sr, &rrs, &rrsPos);
						meanRR = meanRR / rrs.size();
						meanHR = meanHR / rrs.size();
						// cout << "******** minHR = " << minHR <<  ", meanHR = " << meanHR << ", maxHR = " << maxHR << "\n";
						// cout << "******** minRR = " << minRR <<  ", meanRR = " << meanRR << ", maxRR = " << maxRR << "\n";

						for( int i = 0; i < (int)rrs.size(); i++ ) {
							HRsamp = rrs[i]; // rrs[] is measured in "count of samples"
							HRbpm = (60*sr)/HRsamp;
						//	varHR = varHR + (rrs[i] - meanHR) * (rrs[i] - meanHR);
							varHR = varHR + (HRbpm - meanHR) * (HRbpm - meanHR);

							RRmsec = (rrs[i]/sr)*1000; // rrs[] is measured in "count of samples", this converts it to milliseconds.
	//						varRR = varRR + (((sr * (60/rrs[i])) - meanRR) * ((sr * (60/rrs[i])) - meanRR));
	//						varRR = varRR + ((60000/rrs[i] - meanRR) * (60000/rrs[i] - meanRR));
							varRR = varRR + ((RRmsec - meanRR) * (RRmsec - meanRR));
						}
						varRR = varRR / (rrs.size() - 1);
						varHR = varHR / (rrs.size() - 1);
						sdRR = sqrt(varRR);
						sdHR = sqrt(varHR);
						// cout << "******** varHR: " << varHR << ", sdHR: " << sdHR << ", varRR: " << varRR << ", sdRR: " << sdRR << "\n";

						if( ann.GetQTseq( ANN, annNum, sr, &qts, &qtsPos) ) {
							//FILE *fp = _wfopen( L"qts.txt", L"wt" );
							for( int i = 0; i < (int)qts.size(); i++ ) {
								//fwprintf( fp, L"%lf\n", qts[i] );                    
								QTmsec = qts[i]*(1000); // qts[] is measured in Seconds, this converts it to milliseconds.
	//							meanQT = meanQT + qts[i];
	//							if( qts[i] > maxQT) maxQT = qts[i];
	//							if( qts[i] < minQT) minQT = qts[i];
								meanQT = meanQT + QTmsec;
								if( QTmsec > maxQT) maxQT = QTmsec;
								if( QTmsec < minQT) minQT = QTmsec;

								if(i<0){
									// cout << "qts[" << i << "] = " << qts[i] << " : " << QTmsec << " msec, meanQT = " << meanQT << ", min = " << minQT <<  ", max = " << maxQT << "\n";
								}
							}
	//						maxQT = maxQT * sr;
	//						minQT = minQT * sr;
	//						maxQT = maxQT * 1000;
	//						minQT = minQT * 1000;
							//fclose( fp );
						}
						ann.GetQTseq( ANN, annNum, sr, &qts, &qtsPos);
	//					meanQT = meanQT * sr / qts.size();
	//					meanQT = meanQT * 1000 / qts.size();
						meanQT = meanQT / qts.size();
					
						for( int i = 0; i < (int)qts.size(); i++ ){ 
							QTmsec = qts[i]*1000; // qts[] is measured in "count of samples", this converts it to milliseconds.
	//						varQT = varQT + (qts[i] * sr - meanQT) * (qts[i] * sr - meanQT);
	//						varQT = varQT + (qts[i] * 1000 - meanQT) * (qts[i] * 1000 - meanQT);
							varQT = varQT + (QTmsec - meanQT) * (QTmsec - meanQT); 
						}
						varQT = varQT / (qts.size() - 1);
						sdQT = sqrt(varQT);
						std::string bioportalURLStart = "http://purl.bioontology.org/ontology/ECGT/";
						double QTCorrected_Bazett =  meanQT/sqrt(meanRR/1000);  // meanQT/sqrt(meanRR/1000) 
	//					// cout << "******** meanQT = " << meanQT <<  ", meanRR = " << meanRR << ", meanQT/sqrt(meanRR/1000) = " << meanQT/sqrt(meanRR/1000) << ", meanQT/sqrt(60/meanHR) = " << meanQT/sqrt(60/meanHR) << "\n";
						double QTVI_log = log((varQT/(meanQT*meanQT))/(varHR/(meanHR*meanHR))); // log((varQT/(meanQT*meanQT))/(varHR/(meanHR*meanHR)))
						double QT_Dispersion =  maxQT-minQT;

						// cout << "******** minQT = " << minQT <<  ", meanQT = " << meanQT << ", maxQT = " << maxQT << "\n";
						// cout << "******** varQT: " << varQT << ", sdQT: " << sdQT << ", QTCorrected_Bazett: " << QTCorrected_Bazett << ", QTVI_log: " << QTVI_log << ", QT_Dispersion: " << QT_Dispersion << "\n";
						//******** fwprintf the XML result nodes ***********************//
						saveXMLleadResultsHeader(fp_output, (leadNumber+1), leads[leadNumber], QRSnum, ECTnum, 
								qts.size(), rrs.size(), QTCorrected_Bazett, QTVI_log, QT_Dispersion);
						//**************************************
						saveXMLleadResultsRRinterval(fp_output, rrs.size(), meanRR, minRR, maxRR, varRR, sdRR);
						//**************************************
						saveXMLleadResultsHRinterval(fp_output,rrs.size(), meanHR, minHR, maxHR, varHR, sdHR);
						//**************************************
						saveXMLleadResultsQTinterval(fp_output, qts.size(), meanQT, minQT, maxQT, varQT, sdQT);
						//**************************************
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


void saveXMLheader(FILE *fp_output, _TCHAR* dataFile, int LeadsNum, double sr, int UmV, int h, int m, int s, int ms){
	fwprintf( fp_output,L"<?xml version = \"1.0\" encoding = \"UTF-8\"?>\n");
	fwprintf( fp_output,L"<autoQRSResults xmlns=\"http://www.cvrg.org/1/AutoQRSService\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.cvrg.org/1/AutoQRSService\">\n");
	fwprintf( fp_output,L"\t<jobIdentifier/>\n");
	fwprintf( fp_output,L"\t<codeVersion>%.2lf</codeVersion>\n", codeVersion);
	fwprintf( fp_output,L"\t<FileAnalyzed>%s</FileAnalyzed>\n", dataFile);
	fwprintf( fp_output,L"\t<LeadCount>%d</LeadCount>\n", LeadsNum);
	fwprintf( fp_output,L"\t<SR>%.2lf</SR>\n", sr ); // Sample Rate in samples per second (Hz)
	fwprintf( fp_output,L"\t<UmV>%d</UmV>\n", UmV);
	fwprintf( fp_output,L"\t<Length>%02d:%02d:%02d.%03d</Length>\n", h,m,s,ms );
}


void saveXMLleadResultsHeader(FILE *fp_output, int leadID, wchar_t *lead, int totalBeatCount, int ECTnum, 
						double qtsSize, double rrsSize, 
						double QTCorrected_Bazett, double QTVI_log, double QT_Dispersion ){

		fwprintf( fp_output,L"\t\t<Lead>%s</Lead>\n", lead );
		fwprintf( fp_output,L"\t\t<TotalBeatCount>%d</TotalBeatCount>\n", totalBeatCount );
		fwprintf( fp_output,L"\t\t<EctopicBeatCount>%d</EctopicBeatCount>\n", ECTnum);

		if( (qtsSize < 1) || (rrsSize < 1) ) {
			fwprintf( fp_output,L"\t\t<QTCorrected_Bazett bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000000701\">NaN</QTCorrected_Bazett>\n");
			fwprintf( fp_output,L"\t\t<QTVI_log bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001070\">NaN</QTVI_log>\n");
			fwprintf( fp_output,L"\t\t<QT_Dispersion>NaN</QT_Dispersion>\n");
		} else {
			fwprintf( fp_output,L"\t\t<QTCorrected_Bazett bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000000701\">%.5lf</QTCorrected_Bazett>\n", QTCorrected_Bazett);
			fwprintf( fp_output,L"\t\t<QTVI_log bioportalURL=\"http://purl.bioontology.org/ontology/ECGT/ECGTermsv1:ECG_000001070\">%.5lf</QTVI_log>\n", QTVI_log);
			fwprintf( fp_output,L"\t\t<QT_Dispersion>%.5lf</QT_Dispersion>\n", QT_Dispersion);
		}
}


void saveXMLleadResultsRRinterval(FILE *fp_output,  double rrsSize,
			double meanRR,  double minRR, double maxRR, double varRR, double sdRR){

	fwprintf( fp_output,L"\t\t<RRIntervalResults>\n");
		fwprintf( fp_output,L"\t\t\t<RRIntervalCount>%d</RRIntervalCount>\n", (int)rrsSize );
		if( rrsSize < 1 ) {
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
}


void saveXMLleadResultsHRinterval(FILE *fp_output,  double rrsSize,
			double meanHR,  double minHR, double maxHR, double varHR, double sdHR){

	//**************************************
	fwprintf( fp_output,L"\t\t<HRIntervalResults>\n");
		fwprintf( fp_output,L"\t\t\t<HRIntervalCount>%d</HRIntervalCount>\n", (int)rrsSize );
		if( rrsSize < 1 ) {
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
}

void saveXMLleadResultsQTinterval(FILE *fp_output, double qtsSize, 
			double meanQT, double minQT, double maxQT, double varQT, double sdQT){

	fwprintf( fp_output,L"\t\t<QTIntervalResults>\n");
		fwprintf( fp_output,L"\t\t\t<QTIntervalCount>%d</QTIntervalCount>\n", (int)qtsSize );
		if( qtsSize < 1 ) {
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
}
/*
void saveXMLWaveMeanResults(FILE *fp_output, EcgAnnotation &ann){
	int qCount  = ann.GetQcount();
	double qAmp = ann.GetQAmplitudeMean();
	double qDur = ann.GetQDurationMean();
	double qAmpVar = ann.GetQAmplitudeVariance();
	double qDurVar = ann.GetQDurationVariance();


	int rCount  = ann.GetRcount();
	double rAmp = ann.GetRAmplitudeMean();
	double rDur = ann.GetRDurationMean();
	double rAmpVar = ann.GetRAmplitudeVariance();
	double rDurVar = ann.GetRDurationVariance();

	int sCount  = ann.GetScount();
	double sAmp = ann.GetSAmplitudeMean();
	double sDur = ann.GetSDurationMean();
	double sAmpVar = ann.GetSAmplitudeVariance();
	double sDurVar = ann.GetSDurationVariance();

	
//	if(_isnan(meanQT)){ meanQT = 0;}
//	if(_isnan(minQT)){ minQT = 0;}
//	if(_isnan(maxQT)){ maxQT = 0;}
//	if(_isnan(varQT)){ varQT = 0;}
//	if(_isnan(sdQT)){ sdQT = 0;}
	

	fwprintf( fp_output,L"\t\t<WaveMeanResults>\n");
		fwprintf( fp_output,L"\t\t\t<QWaveCount>%i</QWaveCount>\n", qCount );
		if( qCount < 1 ) {
			fwprintf( fp_output,L"\t\t\t<QDurationMean>NaN</QDurationMean>\n");
			fwprintf( fp_output,L"\t\t\t<QDurationVariance>NaN</QDurationVariance>\n");
			fwprintf( fp_output,L"\t\t\t<QDurationStandardDeviation>NaN</QDurationStandardDeviation>\n");

			fwprintf( fp_output,L"\t\t\t<QAmplitudeMean>NaN</QAmplitudeMean>\n");
			fwprintf( fp_output,L"\t\t\t<QAmplitudeVariance>NaN</QAmplitudeVariance>\n");
			fwprintf( fp_output,L"\t\t\t<QAmplitudeStandardDeviation>NaN</QAmplitudeStandardDeviation>\n");
		} else {
			fwprintf( fp_output,L"\t\t\t<QDurationMean>%.4f</QDurationMean>\n", qDur );
			fwprintf( fp_output,L"\t\t\t<QDurationVariance>%.6f</QDurationVariance>\n", qDurVar );
			fwprintf( fp_output,L"\t\t\t<QDurationStandardDeviation>%.6f</QDurationStandardDeviation>\n", sqrt(qDurVar) );

			fwprintf( fp_output,L"\t\t\t<QAmplitudeMean>%.4f</QAmplitudeMean>\n", qAmp );
			fwprintf( fp_output,L"\t\t\t<QAmplitudeVariance>%.6f</QAmplitudeVariance>\n", qAmpVar );
			fwprintf( fp_output,L"\t\t\t<QAmplitudeStandardDeviation>%.6f</QAmplitudeStandardDeviation>\n", sqrt(qAmpVar) );
		}

		fwprintf( fp_output,L"\t\t\t<RWaveCount>%i</RWaveCount>\n", rCount );
		if( rCount < 1 ) {
			fwprintf( fp_output,L"\t\t\t<RDurationMean>NaN</RDurationMean>\n");
			fwprintf( fp_output,L"\t\t\t<RDurationVariance>NaN</RDurationVariance>\n");
			fwprintf( fp_output,L"\t\t\t<RDurationStandardDeviation>NaN</RDurationStandardDeviation>\n");

			fwprintf( fp_output,L"\t\t\t<RAmplitudeMean>NaN</RAmplitudeMean>\n");
			fwprintf( fp_output,L"\t\t\t<RAmplitudeVariance>NaN</RAmplitudeVariance>\n");
			fwprintf( fp_output,L"\t\t\t<RAmplitudeStandardDeviation>NaN</RAmplitudeStandardDeviation>\n");
		} else {
			fwprintf( fp_output,L"\t\t\t<RDurationMean>%.4f</RDurationMean>\n", rDur );
			fwprintf( fp_output,L"\t\t\t<RDurationVariance>%.6f</RDurationVariance>\n", rDurVar );
			fwprintf( fp_output,L"\t\t\t<RDurationStandardDeviation>%.4f</RDurationStandardDeviation>\n", sqrt(rDurVar) );

			fwprintf( fp_output,L"\t\t\t<RAmplitudeMean>%.4f</RAmplitudeMean>\n", rAmp );
			fwprintf( fp_output,L"\t\t\t<RAmplitudeVariance>%.6f</RAmplitudeVariance>\n", rAmpVar );
			fwprintf( fp_output,L"\t\t\t<RAmplitudeStandardDeviation>%.4f</RAmplitudeStandardDeviation>\n", sqrt(rAmpVar) );
		}

		fwprintf( fp_output,L"\t\t\t<SWaveCount>%i</SWaveCount>\n", sCount );
		if( sCount < 1 ) {
			fwprintf( fp_output,L"\t\t\t<SDurationMean>NaN</SDurationMean>\n");
			fwprintf( fp_output,L"\t\t\t<SDurationVariance>NaN</SDurationVariance>\n");
			fwprintf( fp_output,L"\t\t\t<SDurationStandardDeviation>NaN</SDurationStandardDeviation>\n");

			fwprintf( fp_output,L"\t\t\t<SAmplitudeMean>NaN</SAmplitudeMean>\n");
			fwprintf( fp_output,L"\t\t\t<SAmplitudeVariance>NaN</SAmplitudeVariance>\n");
			fwprintf( fp_output,L"\t\t\t<SAmplitudeStandardDeviation>NaN</SAmplitudeStandardDeviation>\n");
		} else {
			fwprintf( fp_output,L"\t\t\t<SDurationMean>%.4f</SDurationMean>\n", sDur );
			fwprintf( fp_output,L"\t\t\t<SDurationVariance>%.6f</SDurationVariance>\n", sDurVar );
			fwprintf( fp_output,L"\t\t\t<SDurationStandardDeviation>%.6f</SDurationStandardDeviation>\n", sqrt(sDurVar) );

			fwprintf( fp_output,L"\t\t\t<SAmplitudeMean>%.4f</SAmplitudeMean>\n", sAmp );
			fwprintf( fp_output,L"\t\t\t<SAmplitudeVariance>%.6lf</SAmplitudeVariance>\n", sAmpVar );
			fwprintf( fp_output,L"\t\t\t<SAmplitudeStandardDeviation>%.6lf</SAmplitudeStandardDeviation>\n", sqrt(sAmpVar) );
		}
	fwprintf( fp_output,L"\t\t</WaveMeanResults>\n");

}
*/
void saveXMLleadResultsFooter(FILE *fp_output){
	fwprintf( fp_output,L"\t</LeadResults>\n");
}


void saveXMLfooter(FILE *fp_output){
	fwprintf( fp_output,L"</autoQRSResults>\n");
}


void help()
{
   // wprintf( L"usage: ecg.exe <filters dir> physionetfile.dat <output filepath> leadnum\n\
       do not forget about \\filters dir to be present." );

	wprintf(L"usage: ecg.exe <filters dir> physionetfile.dat outputFilename.csv [params]\n");
	wprintf(L"   <filters dir>(required) - Do not forget that \\filters dir must be present as a subdirectory of ecg.exe's directory.\n");
	wprintf(L"       physionetfile.dat (required) - name of the data file a Physionet WFDB formatted ECG record.\n");
	wprintf(L"       outputFilename.csv (required) - Path and Name of the CSV output file to save the results into.\n");
	wprintf(L"------------------\n");
}

