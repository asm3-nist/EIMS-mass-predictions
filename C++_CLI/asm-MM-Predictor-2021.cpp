// asm-MM-Predictor-2021.cpp : This file contains the 'main' function. Program execution begins and ends there.
/*Summary: This code will recieve a user indicated input MSP file (query mass spectrum) from file folder,
a user indicated Library (NIST format) from file folder, and compute 3 molecular mass predictions with classifier indicies.
The predictions will be:
(1) Peak Interpretation Method(PIM)
(2) Simple Search Hitlist Method(SS-HM)
(3) iterative Hybrid Search Hitlist  Method(iHS-HM)

Algorithm details are described in manuscript :
Moorthy et.al. "Inferring the nominal molecular mass of an analyte from its electron ionization mass spectrum" (2021)

Date : 2021 - 06 - 23
Version : 0.1
*/
//

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>  
#include <numeric>
#include <fstream>
#include <filesystem>
#include <algorithm>
#include <windows.h>
#include <chrono>
#include <ctime>

using namespace std;

const int DEBUG = 0;
const double minMI_ABUNDANCE = 0.3;
const double epsilon_F = 0.0;
const double beta = 5.0;
const int N = 25;
const int B = 75;
const int mEMF = 700;

class Spectrum
{
    string spectrumFile;
    int numpeaks;
    vector<double> mz;
    vector<double> ab;
public:
    Spectrum(string);
    ~Spectrum();
    string getFN() { return spectrumFile; }
    vector<double> getMZ() { return mz; }
    vector<double> getAB() { return ab; }
    int getNP() { return numpeaks; }
    double PIM_core(double thresholdAB);
    vector<double> PIM(double thresholdAB);
    vector<double> SSHM(string libraryFolder, string library, double thresholdAB);
    vector<double> iHSHM(string libraryFolder, string library,  double thresholdAB);
};

Spectrum::Spectrum(string spectrumFile)
{
    this->spectrumFile = spectrumFile;

    int 	i = 0, j = 0, k = 0; // indices for future use
    int 	numpeaks = 0;

    string line;
    char letter;
    string part_line;
    string empty_string;
    empty_string = '\t';
    int test = 100000;
    int t = 1;

    string s_number;
    double d_number;

    vector<double> mz;
    vector<double> ab;

    ifstream file(spectrumFile);
    if (!file) {
        cout << "Query spectrum not loaded.\n";
    }
    else {
        for (i = 0; !file.eof(); i++) {
            getline(file, line);
            part_line = line.substr(0, 9);
            if (part_line == "Num Peaks") {
                test = i;
            }
            if (i > test) {
                for (k = 0; k < line.length(); k++) {
                    if (line[k] == ' ') {
                        if (t == 1) {
                            d_number = stod(s_number);
                            // cout << s_number;
                            mz.push_back(d_number);
                            s_number.clear();
                        }
                        else {
                            t++;
                        }

                    }
                    else if (line[k] == ';') {
                        d_number = stod(s_number);
                        // cout << ' ' << s_number << endl;
                        ab.push_back(d_number);
                        s_number.clear();
                        t = 0;
                        numpeaks++;
                    }
                    else {
                        s_number += line[k];
                    }
                }
            }
        }
    }

    file.close();


    if (DEBUG == 1) {
        cout << "Iterations: " << i << endl;
        cout << "Num Peaks: " << numpeaks << endl << endl;

        for (i = 0; i < mz.size(); i++) {
            cout << mz[i] << "\t" << ab[i] << endl;
        }
    }

    this->numpeaks = numpeaks;
    this->mz = mz;
    this->ab = ab;

    // return 0;
}

//Spectrum::~Spectrum() { std::cout << spectrumFile <<" spectrum destroyed.\n";}
Spectrum::~Spectrum() { }

double Spectrum::PIM_core(double thresholdAB)
{
    double result=0.0;
    double bp;
    double lowest_ab;
    bp = *max_element(ab.begin(), ab.end());
    lowest_ab = (thresholdAB / 100.0) * bp;

    if (numpeaks < 2) { return(mz[0]);}
    int i = 0;
    double m0 = 0.0;
    double a0 = 0.0;
    double m1 = 0.0;
    double a1 = 0.0;
    double md1 = 0.0;
    double m2 = 0.0;
    double a2 = 0.0;
    double md2 = 0.0;
    double isocalc = 0.0;

    //std::cout << endl;
    for (i = (numpeaks - 1); i >= 2; i--) {
        m0 = mz[i]; a0 = ab[i];
        m1 = mz[(i - 1)]; md1 = m0 - m1; a1 = ab[(i - 1)];

        if ((i - 2) < 1) {
            md2 = 0.0;
        } else {
            m2 = mz[(i - 2)];
            md2 = m0 - m2;
            a2 = ab[(i - 2)];
        }

        if (a0 < lowest_ab) {
            continue;
        }

        if (md1 > 2.0) {
            return(m0); // no peaks within 2 Da
        }

        if ((md1 - 2.0) <= 1e-8) { // possible that compound contains Cl or Br
            if (a1 < (int)(a0 / 2)) {
                return(m0);
            }
        } else if ((md1 - 1.0) <= 1e-8) { // check for c13 isotope peaks

            isocalc = (a1 * 0.011 * m1) / 14.0;
            if ((int)(3 * isocalc) > a0 || (a0 - (int)isocalc) < lowest_ab) {
                if (i == 2) {
                    return(m1);
                }
                continue;
            }

            if (md2 > 2) {
                return(m0);
            } else {
                if (i - 2 >= 1) {
                    if (a2 < (int)(a0 / 2)) {
                        return(m0);
                    }
                }
            }

        }

        //std::cout << i << endl;
    }
   



    //cout << lowest_ab << endl;
    return(mz[numpeaks-1]);
}

vector<double> Spectrum::PIM(double thresholdAB) 
{
    auto t1 = chrono::high_resolution_clock::now();
    cout << "Peak Interpretation Method (PIM): Working ";

    int i,j,k;
    vector<double> results;
    double gamma = PIM_core(thresholdAB);

    results.push_back(gamma); // set results[0] to be the PIM prediction gamma
    
    vector<double> ill_Losses;
    string line;
    ifstream file("data/IllogicalLosses.txt");
    if (!file) {
        cout << "IllogicalLosses.txt not loaded.\n";
    }
    else {
        for (i = 0; !file.eof(); i++) {
            getline(file, line);
            ill_Losses.push_back(stod(line));
        }
    }
    file.close();

    vector<double>::iterator it;
    int index;

    it = find(mz.begin(),mz.end(), gamma);
    index = it - mz.begin();
    double pAb = 0.0;
    pAb = ab[index];

    it = find(ab.begin(), ab.end(), *max_element(ab.begin(),ab.end()) );
    index = it - ab.begin();
    double maxAb = 0.0;
    maxAb = ab[index];

    vector<double> lambda;
    for (i = 0; i < mz.size(); i++) {
        lambda.push_back(gamma - mz[i]);
    }

    vector<int> A;
    for (i = 0; i < lambda.size(); i++) {
        it = find(ill_Losses.begin(), ill_Losses.end(), lambda[i]);
        if (it != ill_Losses.end()) {
            A.push_back(i);
        }
    }

    double r_tau = 0;
    for (i = 0; i < A.size(); i++) {
        r_tau += max((ab[A[i]]-epsilon_F*maxAb),0.0) ;
    }
    r_tau = r_tau / pAb;

    double I = 0.0;
    I = 2 - 2 / (1 + exp(-1.0 * beta * r_tau));

    results.push_back(I);

    auto t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    cout << "\rPeak Interpretation Method (PIM): Completed in " << (duration / 1000000.0) << " seconds.\n";
    return(results);
}

vector<double> Spectrum::SSHM(string libraryFolder,string library,double thresholdAB) 
{
    auto t1 = chrono::high_resolution_clock::now();
    cout << "Simple Search Hitlist Method (SS-HM): Working ";

    vector<double> result;
    int i;
    int j;

    double PIM_pred = PIM_core(thresholdAB);

    ofstream queryWriter("data\\Query_SS.MSP");
    if (!queryWriter) {
        cout << "Error opening file for writing query msp output." << endl;
    }
    queryWriter << "Name: Query" << endl;
    queryWriter << "Num Peaks: " << numpeaks << endl;
    for (i = 0; i < numpeaks; i++) {
        queryWriter << mz[i] << " " << ab[i] << "; ";
    }
    queryWriter.close();

    ofstream writer("batch-files\\asm-SS-HM.bat");
    if (!writer) {
        cout << "Error opening file for writing batch output." << endl;
    }
    writer << "@echo off" << endl;
    writer << "MSPepSearch\\2017_05_15_MSPepSearch\\x64\\MSPepSearch64.exe Sf^" <<endl;
    writer << " /LIB " << libraryFolder+library << "^" << endl;
    writer << " /HITS 25^" << endl;
    writer << " /INP data\\Query_SS.MSP^" << endl;
    writer << " /OutNumMP^" << endl;
    writer << " /OutMW^" << endl;
    writer << " /PROGRESS^" << endl;
    writer << " /LibInMem^" << endl;
    writer << " /OUTTAB data\\HitList_SS.txt" << endl;
    writer.close();

    int pepsearch = system("batch-files\\asm-SS-HM.bat >nul 2>nul");

    // code to read hit lists, identify spectra of library compounds (from msp files), 
    // determine "most probable correction" for query, and then add correction to PIM prediction of query.
    char letter;
    string line;
    string phrase;
    string libraryMSP;
    int i_libraryMW;
    int i_libraryMF;
    int i_libraryMWP;
    int i_libraryMWcorrection;

    vector<int> setCorrections;
    vector<int> setMFs;

    //cout << endl << endl << endl;

    int k = 0;
    ifstream hitlist("data\\HitList_SS.txt");

    if (! hitlist) {
        cout << "Error opening simple search hit list. " << endl;
    }
    else {
        for (i = 0; !hitlist.eof(); i++) {
            getline(hitlist,line);
            letter = line[0];
            if (letter=='Q') {
                k = 0;
                for (j = 0; j < line.length(); j++) {
                    phrase += line[j];
                    if (line[j] == '\t') {
                        k++;
                        if (k == 3) {
                            //cout << "library: " << phrase << endl;
                            phrase.pop_back();
                            libraryMSP = "searchLibrariesMSP\\";
                            libraryMSP += phrase;
                            phrase.clear();
                        }
                        else if (k == 4) {
                            //cout << "id: " << phrase << endl;
                            phrase.pop_back();
                            libraryMSP += "\\";
                            libraryMSP += phrase;
                            libraryMSP += ".MSP";
                            phrase.clear();
                        }
                        else if (k == 6) {
                            //cout << "mw: " << phrase << endl;
                            phrase.pop_back();
                            i_libraryMW = stoi(phrase);
                            phrase.clear();
                        }
                        else if (k == 7) {
                            //cout << "MF: " << phrase << endl;
                            phrase.pop_back();
                            i_libraryMF = stoi(phrase);
                            setMFs.push_back(i_libraryMF);
                            phrase.clear();
                        }
                        else {
                            phrase.clear();
                        }
                    }
                }


                Spectrum lib(libraryMSP);

                i_libraryMWP = lib.PIM_core(thresholdAB);
                i_libraryMWcorrection = i_libraryMW - i_libraryMWP;

                setCorrections.push_back(i_libraryMWcorrection);
                

                //cout << "MW prediction: " << i_libraryMWP << endl;
                //cout << "MW correction: " << i_libraryMWcorrection << endl;

                //cout << endl << endl;

                libraryMSP.clear();

            }  
        }
    }

    hitlist.close();

    vector<int> uniqueCorrections;
    //uniqueCorrections.push_back(0);

    for (i = 0; i < setCorrections.size(); i++) {
        if (find(uniqueCorrections.begin(), uniqueCorrections.end(), setCorrections.at(i)) == uniqueCorrections.end()) {
            uniqueCorrections.push_back(setCorrections.at(i));
        }
        //cout << setCorrections.at(i) << "\t" << setMFs.at(i) << endl;
    }

    //cout << "There are " << setCorrections.size() << " entries values in the hit list" << endl << endl;

    vector<double> correctionProb;

    double expNumerator = 0.0;
    double numerator = 0.0;
    double expDenominator = 0.0;
    double denominator = 0.0;
    int test = 0;

    for (i = 0; i < uniqueCorrections.size(); i++) {
        for (j = 0; j < setCorrections.size(); j++) {
            
            expDenominator = ( (double) setMFs.at(j) - (double) setMFs.at(0)) / (double) B;
            denominator = denominator + pow(2.0, expDenominator);

            if (setCorrections.at(j) == uniqueCorrections.at(i)) {
                expNumerator = ( (double) setMFs.at(j) - (double) setMFs.at(0))/ (double) B;
                numerator = numerator + pow(2.0, expNumerator);
                test = test + 1;
            }

        }
        //cout << numerator << " / " << denominator << endl;

        correctionProb.push_back(numerator / denominator);
        numerator = 0.0;
        denominator = 0.0;

        //cout <<endl << uniqueCorrections.at(i) << "\t" << correctionProb.at(i) << " using " << test << " scores." << endl <<endl;
        test = 0;
    }

    double cstar = 0.0;
    double pstar = 0.0;
    for (i = 0; i < correctionProb.size(); i++) {
        if (correctionProb.at(i) > pstar) {
            pstar = correctionProb.at(i);
            cstar = uniqueCorrections.at(i);
        }
    }

    result.push_back(PIM_pred + cstar);
    result.push_back(pstar);

    auto t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    cout << "\rSimple Search Hitlist Method (SS-HM): Completed in " << (duration / 1000000.0) << " seconds.\n";
    return(result);
}

vector<double> Spectrum::iHSHM(string libraryFolder,string library,double  thresholdAB) 
{
    auto t1 = chrono::high_resolution_clock::now();

    vector<double> result;
    int i = 0;
    int j = 0;
    int k = 0;
    int massIndex = 0;

    double PIM_pred = PIM_core(thresholdAB);
    
    double bp;
    bp = *max_element(ab.begin(), ab.end());
 
    int bp_mz=0;

    for (i = 0; i < ab.size(); i++) {
        if (ab.at(i) == bp) {
            bp_mz = mz.at(i);
        }
    }

    int minMW = PIM_pred - 18;
    int maxMW = PIM_pred + bp_mz;
    int rangeMW = maxMW - minMW;

    // variables to  read  in hit lists and compute Omega (improvement due to hybrid search).
    char letter;
    string line;
    string phrase;
    int i_DeltaMW;
    int i_hMF;
    int i_sMF;
    int Omega1 = 0;
    int Omega2 = 0;
    int OmegaCurrent = 0;
    int n_cognates = 0;
    int n_isomers = 0;
    double cognateOmega = 0.0;
    double isomerOmega = 0.0;
    double theta;

    //cout << "Range of mz values scanned: " << minMW << " to " << maxMW << endl;
    // for progress bar
    int step = 1;
    int displayNext = step;
    int percent = 0;


    for (massIndex = minMW; massIndex <= maxMW; massIndex++) {
        percent = (100 * (massIndex-minMW)) / rangeMW;
        if (percent >= displayNext)
        {
            cout << "\rIterative Hybrid Search Hitlist Method (iHS-HM): " << "[" << std::string(percent / 5, (char)254u) << std::string(100 / 5 - percent / 5, ' ') << "]";
            cout << percent << "%" << " [Iteration " << (massIndex - minMW) << " of " << rangeMW << "]";
            std::cout.flush();
            displayNext += step;
        }

        ofstream queryWriter("data\\Query_HS.MSP");
        if (!queryWriter) {
        cout << "Error opening file for writing query msp output." << endl;
    }
    
            queryWriter << "Name: Query" << endl;
            queryWriter << "MW: " << massIndex << endl;
            queryWriter << "Num Peaks: " << numpeaks << endl;

            for (i = 0; i < numpeaks; i++) {
                queryWriter << mz[i] << " " << ab[i] << "; ";
            }
        
         queryWriter.close();

         ofstream writer("batch-files\\asm-iHS-HM.bat");
         if (!writer) {
             cout << "Error opening file for writing batch output." << endl;
         }
         writer << "@echo off" << endl;
         writer << "MSPepSearch\\2017_05_15_MSPepSearch\\x64\\MSPepSearch64.exe Hf^" << endl;
         writer << " /LIB " << libraryFolder + library << "^" << endl;
         writer << " /HITS 25^" << endl;
         writer << " /INP data\\Query_HS.MSP^" << endl;
         writer << " /OutNumMP^" << endl;
         writer << " /OutMW^" << endl;
         writer << " /OutDeltaMW^" << endl;
         writer << " /PROGRESS^" << endl;
         writer << " /LibInMem^" << endl;
         writer << " /OUTTAB data\\HitList_HS.txt" << endl;
         writer.close();

         int pepsearch = system("batch-files\\asm-iHS-HM.bat >nul 2>nul");

         ifstream hitlist("data\\HitList_HS.txt");

         if (!hitlist) {
             cout << "Error opening Hybrid search hit list. " << endl;
         }
         else {
             for (i = 0; !hitlist.eof(); i++) {
                 getline(hitlist, line);
                 letter = line[0];
                 if (letter == 'Q') {
                     k = 0;
                     for (j = 0; j < line.length(); j++) {
                         phrase += line[j];
                         if (line[j] == '\t') {
                             k++;
                             if (k == 6) {
                                 //cout << "delta  mw: " << phrase << endl;
                                 phrase.pop_back();
                                 i_DeltaMW = stoi(phrase);
                                 phrase.clear();
                             }
                             else if (k == 8) {
                                 //cout << "hMF: " << phrase << endl;
                                 phrase.pop_back();
                                 i_hMF = stoi(phrase);
                                 phrase.clear();
                             }
                             else if (k == 12) {
                                 //cout << "sMF: " << phrase << endl;
                                 phrase.pop_back();
                                 i_sMF = stoi(phrase);
                                 phrase.clear();
                             }
                             else {
                                 phrase.clear();
                             }
                         }
                     }

                     if (i_hMF >= mEMF) {
                         if (i_DeltaMW == 0) {
                             n_isomers++;
                             isomerOmega = isomerOmega + (double)i_hMF * (((double)i_hMF - (double)mEMF) / (1000.0 - (double)mEMF));
                         }
                         else {
                             n_cognates++;
                             cognateOmega = cognateOmega + ((double)i_hMF - (double)i_sMF) * (((double)i_hMF - (double)mEMF) / (1000.0 - (double)mEMF));

                         }
                     }
                     else {
                         break;
                     }

                 }


             }

             hitlist.close();
             if (n_isomers + n_cognates > 0) {
                 OmegaCurrent = (isomerOmega + cognateOmega) / ((double)n_isomers + (double)n_cognates);
                 //cout << massIndex << "\t" << OmegaCurrent << endl;
             }
             else {
                 OmegaCurrent = 0.0;
             }
             

                 if (OmegaCurrent > Omega1) {
                     Omega2 = Omega1; 
                     Omega1 = OmegaCurrent;
                     theta = (double) massIndex;
                 }
                 else if (OmegaCurrent > Omega2) {
                     Omega2 = OmegaCurrent;
                 }
         }


         n_isomers = 0;
         n_cognates = 0;
         cognateOmega = 0.0;
         isomerOmega = 0.0;
    }

    result.push_back(theta);
    result.push_back((Omega1-Omega2)/999.0);

    auto t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    cout << "\rIterative Hybrid Search Hitlist Method (iHS-HM): Completed in " << (duration / 1000000.0) << " seconds.                                                                             \n";
   

    return(result);
}



void Launch() {
    cout << "Molecular Mass Predictor - Console Version\n";
    cout << "=========================================================\n\n";
}

int main()
{
    Launch();

    string queryFolder = "query-spectra-PaperDemonstrations\\";
    string libraryFolder = "searchLibraries\\";
    static const char* queryFolderpath = "query-spectra-PaperDemonstrations\\*.MSP";  
    static const char* libraryFolderpath = "searchLibraries\\*";
    
    vector<string> listofqueries;
    vector<string> listoflibraries;
    int i = 0, max_i = 0;
    int j = 0, max_j = 0; 
    int k = 0, max_k = 0;
    int qSelection = 0;
    int lSelection = 0;
    bool Continue;

    HANDLE hFind1;
    WIN32_FIND_DATAA data1;
    HANDLE hFind2;
    WIN32_FIND_DATAA data2;


    vector<double> PIM_info;
    vector<double> SSHM_info;
    vector<double> iHSHM_info;

    bool run_program = TRUE;
    char continue_program = 'yes';
    int runNum = 0;

    while (run_program) {

    
        hFind1 = FindFirstFileA(queryFolderpath, &data1);
        if (runNum == 0) {
            if (hFind1 != INVALID_HANDLE_VALUE)
            {
                do { listofqueries.push_back(data1.cFileName); } while (FindNextFileA(hFind1, &data1));

                FindClose(hFind1);
            }
        }
       
        Continue = FALSE;
        do {
        cout << "Available example query spectra: \n";
        max_i = listofqueries.size();
        for (i = 0; i < max_i; i++) {
            cout << i + 1 << ". " << listofqueries.at(i) << "\n";
        }
        cout << "\n";
        cout << "Select query spectra for analysis (integer value): ";
        cin >> qSelection;

        if (!cin.fail()) {
            if (qSelection > 0 && qSelection <= max_i) {
                cout << "\n\n";
                Continue = TRUE;
            }
            else {
                cout << "\nInvalid selection. See available integer options.\n\n";
                cin.clear();
            }
        }
        else {
            cout << "\nInvalid entry. Enter an integer value.\n\n";
            cin.clear();
            cin.ignore(10,'\n');
        }

        
    } while (Continue != TRUE);
        

        hFind2 = FindFirstFileA(libraryFolderpath, &data2);
        if(runNum==0){
        i = 0;
            if (hFind2 != INVALID_HANDLE_VALUE)
    {
        do {
            i++;
            if (i>3) { // A hacky way to ignore the current folder, parent folder and README file
                listoflibraries.push_back(data2.cFileName);
            }
        }
        while (FindNextFileA(hFind2, &data2));

        FindClose(hFind2);
    }
        }
        Continue = FALSE;
        do {

        cout << "Search libraries available: \n";
        max_i = listoflibraries.size();
        for (i = 0; i < max_i; i++) {
            cout << i + 1 << ". " << listoflibraries.at(i) << "\n";
        }

        cout << "\n";
        cout << "Select library for searching (integer value): ";
        cin >> lSelection;

        if (!cin.fail()) {
            if (lSelection > 0 && lSelection <= max_i) {
                std::cout << "\n";
                Continue = TRUE;
            }
            else {
                cout << "\nInvalid selection. See available integer options.\n\n";
                cin.clear();
            }
        }
        else {
            cout << "\nInvalid entry. Enter an integer value.\n\n";
            cin.clear();
            cin.ignore(10, '\n');
        }

    } while (Continue != TRUE);


        string queryfile = listofqueries[qSelection - 1];
        string library = listoflibraries[lSelection - 1];

        std::cout << "\nSearching " << queryfile << " against the " << library << " Library. \n\n";
        queryfile = queryFolder + queryfile;
        //cout << libraryFolder + library << endl;

        Spectrum query(queryfile);

        PIM_info = query.PIM(minMI_ABUNDANCE);
        SSHM_info = query.SSHM(libraryFolder,library,minMI_ABUNDANCE);
        iHSHM_info = query.iHSHM(libraryFolder,library,minMI_ABUNDANCE);

        std::cout << endl << "Prediction Results: \n";
        std::cout << "The PIM predicted Molecular Mass is " << PIM_info[0] << " Da with Classification Index " << PIM_info[1] << endl;
        std::cout << "The SS-HM predicted Molecular Mass is " << SSHM_info[0] << " Da with Classification Index " << SSHM_info[1] << endl;
        std::cout << "The iHS-HM predicted Molecular Mass is " << iHSHM_info[0] << " Da with Classification Index " << iHSHM_info[1] << endl;

        std::cout << "\n\nRun another example? (yes/no) ";
        std::cin >> continue_program;

        if (continue_program == 'yes' || continue_program == 'Yes' || continue_program == 'y' || continue_program == 'Y') {
            run_program = TRUE;
            runNum++;
            cout << endl << endl;
            continue_program = 'yes';
        }
        else {
            run_program = FALSE;
        }
        



    }



    return 0;
}


