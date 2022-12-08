#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"

class AziNorm {
    public:
        AziNorm() {}
        AziNorm(GH2D *_source, GH2D *_data){
            source = _source;
            data = _data;
        }

        void SetSourcePeak(double _xlo, double _xhi){
            source_peak = std::make_pair(_xlo,_xhi);
        }

        void AddPeak(double _xlo, double _xhi) {
            if (_xlo > _xhi) std::swap(_xlo,_xhi);

            //check if region is already added
            bool nodupe = true;
            for (auto dp : data_peaks){
                if (_xlo*10000 + _xhi == dp.first*10000 + dp.second) nodupe = false;
            }   
            if (nodupe) data_peaks.push_back(std::make_pair(_xlo,_xhi));
            else printf("%4.0f <> %-4.0f peak already added!\n",_xlo,_xhi);
        }

        void LoadPeaks(std::string filename){
            std::ifstream input(filename.c_str());
            std::string line;
            double xlo, xhi;
            while (getline(input,line)){
                if (line.find('#') != std::string::npos) continue;
                std::stringstream ss(line);
                ss >> xlo >> xhi;
                AddPeak(xlo,xhi);
                printf("Added %4.0f <> %-4.0f\n",xlo,xhi);
            }
        }

        void ListPeaks(){
            for (auto dp : data_peaks){
                printf("%4.0f <> %-4.0f\n",dp.first,dp.second);
            }
        }

        void ClearPeaks(){data_peaks.clear();}

        void GetRatioHists(){
            TFile *fout = new TFile("azimithal_compton_ratios.root","RECREATE");

            int binLo = source->GetYaxis()->FindBin(source_peak.first);
            int binHi = source->GetYaxis()->FindBin(source_peak.second);
            GH1D *s = data->ProjectionX(Form("source_peak_%5.1f-%5.1f",source_peak.first,source_peak.second),binLo,binHi);
            s->Sumw2();

            for (auto dp : data_peaks){
                binLo = data->GetYaxis()->FindBin(dp.first);
                binHi = data->GetYaxis()->FindBin(dp.second);
                GH1D *h = data->ProjectionX(Form("peak_%5.1f-%5.1f_ratio",dp.first,dp.second),binLo,binHi);
                h->Sumw2();
                h->Divide(s);
                h->Write();
            }
            fout->Close();
        }

    private:
        GH2D *source;
        GH2D *data;
        std::pair<double, double> source_peak;
        vector<std::pair<double,double>> data_peaks;
};