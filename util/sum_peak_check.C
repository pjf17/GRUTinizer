#include <string>
#include <cmath>
#include <vector>
#include <iostream>

bool GetPeaks(std::string filename, std::vector<int> &peaks){
  ifstream input(filename.c_str());

  if (!input.is_open()){
    std::cout<<"Error: cannot open file\n";
    return false;
  }

  std::string line;
  while ( getline(input,line) ){
    peaks.push_back(std::stoi(line));
  }
  
  input.close();
  return true;
}

void sum_peak_check(){
    std::vector<int> peaks;
    if (!GetPeaks("peaks.txt",peaks)){
        return;
    }
    int npeaks = (int) peaks.size();
    for (int i=0; i < npeaks; i++){
        for (int j=0; j < npeaks; j++){
            if (i==j) continue;
            for (int k=0; k < npeaks; k++){
                if (k == i || k <= j) continue;
                if (std::abs(peaks[i] - (peaks[j] + peaks[k])) < 5){
                    std::cout<<peaks[k]<<" + "<<peaks[j]<<" = "<<peaks[k]+peaks[j]<<" ~ "<<peaks[i]<<std::endl;
                }
            }
        }
    }
}