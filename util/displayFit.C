#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

bool readPeakFile(std::string list, std::vector<std::string> &filenames, std::vector<double> &scales){
    bool flag = true;
    ifstream input(list);
    if (input.is_open()){
        std::string line;
        std::string word;
        while(getline(input,line)){
            istringstream ss(line);
            ss >> word;
            word = "hist" + word + ".root";
            filenames.push_back(word);
            for (int i=0; i < 2; i++) ss >> word;
            scales.push_back(std::atof(word.c_str()));
        }
    }else{
        std::cout<<"error opening "<<list<<std::endl;
        flag = false;
    }
    return flag;
}

void displayFit(std::string histfile, std::string peakList){
    std::vector <std::string> templateHistFileNames;
    std::vector <double> scales;
    readPeakFile(peakList,templateHistFileNames,scales);
    std::cout<<templateHistFileNames[0]<<" "<<scales[0]<<std::endl;
}