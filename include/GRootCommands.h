#ifndef GROOTCOMMANDS__H
#define GROOTCOMMANDS__H

class TH1;

#include "TDirectory.h"

int  LabelPeaks(TH1*,double,double,Option_t *opt="");
bool ShowPeaks(TH1**,unsigned int);
bool RemovePeaks(TH1**,unsigned int);

//bool PeakFit(TH1*,Double_t,Double_t,Option_t *opt="");

void Help();
void Commands();

#endif