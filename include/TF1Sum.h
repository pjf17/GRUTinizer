#include <vector>
#include <string>

#include <TF1.h>

class TF1Sum : public TNamed {


  public:
    TF1Sum():TNamed("TF1Sum","TF1Sum"),fFit(0),npars(0) { }
    TF1Sum(const TF1Sum& other) 
      : TNamed(other), fFit(0), npars(other.npars), xlow(other.xlow), xhigh(other.xhigh),
        fTF1s(other.fTF1s)  { }

    ~TF1Sum() { if(fFit) delete fFit; } 

    void AddTF1(TF1 *f);
 
    void Print(Option_t *opt="") const;

    //Double_t EvalPar(const Double_t *x,const Double_t *params=0);

    double operator()(double *x,double *params) { 
      int parnum = 0;
      double sum = 0.0;

      if (regions.size() > 0){
        bool reject = !rejectSwitch; //reject switch is default true
        for (auto r : regions){
          if (x[0] >= r.first && x[0] <= r.second) reject = rejectSwitch;
        }
        if (reject) TF1::RejectPoint();
      }

      for(auto fit : fTF1s) {
        if(params==0) sum += fit->EvalPar(x,params);
        else sum += fit->EvalPar(x, params+parnum);

        parnum += fit->GetNpar();
      }
      return sum;
    }

    void SetRange(double l,double h) { xlow =l; xhigh=h; }
    void AddRegion(double l, double h){regions.push_back(std::make_pair(l,h));}
    void FitInRegions(){rejectSwitch = false;}
    void FitOutsideRegions(){rejectSwitch = true;}
    int  GetNpar() const  { return npars; }

    operator TF1*() { return fFit;}
        
    TF1 *GetFunc() { return fFit; }

    virtual void Draw(Option_t *opt="all");


  private:
    TF1 *fFit;
    int npars;
    double xlow;
    double xhigh;
    // double exclude_low;
    // double exclude_high;

    std::vector<double> fParam;
    std::vector<double> fParErr;
    std::vector<double> fParMax;
    std::vector<double> fParMin;
    std::vector<std::string> fParName;
    
    std::vector<TF1*>    fTF1s;

    std::vector<std::pair<double,double>> regions;
    bool rejectSwitch = true;

    ClassDef(TF1Sum,1);
};

