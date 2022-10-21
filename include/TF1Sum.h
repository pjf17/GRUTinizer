#include <vector>
#include <string>

#include <TF1.h>

class TF1Sum : public TNamed {


  public:
    TF1Sum():TNamed("TF1Sum","TF1Sum"),fFit(0),npars(0) { }
    TF1Sum(const TF1Sum& other) 
      : TNamed(other), fFit(0), npars(other.npars), xlow(other.xlow), xhigh(other.xhigh),
        exclude_low(0), exclude_high(0), fTF1s(other.fTF1s)  { }

    ~TF1Sum() { if(fFit) delete fFit; } 

    void AddTF1(TF1 *f);
 
    void Print(Option_t *opt="") const;

    Double_t EvalPar(const Double_t *x,const Double_t *params=0);

    double operator()(double *x,double *params) { 
      int parnum = 0;
      double sum = 0.0;

      if (exclude_low != 0 || exclude_high != 0){
        if (x[0] > exclude_low && x[0] < exclude_high){
          TF1::RejectPoint();
        }
      }
      for(auto fit : fTF1s) {
        //printf("fit->GetNpar() = %i\n",fit->GetNpar()); fflush(stdout);
        if(params==0) sum += fit->EvalPar(x,params);
        else sum += fit->EvalPar(x, params+parnum);

        parnum += fit->GetNpar();
      }
      return sum;
    }

    void SetRange(double l,double h) { xlow =l; xhigh=h; }
    void SetExcludeRange(double l, double h){exclude_low = l; exclude_high = h;}
    int  GetNpar() const  { return npars; }

    operator TF1*() { return fFit;}
        
    TF1 *GetFunc() { return fFit; }

    virtual void Draw(Option_t *opt="all");


  private:
    TF1 *fFit;
    int npars;
    double xlow;
    double xhigh;
    double exclude_low;
    double exclude_high;

    std::vector<double> fParam;
    std::vector<double> fParErr;
    std::vector<double> fParMax;
    std::vector<double> fParMin;
    std::vector<std::string> fParName;
    
    std::vector<TF1*>    fTF1s;

    ClassDef(TF1Sum,0);
};

