

#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include <cassert>
#include <iostream>
#include <math.h>
#include <cmath>
#include <algorithm>

using namespace std;

class GenH1D: public TH1D {
    double pdfWeightRatioMin;
    double pdfWeightRatioMax;
    std::vector<TH1*> hPdfsCut;
    std::vector<TH1*> hPdfs;
    std::vector<TH1*> hPdfEvtCnts;
    TH1* hAlphasUp;
    TH1* hAlphasDwn;
    TH1* hPdfRejRate;
    TH1* hEvtCnt;
    static const unsigned int nScales = 6;
    std::vector<TH1*> hScales;
    bool uncHistBooked;
//    static const int scaleWeightMinIndex = 1; //index of 1st scale weight muR=muF=1 -> same value as nominal weight.
//    static const int pdfWeightMinIndex  =  scaleWeightMinIndex + 9; //weight index of 1st PDF replica
//    static const int alphasIndex = pdfWeightMinIndex + 100;

    static const int scaleWeightMinIndex = 1; // id 1001-1009 for variations of the renormalization and factorization scales
    static const int pdfWeightMinIndex   = 10; // id 2001-2100 (100 replicas) of NNPDF3.0
	static const int  alphasIndex = 110; // id 2101 and 2102 are the central PDF computed with alphas(mZ) = 0.117 and 0.119
	
public:
    
    GenH1D(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup):
    TH1D (name, title,  nbinsx, xlow, xup),
    pdfWeightRatioMin(.5), pdfWeightRatioMax(2), hPdfsCut(100), hPdfs(100),
    hPdfEvtCnts(100), hScales(nScales), uncHistBooked(false)
    {
    }
    
    GenH1D(const char *name, const char *title, Int_t nbinsx, const Float_t *xbins):
    TH1D (name, title,  nbinsx, xbins),
    pdfWeightRatioMin(.5), pdfWeightRatioMax(2), hPdfsCut(100), hPdfs(100),
    hPdfEvtCnts(100), hScales(nScales), uncHistBooked(false)
    {
    }
    
    GenH1D(const char *name, const char *title, Int_t nbinsx, const Double_t *xbins):
    TH1D (name, title,  nbinsx, xbins),
    pdfWeightRatioMin(.5), pdfWeightRatioMax(2), hPdfsCut(100), hPdfs(100),
    hPdfEvtCnts(100), hScales(nScales), uncHistBooked(false)
    {
    }
    
    
    Int_t Fill(Double_t x){ return TH1D::Fill(x); }
    Int_t Fill(Double_t x, Double_t w) { return TH1D::Fill(x, w); }
    Int_t Fill(const char *name, Double_t w) { return TH1D::Fill(name, w); }
    
    
    Int_t Fill(double x, double commonWeight, const std::vector<double>* weights){
        return Fill(x, commonWeight, *weights);
    }
    
    Int_t Fill(double x, double commonWeight, const std::vector<double>& weights){
        if(weights.size() < 6) return Fill(x, commonWeight * weights.at(0)); // for sherpa2 compatibility
        if(!uncHistBooked) bookUncHist();
        if((weights.size() != pdfWeightMinIndex + hPdfs.size())
           && (weights.size() != pdfWeightMinIndex + hPdfs.size() + 2)){
            std::cerr << "Unexpected weights.size() (" << weights.size() << "). "
            << "Expecting " << pdfWeightMinIndex + hPdfs.size()
            << " or " << pdfWeightMinIndex + hPdfs.size() + 2 << "\n";
            abort();
        }
		// fill central value
        //Int_t rc  = TH1D::Fill(x, commonWeight *  weights.at(0));
		Int_t rc  = TH1D::Fill(x, commonWeight *  weights.at(1));
        hEvtCnt->Fill(x, 1.);
        for(unsigned i = 0; i < hPdfsCut.size(); ++i){
            //double wr = weights.at(pdfWeightMinIndex + i) / weights.at(0);
			double wr = weights.at(pdfWeightMinIndex + i) / weights.at(1);
            if(pdfWeightRatioMin < wr && wr < pdfWeightRatioMax){
                hPdfsCut[i]->Fill(x, commonWeight * weights.at(pdfWeightMinIndex + i));
                hPdfEvtCnts[i]->Fill(x, 1.);
            } else{
                //std::cout << "Rejecting weight " << weights.at(pdfWeightMinIndex + i) << "\n";
                hPdfRejRate->Fill(x, 1);
            }
            hPdfs[i]->Fill(x, commonWeight * weights.at(pdfWeightMinIndex + i));
        }
        
        
        if(weights.size() > alphasIndex + 1){
            //Down and up might be swapped, but it does not matter for the unc. calculation
            //std::cout << commonWeight << ", " << weights.at(0) << ", "
            //<< weights.at(alphasIndex) << ", " << weights.at(alphasIndex+1) << "\n";
            hAlphasDwn->Fill(x, commonWeight * weights.at(alphasIndex));
            hAlphasUp->Fill(x, commonWeight * weights.at(alphasIndex+1));
        }
        
        int scaleIds[nScales] = {1, 2, 3, 4, 6, 8};//offset wrt to scaleWeightMinIndex;
        for(unsigned i = 0; i < nScales; ++i){
            hScales[i]->Fill(x, commonWeight * weights.at(scaleWeightMinIndex + scaleIds[i]));
        }
        
        return rc;
    }
    
    void addAlphasUnc(TGraphAsymmErrors* g){
        if(!g || hAlphasUp->GetEntries()==0) return;
        if(g->GetN() != hAlphasUp->GetNbinsX()
           || hAlphasUp->GetNbinsX() != hAlphasUp->GetNbinsX()){
            std::cerr << __FILE__ << ":" << __LINE__ << ":"
            << "Inconsistency in number of histogram bins. "
            << g->GetN() << ", " << hAlphasUp->GetNbinsX()
            << ", " << hAlphasUp->GetNbinsX() << ".\n";
            abort();
        }
        
        for(int i = 0; i < GetNbinsX(); ++i){
            double uncUp = hAlphasUp->GetBinContent(i+1) -  g->GetY()[i];
            double uncDwn = hAlphasDwn->GetBinContent(i+1) -  g->GetY()[i];
            if(uncUp < uncDwn) std::swap(uncUp, uncDwn);
            //std::cout << uncDwn << "\t" << uncUp << "\n";
            g->GetEYhigh()[i] = sqrt(std::pow(g->GetEYhigh()[i], 2) + uncUp*uncUp);
            g->GetEYlow()[i] = sqrt(std::pow(g->GetEYlow()[i], 2) + uncDwn*uncDwn);
        }
        g->SetTitle(TString::Format("%s + #alpha_S", g->GetTitle()));
    }
    
    
    
    
    TGraphAsymmErrors* getPdfUnc(bool withCut){
        std::vector<double> x(GetNbinsX(), 0);
        std::vector<double> y(GetNbinsX(), 0);
        std::vector<double> ey(GetNbinsX(), 0);
        std::vector<double> ex(GetNbinsX(), 0);
        for(int i = 0; i < GetNbinsX(); ++i){
            x[i] = GetXaxis()->GetBinCenter(i+1);
            ex[i] = 0.5*GetXaxis()->GetBinWidth(i+1);
            y[i] = GetBinContent(i+1);
            //calculate pdf uncertainty:
            double sum = 0;
            double sum2 = 0;
            unsigned nPdfs = hPdfsCut.size();
            for(unsigned iPdf = 0; iPdf < nPdfs; ++iPdf){
                double thisPdfY = 0;
                if(hPdfEvtCnts[iPdf]->GetBinContent(i+1) > 0){
                    //correction for the weights excluded by the pdfWeightRatio cut:
                    double corr = withCut ? (hEvtCnt->GetBinContent(i+1) / hPdfEvtCnts[iPdf]->GetBinContent(i+1)) : 1.;
                    std::vector<TH1*> hs = withCut ? hPdfsCut : hPdfs;
                    thisPdfY = hs[iPdf]->GetBinContent(i+1) * corr;
                }
                sum += thisPdfY;
                sum2 += thisPdfY*thisPdfY;
            }
            double s2 = (sum2  - sum*sum / nPdfs) / (nPdfs - 1);
            ey[i] = sqrt(s2);
        }
        
        TGraphAsymmErrors* r = new TGraphAsymmErrors(GetNbinsX(), &x[0], &y[0], &ex[0], &ex[0], &ey[0], &ey[0]);
        if(withCut){
            r->SetName(TString::Format("%s_pdfUnc", GetName()));
            r->SetTitle(TString::Format("%s PDF Unc.", GetTitle()));
        } else{
            r->SetName(TString::Format("%s_pdfUncPlain", GetName()));
            r->SetTitle(TString::Format("%s PDF Unc. (plain)", GetTitle()));
        }
        return r;
    }
    
    TGraphAsymmErrors* getPdfUnc68p(){
        std::vector<double> x(GetNbinsX(), 0);
        std::vector<double> y(GetNbinsX(), 0);
        std::vector<double> ey(GetNbinsX(), 0);
        std::vector<double> ex(GetNbinsX(), 0);
        for(int i = 0; i < GetNbinsX(); ++i){
            x[i] = GetXaxis()->GetBinCenter(i+1);
            ex[i] = 0.5*GetXaxis()->GetBinWidth(i+1);
            y[i] = GetBinContent(i+1);
            //calculate pdf uncertainty:
            unsigned nPdfs = hPdfsCut.size();
            std::vector<double> vals;
            for(unsigned iPdf = 0; iPdf < nPdfs; ++iPdf){
                vals.push_back(hPdfs[iPdf]->GetBinContent(i+1));
            }
            std::sort(vals.begin(), vals.end());
            ey[i] = 0.5 * (vals[83] - vals[15]);
        }
        
        TGraphAsymmErrors* r = new TGraphAsymmErrors(GetNbinsX(), &x[0], &y[0], &ex[0], &ex[0], &ey[0], &ey[0]);
        r->SetName(TString::Format("%s_pdfUnc68p", GetName()));
        r->SetTitle(TString::Format("%s PDF Unc. (68p)", GetTitle()));
        return r;
    }
    
    
    TGraphAsymmErrors* getScaleUnc(){
        std::vector<double> x(GetNbinsX(), 0);
        std::vector<double> y(GetNbinsX(), 0);
        std::vector<double> ehx(GetNbinsX(), 0);
        std::vector<double> elx(GetNbinsX(), 0);
        std::vector<double> ehy(GetNbinsX(), 0);
        std::vector<double> ely(GetNbinsX(), 0);
        for(int i = 0; i < GetNbinsX(); ++i){
            x[i] = GetXaxis()->GetBinCenter(i+1);
            elx[i] = ehx[i] = 0.5 * GetXaxis()->GetBinWidth(i+1);
            y[i] = GetBinContent(i+1);
            //calculate uncertainty with scale variation:
            double down = y[i];
            double up = y[i];
            unsigned nScales = hScales.size();
            for(unsigned iScale = 0; iScale < nScales; ++iScale){
                double thisScaleY = hScales[iScale]->GetBinContent(i+1);
                if(thisScaleY < down) down = thisScaleY;
                if(thisScaleY > up)   up = thisScaleY;
            }
            ehy[i] = up - y[i];
            ely[i] = y[i] - down;
        }
        
        TGraphAsymmErrors* r = new TGraphAsymmErrors(GetNbinsX(), &x[0], &y[0], &elx[0], &ehx[0], &ely[0], &ehy[0]);
        r->SetName(TString::Format("%s_scaleUnc", GetName()));
        r->SetTitle(TString::Format("%s Scale Unc.", GetTitle()));
        return r;
    }
    
    Int_t Write ( const char *  name = 0,
                 Int_t  option = 0,
                 Int_t  bufsize = 0
                 ){
        Int_t rc = TH1::Write(name, option, bufsize);
        
        if(uncHistBooked){
            TGraphAsymmErrors* pdfUnc = getPdfUnc(true);
            addAlphasUnc(pdfUnc);
            pdfUnc->Write();
            delete pdfUnc;
            TGraphAsymmErrors* pdfUncPlain = getPdfUnc(false);
            addAlphasUnc(pdfUncPlain);
            pdfUncPlain->Write();
            delete pdfUncPlain;
            TGraphAsymmErrors* pdfUnc68p = getPdfUnc68p();
            addAlphasUnc(pdfUnc68p);
            pdfUnc68p->Write();
            delete pdfUnc68p;
            TGraphAsymmErrors* scaleUnc = getScaleUnc();
            scaleUnc->Write();
            delete scaleUnc;
            hPdfRejRate->Divide(hEvtCnt);
            hPdfRejRate->Scale(0.01);
            hPdfRejRate->Write();
        }
        
        return rc;
    }
    
    void Scale(Double_t c1 = 1, Option_t *option = ""){
        TH1::Scale(c1, option);
        if(uncHistBooked){
            for(unsigned i = 0; i < hPdfsCut.size(); ++i){
                hPdfsCut[i]->Scale(c1, option);
                hPdfs[i]->Scale(c1, option);
            }
            
            for(unsigned i = 0; i < nScales; ++i){
                hScales[i]->Scale(c1, option);
            }
            hAlphasDwn->Scale(c1, option);
            hAlphasUp->Scale(c1, option);
        }
    }
    
private:
    void bookUncHist(){
        uncHistBooked = true;
        hEvtCnt = (TH1D*) Clone();
        hEvtCnt->SetDirectory(0);
        hEvtCnt->SetName(TString::Format("%s_evtcnt", GetName()));
        hEvtCnt->SetTitle(TString::Format("%s Evt Cnt", GetName()));
        hPdfRejRate = (TH1D*) Clone();
        hPdfRejRate->SetName(TString::Format("%s_rej", GetName()));
        hPdfRejRate->SetTitle(TString::Format("%s Pdf weight rejection rate", GetName()));
        hPdfRejRate->GetYaxis()->SetTitle("PDF weight rejection rate");
        for(unsigned i = 0; i < hPdfsCut.size(); ++i){
            hPdfsCut[i] = (TH1D*) Clone();
            hPdfsCut[i]->SetDirectory(0);
            hPdfsCut[i]->SetName(TString::Format("%s_pdf_cut%d", GetName(), i));
            hPdfEvtCnts[i] = (TH1D*) Clone();
            hPdfEvtCnts[i]->SetDirectory(0);
            hPdfEvtCnts[i]->SetName(TString::Format("%s_pdf%d_evtcnt", GetName(), i));
            hPdfs[i] = (TH1D*) Clone();
            hPdfs[i]->SetDirectory(0);
            hPdfs[i]->SetName(TString::Format("%s_pdf%d", GetName(), i));
        }
        
        for(unsigned i = 0; i < hScales.size(); ++i){
            hScales[i] = (TH1D*) Clone();
            hScales[i]->SetDirectory(0);
            hScales[i]->SetName(TString::Format("%s_scale%d", GetName(), i));
        }
        
        hAlphasUp = (TH1D*) Clone();
        hAlphasUp->SetName(TString::Format("%s_alphasUp", GetName()));
        hAlphasUp->SetTitle(TString::Format("%s #alpha_s +1#sigma", GetTitle()));
        //	hAlphasUp->SetDirectory(0);
        hAlphasDwn = (TH1D*) Clone();
        hAlphasDwn->SetName(TString::Format("%s_alphasDwn", GetName()));
        hAlphasDwn->SetTitle(TString::Format("%s #alpha_s -1#sigma", GetTitle()));
        //	hAlphasDwn->SetDirectory(0);
    }
    
};
