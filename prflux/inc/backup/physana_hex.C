#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "DataFit1D.h"
#include "DataFit1D.C"

void stdfmt(TH1* hist) {
    if (hist == nullptr) return;
    hist->GetXaxis()->CenterTitle();
    hist->GetXaxis()->SetTitleFont(43);
    hist->GetXaxis()->SetTitleSize(20);
    hist->GetXaxis()->SetTitleOffset(2.5);
    hist->GetXaxis()->SetLabelFont(43);
    hist->GetXaxis()->SetLabelSize(15);
    hist->GetYaxis()->SetTitleFont(43);
    hist->GetYaxis()->SetTitleSize(20);
    hist->GetYaxis()->SetLabelFont(43);
    hist->GetYaxis()->SetLabelSize(15);
}

double accp_func(double arig) {
    double crr = 1.03311e+00 + 
                 2.20233e+00 * std::exp(-2.54541e+00 * arig) +
                 1.75704e-01 * std::exp(-1.98416e-01 * arig) +
                 3.46613e-02 * std::exp(-8.62238e-03 * arig);
    return crr;
}

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    std::string subv = "59";
    
    UInt_t cntev = 0;
    
    UInt_t cntpr = 0;
    TFile* fmcpr = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcpr%s/YiMdst.root", subv.c_str()));
    TH1D*  hmcpr = (TH1D*)fmcpr->Get("hHex_cnt_MC");
    TTree* tmcpr = (TTree*)fmcpr->Get("runlist");
    tmcpr->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcpr->GetEntries(); ++it) { tmcpr->GetEntry(it); cntpr+=cntev; }
    
    UInt_t cntap = 0;
    TFile* fmcap = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcap%s/YiMdst.root", subv.c_str()));
    TH1D*  hmcap = (TH1D*)fmcap->Get("hHex_cnt_MC");
    TTree* tmcap = (TTree*)fmcap->Get("runlist");
    tmcap->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcap->GetEntries(); ++it) { tmcap->GetEntry(it); cntap+=cntev; }

    TH1D* haccp = new TH1D("haccp", "", hmcpr->GetXaxis()->GetNbins(), hmcpr->GetXaxis()->GetXbins()->GetArray());
    TH1D* haerr = new TH1D("haerr", "", hmcpr->GetXaxis()->GetNbins(), hmcpr->GetXaxis()->GetXbins()->GetArray());
    for (int ib = 1; ib <= haccp->GetXaxis()->GetNbins(); ++ib) {
        double errpr = hmcpr->GetBinError(ib) / hmcpr->GetBinContent(ib);
        double errap = hmcap->GetBinError(ib) / hmcap->GetBinContent(ib);
        double error = std::sqrt(errpr * errpr + errap * errap);
        double accp = (hmcap->GetBinContent(ib) / static_cast<double>(cntap)) / (hmcpr->GetBinContent(ib) / static_cast<double>(cntpr));
        if (!std::isfinite(accp) || accp <= 0.0) continue;
        haccp->SetBinContent(ib, 1.0/accp);
        haccp->SetBinError  (ib, 1.0/accp * error);
        haerr->SetBinContent(ib, 1.0/accp - accp_func(haccp->GetXaxis()->GetBinCenter(ib)));
        haerr->SetBinError  (ib, 1.0/accp * error);
    }
    
    UInt_t cntapp = 0;
    TFile* fmcapp = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcap_plus%s/YiMdst.root", subv.c_str()));
    TH1D*  hmcapp = (TH1D*)fmcapp->Get("hHex_cnt_MC");
    TTree* tmcapp = (TTree*)fmcapp->Get("runlist");
    tmcapp->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcapp->GetEntries(); ++it) { tmcapp->GetEntry(it); cntapp+=cntev; }
    
    UInt_t cntapm = 0;
    TFile* fmcapm = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcap_minus%s/YiMdst.root", subv.c_str()));
    TH1D*  hmcapm = (TH1D*)fmcapm->Get("hHex_cnt_MC");
    TTree* tmcapm = (TTree*)fmcapm->Get("runlist");
    tmcapm->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcapm->GetEntries(); ++it) { tmcapm->GetEntry(it); cntapm+=cntev; }
    
    TH1D* hcross = new TH1D("hcross", "", hmcapp->GetXaxis()->GetNbins(), hmcapp->GetXaxis()->GetXbins()->GetArray());
    for (int ib = 1; ib <= hcross->GetXaxis()->GetNbins(); ++ib) {
        double errapp = hmcapp->GetBinError(ib) / hmcapp->GetBinContent(ib);
        double errapm = hmcapm->GetBinError(ib) / hmcapm->GetBinContent(ib);
        double error = std::sqrt(errapp * errapp + errapm * errapm);
        double cross = (hmcapp->GetBinContent(ib) / static_cast<double>(cntapp)) / (hmcapm->GetBinContent(ib) / static_cast<double>(cntapm));
        if (!std::isfinite(cross) || cross <= 0.0) continue;
        hcross->SetBinContent(ib, cross);
        hcross->SetBinError  (ib, cross * error);
    }
    
    //Hist* hcc_cc = Hist::New("hHNex_mva_cc", (TH1*)TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcpr_l1o9flux%s/YiMdst.root", subv.c_str()))->Get("hHNex3_mva_MC_FLUX27"));
    Hist* hcc_cc = Hist::New("hHNex_mva_cc", (TH1*)TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcpr_l1o9flux%s/YiMdst.root", subv.c_str()))->Get("hHNex3_mva_MC"));

    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/iss%s", subv.c_str()));

    PdfEditor editor(Window(), "apflux_hex", "out");

    const Axis& AXrig = Hist::Head("hHPex3_mva")->xaxis();
 
    std::vector<THStack*> vhstack;
    std::vector<Hist*> vhist;

    std::vector<Hist*> vhcc_cc = Hist::ProjectAll(HistProj::kY, hcc_cc);
    std::vector<Hist*> vhcc_pos = Hist::ProjectAll(HistProj::kY, Hist::Head("hHPex3_mva"));
    std::vector<Hist*> vhcc_neg = Hist::ProjectAll(HistProj::kY, Hist::Head("hHNex3_mva"));

    Hist* hPcnt = Hist::New("hPcnt", HistAxis(AXrig));
    Hist* hNcnt = Hist::New("hNcnt", HistAxis(AXrig));
    Hist* hStat = Hist::New("hStat", HistAxis(AXrig));
    Hist* hSyst = Hist::New("hSyst", HistAxis(AXrig));
    Hist* hRate = Hist::New("hRate", HistAxis(AXrig));
    Hist* hCrrR = Hist::New("hCrrR", HistAxis(AXrig));
    
    for (int ir = 30; ir <= AXrig.nbin()-4; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
       
        DataFit1D fit1D(Form("R%3d", ir), { (TH1D*)((*vhcc_pos.at(ir))()), (TH1D*)((*vhcc_cc.at(ir))()) }, (TH1D*)((*vhcc_neg.at(ir))()), (TH1D*)((*vhcc_pos.at(ir))()), true, "MVA", "Events/Bin");
        ResultFit1D rlt1D = fit1D.result();

        if (rlt1D.ndof() == 0) continue;
        
        (*hPcnt)()->SetBinContent(ir, rlt1D.num_ref()); 
        (*hNcnt)()->SetBinContent(ir, rlt1D.num_tmps().at(0)); 
        
        (*hStat)()->SetBinContent(ir, rlt1D.num_errs().at(0) / rlt1D.num_ref()); 
        (*hSyst)()->SetBinContent(ir, rlt1D.num_errs().at(0) / rlt1D.num_ref());

        (*hRate)()->SetBinContent(ir, rlt1D.num_tmps().at(0) / rlt1D.num_ref()); 
        (*hRate)()->SetBinError  (ir, rlt1D.num_errs().at(0) / rlt1D.num_ref()); 
        
        (*hCrrR)()->SetBinContent(ir, accp_func(rig) * rlt1D.num_tmps().at(0) / rlt1D.num_ref()); 
        //(*hCrrR)()->SetBinError  (ir, accp_func(rig) * rlt1D.num_errs().at(0) / rlt1D.num_ref()); 
        //(*hCrrR)()->SetBinError  (ir, accp_func(rig) * std::sqrt(rlt1D.num_tmps().at(0)) / rlt1D.num_ref()); 
        (*hCrrR)()->SetBinError  (ir, (*hCrrR)()->GetBinContent(ir) * std::hypot(rlt1D.num_errs().at(0)/rlt1D.num_tmps().at(0), 2*0.02*0.02)); 
        
        Hist* hsmpx = Hist::New(rlt1D.hsmp().get());
        Hist* hsumx = Hist::New(rlt1D.hsum().get());
        Hist* hsigx = Hist::New(rlt1D.htmps().at(0).get());
        Hist* hbkgx = Hist::New(rlt1D.htmps().at(1).get());
        Hist* hrefx = Hist::New(rlt1D.href().get());
		
        hsmpx->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
		hsumx->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
		hbkgx->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
		hsigx->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        
        editor.create();
       
        (*hsmpx)()->GetXaxis()->CenterTitle();
        (*hsmpx)()->SetMaximum( 1.3 * (*hsmpx)()->GetMaximum() );

        (*hsmpx)()->Draw("pe");
        (*hsumx)()->Draw("hist same");
        (*hbkgx)()->Draw("hist same");
        (*hsigx)()->Draw("hist same");
        
        Legend leg_tablex("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.65, 0.85));
        leg_tablex()->SetHeader(Form("Rigidity %.2f - %.2f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
        leg_tablex()->AddEntry((*hsmpx)(), "Data", "lp");
        leg_tablex()->AddEntry((*hsumx)(), "Sum", "lp");
        leg_tablex()->AddEntry((*hbkgx)(), Form("p^{+}  (%.1f #pm %.1f)", rlt1D.num_tmps().at(1), rlt1D.num_errs().at(1)), "lp");
        leg_tablex()->AddEntry((*hsigx)(), Form("p^{-}  (%.1f #pm %.1f)", rlt1D.num_tmps().at(0), rlt1D.num_errs().at(0)), "lp");
        leg_tablex()->SetFillColor(0);
        leg_tablex.draw();
        
        editor.save();
      
        std::cerr << Form("Rigidity %.2f - %.2f [GV/c]\n", AXrig()(ir-1), AXrig()(ir));
        std::cerr << Form("SMP %14.8f TEMPS %14.8f %14.8f REF %14.8f\n", fit1D.data_num_smp(), fit1D.data_num_tmp(0), fit1D.data_num_tmp(1), fit1D.data_num_ref());
        std::cerr << Form("SMP %14.8f SIG %14.8f (%14.8f %14.8f) BKG %14.8f (%14.8f %14.8f) RATE %14.8f\n", rlt1D.num_smp(), rlt1D.num_tmps().at(0), rlt1D.wgt_tmps().at(0), rlt1D.wgt_errs().at(0), rlt1D.num_tmps().at(1), rlt1D.wgt_tmps().at(1), rlt1D.wgt_errs().at(1), rlt1D.num_tmps().at(0)/rlt1D.num_ref());
        std::cerr << Form("NCHI %14.8f\n", rlt1D.nchi());
        std::cerr << "\n";
    }

    editor.close();
    
    TFile * ofle = new TFile("out/apflux_hex.root", "RECREATE");
    ofle->cd();

    haccp->Write();
    haerr->Write();
    hcross->Write();

    (*hPcnt)()->Write();
    (*hNcnt)()->Write();
    (*hStat)()->Write();
    (*hSyst)()->Write();
    (*hRate)()->Write();
    (*hCrrR)()->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
