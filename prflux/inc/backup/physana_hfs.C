#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "DataFit1D.h"
#include "DataFit1D.C"

#include "HistFit1D.h"
#include "HistFit1D.C"

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
    double crr = 1.03171e+00 + 
                 3.95207e+00 * std::exp(-2.70251e+00 * arig) +
                 2.78767e-01 * std::exp(-2.58810e-01 * arig) +
                 4.21401e-02 * std::exp(-1.04794e-02 * arig);
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
    TH1D*  hmcpr = (TH1D*)fmcpr->Get("hHfs_cnt_MC");
    TTree* tmcpr = (TTree*)fmcpr->Get("runlist");
    tmcpr->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcpr->GetEntries(); ++it) { tmcpr->GetEntry(it); cntpr+=cntev; }
    
    UInt_t cntap = 0;
    TFile* fmcap = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcap%s/YiMdst.root", subv.c_str()));
    TH1D*  hmcap = (TH1D*)fmcap->Get("hHfs_cnt_MC");
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
    TH1D*  hmcapp = (TH1D*)fmcapp->Get("hHfs_cnt_MC");
    TTree* tmcapp = (TTree*)fmcapp->Get("runlist");
    tmcapp->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcapp->GetEntries(); ++it) { tmcapp->GetEntry(it); cntapp+=cntev; }
    
    UInt_t cntapm = 0;
    TFile* fmcapm = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcap_minus%s/YiMdst.root", subv.c_str()));
    TH1D*  hmcapm = (TH1D*)fmcapm->Get("hHfs_cnt_MC");
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
    
    //Hist* hcc_cc = Hist::New("hHNfs_mva_cc", (TH1*)TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcpr_l1o9flux%s/YiMdst.root", subv.c_str()))->Get("hHNfs3_mva_MC_FLUX27"));
    Hist* hcc_cc = Hist::New("hHNfs_mva_cc", (TH1*)TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcpr_l1o9flux%s/YiMdst.root", subv.c_str()))->Get("hHNfs_mva_MC"));

    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/iss%s", subv.c_str()));

    PdfEditor editor(Window(), "apflux_hfs", "out");

    const Axis& AXrig = Hist::Head("hHPfs_mva")->xaxis();
 
    std::vector<THStack*> vhstack;
    std::vector<Hist*> vhist;

    std::vector<Hist*> vhcc_cc = Hist::ProjectAll(HistProj::kY, hcc_cc);
    std::vector<Hist*> vhcc_pos = Hist::ProjectAll(HistProj::kY, Hist::Head("hHPfs_mva"));
    std::vector<Hist*> vhcc_neg = Hist::ProjectAll(HistProj::kY, Hist::Head("hHNfs_mva"));

    Hist* hPcnt = Hist::New("hPcnt", HistAxis(AXrig));
    Hist* hNcnt = Hist::New("hNcnt", HistAxis(AXrig));
    Hist* hStat = Hist::New("hStat", HistAxis(AXrig));
    Hist* hSyst = Hist::New("hSyst", HistAxis(AXrig));
    Hist* hRate = Hist::New("hRate", HistAxis(AXrig));
    Hist* hCrrR = Hist::New("hCrrR", HistAxis(AXrig));
    
    HistFit::Axis1D axis(
        "MVA",
        "Events/Bin",
        (*Hist::Head("hHPfs_mva"))()->GetYaxis(), 100, 200);
    
    for (int ir = 40; ir <= AXrig.nbin()-2; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        std::string prefix = Form("BIN%03d_", ir);

        HistFit::Hist1D h1Dref(Form("%sREF", prefix.c_str()), "REF", (TH1D*)((*vhcc_pos.at(ir))()), axis, false);
        HistFit::HistFit1D fit1D(
            (TH1D*)((*vhcc_neg.at(ir))()), 
            { (TH1D*)((*vhcc_pos.at(ir))()), (TH1D*)((*vhcc_cc.at(ir))()) }, 
            axis, prefix);
        if (!fit1D.status()) continue;

        const HistFit::Hist1D& h1Dsmp = fit1D.ref_smp();
        const HistFit::Hist1D& h1Dsum = fit1D.sum_tmps();
        const HistFit::Hist1D& h1Dsig = fit1D.wgt_tmps(0);
        const HistFit::Hist1D& h1Dbkg = fit1D.wgt_tmps(1);

        double num_pr = h1Dref.data().sum();
        double num_ap = fit1D.wgts(0);
        double num_cc = fit1D.wgts(1);

        double err_stat = fit1D.errs(0);
        double err_fluc = fit1D.fluc(0).err;
        double err_syst = err_fluc;

        double app_val  = accp_func(rig) * (num_ap / num_pr);
        double app_stat = app_val * (err_stat / num_ap);
        double app_fluc = app_val * (err_fluc / num_ap);
        double app_syst = app_val * (err_syst / num_ap);

        (*hPcnt)()->SetBinContent(ir, num_pr); 
        (*hNcnt)()->SetBinContent(ir, num_ap); 
        
        (*hStat)()->SetBinContent(ir, app_stat / app_val); 
        (*hSyst)()->SetBinContent(ir, app_syst / app_val);

        (*hCrrR)()->SetBinContent(ir, app_val); 
        (*hCrrR)()->SetBinError  (ir, app_stat); 

        Hist* hsmp = Hist::New(h1Dsmp.get());
        Hist* hsum = Hist::New(h1Dsum.get());
        Hist* hsig = Hist::New(h1Dsig.get());
        Hist* hbkg = Hist::New(h1Dbkg.get());

        hsmp->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
		hsum->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
		hsig->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
		hbkg->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
        
        THStack* hfit = Hist::Collect("hfit", HistList({ hsum, hbkg, hsig }));
        
        editor.create();
   
        hfit->Draw("nostack hist");
        (*hsmp)()->Draw("pe same");

        //(*hsmp)()->GetXaxis()->CenterTitle();
        //(*hsmp)()->SetMaximum( 1.3 * (*hsmp)()->GetMaximum() );
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.65, 0.85));
        leg_table()->SetHeader(Form("Rigidity %.2f - %.2f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
        leg_table()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table()->AddEntry((*hsum)(), "Sum", "lp");
        leg_table()->AddEntry((*hbkg)(), Form("p^{+}  (%.1f #pm %.1f)", fit1D.wgts(1), fit1D.errs(1)), "lp");
        leg_table()->AddEntry((*hsig)(), Form("p^{-}  (%.1f #pm %.1f)", fit1D.wgts(0), fit1D.errs(0)), "lp");
        leg_table()->SetFillColor(0);
        leg_table.draw();
        
        editor.save();
   
        std::cerr << Form("Rigidity %.2f - %.2f [GV/c]\n", AXrig()(ir-1), AXrig()(ir));
        std::cerr << Form("SMP %14.8f SIG %14.8f (%14.8f %14.8f) BKG %14.8f (%14.8f %14.8f) Ratio %14.8f %14.8f\n", 
                          h1Dsmp.data().sum(),
                          fit1D.wgts(0), fit1D.errs(0), fit1D.fluc(0).err,
                          fit1D.wgts(1), fit1D.errs(1), fit1D.fluc(1).err,
                          app_val, app_stat);
        std::cerr << Form("NCHI %14.8f\n", fit1D.nchi());
        std::cerr << "\n";
    }

    editor.close();
    
    TFile * ofle = new TFile("out/apflux_hfs.root", "RECREATE");
    ofle->cd();

    haccp->Write();
    haerr->Write();
    hcross->Write();

    (*hPcnt)()->Write();
    (*hNcnt)()->Write();
    (*hStat)()->Write();
    (*hSyst)()->Write();
    (*hCrrR)()->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
