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
    double crr = 1.03656e+00 + 
                 1.43300e+00 * std::exp(-2.44962e+00 * arig) +
                 1.43388e-01 * std::exp(-1.67525e-01 * arig) +
                 2.73303e-02 * std::exp(-1.18419e-02 * arig);
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
    TH1D*  hmcpr = (TH1D*)fmcpr->Get("hHie_cnt_MC");
    TTree* tmcpr = (TTree*)fmcpr->Get("runlist");
    tmcpr->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcpr->GetEntries(); ++it) { tmcpr->GetEntry(it); cntpr+=cntev; }
    
    UInt_t cntap = 0;
    TFile* fmcap = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcap%s/YiMdst.root", subv.c_str()));
    TH1D*  hmcap = (TH1D*)fmcap->Get("hHie_cnt_MC");
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
    TH1D*  hmcapp = (TH1D*)fmcapp->Get("hHie_cnt_MC");
    TTree* tmcapp = (TTree*)fmcapp->Get("runlist");
    tmcapp->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcapp->GetEntries(); ++it) { tmcapp->GetEntry(it); cntapp+=cntev; }
    
    UInt_t cntapm = 0;
    TFile* fmcapm = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/mcap_minus%s/YiMdst.root", subv.c_str()));
    TH1D*  hmcapm = (TH1D*)fmcapm->Get("hHie_cnt_MC");
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

    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/iss%s", subv.c_str()));

    PdfEditor editor(Window(), "apflux_iex", "out");

    const Axis& AXrig = Hist::Head("hHPie_llr")->xaxis();
 
    std::vector<THStack*> vhstack;
    std::vector<Hist*> vhist;

    std::vector<Hist*> vhllr_pr  = Hist::ProjectAll(HistProj::kY, Hist::Head("hHPie_llr_pr"));
    std::vector<Hist*> vhllr_el  = Hist::ProjectAll(HistProj::kY, Hist::Head("hHNie_llr_el"));
    std::vector<Hist*> vhllr_pos = Hist::ProjectAll(HistProj::kY, Hist::Head("hHPie_llr"));
    std::vector<Hist*> vhllr_neg = Hist::ProjectAll(HistProj::kY, Hist::Head("hHNie_llr"));
    
    Hist* hPcnt = Hist::New("hPcnt", HistAxis(AXrig));
    Hist* hNcnt = Hist::New("hNcnt", HistAxis(AXrig));
    Hist* hStat = Hist::New("hStat", HistAxis(AXrig));
    Hist* hSyst = Hist::New("hSyst", HistAxis(AXrig));
    Hist* hRate = Hist::New("hRate", HistAxis(AXrig));
    Hist* hCrrR = Hist::New("hCrrR", HistAxis(AXrig));

    for (int ir = 20; ir <= 60; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
/*
        
        Fit::RooVar roovar("llr", vhllr_neg.at(ir), HistList({ vhllr_pr.at(ir), vhllr_el.at(ir) }));
        Fit::RooSysResult rlt(roovar, true, 100);
        Fit::RooVar var = rlt.var();
        Fit::RooPar stdpar = rlt.std_par();
        Fit::RooPar syspar = rlt.sys_par();
        
        double int_pos = (*vhllr_pos.at(ir))()->Integral(roovar.min_bin(), roovar.max_bin());
        double int_neg = (*vhllr_neg.at(ir))()->Integral(roovar.min_bin(), roovar.max_bin());

        (*hPcnt)()->SetBinContent(ir, int_pos);
        (*hPcnt)()->SetBinError  (ir, std::sqrt(int_pos));
        (*hNcnt)()->SetBinContent(ir, stdpar.val(0));
        (*hNcnt)()->SetBinError  (ir, stdpar.err(0));

        (*hStat)()->SetBinContent(ir, stdpar.err(0) / int_pos);
        (*hSyst)()->SetBinContent(ir, syspar.err(0) / int_pos);
        
        (*hRate)()->SetBinContent(ir, stdpar.val(0) / int_pos);
        (*hRate)()->SetBinError  (ir, stdpar.err(0) / int_pos);
		
        var.samp() ->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
		var.sumt() ->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
		var.temp(1)->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
		var.temp(0)->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        
        editor.create();
       
        (*var.samp())()->GetXaxis()->CenterTitle();
        (*var.samp())()->GetXaxis()->SetTitle("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]");
        (*var.samp())()->GetYaxis()->SetTitle("Events/Bin");
        (*var.samp())()->SetMaximum( 1.3 * (*var.samp())()->GetMaximum() );

        (*var.samp() )()->Draw("pe");
        (*var.sumt() )()->Draw("hist same");
        (*var.temp(1))()->Draw("hist same");
        (*var.temp(0))()->Draw("hist same");
        
        Legend leg_table("", TextStyle(kBlack, 20, 43), PadWindow(0.58, 0.85, 0.65, 0.85));
        leg_table()->SetHeader(Form("Rigidity %.2f - %.2f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
        leg_table()->AddEntry((*var.samp())(), "Data", "lp");
        leg_table()->AddEntry((*var.sumt())(), "Sum", "lp");
        leg_table()->AddEntry((*var.temp(1))(), Form("e^{-}+#pi^{-}  (%.1f #pm %.1f)", stdpar.val(1), stdpar.err(1)), "lp");
        leg_table()->AddEntry((*var.temp(0))(), Form("p^{-}  (%.1f #pm %.1f)", stdpar.val(0), stdpar.err(0)), "lp");
        //leg_table()->SetTextFont(43);
        //leg_table()->SetTextSize(15);
        leg_table()->SetFillColor(0);
        leg_table.draw();

        editor.save();
      */ 

        DataFit1D fit1D(Form("%3d", ir), { (TH1D*)((*vhllr_pr.at(ir))()), (TH1D*)((*vhllr_el.at(ir))()) }, (TH1D*)((*vhllr_neg.at(ir))()), (TH1D*)((*vhllr_pos.at(ir))()), true);
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
        
        // TEST
        std::cerr << Form("RATE (FIT/STAT) %14.8f\n", rlt1D.num_errs().at(0) / std::sqrt(rlt1D.num_tmps().at(0)) );

        Hist* hsmp = Hist::New(rlt1D.hsmp().get());
        Hist* hsum = Hist::New(rlt1D.hsum().get());
        Hist* hsig = Hist::New(rlt1D.htmps().at(0).get());
        Hist* hbkg = Hist::New(rlt1D.htmps().at(1).get());
        Hist* href = Hist::New(rlt1D.href().get());
		
        hsmp->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
		hsum->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
		hbkg->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
		hsig->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        
        editor.create();
       
        (*hsmp)()->GetXaxis()->CenterTitle();
        (*hsmp)()->GetXaxis()->SetTitle("TRD Estimator");
        (*hsmp)()->GetYaxis()->SetTitle("Events/Bin");
        (*hsmp)()->SetMaximum( 1.3 * (*hsmp)()->GetMaximum() );

        (*hsmp)()->Draw("pe");
        (*hsum)()->Draw("hist same");
        (*hbkg)()->Draw("hist same");
        (*hsig)()->Draw("hist same");
        //(*hbkg)()->Draw("pe same");
        //(*hsig)()->Draw("pe same");
        
        Legend leg_table2("", TextStyle(kBlack, 20, 43), PadWindow(0.60, 0.85, 0.65, 0.85));
        leg_table2()->SetHeader(Form("Rigidity %.2f - %.2f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
        leg_table2()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table2()->AddEntry((*hsum)(), "Sum", "lp");
        //leg_table2()->AddEntry((*hbkg)(), Form("e^{-}+#pi^{-}  (%.1f #pm %.1f)", rlt1D.num_tmps().at(1), rlt1D.num_errs().at(1)), "lp");
        //leg_table2()->AddEntry((*hsig)(), Form("p^{-}  (%.1f #pm %.1f)", rlt1D.num_tmps().at(0), rlt1D.num_errs().at(0)), "lp");
        leg_table2()->AddEntry((*hbkg)(), Form("e^{-}+#pi^{-}  (%.0f)", rlt1D.num_tmps().at(1)), "lp");
        leg_table2()->AddEntry((*hsig)(), Form("p^{-}  (%.0f)", rlt1D.num_tmps().at(0)), "lp");
        //leg_table()->SetTextFont(43);
        //leg_table()->SetTextSize(15);
        leg_table2()->SetFillColor(0);
        leg_table2.draw();
        
        editor.save();
        
        std::cerr << Form("Rigidity %.2f - %.2f [GV/c]\n", AXrig()(ir-1), AXrig()(ir));
        std::cerr << Form("SMP %14.8f SIG %14.8f (%14.8f %14.8f) BKG %14.8f (%14.8f %14.8f) RATE %14.8f\n", rlt1D.num_smp(), rlt1D.num_tmps().at(0), rlt1D.wgt_tmps().at(0), rlt1D.wgt_errs().at(0), rlt1D.num_tmps().at(1), rlt1D.wgt_tmps().at(1), rlt1D.wgt_errs().at(1), rlt1D.num_tmps().at(0)/rlt1D.num_ref());
        std::cerr << Form("NCHI %14.8f\n", rlt1D.nchi());
        std::cerr << "\n";
    }

    editor.close();
    
    TFile * ofle = new TFile("out/apflux_iex.root", "RECREATE");
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
