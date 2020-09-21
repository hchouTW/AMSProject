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
    double crr = 1.04752e+00 + 2.91583e-02 * std::exp(-1.28208e-02 * arig) + 1.60814e-01 * std::exp(-1.81531e-01 * arig);
    return crr;
}

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    std::string subv = "38";
    
    UInt_t cntev = 0;
    
    UInt_t cntpr = 0;
    TFile* fmcpr = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcpr%s/YiMdst.root", subv.c_str()));
    TH1D*  hmcpr = (TH1D*)fmcpr->Get("hI_MC_cnt");
    TH1D*  hTmcpr = (TH1D*)fmcpr->Get("hTI_MC_cnt");
    TTree* tmcpr = (TTree*)fmcpr->Get("ana");
    tmcpr->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcpr->GetEntries(); ++it) { tmcpr->GetEntry(it); cntpr+=cntev; }
    
    UInt_t cntap = 0;
    TFile* fmcap = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcap%s/YiMdst.root", subv.c_str()));
    TH1D*  hmcap = (TH1D*)fmcap->Get("hI_MC_cnt");
    TH1D*  hTmcap = (TH1D*)fmcap->Get("hTI_MC_cnt");
    TTree* tmcap = (TTree*)fmcap->Get("ana");
    tmcap->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcap->GetEntries(); ++it) { tmcap->GetEntry(it); cntap+=cntev; }

    TH1D* haccp = new TH1D("haccp", "", hmcpr->GetXaxis()->GetNbins(), hmcpr->GetXaxis()->GetXbins()->GetArray());
    for (int ib = 1; ib <= haccp->GetXaxis()->GetNbins(); ++ib) {
        double errpr = hmcpr->GetBinError(ib) / hmcpr->GetBinContent(ib);
        double errap = hmcap->GetBinError(ib) / hmcap->GetBinContent(ib);
        double error = std::sqrt(errpr * errpr + errap * errap);
        double accp = (hmcap->GetBinContent(ib) / static_cast<double>(cntap)) / (hmcpr->GetBinContent(ib) / static_cast<double>(cntpr));
        if (!std::isfinite(accp) || accp <= 0.0) continue;
        haccp->SetBinContent(ib, 1.0/accp);
        haccp->SetBinError  (ib, 1.0/accp * error);
    }
    
    TH1D* hTaccp = new TH1D("hTaccp", "", hTmcpr->GetXaxis()->GetNbins(), hTmcpr->GetXaxis()->GetXbins()->GetArray());
    for (int ib = 1; ib <= hTaccp->GetXaxis()->GetNbins(); ++ib) {
        double errpr = hTmcpr->GetBinError(ib) / hTmcpr->GetBinContent(ib);
        double errap = hTmcap->GetBinError(ib) / hTmcap->GetBinContent(ib);
        double error = std::sqrt(errpr * errpr + errap * errap);
        double accp = (hTmcap->GetBinContent(ib) / static_cast<double>(cntap)) / (hTmcpr->GetBinContent(ib) / static_cast<double>(cntpr));
        if (!std::isfinite(accp) || accp <= 0.0) continue;
        hTaccp->SetBinContent(ib, 1.0/accp);
        hTaccp->SetBinError  (ib, 1.0/accp * error);
    }

    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/iss%s", subv.c_str()));

    PdfEditor editor(Window(), "antip_i", "out");

    const Axis& AXrig = Hist::Head("hIP_llr")->xaxis();
    const Axis& AXrig2 = Hist::Head("hTIP_llr")->yaxis();
    const Axis& AXtme = Hist::Head("hTIP_llr")->xaxis();
 
    std::vector<THStack*> vhstack;
    std::vector<Hist*> vhist;

    std::vector<Hist*> vhllr_pr  = Hist::ProjectAll(HistProj::kY, Hist::Head("hIP_llr_pr"));
    std::vector<Hist*> vhllr_el  = Hist::ProjectAll(HistProj::kY, Hist::Head("hIN_llr_el"));
    std::vector<Hist*> vhllr_pos = Hist::ProjectAll(HistProj::kY, Hist::Head("hIP_llr"));
    std::vector<Hist*> vhllr_neg = Hist::ProjectAll(HistProj::kY, Hist::Head("hIN_llr"));
    
    std::vector<Hist*> vhTllr_pr  = Hist::ProjectAll(HistProj::kZY, Hist::Head("hTIP_llr_pr"));
    std::vector<Hist*> vhTllr_el  = Hist::ProjectAll(HistProj::kZY, Hist::Head("hTIN_llr_el"));
    std::vector<Hist*> vhTllr_pos = Hist::ProjectAll(HistProj::kZY, Hist::Head("hTIP_llr"));
    std::vector<Hist*> vhTllr_neg = Hist::ProjectAll(HistProj::kZY, Hist::Head("hTIN_llr"));
        
    Hist* hPcnt = Hist::New("hPcnt", HistAxis(AXrig));
    Hist* hNcnt = Hist::New("hNcnt", HistAxis(AXrig));
    Hist* hStat = Hist::New("hStat", HistAxis(AXrig));
    Hist* hSyst = Hist::New("hSyst", HistAxis(AXrig));
    Hist* hRate = Hist::New("hRate", HistAxis(AXrig));
    Hist* hCrrR = Hist::New("hCrrR", HistAxis(AXrig));

    Hist* hTPcnt = Hist::New("hTPcnt", HistAxis(AXtme, AXrig2));
    Hist* hTNcnt = Hist::New("hTNcnt", HistAxis(AXtme, AXrig2));
    Hist* hTStat = Hist::New("hTStat", HistAxis(AXtme, AXrig2));
    Hist* hTSyst = Hist::New("hTSyst", HistAxis(AXtme, AXrig2));
    Hist* hTRate = Hist::New("hTRate", HistAxis(AXtme, AXrig2));
    Hist* hTCrrR = Hist::New("hTCrrR", HistAxis(AXtme, AXrig2));
    
    (*hTRate)()->GetXaxis()->SetTitle("Date");
    (*hTRate)()->GetXaxis()->SetTimeDisplay(true);
    (*hTRate)()->GetXaxis()->SetTimeOffset(0, "GMT");
    (*hTRate)()->GetXaxis()->SetTimeFormat("%Y/%m/%d");
    (*hTRate)()->GetZaxis()->SetTitle("#bar{p}/p");
    
    (*hTCrrR)()->GetXaxis()->SetTitle("Date");
    (*hTCrrR)()->GetXaxis()->SetTimeDisplay(true);
    (*hTCrrR)()->GetXaxis()->SetTimeOffset(0, "GMT");
    (*hTCrrR)()->GetXaxis()->SetTimeFormat("%Y/%m/%d");
    (*hTCrrR)()->GetZaxis()->SetTitle("#bar{p}/p");

    for (int ir = 11; ir <= 40; ++ir) {
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
        
        Legend leg_table("", PadWindow(0.58, 0.85, 0.65, 0.85));
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
        (*hCrrR)()->SetBinError  (ir, accp_func(rig) * std::sqrt(rlt1D.num_tmps().at(0)) / rlt1D.num_ref()); 
        
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
        (*hsmp)()->GetXaxis()->SetTitle("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]");
        (*hsmp)()->GetYaxis()->SetTitle("Events/Bin");
        (*hsmp)()->SetMaximum( 1.3 * (*hsmp)()->GetMaximum() );

        (*hsmp)()->Draw("pe");
        (*hsum)()->Draw("hist same");
        (*hbkg)()->Draw("hist same");
        (*hsig)()->Draw("hist same");
        //(*hbkg)()->Draw("pe same");
        //(*hsig)()->Draw("pe same");
        
        Legend leg_table2("", PadWindow(0.15, 0.42, 0.65, 0.85));
        leg_table2()->SetHeader(Form("Rigidity %.2f - %.2f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
        leg_table2()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table2()->AddEntry((*hsum)(), "Sum", "lp");
        leg_table2()->AddEntry((*hbkg)(), Form("e^{-}+#pi^{-}  (%.1f #pm %.1f)", rlt1D.num_tmps().at(1), rlt1D.num_errs().at(1)), "lp");
        leg_table2()->AddEntry((*hsig)(), Form("p^{-}  (%.1f #pm %.1f)", rlt1D.num_tmps().at(0), rlt1D.num_errs().at(0)), "lp");
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

    for (int it = 1; it <= AXtme.nbin(); ++it) {
        double tme = AXtme.center(it, AxisScale::kLinear);

        std::vector<Hist*> vhRllr_pr = Hist::ProjectAll(HistProj::kY, vhTllr_pr.at(it));
        std::vector<Hist*> vhRllr_el = Hist::ProjectAll(HistProj::kY, vhTllr_el.at(it));
        std::vector<Hist*> vhRllr_pos = Hist::ProjectAll(HistProj::kY, vhTllr_pos.at(it));
        std::vector<Hist*> vhRllr_neg = Hist::ProjectAll(HistProj::kY, vhTllr_neg.at(it));

        for (int ir = 10; ir <= 20; ++ir) {
            double rig = AXrig2.center(ir, AxisScale::kLog);
       
            DataFit1D fit1D(Form("T%3dR%3d", it, ir), { (TH1D*)((*vhRllr_pr.at(ir))()), (TH1D*)((*vhRllr_el.at(ir))()) }, (TH1D*)((*vhRllr_neg.at(ir))()), (TH1D*)((*vhRllr_pos.at(ir))()), true);
            ResultFit1D rlt1D = fit1D.result();
            if (rlt1D.ndof() == 0) continue; 
            
            (*hTPcnt)()->SetBinContent(it, ir, rlt1D.num_ref()); 
            (*hTNcnt)()->SetBinContent(it, ir, rlt1D.num_tmps().at(0)); 
            
            (*hTStat)()->SetBinContent(it, ir, rlt1D.num_errs().at(0) / rlt1D.num_ref()); 
            (*hTSyst)()->SetBinContent(it, ir, rlt1D.num_errs().at(0) / rlt1D.num_ref());

            (*hTRate)()->SetBinContent(it, ir, rlt1D.num_tmps().at(0) / rlt1D.num_ref()); 
            (*hTRate)()->SetBinError  (it, ir, rlt1D.num_errs().at(0) / rlt1D.num_ref());
            
            (*hTCrrR)()->SetBinContent(it, ir, accp_func(rig) * rlt1D.num_tmps().at(0) / rlt1D.num_ref()); 
            (*hTCrrR)()->SetBinError  (it, ir, accp_func(rig) * rlt1D.num_errs().at(0) / rlt1D.num_ref());

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
            (*hsmp)()->GetXaxis()->SetTitle("Mass^{2}/Z^{2} [(GV/c^{2})^{2}]");
            (*hsmp)()->GetYaxis()->SetTitle("Events/Bin");
            //(*hsmp)()->SetMaximum( 1.3 * (*hsmp)()->GetMaximum() );

            (*hsmp)()->Draw("pe");
            (*hsum)()->Draw("hist same");
            (*hbkg)()->Draw("hist same");
            (*hsig)()->Draw("hist same");
            
            Legend leg_table2("", PadWindow(0.15, 0.42, 0.65, 0.85));
            leg_table2()->SetHeader(Form("Time %d Rigidity %.2f - %.2f [GV/c]", it, AXrig2()(ir-1), AXrig2()(ir)));
            leg_table2()->AddEntry((*hsmp)(), "Data", "lp");
            leg_table2()->AddEntry((*hsum)(), "Sum", "lp");
            leg_table2()->AddEntry((*hbkg)(), Form("e^{-}+#pi^{-}  (%.1f #pm %.1f)", rlt1D.num_tmps().at(0), rlt1D.num_errs().at(0)), "lp");
            leg_table2()->AddEntry((*hsig)(), Form("p^{-}  (%.1f #pm %.1f)", rlt1D.num_tmps().at(1), rlt1D.num_errs().at(1)), "lp");
            leg_table2()->SetFillColor(0);
            leg_table2.draw();
            
            editor.save();

            std::cerr << Form("Rigidity %.2f - %.2f [GV/c]\n", AXrig2()(ir-1), AXrig2()(ir));
            std::cerr << Form("SMP %14.8f SIG %14.8f (%14.8f %14.8f) BKG %14.8f (%14.8f %14.8f) RATE %14.8f\n", rlt1D.num_smp(), rlt1D.num_tmps().at(0), rlt1D.wgt_tmps().at(0), rlt1D.wgt_errs().at(0), rlt1D.num_tmps().at(1), rlt1D.wgt_tmps().at(1), rlt1D.wgt_errs().at(1), rlt1D.num_tmps().at(0)/rlt1D.num_ref());
            std::cerr << Form("NCHI %14.8f\n", rlt1D.nchi());
            std::cerr << "\n";
        }
    }

    editor.close();
    
    TFile * ofle = new TFile("out/antip_i.root", "RECREATE");
    ofle->cd();

    haccp->Write();
    hTaccp->Write();

    (*hPcnt)()->Write();
    (*hNcnt)()->Write();
    (*hStat)()->Write();
    (*hSyst)()->Write();
    (*hRate)()->Write();
    (*hCrrR)()->Write();

    (*hTPcnt)()->Write();
    (*hTNcnt)()->Write();
    (*hTStat)()->Write();
    (*hTSyst)()->Write();
    (*hTRate)()->Write();
    (*hTCrrR)()->Write();

    ofle->Write();
    ofle->Close();

    return 1;
}
