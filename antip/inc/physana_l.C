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
    double crr = 1.21555e+00 + 8.49079e-01 * std::exp(-1.94265e+00 * arig) + 8.49079e-01 * std::exp(-1.94265e+00 * arig);
    return crr;
}

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    std::string subv = "12";

    UInt_t cntev = 0;
    
    UInt_t cntpr = 0;
    TFile* fmcpr = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcpr%s/YiMdst.root", subv.c_str()));
    TH1D*  hmcpr = (TH1D*)fmcpr->Get("hL_MC_cnt");
    TH1D*  hTmcpr = (TH1D*)fmcpr->Get("hTL_MC_cnt");
    TTree* tmcpr = (TTree*)fmcpr->Get("ana");
    tmcpr->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcpr->GetEntries(); ++it) { tmcpr->GetEntry(it); cntpr+=cntev; }
    
    UInt_t cntap = 0;
    TFile* fmcap = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcap%s/YiMdst.root", subv.c_str()));
    TH1D*  hmcap = (TH1D*)fmcap->Get("hL_MC_cnt");
    TH1D*  hTmcap = (TH1D*)fmcap->Get("hTL_MC_cnt");
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

    PdfEditor editor(Window(), "antip_l", "out");

    const Axis& AXrig = Hist::Head("hLP_sqrm")->xaxis();
    const Axis& AXrig2 = Hist::Head("hTLP_sqrm")->yaxis();
    const Axis& AXtme = Hist::Head("hTLP_sqrm")->xaxis();
 
    std::vector<THStack*> vhstack;
    std::vector<Hist*> vhist;

    std::vector<Hist*> vhsqrm_pr  = Hist::ProjectAll(HistProj::kY, Hist::Head("hLP_sqrm_pr"));
    std::vector<Hist*> vhsqrm_ep  = Hist::ProjectAll(HistProj::kY, Hist::Head("hLN_sqrm_el"));
    std::vector<Hist*> vhsqrm_pos = Hist::ProjectAll(HistProj::kY, Hist::Head("hLP_sqrm"));
    std::vector<Hist*> vhsqrm_neg = Hist::ProjectAll(HistProj::kY, Hist::Head("hLN_sqrm"));
    
    std::vector<Hist*> vhTsqrm_pr  = Hist::ProjectAll(HistProj::kZY, Hist::Head("hTLP_sqrm_pr"));
    std::vector<Hist*> vhTsqrm_ep  = Hist::ProjectAll(HistProj::kZY, Hist::Head("hTLN_sqrm_el"));
    std::vector<Hist*> vhTsqrm_pos = Hist::ProjectAll(HistProj::kZY, Hist::Head("hTLP_sqrm"));
    std::vector<Hist*> vhTsqrm_neg = Hist::ProjectAll(HistProj::kZY, Hist::Head("hTLN_sqrm"));
        
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

    for (int ir = 2; ir <= 15; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        bool has_ep_temp = (((TH1D*)((*vhsqrm_ep.at(ir))()))->GetEntries() > 0);
        std::vector<TH1D*> temps = has_ep_temp ? std::vector<TH1D*>({ (TH1D*)((*vhsqrm_pr.at(ir))()), (TH1D*)((*vhsqrm_ep.at(ir))()) }) : std::vector<TH1D*>({ (TH1D*)((*vhsqrm_pr.at(ir))()) });
        
        DataFit1D fit1D(Form("R%3d", ir), temps, (TH1D*)((*vhsqrm_neg.at(ir))()), (TH1D*)((*vhsqrm_pos.at(ir))()), true, "Mass^{2}/Z^{2} [(GV/c^{2})^{2}]", "Events/Bin");
        //DataFit1D fit1D(Form("R%3d", ir), has_ep_temp ? { (TH1D*)((*vhsqrm_pr.at(ir))()), (TH1D*)((*vhsqrm_ep.at(ir))()) } : { (TH1D*)((*vhsqrm_pr.at(ir))()) }, (TH1D*)((*vhsqrm_neg.at(ir))()), (TH1D*)((*vhsqrm_pos.at(ir))()), true, "Mass^{2}/Z^{2} [(GV/c^{2})^{2}]", "Events/Bin");
        ResultFit1D rlt1D = fit1D.result();
        if (rlt1D.ndof() == 0) continue;
        
        (*hPcnt)()->SetBinContent(ir, rlt1D.num_ref()); 
        (*hNcnt)()->SetBinContent(ir, rlt1D.num_tmps().at(0)); 
        
        (*hStat)()->SetBinContent(ir, rlt1D.num_errs().at(0) / rlt1D.num_ref()); 
        (*hSyst)()->SetBinContent(ir, rlt1D.num_errs().at(0) / rlt1D.num_ref());

        (*hRate)()->SetBinContent(ir, rlt1D.num_tmps().at(0) / rlt1D.num_ref()); 
        (*hRate)()->SetBinError  (ir, rlt1D.num_errs().at(0) / rlt1D.num_ref()); 
        
        (*hCrrR)()->SetBinContent(ir, accp_func(rig) * rlt1D.num_tmps().at(0) / rlt1D.num_ref()); 
        (*hCrrR)()->SetBinError  (ir, accp_func(rig) * rlt1D.num_errs().at(0) / rlt1D.num_ref()); 
        
        // TEST
        std::cerr << Form("RATE (FIT/STAT) %14.8f\n", rlt1D.num_errs().at(0) / std::sqrt(rlt1D.num_tmps().at(0)) );

        Hist* hsmp = Hist::New(rlt1D.hsmp().get());
        Hist* hsum = Hist::New(rlt1D.hsum().get());
        Hist* hsig = Hist::New(rlt1D.htmps().at(0).get());
        Hist* hbkg = has_ep_temp ? Hist::New(rlt1D.htmps().at(1).get()) : nullptr;
        Hist* href = Hist::New(rlt1D.href().get());
		
        hsmp->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
		hsum->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
		if (has_ep_temp) hbkg->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
		hsig->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        
        editor.create();
       
        (*hsmp)()->GetXaxis()->CenterTitle();
        (*hsmp)()->SetMaximum( 1.3 * (*hsmp)()->GetMaximum() );

        (*hsmp)()->Draw("pe");
        (*hsum)()->Draw("hist same");
        if (has_ep_temp) (*hbkg)()->Draw("hist same");
        (*hsig)()->Draw("hist same");
        
        Legend leg_table2("", PadWindow(0.15, 0.42, 0.65, 0.85));
        leg_table2()->SetHeader(Form("Rigidity %.2f - %.2f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
        leg_table2()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table2()->AddEntry((*hsum)(), "Sum", "lp");
        if (has_ep_temp) leg_table2()->AddEntry((*hbkg)(), Form("e^{-}+#pi^{-}  (%.1f #pm %.1f)", rlt1D.num_tmps().at(1), rlt1D.num_errs().at(1)), "lp");
        leg_table2()->AddEntry((*hsig)(), Form("p^{-}  (%.1f #pm %.1f)", rlt1D.num_tmps().at(0), rlt1D.num_errs().at(0)), "lp");
        leg_table2()->SetFillColor(0);
        leg_table2.draw();
        
        editor.save();

        std::cerr << Form("Rigidity %.2f - %.2f [GV/c]\n", AXrig()(ir-1), AXrig()(ir));
        if (has_ep_temp) std::cerr << Form("SMP %14.8f SIG %14.8f (%14.8f %14.8f) BKG %14.8f (%14.8f %14.8f) RATE %14.8f\n", rlt1D.num_smp(), rlt1D.num_tmps().at(0), rlt1D.wgt_tmps().at(0), rlt1D.wgt_errs().at(0), rlt1D.num_tmps().at(1), rlt1D.wgt_tmps().at(1), rlt1D.wgt_errs().at(1), rlt1D.num_tmps().at(0)/rlt1D.num_ref());
        else std::cerr << Form("SMP %14.8f SIG %14.8f (%14.8f %14.8f) BKG %14.8f (%14.8f %14.8f) RATE %14.8f\n", rlt1D.num_smp(), rlt1D.num_tmps().at(0), rlt1D.wgt_tmps().at(0), rlt1D.wgt_errs().at(0), 0.0, 0.0, 0.0, rlt1D.num_tmps().at(0)/rlt1D.num_ref());
        std::cerr << Form("NCHI %14.8f\n", rlt1D.nchi());
        std::cerr << "\n";
    }
  
    for (int it = 1; it <= AXtme.nbin(); ++it) {
        double tme = AXtme.center(it, AxisScale::kLinear);
    
        std::vector<Hist*> vhRsqrm_pr  = Hist::ProjectAll(HistProj::kY, vhTsqrm_pr.at(it));
        std::vector<Hist*> vhRsqrm_ep  = Hist::ProjectAll(HistProj::kY, vhTsqrm_ep.at(it));
        std::vector<Hist*> vhRsqrm_pos = Hist::ProjectAll(HistProj::kY, vhTsqrm_pos.at(it));
        std::vector<Hist*> vhRsqrm_neg = Hist::ProjectAll(HistProj::kY, vhTsqrm_neg.at(it));

        for (int ir = 1; ir <= 5; ++ir) {
            double rig = AXrig2.center(ir, AxisScale::kLog);

/*            
            Fit::RooVar roovar("sqrm", vhRsqrm_neg.at(ir), HistList({ vhRsqrm_pr.at(ir), vhRsqrm_ep.at(ir) }));
            Fit::RooSysResult rlt(roovar, true, 100);
            Fit::RooVar var = rlt.var();
            Fit::RooPar stdpar = rlt.std_par();
            Fit::RooPar syspar = rlt.sys_par();

            std::cerr << Form("ROOFIT %14.8f ERR %14.8f %14.8f\n", stdpar.val(0), stdpar.err(0), syspar.err(0));
  */      
            bool has_ep_temp = (((TH1D*)((*vhRsqrm_ep.at(ir))()))->Integral() > 0);
            std::vector<TH1D*> temps = has_ep_temp ? std::vector<TH1D*>({ (TH1D*)((*vhRsqrm_pr.at(ir))()), (TH1D*)((*vhRsqrm_ep.at(ir))()) }) : std::vector<TH1D*>({ (TH1D*)((*vhRsqrm_pr.at(ir))()) });

            DataFit1D fit1D(Form("T%3dR%3d", it, ir), temps, (TH1D*)((*vhRsqrm_neg.at(ir))()), (TH1D*)((*vhRsqrm_pos.at(ir))()), true, "Mass^{2}/Z^{2} [(GV/c^{2})^{2}]", "Events/Bin");
            ResultFit1D rlt1D = fit1D.result();
            if (rlt1D.ndof() == 0) continue; 
           
            (*hTPcnt)()->SetBinContent(it, ir, rlt1D.num_ref()); 
            (*hTNcnt)()->SetBinContent(it, ir, rlt1D.num_tmps().at(0)); 
            
            (*hTStat)()->SetBinContent(it, ir, rlt1D.num_errs().at(0) / rlt1D.num_ref()); 
            (*hTSyst)()->SetBinContent(it, ir, rlt1D.num_errs().at(0) / rlt1D.num_ref());

            (*hTRate)()->SetBinContent(it, ir, rlt1D.num_tmps().at(0) / rlt1D.num_ref()); 
            (*hTRate)()->SetBinError  (it, ir, rlt1D.num_errs().at(0) / rlt1D.num_ref());
            
            (*hTCrrR)()->SetBinContent(it, ir, accp_func(rig) * rlt1D.num_tmps().at(0) / rlt1D.num_ref()); 
            //(*hTCrrR)()->SetBinError  (it, ir, accp_func(rig) * rlt1D.num_errs().at(0) / rlt1D.num_ref());
            (*hTCrrR)()->SetBinError  (it, ir, accp_func(rig) * std::sqrt(rlt1D.num_tmps().at(0)) / rlt1D.num_ref());

            Hist* hsmp = Hist::New(rlt1D.hsmp().get());
            Hist* hsum = Hist::New(rlt1D.hsum().get());
            Hist* hsig = Hist::New(rlt1D.htmps().at(0).get());
            Hist* hbkg = (has_ep_temp ? Hist::New(rlt1D.htmps().at(1).get()) : nullptr);
            Hist* href = Hist::New(rlt1D.href().get());
		    
            hsmp->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
		    hsum->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
		    if (has_ep_temp) hbkg->style(Line(kBlue, 0, 2), Marker(kBlue, MarkerStyle(MarkerShape::kCircle )));
		    hsig->style(Line(kRed    , 0, 2), Marker(kRed, MarkerStyle(MarkerShape::kCircle )));
            
            //editor.create();
            //(*hsmp)()->Draw("hist");
            //editor.save();
            //editor.create();
            //(*href)()->Draw("hist");
            //editor.save();

            //std::cerr << Form("SMP %14.8f REF %14.8f\n", rlt1D.hsmp()->Integral(), rlt1D.href()->Integral());
            //std::cerr << Form("SIG %14.8f BKG %14.8f\n", (*hsig)()->Integral(), (*hbkg)()->Integral());
            //if (rlt1D.hsmp()->Integral() == 0) {
            //    for (auto&& elem : rlt1D.data_smp()) {
            //        std::cerr << Form("%14.8f %14.8f\n", elem[0], elem[1]);
            //    }
            //}

            editor.create();
       
            (*hsmp)()->GetXaxis()->CenterTitle();
            //(*hsmp)()->SetMaximum( 1.3 * (*hsmp)()->GetMaximum() );
            (*hsmp)()->Draw("pe");
            (*hsum)()->Draw("hist same");
            if (has_ep_temp) (*hbkg)()->Draw("hist same");
            (*hsig)()->Draw("hist same");
            
            Legend leg_table2("", PadWindow(0.15, 0.42, 0.65, 0.85));
            leg_table2()->SetHeader(Form("Time %d Rigidity %.2f - %.2f [GV/c]", it, AXrig2()(ir-1), AXrig2()(ir)));
            leg_table2()->AddEntry((*hsmp)(), "Data", "lp");
            leg_table2()->AddEntry((*hsum)(), "Sum", "lp");
            if (has_ep_temp) leg_table2()->AddEntry((*hbkg)(), Form("e^{-}+#pi^{-}  (%.1f #pm %.1f)", rlt1D.num_tmps().at(1), rlt1D.num_errs().at(1)), "lp");
            leg_table2()->AddEntry((*hsig)(), Form("p^{-}  (%.1f #pm %.1f)", rlt1D.num_tmps().at(0), rlt1D.num_errs().at(0)), "lp");
            leg_table2()->SetFillColor(0);
            leg_table2.draw();
            
            editor.save();

            std::cerr << Form("Time %f %f Rigidity %.2f - %.2f [GV/c]\n", AXtme()(it-1)/10000, AXtme()(it)/100000, AXrig2()(ir-1), AXrig2()(ir));
            if (has_ep_temp) std::cerr << Form("SMP %14.8f SIG %14.8f (%14.8f %14.8f) BKG %14.8f (%14.8f %14.8f) RATE %14.8f\n", rlt1D.num_smp(), rlt1D.num_tmps().at(0), rlt1D.wgt_tmps().at(0), rlt1D.wgt_errs().at(0), rlt1D.num_tmps().at(1), rlt1D.wgt_tmps().at(1), rlt1D.wgt_errs().at(1), rlt1D.num_tmps().at(0)/rlt1D.num_ref());
            else std::cerr << Form("SMP %14.8f SIG %14.8f (%14.8f %14.8f) BKG %14.8f (%14.8f %14.8f) RATE %14.8f\n", rlt1D.num_smp(), rlt1D.num_tmps().at(0), rlt1D.wgt_tmps().at(0), rlt1D.wgt_errs().at(0), 0.0, 0.0, 0.0, rlt1D.num_tmps().at(0)/rlt1D.num_ref());
            std::cerr << Form("NCHI %14.8f\n", rlt1D.nchi());
            std::cerr << "\n";
        }
    }

    editor.close();

    TFile * ofle = new TFile("out/antip_l.root", "RECREATE");
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
