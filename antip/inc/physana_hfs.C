#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "DataFit2D.h"
#include "DataFit2D.C"

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
    double crr = 1.05139e+00 + 2.18257e+00 * std::exp(-2.06862e+00 * arig) + 2.00462e-01 * std::exp(-1.52078e-01 * arig);
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
    TH1D*  hmcpr = (TH1D*)fmcpr->Get("hHfs_MC_cnt");
    TTree* tmcpr = (TTree*)fmcpr->Get("ana");
    tmcpr->SetBranchAddress("event", &cntev);
    for (int it = 0; it < tmcpr->GetEntries(); ++it) { tmcpr->GetEntry(it); cntpr+=cntev; }
    
    UInt_t cntap = 0;
    TFile* fmcap = TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcap%s/YiMdst.root", subv.c_str()));
    TH1D*  hmcap = (TH1D*)fmcap->Get("hHfs_MC_cnt");
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
    
    Hist* hcc_cc = Hist::New("hHNfs_lchiy_lrvar_cc", (TH1*)TFile::Open(Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/mcprflux%s/YiMdst.root", subv.c_str()))->Get("hHNfs_lchiy_lrvar"));

    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/antip/19Nov03/iss%s", subv.c_str()));

    PdfEditor editor(Window(), "antip_hfs", "out");

    const Axis& AXrig = Hist::Head("hHPfs_lchiy_lrvar")->xaxis();
 
    std::vector<THStack*> vhstack;
    std::vector<Hist*> vhist;

    std::vector<Hist*> vhcc_cc = Hist::ProjectAll(HistProj::kZY, hcc_cc);
    std::vector<Hist*> vhcc_pos = Hist::ProjectAll(HistProj::kZY, Hist::Head("hHPfs_lchiy_lrvar"));
    std::vector<Hist*> vhcc_neg = Hist::ProjectAll(HistProj::kZY, Hist::Head("hHNfs_lchiy_lrvar"));

    Hist* hPcnt = Hist::New("hPcnt", HistAxis(AXrig));
    Hist* hNcnt = Hist::New("hNcnt", HistAxis(AXrig));
    Hist* hStat = Hist::New("hStat", HistAxis(AXrig));
    Hist* hSyst = Hist::New("hSyst", HistAxis(AXrig));
    Hist* hRate = Hist::New("hRate", HistAxis(AXrig));
    Hist* hCrrR = Hist::New("hCrrR", HistAxis(AXrig));
    
    for (int ir = 30; ir <= AXrig.nbin()-4; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        
        DataFit2D fit2D(Form("R%3d", ir), { (TH2D*)((*vhcc_pos.at(ir))()), (TH2D*)((*vhcc_cc.at(ir))()) }, (TH2D*)((*vhcc_neg.at(ir))()), (TH2D*)((*vhcc_pos.at(ir))()), true, "CC1", "CC2", "Events/Bin");
        ResultFit2D rlt2D = fit2D.result();
        if (rlt2D.ndof() == 0) continue;
        
        (*hPcnt)()->SetBinContent(ir, rlt2D.num_ref()); 
        (*hNcnt)()->SetBinContent(ir, rlt2D.num_tmps().at(0)); 
        
        (*hStat)()->SetBinContent(ir, rlt2D.num_errs().at(0) / rlt2D.num_ref()); 
        (*hSyst)()->SetBinContent(ir, rlt2D.num_errs().at(0) / rlt2D.num_ref());

        (*hRate)()->SetBinContent(ir, rlt2D.num_tmps().at(0) / rlt2D.num_ref()); 
        (*hRate)()->SetBinError  (ir, rlt2D.num_errs().at(0) / rlt2D.num_ref()); 
        
        (*hCrrR)()->SetBinContent(ir, accp_func(rig) * rlt2D.num_tmps().at(0) / rlt2D.num_ref()); 
        //(*hCrrR)()->SetBinError  (ir, accp_func(rig) * rlt2D.num_errs().at(0) / rlt2D.num_ref()); 
        (*hCrrR)()->SetBinError  (ir, accp_func(rig) * std::sqrt(rlt2D.num_tmps().at(0)) / rlt2D.num_ref()); 
       
        // TEST
        std::cerr << Form("RATE (FIT/STAT) %14.8f\n", rlt2D.num_errs().at(0) / std::sqrt(rlt2D.num_tmps().at(0)) );


        Hist* hsmpx = Hist::New(rlt2D.hsmp_x().get());
        Hist* hsumx = Hist::New(rlt2D.hsum_x().get());
        Hist* hsigx = Hist::New(rlt2D.htmps_x().at(0).get());
        Hist* hbkgx = Hist::New(rlt2D.htmps_x().at(1).get());
        Hist* hrefx = Hist::New(rlt2D.href_x().get());
		
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
        
        Legend leg_tablex("", PadWindow(0.15, 0.42, 0.65, 0.85));
        leg_tablex()->SetHeader(Form("Rigidity %.2f - %.2f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
        leg_tablex()->AddEntry((*hsmpx)(), "Data", "lp");
        leg_tablex()->AddEntry((*hsumx)(), "Sum", "lp");
        leg_tablex()->AddEntry((*hbkgx)(), Form("p^{+}  (%.1f #pm %.1f)", rlt2D.num_tmps().at(1), rlt2D.num_errs().at(1)), "lp");
        leg_tablex()->AddEntry((*hsigx)(), Form("p^{-}  (%.1f #pm %.1f)", rlt2D.num_tmps().at(0), rlt2D.num_errs().at(0)), "lp");
        leg_tablex()->SetFillColor(0);
        leg_tablex.draw();
        
        editor.save();
      
        Hist* hsmpy = Hist::New(rlt2D.hsmp_y().get());
        Hist* hsumy = Hist::New(rlt2D.hsum_y().get());
        Hist* hsigy = Hist::New(rlt2D.htmps_y().at(0).get());
        Hist* hbkgy = Hist::New(rlt2D.htmps_y().at(1).get());
        Hist* hrefy = Hist::New(rlt2D.href_y().get());
		
        hsmpy->style(Line(kBlack  , 0, 2), Marker(kBlack  , MarkerStyle(MarkerShape::kCircle )));
		hsumy->style(Line(kGreen+2, 0, 2), Marker(kGreen+2, MarkerStyle(MarkerShape::kCircle )));
		hbkgy->style(Line(kBlue   , 0, 2), Marker(kBlue   , MarkerStyle(MarkerShape::kCircle )));
		hsigy->style(Line(kRed    , 0, 2), Marker(kRed    , MarkerStyle(MarkerShape::kCircle )));
        
        editor.create();
       
        (*hsmpy)()->GetXaxis()->CenterTitle();
        (*hsmpy)()->SetMaximum( 1.3 * (*hsmpy)()->GetMaximum() );

        (*hsmpy)()->Draw("pe");
        (*hsumy)()->Draw("hist same");
        (*hbkgy)()->Draw("hist same");
        (*hsigy)()->Draw("hist same");
        
        Legend leg_tabley("", PadWindow(0.15, 0.42, 0.65, 0.85));
        leg_tabley()->SetHeader(Form("Rigidity %.2f - %.2f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
        leg_tabley()->AddEntry((*hsmpy)(), "Data", "lp");
        leg_tabley()->AddEntry((*hsumy)(), "Sum", "lp");
        leg_tabley()->AddEntry((*hbkgy)(), Form("p^{+}  (%.1f #pm %.1f)", rlt2D.num_tmps().at(1), rlt2D.num_errs().at(1)), "lp");
        leg_tabley()->AddEntry((*hsigy)(), Form("p^{-}  (%.1f #pm %.1f)", rlt2D.num_tmps().at(0), rlt2D.num_errs().at(0)), "lp");
        leg_tabley()->SetFillColor(0);
        leg_tabley.draw();
        
        editor.save();

        std::cerr << Form("Rigidity %.2f - %.2f [GV/c]\n", AXrig()(ir-1), AXrig()(ir));
        std::cerr << Form("SMP %14.8f TEMPS %14.8f %14.8f REF %14.8f\n", fit2D.data_num_smp(), fit2D.data_num_tmp(0), fit2D.data_num_tmp(1), fit2D.data_num_ref());
        std::cerr << Form("SMP %14.8f SIG %14.8f (%14.8f %14.8f) BKG %14.8f (%14.8f %14.8f) RATE %14.8f\n", rlt2D.num_smp(), rlt2D.num_tmps().at(0), rlt2D.wgt_tmps().at(0), rlt2D.wgt_errs().at(0), rlt2D.num_tmps().at(1), rlt2D.wgt_tmps().at(1), rlt2D.wgt_errs().at(1), rlt2D.num_tmps().at(0)/rlt2D.num_ref());
        std::cerr << Form("NCHI %14.8f\n", rlt2D.nchi());
        std::cerr << "\n";
    }

    editor.close();
    
    TFile * ofle = new TFile("out/antip_hfs.root", "RECREATE");
    ofle->cd();

    haccp->Write();

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
