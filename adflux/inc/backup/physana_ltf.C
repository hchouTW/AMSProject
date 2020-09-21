#include <CPPLibs.h>
#include <ROOTLibs.h>

#include "DataFit1D.h"
#include "DataFit1D.C"
#include <TrSys.h>

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
    double crr = 1.04370e+00 + 
                 1.15612e+00 * std::exp(-2.04173e+00 * arig) +
                 2.39598e-01 * std::exp(-2.09577e-01 * arig) +
                 3.34768e-02 * std::exp(-1.29122e-02 * arig);
    return crr;
}

int main(int argc, char* argv[]) {
    using namespace MGROOT;
    MGROOT::LoadDefaultEnvironment();
    Hist::AddDirectory(0);
    std::string subv = "59";

    Hist::Load("YiMdst.root", Form("/eos/ams/user/h/hchou/AMSData/subj/apflux/20Jan15/iss%s", subv.c_str()));

    PdfEditor editor(Window(), "apflux_ltf", "out");

    const Axis& AXrig = Hist::Head("hLP_sqrm")->xaxis();
    std::vector<Hist*> vhsqrm_pos = Hist::ProjectAll(HistProj::kY, Hist::Head("hLP_sqrm"));
    std::vector<Hist*> vhsqrm_neg = Hist::ProjectAll(HistProj::kY, Hist::Head("hLN_sqrm"));
    std::vector<Hist*> vhsqrm_pr  = Hist::ProjectAll(HistProj::kY, Hist::Head("hLP_sqrm_pr"));
    std::vector<Hist*> vhsqrm_ep  = Hist::ProjectAll(HistProj::kY, Hist::Head("hLN_sqrm_el"));
    
    Hist* hPcnt = Hist::New("hPcnt", HistAxis(AXrig));
    Hist* hNcnt = Hist::New("hNcnt", HistAxis(AXrig));
    Hist* hStat = Hist::New("hStat", HistAxis(AXrig));
    Hist* hSyst = Hist::New("hSyst", HistAxis(AXrig));
    Hist* hRate = Hist::New("hRate", HistAxis(AXrig));
    Hist* hCrrR = Hist::New("hCrrR", HistAxis(AXrig));

    const Axis& AXtme  = Hist::Head("hTLP_sqrm")->xaxis();
    const Axis& AXTrig = Hist::Head("hTLP_sqrm")->yaxis();
    std::vector<Hist*> vhTsqrm_pos = Hist::ProjectAll(HistProj::kZY, Hist::Head("hTLP_sqrm"));
    std::vector<Hist*> vhTsqrm_neg = Hist::ProjectAll(HistProj::kZY, Hist::Head("hTLN_sqrm"));
    std::vector<Hist*> vhTsqrm_pr  = Hist::ProjectAll(HistProj::kZY, Hist::Head("hTLP_sqrm_pr"));
    std::vector<Hist*> vhTsqrm_ep  = Hist::ProjectAll(HistProj::kZY, Hist::Head("hTLN_sqrm_el"));
    
    std::vector<Hist*> vhTRsqrm_pr = Hist::ProjectAll(HistProj::kY, Hist::Head("hTRLP_sqrm_pr"));
    std::vector<Hist*> vhTRsqrm_ep = Hist::ProjectAll(HistProj::kY, Hist::Head("hTRLN_sqrm_el"));
        
    Hist* hTPcnt = Hist::New("hTPcnt", HistAxis(AXtme, AXTrig));
    Hist* hTNcnt = Hist::New("hTNcnt", HistAxis(AXtme, AXTrig));
    Hist* hTStat = Hist::New("hTStat", HistAxis(AXtme, AXTrig));
    Hist* hTSyst = Hist::New("hTSyst", HistAxis(AXtme, AXTrig));
    Hist* hTRate = Hist::New("hTRate", HistAxis(AXtme, AXTrig));
    Hist* hTCrrR = Hist::New("hTCrrR", HistAxis(AXtme, AXTrig));

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

    /////////////////////////// TEST ///////////////////////////////////////
    HistFit::Axis1D AX1Dsqrm(
        "Mass^{2}/Z^{2} [(GV/c^{2})^{2}]",
        "Events/Bin",
        (*Hist::Head("hLP_sqrm"))()->GetYaxis());

    for (int ir = 1; ir <= 18; ++ir) {
        double rig = AXrig.center(ir, AxisScale::kLog);
        TH1D* smp = (TH1D*)((*vhsqrm_neg.at(ir))());
        std::vector<TH1D*> tmps = std::vector<TH1D*>({ (TH1D*)((*vhsqrm_pr.at(ir))()), (TH1D*)((*vhsqrm_ep.at(ir))()) });
        
        HistFit::Hist1D hist_pr(Form("BIN%03d_PR", ir), "", (TH1D*)((*vhsqrm_pos.at(ir))()), AX1Dsqrm);
        HistFit::HistFit1D fit1D(smp, tmps, AX1Dsqrm, Form("BIN%03d_", ir));
        if (!fit1D.status()) continue;

        std::cerr << Form("STATUS %d BIN %3d NSMP %8.1f WGT %8.1f %8.1f ERR(%14.8f %14.8f) FUNC(%14.8f %14.8f) NCHI %14.8f APFLUX (%14.8f)\n",
            fit1D.status(),
            ir,
            fit1D.nsmp(),
            fit1D.wgts(0),
            fit1D.wgts(1),
            fit1D.errs(0) / fit1D.wgts(0),
            fit1D.errs(1) / fit1D.wgts(1),
            fit1D.fluc(0).err / fit1D.wgts(0),
            fit1D.fluc(1).err / fit1D.wgts(1),
            fit1D.nchi(),
            1.0e+4 * fit1D.wgts(0) / hist_pr.data().sum());
        
        
        (*hPcnt)()->SetBinContent(ir, hist_pr.data().sum()); 
        (*hNcnt)()->SetBinContent(ir, fit1D.wgts(0)); 
        
        (*hStat)()->SetBinContent(ir, fit1D.errs(0) / hist_pr.data().sum()); 
        (*hSyst)()->SetBinContent(ir, fit1D.errs(0) / hist_pr.data().sum());

        (*hRate)()->SetBinContent(ir, fit1D.wgts(0) / hist_pr.data().sum()); 
        (*hRate)()->SetBinError  (ir, fit1D.wgts(0) / hist_pr.data().sum()); 
        
        (*hCrrR)()->SetBinContent(ir, accp_func(rig) * fit1D.wgts(0) / hist_pr.data().sum()); 
        //(*hCrrR)()->SetBinError  (ir, accp_func(rig) * rlt1D.num_errs().at(0) / rlt1D.num_ref()); 
        (*hCrrR)()->SetBinError  (ir, (*hCrrR)()->GetBinContent(ir) * fit1D.errs(0) / fit1D.wgts(0)); 
        
        Hist* hsmp = Hist::New(fit1D.ref_smp().get());
        Hist* hsum = Hist::New(fit1D.sum_tmps().get());
        Hist* hsig = Hist::New(fit1D.wgt_tmps(0).get());
        Hist* hbkg = Hist::New(fit1D.wgt_tmps(1).get());
        
        editor.create();
       
        (*hsmp)()->GetXaxis()->CenterTitle();
        (*hsmp)()->SetMaximum( 1.3 * (*hsmp)()->GetMaximum() );

        (*hsmp)()->Draw("pe");
        (*hsum)()->Draw("hist same");
        (*hbkg)()->Draw("hist same");
        (*hsig)()->Draw("hist same");
        
        Legend leg_table2("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.65, 0.85));
        leg_table2()->SetHeader(Form("Rigidity %.2f - %.2f [GV/c]", AXrig()(ir-1), AXrig()(ir)));
        leg_table2()->AddEntry((*hsmp)(), "Data", "lp");
        leg_table2()->AddEntry((*hsum)(), "Sum", "lp");
        leg_table2()->AddEntry((*hbkg)(), Form("e^{-}+#pi^{-}  (%.2f #pm %.2f)", fit1D.wgts(1), fit1D.errs(1)), "lp");
        leg_table2()->AddEntry((*hsig)(), Form("p^{-}  (%.2f #pm %.2f)", fit1D.wgts(0), fit1D.errs(0)), "lp");
        leg_table2()->SetFillColor(0);
        leg_table2.draw();
        
        editor.save();
    }

    HistFit::Axis1D AX1DTsqrm(
        "Mass^{2}/Z^{2} [(GV/c^{2})^{2}]",
        "Events/Bin",
        (*Hist::Head("hTLP_sqrm"))()->GetZaxis());

    for (int it = 1; it <= AXtme.nbin(); ++it) {
        double tme = AXtme.center(it, AxisScale::kLinear);
    
        std::vector<Hist*> vhRsqrm_pr  = Hist::ProjectAll(HistProj::kY, vhTsqrm_pr.at(it));
        std::vector<Hist*> vhRsqrm_ep  = Hist::ProjectAll(HistProj::kY, vhTsqrm_ep.at(it));
        std::vector<Hist*> vhRsqrm_pos = Hist::ProjectAll(HistProj::kY, vhTsqrm_pos.at(it));
        std::vector<Hist*> vhRsqrm_neg = Hist::ProjectAll(HistProj::kY, vhTsqrm_neg.at(it));

        for (int ir = 1; ir <= 6; ++ir) {
            double rig = AXTrig.center(ir, AxisScale::kLog);

            HistFit::Hist1D href(Form("BIN%03dT%03d_PR", ir, it), "", (TH1D*)((*vhRsqrm_pos.at(ir))()), AX1DTsqrm);
            HistFit::HistFit1D fit1D(
                    (TH1D*)((*vhRsqrm_neg.at(ir))()), 
                    //{ (TH1D*)((*vhTRsqrm_pr.at(ir))()), (TH1D*)((*vhTRsqrm_ep.at(ir))()) },
                    { (TH1D*)((*vhRsqrm_pr.at(ir))()), (TH1D*)((*vhRsqrm_ep.at(ir))()) },
                    AX1DTsqrm, Form("BIN%03dT%03d_NEG", ir, it));
            
            if (!fit1D.status()) continue;

            std::cerr << Form("STATUS %d BIN %3d %3d NSMP %8.1f WGT %8.1f %8.1f ERR(%14.8f %14.8f) FUNC(%14.8f %14.8f) NCHI %14.8f APFLUX (%14.8f)\n",
                fit1D.status(),
                ir, it,
                fit1D.nsmp(),
                fit1D.wgts(0),
                fit1D.wgts(1),
                fit1D.errs(0) / fit1D.wgts(0),
                fit1D.errs(1) / fit1D.wgts(1),
                fit1D.fluc(0).err / fit1D.wgts(0),
                fit1D.fluc(1).err / fit1D.wgts(1),
                fit1D.nchi(),
                1.0e+4 * fit1D.wgts(0) / href.data().sum());
            
            (*hTPcnt)()->SetBinContent(it, ir, href.data().sum()); 
            (*hTNcnt)()->SetBinContent(it, ir, fit1D.wgts(0)); 
            
            (*hTStat)()->SetBinContent(it, ir, fit1D.errs(0) / href.data().sum()); 
            (*hTSyst)()->SetBinContent(it, ir, fit1D.errs(0) / href.data().sum());

            (*hTRate)()->SetBinContent(it, ir, fit1D.wgts(0) / href.data().sum()); 
            (*hTRate)()->SetBinError  (it, ir, fit1D.wgts(0) / href.data().sum()); 
            
            (*hTCrrR)()->SetBinContent(it, ir, accp_func(rig) * fit1D.wgts(0) / href.data().sum()); 
            (*hTCrrR)()->SetBinError  (it, ir, (*hTCrrR)()->GetBinContent(it, ir) * fit1D.errs(0) / fit1D.wgts(0)); 
            
            Hist* hsmp = Hist::New(fit1D.ref_smp().get());
            Hist* hsum = Hist::New(fit1D.sum_tmps().get());
            Hist* hsig = Hist::New(fit1D.wgt_tmps(0).get());
            Hist* hbkg = Hist::New(fit1D.wgt_tmps(1).get());
            
            editor.create();
            
            (*hsmp)()->GetXaxis()->CenterTitle();
            (*hsmp)()->SetMaximum( 1.3 * (*hsmp)()->GetMaximum() );

            (*hsmp)()->Draw("pe");
            (*hsum)()->Draw("hist same");
            (*hbkg)()->Draw("hist same");
            (*hsig)()->Draw("hist same");
            
            Legend leg_table2("", TextStyle(kBlack, 20, 43), PadWindow(0.15, 0.42, 0.65, 0.85));
            leg_table2()->SetHeader(Form("Rigidity %.2f - %.2f [GV/c]", AXTrig()(ir-1), AXTrig()(ir)));
            leg_table2()->AddEntry((*hsmp)(), "Data", "lp");
            leg_table2()->AddEntry((*hsum)(), "Sum", "lp");
            leg_table2()->AddEntry((*hbkg)(), Form("e^{-}+#pi^{-}  (%.2f #pm %.2f)", fit1D.wgts(1), fit1D.errs(1)), "lp");
            leg_table2()->AddEntry((*hsig)(), Form("p^{-}  (%.2f #pm %.2f)", fit1D.wgts(0), fit1D.errs(0)), "lp");
            leg_table2()->SetFillColor(0);
            leg_table2.draw();
            
            editor.save();
        }
    }
    
    editor.close();

    TFile * ofle = new TFile("out/apflux_ltf.root", "RECREATE");
    ofle->cd();

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
