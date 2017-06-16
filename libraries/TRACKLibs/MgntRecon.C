#ifndef __MgntRecon_C__
#define __MgntRecon_C__

#include "MgntRecon.h"

Bool_t Track::fit(Track::MethodOption method, PhyOption phyOpt) {
	fPhyOpt = phyOpt;
	Bool_t isSuccess = false; 
	switch (method) {
		case ALCARAZ : 
			{ isSuccess = officalFit(method, phyOpt); }
			break;
		case CHOUTKO : 
			{ isSuccess = officalFit(method, phyOpt); }
			break;
		case CHIKANIAN : 
			{ isSuccess = officalFit(method, phyOpt); }
			break;
		case ANALYTIC : 
			{ isSuccess = analyticFit(); }
			break;
		case SIMPLE : 
			{ isSuccess = (analyticFit()) ? (simpleFit()) : false; }
			break;
		case PHYSICAL : 
			{ isSuccess = (analyticFit()) ? ((simpleFit()) ? physicalFit() : false) : false; }
			break;
		default :
			break;
	}
	return isSuccess;
}

Bool_t Track::officalFit(Track::MethodOption method, PhyOption phyOpt) {
	Bool_t isSuccess = false; 
	switch (method) {
		case ALCARAZ : 
			break;
		case CHOUTKO : 
			break;
		case CHIKANIAN : 
			break;
		default :
			break;
	}
	return isSuccess;
}

Bool_t Track::analyticFit() {
  // linear fitting for x or y dir
	Double_t refZCoo = fHits.at(0).CZ();
	SVecD<2> resX;
	SVecD<2> resY;
	SMtxSymD<2> mtxX;
	SMtxSymD<2> mtxY;
	mtxX(0, 0) = ONE;
	mtxY(0, 0) = ONE;
	for (Int_t ih = 0; ih < fHits.size(); ++ih) {
		HitSt & hit = fHits.at(ih);
		Double_t dz  = hit.CZ() - refZCoo;
		Double_t dz2 = dz * dz; 
		if (hit.SideX()) {
			resX(0) += (hit.CX());
			resX(1) += (dz * hit.CX());
			mtxX(0, 0) += ONE;
			mtxX(0, 1) += (dz);
			mtxX(1, 1) += (dz2);
		}
		if (hit.SideY()) {
			resY(0) += (hit.CY());
			resY(1) += (dz * hit.CY());
			mtxY(0, 0) += ONE;
			mtxY(0, 1) += (dz);
			mtxY(1, 1) += (dz2);
		}
	}
	mtxX.Invert();
	mtxY.Invert();
	SVecD<2>&& rslX = mtxX * resX;
	SVecD<2>&& rslY = mtxY * resY;
	if (!MGNumc::Valid(rslX(0))) rslX(0) = fHits.at(0).CX();
	if (!MGNumc::Valid(rslY(0))) rslY(0) = fHits.at(0).CY();
	if (!MGNumc::Valid(rslX(1))) rslX(1) = ZERO;
	if (!MGNumc::Valid(rslY(1))) rslY(1) = ZERO;

	// analytics fitting for curve
	const Double_t MIN_STEP = 20.;
	Double_t resCur = ZERO; 
	Double_t mtxMag = ZERO;
	Double_t magBasedINT1 = MgntMag::GetMagFuncINT1(refZCoo);
	for (Int_t ih = 0; ih < fHits.size(); ++ih) {
		HitSt & hit = fHits.at(ih);
		if (!hit.SideY()) continue;
		Double_t dz  = hit.CZ() - refZCoo;
		Int_t    magNSTP = Int_t(std::ceil(std::fabs(dz / MIN_STEP)));
		Double_t magSTPZ = dz / Double_t(magNSTP);
		Double_t magSIGN = MGNumc::Compare(dz);
		Double_t magINT1_POS = 0.;
		for (Int_t is = 0; is < magNSTP; ++is)
			magINT1_POS += (MgntMag::GetMagFuncINT1(refZCoo + (Double_t(is) + 0.5) * magSTPZ) - magBasedINT1);
		magINT1_POS = magSIGN * std::fabs(magINT1_POS / Double_t(magNSTP));
		resCur += (rslY(0) - hit.CY()) + (dz * rslY(1));
		mtxMag += (dz * magINT1_POS);
	}
	Double_t mtxCurFact = (((fDirOpt) ? NEG : ONE) *
	                       std::sqrt(rslX(1) * rslX(1) + rslY(1) * rslY(1) + ONE) * 
												 (rslY(1) * rslY(1) + ONE));
	Double_t mtxCur = MgntProp::PROP_FACT * std::fabs(fPartSt.ChrgMass()) * mtxCurFact * mtxMag; 

	Double_t rslInvEta = NEG * resCur / mtxCur;
	if (!MGNumc::Valid(rslInvEta)) rslInvEta = 1e-6;
	
	fPartSt.setSpatialWithTan(rslX(0), rslY(0), refZCoo, rslX(1), rslY(1), ((fDirOpt) ? NEG : ONE));
	fPartSt.setInvEta(rslInvEta);

	// analytics fitting
	Int_t curIter = 1;
	Int_t maxIter = NUM_MAX_ITER_FACT_ANALYTIC * 5;
	Bool_t isSuccess = false;
	Double_t lmLambda = LM_LAMBDA_0; 
	Double_t lmChiConvg = ONE / CONVG_EPSILON;
	Double_t lmRigConvg = ONE / CONVG_EPSILON;
	while (curIter <= maxIter && !isSuccess) {
		NLLZero();

		Int_t    ndf   = -5;
		Double_t chisq = ZERO;
		SVecD<5>    gradAna;
		SMtxSymD<5> covAna;
		
		Int_t           hitCnt = 0;
		PhySt           propSt(fPartSt);
		SMtxD<5> propJb = SMtxId();
		for (Int_t ih = 0; ih < fHits.size(); ++ih) {
			HitSt & hit = fHits.at(ih);
			
			SMtxSymD<2> detErr;
			detErr(0, 0) = (!hit.SideX()) ? ZERO : ONE / hit.DetCov(0, 0);
			detErr(1, 1) = (!hit.SideY()) ? ZERO : ONE / hit.DetCov(1, 1);

			SMtxD<5> curJb;
			if (!MgntProp::PropToZ_Analytic(hit.CZ(), propSt, &curJb)) break;
			Double_t cooEst = std::sqrt(propSt.X() * propSt.X() + propSt.Y() * propSt.Y());
			if (cooEst > EPSILON_MAX_COO) break;
			propJb = curJb * propJb;
			SMtxD<2, 5>&& subJb = propJb.Sub<SMtxD<2, 5>>(0, 0);
			
			SVecD<2> detRes;
			detRes(0) = (!hit.SideX()) ? ZERO : detErr(0, 0) * (propSt.X() - hit.CX());
			detRes(1) = (!hit.SideY()) ? ZERO : detErr(1, 1) * (propSt.Y() - hit.CY());

			gradAna += LA::Transpose(subJb) * detRes;
			covAna  += LA::SimilarityT(subJb, detErr);
			
			if (hit.SideX()) { ndf++; chisq += detRes(0) * (propSt.X() - hit.CX()); }
			if (hit.SideY()) { ndf++; chisq += detRes(1) * (propSt.Y() - hit.CY()); }
			hitCnt++;
		}
		Bool_t isFinePath = ((hitCnt == fHits.size()) && (MGNumc::Compare(propSt.Mom(), EPSILON_MIN_MOM) > 0));
		Double_t normChisq = (chisq / ndf);

		Int_t  lmCurIter    = 1;
		Bool_t lmIsSuccess  = false;
		Bool_t lmIsFinePath = false;
		while (lmCurIter <= LM_NUM_MAX_ITER && !lmIsSuccess && isFinePath) {
			SMtxSymD<5> lmCovAna = lmLambda * covAna;
			for (Int_t ipar = 0; ipar < 5; ++ipar)
				lmCovAna(ipar, ipar) = covAna(ipar, ipar);

			if (!lmCovAna.Invert()) break;
			SVecD<5> && rslAna = lmCovAna * gradAna;
			for (Int_t ipar = 0; ipar < 5; ++ipar)
				if (!MGNumc::Valid(rslAna(ipar))) 
					rslAna(ipar) = ZERO;

			PhySt predSt(fPartSt.Part());
			predSt.setSpatialWithCos(
				fPartSt.X() - rslAna(0),
				fPartSt.Y() - rslAna(1),
				fPartSt.Z(),
				fPartSt.TanX() - rslAna(2),
				fPartSt.TanY() - rslAna(3),
				((fDirOpt) ? NEG : ONE)
			);
			predSt.setInvEta(fPartSt.InvEta() - rslAna(4));
			
			Int_t    predNdf   = -5;
			Double_t predChisq = ZERO;
			Int_t predHitCnt = 0;
			PhySt predPropSt(predSt);
			for (Int_t ih = 0; ih < fHits.size(); ++ih) {
				HitSt & hit = fHits.at(ih);
				if (!MgntProp::PropToZ_Analytic(hit.CZ(), predPropSt)) break;
				Double_t cooEst = std::sqrt(predPropSt.X() * predPropSt.X() + predPropSt.Y() * predPropSt.Y());
				if (cooEst > EPSILON_MAX_COO) break;
				
				SVecD<2> detRes;
				detRes(0) = (!hit.SideX()) ? ZERO : (predPropSt.X() - hit.CX()) / hit.EX();
				detRes(1) = (!hit.SideY()) ? ZERO : (predPropSt.Y() - hit.CY()) / hit.EY();

				if (hit.SideX()) { predNdf++; predChisq += detRes(0) * detRes(0); }
				if (hit.SideY()) { predNdf++; predChisq += detRes(1) * detRes(1); }
				predHitCnt++;
			}
			lmIsFinePath = ((predHitCnt == fHits.size()) && (MGNumc::Compare(predPropSt.Mom(), EPSILON_MIN_MOM) > 0));
			Double_t predNormChisq = (predChisq / predNdf);
		
			if (!lmIsFinePath) break;
		
			Double_t chiConvg = (normChisq - predNormChisq) / ((normChisq + predNormChisq) * CONVG_EPSILON + 2 * CONVG_EPSILON * CONVG_EPSILON);
			Double_t rigConvg = (std::fabs(fPartSt.InvRig() - predSt.InvRig()) / (std::fabs(fPartSt.InvRig()) + std::fabs(predSt.InvRig()))) / CONVG_EPSILON;
		
			lmIsSuccess = true;
			if (lmIsSuccess) {
				fNdf       = predNdf;
				fChisq     = predChisq;
				fNormChisq = predNormChisq;
				lmChiConvg = chiConvg; 
				lmRigConvg = rigConvg;
				fPartSt    = predSt;
			}
			else {
				Double_t width = (LM_LAMBDA_UP - LM_LAMBDA_DN) / LM_NUM_MAX_ITER;
				lmLambda = MGRndm::Uniform(LM_LAMBDA_UP - lmCurIter * width, LM_LAMBDA_UP)(); 
				lmCurIter++;
			}
		}
		
		if (!isFinePath || !lmIsFinePath) {
			fPartSt.setInvEta(fPartSt.InvEta() * EPSILON_RAT_MOM);
			curIter++;
			continue;
		}
		
		if (!lmIsSuccess) break;
		
		isSuccess = ((lmChiConvg > (NEG * CONVG_TOLERANCE) && lmChiConvg < ONE) && (lmRigConvg < ONE));
		if (!isSuccess) curIter++;
	}
	
	return isSuccess;
}


Bool_t Track::simpleFit() {
	Int_t curIter = 1;
	Int_t maxIter = NUM_MAX_ITER_FACT_SIMPLE * 5;
	Bool_t isSuccess = false;
	Double_t lmLambda = LM_LAMBDA_0; 
	Double_t lmChiConvg = ONE / CONVG_EPSILON;
	Double_t lmRigConvg = ONE / CONVG_EPSILON;
	while (curIter <= maxIter && !isSuccess) {
		NLLZero();

		Int_t    ndf   = -5;
		Double_t chisq = ZERO;
		SVecD<5>    gradG;
		SMtxSymD<5> covGG;

		Int_t hitCnt = 0;
		PhySt propSt(fPartSt);
		PhyJb propJb(PhyJb::Identity);
		for (Int_t ih = 0; ih < fHits.size(); ++ih) {
			HitSt & hit = fHits.at(ih);
			
			SMtxSymD<2> detErr;
			detErr(0, 0) = (!hit.SideX()) ? ZERO : ONE / hit.DetCov(0, 0);
			detErr(1, 1) = (!hit.SideY()) ? ZERO : ONE / hit.DetCov(1, 1);
		
			PhyJb curJb(PhyJb::Identity);
			if (!MgntProp::PropToZ(hit.CZ(), propSt, &curJb)) break;
			Double_t cooEst = std::sqrt(propSt.X() * propSt.X() + propSt.Y() * propSt.Y());
			if (cooEst > EPSILON_MAX_COO) break;
			propJb.multiplied(curJb);
			SMtxD<2, 5> && subJbG = propJb.SubJbGXY();
			
			SVecD<2> detRes;
			detRes(0) = (!hit.SideX()) ? ZERO : detErr(0, 0) * (propSt.X() - hit.CX());
			detRes(1) = (!hit.SideY()) ? ZERO : detErr(1, 1) * (propSt.Y() - hit.CY());

  		gradG += LA::Transpose(subJbG) * detRes;
  		covGG += LA::SimilarityT(subJbG, detErr);

			if (hit.SideX()) { ndf++; chisq += detRes(0) * (propSt.X() - hit.CX()); }
			if (hit.SideY()) { ndf++; chisq += detRes(1) * (propSt.Y() - hit.CY()); }
			hitCnt++;
		}
		Bool_t isFinePath = ((hitCnt == fHits.size()) && (MGNumc::Compare(propSt.Mom(), EPSILON_MIN_MOM) > 0));
		Double_t normChisq = (chisq / ndf);
		
		Int_t  lmCurIter    = 1;
		Bool_t lmIsSuccess  = false;
		Bool_t lmIsFinePath = false;
		while (lmCurIter <= LM_NUM_MAX_ITER && !lmIsSuccess && isFinePath) {
			SMtxSymD<5> lmCovGG = lmLambda * covGG;
			for (Int_t ipar = 0; ipar < 5; ++ipar)
				lmCovGG(ipar, ipar) = covGG(ipar, ipar);

			if (!lmCovGG.Invert()) break;
			SVecD<5> && rslG = lmCovGG * gradG;
			for (Int_t ipar = 0; ipar < 5; ++ipar)
				if (!MGNumc::Valid(rslG(ipar))) 
					rslG(ipar) = ZERO;
			
			PhySt predSt(fPartSt.Part());
			predSt.setSpatialWithCos(
				fPartSt.X() - rslG(0),
				fPartSt.Y() - rslG(1),
				fPartSt.Z(),
				fPartSt.DirX() - rslG(2),
				fPartSt.DirY() - rslG(3),
				fPartSt.DirZ(),
				false
			);
			predSt.setInvEta(fPartSt.InvEta() - rslG(4));
			
			Int_t predNdf      = -5;
			Double_t predChisq = ZERO;
			Int_t predHitCnt = 0;
			PhySt predPropSt(predSt);
			for (Int_t ih = 0; ih < fHits.size(); ++ih) {
				HitSt & hit = fHits.at(ih);
				if (!MgntProp::PropToZ(hit.CZ(), predPropSt)) break;
				Double_t cooEst = std::sqrt(predPropSt.X() * predPropSt.X() + predPropSt.Y() * predPropSt.Y());
				if (cooEst > EPSILON_MAX_COO) break;
				
				SVecD<2> detRes;
				detRes(0) = (!hit.SideX()) ? ZERO : (predPropSt.X() - hit.CX()) / hit.EX();
				detRes(1) = (!hit.SideY()) ? ZERO : (predPropSt.Y() - hit.CY()) / hit.EY();

				if (hit.SideX()) { predNdf++; predChisq += detRes(0) * detRes(0); }
				if (hit.SideY()) { predNdf++; predChisq += detRes(1) * detRes(1); }
				predHitCnt++;
			}
			lmIsFinePath = ((predHitCnt == fHits.size()) && (MGNumc::Compare(predPropSt.Mom(), EPSILON_MIN_MOM) > 0));
			Double_t predNormChisq = (predChisq / predNdf);
	
			if (!lmIsFinePath) break;

			Double_t chiConvg = (normChisq - predNormChisq) / ((normChisq + predNormChisq) * CONVG_EPSILON + 2 * CONVG_EPSILON * CONVG_EPSILON);
			Double_t rigConvg = (std::fabs(fPartSt.InvRig() - predSt.InvRig()) / (std::fabs(fPartSt.InvRig()) + std::fabs(predSt.InvRig()))) / CONVG_EPSILON;
			
			lmIsSuccess = true;
			if (lmIsSuccess) {
				fNdf       = predNdf;
				fChisq     = predChisq;
				fNormChisq = predNormChisq;
				lmChiConvg = chiConvg; 
				lmRigConvg = rigConvg;
				fPartSt    = predSt;
			}
			else {
				Double_t width = (LM_LAMBDA_UP - LM_LAMBDA_DN) / LM_NUM_MAX_ITER;
				lmLambda = MGRndm::Uniform(LM_LAMBDA_UP - lmCurIter * width, LM_LAMBDA_UP)(); 
				lmCurIter++;
			}
		}
		
		if (!isFinePath || !lmIsFinePath) {
			fPartSt.setInvEta(fPartSt.InvEta() * EPSILON_RAT_MOM);
			curIter++;
			continue;
		}
	
		if (!lmIsSuccess) break;
		
		isSuccess = ((lmChiConvg > (NEG * CONVG_TOLERANCE) && lmChiConvg < ONE) && (lmRigConvg < ONE));
		if (!isSuccess) curIter++;
	}

	return isSuccess;
}


///////////////////////////////////////////////////////////////////////////////////
Bool_t Track::physicalFit() {
	Bool_t optMscat = true, optEngls = false;
	std::vector<MatPropParam> matPars(fHits.size()); // physical parameters between two hit
	PhySt bkSt = fPartSt;
	for (Int_t im = 0; im < matPars.size(); ++im) {
		matPars.at(im).setOption(optMscat, optEngls);
		matPars.at(im).setMscat(0., 0.);
		matPars.at(im).setEngls(0., 0.);
	}
	
	const Int_t dimG = 5;
	const Int_t dimL = 2;

	////////TEST/////////////
	Int_t cntUpdate = 0;
	/////////////////////////

	Int_t curIter = 1;
	Int_t maxIter = NUM_MAX_ITER_FACT_PHYSICAL * ((fHits.size() - 1) * dimL + dimG);
	Bool_t isSuccess = false;
	Double_t lmLambda = LM_LAMBDA_0; 
	Double_t lmChiConvg = ONE / CONVG_EPSILON;
	Double_t lmRigConvg = ONE / CONVG_EPSILON;
	while (curIter <= maxIter && !isSuccess) {
		NLLZero();

		Int_t    ndf   = -5;
		Double_t chisq = ZERO;
		std::vector<SVecD<2>>    detRes(fHits.size());
		std::vector<SMtxSymD<2>> detErr(fHits.size());
		std::vector<PhyJb>              jacb(fHits.size());

		Int_t hitCnt = 0;
		PhySt propSt(fPartSt);
		std::vector<Double_t> nrl(fHits.size());
		for (Int_t ih = 0; ih < fHits.size(); ++ih) {
			HitSt & hit = fHits.at(ih);
			
			detErr.at(ih)(0, 0) = (!hit.SideX()) ? 0. : (1. / hit.EX() / hit.EX());
			detErr.at(ih)(1, 1) = (!hit.SideY()) ? 0. : (1. / hit.EY() / hit.EY());

			PhyJb curJb(PhyJb::Identity);
			if (!MgntProp::PropToZ(hit.CZ(), propSt, &curJb, &matPars.at(ih))) break;
			Double_t cooEst = std::sqrt(propSt.X() * propSt.X() + propSt.Y() * propSt.Y());
			if (cooEst > EPSILON_MAX_COO) break;
			jacb.at(ih) = curJb;
		
			///////////TEST////////
			//COUT("ITER %02d RIG %14.8f\n", curIter, propSt.Rig());
			//////////////////////

			detRes.at(ih)(0) = (!hit.SideX()) ? ZERO : detErr.at(ih)(0, 0) * (propSt.X() - hit.CX());
			detRes.at(ih)(1) = (!hit.SideY()) ? ZERO : detErr.at(ih)(1, 1) * (propSt.Y() - hit.CY());

			if (hit.SideX()) { ndf++; chisq += detRes.at(ih)(0) * (propSt.X() - hit.CX()); }
			if (hit.SideY()) { ndf++; chisq += detRes.at(ih)(1) * (propSt.Y() - hit.CY()); }
			SVecD<4> intchi = matPars.at(ih).Chisq(curJb.NumRadLen());
			chisq += intchi(0) + intchi(1);
			nrl.at(ih) = curJb.NumRadLen();	
			hitCnt++;
		}
		Bool_t isFinePath = ((hitCnt == fHits.size()) && (MGNumc::Compare(propSt.Mom(), EPSILON_MIN_MOM) > 0));
		Double_t normChisq = (chisq / ndf);
		
		//------------ SAT ---------------//
		SVecD<dimG>    gradG;
		SMtxSymD<dimG> covGG;
		std::vector<SVecD<dimL>> gradL(matPars.size());
		std::vector<SMtxD<dimL>> covLL(matPars.size()*matPars.size());
		std::vector<SMtxD<dimG, dimL>> covGL(matPars.size());
		{
			std::vector<SMtxD<2, dimG>> subJbG(fHits.size());
			SMtxD<dimG> jacbG = SMtxId();
			for (Int_t ih = 0; ih < fHits.size(); ++ih) {
				jacbG = jacb.at(ih).G() * jacbG;
				subJbG.at(ih) = jacbG.Sub<SMtxD<2, dimG>>(0, 0);
			}
			
			std::vector<SMtxD<2, dimL>> subJbL(matPars.size()*fHits.size());
			for (Int_t im = 1; im < matPars.size(); ++im) {
    		SMtxD<dimG, dimL> jacbL = jacb.at(im).L().Sub<SMtxD<dimG, dimL>>(0, 0);
				for (Int_t ih = im; ih < fHits.size(); ++ih) {
					if (ih != im) jacbL = jacb.at(ih).G() * jacbL;
					subJbL.at(im*fHits.size()+ih) = jacbL.Sub<SMtxD<2, dimL>>(0, 0);
				}
			}
	
			for (Int_t ih = 0; ih < fHits.size(); ++ih) {
				gradG += LA::Transpose(subJbG.at(ih)) * detRes.at(ih);
				covGG += LA::SimilarityT(subJbG.at(ih), detErr.at(ih));
			}	
		
			for (Int_t im = 1; im < matPars.size(); ++im) {
				gradL.at(im) += matPars.at(im).Grad(jacb.at(im).NumRadLen()).Sub<SVecD<dimL>>(0); 
				for (Int_t ih = im; ih < fHits.size(); ++ih) {
					gradL.at(im) += LA::Transpose(subJbL.at(im*fHits.size()+ih)) * detRes.at(ih);
					covGL.at(im) += LA::Transpose(subJbG.at(ih)) * detErr.at(ih) * subJbL.at(im*fHits.size()+ih);
				}
				covLL.at(im*matPars.size()+im) += matPars.at(im).Cov(jacb.at(im).NumRadLen()).Sub<SMtxSymD<dimL>>(0, 0); 
			}

			for (Int_t ih = 1; ih < fHits.size(); ++ih) {
				for (Int_t im = 1; im <= ih; ++im) {
					for (Int_t jm = 1; jm <= ih; ++jm) {
						covLL.at(im*matPars.size()+jm) += LA::Transpose(subJbL.at(im*fHits.size()+ih)) * detErr.at(ih) * subJbL.at(jm*fHits.size()+ih);
					}
				}
			}
		}

		Int_t sizeF = dimG + (matPars.size() - 1) * dimL;
		TVecD gradF(sizeF);
		TMtxD covFF(sizeF, sizeF);
		for (Int_t ip = 0; ip < dimG; ++ip) {
			gradF(ip) = gradG(ip);
		}
		for (Int_t im = 1; im < matPars.size(); ++im) {
			for (Int_t it = 0; it < dimL; ++it) {
				gradF(dimG+dimL*(im-1)+it) = gradL.at(im)(it);
			}
		}

		for (Int_t ip = 0; ip < dimG; ++ip) {
			for (Int_t jp = 0; jp < dimG; ++jp) {
				covFF(ip, jp) = covGG(ip, jp);
			}
			for (Int_t im = 1; im < matPars.size(); ++im) {
				for (Int_t it = 0; it < dimL; ++it) {
					covFF(ip, dimG+dimL*(im-1)+it) = covGL.at(im)(ip, it);
					covFF(dimG+dimL*(im-1)+it, ip) = covGL.at(im)(ip, it);
				}
			}
		}

		for (Int_t im = 1; im < matPars.size(); ++im) {
			for (Int_t it = 0; it < dimL; ++it) {
				for (Int_t jm = 1; jm < matPars.size(); ++jm) {
					for (Int_t jt = 0; jt < dimL; ++jt) {
						covFF(dimG+dimL*(im-1)+it, dimG+dimL*(jm-1)+jt) = covLL.at(im*matPars.size()+jm)(it, jt);
					}
				}
			}
		}
		//--------------- END -------------//
		
		Int_t  lmCurIter    = 1;
		Bool_t lmIsSuccess  = false;
		Bool_t lmIsFinePath = false;
		while (lmCurIter <= LM_NUM_MAX_ITER && !lmIsSuccess && isFinePath) {
			TMtxD lmCovFF = lmLambda * covFF;
			for (Int_t ipar = 0; ipar < sizeF; ++ipar)
				lmCovFF(ipar, ipar) = covFF(ipar, ipar);
			
			Double_t det = 0;
			lmCovFF.Invert(&det);
			if (!MGNumc::Valid(det) || MGNumc::EqualToZero(det)) break;
			TVecD rslF = lmCovFF * gradF;
			for (Int_t ipar = 0; ipar < sizeF; ++ipar)
				if (!MGNumc::Valid(rslF(ipar))) 
					rslF(ipar) = ZERO;

			//------TEST--------//
			rslF *= (fPartSt.Beta());
			//-------------------//


			SVecD<dimG> rslG;
			for (Int_t ipar = 0; ipar < dimG; ++ipar)
				rslG(ipar) = rslF(ipar);

			std::vector<SVecD<dimL>> rslL(matPars.size());
			for (Int_t im = 1; im < matPars.size(); ++im)
				for (Int_t it = 0; it < dimL; ++it)
					rslL.at(im)(it) = rslF(dimG+(im-1)*dimL+it);
		
			PhySt predSt(fPartSt.Part());
			predSt.setSpatialWithCos(
				fPartSt.X() - rslG(0),
				fPartSt.Y() - rslG(1),
				fPartSt.Z(),
				fPartSt.DirX() - rslG(2),
				fPartSt.DirY() - rslG(3),
				fPartSt.DirZ(),
				false
			);
			predSt.setInvEta(fPartSt.InvEta() - rslG(4));

			std::vector<MatPropParam> predMatPars(fHits.size());
			for (Int_t im = 0; im < predMatPars.size(); ++im) {
				predMatPars.at(im).setOption(optMscat, optEngls);
				predMatPars.at(im).setMscat(
					(dimL >= 1) ? matPars.at(im).MscatDV() - rslL.at(im)(0) : ZERO, 
					(dimL >= 2) ? matPars.at(im).MscatDW() - rslL.at(im)(1) : ZERO
				);
				predMatPars.at(im).setEngls(
					(dimL >= 3) ? matPars.at(im).EnglsIN() - rslL.at(im)(2) : ZERO, 
					(dimL >= 4) ? matPars.at(im).EnglsBR() - rslL.at(im)(3) : ZERO
				);
			}

			Int_t predNdf      = -5;
			Double_t predChisq = ZERO;
			Int_t predHitCnt = 0;
			PhySt predPropSt(predSt);
			//////////TEST////////////
			Double_t _nllx = ZERO, _nlly = ZERO;
			Int_t _ndfx = -2, _ndfy = -3;
			//////////////////////////
			for (Int_t ih = 0; ih < fHits.size(); ++ih) {
				HitSt & hit = fHits.at(ih);
				if (!MgntProp::PropToZ(hit.CZ(), predPropSt, nullptr, &predMatPars.at(ih))) break;
				Double_t cooEst = std::sqrt(predPropSt.X() * predPropSt.X() + predPropSt.Y() * predPropSt.Y());
				if (cooEst > EPSILON_MAX_COO) break;
			
				SVecD<2> detRes;
				detRes(0) = (!hit.SideX()) ? ZERO : ((predPropSt.X() - hit.CX()) / hit.EX());
				detRes(1) = (!hit.SideY()) ? ZERO : ((predPropSt.Y() - hit.CY()) / hit.EY());
			
				if (hit.SideX()) { predNdf++; predChisq += detRes(0) * detRes(0); }
				if (hit.SideY()) { predNdf++; predChisq += detRes(1) * detRes(1); }
				SVecD<4> intchi = predMatPars.at(ih).Chisq(nrl.at(ih));
				predChisq += intchi(0) + intchi(1);
				
				predHitCnt++;
			
				//////////TEST////////////
				if (hit.SideX()) _ndfx += 1;
				if (hit.SideY()) _ndfy += 1;
				_nllx += detRes(0) * detRes(0) + intchi(1), 
				_nlly += detRes(1) * detRes(1) + intchi(0);
				//////////////////////////
			}
			lmIsFinePath = ((predHitCnt == fHits.size()) && (MGNumc::Compare(predPropSt.Mom(), EPSILON_MIN_MOM) > 0));
			Double_t predNormChisq = (predChisq / predNdf);
			
			if (!lmIsFinePath) break;

			Double_t chiConvg = (normChisq - predNormChisq) / ((normChisq + predNormChisq) * CONVG_EPSILON + 2 * CONVG_EPSILON * CONVG_EPSILON);
			Double_t rigConvg = (std::fabs(fPartSt.InvRig() - predSt.InvRig()) / (std::fabs(fPartSt.InvRig()) + std::fabs(predSt.InvRig()))) / CONVG_EPSILON;
			
			//lmIsSuccess = ((curIter <= LM_IGNORE_ITER) || chiConvg > (NEG * CONVG_TOLERANCE));
			lmIsSuccess = true;
			//---------//
			//Double_t stepRSL = std::sqrt(rslF.Norm2Sqr());
			//Bool_t succ = ((chiConvg > (NEG * CONVG_TOLERANCE) && chiConvg < ONE) && (rigConvg < ONE)) || (lmRigConvg < ONE && lmCurIter == LM_NUM_MAX_ITER);
			//if (lmIsSuccess ||  lmCurIter != 1) COUT("PHY ITER %03d(%03d)/%02d(%01d) CHISQ (%10.2f, %10.2f) STEP(%4.2f) %14.8f CONVG %14.8f %14.8f\n", curIter, maxIter, lmCurIter, lmIsFinePath, normChisq, predNormChisq, lmLambda, stepRSL, chiConvg, rigConvg);
			//if (succ || lmCurIter != 1) COUT("PHY ITER %03d(%03d)/%02d(%01d) CHISQ (%10.2f, %10.2f) STEP(%4.2f) %14.8f CONVG %14.8f %14.8f\n", curIter, maxIter, lmCurIter, lmIsFinePath, normChisq, predNormChisq, lmLambda, stepRSL, chiConvg, rigConvg);
			//----------//

			if (lmIsSuccess) {
				fNdf       = predNdf;
				fChisq     = predChisq;
				fNormChisq = predNormChisq;
				lmChiConvg = chiConvg;
				lmRigConvg = rigConvg;
				fPartSt    = predSt;
				matPars    = predMatPars;
				//////////TEST////////////
				fNDFX = _ndfx;
				fNLLX = _nllx;
				fNDFY = _ndfy;
				fNLLY = _nlly;
				fNormNLLX = fNLLX / fNDFX;
				fNormNLLY = fNLLY / fNDFY;
				//////////////////////////
			}
			else {
				Double_t width = (LM_LAMBDA_UP - LM_LAMBDA_DN) / LM_NUM_MAX_ITER;
				lmLambda = MGRndm::Uniform(LM_LAMBDA_UP - lmCurIter * width, LM_LAMBDA_UP)(); 
				lmCurIter++; 
			}
		}
		
		if (!isFinePath || !lmIsFinePath) {
			//////////TEST///////
			cntUpdate++;
			//COUT("UPDATE %03d COO %14.8f %14.8f DIR %14.8f %14.8f RIG %14.8f\n", cntUpdate, fPartSt.X(), fPartSt.Y(), fPartSt.DirX(), fPartSt.DirY(), fPartSt.Rig());
			Double_t inveta = fPartSt.InvEta() * 0.85;
			fPartSt = bkSt;
			fPartSt.setInvEta(inveta);
			//fPartSt.setSpatialCos(fPartSt.DirX() * 0.85, fPartSt.DirY() * 0.85,  fPartSt.DirZ(), false);
			////////////////////
			for (Int_t im = 0; im < matPars.size(); ++im) {
				matPars.at(im).setMscat(0., 0.);
				matPars.at(im).setEngls(0., 0.);
			}
			//fPartSt = bkSt;
			//fPartSt.setInvEta(fPartSt.InvEta() * EPSILON_RAT_MOM);
			curIter++;
			if (curIter >= maxIter) isSuccess = true;
			continue;
		}
		
		if (!lmIsSuccess) break;

		isSuccess = ((lmChiConvg > (NEG * CONVG_TOLERANCE) && lmChiConvg < ONE) && (lmRigConvg < ONE));
		//////////TEST////////////
		if (curIter >= maxIter) isSuccess = true;
		//////////////////////////
		if (!isSuccess) curIter++;
	}
	fMatPars   = matPars;

	//////////TEST////////////
	if      (!isSuccess)     COUT("FAILE.(%03d/%03d) UPDATE %03d CONVG ( %14.8f %14.8f )\n", curIter, maxIter, cntUpdate, lmChiConvg, lmRigConvg);
	if      (cntUpdate >= 1) COUT("UPDAT.(%03d/%03d) UPDATE %03d CONVG ( %14.8f %14.8f )\n", curIter, maxIter, cntUpdate, lmChiConvg, lmRigConvg);
	if (cntUpdate >= 1 || !isSuccess) COUT("COO %14.8f %14.8f DIR %14.8f %14.8f RIG %14.8f\n", fPartSt.X(), fPartSt.Y(), fPartSt.DirX(), fPartSt.DirY(), fPartSt.Rig());
	//std::cout << "--------------------------------------------------------------------------------------------\n";
	//////////////////////////
	
	return isSuccess;
}



//////////////////////////////////////////////////////////////////////////////////////////

/*
// TODO
Bool_t Track::physicalFit() {
	std::vector<MatPropParam> matPars(fHits.size()); // physical parameters between two hit
	for (Int_t im = 0; im < matPars.size(); ++im) {
		matPars.at(im).setOption(true, false);
		matPars.at(im).setMscat(0., 0.);
		matPars.at(im).setEngls(0., 0.);
	}
	std::vector<Double_t> vrig;
	std::vector<Double_t> vchi;
	std::vector<Double_t> vchix;
	std::vector<Double_t> vchiy;
	std::vector<Double_t> vchiv;
	std::vector<Double_t> vchiw;
	std::vector<Double_t> vdix;
	std::vector<Double_t> vovf;
	fNdf       = 0;
	fChisq     = ZERO;
	fNormChisq = ZERO;

	Int_t decayIter = 0;
	const Double_t decayStep = 5.;

	const Int_t dimG = 5;
	const Int_t dimL = 2;

	Int_t curIter = 1;
	const Int_t MaxIter = 200; // 25
	Bool_t isSuccess = false;
	while (curIter <= MaxIter && !isSuccess) {
		Int_t    ndfx  = -2;
		Int_t    ndfy  = -3;
		Double_t chisqx = ZERO;
		Double_t chisqy = ZERO;
		std::vector<SVecD<2>>    detRes(fHits.size());
		std::vector<SMtxSymD<2>> detErr(fHits.size());
		std::vector<PhyJb>              jacb(fHits.size());

		Int_t curIH = 0;
		Double_t numRadLen = 0.0;
		PhySt propSt(fPartSt);
		for (Int_t ih = 0; ih < fHits.size(); ++ih) {
			HitSt & hit = fHits.at(ih);
			detErr.at(ih)(0, 0) = (!hit.SideX()) ? 0. : (1. / hit.EX() / hit.EX());
			detErr.at(ih)(1, 1) = (!hit.SideY()) ? 0. : (1. / hit.EY() / hit.EY());

			PhyJb curJb(PhyJb::Identity);
			if ( ! MgntProp::PropToZ(hit.CZ(), propSt, &curJb, &matPars.at(ih)) ) break;
			jacb.at(ih) = curJb;
			
			detRes.at(ih)(0) = (!hit.SideX()) ? ZERO : detErr.at(ih)(0, 0) * (propSt.X() - hit.CX());
			detRes.at(ih)(1) = (!hit.SideY()) ? ZERO : detErr.at(ih)(1, 1) * (propSt.Y() - hit.CY());

			if (hit.SideX()) { ndfx++; chisqx += detRes.at(ih)(0) * (propSt.X() - hit.CX()); }
			if (hit.SideY()) { ndfy++; chisqy += detRes.at(ih)(1) * (propSt.Y() - hit.CY()); }
			//chisq += matPars.at(ih).Chisq(curJb.NumRadLen());
			numRadLen += curJb.NumRadLen();
			curIH = ih;
		}

		// landau distribution (ion-loss)
		if (MGNumc::Compare(propSt.Mom(), 1e-4) <= 0 && curIH != fHits.size()-1) {
			Double_t rat = 0.85;  // best 0.85
			Int_t type = 0;
			if (curIH <= 1) { fPartSt.setInvEta(MGNumc::Compare(fPartSt.InvEta())); type = 1; }
			if (std::fabs(fPartSt.Eta()) < 4. * numRadLen) { fPartSt.setInvEta(MGNumc::Compare(fPartSt.InvEta()) / (4. * numRadLen)); type = 2; }
			else { fPartSt.setInvEta(fPartSt.InvEta() * rat); type = 3; }
			for (Int_t im = 0; im < matPars.size(); ++im) {
				if (dimL >= 1) matPars.at(im).setMscat(0., 0.);
				if (dimL >= 3) matPars.at(im).setEngls(0., 0.);
			}
			std::cout << "ITER " << curIter << " TYPE " << type << " MIN_MOM " << fPartSt.Rig() << "\n";
			curIter++;
			continue;
		}
		
		Double_t chisqHX = chisqx;
		Double_t chisqHY = chisqy;
		Double_t chisqDV = 0.;
		Double_t chisqDW = 0.;
		Double_t chisqIN = 0.;
		Double_t chisqBR = 0.;
		for (Int_t im = 1; im < matPars.size(); ++im) {
			chisqDV += (matPars.at(im).MscatDV() * matPars.at(im).MscatDV());
			chisqDW += (matPars.at(im).MscatDW() * matPars.at(im).MscatDW());
			chisqIN += (matPars.at(im).EnglsIN() + std::exp(-matPars.at(im).EnglsIN()) - 1.0);
			Double_t bremslen = jacb.at(im).NumRadLen() * 1.44269504088896339e+00;
			Double_t englsBR = (MGNumc::Compare(matPars.at(im).EnglsBR(), 1.0e-7) > 0) ? matPars.at(im).EnglsBR() : 1.0e-7;
			chisqBR += ((englsBR) * bremslen + (ONE - bremslen) * (std::log(englsBR) - std::log(1.0e-7)));
		}
		Double_t normChisq    = (chisqHX+chisqHY+chisqDV+chisqDW)/(ndfx+ndfy);
		Double_t dix = (fNormChisq - normChisq)/((normChisq + fNormChisq) * CONVG_EPSILON + 2 * CONVG_EPSILON * CONVG_EPSILON);
		Double_t ovf = 0.5 * (normChisq + fNormChisq) / dix;

		if (curIter > 5 && fNormChisq > normChisq) { decayIter++; }

		if (dix < ONE && dix > ZERO) { isSuccess = true; break; }
	
		//----------------------------------------------------------------//
		SVecD<dimG>    gradG;
		SMtxSymD<dimG> covGG;
		std::vector<SVecD<dimL>> gradL(matPars.size());
		std::vector<SMtxD<dimL>> covLL(matPars.size()*matPars.size());
		std::vector<SMtxD<dimG, dimL>> covGL(matPars.size());
		{
			std::vector<SMtxD<2, dimG>> subJbG(fHits.size());
			SMtxD<dimG> jacbG = SMtxId();
			for (Int_t ih = 0; ih < fHits.size(); ++ih) {
				jacbG = jacb.at(ih).G() * jacbG;
				subJbG.at(ih) = jacbG.Sub<SMtxD<2, dimG>>(0, 0);
			}
			
			std::vector<SMtxD<2, dimL>> subJbL(matPars.size()*fHits.size());
			for (Int_t im = 1; im < matPars.size(); ++im) {
    		SMtxD<dimG, dimL> jacbL = jacb.at(im).L().Sub<SMtxD<dimG, dimL>>(0, 0);
				for (Int_t ih = im; ih < fHits.size(); ++ih) {
					if (ih != im) jacbL = jacb.at(ih).G() * jacbL;
					subJbL.at(im*fHits.size()+ih) = jacbL.Sub<SMtxD<2, dimL>>(0, 0);
				}
			}
	
			for (Int_t ih = 0; ih < fHits.size(); ++ih) {
				gradG += LA::Transpose(subJbG.at(ih)) * detRes.at(ih);
				covGG += LA::SimilarityT(subJbG.at(ih), detErr.at(ih));
			}	
		
			for (Int_t im = 1; im < matPars.size(); ++im) {
				gradL.at(im) += matPars.at(im).Grad(jacb.at(im).NumRadLen()).Sub<SVecD<dimL>>(0); 
				for (Int_t ih = im; ih < fHits.size(); ++ih) {
					gradL.at(im) += LA::Transpose(subJbL.at(im*fHits.size()+ih)) * detRes.at(ih);
					covGL.at(im) += LA::Transpose(subJbG.at(ih)) * detErr.at(ih) * subJbL.at(im*fHits.size()+ih);
				}
				covLL.at(im*matPars.size()+im) += matPars.at(im).Cov(jacb.at(im).NumRadLen()).Sub<SMtxSymD<dimL>>(0, 0); 
			}

			for (Int_t ih = 1; ih < fHits.size(); ++ih) {
				for (Int_t im = 1; im <= ih; ++im) {
					for (Int_t jm = 1; jm <= ih; ++jm) {
						covLL.at(im*matPars.size()+jm) += LA::Transpose(subJbL.at(im*fHits.size()+ih)) * detErr.at(ih) * subJbL.at(jm*fHits.size()+ih);
					}
				}
			}
		}

		Int_t sizeF = dimG + (matPars.size() - 1) * dimL;
		TVecD gradF(sizeF);
		TMtxD covFF(sizeF, sizeF);
		for (Int_t ip = 0; ip < dimG; ++ip) {
			gradF(ip) = gradG(ip);
		}
		for (Int_t im = 1; im < matPars.size(); ++im) {
			for (Int_t it = 0; it < dimL; ++it) {
				gradF(dimG+dimL*(im-1)+it) = gradL.at(im)(it);
			}
		}

		for (Int_t ip = 0; ip < dimG; ++ip) {
			for (Int_t jp = 0; jp < dimG; ++jp) {
				covFF(ip, jp) = covGG(ip, jp);
			}
			for (Int_t im = 1; im < matPars.size(); ++im) {
				for (Int_t it = 0; it < dimL; ++it) {
					covFF(ip, dimG+dimL*(im-1)+it) = covGL.at(im)(ip, it);
					covFF(dimG+dimL*(im-1)+it, ip) = covGL.at(im)(ip, it);
				}
			}
		}

		for (Int_t im = 1; im < matPars.size(); ++im) {
			for (Int_t it = 0; it < dimL; ++it) {
				for (Int_t jm = 1; jm < matPars.size(); ++jm) {
					for (Int_t jt = 0; jt < dimL; ++jt) {
						covFF(dimG+dimL*(im-1)+it, dimG+dimL*(jm-1)+jt) = covLL.at(im*matPars.size()+jm)(it, jt);
					}
				}
			}
		}

		// testcode
		//for (Int_t i = 0; i < sizeF; ++i) {
		//	covFF(i, i) *= 1.01;
		//}
		///////////


		Double_t det = 0;
		covFF.Invert(&det);
		if (!MGNumc::Valid(det) || MGNumc::EqualToZero(det)) { std::cout << "INVERSE FAIL.\n"; break; }
		TVecD rslF = covFF * gradF;
	
		Double_t decayFT = std::exp(-1. * Double_t(decayIter) / decayStep);
		//rslF *= decayFT;

		SVecD<dimG> rslG;
		for (Int_t ip = 0; ip < dimG; ++ip) {
			rslG(ip) = rslF(ip);
		}

		std::vector<SVecD<dimL>> rslL(matPars.size());
		for (Int_t im = 1; im < matPars.size(); ++im) {
			for (Int_t it = 0; it < dimL; ++it) {
				rslL.at(im)(it) = rslF(dimG+(im-1)*dimL+it);
			}
		}


		for (Int_t im = 1; im < matPars.size(); ++im) {
			if (dimL >= 1) matPars.at(im).setMscat((dimL >= 1) ? matPars.at(im).MscatDV() - rslL.at(im)(0) : 0., (dimL >= 2) ? matPars.at(im).MscatDW() - rslL.at(im)(1) : 0.);
			if (dimL >= 3) matPars.at(im).setEngls((dimL >= 3) ? matPars.at(im).EnglsIN() - rslL.at(im)(2) : 0., (dimL >= 4) ? matPars.at(im).EnglsBR() - rslL.at(im)(3) : 0.);
			//if (curIter%5==0) COUT("ITER %2d  RIG %14.8f  IT %d ION %14.8f BREM %14.8f  %14.8f\n", curIter, fPartSt.Rig(), im, matPars.at(im).EnglsIN(), matPars.at(im).EnglsBR(), rslL.at(im)(3));
		}
		//----------------------------------------------------------------//

		//if (!covGG.Invert()) break;
		//SVecD<dimG> && rslG = covGG * gradG;
		for (Int_t ipar = 0; ipar < dimG; ++ipar)
			if (!MGNumc::Valid(rslG(ipar))) 
				rslG(ipar) = 0.;
		
		fPartSt.setSpatialWithCos(
			fPartSt.X() - rslG(0),
			fPartSt.Y() - rslG(1),
			fPartSt.Z(),
			fPartSt.DirX() - rslG(2),
			fPartSt.DirY() - rslG(3),
			fPartSt.DirZ(),
			false
		);
		fPartSt.setInvEta(fPartSt.InvEta() - rslG(4));


		

		
		//fNdf = ndfx + ndfy + 2 * (matPars.size() - 1);
		fNdf = ndfx + ndfy;
		fChisq = chisqHX+chisqHY+chisqDV+chisqDW;
		fNormChisq = normChisq;
		
		//chisqHX /= ndfx;
		//chisqHY /= ndfy;
		//chisqDV /= (matPars.size() - 1); 
		//chisqDW /= (matPars.size() - 1); 
		//chisqIN /= (matPars.size() - 1); 
		//chisqBR /= (matPars.size() - 1); 

		vrig.push_back(fPartSt.Rig());
		vchi.push_back(fNormChisq);
		vchix.push_back(chisqHX);
		vchiy.push_back(chisqHY);
		vchiv.push_back(chisqDV);
		vchiw.push_back(chisqDW);
		vdix.push_back(dix);
		vovf.push_back(ovf);

		if (dix < ONE && dix > ZERO) { isSuccess = true; break; }
		else curIter++;
		///curIter++;
	}
	
	if (!isSuccess || curIter > MaxIter) {
		//std::cout << "Succ : " << isSuccess << " Iter : " << curIter << "/" << MaxIter << std::endl;
		//for (Int_t it = 0; it < vrig.size(); ++it) {
		//	Double_t rat =  ((it==0) ? 0.0 : std::fabs((1./vrig.at(it)-1./vrig.at(it-1))/(1./vrig.at(it)+1./vrig.at(it-1))));
		//	std::cout << Form("IT %2d RIG %14.8f (%14.8f) CHI %14.8f X %14.8f Y %14.8f V %14.8f W %14.8f DIX %14.8f OVF %14.8f\n", it, vrig.at(it), rat, vchi.at(it), vchix.at(it), vchiy.at(it), vchiv.at(it), vchiw.at(it), vdix.at(it), vovf.at(it));
		//}
		//std::cout << "--------------------------- END " << std::endl;
	}

	fMatPars = matPars;

	return isSuccess;
}

*/














/*
//---- Vertex ----//
Bool_t Vertex::rebuild() {
	Bool_t isSuccess = fit();
  if (!isSuccess) init();
  return isSuccess;
}

Bool_t Vertex::fit() {
	Bool_t isSuccess = true;
  if (NumOfTrack() < 2) return false;
  if (simpleFit() < 0) return false;

	Int_t iter = 1;
	const Int_t MaxIter = 5;
	const Double_t TOL = 1.e-2;
	const Double_t MinChisq = 1.e-2;
	Double_t iterChisq[2] = { -1, -1};

	while (iter <= MaxIter) {
		Double_t chisq = miniParam();
		if (chisq < 0 || std::isnan(chisq) || std::isinf(chisq)) { isSuccess = false; break; }

		Double_t dTOL[2] = {0};
		dTOL[0] = std::fabs(chisq - iterChisq[0]) / (chisq + MinChisq);
		dTOL[1] = std::fabs(chisq - iterChisq[1]) / (chisq + MinChisq);
		iterChisq[1] = iterChisq[0];
		iterChisq[0] = chisq;

		if (iter >= 2 && dTOL[0] < TOL && dTOL[1] < TOL) break;
		iter++;
	}
	if (!isSuccess) init();

	return isSuccess;
}

Double_t Vertex::simpleFit() {
  const Double_t Shift = 3.; // [cm]
  const Double_t ExtLimit = 220.; // [cm]

  // boundary
  Double_t Boundary[2] = { 0, 0 }; // (max, min)
	for (Int_t itr = 0; itr < NumOfTrack(); ++itr) {
    if (TrackStatusRef(itr).Z() > Boundary[0] || Boundary[0] == 0) Boundary[0] = TrackStatusRef(itr).Z();
    if (TrackStatusRef(itr).Z() < Boundary[1] || Boundary[1] == 0) Boundary[1] = TrackStatusRef(itr).Z();
  }
  Boundary[0] -= Shift;
  Boundary[1] += Shift;

  // share cluster (remove this vertex)
  Bool_t IsShareClusterX  = false;
  Bool_t IsShareClusterY  = false;
  Bool_t IsShareClusterXY = false;

  for (Int_t itr = 0; itr < NumOfTrack()-1; ++itr) {
    std::vector<TrHit *> & iTrHits = GetTrack(itr)->TrHits();
	  for (Int_t jtr = itr+1; jtr < NumOfTrack(); ++jtr) {
      std::vector<TrHit *> & jTrHits = GetTrack(jtr)->TrHits();

      for (Int_t ih = 0; ih < GetTrack(itr)->NumOfTrHit(); ++ih) {
        TrHit * iTrHit = iTrHits.at(ih);
        for (Int_t jh = 0; jh < GetTrack(jtr)->NumOfTrHit(); ++jh) {
          TrHit * jTrHit = jTrHits.at(jh);

          Bool_t shareY  = (iTrHit->ClsY() == jTrHit->ClsY());
          Bool_t shareX  = (iTrHit->ClsX() != -1 && jTrHit->ClsX() != -1) ? (iTrHit->ClsX() == jTrHit->ClsX()) : false;
          Bool_t shareXY = (shareX && shareY);

          if (shareX)  IsShareClusterX  = shareX;
          if (shareY)  IsShareClusterY  = shareY;
          if (shareXY) IsShareClusterXY = shareXY;

          if (IsShareClusterX && IsShareClusterY && IsShareClusterXY) break;
        } // j-hits
        if (IsShareClusterX && IsShareClusterY && IsShareClusterXY) break;
      } // i-hits

      if (IsShareClusterX && IsShareClusterY && IsShareClusterXY) break;
    } // j-track
    if (IsShareClusterX && IsShareClusterY && IsShareClusterXY) break;
  } // i-track

  if (IsShareClusterX || IsShareClusterY || IsShareClusterXY) return -1;

	// linear fit for x-dir
  SVecD<2> vtxXZ;
  Bool_t LinearFitX = true;
  if (LinearFitX) {
    const Int_t nIterLnX = 1;
    for (Int_t iter = 1; iter <= nIterLnX; ++iter) {
      SMtxSymD<2> coefXZ;
      SVecD<2>    gradXZ;
      for (Int_t itr = 0; itr < NumOfTrack(); ++itr) {
		  	PhySt & st = TrackStatusRef(itr);
        Double_t tx = st.TanX();
        Double_t invErr2 = 1. / (1. + tx * tx);
        coefXZ(0, 0) += invErr2;
        coefXZ(0, 1) += invErr2 * (-1. * tx);
        coefXZ(1, 1) += invErr2 * (tx * tx);
        gradXZ(0)    += invErr2 * (st.X() - tx * st.Z());
        gradXZ(1)    += invErr2 * (-1. * tx) * (st.X() - tx * st.Z());
      }
      Bool_t isSucs = coefXZ.Invert();
      vtxXZ = coefXZ * gradXZ;
      if (!isSucs) return -11;
      if (std::isinf(vtxXZ(0)) || std::isnan(vtxXZ(0))) return -12;
      if (std::isinf(vtxXZ(1)) || std::isnan(vtxXZ(1))) return -13;

      //if (vtxXZ(1) < Boundary[0] && vtxXZ(1) > Boundary[1]) return -14;
      if (std::fabs(vtxXZ(0)) > ExtLimit || std::fabs(vtxXZ(1)) > ExtLimit) return -15;

      Double_t normDist = 0;
      for (Int_t itr = 0; itr < NumOfTrack(); ++itr) {
		  	PhySt & curSt = TrackStatusRef(itr);
		    MgntProp::Prop(vtxXZ(1) - curSt.Z(), curSt);
        Double_t dist = std::fabs(vtxXZ(0) - curSt.X()) / std::sqrt(1. + curSt.TanX() * curSt.TanX());
        normDist += dist;
      }
      normDist /= NumOfTrack();
      if (normDist > Shift) return -16;
    }
  }

  // linear fit for y-dir
  SVecD<2> vtxYZ;
  Bool_t LinearFitY = true;
  if (LinearFitY) {
    const Int_t nIterLnY = 2;
    for (Int_t iter = 1; iter <= nIterLnY; ++iter) {
      SMtxSymD<2> coefYZ;
      SVecD<2>    gradYZ;
      for (Int_t itr = 0; itr < NumOfTrack(); ++itr) {
		  	PhySt & st = TrackStatusRef(itr);
        Double_t ty = st.TanY();
        Double_t invErr2 = 1. / (1. + ty * ty);
        coefYZ(0, 0) += invErr2;
        coefYZ(0, 1) += invErr2 * (-1. * ty);
        coefYZ(1, 1) += invErr2 * (ty * ty);
        gradYZ(0)    += invErr2 * (st.Y() - ty * st.Z());
        gradYZ(1)    += invErr2 * (-1. * ty) * (st.Y() - ty * st.Z());
      }
      Bool_t isSucs = coefYZ.Invert();
      vtxYZ = coefYZ * gradYZ;
      if (!isSucs) return -21;
      if (std::isinf(vtxYZ(0)) || std::isnan(vtxYZ(0))) return -22;
      if (std::isinf(vtxYZ(1)) || std::isnan(vtxYZ(1))) return -23;

      //if (vtxYZ(1) < Boundary[0] && vtxYZ(1) > Boundary[1]) return -24;
      if (std::fabs(vtxYZ(0)) > ExtLimit || std::fabs(vtxYZ(1)) > ExtLimit) return -25;

      Double_t normDist = 0;
      for (Int_t itr = 0; itr < NumOfTrack(); ++itr) {
		  	PhySt & curSt = TrackStatusRef(itr);
		    MgntProp::Prop(vtxYZ(1) - curSt.Z(), curSt);
        Double_t dist = std::fabs(vtxYZ(0) - curSt.Y()) / std::sqrt(1. + curSt.TanY() * curSt.TanY());
        normDist += dist;
      }
      normDist /= NumOfTrack();
      if (normDist > Shift) return -26;
    }
  }

  // Require XZ & YZ is match
  if (std::fabs(vtxXZ(1) - vtxYZ(1)) > Shift) return -31;

  // linear fit for xy-dir
  MtxLB::SVec3D vtxXYZ;
  Bool_t LinearFitXY = true;
  if (LinearFitY) {
    const Int_t nIterLnXY = 1;
    for (Int_t iter = 1; iter <= nIterLnXY; ++iter) {
      MtxLB::SMtxSym3D coefXYZ;
      MtxLB::SVec3D    gradXYZ;
      for (Int_t itr = 0; itr < NumOfTrack(); ++itr) {
		  	PhySt & st = TrackStatusRef(itr);
        Double_t tx = st.TanX();
        Double_t ty = st.TanY();
        Double_t invXErr2 = 1. / (1. + tx * tx);
        Double_t invYErr2 = 1. / (1. + ty * ty);
        coefXYZ(0, 0) += invXErr2;
        coefXYZ(0, 2) += invXErr2 * (-1. * tx);
        coefXYZ(2, 2) += invXErr2 * (tx * tx);
        coefXYZ(1, 1) += invYErr2;
        coefXYZ(1, 2) += invYErr2 * (-1. * ty);
        coefXYZ(2, 2) += invYErr2 * (ty * ty);
        gradXYZ(0)    += invXErr2 * (st.X() - tx * st.Z());
        gradXYZ(2)    += invXErr2 * (-1. * tx) * (st.X() - tx * st.Z());
        gradXYZ(1)    += invYErr2 * (st.Y() - ty * st.Z());
        gradXYZ(2)    += invYErr2 * (-1. * ty) * (st.Y() - ty * st.Z());
      }
      Bool_t isSucs = coefXYZ.Invert();
      vtxXYZ = coefXYZ * gradXYZ;
      if (!isSucs) return -41;
      if (std::isinf(vtxXYZ(0)) || std::isnan(vtxXYZ(0))) return -42;
      if (std::isinf(vtxXYZ(1)) || std::isnan(vtxXYZ(1))) return -43;
      if (std::isinf(vtxXYZ(2)) || std::isnan(vtxXYZ(2))) return -44;

      //if (vtxXYZ(2) < Boundary[0] && vtxXYZ(2) > Boundary[1]) return -45;
      if (std::fabs(vtxXYZ(0)) > ExtLimit || std::fabs(vtxXYZ(1)) > ExtLimit || std::fabs(vtxXYZ(2)) > ExtLimit) return -46;

      Double_t normDist = 0;
      for (Int_t itr = 0; itr < NumOfTrack(); ++itr) {
		  	PhySt & curSt = TrackStatusRef(itr);
		    MgntProp::Prop(vtxXYZ(2) - curSt.Z(), curSt);
        Double_t distX = std::fabs(vtxXYZ(0) - curSt.X()) / std::sqrt(1. + curSt.TanX() * curSt.TanX());
        Double_t distY = std::fabs(vtxXYZ(1) - curSt.Y()) / std::sqrt(1. + curSt.TanY() * curSt.TanY());
        Double_t dist = std::sqrt(distX * distX + distY * distY);
        normDist += dist;
      }
      normDist /= NumOfTrack();
      if (normDist > Shift) return -47;
    }
  }

  // Geom Vertex
  fGeomVertex = vtxXYZ;

  // Require Down going
  //if (fGeomVertex(2) < Boundary[0]+Shift && fGeomVertex(2) > Boundary[1]-Shift) return -51;

  // simple fit
  fStatus.setSpatialPos(fGeomVertex(0), fGeomVertex(1), fGeomVertex(2));
  for (Int_t itr = 0; itr < NumOfTrack(); ++itr) {
  	PhySt & curSt = TrackStatusRef(itr);
    curSt.setSpatialPos(fGeomVertex(0), fGeomVertex(1), fGeomVertex(2));
  }

  Int_t NPars = NumOfTrack() * 3 + 3;
  const Int_t nIterSimple = 2;
  for (Int_t iter = 1; iter <= nIterSimple; ++iter) {
    //MtxLB::TMtxD Cov(NPars, NPars);
    //MtxLB::TVecD Grad(NPars);
    MtxLB::SMtxSym9D Cov;
    MtxLB::SVec9D    Grad;
    Double_t         Chisq =  0;
    Int_t            NDF   = -3;

    for (Int_t itr = 0; itr < NumOfTrack(); ++itr) {
      SMtxSymD<5> CovSub;
      SVecD<5>    GradSub;
      Double_t         ChisqSub =  0;
      Int_t            NDFSub   = -3;

      Track * track = GetTrack(itr);
      //PhySt curSt = TrackStatusRef(itr);
      //PhyJb curJacb(PhyJb::Identity);
      for (Int_t ih = 0; ih < track->NumOfHit(); ++ih) {
        PhySt curSt = TrackStatusRef(itr);
      	HitSt * hit = track->HitPtr(ih);

        PhyJb curJacb;
      	MgntProp::PropWithJacb(hit->Z() - curSt.Z(), curSt, curJacb);
        //PhyJb tmpJacb;
      	//MgntProp::PropWithJacb(hit->Z() - curSt.Z(), curSt, tmpJacb);
		    //MgntProp::PropWithJacb(hit->Z() - curSt.Z(), curSt, tmpJacb, MgntProp::kPropWithMaterial, &fMats);
      	//curJacb() = tmpJacb() * curJacb();
      	SMtxD<2, 5> subJacb = curJacb.getSubJacb();

      	SVecD<2> dst;
        dst(0) = (hit->SideX() ? curSt.X() - hit->X() : 0);
        dst(1) = (hit->SideY() ? curSt.Y() - hit->Y() : 0);

		    SMtxSymD<2> det;
		    det(0, 0) = (hit->ErrX()*hit->ErrX());
		    det(1, 1) = (hit->ErrY()*hit->ErrY());

        MaterialPhy matPhy = MgntMat::Calculate(track->Mat(), TrackStatusRef(itr), TrackStatusRef(itr).Z(), hit->Z());
        SMtxSymD<2> idetmat = (det + matPhy.fMatCov.Sub<SMtxSymD<2>>(0, 0));
		    idetmat.Invert();

        //SMtxSymD<2> idet;
      	//idet(0, 0) = 1. / (hit->ErrX()*hit->ErrX());
      	//idet(1, 1) = 1. / (hit->ErrY()*hit->ErrY());

        CovSub   += MtxLB::SimilarityT(subJacb, idetmat);
        GradSub  += MtxLB::Transpose(subJacb) * idetmat * dst;
        ChisqSub += MtxLB::Similarity(dst, idetmat);
      	NDFSub   += (hit->SideX() + hit->SideY());
      }
      //MtxLB::SMtxR5C6D trans = PhyJb::GetTransJacbDzDs(curSt); // problem : S6 = trans * S5
      MtxLB::SMtxR5C6D trans = PhyJb::GetTransJacbDzDs(TrackStatusRef(itr)); // problem : S6 = trans * S5
      MtxLB::SMtxSym6D CovTrans  = MtxLB::SimilarityT(trans, CovSub);
      MtxLB::SVec6D    GradTrans = MtxLB::Transpose(trans) * GradSub;

      for (Int_t it = 0; it < 3; ++it) {
        for (Int_t jt = 0; jt < 3; ++jt) {
          Cov(it, jt) += CovTrans(it, jt);
          Cov(it+itr*3, jt+itr*3+3) += CovTrans(it, jt+3);
          Cov(it+itr*3+3, jt+itr*3+3) += CovTrans(it+3, jt+3);
        }
        Grad(it) += GradTrans(it);
        Grad(it+3) += GradTrans(it+3);
      }
      Chisq += ChisqSub;
      NDF   += NDFSub;
    }
    Bool_t isSucs = Cov.Invert();
    if (!isSucs) return -61;
    MtxLB::SVec9D Sol = Cov * Grad;

    MtxLB::SVec4D momVec;
    fStatus.setSpatialPos(fStatus.X() - Sol(0), fStatus.Y() - Sol(1), fStatus.Z() - Sol(2));
    for (Int_t itr = 0; itr < NumOfTrack(); ++itr) {
      PhySt & curSt = TrackStatusRef(itr);
      curSt.setSpatialWithTan(
        fStatus.X(),
        fStatus.Y(),
        fStatus.Z(),
        curSt.TanX() - Sol(itr*3+3),
        curSt.TanY() - Sol(itr*3+4),
        curSt.DirZ()
      );
      Double_t InvR = curSt.InvRigWithPropFact() - Sol(itr*3+5);
      curSt.setInvRigWithPropFact(InvR);
      momVec += curSt.MomVec();
    }
    fStatus.setMomVec(momVec);

    fCov   = Cov;
    fChisq = Chisq;
    fNDF   = NDF;
  }

	return (fNDF > 0) ? (fChisq / fNDF) : 0;
}
*/
#endif // __MgntRecon_C__
