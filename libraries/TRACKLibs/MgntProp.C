#ifndef __MgntProp_C__
#define __MgntProp_C__
#include "MgntProp.h"


/********************************
 **    MgntProp based on DS    **
 ********************************/

// Step Length
// step < 0, backward trace
// step > 0, forward trace
Double_t MgntProp::GetPropStep(PhySt & phySt, Double_t sign) {
	MtxLB::SVecD<3> && coo    = phySt.Pos() + sign * MgntProp::INT_STEP * phySt.Dir();
	MtxLB::SVecD<3> &  magVec = MgntMag::GetMagFast(coo);
	Double_t           magMag = MtxLB::Mag(magVec);

	Double_t prp = MgntProp::PROP_FACT * std::fabs(phySt.ChrgMass()) / (std::fabs(phySt.Eta()) + TUNE_ETA);
	Double_t cos = MtxLB::Dot(phySt.Dir(), magVec);
	Double_t sin = std::sqrt(1.001 - cos * cos / magMag / magMag);
	Double_t chr = prp * (0.2 + 0.8 * sin); 

	Double_t length = (MgntProp::TUNE_STEP / chr /  (magMag + MgntProp::TUNE_MAGF));
	if (!MgntNum::Valid(length))     length = MgntProp::MAX_STEP;
	if (length > MgntProp::MAX_STEP) length = MgntProp::MAX_STEP;
	if (length < MgntProp::MIN_STEP) length = MgntProp::MIN_STEP;

	if (MgntProp::CheckEnvironment(MgntProp::kMaterial)) {
		MtxLB::SVecD<3> finCoo = phySt.Pos() + (sign * length) * phySt.Dir();
		Double_t rat = MgntMatPhyCal::GetSimplicityNumRadLen(phySt.Pos(), finCoo) / MgntProp::TUNE_MAT;
		if (MgntNum::Valid(rat) && rat > ONE) length /= rat;
	}

	return length;
}

Double_t MgntProp::GetStep(PhySt & phySt, Double_t resStep) {
	Double_t sign = MgntNum::Compare(resStep);
	Double_t len = MgntProp::GetPropStep(phySt, sign);
	Double_t res = std::fabs(resStep);

	Double_t length = ZERO;
	if      (res < 1.2 * len) length = res;
	else if (res < 1.7 * len) length = 0.7 * len;
	else                      length = len;

	Double_t step = sign * length;

	if (MgntProp::CheckEnvironment(MgntProp::kMaterial)) {
		MtxLB::SVecD<3>    satCoo = phySt.Pos();
		MtxLB::SVecD<3> && endCoo = satCoo + step * phySt.Dir();
		MatCal = MgntMatPhyCal::GetMatPhyCalParam(phySt, satCoo, endCoo);
	}

	return step;
}

Double_t MgntProp::GetStepToZ(PhySt & phySt, Double_t resStepZ) {
	Double_t signDz = MgntNum::Compare(resStepZ);
	Double_t sign = MgntNum::Compare(phySt.DirZ() * resStepZ);
	Double_t lens = sign * MgntProp::GetPropStep(phySt, sign);
	Double_t resz = std::fabs(resStepZ);

	MtxLB::SVecD<3> & magVec = MgntMag::GetMagFast();
	Double_t propFT = MgntProp::PROP_FACT * std::fabs(phySt.ChrgMass());
	Double_t crsUBZ = (phySt.DirX() * magVec(1) - phySt.DirY() * magVec(0));
	Double_t termS0 = 0.5 * propFT * phySt.InvEta() * crsUBZ;
	Double_t termS1 = phySt.DirZ();
	Double_t lenz = std::fabs(termS0 * lens * lens + termS1 * lens);

	Double_t lengthz = 0;
	if      (resz < 1.2 * lenz) lengthz = resz;
	else if (resz < 1.7 * lenz) lengthz = 0.7 * lenz;
	else                        lengthz = lenz;
	Double_t termS2 = NEG * signDz * lengthz;

	// solve step
	Double_t step = ZERO;
	if (MgntNum::EqualToZero(termS0)) {
		step = NEG * termS2 / termS1;
	}
	else {
		Double_t discriminant = std::sqrt(termS1 * termS1 - 4. * termS0 * termS2); 
		Double_t solve1 = HALF * (NEG * termS1 + discriminant) / termS0;
		Double_t solve2 = HALF * (NEG * termS1 - discriminant) / termS0;
		Bool_t isSol1 = (!MgntNum::EqualToZero(solve1) && MgntNum::Valid<Double_t>(solve1 / lens) && std::fabs(solve1 / lens) < 10.);
		Bool_t isSol2 = (!MgntNum::EqualToZero(solve2) && MgntNum::Valid<Double_t>(solve2 / lens) && std::fabs(solve2 / lens) < 10.);
		if      (isSol1 && MgntNum::Compare(solve1 * sign) > 0) step = solve1;
		else if (isSol2 && MgntNum::Compare(solve2 * sign) > 0) step = solve2;
		else                                                    step = NEG * termS2 / termS1;
	}
	if (!MgntNum::Valid(step)) step = lens;

	if (MgntProp::CheckEnvironment(MgntProp::kMaterial)) {
		MtxLB::SVecD<3> satCoo = phySt.Pos();
		MtxLB::SVecD<3> endCoo = satCoo + step * phySt.Dir();
		MatCal = MgntMatPhyCal::GetMatPhyCalParam(phySt, satCoo, endCoo);
	}

	return step;
}


// Prop with Numeric Method
Bool_t MgntProp::Prop(Double_t step, PhySt & phySt, PhyJb * phyJb, MatPropParam * matPar) {
	if (matPar == nullptr || matPar->Vacuum()) { MgntProp::SetEnvironment(MgntProp::kVacuum); MatPar.init(); MatCal.init(); }
	else { MgntProp::SetEnvironment(MgntProp::kMaterial); MatPar = *matPar; MatCal.init(); }
	Bool_t WithJacb = (phyJb != nullptr);
	if (WithJacb) phyJb->init(PhyJb::Identity);

	Long64_t iter = 0;
	Bool_t   isSuccess = true;
	Double_t currentStep = 0.;
	Double_t chrgSign = MgntNum::Compare(phySt.Eta());
	while (MgntNum::Valid(step - currentStep) && std::fabs(step - currentStep) > MgntProp::MIN_CONT && iter < MgntProp::MAX_ITER) {
		Double_t residualStep = step - currentStep;
		Double_t subStep = GetStep(phySt, residualStep);

		Bool_t vaild = true;;
		PhyJb * subPhyJb = (!WithJacb) ? 0 : (new PhyJb());
		switch (Method) {
			case kEuler :
				vaild = Prop_EulerMethod(subStep, phySt, subPhyJb); break;
			case kEulerHeun :
				vaild = Prop_EulerHeunMethod(subStep, phySt, subPhyJb); break;
			case kRungeKuttaNystrom :
				vaild = Prop_RungeKuttaNystromMethod(subStep, phySt, subPhyJb); break;
			default :
				std::cerr << "\nERROR : MgntProp method is error.\n";
		}

		iter++;
		currentStep += subStep;
		if (WithJacb) { phyJb->multiplied(*subPhyJb); delete subPhyJb; subPhyJb = nullptr; }
		if (MgntNum::Compare(chrgSign * phySt.Eta()) < 0 || phySt.Mom() < MgntProp::MIN_MOM) { isSuccess = false; break; }
		if (!vaild) { isSuccess = false; break; }
	}
	
	if (iter >= MgntProp::MAX_ITER) isSuccess = false;
	if (!isSuccess) {
		phySt.init(phySt.Part()); 
		if (WithJacb) phyJb->init(PhyJb::Identity);
	}
	return isSuccess;
}
		
Bool_t MgntProp::PropToZ(Double_t zcoo, PhySt & phySt, PhyJb * phyJb, MatPropParam * matPar) {
	if (MgntNum::EqualToZero(phySt.DirZ())) return false;
	if (matPar == nullptr || matPar->Vacuum()) { MgntProp::SetEnvironment(MgntProp::kVacuum); MatPar.init(); MatCal.init(); }
	else { MgntProp::SetEnvironment(MgntProp::kMaterial); MatPar = *matPar; MatCal.init(); }
	Bool_t WithJacb = (phyJb != nullptr);
	if (WithJacb) phyJb->init(PhyJb::Identity);

	Long64_t iter = 0;
	Bool_t   isSuccess = true;
	Double_t chrgSign = MgntNum::Compare(phySt.Eta());
	while (MgntNum::Valid(zcoo - phySt.Z()) && std::fabs(zcoo - phySt.Z()) > MgntProp::MIN_CONT && iter < MgntProp::MAX_ITER) {
		Double_t residualStepZ = zcoo - phySt.Z();
		Double_t subStep = GetStepToZ(phySt, residualStepZ);

		Bool_t vaild = true;;
		PhyJb * subPhyJb = (!WithJacb) ? 0 : (new PhyJb());
		switch (Method) {
			case kEuler :
				vaild = Prop_EulerMethod(subStep, phySt, subPhyJb); break;
			case kEulerHeun :
				vaild = Prop_EulerHeunMethod(subStep, phySt, subPhyJb); break;
			case kRungeKuttaNystrom :
				vaild = Prop_RungeKuttaNystromMethod(subStep, phySt, subPhyJb); break;
			default :
				std::cerr << "\nERROR : MgntProp method is error.\n";
		}

		iter++;
		if (WithJacb) { phyJb->multiplied(*subPhyJb); delete subPhyJb; subPhyJb = nullptr; }
		if (MgntNum::Compare(chrgSign * phySt.Eta()) < 0 || phySt.Mom() < MgntProp::MIN_MOM) { isSuccess = false; break; }
		if (!vaild) { isSuccess = false; break; }
	}
	
	if (iter >= MgntProp::MAX_ITER) isSuccess = false;
	if (!isSuccess) {
		phySt.print(); // test
		phySt.init(phySt.Part()); 
		if (WithJacb) phyJb->init(PhyJb::Identity);
	}
	return isSuccess;
}

// Prop with Monte Carlo
Bool_t MgntProp::Prop_MonteCarlo(Double_t step, PhySt & phySt, Bool_t optMscat, Bool_t optEngls) { 
	if (optMscat || optEngls) MgntProp::SetEnvironment(MgntProp::kMaterial);
	else MgntProp::SetEnvironment(MgntProp::kVacuum);
	MatCal.init(); MatPar.init();
	
	Long64_t iter = 0;
	Bool_t   isSuccess = true;
	Double_t currentStep = 0.;
	Double_t chrgSign = MgntNum::Compare(phySt.Eta());
	while (MgntNum::Valid(step - currentStep) && std::fabs(step - currentStep) > MgntProp::MIN_CONT && iter < MgntProp::MAX_ITER) {
		Double_t residualStep = step - currentStep;
		Double_t subStep = GetStep(phySt, residualStep);
		MatPar.random(MatCal, optMscat, optEngls);

		Bool_t vaild = true;;
		switch (Method) {
			case kEuler :
				vaild = Prop_EulerMethod(subStep, phySt); break;
			case kEulerHeun :
				vaild = Prop_EulerHeunMethod(subStep, phySt); break;
			case kRungeKuttaNystrom :
				vaild = Prop_RungeKuttaNystromMethod(subStep, phySt); break;
			default :
				std::cerr << "\nERROR : MgntProp method is error.\n";
		}

		iter++;
		currentStep += subStep;
		if (MgntNum::Compare(chrgSign * phySt.Eta()) < 0 || phySt.Mom() < MgntProp::MIN_MOM) { isSuccess = false; break; }
		if (!vaild) { isSuccess = false; break; }
	}
	
	if (iter >= MgntProp::MAX_ITER) isSuccess = false;
	if (!isSuccess) phySt.init(phySt.Part()); 
	return isSuccess;
}

Bool_t MgntProp::PropToZ_MonteCarlo(Double_t zcoo, PhySt & phySt, Bool_t optMscat, Bool_t optEngls) {
	if (MgntNum::EqualToZero(phySt.DirZ())) return false;
	if (optMscat || optEngls) MgntProp::SetEnvironment(MgntProp::kMaterial);
	else MgntProp::SetEnvironment(MgntProp::kVacuum);
	MatCal.init(); MatPar.init();
	
	Long64_t iter = 0;
	Bool_t   isSuccess = true;
	Double_t chrgSign = MgntNum::Compare(phySt.Eta());
	while (MgntNum::Valid(zcoo - phySt.Z()) && std::fabs(zcoo - phySt.Z()) > MgntProp::MIN_CONT && iter < MgntProp::MAX_ITER) {
		Double_t residualStepZ = zcoo - phySt.Z();
		Double_t subStep = GetStepToZ(phySt, residualStepZ);
		MatPar.random(MatCal, optMscat, optEngls);

		Bool_t vaild = true;;
		switch (Method) {
			case kEuler :
				vaild = Prop_EulerMethod(subStep, phySt); break;
			case kEulerHeun :
				vaild = Prop_EulerHeunMethod(subStep, phySt); break;
			case kRungeKuttaNystrom :
				vaild = Prop_RungeKuttaNystromMethod(subStep, phySt); break;
			default :
				std::cerr << "\nERROR : MgntProp method is error.\n";
		}
		
		iter++;
		if (MgntNum::Compare(chrgSign * phySt.Eta()) < 0 || phySt.Mom() < MgntProp::MIN_MOM) { isSuccess = false; break; }
		if (!vaild) { isSuccess = false; break; }
	}
	
	if (iter >= MgntProp::MAX_ITER) isSuccess = false;
	if (!isSuccess) phySt.init(phySt.Part()); 
	return isSuccess;
}


// Official
Bool_t MgntProp::PropToZ_Official(Double_t zcoo, PhySt & phySt) {
	if (MgntNum::EqualToZero(phySt.DirZ())) return false;
	if (MgntNum::Valid(zcoo - phySt.Z()) && std::fabs(zcoo - phySt.Z()) < MgntProp::MIN_CONT) return true;
	MgntProp::SetEnvironment(MgntProp::kVacuum);
	MatCal.init(); MatPar.init();

	Double_t sign = MgntNum::Compare(phySt.DirZ());
	AMSPoint pos(phySt.X(),    phySt.Y(),    phySt.Z());
	AMSDir   dir(phySt.DirX(), phySt.DirY(), phySt.DirZ());
	Double_t rig = (-sign) * phySt.Rig();
	
	TrProp trProp(pos, dir, rig);
	Double_t len = trProp.Propagate(zcoo);
	if (len < 0) { phySt.init(phySt.Part()); return false; }

	phySt.setSpatialWithCos(
		trProp.GetP0x(),
		trProp.GetP0y(),
		trProp.GetP0z(),
		trProp.GetD0x() * sign,
		trProp.GetD0y() * sign,
		trProp.GetD0z() * sign
	);
	Bool_t vaild = phySt.vaild();

	return vaild;
}
	

// Analytic
Bool_t MgntProp::PropToZ_Analytic(Double_t zcoo, PhySt & phySt, MtxLB::SMtxD<5> * phyJb) {
	if (MgntNum::EqualToZero(phySt.DirZ())) return false;
	if (MgntNum::EqualToZero(phySt.InvEta())) return false;
	if (MgntNum::Valid(zcoo - phySt.Z()) && MgntNum::Compare(std::fabs(zcoo - phySt.Z()), MgntProp::MIN_CONT) < 0) {
		if (phyJb != nullptr) (*phyJb) = MtxLB::SIdMtx();
		return true;
	}
	MgntProp::SetEnvironment(MgntProp::kVacuum);
	MatCal.init(); MatPar.init();

	Double_t stpz = zcoo - phySt.Z();
	Double_t chms = std::fabs(phySt.ChrgMass());
	Double_t ieta = phySt.InvEta();
	Double_t dsdz = std::fabs(ONE / phySt.DirZ());
	Double_t sign = MgntNum::Compare(phySt.DirZ());
	Double_t prop = MgntProp::PROP_FACT * chms * sign;
	Double_t posx = phySt.X();
	Double_t posy = phySt.Y();
	Double_t tanx = phySt.TanX();
	Double_t tany = phySt.TanY();
	Double_t prpx = prop * dsdz * (tanx * tany);
	Double_t prpy = prop * dsdz * (tany * tany + ONE);

	// Integral Magnetic Field
	const Double_t stableFT = 1.0;
	Int_t nStpz = std::ceil(std::fabs(stableFT * stpz / MgntProp::MAX_STEP));
	Double_t lStpz = stpz / Double_t(nStpz); 

	Double_t stpzSIGN = MgntNum::Compare(stpz);
	Double_t magxBasedINT1 = MgntMag::GetMagFuncINT1(phySt.Z());
	Double_t magxINT1 = ZERO;
	for (Int_t iStpz = 0; iStpz < nStpz; ++iStpz) {
		Double_t ztmp = phySt.Z() + lStpz * (Double_t(iStpz) + 0.5);
		magxINT1 += (MgntMag::GetMagFuncINT1(ztmp) - magxBasedINT1);
	}

	Double_t magxINT1_POS = stpzSIGN * std::fabs(magxINT1 / Double_t(nStpz));
	Double_t magxINT1_TAN = stpzSIGN * std::fabs(MgntMag::GetMagFuncINT1(zcoo) - magxBasedINT1);
	
	// Final Status
	Double_t rslPosX = posx + stpz * (tanx + prpx * magxINT1_POS * ieta);
	Double_t rslPosY = posy + stpz * (tany + prpy * magxINT1_POS * ieta);
	Double_t rslTanX = tanx + prpx * magxINT1_TAN * ieta;
	Double_t rslTanY = tany + prpy * magxINT1_TAN * ieta;

	phySt.setSpatialWithTan(rslPosX, rslPosY, zcoo, rslTanX, rslTanY, phySt.DirZ());
	Bool_t vaild = phySt.vaild();

	if (phyJb != nullptr && vaild) {
		Double_t tanx2 = tanx * tanx;
		Double_t tany2 = tany * tany;
		Double_t xtanDevX = (tany * tany2 + TWO * tany * tanx2 + tany) / dsdz;
		Double_t xtanDevY = (tanx * tanx2 + TWO * tanx * tany2 + tanx) / dsdz;
		Double_t ytanDevX = tanx * (tany2 + ONE) / dsdz;
		Double_t ytanDevY = tany * (TWO * tanx2 + THREE * tany2 + THREE) / dsdz;

		(*phyJb) = MtxLB::SIdMtx();
		
		(*phyJb)(0, 2) += stpz * (ONE + prop * xtanDevX * magxINT1_POS * ieta); 
		(*phyJb)(0, 3) += stpz * (prop * xtanDevY * magxINT1_POS * ieta);
		(*phyJb)(0, 4) += stpz * prpx * magxINT1_POS;

		(*phyJb)(1, 2) += stpz * (prop * ytanDevX * magxINT1_POS * ieta);
		(*phyJb)(1, 3) += stpz * (ONE + prop * ytanDevY * magxINT1_POS * ieta);
		(*phyJb)(1, 4) += stpz * prpy * magxINT1_POS;

		(*phyJb)(2, 2) += prop * xtanDevX * magxINT1_TAN * ieta;
		(*phyJb)(2, 3) += prop * xtanDevY * magxINT1_TAN * ieta;
		(*phyJb)(2, 4) += prpx * magxINT1_TAN;

		(*phyJb)(3, 2) += prop * ytanDevX * magxINT1_TAN * ieta;
		(*phyJb)(3, 3) += prop * ytanDevY * magxINT1_TAN * ieta;
		(*phyJb)(3, 4) += prpy * magxINT1_TAN;
	}

	return vaild;
}

// Euler Method
Bool_t MgntProp::Prop_EulerMethod(Double_t step, PhySt & phySt, PhyJb * phyJb) {
	Double_t halfSquareStep = MgntProp::HALF * step * step;
	Bool_t WithJacb = (phyJb != nullptr);

	PhySt St0(phySt);
	DevStatusFunc Zeta0(St0, &MatPar, &MatCal);

	phySt.setSpatialWithCos(
		phySt.X()    + step * Zeta0(PX) + halfSquareStep * Zeta0(UX),
		phySt.Y()    + step * Zeta0(PY) + halfSquareStep * Zeta0(UY),
		phySt.Z()    + step * Zeta0(PZ) + halfSquareStep * Zeta0(UZ),
		phySt.DirX() + step * Zeta0(UX),
		phySt.DirY() + step * Zeta0(UY),
    phySt.DirZ() + step * Zeta0(UZ)
	);
	phySt.setInvEta( phySt.InvEta() + step * Zeta0(GB) );
	Bool_t vaild = phySt.vaild();

	if (WithJacb && vaild) {
		DevParamFunc Ki0(St0, &MatPar, &MatCal);
		
		phyJb->init(PhyJb::Identity, MatCal.Vacuum());
		phyJb->DeltaZ() = (phySt.Z() - St0.Z());
		phyJb->NumRadLen() = MatCal.NumRadLen();
		phyJb->G(JPX, JUX) += step * Ki0.G(JPX, JUX);
		phyJb->G(JPY, JUY) += step * Ki0.G(JPY, JUY);
		for (Int_t it = JUX; it < PhyJb::GDim; ++it) {
			phyJb->G(JPX, it) += halfSquareStep * Ki0.G(JUX, it);
			phyJb->G(JPY, it) += halfSquareStep * Ki0.G(JUY, it);
			phyJb->G(JUX, it) += step * Ki0.G(JUX, it);
			phyJb->G(JUY, it) += step * Ki0.G(JUY, it);
			phyJb->G(JGB, it) += step * Ki0.G(JGB, it);
		}
		for (Int_t it = 0; it < PhyJb::LDim && !phyJb->Vacuum(); ++it) {
			phyJb->L(JPX, it) += step * Ki0.L(JPX, it) + halfSquareStep * Ki0.L(JUX, it);
			phyJb->L(JPY, it) += step * Ki0.L(JPY, it) + halfSquareStep * Ki0.L(JUY, it);
			phyJb->L(JUX, it) += step * Ki0.L(JUX, it);
			phyJb->L(JUY, it) += step * Ki0.L(JUY, it);
			phyJb->L(JGB, it) += step * Ki0.L(JGB, it);
		}
	}
	
	return vaild;
}


// Euler-Heun Method
Bool_t MgntProp::Prop_EulerHeunMethod(Double_t step, PhySt & phySt, PhyJb * phyJb) {
	Double_t halfStep = MgntProp::HALF * step;
	Double_t oneSixthStep = MgntProp::ONE_SIXTH * step;
	Double_t halfSquareStep = MgntProp::HALF * step * step;
	Double_t oneSixthSquareStep = MgntProp::ONE_SIXTH * step * step;
	Bool_t WithJacb = (phyJb != nullptr);
	
	PhySt St0(phySt);
	DevStatusFunc Zeta0(St0, &MatPar, &MatCal);

	PhySt St1(phySt);
	St1.setSpatialWithCos(
		phySt.X()    + step * Zeta0(PX) + halfSquareStep * Zeta0(UX),
		phySt.Y()    + step * Zeta0(PY) + halfSquareStep * Zeta0(UY),
		phySt.Z()    + step * Zeta0(PZ) + halfSquareStep * Zeta0(UZ),
		phySt.DirX() + step * Zeta0(UX),
		phySt.DirY() + step * Zeta0(UY),
    phySt.DirZ() + step * Zeta0(UZ)
	);
	St1.setInvEta( phySt.InvEta() + step * Zeta0(GB) );

	DevStatusFunc Zeta1(St1, &MatPar, &MatCal);

	phySt.setSpatialWithCos(
		phySt.X()    + step * Zeta0(PX) + oneSixthSquareStep * (TWO * Zeta0(UX) + Zeta1(UX)),
		phySt.Y()    + step * Zeta0(PY) + oneSixthSquareStep * (TWO * Zeta0(UY) + Zeta1(UY)),
		phySt.Z()    + step * Zeta0(PZ) + oneSixthSquareStep * (TWO * Zeta0(UZ) + Zeta1(UZ)),
		phySt.DirX() + halfStep * (Zeta0(UX) + Zeta1(UX)),
		phySt.DirY() + halfStep * (Zeta0(UY) + Zeta1(UY)),
    phySt.DirZ() + halfStep * (Zeta0(UZ) + Zeta1(UZ))
	);
	phySt.setInvEta( phySt.InvEta() + halfStep * (Zeta0(GB) + Zeta1(GB)) );
	Bool_t vaild = phySt.vaild();

	if (WithJacb && vaild) {
		DevParamFunc  Ki0(St0, &MatPar, &MatCal);
		DevParamFunc  Ki1(St1, &MatPar, &MatCal);

		PhyJb Jb1(PhyJb::Identity, MatCal.Vacuum());
		Jb1.DeltaZ() = (St1.Z() - St0.Z());
		Jb1.G(JPX, JUX) += step * Ki0.G(JPX, JUX);
		Jb1.G(JPY, JUY) += step * Ki0.G(JPY, JUY);
		for (Int_t it = JUX; it < PhyJb::GDim; ++it) {
			Jb1.G(JPX, it) += halfSquareStep * Ki0.G(JUX, it);
			Jb1.G(JPY, it) += halfSquareStep * Ki0.G(JUY, it);
			Jb1.G(JUX, it) += step * Ki0.G(JUX, it);
			Jb1.G(JUY, it) += step * Ki0.G(JUY, it);
			Jb1.G(JGB, it) += step * Ki0.G(JGB, it);
		}
		for (Int_t it = 0; it < PhyJb::LDim && !Jb1.Vacuum(); ++it) {
			Jb1.L(JPX, it) += step * Ki0.L(JPX, it) + halfSquareStep * Ki0.L(JUX, it);
			Jb1.L(JPY, it) += step * Ki0.L(JPY, it) + halfSquareStep * Ki0.L(JUY, it);
			Jb1.L(JUX, it) += step * Ki0.L(JUX, it);
			Jb1.L(JUY, it) += step * Ki0.L(JUY, it);
			Jb1.L(JGB, it) += step * Ki0.L(JGB, it);
		}
		MtxLB::SMtxD<5>    && KJb1G = Ki1.G() * Jb1.G();
		MtxLB::SMtxD<5, 4> && KJb1L = (MatCal.Vacuum()) ? MtxLB::SMtxD<5, 4>() : (Ki1.G() * Jb1.L() + Ki1.L());

		phyJb->init(PhyJb::Identity, MatCal.Vacuum());
		phyJb->DeltaZ() = (phySt.Z() - St0.Z());
		phyJb->NumRadLen() = MatCal.NumRadLen();
		phyJb->G(JPX, JUX) += step * Ki0.G(JPX, JUX);
		phyJb->G(JPY, JUY) += step * Ki0.G(JPY, JUY);
		for (Int_t it = JUX; it < PhyJb::GDim; ++it) {
			phyJb->G(JPX, it) += oneSixthSquareStep * (MgntProp::TWO * Ki0.G(JUX, it) + KJb1G(JUX, it));
			phyJb->G(JPY, it) += oneSixthSquareStep * (MgntProp::TWO * Ki0.G(JUY, it) + KJb1G(JUY, it));
			phyJb->G(JUX, it) += halfStep * (Ki0.G(JUX, it) + KJb1G(JUX, it));
			phyJb->G(JUY, it) += halfStep * (Ki0.G(JUY, it) + KJb1G(JUY, it));
			phyJb->G(JGB, it) += halfStep * (Ki0.G(JGB, it) + KJb1G(JGB, it));
		}	
		for (Int_t it = 0; it < PhyJb::LDim && !phyJb->Vacuum(); ++it) {
			phyJb->L(JPX, it) += step * Ki0.L(JPX, it) + oneSixthSquareStep * (MgntProp::TWO * Ki0.L(JUX, it) + KJb1L(JUX, it));
			phyJb->L(JPY, it) += step * Ki0.L(JPY, it) + oneSixthSquareStep * (MgntProp::TWO * Ki0.L(JUY, it) + KJb1L(JUY, it));
			phyJb->L(JUX, it) += halfStep * (Ki0.L(JUX, it) + KJb1L(JUX, it));
			phyJb->L(JUY, it) += halfStep * (Ki0.L(JUY, it) + KJb1L(JUY, it));
			phyJb->L(JGB, it) += halfStep * (Ki0.L(JGB, it) + KJb1L(JGB, it));
		}
	}

	return vaild;
}


// Runge-Kutta-Nystrom Method
Bool_t MgntProp::Prop_RungeKuttaNystromMethod(Double_t step, PhySt & phySt, PhyJb * phyJb) {
	Double_t halfStep = MgntProp::HALF * step;
	Double_t oneSixthStep = MgntProp::ONE_SIXTH * step;
	Double_t halfSquareStep = MgntProp::HALF * step * step;
	Double_t oneSixthSquareStep = MgntProp::ONE_SIXTH * step * step;
	Double_t oneEighthSquareStep = MgntProp::ONE_EIGHTH * step * step;
	Bool_t WithJacb = (phyJb != nullptr);

	PhySt St0(phySt);
	DevStatusFunc Zeta0(St0, &MatPar, &MatCal);

	PhySt St1(phySt);
	St1.setSpatialWithCos(
		phySt.X()    + halfStep * Zeta0(PX) + oneEighthSquareStep * Zeta0(UX),
		phySt.Y()    + halfStep * Zeta0(PY) + oneEighthSquareStep * Zeta0(UY),
		phySt.Z()    + halfStep * Zeta0(PZ) + oneEighthSquareStep * Zeta0(UZ),
		phySt.DirX() + halfStep * Zeta0(UX),
		phySt.DirY() + halfStep * Zeta0(UY),
    phySt.DirZ() + halfStep * Zeta0(UZ)
	);
	St1.setInvEta( phySt.InvEta() + halfStep * Zeta0(GB) );
	
	DevStatusFunc Zeta1(St1, &MatPar, &MatCal);

	PhySt St2(phySt);
	St2.setSpatialWithCos(
		phySt.X()    + halfStep * Zeta0(PX) + oneEighthSquareStep * Zeta0(UX),
		phySt.Y()    + halfStep * Zeta0(PY) + oneEighthSquareStep * Zeta0(UY),
		phySt.Z()    + halfStep * Zeta0(PZ) + oneEighthSquareStep * Zeta0(UZ),
		phySt.DirX() + halfStep * Zeta1(UX),
		phySt.DirY() + halfStep * Zeta1(UY),
    phySt.DirZ() + halfStep * Zeta1(UZ)
	);
	St2.setInvEta( phySt.InvEta() + halfStep * Zeta1(GB) );
	
	DevStatusFunc Zeta2(St2, &MatPar, &MatCal);
	
	PhySt St3(phySt);
	St3.setSpatialWithCos(
		phySt.X()    + step * Zeta0(PX) + halfSquareStep * Zeta2(UX),
		phySt.Y()    + step * Zeta0(PY) + halfSquareStep * Zeta2(UY),
		phySt.Z()    + step * Zeta0(PZ) + halfSquareStep * Zeta2(UZ),
		phySt.DirX() + step * Zeta2(UX),
		phySt.DirY() + step * Zeta2(UY),
    phySt.DirZ() + step * Zeta2(UZ)
	);
	St3.setInvEta( phySt.InvEta() + step * Zeta2(GB) );
	
	DevStatusFunc Zeta3(St3, &MatPar, &MatCal);
	
	phySt.setSpatialWithCos(
		phySt.X()    + step * Zeta0(PX) + oneSixthSquareStep * (Zeta0(UX) + Zeta1(UX) + Zeta2(UX)),
		phySt.Y()    + step * Zeta0(PY) + oneSixthSquareStep * (Zeta0(UY) + Zeta1(UY) + Zeta2(UY)),
		phySt.Z()    + step * Zeta0(PZ) + oneSixthSquareStep * (Zeta0(UZ) + Zeta1(UZ) + Zeta2(UZ)),
		phySt.DirX() + oneSixthStep * (Zeta0(UX) + MgntProp::TWO * Zeta1(UX) + MgntProp::TWO * Zeta2(UX) + Zeta3(UX)),
		phySt.DirY() + oneSixthStep * (Zeta0(UY) + MgntProp::TWO * Zeta1(UY) + MgntProp::TWO * Zeta2(UY) + Zeta3(UY)),
    phySt.DirZ() + oneSixthStep * (Zeta0(UZ) + MgntProp::TWO * Zeta1(UZ) + MgntProp::TWO * Zeta2(UZ) + Zeta3(UZ))
	);
	phySt.setInvEta( phySt.InvEta() + oneSixthStep * (Zeta0(GB) + MgntProp::TWO * Zeta1(GB) + MgntProp::TWO * Zeta2(GB) + Zeta3(GB)) );
	Bool_t vaild = phySt.vaild();

	if (WithJacb && vaild) {
		DevParamFunc  Ki0(St0, &MatPar, &MatCal);
		DevParamFunc  Ki1(St1, &MatPar, &MatCal);
		DevParamFunc  Ki2(St2, &MatPar, &MatCal);
		DevParamFunc  Ki3(St3, &MatPar, &MatCal);

		PhyJb Jb1(PhyJb::Identity, MatCal.Vacuum());
		Jb1.DeltaZ() = (St1.Z() - St0.Z());
		Jb1.G(JPX, JUX) += halfStep * Ki0.G(JPX, JUX);
		Jb1.G(JPY, JUY) += halfStep * Ki0.G(JPY, JUY);
		for (Int_t it = JUX; it < PhyJb::GDim; ++it) {
			Jb1.G(JPX, it) += oneEighthSquareStep * Ki0.G(JUX, it);
			Jb1.G(JPY, it) += oneEighthSquareStep * Ki0.G(JUY, it);
			Jb1.G(JUX, it) += halfStep * Ki0.G(JUX, it);
			Jb1.G(JUY, it) += halfStep * Ki0.G(JUY, it);
			Jb1.G(JGB, it) += halfStep * Ki0.G(JGB, it);
		}
		for (Int_t it = 0; it < PhyJb::LDim && !Jb1.Vacuum(); ++it) {
			Jb1.L(JPX, it) += halfStep * Ki0.L(JPX, it) + oneEighthSquareStep * Ki0.L(JUX, it);
			Jb1.L(JPY, it) += halfStep * Ki0.L(JPY, it) + oneEighthSquareStep * Ki0.L(JUY, it);
			Jb1.L(JUX, it) += halfStep * Ki0.L(JUX, it);
			Jb1.L(JUY, it) += halfStep * Ki0.L(JUY, it);
			Jb1.L(JGB, it) += halfStep * Ki0.L(JGB, it);
		}
		MtxLB::SMtxD<5>    && KJb1G = Ki1.G() * Jb1.G();
		MtxLB::SMtxD<5, 4> && KJb1L = (MatCal.Vacuum()) ? MtxLB::SMtxD<5, 4>() : (Ki1.G() * Jb1.L() + Ki1.L());
		
		PhyJb Jb2(PhyJb::Identity, MatCal.Vacuum());
		Jb2.DeltaZ() = (St2.Z() - St0.Z());
		Jb2.G(JPX, JUX) += halfStep * Ki0.G(JPX, JUX);
		Jb2.G(JPY, JUY) += halfStep * Ki0.G(JPY, JUY);
		for (Int_t it = JUX; it < PhyJb::GDim; ++it) {
			Jb2.G(JPX, it) += oneEighthSquareStep * Ki0.G(JUX, it);
			Jb2.G(JPY, it) += oneEighthSquareStep * Ki0.G(JUY, it);
			Jb2.G(JUX, it) += halfStep * KJb1G(JUX, it);
			Jb2.G(JUY, it) += halfStep * KJb1G(JUY, it);
			Jb2.G(JGB, it) += halfStep * KJb1G(JGB, it);
		}
		for (Int_t it = 0; it < PhyJb::LDim && !Jb2.Vacuum(); ++it) {
			Jb2.L(JPX, it) += halfStep * Ki0.L(JPX, it) + oneEighthSquareStep * Ki0.L(JUX, it);
			Jb2.L(JPY, it) += halfStep * Ki0.L(JPY, it) + oneEighthSquareStep * Ki0.L(JUY, it);
			Jb2.L(JUX, it) += halfStep * KJb1L(JUX, it);
			Jb2.L(JUY, it) += halfStep * KJb1L(JUY, it);
			Jb2.L(JGB, it) += halfStep * KJb1L(JGB, it);
		}
		MtxLB::SMtxD<5>    && KJb2G = Ki2.G() * Jb2.G();
		MtxLB::SMtxD<5, 4> && KJb2L = (MatCal.Vacuum()) ? MtxLB::SMtxD<5, 4>() : (Ki2.G() * Jb2.L() + Ki2.L());
		
		PhyJb Jb3(PhyJb::Identity, MatCal.Vacuum());
		Jb3.DeltaZ() = (St3.Z() - St0.Z());
		Jb3.G(JPX, JUX) += step * Ki0.G(JPX, JUX);
		Jb3.G(JPY, JUY) += step * Ki0.G(JPY, JUY);
		for (Int_t it = JUX; it < PhyJb::GDim; ++it) {
			Jb3.G(JPX, it) += halfSquareStep * KJb2G(JUX, it);
			Jb3.G(JPY, it) += halfSquareStep * KJb2G(JUY, it);
			Jb3.G(JUX, it) += step * KJb2G(JUX, it);
			Jb3.G(JUY, it) += step * KJb2G(JUY, it);
			Jb3.G(JGB, it) += step * KJb2G(JGB, it);
		}
		for (Int_t it = 0; it < PhyJb::LDim && !Jb3.Vacuum(); ++it) {
			Jb3.L(JPX, it) += step * Ki0.L(JPX, it) + halfSquareStep * KJb2L(JUX, it);
			Jb3.L(JPY, it) += step * Ki0.L(JPY, it) + halfSquareStep * KJb2L(JUY, it);
			Jb3.L(JUX, it) += step * KJb2L(JUX, it);
			Jb3.L(JUY, it) += step * KJb2L(JUY, it);
			Jb3.L(JGB, it) += step * KJb2L(JGB, it);
		}
		MtxLB::SMtxD<5>    && KJb3G = Ki3.G() * Jb3.G();
		MtxLB::SMtxD<5, 4> && KJb3L = (MatCal.Vacuum()) ? MtxLB::SMtxD<5, 4>() : (Ki3.G() * Jb3.L() + Ki3.L());
		
		phyJb->init(PhyJb::Identity, MatCal.Vacuum());
		phyJb->DeltaZ() = (phySt.Z() - St0.Z());
		phyJb->NumRadLen() = MatCal.NumRadLen();
		phyJb->G(JPX, JUX) += step * Ki0.G(JPX, JUX);
		phyJb->G(JPY, JUY) += step * Ki0.G(JPY, JUY);
		for (Int_t it = JUX; it < PhyJb::GDim; ++it) {
			phyJb->G(JPX, it) += oneSixthSquareStep * (Ki0.G(JUX, it) + KJb1G(JUX, it) + KJb2G(JUX, it));
			phyJb->G(JPY, it) += oneSixthSquareStep * (Ki0.G(JUY, it) + KJb1G(JUY, it) + KJb2G(JUY, it));
			phyJb->G(JUX, it) += oneSixthStep * (Ki0.G(JUX, it) + MgntProp::TWO * KJb1G(JUX, it) + MgntProp::TWO * KJb2G(JUX, it) + KJb3G(JUX, it));
			phyJb->G(JUY, it) += oneSixthStep * (Ki0.G(JUY, it) + MgntProp::TWO * KJb1G(JUY, it) + MgntProp::TWO * KJb2G(JUY, it) + KJb3G(JUY, it));
			phyJb->G(JGB, it) += oneSixthStep * (Ki0.G(JGB, it) + MgntProp::TWO * KJb1G(JGB, it) + MgntProp::TWO * KJb2G(JGB, it) + KJb3G(JGB, it));
		}	
		for (Int_t it = 0; it < PhyJb::LDim && !phyJb->Vacuum(); ++it) {
			phyJb->L(JPX, it) += step * Ki0.L(JPX, it) + oneSixthSquareStep * (Ki0.L(JUX, it) + KJb1L(JUX, it) + KJb2L(JUX, it));
			phyJb->L(JPY, it) += step * Ki0.L(JPY, it) + oneSixthSquareStep * (Ki0.L(JUY, it) + KJb1L(JUY, it) + KJb2L(JUY, it));
			phyJb->L(JUX, it) += oneSixthStep * (Ki0.L(JUX, it) + MgntProp::TWO * KJb1L(JUX, it) + MgntProp::TWO * KJb2L(JUX, it) + KJb3L(JUX, it));
			phyJb->L(JUY, it) += oneSixthStep * (Ki0.L(JUY, it) + MgntProp::TWO * KJb1L(JUY, it) + MgntProp::TWO * KJb2L(JUY, it) + KJb3L(JUY, it));
			phyJb->L(JGB, it) += oneSixthStep * (Ki0.L(JGB, it) + MgntProp::TWO * KJb1L(JGB, it) + MgntProp::TWO * KJb2L(JGB, it) + KJb3L(JGB, it));
		}
	}

	return vaild;
}

//---- OrthCoord ----//
// Orthogonal Coordinate based on seed vector (1, 0, 0)
// Taget is particle direction
OrthCoord::OrthCoord(PhySt & phySt, Bool_t isWithDev) {
	Double_t parallel = std::fabs(MtxLB::Dot(phySt.Dir(), OrthCoord::ORTH_SEED));
	if (MgntNum::Compare(parallel, ONE) == 0) return;
	else { init(); fTaget = phySt.Dir(); }
	Double_t idot = ONE / std::sqrt(ONE - fTaget(X) * fTaget(X));
	
	fAxisV(X) = ZERO;
	fAxisV(Y) = idot * (fTaget(Z)); 
	fAxisV(Z) = idot * (NEG * fTaget(Y));

	fAxisW(X) = idot * (ONE - fTaget(X) * fTaget(X));
	fAxisW(Y) = idot * (NEG * fTaget(X) * fTaget(Y)); 
	fAxisW(Z) = idot * (NEG * fTaget(X) * fTaget(Z));
	
	
	if (isWithDev) {
		Double_t idot1 = NEG * idot;
		Double_t idot3 = NEG * idot * idot * idot;
		
		fAxisVDev(X, X) = ZERO;
		fAxisVDev(X, Y) = ZERO;
		
		fAxisVDev(Y, X) = idot3 * (fTaget(X) * fTaget(Y) * fTaget(Y) / fTaget(Z));
		fAxisVDev(Y, Y) = idot1 * (fTaget(Y) / fTaget(Z)); 

		fAxisWDev(X, X) = idot1 * (fTaget(X));
		fAxisWDev(X, Y) = ZERO;

		fAxisWDev(Y, X) = idot3 * (fTaget(Y));
		fAxisWDev(Y, Y) = idot1 * (fTaget(X));
	}
}


//---- DevStatusFunc ----//
// Status (x, y, z, 1/eta) --- eta := ChrgSign * GammaBeta
// 1st Derivative (ds) of Status Function
// 2st Derivative (ds) of Status Function
// devFunc(x, y, z, ux, uy, uz, 1/eta)
DevStatusFunc::DevStatusFunc(PhySt & phySt, MatPropParam * matPar, MatPhyCalParam * matCal) {
	init();
	if (matPar == nullptr || matCal == nullptr) fVacuum = true;
	else { fVacuum = (matPar->Vacuum() || matCal->Vacuum()); }

	MtxLB::SVecD<3> & magVec = MgntMag::GetMagFast(phySt.Pos());
	Double_t        propFT = MgntProp::PROP_FACT * std::fabs(phySt.ChrgMass());

	Double_t sign = MgntNum::Compare(phySt.Eta());
	Double_t ieta = std::fabs(phySt.InvEta());
	Double_t propWithIETA = propFT * (sign * ieta);
	MtxLB::SVecD<3> & dir = phySt.Dir();
	
	Double_t crsUB[3] = { ZERO };
	crsUB[X] = (dir(Y) * magVec(Z) - dir(Z) * magVec(Y));
	crsUB[Y] = (dir(Z) * magVec(X) - dir(X) * magVec(Z));
	crsUB[Z] = (dir(X) * magVec(Y) - dir(Y) * magVec(X));
	
	fDevFunc(0) = dir(X); 
	fDevFunc(1) = dir(Y);
	fDevFunc(2) = dir(Z);
	fDevFunc(3) = propWithIETA * crsUB[X];
	fDevFunc(4) = propWithIETA * crsUB[Y];
	fDevFunc(5) = propWithIETA * crsUB[Z];

	if (!fVacuum) {
		Bool_t isLowEngLimit = (phySt.Beta() < MgntProp::BetaLimit);
		Double_t ieta2 = (isLowEngLimit) ? (MgntProp::IEtaLimit *  MgntProp::IEtaLimit) : (ieta * ieta);
		Double_t eng   = std::sqrt(ONE / ieta2 + ONE);
		Double_t ibeta = std::sqrt(ieta2 + ONE);
		Double_t ietab = (isLowEngLimit) ? (MgntProp::IEtaLimit / MgntProp::BetaLimit) : (ieta * ibeta);
		
		if (matPar->OptMscat()) {
			OrthCoord orth(phySt, 0);
			Double_t mscatDV = (matPar->MscatDV() * matCal->fMscatD);
			Double_t mscatDW = (matPar->MscatDW() * matCal->fMscatD);

			Double_t mscatE_Psi = ietab; 
			Double_t mscatDVPsi = mscatDV * mscatE_Psi;
			Double_t mscatDWPsi = mscatDW * mscatE_Psi;
			
			fDevFunc(3) += (mscatDVPsi * orth.V(X) + mscatDWPsi * orth.W(X));
			fDevFunc(4) += (mscatDVPsi * orth.V(Y) + mscatDWPsi * orth.W(Y));
			fDevFunc(5) += (mscatDVPsi * orth.V(Z) + mscatDWPsi * orth.W(Z));
		}	
		
		if (matPar->OptEngls()) {
			Double_t englsIN = (matPar->EnglsIN() * matCal->fEnglsISGM + matCal->fEnglsIMPV);	
			Double_t englsBR = (matPar->EnglsBR() * matCal->fEnglsBMEN);
			if (MgntNum::Compare(englsIN) <= 0) englsIN = ZERO;
			
			Double_t englsI_Psi = (ieta2 * eng) * ((isLowEngLimit) ? (ONE/MgntProp::BetaLimit/MgntProp::BetaLimit + ONE) : (ibeta * ibeta + ONE)); 
			Double_t englsB_Psi = (ieta2 * eng) * (eng - ONE); 
			Double_t englsIBPsi = (englsIN * englsI_Psi + englsBR * englsB_Psi);
			
			//fDevFunc(3) += englsIBPsi * (-dir(X));
			//fDevFunc(4) += englsIBPsi * (-dir(Y));
			//fDevFunc(5) += englsIBPsi * (-dir(Z));
			fDevFunc(6) += englsIBPsi * (sign * ieta);
		}
	}
}


//---- DevParamFunc ----//
// Status (x, y, z, 1/eta) --- eta := ChrgSign * GammaBeta
// 1st Derivative of Status Function Relation with Parameters
// 2st Derivative of Status Function Relation with Parameters
// devMtxG(x, y, ux, uy, 1/eta)
// devMtxL(mscatL, mscatD, mscatC, englsI)
DevParamFunc::DevParamFunc(PhySt & phySt, MatPropParam * matPar, MatPhyCalParam * matCal) {
	init();
	if (matPar == nullptr || matCal == nullptr) fVacuum = true;
	else { fVacuum = (matPar->Vacuum() || matCal->Vacuum()); }

	MtxLB::SVecD<3> & magVec = MgntMag::GetMagFast(phySt.Pos());
	Double_t        propFT = MgntProp::PROP_FACT * std::fabs(phySt.ChrgMass());
	
	Double_t sign = MgntNum::Compare(phySt.InvEta());
	Double_t ieta = std::fabs(phySt.InvEta());
	Double_t propWithIETA = propFT * (sign * ieta);
	MtxLB::SVecD<3> & dir = phySt.Dir();
	Double_t tanx = phySt.TanX();
	Double_t tany = phySt.TanY();
	
	Double_t crsUB[2] = { ZERO };
	crsUB[X] = (dir(Y) * magVec(Z) - dir(Z) * magVec(Y));
	crsUB[Y] = (dir(Z) * magVec(X) - dir(X) * magVec(Z));
	
	Double_t crsUBX_DevU[2] = { ZERO };
	Double_t crsUBY_DevU[2] = { ZERO };
	crsUBX_DevU[X] = (tanx * magVec(Y));
	crsUBX_DevU[Y] = (+magVec(Z) + tany * magVec(Y));
	crsUBY_DevU[X] = (-tanx * magVec(X) - magVec(Z));
	crsUBY_DevU[Y] = (-tany * magVec(X));

	fDevMtxG(0, 2) = ONE;
	fDevMtxG(1, 3) = ONE;

	fDevMtxG(2, 2) = propWithIETA * crsUBX_DevU[X];
	fDevMtxG(2, 3) = propWithIETA * crsUBX_DevU[Y];
	fDevMtxG(2, 4) = propFT * crsUB[X];

	fDevMtxG(3, 2) = propWithIETA * crsUBY_DevU[X];
	fDevMtxG(3, 3) = propWithIETA * crsUBY_DevU[Y];
	fDevMtxG(3, 4) = propFT * crsUB[Y];

	if (!fVacuum) {
		Bool_t isLowEngLimit = (phySt.Beta() < MgntProp::BetaLimit);
		Double_t ieta2 = (isLowEngLimit) ? (MgntProp::IEtaLimit * MgntProp::IEtaLimit) : (ieta * ieta);
		Double_t eng   = std::sqrt(ONE / ieta2 + ONE);
		Double_t ibeta = std::sqrt(ieta2 + ONE);
		Double_t ietab = (isLowEngLimit) ? (MgntProp::IEtaLimit / MgntProp::BetaLimit) : (ieta * ibeta);

		if (matPar->OptMscat()) {
			OrthCoord orth(phySt, 1);

			Double_t mscatE_Psi  = ietab;
			Double_t mscatD_DevD = (matCal->fMscatD * mscatE_Psi); 

			fDevMtxL(2, 0) += mscatD_DevD * orth.V(X);
			fDevMtxL(2, 1) += mscatD_DevD * orth.W(X);

			fDevMtxL(3, 0) += mscatD_DevD * orth.V(Y);
			fDevMtxL(3, 1) += mscatD_DevD * orth.W(Y);
		}
		
		if (matPar->OptEngls()) {
			Double_t englsI_Psi = (ieta2 * eng) * (ibeta * ibeta + ONE); 
			Double_t englsB_Psi  = (ieta2 * eng) * (eng - ONE); 

			Double_t englsI_DevI = (matCal->fEnglsISGM);
			Double_t englsB_DevB = (matCal->fEnglsBMEN);

			//fDevMtxL(2, 2) += englsI_DevI * englsI_Psi * (-dir(X));
			//fDevMtxL(2, 3) += englsB_DevB * englsB_Psi * (-dir(X));

			//fDevMtxL(3, 2) += englsI_DevI * englsI_Psi * (-dir(Y));
			//fDevMtxL(3, 3) += englsB_DevB * englsB_Psi * (-dir(Y));

			fDevMtxL(4, 2) += englsI_DevI * englsI_Psi * (sign * ieta);
			fDevMtxL(4, 3) += englsB_DevB * englsB_Psi * (sign * ieta);
		}
	}
}


#endif // __MgntProp_C__
