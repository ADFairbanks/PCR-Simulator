//Species indecies as an object
const species_index = {P1: 0, P2: 1, S1: 2, S2: 3, S1S2: 4, E: 5, S1P2: 6, S2P1: 7, S1P2E: 8, S2P1E: 9, dNTP: 10, S1Q2E: 11, S2Q1E: 12, S1Q2: 13, S2Q1: 14, Q1: 15, Q2: 16, M2: 17, M1: 18, S1M2: 19, S2M1: 20, S1M2E: 21, S2M1E: 22, S1N2E: 23, S2N1E: 24, S1N2: 25, S2N1: 26, S1L2: 27, S2L1: 28, L2: 29, L1: 30, L2P1: 31, L1P2: 32, L2P1E: 33, L1P2E: 34, L2Q1E: 35, L1Q2E: 36, L2Q1: 37, L1Q2: 38, D_E : 39, SPECIES_COUNT: 40}; 
//Labels used in model in concentration:
//P1 and P2: Primers
//S1 and S2: DNA/Plasmid single strands
//S1S2: DNA/Plasmid double strand
//E: Polymerase
//S1P2 and S2P1: DNA + Primer bound
//S1P2E and S2P1E: DNA + Primer bound + Polymerase
//dNTP: Free Nucleotides
//S1Q2E and S2Q1E: DNA + Extended Primer bound + Polymerase
//S1Q2 and S2Q1: DNA + Extended Primer bound
//Q1 and Q2: Extended Primer
//M2 to L1Q2. Unused species. Intended for Szabo's mismatch model.
//Denatured polymerase (D_E) was added kind of in the last second because I thought it was needed for purity calculations but it wasn't.

//PCR Constants
const kh = 2.5e5;
const kPolym = 1.9 * 1e7; // 1/Ms
const ke_curveH = 1e5;
const kd_A = 0.0008; //Polymerase Denaturation kinetic constant
const kd_Ea = 800; // cal/mol //Polymerase Denaturation kinetic constant
const P1_stock_conc = 10e-6;// M
const P2_stock_conc = 10e-6;// M
const dNTP_stock_conc = 10e-3;// M
const total_reaction_volume = 50; // uL -- total reaction volume
const plasmid_length = 5000; //bp
const amplicon_length = 1000; //bp
const primer_length = 20; //nt
const deltaS = -300; // cal/mol
const deltaG_polymerase = -11000; // cal/mol
const R = 1.9872036; // cal/molK | Gas Constant
const CT_melt_curves = 500e-9; 
const polymerase_nucleotide_extension_step = 100;
const extension_length = amplicon_length - polymerase_nucleotide_extension_step - primer_length;
const Amplicon_Tm = NEB_Tm(amplicon_length);
const deltaH_Amplicon = calc_deltaH(Amplicon_Tm, deltaS, R, CT_melt_curves);
const Primer1_Tm = Primer2_Tm = NEB_Tm(primer_length);
const Long_Primer1_Tm = LongPrimer2_Tm = NEB_Tm(primer_length + polymerase_nucleotide_extension_step);
const deltaH_Primer1 = deltaH_Primer2 = calc_deltaH(Primer1_Tm, deltaS, R, CT_melt_curves);
const deltaH_Ext_Primer1 = deltaH_Ext_Primer2 = calc_deltaH(Long_Primer1_Tm, deltaS, R, CT_melt_curves);

//The hold temperature and time after all the cycles are completed. Currently set as constants.
const HOLD_TEMP = 5; //Celcius
const HOLD_TIME = 600; //Seconds

function calc_polymerase_denaturation_rate(Tc){
    return kd_A * Math.exp((-1*kd_Ea) / (R * (Tc + 273.15)));
}

function lognormal(Tc,M,S){
	//The log-normal curve. The horizontally flipped version of this 
	//is used to work out polymerase extension kinetic forward rate, given temperature
	//Tc:		Celsius
    var epower = (-1 * Math.pow((Math.log(Tc) - M),2))/(2*S*S);
    var coeff = 1/(S*Math.sqrt(6.28)*Tc)
    return coeff * Math.exp(epower);
}

const ke_curveV_TAQ = 1;
const ke_curveM_TAQ = 3.53;
const ke_curveS_TAQ = 0.27;
const ke_curveV_Phusion = 2.6;
const ke_curveM_Phusion = 3.4;
const ke_curveS_Phusion = 0.27;

function calc_polymerase_extension_rate(temperature_C, polymerase_type="Taq"){
    var ke_curveV, ke_curveM, ke_curveS;
    if (polymerase_type == "Taq"){
        ke_curveV = ke_curveV_TAQ;
        ke_curveM = ke_curveM_TAQ;
        ke_curveS = ke_curveS_TAQ;
    }else if (polymerase_type == "Phusion"){
        ke_curveV = ke_curveV_Phusion;
        ke_curveM = ke_curveM_Phusion;
        ke_curveS = ke_curveS_Phusion;
    }
    return ke_curveH * ke_curveV * lognormal((100.0001-temperature_C),ke_curveM,ke_curveS);
}
function Keq_vant_hoff(Tc, deltaH, deltaS, R){
	T = Tc + 273.15;
	frac1 = -1 * (deltaH / (R * T));
	frac2 = deltaS / R;

	return Math.exp(frac1 + frac2);
}
function calc_hybridisation_reverse_rate_constant(kf, Tc, deltaH, deltaS, R){
    return kf / Keq_vant_hoff(Tc, deltaH, deltaS, R);
}
function NEB_Tm(basepairs){
    //Taq and Phusion have the same params currently
    param = [0.00221, -0.003396, 83.1336];
    return -1 * (1/((param[0] * basepairs) + param[1])) + param[2];
}
function calc_deltaH(Tm, deltaS, R, CT){
	T = Tm + 273.15;
	return T * (deltaS - (R * Math.log(4.0/CT)));
}
function calc_polymerase_reverse_rate_constant(kf, Tc, deltaG, R){
	T = Tc + 273.15;
	Keq = Math.exp((-1 * deltaG) / (R * T));
	return kf / Keq;
}
function ng_mass_to_nanomolar(ng, bp, V_ul){
    mol = (ng / 1e9) / (660 * bp);
    M = mol / (V_ul / 1e6);
    return M * 1e9;
}
function polymerase_units_to_nanomolar(U){
    return U * (200/1.25);
}