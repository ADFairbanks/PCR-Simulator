//The only function in this file is PCR_ODE

//Required kinetic constants for this ODE. 
//("*" denotation means it's effected by temperature)
//(1 and 2 pairs are currently the same numbers because the primers are assumed to be the same length and Tm)
//  kh
//  kh_amplicon_hybridize*
//  kh_primer1_bind_reverse*
//  kh_primer2_bind_reverse*
//  kPolym
//  kh_polymerase_reverse*
//  polymerase_extension_rate*
//  kh_ext_primer1_bind_reverse*
//  kh_ext_primer2_bind_reverse*
//  polymerase_denaturation_rate*
//Required pcr constants for this ODE
//  extension_length (It's the amplicon - long primer)
//  polymerase_nucleotide_extension_step (long primer - primer)

Array.prototype.vectorAdd = function(other) {
    return this.map((x,i) => x + other[i]);
}

function PCR_ODE(state, t, rates){
    // teach Array some vector arithmetics
    Array.prototype.vectorAdd = function(other) {
        return this.map((x,i) => x + other[i]);
    }
    var s = species_index; //Shorten the variable name since it's used a lot

    //Note: The forward hybridization rate constant is "kh" which is often used in these reactions and does not change by temperature. When "kh" is used, the other rate should change by temperature.

    var SPECIES_COUNT = s.SPECIES_COUNT; //another shortened variable name
    
    //Denaturation and Hybridization #1
    function hybridization_reaction(){
        // S1 S2 <-> S1S2
        const kh = rates.get_rate("kh");
        const kh_amplicon_hybridize = rates.get_rate("kh_amplicon_hybridize");
        const flow_hybridization = kh * state[s.S1] * state[s.S2];
        const flow_dehybridization = kh_amplicon_hybridize * state[s.S1S2];
        var ret = new Array(SPECIES_COUNT).fill(0);
        ret[s.S1S2] = flow_hybridization - flow_dehybridization;
        ret[s.S1] = -flow_hybridization + flow_dehybridization;
        ret[s.S2] = -flow_hybridization + flow_dehybridization;
        return ret;
    }
    
    //Primer Binding [1st primer] #2
    function primer_binding_reaction1(){
        // S1 P2 <-> S1P2
        //I have no idea why there needs to be a 0.5 here but the model is very inaccurate without it.
        //const flow_binding = 0.5*rates.get_rate("kh") * state[s.S1] * state[s.P2];
        //const flow_binding = rates.get_rate("binding_multiplier")*rates.get_rate("kh") * state[s.S1] * state[s.P2];
        const flow_binding = rates.get_rate("kh") * state[s.S1] * state[s.P2];
        const flow_unbinding = rates.get_rate("kh_primer1_bind_reverse") * state[s.S1P2];
        var ret = new Array(SPECIES_COUNT).fill(0);
        ret[s.S1P2] = flow_binding - flow_unbinding;
        ret[s.S1] = -flow_binding + flow_unbinding;
        ret[s.P2] = -flow_binding + flow_unbinding;
        return ret;
    }
    
    //Primer Binding [2nd primer] #2
    function primer_binding_reaction2(){
        // S2 P1 <-> S2P1
        //const flow_binding = rates.get_rate("binding_multiplier")*rates.get_rate("kh") * state[s.S2] * state[s.P1];
        const flow_binding = rates.get_rate("kh") * state[s.S2] * state[s.P1];
        const flow_unbinding = rates.get_rate("kh_primer2_bind_reverse") * state[s.S2P1];
        var ret = new Array(SPECIES_COUNT).fill(0);
        ret[s.S2P1] = flow_binding - flow_unbinding;
        ret[s.S2] = -flow_binding + flow_unbinding;
        ret[s.P1] = -flow_binding + flow_unbinding;
        return ret;
    }

    //Polymerase Binding [1st primer] #3
    function polymerase_binding_reaction1(){
        //S1P2 <-> S1P2E
        const flow_binding = rates.get_rate("kPolym") * state[s.S1P2] * state[s.E];
        const flow_unbinding = rates.get_rate("polymerase_reverse") * state[s.S1P2E];
        var ret = new Array(SPECIES_COUNT).fill(0);
        ret[s.S1P2] = -flow_binding + flow_unbinding;
        ret[s.E] = -flow_binding + flow_unbinding;
        ret[s.S1P2E] = flow_binding - flow_unbinding;
        return ret;
    }
    
    //Polymerase Binding [2nd primer] #3
    function polymerase_binding_reaction2(){
        //S2P1 E <-> S2P1E
        const flow_binding = rates.get_rate("kPolym") * state[s.S2P1] * state[s.E];
        const flow_unbinding = rates.get_rate("polymerase_reverse") * state[s.S2P1E];
        var ret = new Array(SPECIES_COUNT).fill(0);
        ret[s.S2P1] = -flow_binding + flow_unbinding;
        ret[s.E] = -flow_binding + flow_unbinding;
        ret[s.S2P1E] = flow_binding - flow_unbinding;
        return ret;
    }
    
    //Primer Extension [1st Primer] #4
    function primer_short_extension_reaction1(){
        //S1P2E -> S1Q2E
        const Amplicon_Ratio = 1 - (polymerase_nucleotide_extension_step / amplicon_length);
        const Extension_Rate = state[s.S1P2E] * (state[s.dNTP]/4) * rates.get_rate("polymerase_extension_rate") * Amplicon_Ratio;
        var ret = new Array(SPECIES_COUNT).fill(0);
        ret[s.S1P2E] = -Extension_Rate;
        ret[s.S1Q2E] = Extension_Rate;
        ret[s.dNTP] = -Extension_Rate * polymerase_nucleotide_extension_step;
        return ret;
    }
    
    //Primer Extension [2nd Primer] #4
    function primer_short_extension_reaction2(){
        //S2P1E -> S2Q1E
        const Amplicon_Ratio = 1 - (polymerase_nucleotide_extension_step / amplicon_length);
        const Extension_Rate = state[s.S2P1E] * (state[s.dNTP]/4) * rates.get_rate("polymerase_extension_rate") * Amplicon_Ratio;
        var ret = new Array(SPECIES_COUNT).fill(0);
        ret[s.S2P1E] = -Extension_Rate;
        ret[s.S2Q1E] = Extension_Rate;
        ret[s.dNTP] = -Extension_Rate * polymerase_nucleotide_extension_step;
        return ret;
    }
    
    //Long Primer Polymerase Binding [1st Primer] #5 (Note: Works the same as non-extended primer)
    function polymerase_binding_reaction_long1(){
        //S1Q2E <-> S1Q2 E
        const flow_binding = rates.get_rate("kPolym") * state[s.S1Q2] * state[s.E];
        const flow_unbinding = rates.get_rate("polymerase_reverse") * state[s.S1Q2E];//*10
        var ret = new Array(SPECIES_COUNT).fill(0);
        ret[s.S1Q2] = -flow_binding + flow_unbinding;
        ret[s.E] = -flow_binding + flow_unbinding;
        ret[s.S1Q2E] = flow_binding - flow_unbinding;
        return ret;
    }
    
    //Long Primer Polymerase Binding [2nd Primer] #5 (Note: Works the same as non-extended primer)
    function polymerase_binding_reaction_long2(){
        //S2Q1E <-> S2Q1 E
        const flow_binding = rates.get_rate("kPolym") * state[s.S2Q1] * state[s.E];
        const flow_unbinding = rates.get_rate("polymerase_reverse") * state[s.S2Q1E];//*10
        var ret = new Array(SPECIES_COUNT).fill(0);
        ret[s.S2Q1] = -flow_binding + flow_unbinding;
        ret[s.E] = -flow_binding + flow_unbinding;
        ret[s.S2Q1E] = flow_binding - flow_unbinding;
        return ret;
    }
    
    //Long Primer Binding [1st Primer] #6
    function primer_binding_reaction_long1(){
        //S1Q2 <-> S1 Q2
        const flow_binding = rates.get_rate("kh") * state[s.S1] * state[s.Q2];
        const flow_unbinding = rates.get_rate("kh_ext_primer1_bind_reverse")  * state[s.S1Q2];
        var ret = new Array(SPECIES_COUNT).fill(0);
        ret[s.S1Q2] = flow_binding - flow_unbinding;
        ret[s.S1] = -flow_binding + flow_unbinding;
        ret[s.Q2] = -flow_binding + flow_unbinding;
        return ret;
    }
    
    //Long Primer Binding [2nd Primer] #6
    function primer_binding_reaction_long2(){
        //S2Q1 <-> S2 Q1
        const flow_binding = rates.get_rate("kh") * state[s.S2] * state[s.Q1];
        const flow_unbinding = rates.get_rate("kh_ext_primer2_bind_reverse") * state[s.S2Q1];
        var ret = new Array(SPECIES_COUNT).fill(0);
        ret[s.S2Q1] = flow_binding - flow_unbinding;
        ret[s.S2] = -flow_binding + flow_unbinding;
        ret[s.Q1] = -flow_binding + flow_unbinding;
        return ret;
    }
    
    //Annealing [1st Primer] #7
    function full_extension_reaction1(){
        //S1Q2E dNTP <-> S1S2 E
        const Amplicon_Ratio = (polymerase_nucleotide_extension_step / amplicon_length);
        const Extension_Rate = state[s.S1Q2E] * (state[s.dNTP]/4) * rates.get_rate("polymerase_extension_rate") * Amplicon_Ratio; //(amplicon_length/polymerase_nucleotide_extension_step);
        var ret = new Array(SPECIES_COUNT).fill(0);
        ret[s.dNTP] = -Extension_Rate * extension_length;
        ret[s.S1Q2E] = -Extension_Rate;
        ret[s.S1S2] = Extension_Rate;
        ret[s.E] = Extension_Rate;
        return ret;
    }

    //Annealing [2nd Primer] #7
    function full_extension_reaction2(){
        //S2Q1E dNTP <-> S1S2 E
        const Amplicon_Ratio = (polymerase_nucleotide_extension_step / amplicon_length);
        const Extension_Rate = state[s.S2Q1E] * (state[s.dNTP]/4) * rates.get_rate("polymerase_extension_rate") * Amplicon_Ratio;//(amplicon_length/polymerase_nucleotide_extension_step);
        var ret = new Array(SPECIES_COUNT).fill(0);
        ret[s.dNTP] = -Extension_Rate * extension_length;
        ret[s.S2Q1E] = -Extension_Rate;
        ret[s.S1S2] = Extension_Rate;
        ret[s.E] = Extension_Rate;
        return ret;
    }
    
    //Polymerase Denaturation #8
    function polymerase_denaturation(){
        var ret = new Array(SPECIES_COUNT).fill(0);
        denatured_polymerase = rates.get_rate("polymerase_denaturation_rate") * state[s.E];
        ret[s.E] = -denatured_polymerase;
        ret[s.D_E] = denatured_polymerase;
        return ret;
    }
    
    //Incomplete model notes using mismatched primers as defined by Szabo.
    //S1 + P2 <-> S1M2
    //S2 + P1 <-> S2M1
    //S1M2 + E <-> S1M2E
    //S2M1 + E <-> S2M1E
    //S1M2E + dNTP -> S1N2E
    //S2M1E + dNTP -> S2N1E
    //S1N2E <-> S1N2 + E
    //S2N1E <-> S2N1 + E
    //S1N2E + dNTP -> S1L2 + E
    //S2N1E + dNTP -> S2L1 + E
    //S1L2 <-> S1 + L2
    //S2L1 <-> S2 + L1
    //L1 + P2 <-> L1P2
    //L2 + P1 <-> L2P1
    //L1P2 + E <-> L1P2E
    //L2P1 + E <-> L2P1E
    //L1P2E + dNTP -> L1Q2E
    //L2P1E + dNTP -> L2Q1E
    //L1Q2E + dNTP -> S2L1
    //L2Q1E + dNTP -> S1L2
    
    var ODE_LIST = [
        hybridization_reaction,// (#1)
        primer_binding_reaction1,// (#2)
        primer_binding_reaction2,
        polymerase_binding_reaction1,// (#3)
        polymerase_binding_reaction2,
        primer_short_extension_reaction1,// (#4)
        primer_short_extension_reaction2,
        polymerase_binding_reaction_long1,// (#5)
        polymerase_binding_reaction_long2,
        primer_binding_reaction_long1,// (#6)
        primer_binding_reaction_long2,
        full_extension_reaction1,// (#7)
        full_extension_reaction2,
        polymerase_denaturation //(#8)
    ];
    var state_vector = new Array(SPECIES_COUNT).fill(0);
    for (var i=0; i<ODE_LIST.length; i++){
        var vector = ODE_LIST[i]();
        for (var j=0; j<vector.length; j++){
            //Sanity checks with the  ODE solver and input state.
            //If you want the error messages to use state labels instead of index numbers, use state_labels[j] instead of j
            if (isNaN(vector[j])){
                console.log("NaN ERROR of equation #"+i+" on value "+j);
                console.log("Current state: "+state);
                throw "Equation "+i+" in the ODE set on state index "+j+" is not a number";
            }
            if (!isFinite(vector[j])){
                console.log("Infinite ERROR of equation #"+i+" on value "+j);
                console.log("Current state: "+state);
                throw "Equation "+i+" in the ODE set on state index "+j+" approaches infinity";
            }
            
            //Add reaction to state
            state_vector[j] += vector[j];

            //Another sanity check
            if (isNaN(state_vector[j])){
                console.log("NaN ERROR in accumulated state vector value "+j+" of equation "+i);
                console.log("Current accumulated vector state: "+state_vector);
                throw "Equation "+i+" in the ODE set caused the state index "+j+" number to not be a number";
                //Maybe that only throws when the initial state of that index is not a number, which is not the ODE's fault but the fault of the initial state given.
            }      
        }
    }
    return state_vector;
}
