//This file is the entry point of the simulator and the callback function that exits the simulator
//Main functions:
//  pcr_main() - Entry point from UI
//  yield_purity_callback() - Exit point
//  calculate_yield()
//  calculate_purity()
//      dna_concentration_to_mass() - Used to calculate yield and purity
//  start_pcr() - Calculates needed constants from params
//  build_PCR_schedule()
//  call_pcr_ODE() - Calls the ODE recursively

function pcr_main(){
    const params = {};
    params.cycles = document.getElementById("cycles_input").value*1.0;
    params.time_denat = document.getElementById("denat_time_input").value*1.0;
    params.temp_denat = document.getElementById("denat_temp_input").value*1.0;
    params.time_anneal = document.getElementById("anneal_time_input").value*1.0;
    params.temp_anneal = document.getElementById("anneal_temp_input").value*1.0;
    params.time_extend = document.getElementById("extend_time_input").value*1.0;
    params.temp_extend = document.getElementById("extend_temp_input").value*1.0;
    params.mass_dna = document.getElementById("mass_dna_input").value*1.0;
    params.vol_primers = document.getElementById("vol_primers_input").value*1.0;
    params.vol_dntp = document.getElementById("vol_dntp_input").value*1.0;
    params.polymerase_units = document.getElementById("polymerase_units_input").value*1.0;
    params.polymerase_type = document.getElementById("polymerase_type_input").value;
    const callback_function = yield_purity_callback;
    start_pcr(params, callback_function);
}


function yield_purity_callback(init_state, final_state){
    const yield_mass = calculate_yield(init_state, final_state);
    const purity = calculate_purity(final_state);
    document.getElementById("outDiv").innerHTML = "<p>Yeild: "+yield_mass.toFixed(2)+" ng/uL</p><p>Purity: "+purity.toFixed(2)+"%</p><p>Init DNA: "+init_state[species_index.S1S2]+"<br>Final DNA: "+final_state[species_index.S1S2]+"</p>";
    const a = [yield_mass,purity];
    return a;
}

function dna_concentration_to_mass(dna_conc,base_pairs,nucleotides=0){
    //2x to convert from double to single strand
    const dna_nucleotides = base_pairs*2 + nucleotides;
    //308.97 g/mol is the average molecular weight of each nucleotide
    //18.02 g/mol is the -OH and -H added to the ends
    //1e9 is the converstion from nanograms to grams
    //1e6 is the conversion from L to uL
    const dna_mass = ((dna_nucleotides * 1e9 * 308.97)+18.02)/1e6;
    const conc_mass = dna_conc * dna_mass;
    return conc_mass;
}

function calculate_yield(init_state, final_state){
    //Yield in ng/ul
    //The unit of the state is mol/liter
    const amplicon_index = species_index.S1S2;
    const init_amplicon_conc = init_state[amplicon_index];
    const final_amplicon_conc = final_state[amplicon_index];
    const diff_conc = final_amplicon_conc - init_amplicon_conc;
    var yield_mass = 0;
    if (diff_conc > 0){
        yield_mass = dna_concentration_to_mass(diff_conc,amplicon_length);
    }
    return yield_mass;
}

function calculate_purity(final_state){
    //Purity as percent of DNA amplicon / total mass 
    const dna_mass = dna_concentration_to_mass(final_state[species_index.S1S2],amplicon_length);
    var all_mass = 0
    for (i=0; i<final_state.length; i++){
        var nM_conc = final_state[i];
        var nt = 0
        var bp = 0
        //Single Strand
        if (i == species_index.S1 || i == species_index.S2){
            nt = amplicon_length;
        //Double Strand (DNA Amplicon)
        }else if (i == species_index.S1S2){
            bp = amplicon_length;
        //Primer
        }else if (i == species_index.P1 || i == species_index.P2){
            nt = primer_length;
        //Elongated Primer shorter than amplicon
        }else if (i == species_index.Q1 || i == species_index.Q2){
            nt = primer_length + polymerase_nucleotide_extension_step;
        //Single Strand with Primer
        }else if (i == species_index.S1P2 || i == species_index.S2P1 || i == species_index.S1P2E || i == species_index.S2P1E){
            nt = amplicon_length - primer_length;
            bp = primer_length;
        //Single Strand with Elongated Primer
        }else if (i == species_index.S1Q2 || i == species_index.S2Q1 || i == species_index.S1Q2E || i == species_index.S2Q1E){
            nt = amplicon_length - primer_length - polymerase_nucleotide_extension_step;
            bp = primer_length + polymerase_nucleotide_extension_step;
        }
        //if (nt!=0 || bp!=0){
        if (i!=species_index.dNTP){
            all_mass += dna_concentration_to_mass(nM_conc,bp,nt);
        }
    }
    const purity = 100*dna_mass / all_mass;
    return purity;
}

/*
//Display the final result of each species. Useful for testing.
const state_labels = [
    "P1", "P2", "S1", "S2", "S1S2", "E", "S1P2", "S2P1", "S1P2E", "S2P1E",
    "dNTP", "S1Q2E", "S2Q1E", "S1Q2", "S2Q1", "Q1", "Q2", "M2", "M1", "S1M2",
    "S2M1", "S1M2E", "S2M1E", "S1N2E", "S2N1E", "S1N2", "S2N1", "S1L2", "S2L1", "L2",
    "L1", "L2P1", "L1P2", "L2P1E", "L1P2E", "L2Q1E", "L1Q2E", "L2Q1", "L1Q2"
];
function state_pcr_callback(init_state, final_state){
    let s='';
    for (var i=0; i<state_labels.length; i++){
        s+=state_labels[i] + ": " + final_state[i] + "<br>";
    }
    document.getElementById("outDiv").innerHTML = s;
}
*/

function start_pcr(user_params, callback_function){
    //Initialization
    const P1_vol = user_params.vol_primers;
    const P2_vol = user_params.vol_primers;
    const dNTP_vol = 4*user_params.vol_dntp;
    const E_Units = user_params.polymerase_units;
    const Plasmid_mass_ng = user_params.mass_dna;
    
    var init_state = new Array(species_index.SPECIES_COUNT).fill(0);
    init_state[species_index.P1] = P1_vol * P1_stock_conc / total_reaction_volume; //Primer 1 concentration
    init_state[species_index.P2] = P2_vol * P2_stock_conc / total_reaction_volume; //Primer 2 concentration
    init_state[species_index.S1S2] = ng_mass_to_nanomolar (Plasmid_mass_ng, plasmid_length, total_reaction_volume) / 1e9; //DNA
    init_state[species_index.E] = polymerase_units_to_nanomolar (E_Units) / 1e9; //Polymerase
    init_state[species_index.dNTP] = dNTP_vol * dNTP_stock_conc / (total_reaction_volume); //Nucleotides
    const PCR_Schedule = build_PCR_schedule(user_params);
    //console.log("starting pcr at time 0 up to " + PCR_Schedule[0].time + " at temperature "+PCR_Schedule[0].temperature);
    //console.log("initial state: "+init_state);
    const final_state = call_pcr_ODE(init_state, PCR_Schedule, user_params.polymerase_type);
    return callback_function(init_state, final_state);
}

function build_PCR_schedule(user_parameters){
    const schedule = [];
    for (var i=0; i<user_parameters.cycles; i++){
        const denat_schedule_obj = {temperature:user_parameters.temp_denat, time:user_parameters.time_denat};
        schedule.push(denat_schedule_obj);
        const anneal_schedule_obj = {temperature:user_parameters.temp_anneal, time:user_parameters.time_anneal};
        schedule.push(anneal_schedule_obj);
        const extend_schedule_obj = {temperature: user_parameters.temp_extend, time:user_parameters.time_extend};
        schedule.push(extend_schedule_obj);
    }
    //Hold is added to hybridize both strands after PCR
    const hold_schedule_obj = {temperature:HOLD_TEMP, time:HOLD_TIME};
    schedule.push(hold_schedule_obj);
    return schedule;
}

//Initial call:
//  State: Initial y0 vector
//  Schedule: List of tuple/object: (temperature and time)
//  Time_Passed: Default of t0 = 0 can be left as is. Used in the recursive calls
function call_pcr_ODE(state, schedule, polymerase_type="Taq", time_passed=0){
    if (schedule.length == 0){
        return state;
    }
    //Determine temperature based constants for the ODE
    const temperature = schedule[0].temperature;
    const kh_amplicon_hybridize = calc_hybridisation_reverse_rate_constant(kh, temperature, deltaH_Amplicon, deltaS, R);
    const kh_polymerase_reverse = calc_polymerase_reverse_rate_constant(kPolym, temperature, deltaG_polymerase, R); 
    const kh_primer1_bind_reverse = kh_primer2_bind_reverse = calc_hybridisation_reverse_rate_constant(kh, temperature, deltaH_Primer1, deltaS, R);
    const kh_ext_primer1_bind_reverse = kh_ext_primer2_bind_reverse = calc_hybridisation_reverse_rate_constant(kh, temperature, deltaH_Ext_Primer1, deltaS, R);
    const polymerase_denaturation_rate = calc_polymerase_denaturation_rate(temperature);
    
    //Need to push the rates into a class for this ODE solver
    const kinetic_rates = new Rates_Class();
    //non-temperature based
    kinetic_rates.push_rate("kh",kh);
    kinetic_rates.push_rate("kPolym",kPolym);
    //temperature based
    kinetic_rates.push_rate("kh_primer1_bind_reverse",kh_primer1_bind_reverse);
    kinetic_rates.push_rate("kh_primer2_bind_reverse",kh_primer2_bind_reverse);
    kinetic_rates.push_rate("polymerase_reverse",kh_polymerase_reverse);
    kinetic_rates.push_rate("polymerase_extension_rate",calc_polymerase_extension_rate(temperature,polymerase_type));
    kinetic_rates.push_rate("kh_ext_primer1_bind_reverse",kh_ext_primer1_bind_reverse);
    kinetic_rates.push_rate("kh_ext_primer2_bind_reverse",kh_ext_primer2_bind_reverse);
    kinetic_rates.push_rate("polymerase_denaturation_rate",polymerase_denaturation_rate);
    kinetic_rates.push_rate("kh_amplicon_hybridize",kh_amplicon_hybridize);
    
    //Debug parameter checking. Convenient to see if you want to modify or verify the ODE model.
    /*
    console.log("Temperature: "+temperature);
    console.log("kh: "+kh);
    console.log("kPolym: "+kPolym);
    console.log("kh_primer_bind_reverse: "+kh_primer1_bind_reverse);
    console.log("kh_polymerase_reverse: "+kh_polymerase_reverse);
    console.log("polymerase_extension_rate: "+calc_polymerase_extension_rate(temperature,polymerase_type));
    console.log("kh_ext_primer_bind_reverse: "+kh_ext_primer1_bind_reverse);
    console.log("polymerase_denaturation_rate: " +polymerase_denaturation_rate);
    console.log("kh_amplicon_hybridize: "+kh_amplicon_hybridize);
    console.log("extension_length: "+extension_length);
    */
    let current_state = state; //This would have to be a bit different if you want all state steps returned
    //Set parameters
    current_state = Mass_Action_Kinetics_ODE_Solve(current_state, schedule[0].time, PCR_ODE, kinetic_rates, dt=0.05);
    //console.log("Resulting state at time: ",(schedule[0].time + time_passed)," - ",current_state);
    return call_pcr_ODE(current_state, schedule.slice(1), polymerase_type, schedule[0].time + time_passed);
}
