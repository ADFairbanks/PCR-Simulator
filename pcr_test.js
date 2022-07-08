//Not needed at all for final sim

function big_test(){
    out_str = "Cycle,TD_denat,Ti_denat,TD_anneal,Ti_anneal,TD_extend,Ti_extend,dNTP,primer,plasmid,polymerase,yield,purity,polymerase_type\r\n";
    i=0;
    j=0;
    for (cycle_i = 10; cycle_i <= 15; cycle_i+=5){
        for (td_denat_i = 80; td_denat_i < 100; td_denat_i += 5){
            for (ti_denat_i = 10; ti_denat_i <= 30; ti_denat_i += 10){
                for (td_anneal_i = 40; td_anneal_i <= 70; td_anneal_i +=10){
                    for (ti_anneal_i = 30; ti_anneal_i <= 90; ti_anneal_i += 10){
                        for (td_extend_i = 40; td_extend_i <= 70; td_extend_i += 10){
                            for (ti_extend_i = 20; ti_extend_i <= 90; ti_extend_i += 10){
                                for (dntp_i = 0.5; dntp_i <= 3; dntp_i += 0.5){
                                    for (primer_i = 0.5; primer_i <= 3; primer_i += 0.5){
                                        for (plasmid_i = 10; plasmid_i <= 50; plasmid_i += 10){
                                            for (polymerase_i = 1; polymerase_i <= 3; polymerase_i++){
                                                j+=1;
                                                if (j%57187==0){
                                                    const params = {};
                                                    params.cycles = cycle_i;
                                                    params.time_denat = ti_denat_i;
                                                    params.temp_denat = td_denat_i;
                                                    params.time_anneal = ti_anneal_i;
                                                    params.temp_anneal = td_anneal_i;
                                                    params.time_extend = ti_extend_i;
                                                    params.temp_extend = td_extend_i;
                                                    params.mass_dna = plasmid_i;
                                                    params.vol_primers = primer_i;
                                                    params.vol_dntp = dntp_i;
                                                    params.polymerase_units = polymerase_i;
                                                    if (j%2==0){
                                                        params.polymerase_type = "Taq";
                                                    }else{
                                                        params.polymerase_type = "Phusion";
                                                    }
                                                    const callback_function = yield_purity_callback;
                                                    ret = start_pcr(params,callback_function);
                                                    out_str += cycle_i+","+td_denat_i+","+ti_denat_i+","+td_anneal_i+","+ti_anneal_i+","+td_extend_i+","+ti_extend_i+","+dntp_i+","+primer_i+","+plasmid_i+","+polymerase_i+","+ret[0].toFixed(2)+","+ret[1].toFixed(2)+","+params.polymerase_type+"\r\n";
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    }
    download_text(out_str,"big_test.csv");
    //download_text(out_str,"big_test.csv");
}
function download_text(text, filename) {
    var file = new Blob([text], {type: "text/plain;charset=utf-8"});
    var a = document.createElement("a"),
            url = URL.createObjectURL(file);
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    setTimeout(function() {
        document.body.removeChild(a);
        window.URL.revokeObjectURL(url);  
    }, 0); 

}