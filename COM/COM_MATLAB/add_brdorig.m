function [ s11out, s12out, s21out, s22out ] = add_brdorig(chdata, param, OP)
% Used in Clause 92 for adding board trace between TP0 and TP2
switch chdata.type
    case 'THRU'
        z_bp_tx = param.z_bp_tx;
    case 'NEXT'
        z_bp_tx = param.z_bp_next;
    case 'FEXT'
        z_bp_tx = param.z_bp_fext;
end

[ s11tx, s12tx, s21tx, s22tx ] = synth_tline(chdata.faxis, param.brd_Z_c(1), param.Z0, param.brd_gamma0_a1_a2, param.brd_tau, z_bp_tx);
[ s11rx, s12rx, s21rx, s22rx ] = synth_tline(chdata.faxis, param.brd_Z_c(2), param.Z0, param.brd_gamma0_a1_a2, param.brd_tau, param.z_bp_rx);

switch OP.include_pcb
    case 1
        [ s11out1, s12out1, s21out1, s22out1 ]=combines4p(  s11tx, s12tx, s21tx, s22tx, chdata.sdd11_raw, chdata.sdd12_raw, chdata.sdd21_raw, chdata.sdd22_raw );
        [ s11out, s12out, s21out, s22out ]=combines4p(  s11out1, s12out1, s21out1, s22out1,  s11rx, s12rx, s21rx, s22rx);
    case 2
        [ s11out, s12out, s21out, s22out ]=combines4p( chdata.sdd11_raw, chdata.sdd12_raw, chdata.sdd21_raw, chdata.sdd22_raw, s11rx, s12rx, s21rx, s22rx);
end