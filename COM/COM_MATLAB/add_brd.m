function [ s11out, s12out, s21out, s22out ] = add_brd(chdata, param, OP)
%% Used in Clause 92 for adding board trace between TP0 and TP2

switch chdata.type
    case 'THRU'
        z_bp_tx = param.z_bp_tx;
        z_bp_rx = param.z_bp_rx;
    case 'NEXT'
        z_bp_tx = param.z_bp_rx;
        z_bp_rx = param.z_bp_next;
    case 'FEXT'
        z_bp_tx = param.z_bp_fext;
        z_bp_rx = param.z_bp_rx;
end
% Same cap on each tx and rx three is a data stratue for bifrucation but
% logic no implemented here RIM 06/28/2019
zref=param.Z0;
c1=param.C_0;
c2=param.C_1;
f=chdata.faxis;
f(f<eps)=eps;
s11pad1t= -1i*2*pi.*f*c1(1)*zref./(2+1i*2*pi.*f*c1(1)*zref);
s21pad1t= 2./(2+1i*2*pi.*f*c1(1)*zref);
s11pad2t= -1i*2*pi.*f*c2(1)*zref./(2+1i*2*pi.*f*c2(1)*zref);
s21pad2t= 2./(2+1i*2*pi.*f*c2(1)*zref);
[ s11tx, s12tx, s21tx, s22tx ] = synth_tline(chdata.faxis, param.brd_Z_c(1), param.Z0, param.brd_gamma0_a1_a2, param.brd_tau, z_bp_tx);
% add Tx caps
[s11tx, s12tx, s21tx, s22tx  ]= ...
    combines4p(  s11pad1t, s21pad1t, s21pad1t, s11pad1t, s11tx, s12tx, s21tx, s22tx   );
[s11tx, s12tx, s21tx, s22tx  ]= ...
    combines4p(   s11tx, s12tx, s21tx, s22tx,s11pad2t, s21pad2t, s21pad2t, s11pad2t  );


s11pad1r= -1i*2*pi.*f*c1(2)*zref./(2+1i*2*pi.*f*c1(2)*zref);
s21pad1r= 2./(2+1i*2*pi.*f*c1(2)*zref);
s11pad2r= -1i*2*pi.*f*c2(2)*zref./(2+1i*2*pi.*f*c2(2)*zref);
s21pad2r= 2./(2+1i*2*pi.*f*c2(2)*zref);
[ s11rx, s12rx, s21rx, s22rx ] = synth_tline(chdata.faxis, param.brd_Z_c(2), param.Z0, param.brd_gamma0_a1_a2, param.brd_tau, z_bp_rx);
% add Rx caps
[s11rx, s12rx, s21rx, s22rx  ]= ...
    combines4p(  s11pad2r, s21pad2r, s21pad2r, s11pad2r, s11rx, s12rx, s21rx, s22rx   );
[s11rx, s12rx, s21rx, s22rx  ]= ...
    combines4p(   s11rx, s12rx, s21rx, s22rx,s11pad1r, s21pad1r, s21pad1r, s11pad1r  );


switch OP.include_pcb
    case 1
        [ s11out1, s12out1, s21out1, s22out1 ]=combines4p(  s11tx, s12tx, s21tx, s22tx, chdata.sdd11_raw, chdata.sdd12_raw, chdata.sdd21_raw, chdata.sdd22_raw );
        [ s11out, s12out, s21out, s22out ]=combines4p(  s11out1, s12out1, s21out1, s22out1,  s11rx, s12rx, s21rx, s22rx);
    case 2
        [ s11out, s12out, s21out, s22out ]=combines4p( chdata.sdd11_raw, chdata.sdd12_raw, chdata.sdd21_raw, chdata.sdd22_raw, s11rx, s12rx, s21rx, s22rx);
end
