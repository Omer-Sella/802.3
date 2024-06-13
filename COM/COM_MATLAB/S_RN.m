function S_RN_of_f = S_RN(f,G_DC,G_DC2,param)
p1=param.CTLE_fp1(1);
z1=param.CTLE_fz(1);
p2=param.CTLE_fp2(1);
zlf=param.f_HP(1);
plf=param.f_HP(1);
f_b=param.fb;
f_r=param.f_r;
eta_0=param.eta_0;
H_CTF =   ( 10^(G_DC/20)+ 1i*f/z1 ) .*( 10^(G_DC2/20) + 1i*f/zlf )./ ( ( 1+1i*f/p1) .* ( 1+1i*f/p2)  .* ( 1+1i*f/plf));
H_R =1./polyval([1 2.613126 3.414214 2.613126 1], 1i*f./(f_r*f_b));
S_RN_of_f = eta_0/2.*abs( H_CTF.*H_R).^2; %EQ healey_3dj_01_2401 slide 5
if 0
    figure
    set(gcf, 'tag', 'COM');movegui(gcf,'southeast');
    % see if it looks correct
    semilogx(f/1e9, 20*log10( abs( H_CTF.*H_R)) );
    ylabel('dB');
    xlabel('GHz');
    title( 'H_ctf with H_r')
    grid on
    ylim([-30 0])
end