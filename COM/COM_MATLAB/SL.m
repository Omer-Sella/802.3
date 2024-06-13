function SLD=SL(S,f,R)
% source load impact return loss add to S21
% S and SLD are the same structure
% S.Parameters
% S.Impedance
% S.NumPorts
% S.Frequencies
SLD=S; % assign the fields
zref=100;
if R==0
    warndlg('Termination should not be set to zero');
    SLD=S;
    return
end

if R > zref
    spr =R_series2(zref,f,(R-zref)); % make series sparameter
    %     SLD=sparameters(cascadesparams(S.Parameters,spr.Parameters),f,S.Impedance); % Casdeade
    [ SLD.Parameters(1,1,:), SLD.Parameters(1,2,:), SLD.Parameters(2,1,:), SLD.Parameters(2,2,:)]  = ...
        combines4p( ...
        spr.Parameters(1,1,:), spr.Parameters(1,2,:), spr.Parameters(2,1,:), spr.Parameters(2,2,:),...
        S.Parameters(1,1,:),   S.Parameters(1,2,:),   S.Parameters(2,1,:),   S.Parameters(2,2,:) ...
        );
elseif R < zref
    spr =r_parrelell2(zref,f,-R*zref/(R-zref)); % make parrellel sparameter
    %     SLD=sparameters(cascadesparams(S.Parameters,spr.Parameters),f,S.Impedance);
    [ SLD.Parameters(1,1,:), SLD.Parameters(1,2,:), SLD.Parameters(2,1,:), SLD.Parameters(2,2,:)]  = ...
        combines4p( ...
        spr.Parameters(1,1,:),  spr.Parameters(1,2,:), spr.Parameters(2,1,:), spr.Parameters(2,2,:),...
        S.Parameters(1,1,:),   S.Parameters(1,2,:),   S.Parameters(2,1,:),   S.Parameters(2,2,:) ...
        );
else
    SLD=S;
end

%%