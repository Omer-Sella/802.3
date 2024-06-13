%Make RX Package
%Important:  RX pkg doesn't change based on mode being dd or dc.  So 'dd' is always passed
[ s11in, s12in, s21in, s22in]=make_full_pkg('RX',faxis,param,channel_type,'dd',include_die);


%     p(1 ,1, :)=s11in;
%     p(2 ,2, :)=s22in;
%     p(1 ,2, :)=s12in;
%     p(2 ,1, :)=s21in;
%
%     S=sparameters(p,faxis);
%     rfwrite(S,'temp.s4p');

if strcmpi(mode,'dc')
    RTX=R_diepad(param.Tx_rd_sel)/2;
    RRX=R_diepad(param.Rx_rd_sel)/2;
    Z0gamma=Z0/2;
else
    RTX=R_diepad(param.Tx_rd_sel);
    RRX=R_diepad(param.Rx_rd_sel);
    Z0gamma=Z0;
end
if OP.IDEAL_TX_TERM || (OP.RX_CALIBRATION == 1 && channel_number == 2) || OP.include_pcb == 2
    gamma_tx=0;
else
    gamma_tx=(RTX-Z0gamma)/(RTX+Z0gamma);% equation 93A-17
end
if OP.IDEAL_RX_TERM
    gamma_rx=0;
else
    gamma_rx=(RRX-Z0gamma)/(RRX+Z0gamma);% equation 93A-17
end

if OP.INC_PACKAGE==0
    s21p= s21;
    warning('do not use INC_PACKAGE = 0. Instead use package parameters)');
else
    if OP.RX_CALIBRATION == 1 && channel_number == 2
        %   for calibration do not include the transmitter package
        [s11out_rx, s12out_rx, s21out_rx, s22out_rx ] = combines4p( s11, s12, s21, s22, s22out, s12out, s21out, s11out ); %#ok<ASGLU> % s22 is ball side of package
        SCH.Frequencies=faxis;
        SCH.Parameters(1,1,:)=s11out_rx;
        SCH.Parameters(2,2,:)=s22out_rx;
        SCH.Parameters(1,2,:)=s12out_rx;
        SCH.Parameters(2,1,:)=s21out_rx;
        SCH.NumPorts=2;
        SCH.Impedance=100;
        %% Equation 93A-18
        if include_die
            s21p= s21out_rx.*(1-gamma_tx).*(1+gamma_rx)./(1.- s11out_rx.*gamma_tx - s22out_rx.*gamma_rx  -s21out_rx.^2.*gamma_tx.*gamma_rx +s11out_rx.*s22out_rx.*gamma_tx.*gamma_rx);
        else
            s21p=s21out_rx; % if no die we do not want a VTF
        end
    else
        %% Equations 93A-4 to 93A-7
        if ~OP.IDEAL_TX_TERM
            [s11, s12, s21, s22] = combines4p( s11out, s12out, s21out, s22out, s11, s12, s21, s22 ); %#ok<ASGLU>
        end
        H_t=ones(1,length(faxis)); % .3bj compatibility
        if OP.IDEAL_TX_TERM ||  OP.T_r_filter_type == 1
            % for RITT testing with good termination as in some instruments
            % and tx filter when required
            if OP.T_r_filter_type==0
                H_t = exp(-(pi*faxis/1e9*OP.transmitter_transition_time/1.6832).^2); %% Equation 93A-46 %%
            else
                tr=OP.transmitter_transition_time;
                f9=faxis/1e9;
                if OP.T_r_meas_point == 1
                    k=1.9466+7.12*sqrt(1-6.51e-3/tr);
                    H_t=105./(f9.^4*(k*tr)^4 - f9.^3*(k*tr)^3*10i - 45*f9.^2*(k*tr)^2 + f9*(k*tr)*105i + 105);
                else
                    H_t = exp( -2*(pi*f9*tr/1.6832).^2 ).*exp(-1j*2*pi*f9*tr*3);
                end
                
            end
        end
        if strcmpi(mode,'dc')
%             H_t=ones(1,length(faxis)); % not sure if we need a H_t or not. Is the CM noise correlated to the edge rate?
        end
        if ~OP.IDEAL_RX_TERM
            [s11, s12, s21, s22] = combines4p( s11, s12, s21, s22, s22in, s21in, s12in, s11in ); %#ok<ASGLU> % s22 is ball side of package
        else
            warning('do not use IDEAL_RX_TERM. Instead hard code package and TR parameters')
        end
        %% Equation 93A-18 and part of 93A-1: Ht fix in V290 identified by Bill Kirkland and Ed Frlan ( s21^2 changed to s12*s21 )
        if include_die
            s21p= H_t.*s21.*(1-gamma_tx).*(1+gamma_rx)./(1.- s11.*gamma_tx - s22.*gamma_rx  -s21.*s12.*gamma_tx.*gamma_rx +s11.*s22.*gamma_tx.*gamma_rx);
        else
            s21p=s21; % if no die we do not want a VTF
        end
    end
    
    if strcmpi(mode,'dc')
        % compute AC_CM_RMS at tp0
        OP.TX_BesselThomson=1; %AC CM is measured in scopes with at BT filter
        H_bt=Bessel_Thomson_Filter(param,faxis,OP.TX_BesselThomson);
        if channel_number == 1
            f_int= faxis( faxis<=param.ACCM_MAX_Freq );
            H_cc= s21out.*(1-gamma_tx)./(1.- s11out.*gamma_tx ).*H_bt;
            sigma_ACCM_at_tp0= sqrt(2*param.AC_CM_RMS_TX^2*sum( abs( H_cc(2:length(f_int)) ).^2 .* diff(f_int))/f_int(end)) ;
            %     S=sparameters(p,faxis);
            %     rfwrite(S,'temp.s4p');
        end
    end
    
    SCH.Frequencies=faxis;
    SCH.Parameters(1,1,:)=s11;
    SCH.Parameters(2,2,:)=s22;
    SCH.Parameters(1,2,:)=s12;
    SCH.Parameters(2,1,:)=s21;
    SCH.NumPorts=2;
    if strcmpi(mode,'dc')
        SCH.Impedance=25;
    else
        SCH.Impedance=100;
    end
    
end