function [data, SDD, SDC, SCC ] = read_p4_s4params(infile, plot_ini_s_params, plot_dif_s_params, ports,OP,param)
%% FUNCTION :: read_sp4_sparams
%
% Description
%   Read the fid of single-ended 4-port complex S-parameters
%   in Touchstone format 'file' and convert to the internal
%   format using the port transform 'ports'
%
%   Created by Mike Y. He
%   April 22, 2005
%
%   Reused some code from
%   Anthony Sanders, Alex Deas, Bob Davidov (24 January 2005)
%   for touchstone 4-port S-matrix import.
%
%   Modified (2012-July-27) by Ken Young to match current indexing scheme and
%   optimized for quicker parameter matching and parsing. also, separated out
%   the plotting algorithms into their own sub-function routines
%
%  Modified December 2021 to use read_Nport_touchstone
%  This is faster reader that is capable of reading touchstone with any number of ports
%
% Input Variables (required)
%   infile              -- The s4p file to be read and converted
%   plot_ini_sparams    -- Plot the initial s-parameter information. For debugging purposes
%   plot_dif_s_params   -- Plot the differential s-parameter information. For debugging purposes
%   ports               -- Re-order the port layout
%   OP
%   param
% Output/Return Variables
%   data        -- structure containing network parameter data points and frequency axis
%   sdd         -- the differential in/differential out s-parameter data matrix
%   sdc         -- the differential in/common-mode out s-parameter data matrix
%   scc         -- the common mode in/common-mode out s-parameter data matrix
%
%


% backwards compatibility settings. can be removed in updated code.
if ~exist('OP', 'var'); OP.DISPLAY_WINDOW = true; end
if isempty(ports); ports = [1 3 2 4]; end % default order normally used.

% adjust ports to maintain the meaning [in1, in2 , out1, out2] when one
% pair is reversed.
%ports_adj=ports; for k=1:4, ports_adj(k)=find(ports==k); end; ports=ports_adj;

if OP.DISPLAY_WINDOW
    hMsgBox = msgbox(infile, 'Reading S-Parameter File'); % display a progress bar for reading the s-parameter file(s)
end

%AJG:  fast touchstone read for any number of ports
[sch,schFreqAxis]=read_Nport_touchstone(infile,ports);
% matrix to introduce p or n skew on Tx or Rx RIM 12/29/2023
% Sigma's will be form exp(2i*pi*f*skew*1e-12). i.e. if skew = 0 sigma = 1
% need to swap sigma for 1 and 3 and 2 and 4 not RIM 12/29/2023
Sigfct = ...
    @(sigma2,sigma1,sigma4,sigma3)reshape([sigma1.^2,sigma1.*sigma2,sigma1.*sigma3,sigma1.*sigma4,sigma1.*sigma2,sigma2.^2,sigma2.*sigma3,sigma2.*sigma4,sigma1.*sigma3,sigma2.*sigma3,sigma3.^2,sigma3.*sigma4,sigma1.*sigma4,sigma2.*sigma4,sigma3.*sigma4,sigma4.^2],[4,4]);
D=NaN(size(sch));
% calculate differential s parameter matrix from single ended
% skew added RIM 12/29/2023
for i=1:size(sch,1)
    f=schFreqAxis(i);
    sigma_matrix=Sigfct(exp(2i*pi*f*param.Txpskew*1e-12),exp(2i*pi*f*param.Txnskew*1e-12),exp(2i*pi*f*param.Rxpskew*1e-12),exp(2i*pi*f*param.Rxnskew*1e-12) );
    S(:,:) = sch(i,:,:);
    Snew=sigma_matrix.*S;
    T = [1 1 0 0 ; 1 -1 0 0 ; 0 0 1 1 ; 0 0 1 -1];
    W = T * (Snew / T);
    D(i,:,:) = W(:,:);
end

% D matrix should be
% Scc11 Scd11 Scc12 Scd21
% Sdc11 Sdd11 Sdc12 Sdd12
% Scc21 Scd21 Scc22 Scd22
% Sdc21 Sdd21 Sdc22 Sdd22

% proper values
SDD(:,1,1) = D(:,2,2);
SDD(:,2,2) = D(:,4,4);
SDD(:,1,2) = D(:,2,4);
SDD(:,2,1) = D(:,4,2);

SDC(:,1,1) = D(:,2,1);
SDC(:,2,2) = D(:,4,3);
SDC(:,1,2) = D(:,2,3);
SDC(:,2,1) = D(:,4,1);

SCC(:,1,1) = D(:,1,1);
SCC(:,2,2) = D(:,3,3);
SCC(:,1,2) = D(:,1,3);
SCC(:,2,1) = D(:,3,1);

% backwards compatibility output variables
data.m = sch;
%schFreqAxis=schFreqAxis(1:freqCounter); % truncating preallocated array to number of freq points.
data.freq = schFreqAxis;
colors = 'rgbk';

if (plot_ini_s_params == 1)
    figure('name', 'Single-ended s-parameters');set(gcf,'Tag','COM');
    for mj=1:4
        subplot(2,2,mj);
        for mi=1:4
            plot(data.freq, 20*log10(abs(data.m(:,mj,mi)+1.0e-15)), ...
                colors(mi), 'linewidth', 2, 'disp', sprintf('S%d%d', mj, mi));
            hold on
        end
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        legend show
        grid on
        title(sprintf('Output port %d', mj));
    end
end
plot_dif_s_params =0;
if (plot_dif_s_params == 1)
    figure('name', 'Mixed-mode s-parameters');set(gcf,'Tag','COM');
    %     subplot(2,1,1);
    for mj=1:2
        for mi=1:2
            plot(data.freq, 20*log10(abs(SDD(:,mj,mi))), ...
                colors((mj-1)*2+mi), 'linewidth',2, 'disp', sprintf('SDD%d%d', mj, mi));
            hold on
        end
    end
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    legend show
    grid on
    title(infile);
    %
    %     subplot(2,1,2);
    %     for mj=1:2
    %         for mi=1:2
    %             plot(data.freq, 20*log10(abs(SDC(:,mj,mi))+1.0e-15), ...
    %                 colors((mj-1)*2+mi), 'linewidth',2, 'disp', sprintf('SDC%d%d', mj, mi));
    %             hold on
    %         end
    %     end
    %     xlabel('Frequency (Hz)');
    %     ylabel('Magnitude (dB)');
    %     legend show
    %     grid on
end

if OP.DISPLAY_WINDOW, close(hMsgBox); end