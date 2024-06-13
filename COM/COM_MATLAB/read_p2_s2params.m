function [data, SDD, SDC] = read_p2_s2params(infile, plot_ini_s_params, plot_dif_s_params, ports,OP)
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
%   Modified December 2021 to use read_Nport_touchstone
%   This is faster reader that is capable of reading touchstone with any number of ports
%
% Input Variables (required)
%   infile              -- The s4p file to be read and converted
%   plot_ini_sparams    -- Plot the initial s-parameter information. For debugging purposes
%   plot_dif_s_params   -- Plot the differential s-parameter information. For debugging purposes
%   ports               -- Re-order the port layout
%
% Output/Return Variables
%   data        -- structure containing network parameter data points and frequency axis
%   sdc         -- the differential in/common-mode out s-parameter data matrix
%   sdd         -- the differential in/differential out s-parameter data matrix
%


% backwards compatibility settings. can be removed in updated code.
if ~exist('OP', 'var'); OP.DISPLAY_WINDOW = true; end
if isempty(ports); ports = [1 2]; end % default order normally used.
ports = [1 2];


if OP.DISPLAY_WINDOW
    set(0,'defaulttextinterpreter','none') % prevents subscripting character in displayed messages
    hMsgBox = msgbox(infile, 'Reading S-Parameter File'); % display a progress bar for reading the s-parameter file(s)
end

%AJG:  fast touchstone read for any number of ports
[sch,schFreqAxis]=read_Nport_touchstone(infile,ports);



D=NaN(size(sch));
% calculate differential s parameter matrix from single ended
for i=1:size(sch,1)
    S(:,:) = sch(i,:,:);
    T = [1 1  ; 1 -1  ];
    W = T * (S / T);
    D(i,:,:) = W(:,:);
end

% D matrix should be
% Scc11 Scd11 Scc12 Scd21
% Sdc11 Sdd11 Sdc12 Sdd12
% Scc21 Scd21 Scc22 Scd22
% Sdc21 Sdd21 Sdc22 Sdd22

% proper values
%AJG:  matrix can be properly referenced after fixing mapping
SDD(:,1,1) = D(:,2,2);
SDC(:,1,1)=  D(:,2,1);
SCC(:,1,1)=  D(:,1,1);
SCD(:,1,1)=  D(:,1,2);



% backwards compatibility output variables
data.m = sch;
%schFreqAxis=schFreqAxis(1:freqCounter); % truncating preallocated array to number of freq points.
data.freq = schFreqAxis;
colors = 'rgbk';

if (plot_ini_s_params == 1)
    figure('name', 'Single-ended s-parameters');set(gcf,'Tag','COM');
    for mj=1:4
        %         subplot(2,2,mj);
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
    for mj=1:1
        for mi=1:1
            plot(data.freq, 20*log10(abs(squeeze(SDD(:,mj,mi)))), ...
                colors((mj-1)*2+mi), 'linewidth',2, 'disp', sprintf('SDD%d%d', mj, mi));
            hold on
        end
    end
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    legend show
    grid on
    title(infile);
    
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
% end read_sp2_sparam