function result = readdataSnPx(filename, nport)
%function [freq, cs] = readdataSnPx(filename, nport)
% [freq, cs] = readdataSnP(filename, nport, format, nheader)
%
% Read Touchstone file with frequencies in units of Hertz
%
% Input:
% ======
% filename: Name of the Touchstone/SnP file
% nport: Number of ports
% format: 'RI' for real/imag, 'MA' for mag/angle (check option line in the
%         Touchstone file)
% nheader: Number of header lines (comment lines plus option line in the
%          Touchstone file)
%
% Output:
% =======
% freq: Vector of frequencies [Hz]
% cs: 3-D array of complex-valued S parameters where cs(i,j,k) is S(i,j)
%     at frequency freq(k)
%
% Note: If frequency unit is not Hertz (but GHz, MHz etc.) simply scale
% frequencies appropriately after reading the data.
%
% Ref.: Touchstone(R) File Format Specification, Rev.1.1,
% EIA/IBIS Open Forum, 2002.
%
% Written by Henning Braunisch, September 2004.
% Updated by Steven Krooswyk, April 2006.


fid = fopen(filename, 'r');


% Skip header lines
str = ' ';
n = 0;
while ~strcmp(str(1),'#')
    str = fgetl(fid);
    if isempty(str)
        str=' ' ;
        if n > 1000
            display('error: could not find config line (#)')
            break
        end
    end
    n = n + 1;
end

% parse configuration line
A=sscanf(str,'%1s %2s %1s %2s %1s %2s',[1,inf]);
p = find(A=='S');           %position of 'S'
units = lower(A(2:p-1));    %units before 'S'
format = A(p+1:p+2);        %format after 'S'

% skip any more header lines
%while ~str

nk = 0; % frequency counter
while 1
    
    [temp, count] = fscanf(fid, '%f', 1);
    if count == 0
        temp2 = fscanf(fid, '%s', 1);
        if ~isempty(temp2), fgetl(fid); continue, end;
        break
    end
    nk = nk+1; freq(1,nk) = temp;  %#ok<AGROW>
    for ni = 1:nport
        for nj = 1:nport
            switch lower(format)
                case 'ma'
                    mag = fscanf(fid, '%f', 1); ang = fscanf(fid, '%f', 1);
                    cs(ni,nj,nk) = mag * exp(1i*ang*pi/180); %#ok<AGROW>
                case 'ri'
                    re = fscanf(fid, '%f', 1); im = fscanf(fid, '%f', 1);
                    cs(ni,nj,nk) = complex(re, im); %#ok<AGROW>
                case 'db'
                    db = fscanf(fid, '%f', 1); ang = fscanf(fid, '%f', 1);
                    M = 10^(db/20);
                    %re = M*cos(ang);
                    %im = M*sin(ang);
                    re =  M*cos(ang * pi / 180);
                    im =  M*sin(ang * pi / 180);
                    cs(ni,nj,nk) = complex(re, im); %#ok<AGROW>
                otherwise
                    error('readdataSnP: Unknown data format');
            end
        end
    end
end

fclose(fid);

% If 2-port then swap S_12 and S_21 per Touchstone spec
if nport == 2
    temp = cs(2,1,:);
    cs(2,1,:) = cs(1,2,:);
    cs(1,2,:) = temp;
end

% Update freq units to Hz
switch lower(units)
    case 'hz'
        
    case 'khz'
        freq=freq.*1e3;
    case 'mhz'
        freq=freq.*1e6;
    case 'ghz'
        freq=freq.*1e9;
end

% passivity check
result.freq = freq;
result.cs    = cs;