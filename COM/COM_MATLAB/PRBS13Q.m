function [ seq syms syms_nrz ] = PRBS13Q( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


taps = ([13 12 2  1]);
seed =([0 0 0 0 0 1 0 1 0 1 0 1 1]);
[seq_nrz c] =LFSR(seed,taps);
seq_nrz=2*(seq_nrz-0.5);
seq=pam(seq_nrz);
% syms=round(2*(seq+1));
syms((round(2*(seq+1))/2==2))=3;
syms((round(2*(seq+1))/2==1.5))=2;
syms((round(2*(seq+1))/2==.5))=1;
syms((round(2*(seq+1))/2==0))=0;

% syms_nrz=((seq_nrz+1)/2);

syms_nrz=seq_nrz;
