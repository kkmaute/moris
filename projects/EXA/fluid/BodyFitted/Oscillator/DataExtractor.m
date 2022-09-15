function [Fedata, Eledata_tri6, Eledata_quad8 ] = Input_data()
%This function extracts Grid Data, Triangular Element Data and Quad element
%data from a *.nas file (NASTRAN format) or *.txt file. 
% Input File: *.nas or *.txt
% Output: Matrix of [Fedata, Eledata_tri6, Eledata_quad8]
%   Detailed explanation goes here

clc;
clear all;
format long;
fid = fopen('sharpGeom_studentMesh.nas','rt');
A = textscan(fid, '%s','HeaderLines',4);  % After running the code until here, check A to see each cell entry
Ax = [A{:}];

%% FEDATA Extraction
% Input file consists of 3 data set. idx1 find total entries of "GRID".
idx1 = find(contains(Ax,'GRID*'));
% n1 is a counter that allows the for loop to iterate over first set of
% data to extract Fedata.
n1 = length(idx1);

for i=1:n1
    A1(:,i) = string(A{1}{5*i-3});
    B(:,i) = string(A{1}{5*i-2});
%     if strlength(B(:,i)) == 36
%        B(:,i) = insertAfter(B(:,i),15,'0');
%     else
%         B(:,i) = B(:,i);
%     end
%     Bprime(:,i) = strtrim(regexprep(B(:,i),'.{16}',' $0 '));
    C(:,i) = sscanf(B(:,i),'%16e %f*CONST'); % If it shows error here due to size on each side, make sure to check "B". If some entries in *.nas file are not properly aligned it may cause problem.
    z(:,i) = str2double(string(A{1}{5*i}));
    D(i,:) = transpose(C(:,i));
 end
 
 
 %Zp = sscanf(char(z),'%f');
 Z = transpose(z);
 
 %C = sscanf(char(B),'%e');
 
 Fedata = [str2double(transpose(A1)),D, Z];

%% Tri6 Extraction
% Tri6 Input data has 8 columns. Element-id, Block-id, node-1,
% node-2,....,node-6
cp1 = 8; % 8 columns of data, look at the output of Ay to identify this number
idx2 = find(contains(Ax,'CTRIA6'));
p1 = idx2(1);p2=idx2(length(idx2));
Ay = Ax(p1:p2+cp1);

for j=1:length(idx2)
    ele_idx(j,1)= Ay(9*j-7,1);
    block_id(j,1) = Ay(9*j-6,1);
    N1(j,1) = Ay(9*j-5,1);
    N2(j,1)=Ay(9*j-4,1);
    N3(j,1)=Ay(9*j-3,1);
    N4(j,1)=Ay(9*j-2,1);
    N5(j,1)=Ay(9*j-1,1);
    N6(j,1)=Ay(9*j,1);
end

Eledata_tri6(:,:) = str2double([ele_idx, block_id,N1, N2,N3,N4,N5,N6]);



%% Quad8 Extraction
% Quad8 input file has 11 columns. Element Index, Block-id, Node-1, Node-2,
% Node-3,...,Node-6, "CONST",..,Node-8.
cp2 = 11 ; %Look at the output of Az to identify this number
idx3 = find(contains(Ax,'QUAD8'));
p3 = idx3(1); p4 = idx3(length(idx3));
Az = Ax(p3:p4+cp2);

for k=1:length(idx3)
    ele_id2(k,1) = str2double(Az(12*k-10,1));
    block_id2(k,1) = str2double(Az(12*k-9,1));
    N11(k,1) = str2double(Az(12*k-8,1));
    N21(k,1) = str2double(Az(12*k-7,1));
    N31(k,1) = str2double(Az(12*k-6,1));
    N41(k,1) = str2double(Az(12*k-5,1));
    N51(k,1) = str2double(Az(12*k-4,1));
    N61(k,1) = sscanf(char(Az(12*k-3,1)),'%e');
    N71(k,1) = str2double(Az(12*k-1,1));
    N81(k,1) = str2double(Az(12*k,1));
end

Eledata_quad8(:,:) = [ele_id2, block_id2, N11,N21,N31,N41,N51,N61,N71,N81];


end