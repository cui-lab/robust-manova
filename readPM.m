% extract PM values from CEL files. 
%To run it, it needs the bioinformatics toolbox installed. 
%Following the example and then assemble PM values from multiple CEL files into the 3-D format for rMANOVA package.


function [PM,PSname] = readPM(CelName,libDir);
% return PM (NumGene by NumProbePair) and ProbeSetName of probe sets 
%   with regular probe pair number (ususally 11)
%
% e.g., [PM,PSname]=readPM(CelName,libDir);
%       CelName,libDir are in character arrays
%  CelName='38347-38341.CEL';
%  libDir = 'P:\AG\Affytools\CDF Files';

CELStruct=affyread(char(CelName));
CDFStruct=affyread([CELStruct.ChipType '.CDF'],libDir);

% Probe pair number detection 
NumProbeSets=length(CDFStruct.ProbeSets);
mNumPairs(1)=CDFStruct.ProbeSets(round(NumProbeSets*0.15)).NumPairs;
mNumPairs(2)=CDFStruct.ProbeSets(round(NumProbeSets*0.35)).NumPairs;
mNumPairs(3)=CDFStruct.ProbeSets(round(NumProbeSets*0.55)).NumPairs;
mNumPairs(4)=CDFStruct.ProbeSets(round(NumProbeSets*0.75)).NumPairs;
mNumPairs(5)=CDFStruct.ProbeSets(round(NumProbeSets*0.95)).NumPairs;
NumPairs=median(mNumPairs);
%
numCols = CDFStruct.Cols;
%
RegPairID=find([CDFStruct.ProbeSets.NumPairs]==NumPairs);
NumRegPair=length(RegPairID);
PSname={CDFStruct.ProbeSets(RegPairID).Name}';
%
thePairs = [CDFStruct.ProbeSets(RegPairID).ProbePairs];
PMX = thePairs(:,[3:6:end]);
PMY = thePairs(:,[4:6:end]);
PMRow = reshape(PMY*numCols+PMX+1,NumPairs*NumRegPair,1);
PM = reshape([CELStruct.Probes(PMRow,3)],NumPairs,NumRegPair)';

