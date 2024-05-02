
# TODO: compare my version to the below for consistency
# From DL's matlab code
%% CREATES EXTERNAL ADJ FILE FOR ONLINE PPD
% INPUTS:
%   DICOM OF 1 GRE VOLUME (EXTRACT AFFINE MATRIX FROM HEADERS)
%   NIFTI OF ALL GRE VOLUMES OR RAW DATA OF ALL GRE VOLUMES
%   NIFTI OF AFI VOLUME(S)
% OUTPUTS:
%   EXTERNAL ADJ STRUCTURE
%
% PIPELINE:
% 1.CALCULATE AFFINE MATRIX FROM DICOM
% 2.MAP {GRE, AFI, BET} TO DICOM SPACE
% 3.CALCULATE TX MAPS AND COORDINATES (IN DICOM SPACE)
% 4.SAVE EXTERNAL ADJ STRUCTURE
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%    
    
%%% CALCULATE AFFINE TRANSFORMATION MATRIX FROM DICOM HEADERS
dcm_files = dir(DCMdir); dcm_files(1:2) = [];

dcmhdr1 = dicominfo([dcm_files(1).folder,'/',dcm_files(1).name]);
dcmhdr2 = dicominfo([dcm_files(end).folder,'/',dcm_files(end).name]);

nRow = double(dcmhdr1.Rows); 
nCol = double(dcmhdr1.Columns); 
nSli = numel(dcm_files);   

%CALCULATE SCALING COMPONENT
RowColSpacing  = dcmhdr1.PixelSpacing;
dRow = double(RowColSpacing(1));       dX = [1; 1; 1].*dRow;
dCol = double(RowColSpacing(2));       dY = [1; 1; 1].*dCol;
dSli = double(dcmhdr1.SliceThickness); dZ = [1; 1; 1].*dSli;
Adcm_scaling = [dX dY dZ];

%CALCULATE ORIENTATION COMPONENT
dircosRowCol = double(dcmhdr1.ImageOrientationPatient);
dircosRow = dircosRowCol(1:3); 
dircosCol = dircosRowCol(4:6);
T1 = double(dcmhdr1.ImagePositionPatient); 
TN = double(dcmhdr2.ImagePositionPatient);
dircosSli = ((TN-T1)./(nSli-1))./dZ;
Adcm_orientation = [dircosCol dircosRow dircosSli]; %why not Row Col Sli?

Adcm = [[Adcm_orientation.*Adcm_scaling, T1];[0 0 0 1]];

