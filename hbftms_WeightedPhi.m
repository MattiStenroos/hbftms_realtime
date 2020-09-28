function Phiw=hbftms_WeightedPhi(bmeshes,Phi,ci,co)
% HBFTMS_WEIGHTEDPHI weights the Phi matrix
%
% Phiw=HBFTMS_WEIGHTEDPHI(bmeshes,Phi,ci,co)
% implements \phi in Eq. 5 of Stenroos--Koponen, NeuroImage 2019
%
% bmeshes:  boundary meshes, {Nmeshes x 1} cell array of hbg mesh structs
% Phi:      Potential matrix
% ci:       conductivities inside each boundary mesh, [Nmeshes x 1]
% co:       conductivities outside each boundary mesh, [Nmeshes x 1]
%
% Phiw:     weighted potentials (single precision, transposed)
% v200109 (c) Matti Stenroos, matti.stenroos@aalto.fi
Nm=numel(bmeshes);                  
[startind,endind]=NodeIndices(bmeshes);   %indexing of meshes in Phi matrix
mu0over4pi=1e-7;

Phiw=single(zeros(size(Phi)));
for M=1:Nm
    Phiw(startind(M):endind(M),:)=mu0over4pi*(co(M)-ci(M))*Phi(startind(M):endind(M),:);
  
    %Remove mean to correct for the numerical inaccuracy due to approximations in hbftms_BetaQ
    mp=mean(Phiw(startind(M):endind(M),:),1);
    Phiw(startind(M):endind(M),:)=bsxfun(@minus,Phiw(startind(M):endind(M),:),mp);
end
Phiw=Phiw.';

function [startinds,endinds]=NodeIndices(meshes)
Nsurf=length(meshes);
startinds=zeros(Nsurf,1);
endinds=zeros(Nsurf,1);
Nop=0;
for I=1:Nsurf
    startinds(I)=Nop+1;
    endinds(I)=startinds(I)+size(meshes{I}.p,1)-1;
    Nop=endinds(I);
end