function Beta=hbftms_BetaQ(QP,QN,QW,dpos,dmom)
% HBFTMS_BETAQ computes fast beta integrals for a TMS coil
%
% Beta=HBFTMS_BETAQ_(QP,QN,QW,dpos,dmom)
%   QP = coil quadrature points, [N x 3]
%   QN = coil normals, [N x 3]
%   QW = quadrature weights, [1 x N]
%   dpos  = Beta dipole positions, [M x 3]
%   dmom =  Beta dipole moments, [M x 3]
%
%   Beta  = resulting magnetic flux for each dipole, [M x 1]
%
%   Computes the beta integral of Eq. 8 in Stenroos--Koponen, 2019.
%   The needed 'dpos' and 'dmom' are computed with
%   HBFTMS_BETAQDIPOLES.
%
% v191119 (c) Matti Stenroos, matti.stenroos@aalto.fi

%original: hbf_Binf_dir_TMS_v4
Nqp=size(QP,1);
Nsp=size(dpos,1);
Beta=zeros(Nsp,1);

sd1=dmom(:,1);
sd2=dmom(:,2);
sd3=dmom(:,3);

for F=1:Nqp
    R=bsxfun(@minus,QP(F,:),dpos);
    absRsq=R(:,1).^2+R(:,2).^2+R(:,3).^2;
    
    Beta=Beta+QW(F)*...
    (R(:,1).*(QN(F,2).*sd3-QN(F,3).*sd2)+...
    R(:,2).*(QN(F,3).*sd1-QN(F,1).*sd3)+...
    R(:,3).*(QN(F,1).*sd2-QN(F,2).*sd1))./(absRsq.*sqrt(absRsq));   
end
