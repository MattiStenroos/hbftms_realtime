function Bfinf=hbftms_BpFlux_xyz(QP,QN,QW,spos)
% HBFTMS_BPFLUX_XYZ computes magnetic flux thru a coil
%   (due to a unit dipole triplets in infinite non-magnetic medium) 
%
% Binf=HBFTMS_BPFLUX_XYZ(QP,QN,QW,spos)
%   QP = coil quadrature points, [N x 3]
%   QN = coil normals, [N x 3]
%   QW = quadrature weights, [1 x N]
%   spos  = source positions, [M x 3]
%
%   Bfinf = resulting magnetic flux for each dipole component-wise,
%       [Bfinfx Bfinfy Bfinfz], [M x 3]
% v191119 (c) Matti Stenroos, matti.stenroos@aalto.fi

Nfp=size(QP,1);
Nsp=size(spos,1);

Bfinf=zeros(Nsp,3);
KQW=1e-7*QW; %mu0/(4pi)*QW
zerov=zeros(Nsp,1);
for F=1:Nfp
    R=bsxfun(@minus,QP(F,:),spos);
    absRsq=sum(R.*R,2);
    absRm3 = 1./(absRsq.*sqrt(absRsq));
    sourcex=[zerov -R(:,3) R(:,2)]*QN(F,:)';
    sourcey=[R(:,3) zerov -R(:,1)]*QN(F,:)';
    sourcez=[-R(:,2) R(:,1) zerov]*QN(F,:)';
    Btemp=bsxfun(@times,[sourcex sourcey sourcez],absRm3);
    
    Bfinf=Bfinf+KQW(F)*Btemp;
end
Bfinf=reshape(Bfinf',1,Nsp*3)';

