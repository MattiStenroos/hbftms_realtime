function  coil  = make_coil_reference()
%MAKE_COIL_REFERENCE makes the coil model used as reference in (Stenroos & Koponen 2019)
%
%   coil  = MAKE_COIL_REFERENCE()
%
%   The function can, by changing the hard-coded coil parameters, be also
%   used for making other dipole-grid-based models of circular
%   figure-of-eight coils.
%   
%   DISCLAIMER: the function has been tested and verified only with the default parameters.
%
%   v210819 (c) Matti Stenroos (matti.stenroos@aalto.fi)

% changes:
% v210819: bugfix ('Nz' and 'el' are the same thing, changed to 'Nz')

%specs of the coil model (in mm)
Ri=26; %inner radius
Ro=44; %outer radius
N=9;   %number of windings
ww=1;  %wire width
wh=7;  %wire height
wo=3;  %coil casing thickness = wire offset from the bottom of the coil

%discretization parameters (lengths in mm)
ew=1;  %element side length
Nz=7;  %number of layers in z direction
sd=10; %sub-discretization of each element
%DISCLAIMER: the function has been tested & validated with only these
%parameters.

dr=(Ro-ww-Ri)/(N-1); %wire spacing
ri=(Ri:dr:Ro-1)';    %inner radius of each wire loop
ro=(Ri+ww:dr:Ro)';   %outer radius of each wire loop
k=1/ww; %slope of the moment density function within a loop

%sub-discretization of each element
sdsteps=linspace(-ew/2+ew/(2*sd),ew/2-ew/(2*sd),10);
[X0,Y0]=meshgrid(sdsteps,sdsteps);
Nsd=numel(X0);

%centers of all elements --- one layer , centered at [0 0 0]
dsteps=ew/2:ew:Ro;
dsteps=[-fliplr(dsteps),dsteps];
[PCX,PCY]=meshgrid(dsteps,dsteps);
Npc=numel(PCX);
pcset=[reshape(PCX,Npc,1),reshape(PCY,Npc,1)];

%resulting densities
Qset=zeros(Npc,1);

%loop over elements
for P=1:Npc
    pc=pcset(P,:);
    
    %element corners and the max & min radius from the coil centre
    ex=pc(1)+[-1 1]*ew/2;
    ey=pc(2)+[-1 1]*ew/2;
    pe=[ex(1) ey(1);ex(2) ey(1);ex(2) ey(2);ex(1) ey(2)];
    Re=sqrt(sum(pe.*pe,2));
    Rmax=max(Re);Rmin=min(Re);
    
    %coil rings that are fully outside or inside the element
    fullyout=ri>=Rmax;
    fullyin=ro<=Rmin;
    Nencircling=nnz(fullyout);
    
    %these rings cross the element
    targetrings=find(fullyout==fullyin);
    
    Ntarget=numel(targetrings);
    if Ntarget==0
        Qset(P)=Nencircling;
    else
        X=X0+pc(1);
        Y=Y0+pc(2);
        R=sqrt(X.^2+Y.^2);
        if Ntarget==1
            rothis=ro(targetrings);
            Rdiff=rothis-R;
            Rdiff(Rdiff>ww)=1; %the whole ring is outside this point
            Rdiff(Rdiff<0)=0; %the whole ring is inside this point
            betweenind=Rdiff<1&Rdiff>0;
            fdiff=Rdiff;
            fdiff(betweenind)=Rdiff(betweenind)*k;
            Q=Nencircling+sum(fdiff(:))/Nsd;
            Qset(P)=Q;
        else
            Q=0;
            for I=1:Ntarget
                rothis=ro(targetrings(I));
                Rdiff=rothis-R;
                Rdiff(Rdiff>ww)=1;
                Rdiff(Rdiff<0)=0;
                betweenind=Rdiff<1&Rdiff>0;
                fdiff=Rdiff;
                fdiff(betweenind)=Rdiff(betweenind)*k;
                Q=Q+sum(fdiff(:))/Nsd;
            end
            Qset(P)=Q+Nencircling;
        end
    end
end
% make full coil by copying & shifting & scaling
ol=[-Ro 0];
or=[Ro 0];
indnz=find(Qset>0);
pcsetnz=pcset(indnz,:);
Qsetnz=Qset(indnz);
ptemp=[bsxfun(@plus,pcsetnz,ol);bsxfun(@plus,pcsetnz,or)];
Qtemp=[-Qsetnz; Qsetnz];
% now make layers & change to meters
eh=wh/Nz;
zlevels=wo+(eh/2:eh:wh);
Wtemp=Qtemp*(ew*1e-3)^2/Nz; %1 mm2 area, Nz layers with the same weight

% and build the pointset...
Npl=numel(Wtemp);
Np=Nz*Npl;
coil.QP=[];
coil.QW=[];
coil.QW=[];
for I=1:Nz
    coil.QP=[coil.QP;[ptemp ones(Npl,1)*zlevels(I)]];
    coil.QW=[coil.QW;Wtemp];
end
coil.QP=coil.QP/1000;
coil.QN=ones(Np,1)*[0 0 1];
coil.QPinds=[1 Np];