function [betaQpos,betaQmom]=hbftms_BetaQDipoles(bmeshes)
% HBFTMS_BETAQDIPOLES implements Eqs. 6 and 7 of Stenroos--Koponen, NeuroImage 2019
%
% [betaQpos,betaQmom]=HBFTMS_BETAQDIPOLES(bmeshes)
% [betaQpos,betaQmom]=HBFTMS_BETAQDIPOLES(mesh)
% computes dipoles to be used in quick beta computation as defined in 
% Eqs. 6 and 7 of Stenroos--Koponen, 2019
%   bmeshes = boundary meshes, a cell array of hbf mesh structs
%   mesh = boundary mesh, a hbf mesh struct
%
%   betaQpos = dipole locations (pooled across meshes)
%   betaQmom = dipole moments (pooled across meshes)
%
% v191113 (c) Matti Stenroos, matti.stenroos@aalto.fi

if iscell(bmeshes)
[startind,endind,Nop]=NodeIndices(bmeshes); 
betaQpos=zeros(Nop,3);
betaQmom=zeros(Nop,3);
for M=1:numel(bmeshes)
    [betaQpos(startind(M):endind(M),:),betaQmom(startind(M):endind(M),:)]=...
        BetaQDipoles(bmeshes{M});
end
else
    [betaQpos,betaQmom]=BetaQDipoles(bmeshes);
end

function [Qpos,Qmom]=BetaQDipoles(mesh)
% BETAQDIPOLES computes dipoles to be used in quick beta computation
%
% [Qpos,Qmom]=BETAQDIPOLES(mesh)
%   mesh = triangle mesh structure
%
%   Qpos = dipole locations
%   Qmom = dipole moments
%
% Implements equations (6) and (7) in Stenroos--Koponen, 2019.
% v191007 (c) Matti Stenroos, matti.stenroos@aalto.fi

nop=size(mesh.p,1);
[un,area]=TriangleNormals(mesh.p,mesh.e);
[ntri,ntri_s,ntri_n]=TrianglesForNodes(mesh.e,nop);

Qmom=zeros(nop,3);
for I=1:nop
    Qtemp=[0 0 0];
    for J=1:ntri_n(I)
        tri=ntri(I,J);
        Qtemp = Qtemp +un(tri,:)*area(tri);
    end
    Qmom(I,:)=Qtemp/3;
end

Qpos=zeros(nop,3);
for I=1:nop
    Qpostemp=[0 0 0];
    atot=0;
    for J=1:ntri_n(I)
        tri=ntri(I,J);
        triptemp=[mesh.e(tri,:) mesh.e(tri,:)];
        tipind=ntri_s(I,J);
        trip=triptemp(tipind+(0:2));
        
        p1=mesh.p(trip(1),:)';
        p2=mesh.p(trip(2),:)';
        p3=mesh.p(trip(3),:)';
        
        a=p2-p1;
        b=p3-p1;
        pwthis=p1+[a,b]*[1/4; 1/4];
        Qpostemp=Qpostemp+pwthis'*area(tri);
        atot=atot+area(tri);
    end
    Qpos(I,:)=Qpostemp/atot;
end

function [startinds,endinds,Nop]=NodeIndices(meshes)
Nsurf=length(meshes);
startinds=zeros(Nsurf,1);
endinds=zeros(Nsurf,1);
Nop=0;
for I=1:Nsurf
    startinds(I)=Nop+1;
    endinds(I)=startinds(I)+size(meshes{I}.p,1)-1;
    Nop=endinds(I);
end

function [unormals,areas]=TriangleNormals(nodes,elements)
% function [unormals,areas]=TriangleNormals(nodes,elements)
% function [unormals,areas]=TriangleNormals(mesh)
%
% unormals = unit normals
% areas = areas of the triangles
% v191007 (c) Matti Stenroos

if nargin==1 && isstruct(nodes)
    elements=nodes.e;
    nodes=nodes.p;
end
p1=nodes(elements(:,1),:);
p2=nodes(elements(:,2),:);
p3=nodes(elements(:,3),:);
normals=crossp(p2-p1,p3-p1);
areas=sqrt(sum(normals.*normals,2))/2;
unormals=normals./(2*areas*[1 1 1]);

function [tr,sides,n]=TrianglesForNodes(triangles,maxpind)
% function [tr,sides,n]=TrianglesForNodes(triangles,nop)
% Finds, which triangles belong to the neighborhood of each node
% triangles = triangle description
% maxpin    = maximum vertex index (optional) 
% tr = list of triangles for each node
% sides = index of the node in the triangle (1, 2, or 3)
% n = number of triangles for each node
%
% v191007 (c) Matti Stenroos
if nargin==1
    maxpind=max(triangles(:));
end
maxt=20;%maximum number of triangles for one node (for memory allocation)
tr=zeros(maxpind,maxt);
sides=zeros(maxpind,maxt);
n=zeros(maxpind,1);
not=size(triangles,1);
% count=zeros(non,1);
e1=triangles(:,1);
e2=triangles(:,2);
e3=triangles(:,3);

for I=1:not
    n(e1(I))=n(e1(I))+1;
    tr(e1(I),n(e1(I)))=I;
    sides(e1(I),n(e1(I)))=1;
    
    n(e2(I))=n(e2(I))+1;
    tr(e2(I),n(e2(I)))=I;
    sides(e2(I),n(e2(I)))=2;
    
    n(e3(I))=n(e3(I))+1;
    tr(e3(I),n(e3(I)))=I;
    sides(e3(I),n(e3(I)))=3;
    
end
nmax=max(n);
tr=tr(:,1:nmax);
sides=sides(:,1:nmax);

function R=crossp(R1,R2)
R(:,1)=R1(:,2).*R2(:,3)-R1(:,3).*R2(:,2);
R(:,2)=R1(:,3).*R2(:,1)-R1(:,1).*R2(:,3);
R(:,3)=R1(:,1).*R2(:,2)-R1(:,2).*R2(:,1);
