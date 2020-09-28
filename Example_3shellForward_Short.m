% This example demonstrates, how to compute TMS-induced E-fields fast using
% the approach presented in Stenroos & Koponen, NeuroImage 2019.
%
% Before going further, please download sample data and an example BEM solver from
% github.com/MattiStenroos/hbf_sampledata
% github.com/MattiStenroos/hbf_lc_p

% v200928 Matti Stenroos
%
% (c) Matti Stenroos and Lari M. Koponen, Aalto University, 2019.
% Contact:  matti.stenroos@aalto.fi (computational routines, plot functions)
%           lari.koponen@aalto.fi   (coil optimization, optimized coil model)
%

clear

% Load the model geometry.
% The example uses the head geometry of Subject 1 of the data set described in
% Wakeman, D.G. & Henson, R.N. (2015). A multi-subject, multi-modal
% human neuroimaging dataset. Sci. Data 2:150001 doi: 10.1038/sdata.2015.1,
% available through 
% ftp://ftp.mrc-cbu.cam.ac.uk/personal/rik.henson/wakemandg_hensonrn/
%
% The meshes were constructed with SPM, using scripts provided with the
% data. If you use the head geometry in a publication,
% please cite the Wakeman--Henson publication and include the above
% information.
% 
% Please download the head geometry data from
% github.com/MattiStenroos/hbf_sampledata

% path to the head geometry data --- put your path here:
sampledatapath = '~/gitwork/hbf_sampledata/hbf_samplehead_3shell_wh';
load(sampledatapath,'bmeshes','cortex');
% -> bmeshes, cortex

scalp = bmeshes{end}; 

% set conductivities for brain, skull and scalp
% (conductivity inside and outside each boundary)
ci = [1 1/50 1]*.33; 
co = [ci(2:3) 0]; % remember, meshes are ordered in -> out

%% Compute boundary potentials for dipoles placed in nodes of the cortex mesh
% In this example, we use simple and fast LCISA BEM as shared in 
% github.com/MattiStenroos/hbf_lc_p
%
% The more accurate but slower LGISA BEM that was used in the paper is 
% not in release state. It may, however, be available for collaborative use.

% Set here your path to hbf_lc_isa_p
hbfpath = '~/gitwork/hbf_lc_p';

thisdir = cd;
cd(hbfpath);
hbf_SetPaths();
cd(thisdir);

% Hey! Ho! Let's Go!
D = hbf_BEMOperatorsPhi_LC(bmeshes);
Tphi_full = hbf_TM_Phi_LC_ISA2(D,ci,co,1);
Phi = single(hbf_LFM_Phi_LC(bmeshes,Tphi_full,cortex.p));

clear D Tphi_full

%% Related to Eq. 4 & 5: weighted potentials

Phiw = hbftms_WeightedPhi(bmeshes,Phi,ci,co);

%% Eq. 6 & 7: dipoles for fast beta integrals

[betaQpos,betaQmom] = hbftms_BetaQDipoles(bmeshes);

%% Load a coil model

coiltemplate = make_coil_moment(9,1); %the 42-dipole coil used in the paper

%% How to define a coil position & mapping?
% coil position cp : the center point at the bottom of the coil casing,
%   corresponds to the origin of the coil template.
% coil normal cn: outer normal of the coil, i.e. away from the head
% coil tangent 1 ct1: long axis
% coil tangent 2 ct2: short axis ~ direction of primary E under the coil 

%% When everything has been set up, the GPU computation is initialized by
% calling function 'gpuDevice'.
%   gpustat=gpuDevice(index);
% where the index is an integer (usually 1) to the desired GPU device.

gpustat = gpuDevice(1);

% In the thorough script, we give a few pointers to getting to this stage.

%% Speedtest: rotate the coil at 2 deg angles and compute the fields, with and without GPU

minusdIPerdt = -6600*10000; 

% define rotations
rot = 0:2:360;
Nrot = numel(rot);

% pick the ROI
fpind = 1:9000;
fp = double(cortex.p(fpind,:));
Nfp = size(fp,1);

matrixind = reshape([3*fpind-2; 3*fpind-1;3*fpind],Nfp*3,1);
Phiw_gpu = gpuArray(Phiw(matrixind,:));
Eset_gpu = zeros(size(Phiw_gpu,1),Nrot);

% set the coil in a location on the scalp
pind = 184;
cp = scalp.p(pind,:); 
cn   =  [-0.670932428422181  -0.087159683501564   0.736378208574522];
ct1  =  [0.472266524180583  -0.815836492774135   0.333729152452088];
ct2  =  [0.571676487735770   0.571676487735770   0.588533760068349];
coil = coiltemplate;

%with GPU:
tstart = tic;
for I = 1:Nrot
    %tangents for the rotated coil
    ctrot1 = cosd(rot(I))*ct1+sind(rot(I))*ct2;
    ctrot2 = cosd(rot(I)+90)*ct1+sind(rot(I)+90)*ct2;
    %rotation matrix
    Trot = [ctrot1;ctrot2;cn];
    %apply rotation and translation
    coil.QP = coiltemplate.QP*Trot+cp;
    coil.QN = coiltemplate.QN*Trot;
    
    %the actual field computation
    Bs = Phiw_gpu*gpuArray(hbftms_BetaQ(coil.QP,coil.QN,coil.QW,betaQpos,betaQmom));
    Bp = hbftms_BpFlux_xyz(coil.QP,coil.QN,coil.QW,fp);
    Eset_gpu(:,I) = minusdIPerdt*(Bp+gather(Bs));
end
fprintf('With GPU (''%s''): solved 3D field in %d points for %d coil positions in %.1f s\n',gpustat.Name,Nfp,Nrot,toc(tstart));

%% Without GPU:
% define rotations
rot = 0:2:360;
Nrot = numel(rot);

% pick the ROI
fpind = 1:9000;
fp = double(cortex.p(fpind,:));
Nfp = size(fp,1);

matrixind = reshape([3*fpind-2; 3*fpind-1;3*fpind],Nfp*3,1);
Phiw_ROI = Phiw(matrixind,:);
Eset = zeros(size(Phiw_ROI,1),Nrot);

% set the coil in a location on the scalp
pind = 184;
cp = scalp.p(pind,:); 
cn   =  [-0.670932428422181  -0.087159683501564   0.736378208574522];
ct1  =  [0.472266524180583  -0.815836492774135   0.333729152452088];
ct2  =  [0.571676487735770   0.571676487735770   0.588533760068349];
coil = coiltemplate;

tstart = tic;
for I = 1:Nrot
    % tangents for the rotated coil
    ctrot1 = cosd(rot(I))*ct1+sind(rot(I))*ct2;
    ctrot2 = cosd(rot(I)+90)*ct1+sind(rot(I)+90)*ct2;
    Trot = [ctrot1;ctrot2;cn];
    % apply rotation and translation
    coil.QP = coiltemplate.QP*Trot+cp;
    coil.QN = coiltemplate.QN*Trot;
    
    % the actual field computation
    Bs = Phiw_ROI*hbftms_BetaQ(coil.QP,coil.QN,coil.QW,betaQpos,betaQmom);
    Bp = hbftms_BpFlux_xyz(coil.QP,coil.QN,coil.QW,fp);
    Eset(:,I) = minusdIPerdt*(Bp+Bs);
    
end
fprintf('No GPU: solved 3D field in %d points for %d coil positions in %.1f s\n',Nfp,Nrot,toc(tstart));
