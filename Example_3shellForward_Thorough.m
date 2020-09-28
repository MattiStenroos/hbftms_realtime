% This example demonstrates, how to compute TMS-induced E-fields fast using
% the approach presented in Stenroos & Koponen, NeuroImage 2019.
%
% Before going further, please download sample data and an example BEM solver from
% github.com/MattiStenroos/hbf_sampledata
% github.com/MattiStenroos/hbf_lc_p
%
% v200928 Matti Stenroos
%
% (c) Matti Stenroos and Lari M. Koponen, Aalto University, 2020.
% Contact:  matti.stenroos@aalto.fi (computational routines, plot functions)
%           lari.koponen@aalto.fi   (coil optimization, optimized coil model)
%

clear

%% Load the model geometry
% Please download the head geometry data from
% github.com/MattiStenroos/hbf_sampledata
%
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

% This 3-shell model has much less vertices than the 4-compartment
% and 5-compartment models described and used in the paper. Thus, the
% resulting system matrix is smaller, and the use of GPU in the big
% matrix--vector computation does not give as much advantage as in the paper.

% path to the head geometry data --- put your path here:
sampledatapath = '~/gitwork/hbf_sampledata/hbf_samplehead_3shell_wh';
load(sampledatapath,'bmeshes','cortex');
% -> bmeshes, cortex

scalp = bmeshes{end};

% bmeshes: 3 x 1 cell array that contains mesh structs for inner skull,
%          outer skull and scalp
% cortex: mesh struct that contains a cortical mesh
% Helsinki BEM framework mesh struct contains at least
%          mesh.p: vertices, [number of vertices x 3]
%          mesh.e: element description, [number of faces x 3], starts from 1
% The meshes need to have CCW orientation.

% set conductivities for brain, skull and scalp
% (conductivity inside and outside each boundary)
ci = [1 1/50 1]*.33; 
co = [ci(2:3) 0]; % remember, meshes are ordered in -> out

% Load a coil model
coiltemplate = make_coil_moment(9,1); %the 42-dipole coil used in the paper

% Plot head and coil models
    set(figure(1),'outerposition',[0 50 800 700],'name','headmodel');clf;
    patch('faces',cortex.e,'vertices',cortex.p,'facecolor',[1 .7 .7],...
        'facealpha',1, 'edgecolor','k','edgealpha',.3);
    fcols = [1 .7 .7;.7 .7 1;.9 .9 .9];
    for I = 1:3
        patch('faces',bmeshes{I}.e,'vertices',bmeshes{I}.p,'facecolor',fcols(I,:),...
            'facealpha',.3,'edgecolor','none');
    end
    view([-90 0]);
    axis tight equal off; material dull;lighting gouraud;camlight
    
    set(figure(2),'outerposition',[800 50 800 700],'name','coiltemplate');clf;hold on
    plot3(coiltemplate.QP(:,1),coiltemplate.QP(:,2),coiltemplate.QP(:,3),...
        'k.','MarkerSize',10);
    quiver3(coiltemplate.QP(:,1),coiltemplate.QP(:,2),coiltemplate.QP(:,3),...
        coiltemplate.QW.*coiltemplate.QN(:,1),coiltemplate.QW.*coiltemplate.QN(:,2),...
        coiltemplate.QW.*coiltemplate.QN(:,3));
    axis([-.1 .1 -.05 .05 -.05 .05]);axis equal;view([120 30]);
    xlabel('x');ylabel('y');zlabel('z');

%% Compute boundary potentials for dipoles placed in nodes of the cortex mesh
% In the paper, the potentials were computed using ISA-formulated linear
% Galerkin BEM (LGISA), as described in Stenroos & Sarvas, PMB 2012. For
% simplicity & to provide a complete example, we use here a much simpler
% collocation BEM implemented in the Helsinki BEM Framework LCISA solver
% (Stenroos et al., CPMB 2007; Stenroos et al., NeuroImage 2014; Stenroos &
% Nummenmaa, PLoS ONE 2016), shared at 
%
% www.github.com/MattiStenroos/hbf_lc_p
%%
% We are not saying that LCISA BEM would be optimal or even adequate for
% building high-detail TMS E-field models. It is, however, good enough for
% 3-shell models, when the field points are not too close to the meshes
% (see also Nummenmaa, Stenroos et al., Clin Neurophys 2013).
% 
% The more accurate but slower LGISA BEM that was used in the paper is 
% not in release state. It may, however, be available for collaborative use.

% Set here your path to hbf_lc_isa_p
hbfpath = '~/gitwork/hbf_lc_p';

thisdir = cd;
cd(hbfpath);
hbf_SetPaths();
cd(thisdir);

% Hey! Ho! Let's go!
% double-layer operators for BEM
D = hbf_BEMOperatorsPhi_LC(bmeshes);
% BEM transfer matrix, setting ISA to surface 1.
Tphi_full = hbf_TM_Phi_LC_ISA2(D,ci,co,1);
% potential matrix
Phi = single(hbf_LFM_Phi_LC(bmeshes,Tphi_full,cortex.p));
clear D Tphi_full

%% Related to Eq. 4 & 5: weighted potentials

Phiw = hbftms_WeightedPhi(bmeshes,Phi,ci,co);

%% Eq. 6 & 7: dipoles for fast beta integrals

[betaQpos,betaQmom] = hbftms_BetaQDipoles(bmeshes);

%% How to define a coil position & mapping?
% coil position cp : the center point at the bottom of the coil casing,
%                    corresponds to the origin of the coil template.
% coil normal cn: outer normal of the coil, i.e. away from the head
% coil tangent 1 ct1: long axis
% coil tangent 2 ct2: short axis ~ direction of primary E under the coil 

pind = 184;
cp = scalp.p(pind,:); % set the coil in a location on the scalp
cn   =  [-0.670932428422181  -0.087159683501564   0.736378208574522];
ct1  =  [0.472266524180583  -0.815836492774135   0.333729152452088];
ct2  =  [0.571676487735770   0.571676487735770   0.588533760068349];

coil = coiltemplate;
% rotation matrix for the coil coordinates
T = [ct1;ct2;cn];
% apply rotation and translation
coil.QP = coiltemplate.QP*T+cp;
coil.QN = coiltemplate.QN*T;

% plot geometry
    set(figure(3),'outerposition',[0 750 800 700],'name','geometry');clf;hold on;
     hp  =  patch('faces',scalp.e,'vertices',scalp.p,'facecolor',[1 .7 .7],...
        'facealpha',.3, 'edgecolor','k','edgealpha',.3);
    hp  =  patch('faces',cortex.e,'vertices',cortex.p,'facecolor',[.7 .7 1],...
        'facealpha',1, 'edgecolor','none');
    h = plot3(coil.QP(:,1),coil.QP(:,2),coil.QP(:,3),...
        'k.','MarkerSize',15);
    h = quiver3(coil.QP(:,1),coil.QP(:,2),coil.QP(:,3),...
        coil.QW.*coil.QN(:,1),coil.QW.*coil.QN(:,2),...
        coil.QW.*coil.QN(:,3),'linewidth',1.5);
    ht1 = quiver3(cp(1),cp(2),cp(3),ct1(1),ct1(2),ct1(3),.05,'linewidth',2);
    ht2 = quiver3(cp(1),cp(2),cp(3),ct2(1),ct2(2),ct2(3),.05,'linewidth',2);
    hn = quiver3(cp(1),cp(2),cp(3),cn(1),cn(2),cn(3),.05,'linewidth',2);
       view([-90,30]);
     axis tight equal off; material dull;lighting gouraud;camlight
    camzoom(1.5);
    legend([ht1,ht2,hn],{'Tangent 1','Tangent 2','Normal'});

%% Eqs. 5 & 2: how to compute the E-field
% dI/dt  =  6600 A/(100 \mu s), sign comes from reciprocity (in sphere produces ~ 100 V/m)

minusdIPerdt = -6600*10000; 
BetaQ = hbftms_BetaQ(coil.QP,coil.QN,coil.QW,betaQpos,betaQmom);
Bp = hbftms_BpFlux_xyz(coil.QP,coil.QN,coil.QW,cortex.p); %flux of primary B thru the coil
Bs = Phiw*BetaQ; %flux of secondary B thru the coil
Etms = minusdIPerdt*(Bp+Bs); %use reciprocity to convert flux of B to cortical E.

%% How to use GPU for matrix--vector multiplication?
%
% Init GPU
%
% When everything has been set up, the GPU computation is initialized by
% calling function 'gpuDevice'.
%   gpustat = gpuDevice(index);
% where the index is an integer (usually 1) to the desired GPU device.
%
% We give a few pointers to getting to this stage.
%
% First, the script checks if 'Parallel Computing Toolbox' is installed.
%
% Then, it either sets up the computation with the first applicable GPU or,
% if no devices are available, simulates GPU computation with CPU.
%

gpu_device_count = 0;

% Check toolbox version, with function 'ver':
%   Depending on MATLAB version, the toolbox is known either as
%   'parallel' or 'distcomp'.
%   We have succesfully used, for example, versions 6.12 in R2018a and 7.0 in R2019a.
if isempty(ver('parallel'))
    fprintf('MATLAB Parallel Computing Toolbox is not installed.\n');
else
    fprintf('Found MATLAB Parallel Computing Toolbox.\n');
    gpu_device_count = gpuDeviceCount();
    if gpu_device_count == 0
        fprintf('No suitable GPU devices were found.\n');
        fprintf('If you have an nVidia GPU try updating its drivers.\n');
        fprintf('If not, you have to insert one into your computer.\n');
    else
        fprintf('Found %d GPU device(s).\n',gpu_device_count);
    end
end

if gpu_device_count>0
    % Clear GPU simulation, if any.
    clear('gpuArray');
    % Init the real GPU computation with first GPU.
    %   If more than one are installed, you might want to use the others.
    gpustat = gpuDevice(1);
    gpu_device_name = gpustat.Name;
    fprintf('Initialized GPU with model ''%s''.\n',gpu_device_name);
else
    % Simulate GPU computation:
    %   The GPU simulation is realized with anonymous function 'gpuArray' 
    %   that replaces its namesake from 'Parallel Processing Toolbox'.
    %   Instead of copying data to GPU VRAM, it creates copies in RAM.
    gpuArray = @(x) x;
    gpu_device_name = 'CPU simulation';
    fprintf('No GPU acceleration available, "simulating" GPU with CPU!\n');
end
% END initialize GPU computation.

% Compute the E-field with GPU
Phiw_gpu = gpuArray(Phiw);
% dI/dt  =  6600 A/(100 \mu s), sign comes from reciprocity (in sphere produces ~ 100 V/m)
minusdIPerdt = -6600*10000; 
Bs = Phiw_gpu*gpuArray(hbftms_BetaQ(coil.QP,coil.QN,coil.QW,betaQpos,betaQmom)); 
Bp = hbftms_BpFlux_xyz(coil.QP,coil.QN,coil.QW,cortex.p); %flux of primary B thru the coil
Etms_gpu = minusdIPerdt*(Bp+Bs); %use reciprocity to convert flux of B to cortical E.

%% Example: Plot the E-field norm

    Earray = reshape(Etms,3,numel(Etms)/3)';
    Enorm = sqrt(sum(Earray.^2,2));
    cscale = floor(max(Enorm)/10)*10;
    cmap = jet(18);
    cmap = cmap(10:end,:);

    set(figure(4),'outerposition',[850 750 800 700],'name','E-field magnitude');clf;hold on;
      hp  =  patch('faces',cortex.e,'vertices',cortex.p,'facevertexcdata',Enorm,...
          'facecolor','interp','edgecolor','none','facealpha',1);
    colormap(cmap);caxis([0 cscale]);
     h = plot3(coil.QP(:,1),coil.QP(:,2),coil.QP(:,3),...
        'k.','MarkerSize',15);
     ht1 = quiver3(cp(1),cp(2),cp(3),ct1(1),ct1(2),ct1(3),.05,'color','blue','linewidth',2);
    ht2 = quiver3(cp(1),cp(2),cp(3),ct2(1),ct2(2),ct2(3),.05,'color','red','linewidth',2);
    hn = quiver3(cp(1),cp(2),cp(3),cn(1),cn(2),cn(3),.05,'color','blue','linewidth',2);
    view([-90,30]);
    axis tight equal off; material dull;lighting gouraud;camlight
    colorbar;
    camzoom(1.5);

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

% with GPU:
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
    Bs = Phiw_gpu*gpuArray(hbftms_BetaQ(coil.QP,coil.QN,coil.QW,betaQpos,betaQmom));
    Bp = hbftms_BpFlux_xyz(coil.QP,coil.QN,coil.QW,fp);
    Eset_gpu(:,I) = minusdIPerdt*(Bp+gather(Bs));
    
end
fprintf('With GPU (''%s''): solved 3D field in %d points for %d coil positions in %.1f s\n',gpu_device_name,Nfp,Nrot,toc(tstart));

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
