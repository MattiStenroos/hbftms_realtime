**Helsinki BEM Framework TMS realtime library**

This hbftms library contains Matlab routines for fast computation of TMS-induced electric field in the cortex. The methods were presented in

Stenroos & Koponen, Real-time computation of the TMS-induced electric field in a realistic head model, NeuroImage, 2019 Dec; 203:116159, doi: 10.1016/j.neuroimage.2019.116159.

The routines contained in this library / distribution package are copyrighted by Matti Stenroos and Lari Koponen. If you produce any results with the help of this library, we kindly ask to cite our publication(s) and acknowledge the use of this library in an appropriate way. The library is shared under GPL v3.0 license.

The methods presented in the above article are implemented in routines
*  hbftms_BetaQ.m
*  hbftms_BetaQDipoles.m
*  hbftms_BpFlux_xyz.m
*  hbftms_WeightedPhi.m

The coil models used in the article can be built with routines

*  make_coil_moment.m
*  make_coil_reference.m

To demonstrate the use of the routines, with or without the GPU, the library contains example scripts 
*  Example_3shellForward_Thorough.m
*  Example_3shellForward_Short.m

that use a three-shell head model and a basic previously-described BEM solver for computing the needed boundary potentials. The head model and the BEM solver are available at 

https://github.com/MattiStenroos/hbf_sampledata

https://github.com/MattiStenroos/hbf_lc_p

In the article, the potentials were computed using ISA-formulated linear Galerkin BEM (LGISA), as described in Stenroos & Sarvas, PMB 2012. For simplicity & to provide a complete example, we use here a much simpler and faster collocation BEM implemented in the Helsinki BEM Framework LCISA solver (Stenroos et al., CPMB 2007; Stenroos & Sarvas, PMB 2012; Stenroos et al., NeuroImage 2014; Stenroos & Nummenmaa, PLoS ONE 2016). 

We are not saying that the LCISA BEM would be optimal or even adequate for building high-detail TMS E-field models --- it will not produce as accurate results as a proper LGISA solver, and it is not suitable, for example, for evaluating the accuracy of our fast TMS computation methods. It is, however, good enough for 3-shell models, when the field points are not too close to the boundary meshes (see also Nummenmaa, Stenroos et al., Clin Neurophys 2013). The more accurate but slower LGISA BEM that was used in the article is not publicly distributed. It may, however, be available for collaborative use.

v200928 (c) Matti Stenroos and Lari Koponen
contact: matti.stenroos@aalto.fi

**Example data description**

The example data, 'hbf_samplehead_3shell_wh' is to be downloaded from another repository,

https://github.com/MattiStenroos/hbf_sampledata
  
It contains meshes of the head geometry of Subject 1 of the data set described in

Wakeman, D.G. & Henson, R.N. (2015). A multi-subject, multi-modal
human neuroimaging dataset. Sci. Data 2:150001 doi: 10.1038/sdata.2015.1,

available through 
ftp://ftp.mrc-cbu.cam.ac.uk/personal/rik.henson/wakemandg_hensonrn/

If you use the head geometry in a publication,
please cite the Wakeman--Henson publication and include the above
information.

The meshes were constructed with SPM, using scripts provided with the
data. Sensor geometries of a 306-channel Neuromag MEG system (as described in
MNE-C software, http://mne.tools/) and 256-channel ABC EEG-layout (as
described by Biosemi, http://www.biosemi.com/download/Cap_coords_all.xls,
were added manually by the author. The sensors and their coregistrations
do not correspond to the sensors used in the of the Wakeman--Henson
dataset.


****