# MSM_HOCR
Multimodal Surface Matching with Higher order Clique Reduction. Source code. Version 3.00 (Dec 2017)

Copyright 2016 Emma C Robinson. This software may only be used or distributed for non-commercial, research purposes. If you have any questions please contact emma.robinson[at]kcl.ac.uk

THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Files in this distribution:

src/newmesh - folding containing surface mesh class, depending on FSL's gifti reader class

src/DiscreteOpt - class for representing the surface mesh registration problem as energies for discrete optimisation

src/MSMRegLib - the main registration class

src/MSM - binaries

extras/ELC1.04 - the ELC library. This is needed for optimisation

Licenses - including copy of the ELC and FastPD licence info.  (for information only; I supply the ELC library with permission from the author (for research use only), FSL have a licence agreement with FastPD and it is released with FSL)

allparameterssulcDRconfStage1to4 - config file for run MSM cortical folding alignment (HCP multimodal parcellation compatible)

config_strain_STRAIN_SPHERE_shear20_within -config file to run MSM cortical folding alignment for optimised alignment of folds (optimised on neonatal data)

config_MSMsulc_pairwise - this is a config file compatible with the original (pairwise) MSM framework. It will run much faster than the new forms, but the methos is less robust to noise and topological variations in the data

Description:

This MSM github repository has been moved from CVS. It provides the source code for MSM specific functions. To run the code will require external libraries:

FSL - you will the source code (available here: https://fsl.fmrib.ox.ac.uk/fsldownloads/fsldownloadmain.html)

ELC - packaged with this repository (extras ELC1.04)

FastPD - FSL have obtained a licence for this and the necesary version of this library will be released in the next version of FSL

To Compile:

see https://github.com/ecr05/MSM_HOCR/blob/master/compilation%20instructions

Overview of the method.

The MSM software series are cortical surface registration tools designed for flexible alignment of multiple different types of data on the cortical surface. MSM_HOCR represents an adaption to the 2014 Neuroimage version of the MSM software that allows improved regularisation through use of the Higher Order Clique Reduction libraries (HOCR and ELC) written by Hiroshi Ishikawa. These improvements have proved essential to achieve the alignments of areal features as described in the HCP's "A Multimodal parcellation of the Human Cerebral Cortex" Glasser, et al. Nature 2016. It is important to note this version MSM performs discrete optimisation implemented using fastPD not QPBO (as quoted in the upcoming MSM paper, and in the ELC READMetxt). The licence for QPBO is not compatible with FSL or this licence. FastPD and ELC are restricted licence libraries which we have gained permission to use as part of the MSM software (see below). For this reason, we are releasing only binaries at this time (for Ubuntu and Centos Linux, and MacOs). Source code will be released as part of the next FSL software library release. We supply an example config file for MSMsulc (MSMall parameters are on the way). Please contact if you have any problems

Usage:

See http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/msm for general usage

Additional config parameters for MSM_HOCR

--regoption can take vales from 1 (original pairwise MSM) to 3. Options reflect different regulariser choices. 2-5 are higher order methods: 2) Angular Deviation Penalty (ADP); used by the HCP parcellation paper; 3) Strain penalty --triclique estimate a triclique data likelihood --shearmod shear modulus parameter for strain likelihood --bulkmod bulk modulus parameter for strain likelihood

A new paper on the methods used in MSM_HOCR is under review. In the meantime please refer to Glocker et al ECCV 2010 (below) for more details on the Higher Order Clique Reduction and the motivations for the ADP and triclique data terms

Contributors:

Emma C. Robinson, Tim Coalson, Kara Garcia Matthew Webster Saad Jbabdi,

Licence information:

To build MSM_HOCR requires libraries with restricted licences (ELC and FastPD). I have provided information on these licences in src/Licences

When using please cite:

Robinson, Emma C., Saad Jbabdi, Matthew F. Glasser, Jesper Andersson, Gregory C. Burgess, Michael P. Harms, Stephen M. Smith, David C. Van Essen, and Mark Jenkinson. "MSM: A new flexible framework for Multimodal Surface Matching." Neuroimage 100 (2014): 414-426.

Robinson, E.C., Garcia, K., Glasser, M.F., Chen, Z., Coalson, T.S., Makropoulos, A., Bozek, J., Wright, R., Schuh, A., Webster, M. and Hutter, J., 2017. Multimodal surface matching with higher-order smoothness constraints. NeuroImage.

Ishikawa, Hiroshi. "Higher-order clique reduction without auxiliary variables." Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition. 2014.

N. Komodakis and G. Tziritas "Approximate Labeling via Graph-Cuts Based on Linear Programming". IEEE Transactions on Pattern Analysis and Machine Intelligence, 2007.

N. Komodakis, G. Tziritas and N. Paragios, "Performance vs Computational Efficiency for Optimizing Single and Dynamic MRFs: Setting the State of the Art with Primal Dual Strategies". Computer Vision and Image Understanding, 2008 (Special Issue on Discrete Optimization in Computer Vision).

Glocker, Ben, et al. "Triangleflow: Optical flow with triangulation-based higher-order likelihoods." European Conference on Computer Vision. Springer Berlin Heidelberg, 2010.
