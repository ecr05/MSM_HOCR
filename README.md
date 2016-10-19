# MSM_HOCR
Multimodal Surface Matching with Higher order Clique Reduction
Version 1.00 (October 20th, 2016)

Copyright 2016 Emma C Robinson. All rights reserved.
This software is for research purposes only.
If you have any questions please contact emma.robinsom01@gmail.com

THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Description:

The MSM software series are cortical surface registration tools designed for flexible alignment of multiple different types of data on the cortical surface. MSM_HOCR represents an adaption to the 2014 Neuroimage version of the MSM software that allows improved regularisation 
through use of the Higher Order Clique Reduction libraries (HOCR and ELC) written by Hiroshi Ishikawa. These improvements have proved essential to achieve the alignments of areal features as described in the HCP's "A Multimodal parcellation of the Human Cerebral Cortex" Glasser, et al. Nature 2016. MSM uses discrete optimisation implemented using fastPD. FastPD and ELC are restricted licence libraries which we have gained permission to use as part of the MSM software (see below). For this reason, we are releasing only binaries at this time (for Ubuntu and Centos Linux, and MacOs). Source code will be released as part of the next FSL software library release. PLease contact if you have any problems 


Usage:

See http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/msm

Contributors:

Emma C. Robinson
Matthew Webster
Saad Jbabdi

Licence information:

MSM_HOCR uses libraries with restricted licences (ELC and FastPD) for which we have permission for research only distribution:

"I hereby give permission to redistribute the ELC and HOCR library as a complete source code package as well as
 in a binary form, for non-commercial use, as part of the Multimodal Surface Matching (MSM) software. 
 
Hiroshi Ishikawa 
14/10/2016
"
FSL have entered into a licence agreement with the licencees of FastPD. F

When using please cite 

Robinson, Emma C., Saad Jbabdi, Matthew F. Glasser, Jesper Andersson, Gregory C. Burgess, Michael P. Harms, Stephen M. Smith, David C. Van Essen, and Mark Jenkinson. "MSM: A new flexible framework for Multimodal Surface Matching." Neuroimage 100 (2014): 414-426.
Emma C. Robinson,, Ben Glocker, Kara Garcia, Matthew F. Glasser, Antonios Makropoulos, Jelena Bozek, Robert Wright, Andreas Schuh, Matthew Webster, Jana Hutter, Anthony Price, Lucilio Cordero Grand, Emer Hughes, Nora Tusor, Timothy S. Coalson, Philip V. Bayly, David C. Van Essen, Stephen M. Smith, A. David Edwards, Joseph Hajnal, Mark Jenkinson, Daniel Rueckert,. "Multimodal Surface Matching with Higher-Order Smoothness Constraints."  (under review)
Ishikawa, Hiroshi. "Higher-order clique reduction without auxiliary variables." Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition. 2014.
N. Komodakis and G. Tziritas "Approximate Labeling via Graph-Cuts Based on Linear Programming". IEEE Transactions on Pattern Analysis and Machine Intelligence, 2007.
N. Komodakis, G. Tziritas and N. Paragios, "Performance vs Computational Efficiency for Optimizing Single and Dynamic MRFs: Setting the State of the Art with Primal Dual Strategies". Computer Vision and Image Understanding, 2008 (Special Issue on Discrete Optimization in Computer Vision).

