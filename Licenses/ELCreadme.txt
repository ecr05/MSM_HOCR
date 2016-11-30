Clique Reduction by Excludable Local Configuration software
Version 1.04 (September 12th, 2014)

Copyright 2014 Hiroshi Ishikawa. All rights reserved.
This software can be used for research purposes only.
This software or its derivatives must not be publicly distributed
without a prior consent from the author (Hiroshi Ishikawa).

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

For the latest version, check: http://www.f.waseda.jp/hfs/indexE.html
----------------------------------------------------------------------------------------------

OVERVIEW

Files in this distribution.

readme.txt        This file.
ELC/ELC.h         The main ELC software. Include this to use the ELC software.
ELC/ELC0.h        Helper classes and functions used in ELC.h.
ELCRdemo.cpp      A denoisng demo program using ELC.
Image.h           A minimal implementation of a PGM image class.
ELCdemo.vcxproj   A Microsoft Visual C++ 2010 project file.
ELCdemo.sln       A Microsoft Visual C++ 2010 solution file.
test001.pgm       A test image.
test001_020.pgm   Noise-added version of the test image (sigma=20).


Software for minimizing a higher-order function of binary variables x_1,...,x_n.
What it actually does is to reduce the function into a first-order MRF, or a 
Quadratic Pseudo-Boolean function, i.e., a function of the form
E(x_1, ..., x_n, ..., x_m) = \sum_i Ei(x_i) + \sum_{ij} Eij(x_i,x_j)
possibly with additional variables.
Once the function is reduced, minimization itself can be done with known first-order methods,
such as the QPBO (roof dual) algorithm below.

In this software, there are two ways to reduce a higher-order monomial to first order (second degree).


A) Finding Excludable Local Configuration (ELC)
Reduces the monomial without adding any variable. For some monomials, it cannot be done.
So the B) below must be used after A). The technique is described in the following paper:

[1] Hiroshi Ishikawa. "Higher-Order Clique Reduction without Auxiliary Variables."
In CVPR2014, Columbus, Ohio. June 23-28, 2014.

Additionally, this software provides an approximation that always reduces the term
without additional variables. It does so correctly most of the time, but there is no guarantee.
It is much faster, though. This is also described in [1].


B) Higher-Order Clique Reduction (HOCR)
Additional variables are added to reduce the order of the energy.
The number of variables increases exponentially as the order of the given energy increases.
The technique is described in the following papers:

[2] Hiroshi Ishikawa. "Transformation of General Binary MRF Minimization to the First Order Case."
IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 33, no. 6, pp. 1234-1249,
June 2011.

[3] Hiroshi Ishikawa. "Higher-Order Clique Reduction in Binary Graph Cut."
In CVPR2009: IEEE Computer Society Conference on Computer Vision and Pattern Recognition,
Miami Beach, Florida. June 20-25, 2009.

This software is implemented so that it can be used most conveniently 
in combination with the QPBO software by Vladimir Kolmogorov available
at http://pub.ist.ac.at/~vnk/software.html

This software has been tested on Windows 7 (x64) with Visual Studio 2010,
Ubuntu 12.04 with g++ 4.6.3, and Ubuntu 12.10 with g++ 4.8.1.
Any report on bugs and results of trying on other platforms is appreciated.


CHANGE LOG

Version 1.00 (June 23rd, 2014)  First release.
Version 1.01 (July  3rd, 2014)  Modified so that it runs on Linux.
Version 1.02 (August 7th, 2014) Modified the interface for toQuadratic and convert function so that
   the user now has to provide the variable ID used so far. (In the rare case of some variables
   dissapearing from the PBF, using maxID() to determine the ID for new variable (in the case of
   toQuadratic) or the largest ID that is used (in the case of convert) can lead to an incorrect
   behavior.
Version 1.03 (September 5th, 2014) Bug fix.
Version 1.04 (September 12th, 2014) Bug fix.


