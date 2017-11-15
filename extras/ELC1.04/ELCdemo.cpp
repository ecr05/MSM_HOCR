/*//////////////////////////////////////////////////////////////////////////////////////////////////
///  ELCdemo.cpp    Clique Reduction by Excludable Local Configuration: Denoising Demo
///  Version 1.04         September 12th, 2014
////////////////////////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////////////////////////

Loads the noise-added image and denoise it with a third-order FoE prior.
The experiment in this demo is described in the following paper:

Hiroshi Ishikawa, "Higher-Order Clique Reduction without Auxiliary Variables,"
In CVPR2014, Columbus, Ohio, June 23-28, 2014.

It also uses techniques described in the following papers:

Hiroshi Ishikawa, "Transformation of General Binary MRF Minimization to the First Order Case,"
IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 33, no. 6, pp. 1234-1249,
June 2011.

Hiroshi Ishikawa, "Higher-Order Clique Reduction in Binary Graph Cut,"
In CVPR2009, Miami Beach, Florida. June 20-25, 2009.

FoE is described in the following paper:

Stefan Roth and Michael J. Black, "Fields of Experts: A Framework for Learning Image Priors,"
In CVPR2005, San Diego, California, June 20-25, 2005, p.p. II:860-867.

This software requires the QPBO software by Vladimir Kolmogorov available
at http://pub.ist.ac.at/~vnk/software.html

This software has been tested on Windows 7 (x64) with Visual Studio 2010,
Ubuntu 12.04 with g++ 4.6.3, and Ubuntu 12.10 with g++ 4.8.1.
Any report on bugs and results of trying on other platforms is appreciated.

//////////////////////////////////////////////////////////////////////////////////////////////////*/

#include "ELC/ELC.h"
#include "QPBO/QPBO.h"
#include "Image.h"
#include <string>
#include <time.h>

using namespace ELCReduce;

// Image filenames
const char* NOISE_ADDED = "test001_020.pgm"; // The given noisy image.
const char* ORIGINAL = "test001.pgm"; // For PSNR calculation
const char* RESULT = "test001_020_result.pgm";  // To save the denoised image.

// Parameters
typedef int REAL; // Type for the energy value (only integer has been tested)
const REAL sigma = 20;	// Sigma of the noise added
const int mult = 10000;	// Multiplier for representing real values by integers
const int Ecycle = 20;	// Number of iterations to check energy decrease.
const int maxIt = 300;	// Maximum number of iteration.
const double stopThreashold = 100.0; // Stop if the energy has changed less than this after Ecycle iterations. 
const int Bcycle = 30;	// Number of iterations to renew the blurred image.

// direction vectors
const int vx[4] = {0, 1, 0, 1};
const int vy[4] = {0, 0, 1, 1};

// FoE experts provided by Stefan Roth 
double alpha[3] = {0.586612685392731, 1.157638405566669, 0.846059486257292};
double expert[3][4] = {
	{-0.0582774013402734, 0.0339010363051084, -0.0501593018104054, 0.0745568557931712},
	{0.0492112815304123, -0.0307820846538285, -0.123247230948424, 0.104812330861557},
	{0.0562633568728865, 0.0152832583489560, -0.0576215592718086, -0.0139673758425540}
};


template<typename T> T square(const T& t) {return t * t;}


double getFoELocal(unsigned int h[4])
{
	double e = 0;
	for (int j = 0; j < 3; j++)
	{
		double s = 0;
		for (int k = 0; k < 4; k++)
			s += expert[j][k] * h[k];
		e += alpha[j] * log(1 + 0.5 * square(s));
	}
	return e;
}


double unaryEnergy(int x, int data, int sigma)
{
	return (double)square(x - data) / (square(sigma) * 2);
}


template<typename F>
void addEnergy(F& f, const image& cur, const image& im, const image& proposal, int W, int H, int N, REAL mult, REAL sigma)
{
	for (int y = 0; y + 1 < H; y++)	// FoE prior
		for (int x = 0; x + 1 < W; x++) 
		{
			int nds[4];
			for (int j = 0; j < 4; j++)
				nds[j] = x + vx[j] + (y + vy[j]) * W;
			REAL E[16];
			for (int i = 0; i < 16; i++) // i = 0000 (binary) means all current, i = 1000 means only node 0 is proposed, 
			{                                      // i = 1001 means node 0 and 3 are proposed, etc.
				unsigned int h[4];
				int b = 8;
				for (int j = 0; j < 4; j++)
				{
					h[j] = (i & b) ? proposal.buf[nds[j]] : cur.buf[nds[j]];
					b >>= 1;
				}
				E[i] = (REAL)(getFoELocal(h) * mult);
			}
			f.AddHigherTerm(4, nds, E);
		}
	for (int j = 0; j < N; j++) // Data term
	{
		double e0 = unaryEnergy(cur.buf[j], im.buf[j], sigma);
		double e1 = unaryEnergy(proposal.buf[j], im.buf[j], sigma);
		f.AddUnaryTerm(j, (REAL)(e0 * mult), (REAL)(e1 * mult));
	}
}


void buildQPBF(QPBO<REAL>& qpbo, const image& cur, const image& im, const image& proposal, int W, int H, int N, REAL mult, REAL sigma, int mode)
{
	PBF<REAL> pbf(W * H * 10);
	addEnergy(pbf, cur, im, proposal, W, H, N, mult, sigma);
	if (mode == 0)
	{	// mode 0: reduce only ELCs, convert the rest with HOCR
		pbf.reduceHigher(); // Reduce ELCs
		PBF<REAL> qpbf(W * H * 10);
		pbf.toQuadratic(qpbf, N); // Reduce what is left with HOCR. N is the ID for new variable.
		pbf.clear(); // free memory
		qpbf.convert(qpbo, N); // copy to QPBO object by V. Kolmogorov. The QPBO object needs to know that at least N variables exist.
	}
	else if (mode == 1) 
	{	// mode 1: reduce all higher-order terms using the approximation
		pbf.reduceHigherApprox();
		pbf.convert(qpbo, N); // copy to QPBO object by V. Kolmogorov. The QPBO object needs to know that at least N variables exist.
	}
	else if (mode == 2) 
	{	// mode 2: use only HOCR
		PBF<REAL> qpbf(W * H * 10);
		pbf.toQuadratic(qpbf, N); // Reduce to Quadratic pseudo-Boolean function. N is the ID for new variable.
		pbf.clear(); // free memory
		qpbf.convert(qpbo, N); // copy to QPBO object by V. Kolmogorov. The QPBO object needs to know that at least N variables exist.
	}
}


double getEnergy(image& f, const image& data, int W, int H, REAL sigma)
{
	double E = 0;
	unsigned int h[4];
	for (int y = 0; y + 1 < H; y++)
		for (int x = 0; x + 1 < W; x++)
		{
			for (int j = 0; j < 4; j++)
				h[j] = f(x + vx[j], y + vy[j]);
			E += getFoELocal(h);
		}
	for (int j = 0; j < H * W; j++)
		E += unaryEnergy(f.buf[j], data.buf[j], sigma);
	return E;
}


double getPSNR(const image& d, const image& o, int N)
{
	int s = 0;
	for (int i = 0; i < N; i++)
		s += square(d.buf[i] - o.buf[i]);
	return 20 * log10(255.0 / sqrt((double)s / N));
}


double rand01()
{
	return ((double)rand()-0.00000000001) / RAND_MAX;
}


const char* modename[] = {"ELC+HOCR", "Approx. ELC", "HOCR"};
int main(int argc, char *argv[])
{
	int mode = 1;
	if (argc == 2)
		mode = atoi(argv[1]);
	double Erec[Ecycle];
	// mode 0: reduce only ELC, convert the rest with HOCR
	// mode 1: reduce all higher-order terms using the approximation
	// mode 2: use only HOCR
	printf("Mode\t%d (%s)\tOriginal image:\t%s\tNoise-added image:\t%s\n", mode, modename[mode], ORIGINAL, NOISE_ADDED);
	image im, org, blr;
	im.readPGM(NOISE_ADDED);
	org.readPGM(ORIGINAL);
	if (im.empty() || org.empty())
	{
		printf("Error. Cannot load the image.\n");
		return 0;
	}
	int W = im.W, H = im.H, N = W * H;
	image cur(im);
	image proposal(W, H);

	printf("Initial Energy =\t%.1f\tInitial PSNR =\t%.4f\n", (double)getEnergy(im, im, W, H, sigma), getPSNR(im, org, N));
	printf("\nIter\tt(sec.)\tlabeled\tEnergy\tPSNR\t#vars\n");
	clock_t begin = clock();
	for (int i = 0; i < maxIt; i++)
	{
		if ((i % Bcycle) == 0)
			cur.gaussianblur(blr, .5625);
		if ((i % 2) == 0)
			for (int j = 0; j < N; j++)
				proposal.buf[j] = std::min((REAL)255, std::max((REAL)0, (REAL)blr.buf[j] + (REAL)((rand01()-0.5) * sigma * 3)));
		else
			for (int j = 0; j < N; j++)
				proposal.buf[j] = (int)(rand01() * 256);
		QPBO<REAL> qpbo(N * 10, N * 20);
		buildQPBF(qpbo, cur, im, proposal, W, H, N, mult, sigma, mode);
		qpbo.MergeParallelEdges();
		qpbo.Solve();
		qpbo.ComputeWeakPersistencies();
		int labeled = 0;
		for (int j = 0; j < N; j++)
		{
			int res = qpbo.GetLabel(j);
			if (res == 1)
				cur.buf[j] = proposal.buf[j];
			if (res >= 0)
				labeled++;
		}
		double psnr = getPSNR(cur, org, N);
		double te = getEnergy(cur, im, W, H, sigma);
		int elapsed = clock() - begin;
		printf("%d\t%.3f\t", i, (float)elapsed / CLOCKS_PER_SEC);
		printf("%.3f\t", labeled * 100.0 / N);
		printf("%.1f\t%.4f\t%d\n", te, psnr, qpbo.GetNodeNum());
		fflush(stdout);
		if (i >= Ecycle && Erec[i % Ecycle] - te < stopThreashold)
			break;
		Erec[i % Ecycle] = te;
	}
	cur.writePGM(RESULT);
	printf("Denoised image written to\t%s\n", RESULT);
	return 0;
}



