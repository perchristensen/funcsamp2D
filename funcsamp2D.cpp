//
// funcsamp2D.cpp
// Sample a 2D function on the unit square domain with various sample sequences.
// The result is a table of sampling errors for increasing sample numbers.
// Per Christensen, 2015-2019.
// 
// -----------------------------------------------------------------------------
// Copyright (c) 2019 Disney/Pixar
// 
// Permission is hereby granted to use this software solely for non-commercial applications
// and purposes including academic or industrial research, evaluation and not-for-profit media
// production. All other rights are retained by Pixar. For use for or in connection with
// commercial applications and purposes, including without limitation in or in connection with
// software products offered for sale or for-profit media production, please contact Pixar at
// tech-licensing@pixar.com.
// 
// THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
// NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PIXAR OR ITS AFFILIATES BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// -----------------------------------------------------------------------------
//
// To compile (debug or optimized):
// g++ -Wall -o funcsamp2D funcsamp2D.cpp
// g++ -O3 -o funcsamp2D funcsamp2D.cpp
//
// To run:
// funcsamp2D functionName samplesFilename [numSamples numSequences]
// For example:
// funcsamp2D quarterdisk random_1024samples_100sequences.data 1024 100
// funcsamp2D quartergaussian halton_base23_owen_1024samples_100sequences.data 1024 100 >
//                                                                 errors_gaussian_halton23.data
//
// The file format for the input sample file is:
// - A two-line text header
// - A sequence number
//   ... numSamples pairs of floating-point numbers ...
// - Next sequence number 
//   ... numSamples pairs of floating-point numbers ...
// ... etc.
//
// For example:
// 
//   // Table of 100 sequences of 1024 uniform random 2D samples
//   // Each sample is generated with drand48().
//   // Sequence 0:
//   0.000000000000 0.000985394675
//   0.041631001595 0.176642642543
//   0.364602248391 0.091330612112
//   ...
//   // Sequence 99:
//   0.380256594000 0.504357940608
//   0.202482347415 0.957024735721
//   0.407707287127 0.166596960495
//   ...
//
// The output is a table of the calculated error for sample counts 4, 8, 12, ... numSamples.
// For example:
//   4 0.170000
//   8 0.127500
//   12 0.101667
//   16 0.086875
//   20 0.086000
//   ...
//   1020 0.012088
//   1024 0.012021
//
// These errors can then be plotted with a plotting program such as Gnuplot or similar.
//
// Feel free to modify this program in any way you want!
//


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>


#define MIN(a,b) ((a < b) ? (a) : (b))
#define MAX(a,b) ((a > b) ? (a) : (b))

#define MAXSAMPLES 4096
#define MAXTABLES 10000
#define NUMFUNCTIONS 18

typedef struct Point { double x, y; } Point;


Point samplePoints[MAXTABLES][MAXSAMPLES];   // sample points read from file

// Known functions and their reference values:
typedef struct Functions {
    const char* name;
    double refValue;
} Functions;

Functions functionTable[NUMFUNCTIONS] =
{
    // Discontinuous 2D functions:
    {"quarterdisk", 0.5},
    {"fulldisk", 0.5},
    {"triangle", 0.5},
    // Piece-wise linear 2D functions:
    {"quarterdiskramp", 0.505273},   // ref value ~ 0.505273
    {"fulldiskramp", 0.505273},   // ref value ~ 0.505273
    {"triangleramp", 0.5},
    // Smooth 2D functions:
    {"quartergaussian", 0.55774629},   // = pi/4 erf(1)^2
    {"fullgaussian", 0.851121},   // = 4 * pi/4 * erf(0.5)^2 = pi * erf(0.5)^2
    {"bilinear", 0.25},
    {"biquadratic", 1.0/9.0},   // ref value = 1/9
    {"sinxy", 0.0},
    {"sininvr", -0.220242},
    // Discontinuous 1D function:
    {"stepx", 1.0/M_PI},
    // Piece-wise linear 1D function:
    {"rampx", 0.3},
    // Smooth 1D functions:
    {"lineary", 0.5},
    {"gaussianx", 0.74682413},   // = sqrt(pi)/2 erf(1)
    {"siny", 2.0/M_PI},
    {"sin2x", 0.0},
};


// Skip a line in file (usually a comment)
static void
skipLine(FILE* fd) {
    char c;
    do {
        fscanf(fd, "%c", &c);
    } while (c != '\n');
}


// Evaluate quarter-disk function at (x,y).  Discontinuous.  The quarter-disk is centered
// at (0,0) and has radius sqrt(2/pi).  The area of the quarter-disk is 0.5.
// Returns 1 if point radius < sqrt(2/pi), 0 otherwise.
double
quarterdisk(double x, double y)
{
    const double radius2 = 2.0 / M_PI;
    double r2;

    r2 = x*x + y*y;
    return (r2 < radius2) ? 1.0 : 0.0;
}


// Evaluate disk function at (x,y).  Discontinuous.  The disk is centered at (0.5,0.5) and
// has radius 1/sqrt(2pi).  The area of the disk is 0.5.
// Returns 1 if point radius < 1/sqrt(2 pi), 0 otherwise.
double
fulldisk(double x, double y)
{
    const double radius2 = 1.0 / (2.0 * M_PI);
    double r2;

    x -= 0.5;
    y -= 0.5;
    r2 = x*x + y*y;
    return (r2 < radius2) ? 1.0 : 0.0;
}


// Evaluate triangle function f(x,y) = (y > x).  Discontinuous.
double
triangle(double x, double y)
{
    return (x + y < 1.0) ? 1.0 : 0.0;
}


// Evaluate quarter-disk ramp function at (x,y).  Piece-wise linear.
// The quarter-disk  is centered
// at (0,0) and has a linear ramp fall-off between 0.7 and 0.9.  The area of the quarter-disk is ???.
double
quarterdiskramp(double x, double y)
{
    const double innerRadius = 0.7, outerRadius = 0.9;
    double r = sqrt(x*x + y*y);

    if (r <= innerRadius)
        return 1.0;
    else if (r >= outerRadius)
        return 0.0;
    else
        return 1.0 - (r - innerRadius) / (outerRadius - innerRadius);
}


// Evaluate disk ramp function at (x,y).  Piece-wise linear.  The disk is centered at (0.5,0.5) and
// has a linear ramp fall-off between 0.35 and 0.45.  The area of the disk is ???.
double
fulldiskramp(double x, double y)
{
    const double innerRadius = 0.35, outerRadius = 0.45;

    x -= 0.5;
    y -= 0.5;
    double r = sqrt(x*x + y*y);
   if (r <= innerRadius)
        return 1.0;
    else if (r >= outerRadius)
        return 0.0;
    else
        return 1.0 - (r - innerRadius) / (outerRadius - innerRadius);
}


// Evaluate triangle ramp function.  Piece-wise linear.
double
triangleramp(double x, double y)
{
    double ymx = 5.0 * (y - x);

    if (ymx >= 0.5)
        ymx = 0.5;
    else if (ymx <= -0.5)
        ymx = -0.5;

    return ymx + 0.5;
}


// Evaluate 2D Gaussian function e^(-x^2-y^2) at (x,y).
double
quartergaussian2D(double x, double y)
{
    return exp(-x*x - y*y);
}


// Evaluate 2D Gaussian function centered at (0.5,0.5): e^(-(x-0.5)^2-(y-0.5)^2)
double
fullgaussian2D(double x, double y)
{
    x -= 0.5;
    y -= 0.5;
    return exp(-x*x - y*y);
}


// Evaluate smooth bilinear function f(x,y) = xy.
double
bilinear(double x, double y)
{
    return x*y;
}


// Evaluate smooth biquadratic function f(x,y) = x^2 * y^2.
double
biquadratic(double x, double y)
{
    return x*x*y*y;
}


// Evaluate smooth function f(x,y) = sin(pi*(x+y)).
double
sinxy(double x, double y)
{
    return sin(M_PI * (x+y));
}


// Evaluate sin(pi/r) function.  Mostly smooth but has very large derivatives near (0,0).
double
sininvr(double x, double y)
{
    double r = sqrt(x*x + y*y);

    return (r > 0.0) ? sin(M_PI/r) : 1.0;
}


// Evaluate 1D step function f(x,y) = 1 if x < 1/pi, 0 otherwise.
double
step(double x)
{
    return (x < 1.0/M_PI) ? 1.0 : 0.0;
}


// Evaluate 1D ramp function f(x,y) = 1 if x < 1/pi, 0 otherwise.
double
ramp(double x)
{
    const double left = 0.2, right = 0.4;

    if (x <= left)
        return 1.0;
    else if (x >= right)
        return 0.0;
    else
        return 1.0 - (x - left) / (right - left);
}


// Evaluate linear function f(x,y) = x.
double
linear(double x)
{
    return x;
}


// Evaluate 1D Gaussian function e^(-x^2) at x.
double
gaussian1D(double x)
{
    return exp(-x*x);
}


// Evaluate sin(pi*x) function at x.
double
sinx(double x)
{
    return sin(M_PI * x);
}


// Evaluate sin(2*pi*x) function at x.
double
sin2x(double x)
{
    return sin(2.0 * M_PI * x);
}


// Random value between 0 and 1
double
uniformrandom()
{
    return drand48();
}


// Evaluate function at sample point s from table t
double
evaluateFunction(int functionNum, int t, int s)
{
    Point sample;
    double result = 0.0;

    sample = samplePoints[t][s];

    switch (functionNum) {
    // 2D:
    case 0: result = quarterdisk(sample.x, sample.y); break;
    case 1: result = fulldisk(sample.x, sample.y); break;
    case 2: result = triangle(sample.x, sample.y); break;
    case 3: result = quarterdiskramp(sample.x, sample.y); break;
    case 4: result = fulldiskramp(sample.x, sample.y); break;
    case 5: result = triangleramp(sample.x, sample.y); break;
    case 6: result = quartergaussian2D(sample.x, sample.y); break;
    case 7: result = fullgaussian2D(sample.x, sample.y); break;
    case 8: result = bilinear(sample.x, sample.y); break;
    case 9: result = biquadratic(sample.x, sample.y); break;
    case 10: result = sinxy(sample.x, sample.y); break;
    case 11: result = sininvr(sample.x, sample.y); break;
    // 1D:
    case 12: result = step(sample.x); break;
    case 13: result = ramp(sample.x); break;
    case 14: result = linear(sample.y); break;
    case 15: result = gaussian1D(sample.x); break;
    case 16: result = sinx(sample.y); break;
    case 17: result = sin2x(sample.x); break;
    }

    return result;
}


int
main(int argc, char *argv[]) {
    FILE *fd;
    double* sumresults;
    double reference, result, estimate, error, sumerror, aveerror, maxerror;
    double x, y;
    int functionNumber;
    int numSamples = 1024, numSequences = 100;
    int s, t, i;
    int ok, line;
    char *functionName = NULL, *samplesFilename = NULL;

    if (argc < 3 || argc > 5) {
	printf("Usage: funcsamp2D functionName samplesFilename [numSamples numSequences]\n");
	return 1;
    }

    // Find function name in table of known functions
    functionName = argv[1];
    bool match = false;
    for (i = 0; i < NUMFUNCTIONS; i++) {
        match = (strcmp(functionName, functionTable[i].name) == 0);
        if (match) break;
    }

    if (i == NUMFUNCTIONS) {
        printf("Unknown function: '%s'\n", functionName);
        exit(1);
    }

    functionNumber = i;
    reference = functionTable[i].refValue;

    samplesFilename = argv[2];

    if (argc > 3)
        numSamples = atoi(argv[3]);   // number of sample points

    if (argc > 4)
        numSequences = atoi(argv[4]);   // number of sequences (trials)

    // Read tables

    // Open file with tables of sample points
    fd = fopen(samplesFilename, "r");
    if (!fd) {
        printf("cannot open file '%s'\n", samplesFilename);
        exit(1);
    }
    // Skip comments on first 3 lines
    for (line = 0; line < 3; line++) skipLine(fd);

    // Read numSequences sequences with numSamples sample points in each
    for (t = 0; t < numSequences; t++) {
        for (s = 0; s < numSamples; s++) {
            ok = fscanf(fd, "%lf %lf", &x, &y);
            if (ok == -1) break;   // too few sample points?
            samplePoints[t][s].x = x;
            samplePoints[t][s].y = y;
        }
        // Skip newline and one-line comment (Sequence number)
        for (line = 0; line < 2; line++) skipLine(fd);
    }

    // Allocate and init
    sumresults = (double *) malloc(numSequences * sizeof(double));
    for (t = 0; t < numSequences; t++)
        sumresults[t] = 0.0;   // memset?

    // Loop over sample counts 0 .. numSamples-1
    for (s = 0; s < numSamples; s++) {
        // Loop over sequences (aka. "trials")
        sumerror = 0.0;
        maxerror = 0.0;
        for (t = 0; t < numSequences; t++) {
            result = evaluateFunction(functionNumber, t, s);

            sumresults[t] += result;
            estimate = sumresults[t] / (s+1);

            error = fabs(estimate - reference);
            sumerror += error;
            maxerror = MAX(error, maxerror);
        }
        aveerror = sumerror / numSequences; 

        // Print error for 4, 8, 12, 16, ... samples
        if ((s+1) % 4 == 0) {
            printf("%i %f\n", s+1, aveerror);
            fflush(stdout);
        } 
    }

    return 0; // ok
}
