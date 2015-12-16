// Chapel's data-parallel block-distributed stencil

// TODO - probably should use the dataParTasksPerLocale arg for Block() so that
// I can specify something similar to the idea of "number of threads"

use BlockDist;
use Time;

// Note: Defaulting to STAR stencil (defines weight)
// Configurable runtime constants
config const n: int = 100;
config const iterations: int = 100;
config const debug: bool = false;
config const verbose: bool = false;

// Compile constants
param R = 2;
param coefx = 1.0;
param coefy = 1.0;
param PRKVERSION = "2.15";
param epsilon = 1.e-8;

// Runtime constants
const activePoints = (n-2*R)*(n-2*R);
const stencilSize = 4*R + 1;

// Timer
var timer: Timer;

// Domains
const LocDom = {0.. # n, 0.. # n},
         Dom = LocDom dmapped Block(LocDom),
      BigDom = {-1..n, -1..n} dmapped Block(LocDom),
           W = {-R..R, -R..R},
    InnerDom = {R.. # (n-2*R), R.. # (n-2*R)};

// Arrays (initialized to zeros)
var input:  [BigDom] real;
var output: [BigDom] real;
var weight: [W] real;

// Lots of arg error checking (TODO)

// Setup weight array
for i in 1..R do {
  weight[0, i] = 1.0 / (2.0*i*R);
  weight[i, 0] = 1.0 / (2.0*i*R);
  weight[-i, 0] = -1.0 / (2.0*i*R);
  weight[0, -i] = -1.0 / (2.0*i*R);
}

// Print info
writeln("Parallel Research Kernels Version ", PRKVERSION);
writeln("Serial stencil execution on 2D grid");
writeln("Grid size            = ", n);
writeln("Radius of stencil    = ", R);
writeln("Type of stencil      = star"); // Temporarily hard-coded
writeln("Data type            = double precision");
writeln("Untiled");                     // Temporarily hard-coded
writeln("Number of iterations = ", iterations);

// Setup input array
for j in 0.. # n do {
  for i in 0.. # n do {
    input[i,j] = coefx*i+coefy*j;
  }
}

for iteration in 0..iterations do {

  // Start timer after warmup iteration
  if (iteration == 1) {
    timer.start();
  }

  forall (i,j) in InnerDom do {
    for jj in -R..R do {
      output[i, j] += weight[0, jj]*input[i, j+jj];
    }
    for ii in -R..-1 do {
      output[i, j] += weight[ii, 0]*input[i+ii, j];
    }
    for ii in 1..R do {
      output[i, j] += weight[ii, 0]*input[i+ii, j];
    }
  }

  // Add constant to solution to force refresh of neighbor data, if any
  for (i,j) in LocDom do {
      input[i, j] += 1.0;
  }

} // end of iterations

timer.stop();

// Timings
var stencilTime = timer.elapsed();
writeln("stencil_time: ", stencilTime);

var norm = 0.0;

// compute L1 norm in parallel
for (i,j) in InnerDom do {
  norm += abs(output[i, j]);
}

norm /= activePoints;

/*******************************************************************************
** Analyze and output results.
********************************************************************************/

// Verify correctness
var referenceNorm = (iterations + 1) * (coefx + coefy);

if abs(norm-referenceNorm) > epsilon then {
  writeln("ERROR: L1 norm = ", norm, ", Reference L1 norm = ", referenceNorm);
  exit(1);
} else {
  writeln("Solution validates");
  if verbose then {
    writeln("L1 norm = ", norm, ", Reference L1 norm = ", referenceNorm);
  }
}

var flops = (2*stencilSize + 1) * activePoints;
var avgTime = stencilTime / iterations;
writeln("Rate (MFlops/s): ", 1.0E-06 * flops/avgTime,
        "  Avg time (s): ", avgTime);
