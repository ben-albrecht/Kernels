// Chapel's serial stencil implementation
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

// Timeer
var timer: Timer;

// Arrays and Domains
const Dom: domain(2) = {0.. # n, 0.. # n};
var input: [Dom] real;
var output: [Dom] real;

const W: domain(2) = {-R..R, -R..R};
var weight: [W] real;

// Lots of error checking (TODO)

for j in -R..R do {
  for i in -R..R do {
    weight[i, j] = 0.0;
  }
}

for i in 1..R do {
  weight[0, i] = 1.0 / (2.0*i*R);
  weight[i, 0] = 1.0 / (2.0*i*R);
  weight[-i, 0] = -1.0 / (2.0*i*R);
  weight[0, -i] = -1.0 / (2.0*i*R);
}


writeln("Parallel Research Kernels Version ", PRKVERSION);
writeln("Serial stencil execution on 2D grid");
writeln("Grid size            = ", n);
writeln("Radius of stencil    = ", R);
// Hard-coded for now...
writeln("Type of stencil      = star");
writeln("Data type            = double precision");
writeln("Untiled");
writeln("Number of iterations = ", iterations);

// Initialize the input and output arrays
for j in 0.. # n do {
  for i in 0.. # n do {
    input[i,j] = coefx*i+coefy*j;
  }
}

for j in R.. # n-R do {
  for i in R.. # n-R do {
    output[i, j] = 0.0;
  }
}

var broadrange = R.. # (n-2*R);

var works = 0;
for iteration in 0..iterations do {
  // Start timer after warmup iteration
  if (iteration == 1) {
    timer.start();
  }

  for j in broadrange do {
    for i in broadrange do {
      for jj in -R..R do {
        output[i, j] += weight[0, jj]*input[i, j+jj];
        works += 1; }
      for ii in -R..-1 do {
        output[i, j] += weight[ii, 0]*input[i+ii, j];
        works += 1; }
      for ii in 1..R do {
        output[i, j] += weight[ii, 0]*input[i+ii, j];
        works += 1; }
    }
  }

  // add constant to solution to force refresh of neighbor data, if any
  for j in 0.. # n do {
    for i in 0.. # n do {
      input(i, j) += 1.0;
      works += 1;
    }
  }

} // end of iterations

timer.stop();

// Timings
var stencilTime = timer.elapsed();
writeln("stencil_time: ", stencilTime);

var norm = 0.0;

// compute L1 norm in parallel
for j in broadrange do {
  for i in broadrange do {
    norm += abs(output[i, j]);
  }
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

if debug then
  writeln("Units of work executed: ", works);

