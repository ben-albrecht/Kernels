// Chapel's serial implementation of transpose

use PRK;

config const iterations : int = 100;
config const order : int = 100;
config const tileSize : int = 32;
config const debug: bool = false;

// Timer
var timer: Timer;

var tiled : bool;
var bytes = 2.0 * numBytes(real) * order * order;

const Dom = {0.. # order, 0.. # order},
 tiledDom = {0.. # order by tileSize, 0.. # order by tileSize};
var A, B : [Dom] real;


writeln("Parallel Research Kernels version ", PRKVERSION);
writeln("Serial Matrix transpose: B = A^T");
writeln("Matrix order          = ", order);
if (tileSize < order && tileSize > 0) {
  writeln("Tile size              = ", tileSize);
  tiled = true;
} else {
  writeln("Untiled");
  tiled = false;
}
writeln("Number of iterations = ", iterations);


// Fill original column matrix
for (i,j) in Dom {
  A[i,j] = order*j + i;
}

// Set transpose matrix to known garbage value
B = -1.0;

var work = 0;
// Main loop
for iteration in 0..iterations {
  // Start timer after a warmup lap
  if (iteration == 1) then timer.start();

  if (tiled) {
    for (i,j) in tiledDom {
      for it in i .. # min(order - i, tileSize) {
        for jt in j .. # min(order - j, tileSize) {
          B[jt,it] = A[it,jt];
          work += 1;
        }
      }
    }
  }
  else {
    for (i,j) in Dom {
      B[j,i] = A[i,j];
      work += 1;
    }
  }

} // end of main loop

timer.stop();

// Analyze and output results

var transposeTime = timer.elapsed();


var abserr = 0.0;
for (i,j) in Dom {
  abserr += abs(B[i,j] - (order*i + j));
}

if (debug) {
  writeln("transposeTime = ", transposeTime);
  writeln("work = ", work);
  writeln("Sum of absolute differences: ", abserr);
}

// Error tolerance
const epsilon = 1.e-8;

if (abserr < epsilon) {
  writeln("Solution validates\n");
  var avgtime = transposeTime / iterations;
  writeln("Rate (MB/s): ", 1.0E-06 * bytes / avgtime,
          " Avg time (s): ", avgtime);
} else {
  writeln("ERROR: Aggregate squared error", abserr,
          " exceeds threshold ", epsilon);
}
