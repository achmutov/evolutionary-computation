# Problem description

We are given three columns of integers with a row for each node. The first two
columns contain x and y coordinates of the node positions in a plane. The third
column contains node costs. The goal is to select exactly 50% of the nodes (if
the number of nodes is odd we round the number of nodes to be selected up) and
form a Hamiltonian cycle (closed path) through this set of nodes such that the
sum of the total length of the path plus the total cost of the selected nodes
is minimized.

The distances between nodes are calculated as Euclidean distances rounded
mathematically to integer values. The distance matrix should be calculated just
after reading an instance and then only the distance matrix (no nodes
coordinates) should be accessed by optimization methods to allow instances
defined only by distance matrices.


## How to run

### Quick way

`cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`

`cmake --build build`


### Build the project

`mkdir build && cd build && cmake .. && make`

### Run the solver

`./build/samples/greedy/greedy 200 ./data/TSPA.csv ./data/TSPB.csv`

### Save results to file

`./build/samples/msls_ils/msls_ils 200 ./data/TSPA.csv ./data/TSPB.csv > results.csv`

### Generate visualizations

`python plot.py results.csv`