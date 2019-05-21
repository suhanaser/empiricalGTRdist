# non-SRH simulations pipeline
_____________________________________
*All scripts have hard-coded input and output destinations.*
*If you want to run them for yourself, you will need to adjust these destinations in each script as you go.*

##### 1. run birth-death.R
This script simulates birth-death trees under the constant rate birth-death process with incomplete sampling.
The speciation rate, extinction rate, and sampling fraction are sampled from a unifrom distribution (0,1).
However, we force the extinction rate to be smaller than the speciation rate.
