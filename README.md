# non-SRH simulations pipeline
_____________________________________
*All scripts have hard-coded input and output destinations.*
*If you want to run them for yourself, you will need to adjust these destinations in each script as you go.*

#### 1. run birth-death.R
This script simulates birth-death trees under the constant rate birth-death process with incomplete sampling. The speciation rate, extinction rate, and sampling fraction are sampled from a unifrom distribution (0,1). However, we force the extinction rate to be smaller than the speciation rate.

#### 2. run run_iqtree.sh
This script will run iqtree on all alignment files to find the best-fit GTR models, base frequencies, and branch lengths for each dataset.
(note there is a threads argument at the top of the script which you should change as appropriate, it also relies on GNU parallel)

#### 3. run extract_GTR_parameters_dist_from_empirical_data.py
This script will create a table with the GTR parameters, base frequencies and proportion of invariant sites from all partitions as have been estimated by IQ-TREE. the output file will be saved in the root directory as **GTRparam.csv**

#### 3. run extract_brach_lengths_dist_from_empirical_data.py
This script will create a table with the branch lengths from all datasets as have been estimated by IQ-TREE. the output file will be saved in the root directory as **BranchLen.csv**

#### 3. run Hetero2_Simulations.py
This script will generate non-SRH GTR matrices that use to simulate alignments under non-SRH conditions. In order to run this script you need the software **Hetero2** that is available on [Github](https://github.com/thomaskf/Hetero/releases). In addition to the non-SRH parameters, this script will assign new branch lengths from the empirical distribution to the trees that were simulated under the birth-death process > run birth-death.R. The simulated alignments will be saved in the root directory.
