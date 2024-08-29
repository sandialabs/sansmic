![Lines of code](https://sloc.xyz/github/sandialabs/sansmic/?category=code)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)


## SANSMIC

SANSMIC (or sansmic) is research software developed to simulate the 
leaching of salt caverns.
Its primary use has been modeling the leaching for the caverns at the 
U.S. Strategic Petroleum Reserve (SPR). 
SANSMIC differs from other leaching software as it implements a 
simultaneous leach+fill simulation which was used in the 1980s during 
the construction and original fill of the SPR Bryan Mound site. The 
primary use for SANSMIC is for modeling liquid petroleum product 
storage caverns that use raw water for product withdrawals, and as 
a comparison point for newer salt disolution models.

All lower-case 'sansmic' is used as the repository name and remaining
documentation to differentiate this package from the older SANSMIC
program written in FORTRAN - and because all-caps feels very loud. 
The sansmic package provided here is a re-write of the original program 
using the C++ and Python programming languages.


### Installation
The sansmic package requires Python 3.9 or greater and at least the 
numpy and pandas packages.

Installation can be accomplished most easily by using the PyPI. 
It can also be installed by downloading a wheel from the [releases] 
in this repository, or by cloning this repository and building it
yourself. (Because there is C++ code at its core, building it from
scratch requires a C++ compiler and the pybind11 package - but 
if you are building this from scratch, the authors assume you will 
be okay with that :grin: )


#### PyPI wheels
To install a pre-compiled version of sansmic, use the pip command

    python -m pip install sansmic


#### Build from source
To download and build from source, you should clone the repository
(or, fork sansmic and clone your repository) and then install using
the editable (-e) flag.

    git clone https://github.com/sandialabs/sansmic.git
    cd sansmic
    python -m pip install -e .


### Usage
Once installed, you can use

    sansmic --help

and

    sansmic-convert --help

to get help on how to run sansmic commands from the command line.
For more detailed usage and API information, please see
[our documentation][docs].


### Acknowledgements and Credits
There are undoubtedly more people who deserve acknowledgement for
their contributions to the science, and if the current maintainers
have missed anyone, please let us know. 

Much of the research into numerical modeling of leaching has come 
from the solution mining industry or from nuclear waste disposal 
research in the 1970s and 1980s. The Solution Mining Research 
Institute (SMRI) funded research that was done by Ahmad Saberian 
and A. J. Podio in the early 1970s. Saberian wrote the SALT code 
which became the SMRI software product SalGas.

A. J. Russo from Sandia National Laboratories released the first 
version of SANSMIC in 1983. Modifications were made to SANSMIC in 
the 1990s by P. S. Kuhlman. From 2013-2015 the code underwent a 
code review by Paula Weber and David Rudeen, who found that other, 
undocumented changes had been made at some point and documented 
those changes.

Thomas Eyerman performed an independent evaluation of SANSMIC in 
1984, and also generously provided a technical review of a "SANSMIC 
Validation" by Weber, Rudeen and Lord from 2014.

The U.S. Strategic Petroleum Reserve (DOE) has been the primary 
sponsor for sansmic development and current use, and the Cavern 
Integiry group at the SPR Project Management Office in New Orleans 
has been very supportive.


### License & Copyright
See [LICENSE.md](LICENSE.md) and [COPYRIGHT.md](COPYRIGHT.md).

Sandia National Laboratories is a multimission laboratory managed 
and operated by National Technology & Engineering Solutions of Sandia, 
LLC, a wholly owned subsidiary of Honeywell International Inc., for 
the U.S. Department of Energy's National Nuclear Security 
Administration under contract DE-NA0003525.
