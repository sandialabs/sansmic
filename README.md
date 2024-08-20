![Lines of code](https://sloc.xyz/github/sandialabs/sansmic/?category=code)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)

## SANSMIC

SANSMIC (or sansmic) is research software developed to simulate the 
leaching of salt caverns.
Its primary use has been modeling the leaching for the caverns at the 
[U.S. Strategic Petroleum Reserve][SPR]
(SPR). SANSMIC differs from other leaching software as it implements a 
simultaneous leach+fill simulation which was used in the 1980s during the
construction and original fill of the SPR Bryan Mound site. The primary 
use for SANSMIC is for modeling liquid petroleum product storage caverns
that use raw water for product withdrawals, and as a comparison point for
newer salt disolution models.

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
if you are building this from scratch, you're probably okay with that
:grin: )


#### PyPI
To install a pre-compiled version of sansmic, use the pip command

    python -m pip install sansmic


#### Download a wheel


#### Build from source


### Usage

Once installed, you can use

    sansmic --help

to get help on how to run sansmic from the command line.
For more detailed usage and API information, please see
[our documentation][docs].


### License & Copyright
See [LICENSE.md](LICENSE.md) and [COPYRIGHT.md](COPYRIGHT.md).


### Credits
The version of **sansmic** provided here is a rewrite of a FORTRAN
program[^5] of the same name that was a significant modification
from the original **SANSMIC** program that was written from
approximately 1980 to 1983 by J. Russo[^4] to meet specific needs
during the SPR's development by the U.S. DOE. *That* code was
at-least-partially based on the **SALT** code written by A. Saberian[^1]
and A. Saberian and A. L. Podio[^2] [^3]. SALT became the [SMRI][SMRI]
software **SALGAS**, the current version of which is available from
[SMRI](https://www.solutionmining.org/software).


[^1]: Saberian, A. (1974). "Numerical Simulation of Development of Solution Mined Salt Cavities". Solution Mining Research Institute. [SMRI:74-0002](https://www.solutionmining.org/historical-research-projects)
[^2]: Saberian, A., & Podio, A.L. (1976). "A numerical model for development of solution mined cavities". Solution Mining Research Institute. [SMRI:76-0004](https://www.solutionmining.org/historical-research-projects)
[^3]: Saberian, A., & Podio, A.L. (1977). "Computer model for describing the development of solution-mined cavities". *In Situ* 1(1). [OSTI:7278331](https://www.osti.gov/biblio/7278331)
[^4]: Russo, J. (1983). "User's manual for the salt solution mining code, SANSMIC". Sandia National Laboratories. [OSTI:5680228](https://www.osti.gov/biblio/5680228)
[^5]: Weber, P.D., Rudeen, D.K., & Lord, D.L. (2014). "SANSMIC Validation". Sandia National Laboratories. [doi:10.2172/1150848](https://doi.org/10.2172/1150848)

[SMRI]: https://www.solutionmining.org/
[SPR]: https://www.energy.gov/ceser/strategic-petroleum-reserve
