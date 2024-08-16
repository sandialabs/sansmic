![Lines of code](https://sloc.xyz/github/sandialabs/sansmic/?category=code)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
[![GitHub contributors](https://img.shields.io/github/contributors/sandialabs/sansmic.svg)](https://github.com/sandialabs/sansmic/graphs/contributors)

## sansmic

sansmic (or SANSMIC) is research software developed to simulate the 
leaching of salt caverns. Its primary use has been modeling the leaching
for the caverns at the 
[U.S. Strategic Petroleum Reserve](https://www.energy.gov/ceser/strategic-petroleum-reserve)
(SPR). sansmic differs from other leaching software as it implements a 
simultaneous leach+fill simulation which was used in the 1980s during the
construction and original fill of the SPR Bryan Mound site. The primary 
use for sansmic is for modeling liquid petroleum product storage caverns
that use raw water for product withdrawals, and as a comparison point for
newer salt disolution models.

All lower-case 'sansmic' is used as the repository name to differentiate
the software from the program written in FORTRAN;
sansmic is a rewrite of the older code using the C++ and python
programming languages.

### Installation


#### PyPI
To start using sansmic, simply:
``bash
python -m pip install sansmic
``

#### Download a wheel


#### Build from source


### Usage

Once installed, you can use
``bash
sansmic myInputFile.ext
``
to get started.

For more detailed usage and API information, please see
[our documentation][docs].

[docs]: https://dbhart.github.io/

% ### Contributing

% If you're interested in contributing to the development of sansmic,
% please see
% [contributing guidelines](CONTRIBUTING.md) for how to get started.

% [Contributors][contributors] include:
%
% [contributors]: https://github.com/sandialabs/sansmic/graphs/contributors

### License & Copyright

See [LICENSE.md](LICENSE.md) and [COPYRIGHT.md](COPYRIGHT.md).


### Credits

The version of sansmic provided here is a rewrite of a program of the 
same name that was in turn a modification of code that was written in the 
early 1980s by J. Russo to meet specific SPR needs.
That code was at least partially based on the [Solution Mining Research 
Institute (SMRI)](https://www.solutionmining.org/) code SALT77,
written by A. Saberian and A.L. Podio[^1].

SALT77 became SALGAS, version 4 of which is currently available from SMRI.

[^1]: Saberian, A. and A.L. Podio. (1977) "Computer model for describing the development of solution-mined cavities". *In Situ* 1(1). OSTI:7278331
