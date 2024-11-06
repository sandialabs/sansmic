# CHANGELOG


## v1.0.7 (2024-11-06)

### Bug Fixes

* fix(dat-convert): change default conversion behavior to match UM (#94)

Add additional comments to scenarios converted from .DAT files. Change subsequent-stage cavern SG initialization to match UM-documented behavior rather than the undocumented reset-the-SG behavior that was occurring in the FORTRAN version. ([`932b2d1`](https://github.com/sandialabs/sansmic/commit/932b2d1f09bc8fe600f3e0e6ddeb3e43cd293676))

* fix(dat-convert): change default conversion behavior to match sansmic UM (2015) ([`b3fea1c`](https://github.com/sandialabs/sansmic/commit/b3fea1c7eaa7058d91d717c12889c217542c413b))

### Documentation

* docs: Update scenario documentation

Update to include both -depth and -height configuraiton options for brine injection, production, and interface. ([`d28bdf4`](https://github.com/sandialabs/sansmic/commit/d28bdf4254af3d437b1dfd6bcce2372d5ab37d7a))

* docs: update badges ([`0300c09`](https://github.com/sandialabs/sansmic/commit/0300c0959816d5737c35bbddf26ffb035c06989b))

### Testing

* test: rename regression outputs as txt

Rename regression testing output files to be explicitly .txt to avoid filetype confusion. ([`785716d`](https://github.com/sandialabs/sansmic/commit/785716d985f05a47c81d3ac226dce9916dc43f48))


## v1.0.6-rc.2 (2024-11-04)


## v1.0.6 (2024-11-01)

### Bug Fixes

* fix(StageDefinition): handle null values in heights ([`fde972f`](https://github.com/sandialabs/sansmic/commit/fde972fba3bfb87d20899ac81cde11ccb3c353b8))

* fix(StageDefinition): property error in class definition due to dataclass use ([`cc194c4`](https://github.com/sandialabs/sansmic/commit/cc194c408d6319df2c0e428064a51d282c08b735))

* fix(StageDefinition): test the lack of stage title in basic example ([`f2fa701`](https://github.com/sandialabs/sansmic/commit/f2fa7014addf26ba3e1fa015b96e52a0a1037a05))

* fix(StageDefinition): provide default name for stages when creating a CStage object ([`399ea65`](https://github.com/sandialabs/sansmic/commit/399ea658cfda9cee0ba630573c341803e5090747))

* fix(StageDefinition): provide default name for stages when creating a CStage object ([`4eae2b5`](https://github.com/sandialabs/sansmic/commit/4eae2b58a1c4a02cd47615abab75aa57e16bca68))

### Documentation

* docs: update example with depths not heights ([`4426fc7`](https://github.com/sandialabs/sansmic/commit/4426fc7495f5a4dd1f3ba18b72c8367e67c1b01e))

### Testing

* test(python): improve test coverage for the python code ([`b6d1ad9`](https://github.com/sandialabs/sansmic/commit/b6d1ad9fe05019ba788054984dd689592bbe6ffe))

* test(model): test null values ([`8bc3d61`](https://github.com/sandialabs/sansmic/commit/8bc3d6109a2af404a00e1b6f64972f369d3b72c7))

* test: add tests for StageDefinition class ([`93cc9d5`](https://github.com/sandialabs/sansmic/commit/93cc9d5606440555b972fa2b358fdd0eab7d1e9f))


## v1.0.5-post.1 (2024-10-31)


## v1.0.5 (2024-10-31)


## v1.0.4 (2024-10-29)

### Refactoring

* refactor: update libsansmic cpp code (#79)

* refactor: change the way logging is done, make output directories more explicit (not enforced yet)

* style: apply black formatting

* ci: change testing output to use 'tee' instead of redirect

* test: no subprocess test in ipynb on linux

* test: strip metadata by hand from ipynb

* refactor: modify logging output formats

* refactor: move certain elements out of Model into a BaseModel

* refactor(version): move location of version number to avoid circular imports

* refactor(test): fix test to match refactor of license and copyright text names ([`ed6ddd8`](https://github.com/sandialabs/sansmic/commit/ed6ddd896cffab1e3c6572ecde2a04b0b8b97ca9))


## v1.0.3 (2024-10-23)


## v1.0.2 (2024-10-17)


## v1.0.1 (2024-10-12)

### Bug Fixes

* fix: change workflows for proper execution ([`c4ea2f4`](https://github.com/sandialabs/sansmic/commit/c4ea2f4464cbcbdf8d16a94e392fca404348e48e))

* fix: semantic release workflow syntax errors ([`9f7d215`](https://github.com/sandialabs/sansmic/commit/9f7d215439feb7348068bcae8748452b3737b2d5))

* fix: test semantic release and conventional commits ([`54423b1`](https://github.com/sandialabs/sansmic/commit/54423b194057e493def5456b290a3c4ec2f62098))

### Documentation

* docs: build script failed to check out ([`a39a03f`](https://github.com/sandialabs/sansmic/commit/a39a03faa171c1fa36b955c68db1ebfc17ad7813))

* docs: build script failed to check out ([`df70c39`](https://github.com/sandialabs/sansmic/commit/df70c392ee5732daa6c18ed3ab9a88e3220bb4c4))

* docs: build docs in correct order ([`42dcbb7`](https://github.com/sandialabs/sansmic/commit/42dcbb7f1d062cf04bdcda039c98ff3b2a1cc932))

* docs: build docs with preferred latest version ([`f1a71f1`](https://github.com/sandialabs/sansmic/commit/f1a71f1d84c17a0766fb51d7749fc0a7cf1fdd60))

* docs: get correct documentation when building older versions ([`1ea87da`](https://github.com/sandialabs/sansmic/commit/1ea87daae47fa0a96c8d0256f76a007b4de4021a))

* docs: Update README.md

Add badge for gh-pages documentation status

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`0ef7f87`](https://github.com/sandialabs/sansmic/commit/0ef7f87546c84874473a3f20fedd04a7acf69b3a))

* docs: Update conf.py

update analytics

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`c183738`](https://github.com/sandialabs/sansmic/commit/c183738252625027efadbc645e03781d7cb7d5f0))

* docs: try executing multi-version build process ([`669f1be`](https://github.com/sandialabs/sansmic/commit/669f1be23c24c212daa42265a2a46ca38960f770))

* docs: use current sphinx config even for old documentation ([`8c45da3`](https://github.com/sandialabs/sansmic/commit/8c45da31832faa246f31d0903f2e44c179662daa))

* docs: clean apidocs before building each tag ([`cbadc37`](https://github.com/sandialabs/sansmic/commit/cbadc372e9721d88bd10c91cd9d95c8cf177ba00))

* docs: clean apidocs before building each tag ([`f2a1a82`](https://github.com/sandialabs/sansmic/commit/f2a1a826192537408f1a776884dd5a5f6ba8e910))

* docs: correct subprocess syntax

Signed-off-by: David Hart <dbhart@sandia.gov> ([`f2a0568`](https://github.com/sandialabs/sansmic/commit/f2a056892baac2c1f8c2d9165d5502e70dd707e1))

* docs: set up multi-version documentation build

Signed-off-by: David Hart <dbhart@sandia.gov> ([`a181f3a`](https://github.com/sandialabs/sansmic/commit/a181f3a4a0133ffd62ea4792a66db6b77e1e39b7))

* docs: Update README.md

Update the README on github.

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`11a5215`](https://github.com/sandialabs/sansmic/commit/11a521545d994ea8777ed654fd7c1f40e287e683))

### Refactoring

* refactor: Convert command-line interface to click instead of argparse (#54)

* refactor: convert to click from argparse

* test: correct the tests so that results are returned

* refactor: put back version and copyright cli options and test them

* test: change path to join abspath ([`323c99b`](https://github.com/sandialabs/sansmic/commit/323c99b5bdf3685c0b9af28d50ad6d1dffecb376))

### Unknown

* adding logo images ([`eb66068`](https://github.com/sandialabs/sansmic/commit/eb6606804a19a3324ebedf77ed7638151a477611))

* [StepSecurity] ci: Harden GitHub Actions (#47)

Signed-off-by: StepSecurity Bot <bot@stepsecurity.io> ([`aff8932`](https://github.com/sandialabs/sansmic/commit/aff89323e841a360dcca7220847c26069ca0403e))

* Revert "chore: exclude dependabot action bumps from change log"

This reverts commit 427115349b5d17a22dfe2387623e9295f67e0a50. ([`2433508`](https://github.com/sandialabs/sansmic/commit/2433508b94b379df16507b4925e8905450543a6b))

* Fix example files (#42)

* build: update and test example jupyter/ipython notebooks

Resolves: Missing class call in example files #41

The ipython notebooks were not being tested. This led to an error in the example notebooks which are now fixed. The testing of any ipynb notebooks checked into the repo has also been added.

* build: add nbstripout to pre-commit actions to clean notebooks files

* build(examples): need matplotlib to test the example notebooks

* test: add missing regression testing data file

* build: requirements.txt

matplotlib missing ([`edb6d70`](https://github.com/sandialabs/sansmic/commit/edb6d707fb60827561e8a439bbd1074fa7455fb5))

* Fix pypi upload ([`abd9b67`](https://github.com/sandialabs/sansmic/commit/abd9b6784ee3150ff11a3964ed73ae330c4fbf08))

* Fix pypi upload ([`5b10c75`](https://github.com/sandialabs/sansmic/commit/5b10c75b589b92fd2a289480a8faa7d0eaca7b6b))


## v1.0.0 (2024-10-01)

### Unknown

* 1.0.0

Automatically generated by python-semantic-release ([`b02f15e`](https://github.com/sandialabs/sansmic/commit/b02f15ec018ca9aeeca6c25748f544862f65ac1e))

* Harden runner updates ([`f42c300`](https://github.com/sandialabs/sansmic/commit/f42c30023ea0d616dc7849b42b8920051a939158))

* Patch reset cavern sg on new stage start (#25)

* Bring set-cavern-sg back to old behavior

---------

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`fdb71c2`](https://github.com/sandialabs/sansmic/commit/fdb71c25c1d4534edf858707803e7caf951999d5))

* Update .pre-commit-config.yaml

Pylint and cpplint need way more configuration before being added to pre-commit

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`2a98371`](https://github.com/sandialabs/sansmic/commit/2a9837129e1fe70bec4a626086398db160ceaa05))

* [StepSecurity] Apply security best practices (#26)

Signed-off-by: StepSecurity Bot <bot@stepsecurity.io> ([`a7c0e2c`](https://github.com/sandialabs/sansmic/commit/a7c0e2cbd92e1cdba247ee4e8fc52d7edd3c19aa))

* Update release.yml

Update so that the PyPI installations only occur on a new tagged release.

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`ca2f5a1`](https://github.com/sandialabs/sansmic/commit/ca2f5a1c132a44df747a36f7a20602cbd761c20f))

* Update gh-pages.yml

Publish GHPages action had a condition that would always evaluate as false. This commit fixes that.

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`28614d0`](https://github.com/sandialabs/sansmic/commit/28614d07aa5532fb208c696a970d2b0481ca9251))

* Continuous integration updates (#24)

* Update release.yml

* Add Documentation URL to pyproject.

* Trying the CI build and release

---------

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`5031b1f`](https://github.com/sandialabs/sansmic/commit/5031b1f7d71e4f00fe7c52dd38f7bf49867bd99a))

* Update release.yml ([`c27a0fe`](https://github.com/sandialabs/sansmic/commit/c27a0fef643942e2ff1986d2536c875df5f873db))

* Update README.md

Adding link to documentation.

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`df46a5d`](https://github.com/sandialabs/sansmic/commit/df46a5d4feaaf97b4decc115dabe8b89fdc3b3e6))

* Fixing workflows

Adding license and copyrights ([`964ecd4`](https://github.com/sandialabs/sansmic/commit/964ecd4b94a054af6a890de6ede7750fe8c4567e))

* GitHub pages (#23)

* Create and upload sphinx docs.
* Still needs to fix version switcher.

---------

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`d1905da`](https://github.com/sandialabs/sansmic/commit/d1905da3c15a239897dbda050ce6528d9bcd56a5))

* Update release.yml

Typo that was not caught by github-actions validation, python 3.12 didn't build. This commit does not change any of the source code.

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`ede0514`](https://github.com/sandialabs/sansmic/commit/ede0514508ae1d3103bc562114689891f69b86a4))


## v1.0.0-rc.1 (2024-09-14)

### Unknown

* Adding release action ([`ccad071`](https://github.com/sandialabs/sansmic/commit/ccad0714e34b01c43fceec414ebd9c4f8f9dd5f4))

* Update CI workflow to use all OSes and flag coverage by OS ([`3be7f74`](https://github.com/sandialabs/sansmic/commit/3be7f7430544aa13a5d9b5e0c7e56bb1a06130a4))

* Pull in updates from development branch (#10)

Final initial import for release candidate 1 for v1.0.0. Includes libsansmic and sansmic modules. Includes tests and code coverage configuration. Includes sphinx documentation for UM and API docs, though github pages still needs to be configured.

---------

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`8dc4f33`](https://github.com/sandialabs/sansmic/commit/8dc4f33a21ea038a70b1e5ca04c2c682a53e590f))

* Create codecov.yml

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`cdaf4dc`](https://github.com/sandialabs/sansmic/commit/cdaf4dc80a136f3c251a7493fd2cde1e2172a5c5))

* Create dependabot.yml

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`0ed0a33`](https://github.com/sandialabs/sansmic/commit/0ed0a331743e08a0e77c5a00ecba328ef8297d2c))

* Update README.md

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`e0d6360`](https://github.com/sandialabs/sansmic/commit/e0d6360a9ddc6463732a068c54157119d0487064))

* Update LICENSE

Trying to figure out why GitHub and OpenCSS don't recognize the license

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`f602fc2`](https://github.com/sandialabs/sansmic/commit/f602fc2aa0097b66c3df7aba6e6ff65f415b03f3))

* Pull in smoke tests for build_ext (#8)

* Update to main (#7)

* Remove testing output files

* Fix license text formatting for github to recognize it

* Adding smoke test CXX file, header, and pybind files ([`5789690`](https://github.com/sandialabs/sansmic/commit/57896907544acca6ccabc5d8ef50fc995e4761ca))

* Fix license text formatting for github to recognize it ([`59d6e32`](https://github.com/sandialabs/sansmic/commit/59d6e32224fca56364a3178614f8d3ebe17db0f7))

* Remove testing output files ([`9cffcf1`](https://github.com/sandialabs/sansmic/commit/9cffcf1f80be875f8a3e2b2bbc064598f71907a5))

* First pre-commit hooks update ([`743cace`](https://github.com/sandialabs/sansmic/commit/743cace17bb4c4e1e6b5ac2c5d13b6f99986b3b1))

* Add openssf scorecard badge

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`20040d6`](https://github.com/sandialabs/sansmic/commit/20040d6cd98a63bed9049540a235f89a7de823d2))

* Create scorecard.yml

Add in security scorecard scanning

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`212f5f3`](https://github.com/sandialabs/sansmic/commit/212f5f34985d406284d8d92472914dffaf3e7e6a))

* Adding OpenSSF badge

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`4236b9d`](https://github.com/sandialabs/sansmic/commit/4236b9d8ce1543c0c8d9ea5e41ff07b46e07da6e))

* Documentation updates ([`a461b58`](https://github.com/sandialabs/sansmic/commit/a461b58256042d6e69ed5d6a5597fa89c103472f))

* Fix the readme documentation ([`73c179a`](https://github.com/sandialabs/sansmic/commit/73c179a69e01f4392786f2b8832e936ee775dc88))

* Updates to the app and the model classes ([`35619ea`](https://github.com/sandialabs/sansmic/commit/35619eaa59bb4d2688fd534b188872d9f68581a9))

* published ghost-test wheel to pypi ([`e58e129`](https://github.com/sandialabs/sansmic/commit/e58e129ee3984d9a3d4f1ab8a03e78b9ae9aaab8))

* Initial python classes import and attempt at pypi test ([`d8e8401`](https://github.com/sandialabs/sansmic/commit/d8e8401719f9f24dbc527018bdd1a4985005c458))

* Committing the setup and configuration files. ([`560c6bd`](https://github.com/sandialabs/sansmic/commit/560c6bd158049404a6a271171447f0f9b83e9b0f))

* Update README.md

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`d8ecb09`](https://github.com/sandialabs/sansmic/commit/d8ecb09a007b39bd3f7cef4a71738b250b9fa38c))

* Update README.md

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`d1324fe`](https://github.com/sandialabs/sansmic/commit/d1324feb077dbaf5d77450bfa6976cb91178d90f))

* Update README.md

Fix broken hyperlinks

Signed-off-by: David Hart <dbhart@users.noreply.github.com> ([`4a5613c`](https://github.com/sandialabs/sansmic/commit/4a5613cf949f7930d1794b029afc2c20e6f01ec3))

* Update to base documentation ([`5db6234`](https://github.com/sandialabs/sansmic/commit/5db6234b3d9e5623b41b345d3d4945b530e908c9))

* Initial commit ([`4cf7cf7`](https://github.com/sandialabs/sansmic/commit/4cf7cf73416b10de77640b189e61fced932f15e0))
