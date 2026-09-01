# Changelog

All notable changes to the **Trajectory Analysis for LAMMPS** (`trj_analysis`) project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased] - 2026-09-01
### Fixed
- **FFT 1D Transform & High-Frequency Spectrum (`src/core/fftwlib.f90`):**
  - Fixed scaling in `fftw1` for positive frequency bins up to and including the Nyquist frequency ($i = n/2 + 1$).
  - Corrected FFT zero-padding sizing ($n$) in `fftw1d` to ensure $n/2 \ge nin$, guaranteeing all extracted points represent physical positive frequencies strictly below the Nyquist limit.
  - Resolved spurious high-frequency jumps/divergence in $Z(\omega)$ and $S(q,\omega)$ caused by reading unscaled negative frequency bins when $nin > n/2$.
  - Fixed deallocation in `fftw1d` to include the temporary time grid array `tx`.

## [1.7.1] - 2026-08-25
### Added
- Add structural order parameter analysis examples for BCC, FCC, and SC lattices. debug: Corrected sign error in order params (`6b6e0b0`)
- Add literature references and cubic phase benchmark values to order module header (`d2a2cd3`)

## [1.7.0] - 2026-07-30
### Added
- Upgrade to v1.7 and implement geometric asymmetry-based border atom identification for cluster analysis (`a967a0a`)
- Add asym_threshold parameter to input configuration for geometric surface detection (`cbf9e25`)

### Refactored
- Track and report ASCII I/O performance metrics in simulation output (`7088948`)

## [1.6.9] - 2026-07-28
### Added
- Add print_last_brdconf subroutine to output cluster boundary configurations in LAMMPS format (`8d96e84`)

### Changed
- Increase float precision for rcl and mass output fields (`99bf083`)
- Increased max no. of neighbors (`074783a`)

### Fixed
- Correct units treatment in rdfcomp, and handling of idir (`b81cb1a`)
- Update format string to prevent integer overflow in util.f90 (`9945850`)
- Ensure cluster grid distribution is at least one particle (`d45a613`)
- Update cell error messages and netcdf periodic (`e295c1d`)
- Fix segmentation fault in cluster_analysis when velocities are absent (`bedbc6e`)
- Fix NUMBER OF ATOMS in last_conf.lammpstrj header when filtering species (`340a738`)

### Refactored
- Simplify normalization and output logic in cluster log (`d1945dd`)

## [1.6.8] - 2026-05-19
### Added
- Add molecular mapping support to netcdf input (`3d66ac4`)
- Add ex_mol to mod_common imports in util.f90 (`0cc73da`)

### Changed
- Added charges and energies in last configuration printout (when present) (`767c504`)
- Replace Natoms by Nsites in last conf printout (`b63f6b4`)

### Fixed
- Use mapped atom types when writing configuration output (`f7c8f87`)

### Refactored
- Extract print_last_conf into a separate subroutine (`fe51741`)

## [1.6.7] - 2026-05-13
### Added
- Wrap cluster correlation analysis in run_order check (`ab3244d`)
- Add worm analysis and spine trajectory output (`e6e1d9f`)
- Add worms parameter to INPUT_CL namelist (`6508e9b`)
- Implement multi-worm spine tracking and export (`0de649a`)
- Add debug print statements for spine atom tracking (`e7a7038`)
- Implement smoothed spine coordinates via local COM calculation (`8b37995`)
- Rename folded_worm_analysis to worm_analysis and add branched worm analysis (`c845bed`)
- Implement branched worm analysis and skeletonization (`b343594`)
- Implement centered worm analysis with active contours (`804a615`)
- Implement pruned branched worm analysis subroutine (`b091660`)
- Replace branched pruning with perfect worm analysis (`a9ee3eb`)
- Replace robust analysis with iterative branch extraction (`32f4e3c`)
- Add branch length and width calculation (`eaf44ee`)
- Implement faster branched skeletonization (`cab7c3b`)
- Replace branched analysis with single spine routine (`212d718`)
- Add molecule ID to cluster configuration output formats (`279f6e1`)
- Add missing index field to lastconf output format (`eb6e78b`)

### Changed
- Remove debug print for ntypes and nsp (`062ad3c`)
- Import write_spine_lammpstrj from mod_util (`d2f6c54`)
- Update worm analysis (`0721405`)

### Fixed
- Apply PBC folding based on periodicity instead of dimension index (`f6f6efa`)
- Correct variable names and types in worm analysis loop (`93c391c`)
- Correct cluster indexing and variable type declarations (`477f6b2`)
- Update arguments passed to write_all_spines_lammpstrj (`ff880cf`)
- Correct variable name from indim to ndim in allocation (`d27d0de`)
- Rename head variable to headl to avoid shadowing (`5e32829`)
- Use correct queue index in BFS traversal (`354df50`)
- Rename variables to lengths and widths and remove unused spine_len (`bf8f145`)
- Use cluster index instead of undefined variable in write loop (`67d0c17`)
- Correct format string for cluster configuration output (`0e570da`)
- Update write format strings to match output arguments in util.f90 (`df2f1d6`)
- Update format string widths for cluster configuration output (`250f944`)
- Adjust field widths in cluster and configuration output formats (`c1f5d2e`)

### Performance
- Optimize atom type identification using boolean tracking array (`d619e0d`)

### Refactored
- Move write_spine_lammpstrj from util to clusters module (`ae231e1`)
- Move spine output and deallocation outside cluster loop (`35fd947`)
- Rename head import to head_common in analyze_worm_robust (`cced5c8`)
- Rename head import to head_common in iterative analysis (`acd8905`)
- Rename head variable to head_common in BFS helper (`7cc0cbb`)
- Remove worm_analysis module and associated logic (`d3df805`)

## [1.6.6] - 2026-05-11
### Added
- Add folded worm analysis and clean up neighbor search. Testing phase (`c2d84a2`)
- Add debug print for ntypes and nsp variables (`6467864`)
- Add debug print for tempty in read_nc_cfg (`3deeb6b`)

## [1.6.5] - 2026-05-07
### Added
- Map atom data using ID from LAMMPS file (`a932e2b`)

### Changed
- Minor revision (`5200485`)
- Major directory restructuring performed (`15a2b1b`)

### Refactored
- Move ex-scan.cu to src/core directory (`4cf9ada`)
- Remove unused clusters module and gemini config (`4f2722a`)

## [1.6.4] - 2026-04-30
### Added
- Remove debug print and add comment for Tfact calculation (`fdf68b4`)

### Changed
- Working version with topology analysis and rotational dof for rigid molecules (`bf75a27`)
- Remove debug print statements from log and rdf modules (`32a4c2c`)
- Changes in examples folder (`9d44e5e`)

### Refactored
- Remove moltool.f90 module and associated subroutines (`c636235`)

### Documentation
- Update README with molecular topology and rigid body parameters (`b91de69`)
- Update README limitations and formatting (`86fe516`)
- Update READMEs with rigid molecules example and formatting fixes (`e01909d`)

## [1.6.3] - 2026-04-29
### Added
- Add support for rigid molecule definitions in input and common modules (`7673b87`)

### Changed
- Nmol changed to Nsites, variable rigid introduced (`ecf4e40`)
- Testing species selection (`16f501c`)
- Corrected version by Gemini (`aac94bd`)

## [1.6.2] - 2026-04-29
### Changed
- Restored version with Nmol back (`135b3c8`)

## [1.6.1] - 2026-04-28
### Added
- Add molecule identification and topology support (`3ab554d`)
- Add bounds check for atom type remapping and debug print (`6f403b4`)
- Add debug print statements for atom and species types (`e20d428`)
- Add bounds check for atom type remapping (`9863ffa`)
- Remove dead code and add debug prints for species types (`22f9ca5`)
- Add debug print for atom counts in select_ncdfinput (`124af43`)
- Uncomment call to dprof in densprof.cuf (`c6114b9`)
- Add debug prints to read_nc_cfg in netcdf.f90 (`11de3a3`)
- Add debug print to read_nc_cfg in netcdf.f90 (`901d946`)

### Fixed
- Use correct atom count in netcdf reset_natoms subroutine (`9500bf5`)
- Comment out dprof call in densprof.cuf (`91ebcc6`)
- Remove trailing parenthesis in read_nc_cfg subroutine (`af094e4`)
- Assign original atom types to wtypes in read_nc_cfg (`a222e65`)
- Use natoms_in instead of natoms in select_ncdfinput loop (`cd03800`)
- Use constant unit number for atoms output in netcdf.f90 (`b081356`)
- Correct atom type remapping and tracking in read_nc_cfg (`8c73839`)
- Correct atom type indexing and remapping in read_nc_cfg (`3f79b20`)
- Remove syntax error in read_nc_cfg assignment (`03b0684`)
- Correct variable name in print statement inside read_nc_cfg (`90ffa52`)

### Refactored
- Import natoms_in from mod_common in netcdf.f90 (`12de98a`)

## [1.6.0] - 2026-04-27
### Changed
- Major revision. Molecule identification implemented, variable Nmol now refers to molecules exclusively (`faf5454`)

## [1.5.6] - 2026-04-26
### Fixed
- Bug in atom remapping removed (`b77ba35`)

## [1.5.5] - 2026-04-24
### Added
- Add validation for zslice input values (`7e3799f`)
- Add limitations, usage instructions, and supported units to README (`0ec38a2`)

### Fixed
- Rename gxy_cc to gxy_qq and fix normalization in 2D RDFs (`d17f120`)

### Refactored
- Remove idir parameter and hardcode z-axis confinement (`e05ea34`)

### Documentation
- Update module execution log messages for clarity (`b6d2e21`)
- Update README.md (`dac3864`)
- Updated Readme (`9fb8385`)

## [1.5.4] - 2026-04-24
### Changed
- Bound controls for zslice (`e8af68b`)
- Test ... (`f4deb43`)

### Documentation
- Update features heading title (`9858277`)

## [1.5.3b] - 2026-04-23
### Documentation
- Rewrite features and thermodynamic properties sections in README (`12a0bc2`)
- Update analysis capabilities and output file list in header (`db32418`)

## [1.5.3] - 2026-04-22
### Changed
- More tests to go, without water (`0f3f2d2`)

## [1.5.2r1] - 2026-04-22
### Changed
- All charge-charge correlations incorporated. g_qq(r<sigma) has values, which is incorrect (`bbaace9`)

## [1.5.1] - 2026-04-20
### Changed
- Checked first version of q-q correlations under confinement (`896a543`)

## [1.5.0] - 2026-04-20
### Changed
- Main update g_xy(r) and S_xy(q) implemented for z-confined systems (`917ac0c`)

## [1.4.3] - 2026-04-20
### Fixed
- Error in density profile removed. Refactor density profile computation by removing unnecessary parameters and improving error handling in trajectory analysis. (`b3d33cd`)

### Refactored
- Refactor output formatting and improve thread handling in GPU routines. Adjusted time unit messages for better readability, modified thread count initialization, and streamlined S(Q) calculation logic for clarity and efficiency. (`ef4aea3`)

## [1.4.2] - 2026-04-17
### Added
- Add countpslice and countpslice_d arrays for improved histogram normalization in RDF calculations. Working version for rdfxy (`7cb4fcd`)

### Changed
- Enhance S(Q) calculations by adding countpslice and countsliced parameters in sqf3xy and sqf2xy subroutines for improved histogram normalization. (`061e655`)

## [1.4.1] - 2026-04-17
### Refactored
- Refactor Makefile and clean up debug statements in multiple modules. Partially working version for rdf_xy under confinement (`e7a2e4a`)

## [1.4.0] - 2026-04-17
### Added
- Add twoDsq_in_3D flag to mod_common and update input processing for 2D structure factor in 3D confined systems (`54602e3`)
- Add nslice and zslice variables to mod_common and update input processing for confinement parameters (`feb3fb1`)
- Add support for 2D directional sampling in S(Q) calculations and allocate necessary arrays (`305843b`)
- Add zsliced variable to mod_common and update S(Q) calculation subroutines for 2D directional sampling (`b5e446d`)
- Refactor S(Q) calculation variables for improved 2D directional sampling support. Missing correction of output ... (`e8822b5`)
- New check of zslice (`42cb916`)
- Add gtot variable to mod_common for enhanced data handling (`921ebad`)
- Add gtot and islice variables to mod_common and mod_rdf for enhanced data handling (`fda8e45`)
- Add debug print statement for slice counts in 3D RDF kernel execution (`fc7ff42`)
- Add shared memory usage print statement in 3D RDFxy kernel for debugging (`92e49bd`)
- Add debug print statements for histogram gathering in s2D3_sh kernel (`689d43e`)
- Add debug print statements for histogram normalization and contributions in s2D3_sh kernel (`8cb7679`)
- Add debug print statement for histogram accumulation in s2D3_sh kernel (`afa80cc`)
- Add lsmax information to debug print statements in s2D3_sh kernel (`e090905`)
- Add error handling and debug print statements for initialization in multiple modules (`29ceab9`)

### Changed
- Prepare output of directional Sqs (`a48fda6`)
- Minor corrections (`6859596`)
- Minor chage (`5601f20`)
- Reorder input (`17fccaa`)
- Remove extra allocate sqf (`013d19f`)
- Reformat output for directional Sqs (`b075430`)
- Minor correction (`c3619f9`)
- Minor correction bis (`76d2589`)
- Control of confinement direction (`2ecbfdb`)
- Partly modified rdfxy (`3cd111c`)
- Update rdf_sh subroutine to use device array 'r_d' for particle coordinates (`41775c1`)
- Update condition in print_order subroutine to use twoDstruc_3D for 2D structure factor computation (`48c2b34`)
- Update print_order subroutine to use twoDstruc_3D for 2D structure factor calculations and improve file naming conventions (`9aacce9`)
- Enhance print_order subroutine by adding 'xfi' variable and updating integer declarations for improved clarity (`a82a6ac`)

### Fixed
- Fix array dimension in histomix_xyi assignment for 3D RDF handling (`293cafb`)
- Fix slice index validation in rdfxy_sh kernel to prevent out-of-bounds errors (`5b2ba39`)
- Fix histogram declaration in s2D3_sh kernel to ensure proper memory allocation (`6623d43`)
- Fix histogram declaration in s2D3_sh kernel to ensure shared memory allocation (`9df152d`)
- Fix kernel launch parameters in rdfxy_sh for proper execution (`be67acc`)
- Fix memory allocation and deallocation for 2D RDF histogram in s2D3_sh kernel (`9fb35b9`)

### Refactored
- Refactor rdfxy_sh subroutine parameters for improved clarity and functionality (`f32375a`)
- Remove unused input variable 'r' from rdfxy_sh subroutine to clean up code (`f25fa3b`)
- Refactor rdfxy_sh subroutine parameters for improved clarity and remove unused variable (`13026d7`)
- Refactor rdfxy_sh and rdf_sh subroutines to improve parameter handling and remove unused variables (`c5bbcf8`)
- Refactor print_order subroutine to improve clarity by removing unused variables and simplifying array indexing (`d2d445a`)
- Refactor print_order subroutine to improve file handling and output formatting for 2D and 3D systems (`a3a67d4`)

### Documentation
- Comment out debug print statements for histogram gathering in s2D3_sh kernel (`abc772a`)
- Comment out debug print statement in s2D3_sh kernel (`64ca5f8`)

## [1.3.27] - 2026-03-27
### Added
- Add directory for EAM dynamics examples in README (`f3fbaa7`)
- Update README to include support for "metal" units and add MEAM dynamics example (`fb7acda`)

### Changed
- Final version with metal units. Improved output and a MEAM example (`87d5239`)

## [1.3.26] - 2026-03-27
### Added
- Enhance thermodynamic output formatting by adding support for eV and LJ units, and refactor related calculations for clarity (`91db6f9`)

## [1.3.25] - 2026-03-27
### Changed
- Enhance logging output by adding active analysis modules header and formatting for clarity (`7dc1326`)

### Fixed
- Fix formatting in active modules logging output by correcting underscore character (`21cc8bf`)
- Fix formatting in flow control logging output by removing extra newline character (`a42c415`)

### Refactored
- Refactor cluster handling by removing outlier purge option and updating related logic in input and log modules. Important: cluster merge incorporated (`9fb2324`)

## [1.3.24] - 2026-03-26
### Added
- Add debug print statement for new cluster initiation in BFS expansion (`20a35e6`)
- Remove debug print statement for new cluster initiation in BFS expansion (`51ba855`)
- Add gborder_next for improved cluster expansion and memory management (`6c07487`)
- Add asymmetry check and neighbor verification in cluster analysis (`c7dbb82`)
- Add additional integer variables for enhanced functionality in cluster analysis (`adc90ed`)
- Add debug print statement for border node processing in cluster analysis (`8b6dd64`)
- Add debug print statements for neighbor and offset values in cluster analysis (`d485caa`)
- Add debug print statement for count mismatch in GPU adjacency cell processing (`685f5b1`)
- Add debug print statements for count validation in GPU adjacency and cluster analysis (`7bb4754`)
- Add G-DBSCAN algorithm implementation for cluster identification (`dcde67a`)
- Add border_next array for neighbor tracking in cluster analysis (`13f783f`)
- Add iteration tracking for cluster expansion in GPU BFS (`bc02b00`)
- Add debug print statement for labeled particles in cluster identification (`f974e71`)
- Add check for already labeled particles during cluster expansion (`1afbd09`)
- Add debug output for seed particles during cluster expansion (`0a0ae5f`)
- Add variable 'k' to cluster analysis module for enhanced iteration control (`33026f5`)
- Implement alternative GPU BFS routine for cluster expansion and update active node tracking (`3d050fd`)

### Changed
- Initialize gborder array to false before BFS search in cluster identification (`c31f06e`)
- Use atomicCAS for thread-safe border marking in cluster expansion (`8a9ca1d`)
- < by <= in neighbors (`82c7601`)
- Check adjacency (`dfcefe9`)
- Offset (`632fd74`)
- Adjacency check (`52e5bdd`)
- Adj (`1831d45`)
- Remove deprecated cluster_search subroutine and streamline GPU BFS logic (`644a3d7`)
- Remove unnecessary blank line before cluster_analysis subroutine (`81f795b`)
- Update cluster particle count threshold to use minPts variable for warnings (`0ec040c`)
- Remove unused border and label arrays from cluster expansion logic (`f0d765a`)
- T (`ec43e69`)
- Enhance BFS cluster expansion by inheriting cluster IDs from labeled neighbors (`59a7f76`)
- Update GPU BFS routine to correctly initialize active flags and adjust border node conditions (`4c7ab0f`)

### Fixed
- Fix cluster labeling logic by using glabel for core particle identification (`3e2ae57`)
- Fix cluster labeling logic by initializing label array and updating condition for core particle identification (`4107951`)
- Fix parameter name in gbfs subroutine for correct border processing (`6e49545`)
- Fix neighbor condition checks for cluster expansion logic (`dca8ef1`)
- Correct minor error in scan (`4761d7f`)
- Fix active variable handling in GPU BFS for improved neighbor search (`f7c8034`)
- Fix border node handling in GPU BFS to improve cluster marking logic (`5cd9211`)
- Fix atomicCAS usage for neighbor marking in BFS expansion (`4938854`)
- Fix typo in variable declaration for neighbor tracking in cluster analysis (`22eb277`)
- Fix typo in border_next assignment for cluster analysis (`66f83b0`)
- Fix active_old flag assignment for cluster expansion in BFS (`c7fb3d4`)
- Fix label assignment after BFS completion for accurate cluster labeling (`534034f`)
- Fix variable name in cluster label inheritance for accurate cluster identification (`16e8748`)
- Fix variable names for border and label tracking in BFS expansion (`8f649ba`)

### Performance
- Optimize border marking logic during active cluster expansion (`3129d5e`)

### Refactored
- Refactor memory allocation by removing core arrays and adjusting related logic in cluster and GPU routines (`de08fe8`)
- Refactor cluster size logging and enhance warning for small clusters in analysis (`0f9ee57`)
- Refactor neighbor determination logic and restore atomicCAS check for border marking in cluster expansion (`13284d8`)
- Refactor neighbor processing logic for cluster expansion and improve offset calculation (`1860a1e`)
- Refactor adjacency handling in cluster analysis for improved neighbor verification (`1106336`)
- Refactor BFS expansion logic to simplify neighbor marking in GPU subroutine (`22d2bf4`)
- Refactor border variable types from logical to integer for consistency in cluster analysis and GPU routines (`a800250`)
- Refactor border condition check in BFS expansion for improved clarity and performance (`c874a97`)
- Refactor variable name in cluster processing for clarity and consistency (`65aed7b`)
- Refactor BFS expansion logic to check for unvisited neighbors before marking as border (`088fe52`)
- Refactor gbfs subroutine to streamline parameters and improve border node handling (`c7a2358`)
- Refactor BFS cluster expansion logic to improve clarity and remove debug prints (`490ebaf`)

### Documentation
- Comment out atomicCAS check for border marking in cluster expansion (`c5a46d0`)
- Comment out GPU graph cell method call for neighbor determination (`8c21538`)
- Comment out atomicCAS check for border marking during cluster expansion (`85e8c6f`)
- Comment out label assignment in neighbor marking logic for border processing (`2194984`)
- Comment out unused border tracking assignments in GPU BFS routine (`1912ddc`)

## [1.3.23] - 2026-03-26
### Added
- Enhance logging and output messages for clarity, add total time reporting, and improve GPU timing initialization (`3386721`)
- Add trj_analysis to .gitignore to exclude analysis results from version control (`539ec59`)
- Add trj_analysis to .gitignore (`6ab67b0`)
- Add debug print statement for core particle identification in cluster expansion (`38ba54f`)
- Enhance BFS expansion logic in cluster identification and add debug output for core particle neighbors (`ed9a714`)

### Changed
- Trj_analysis (`80f74fc`)
- Remove trj_analysis binary file from the repository (`a7169e7`)

### Fixed
- Fix formatting by removing unnecessary newline before end of module in common.cuf (`7a31d54`)

### Refactored
- Refactor logging and device property functions for improved clarity and functionality (`5975920`)

## [1.3.22] - 2026-03-23
### Changed
- Enhance input validation and default settings in trajectory analysis module (`3b87606`)

## [1.3.21] - 2026-03-23
### Fixed
- Update configuration description and enhance error handling in NetCDF reading (`1ae2970`)

## [1.3.20] - 2026-03-18
### Added
- Add density variable and update input parameters for enhanced cluster analysis (`39db8f8`)
- Add printout frequency parameter to input file documentation (`e491f53`)

### Changed
- Update output formatting for cluster density to clarify terminology (`b3d8cb7`)

## [1.3.19] - 2026-03-18
### Added
- Add rclcls variable for first minimum determination in cluster analysis (`01af08a`)
- Add progress update separators in log output for better readability (`8015e6d`)
- Add contributing guidelines to enhance collaboration and project clarity (`aaf6282`)

### Changed
- Enhance output message for suggested rcl-cl in cluster analysis to include the used rcl-cl value (`1589034`)
- Improve progress update formatting in log output for better visibility (`8728f69`)
- Improve output formatting for species summary with color coding for better visibility (`d8aff93`)

### Fixed
- Fix data type declaration for cluster properties in clusters.cuf (`0352f43`)
- Fix output formatting in S(Q) calculation messages for improved clarity (`6ef8068`)
- Fix output formatting in progress updates and program header for improved visibility (`4dc6b9f`)
- Fix GPU device property reporting by passing the correct device identifier (`ec59525`)

### Refactored
- Refactor output formatting in cluster analysis and RDF calculations for improved readability (`f4968a5`)
- Refactor code structure for improved readability and maintainability (`c0a3b8d`)

## [1.3.18] - 2026-03-17
### Added
- Add units parameter to RDFcomp subroutine for enhanced functionality (`ff00aa9`)
- Add first minimum determination for cluster analysis cutoff in RDF calculations (`a338422`)
- Add threshold check for peak detection in find_first_minimum_interp subroutine (`a30ae47`)

### Changed
- Enhance RDF computation output with suggested connectivity distance based on units (`5d48aba`)
- Update RDFcomp call to use units from configuration for improved accuracy (`4d01cda`)
- Enhance cluster analysis by adding gclt array for improved density calculations (`7153275`)

### Fixed
- Fix threshold condition in find_first_minimum_interp subroutine for peak detection (`8ed33a2`)

## [1.3.17] - 2026-03-17
### Changed
- Added gtot calculation (`75f9f18`)
- Added rmingr (`b3371bd`)

## [1.3.16] - 2026-03-17
### Changed
- Check for too large rcl (`e37b143`)

## [1.3.15] - 2026-03-17
### Changed
- Reorganizad sample & tools directory (`877d379`)
- Testing cells for large rc (`09a47c7`)

## [1.3.14] - 2026-03-16
### Fixed
- Fix directory structure in examples README (`50fb60e`)

### Documentation
- Improved README.md files and examples/tools added (`d846f90`)
- Update README with SALR system examples (`ca18461`)

## [1.3.13] - 2026-03-16
### Added
- New Makefiles (`0c08459`)
- New examples added (`c8bf8a8`)

### Changed
- Corrected heading for order.dat in 3D (`d9cc682`)
- Added examples for electrolyte within two charged plates (`699b7ac`)
- Added missing result files (`fd5f595`)
- Change .gitignore (`14b9a5d`)
- Added dynamics example and scripts to process files (`9b378ca`)

### Documentation
- Added README.md for examples (`af0884c`)

## [1.3.12] - 2026-03-13
### Changed
- Print message on outliers (`b4ba9ba`)

## [1.3.11] - 2026-03-13
### Added
- Added new examples (`c55a1fa`)

### Changed
- Optional purge of outliers (`c108f28`)

### Fixed
- Corrected error in purge (`3200359`)
- Syntax error corrected (`82803b5`)

## [1.3.10] - 2026-03-13
### Changed
- Sp_labels removed (`f781dcd`)

## [1.3.9] - 2026-03-13
### Added
- Add nstoms to util.f90 (`fa08493`)

### Changed
- Added printing of last conf in vmd compatible format (`749c4c4`)
- 1.3.9.1 removed velocity headers from last conf (`19efbd0`)
- Adding examples in 3D (`ccb8bff`)

### Documentation
- Modified README.md (`3460b3b`)

## [1.3.8] - 2026-03-10
### Changed
- Multiple defaults defined, kmin changed to minPts for coherence (`86ef61f`)

## [1.3.7.1] - 2026-03-10
### Changed
- Use tunits (`d5b28a9`)
- Assing tunits (`66e0210`)
- Added us of module nc_conf (`28d6173`)
- Reorder Makefile (`a0bfe66`)
- Reset rcl infor (`9a91e36`)
- Lines commented (`bad9b4a`)
- Reorder rcl (`ace9337`)
- Nc_conf removed from input (`bcb6c24`)

### Fixed
- Correct output units in rcl (`a9f708d`)
- Compilation error removed ? (`a33726c`)

## [1.3.7] - 2026-03-10
### Changed
- Corrected run_clusters condition in input.f90 (`8a86591`)

## [1.3.6] - 2026-03-10
### Changed
- Define default for minPts for cluster definitions (`53a90b1`)
- Compilation errors removed (`b95a139`)
- More comp errors (`a4566c8`)
- Remove initialization of minPts (`d60b116`)
- Test minPts (`fdf4a2f`)

## [1.3.5] - 2026-03-09
### Changed
- First version without outliers (`9319244`)

## [1.3.4.4] - 2026-03-09
### Changed
- First trial of outlier removal (`6de5879`)
- 1st debug (`b8f8c3d`)

## [1.3.4.3] - 2026-03-09
### Changed
- Errors in density profiles corrected (`53e4e94`)
- Format in fshape.dat (`c0fe15e`)
- Test cldens (`a748044`)
- Test dens (`d44aef9`)

## [1.3.4.2] - 2026-03-09
### Changed
- Minimal changes (`c2a493f`)
- Control superdense clusters (`2152412`)

### Fixed
- Correct for fort.112 error (`0646cf6`)

## [1.3.4.1] - 2026-03-06
### Changed
- Corrections in average radius (`364db77`)
- Removed Makefile for debuging purposes (`32b8be8`)

## [1.3.4] - 2026-03-06
### Changed
- Introduced variable to control computation of cluster geometry and correlations (`abe0a83`)

### Documentation
- Minor change in README (`6d576cb`)

## [1.3.3.2] - 2026-03-06
### Added
- Add compute capability 7.5 yo make file (`cac4f10`)

### Changed
- Added path to libfftw in Makefiles (`483a19d`)
- Test change (`0012cef`)
- More tests (`879dcfb`)
- Mas (`2e35160`)
- More .. (`7b59cc1`)
- Check pwall (`6393970`)
- More (`f9d8f69`)
- Mmore (`2b05b34`)
- K (`f7b040e`)
- Ne (`8a62ab5`)
- Variable geometry controls whether cluster shape is analyzed (`764761c`)
- Corrected allocates (`8646ffd`)
- More additions to geometry (`979c26f`)
- Added conditional variables in Makefile (`6574bf2`)

## [1.3.3.1] - 2026-03-05
### Changed
- Minor change (`635ffcb`)
- More minor changes (`ace8a15`)
- Reshaping of cluster geometry analysis (`660b201`)
- Corr 1 (`9f83469`)
- Rclus (`ed94548`)

### Fixed
- Correct for rcluster (`f1b6895`)
- Correct cluster radius (`b6b4dbd`)

### Documentation
- Added info on environment var in README.md (`b3c24e5`)

## [1.3.3] - 2026-03-02
### Added
- New README updates (`e5ec280`)

### Changed
- Makefile includes path to fftw includes if needed (`9084582`)
- Tools removed from tree (`f47bd0c`)

### Fixed
- Error in density profile calculation apparently corrected (`f53e521`)

### Documentation
- Updated README.md (`c83ba16`)
- Readme update (`b44f1ce`)

## [1.3.2] - 2026-02-17
### Changed
- Comments added to the Makefiles (`db049b7`)

### Refactored
- Minor comments afeter major refactor (`02819b9`)

## [1.3.1] - 2026-02-13
### Changed
- Makefile changed to use environment variables (`573259e`)
- Minor change in gpu (`b301978`)
- Minor change (`41a42b7`)

## [1.3.0] - 2026-03-31
### Added
- Add MEAM SiC example description to README (`982aba9`)

### Changed
- Merge remote-tracking branch 'github' into HEAD (`ed915cb`)
- Enhance documentation across multiple modules with detailed purpose, functionality, and usage descriptions for improved clarity and maintainability (`ee059d6`)
- Enhance output formatting in print_active_modules subroutine to conditionally apply color based on unit (`dc88f11`)

### Refactored
- Refactor code structure for improved readability and maintainability (`834f919`)

### Documentation
- Actualizar README.md (`436d52a`)

## [1.2.5] - 2026-02-13
### Changed
- Improved NetCDF input and examples added to file structure (`6338b9c`)
- Changed to GPL Licensse (`0f2335c`)
- Further comments added (`b153552`)
- Update program header (`f6a9077`)

### Documentation
- Updated README (`9b6254f`)

## [1.2.4] - 2026-02-11
### Added
- Add installation instructions for Nvidia Fortran and NetCDF (`192513b`)

### Changed
- Added control to prevent analysis of configuration with changing N (GCMC) (`99f8c60`)
- Minor changes in comments (`f229727`)

## [1.2.3] - 2026-02-10
### Changed
- More code rearragements (`96538cf`)

## [1.2.2] - 2026-02-10
### Changed
- Basic code reorganization and comments added (`99f1462`)

## [1.2.1] - 2026-02-10
### Added
- Add extra command line arg to specify gpu (`202488e`)

### Changed
- Replace arrays done and gdone by label & glabel (`880dbcd`)

## [1.2.0] - 2026-02-09
### Changed
- Removed unneeded matrices. BFS optimized following R. Lomba implementation (`b8333f0`)

## [1.1.1] - 2026-02-09
### Changed
- Improved BFS search (`ed66bb3`)

## [1.1.0] - 2026-02-05
### Changed
- Implements corrections when starting step==0 (`04e891e`)

## [1.0.4] - 2026-02-05
### Fixed
- Correct output the kmin, only if clusters are computed (`54a3d54`)

## [1.0.3] - 2026-02-03
### Changed
- Errors in densprof corrected (`3feb52a`)

## [1.0.1] - 2026-02-02
### Changed
- Minor corrections.. ongoing (`184f503`)
- Changes in Makefiles (`206bc98`)
- Modify net charge (`cdf3323`)

### Fixed
- V1..0.2 Final clean up. Correct implementation of cluster numbre threshold for colective properties (`d5260a2`)

## [0.7.7] - 2026-01-26
### Changed
- First version with 3D cluster order params (`1363cf6`)

## [0.7.6] - 2026-01-26
### Changed
- Output formats corrected (`2d76652`)

## [0.7.5] - 2026-01-26
### Added
- New version implementing neighbor lists for 2D order params computation (`3362d64`)

### Changed
- Large improvement, order parameters used neighbor lists precalculated (`4eba59c`)
- Minor corrections (`749497c`)

## [0.7.4] - 2026-01-25
### Changed
- Version implementing maximum no. of neighbors checking (`954c47c`)

## [0.7.3] - 2026-01-25
### Changed
- Corrected problem with nnbnd=0 (`ca41ccf`)

## [0.7.2] - 2026-01-25
### Changed
- Random segfault due to uninitialiazed im solved (`4b692e4`)

### Fixed
- This version reproduces order parms for sc, bcc, and fcc lattices. Error when nnbnd=0 (`7fd0ccb`)

## [0.7.1.5] - 2026-01-23
### Changed
- First version producing sensible numbers (still wrong) (`f8d4108`)
- Last changes (`6100351`)

## [0.7.1.4] - 2026-01-21
### Fixed
- Big error discovered ... yll cosmphi cannot be global ! (`c58461f`)

## [0.7.1.3] - 2026-01-21
### Changed
- Corrected dimensions in yll, cosmphi, sinmphi (`fce7fd8`)

## [0.7.1.2] - 2026-01-21
### Changed
- Corrected strange modification of N! (`8fb516e`)

## [0.7.1.1] - 2026-01-20
- Maintenance updates and minor code cleanups.

## [0.7.1] - 2026-01-20
### Changed
- Untested version of 3D order parameters complete (`cc9dd92`)

## [0.7.0] - 2026-01-20
### Changed
- Unstable version with 3D order params (`dc816ca`)

## [0.6.4] - 2026-01-16
### Changed
- Added comments concerning periodic info (`0be9fb9`)
- Starting coding 3D order parameters (`1c638b2`)

### Fixed
- Major changes: correct centers and cluster counting in non-periodic conditions, proper handling of single configuration analysis (`c09ba0e`)

## [0.6.3] - 2026-01-13
### Changed
- Various problems with periodic and non-periodic boundaries corrected (`e6582aa`)

## [0.6.2] - 2025-12-01
### Changed
- Version complete except for 3D order parameters (`9eb73cb`)
- Minor corrections in output (`6f6ad78`)
- Minimal changes adding comments (`f74f936`)
- Changes in Makefile (`b98dfff`)

## [0.6.1.8] - 2025-10-15
### Documentation
- Added changes to README.md (`520fc2c`)

## [0.6.1.7] - 2025-10-09
### Documentation
- Changes in README.MD (`75b2c68`)

## [0.6.1.6] - 2025-10-15
### Added
- Add License File (`7ffc5a2`)

### Changed
- Removed memory leak in cluster profile calculation (`557f9e7`)
- Added cumulative order parameter profiles (`2cef12e`)

### Fixed
- Sintax error corrected in order.cuf (`89d2e84`)
- Sintax error corrected (`2396d9e`)

## [0.6.1.5] - 2025-10-06
### Changed
- Adapt cluster ids in last conf to vmd colors (`d320c14`)
- Minor corrections to output (`b8cd46a`)
- Close lastclconf (`cd403b0`)
- Test credential helper (`c9dfb19`)
- Test again (`bf0b244`)
- Final test (`ee0de11`)

## [0.6.1.4] - 2025-10-02
### Changed
- Corrected errors in gcl(r) when cluster search included very small clusers (`8c490d7`)
- Added final cluster configuration printout in lammps traj format. CLusters id generated by cluster size (`ed6d4ee`)
- Minor changes (`38c8f13`)
- Last cluster conf to be written (`52f80b2`)

### Documentation
- Changes in README.md (`b821261`)

## [0.6.1.2] - 2025-09-19
### Added
- Add checks in dynamics (`8d8c332`)

### Changed
- Changes in Makefiles (`a60a0ae`)
- Print removed (`b468778`)

### Fixed
- Removed bug in the histogram of cluster radii (`03a505f`)

### Documentation
- Changes in README.md and comments in main prog (`86ccbd8`)

## [0.6.1.1] - 2025-08-31
### Changed
- Added consistency checks in input file (`5a8af96`)

## [0.6.1] - 2025-09-23
### Changed
- Corrected minor mistakes and dependency of order parameters on link cell (`6699ad9`)
- 3 Added control on no. of buffers (`4d51eaa`)
- Minor correction in output format (`6c9e2de`)

## [0.6.0] - 2025-08-26
### Changed
- Working version implementing order parameters and their profiles in clusters in 2D systems (`f40b3c4`)

## [5.5.2] - 2025-08-25
### Changed
- First version implementing cluster order params profile. Not working (`25b85d1`)

## [5.5.1] - 2025-08-23
### Changed
- Agrees with lammps. nnbond introduced, as the minimum number of neighbors to compute contribution to phi_m (`64beeda`)

## [0.5.5] - 2025-08-23
### Changed
- First working version with hexatic order analysis in 2D (`b29a963`)

## [0.5.3] - 2025-08-21
### Added
- First changes to implement order params (`09e0a7f`)

## [0.5.2] - 2025-07-04
### Changed
- Some comments added (`54839f1`)

### Fixed
- Corrected error in non-periodic directions (`7533bef`)

### Documentation
- Changes in README file (`6b82267`)

## [0.5.1] - 2025-07-03
### Changed
- Improved species selection (`53c441c`)

## [0.5.0] - 2025-07-01
### Changed
- Major corrections in form factors. Origin repositioning removed in netcdf.f90 (`9f23523`)

## [0.4.8] - 2025-06-24
### Changed
- Distribution of internal cluster densities corrected, no input parameter needed (`92fed26`)

## [0.4.7] - 2025-06-18
### Changed
- Corrected major memory leak in shape distros (`7f829a4`)
- Minor changes (`81b54d5`)

## [0.4.6] - 2025-06-16
### Changed
- First version with analysis of inertia moments (`acd2602`)
- Minor corrections (`41afaa9`)
- Checking problems with cldens. Likely corruption of netcdf file (`0dff4f5`)

## [04.5] - 2025-06-10
### Changed
- Added composition analysis of clusters (`bc8c59c`)

## [0.4.4.1] - 2025-06-09
### Changed
- Expand format of cluster evolution file (`ca9aa78`)

## [0.4.4] - 2025-05-29
### Changed
- Removed sigma, added drho for cluster density analysis, corrected major bugs for cluster analysis for sizes > minclsize (`6e21e01`)
- Changes in output messages (`911d16a`)
- Minor changes (`cd4a25b`)

## [0.4.3.3] - 2025-05-28
### Changed
- Check memory freed and errors in bounds of cl-cl g(r) (`53d0093`)

## [0.4.3.2] - 2025-05-22
### Changed
- Removed errors in computation of cluster kinetic energies and temperatures (`175015d`)

## [0.4.3.1] - 2025-05-22
### Fixed
- Corrected error concerning Tfact (`ce7c633`)

## [0.4.3] - 2025-05-14
### Changed
- Internal cluster g(r) only computed for clusters >= 4*minclsize to avoid singularities (`4fe0672`)

## [0.4.2] - 2025-05-11
### Changed
- Remove warning in window function using result() (`a57bed4`)
- Minor corrections in functions (`b0a8cf5`)

### Fixed
- Tested for correct units. Metal units still missing (`fc1bcaa`)

## [0.4.1.1] - 2025-05-09
### Changed
- Corrected analysis of units lj and other (`36a79e0`)

## [0.4.1] - 2025-05-09
### Changed
- LJ units implemented (`4b020d8`)
- Comments added (`f0a418c`)

## [0.4.01rev] - 2025-05-08
### Changed
- Starting revision for lj units (`87dae68`)

## [0.40.1] - 2025-05-07
### Added
- New correction (`3187d7f`)

### Changed
- Implementation of LJ units warning (`a21b489`)
- Shift coordinate origin to (0,0,0) (`d34be92`)
- Origin in (0,0,0) (`4735cdc`)
- Nspmax increased to 10 (`f6ab04d`)

### Fixed
- Correct last conf in 3D (`962f047`)

## [0.40] - 2025-04-11
### Changed
- Added output for final cluster conf filtered in atom_style charge (`5522cdb`)

### Fixed
- Corrected error in write (`3daf957`)

## [0.39] - 2025-04-10
### Changed
- Charge densities, average charges and net charges finally consistent (`523d9bc`)
- First preliminary version with last conf out filtered for clusters > minclsize (needs to be adapted to atom_style charge) (`67d4e70`)

## [0.38] - 2025-02-14
### Fixed
- Correct output of net charges (`27b45e0`)

## [0.37] - 2025-02-14
### Fixed
- Correct handling of atom indices and charges (`d61c9c9`)

## [0.36.1] - 2025-02-12
### Changed
- Output corrected for large no. of species (`7a5e58e`)
- Set all types selected to contiguous IDs by default (`a42d8b4`)

## [0.36] - 2025-02-12
### Changed
- Version using non consecutive atom IDs (`1510a48`)

## [0.35] - 2025-01-15
### Changed
- To date optimized version, viscosity calculation transferred to CPU (`4b8301d`)
- Added commentaries (`df00faa`)
- Changes to makefile no avx2 (`0b203a7`)

## [0.34] - 2024-11-26
### Changed
- Code with corrected stress correlation+DP version commented out (`7147f43`)
- Change nblock to 100 in dynamics (`bebc192`)

## [0.33] - 2024-11-20
### Changed
- Removed Evans&Davies algorithm, only off diag stress. LAMMPS <Pxx(t)Pxx(0)> is intensive <Pxy(t)(Pxy(0)> \propto 1/V (`08fce38`)
- Double precision test for viscosity (`01b21c4`)
- Corrected multiple errors (`17f0610`)

## [0.32] - 2024-10-29
### Changed
- Use Daivies& Evans viscosity formula: all tensor components (`f9e834d`)

## [0.31.3] - 2024-10-29
### Added
- New origin handling algorithm working (`126af11`)

## [0.31.2] - 2024-10-27
### Added
- New buffer handling algorithm (`238607d`)

## [0.31.1] - 2024-10-25
### Changed
- Modify F(Q,t) to include proper buffer origins (`5056465`)

## [0.31] - 2024-10-24
### Added
- Version with new buffer origins (`8398dc7`)

## [0.2.9] - 2024-10-23
### Changed
- Not working origins of buffers shifted ..Nconf=nbuffer**2 ??? (`5e26db1`)

## [0.2.8] - 2024-10-23
### Changed
- First version with stress self-correlation (`51f55ce`)
- Compiles & runs (`f98c766`)

## [0.2.2.30] - 2024-10-10
### Changed
- Removed pressure calculations, retains input of forces (`4503cc9`)

## [0.2.29] - 2024-10-10
### Changed
- Non working version of pressure calculation (`8d81a4a`)
- Updated values of constants (`8e00390`)

### Fixed
- Error in forces transfer corrected. Virial values still too large (`874f27b`)

## [0.28.1] - 2024-10-09
### Changed
- First changes to compute pressure (`2ea17dc`)

### Documentation
- Readme modified to include pressure (`a7632a1`)

## [0.2.28] - 2024-10-09
### Changed
- Changes in Makefile (`4c4f16f`)
- Adapt Makefile for non avx2 cpus (`858d07c`)

## [0.2.27] - 2024-10-08
### Added
- New comments added. Missing current correlations (`878bb84`)

## [0.2.26] - 2024-10-08
### Changed
- Errors in sp_selected corrected (`b1a5efb`)

## [0.2.25] - 2024-10-08
### Changed
- Version with species selection implemented (`07e6e13`)

## [0.2.24] - 2024-10-07
### Changed
- Added FT's to calculate sqw (`c789c3f`)

## [0.2.22] - 2024-10-06
### Changed
- Full working version with F(Q,t) (`f2ecf1c`)
- V 0.2.23 2D and 3D calculation of S(Q)'s unified (`c13ba6f`)

## [0.2.21] - 2024-10-05
### Added
- Add fkt.dat (`7f401c1`)

### Changed
- F(Q,t) added (`843b7b7`)

## [1.0.0] - 2026-01-27
### Changed
- Remove use of bigcl when no clusters are computed (`901afd4`)
- Delete input.json (`256e090`)
- Fixed exiting when table potentials do not exist (`ac497b5`)
- Use >= minclsize in all code (`c727d8c`)
- Added radius of gyration (`cc9ac05`)
- V 0.2.10 Radius of gyration working (`babf0c4`)
- Radius of gyration distribution instead of Rcl (`60058fc`)
- Added scattering lengths to S(Q) calculation (`4791af5`)
- V 0.2.11 moved bsc to namelist INPUT_SQ (`a02df44`)
- First complete version. Added 3D cluster order parameters, minimum cluster size minclsize eliminated in favor or kmin (`d9eaf34`)
- Example input files modified (`c1a1c60`)

### Fixed
- Correct README, remove trash Makefiles (`7113984`)
- Correct printout of S(Q) when b.eq.1 (`b774c3c`)
- V 0.2.20 Corrected error in allocation of dynamics, implemented local number fluctuations (`ada5f22`)

### Performance
- Final fully parallel version implementing thermodynamics as well (`8869127`)

### Documentation
- Update README.md (`38ed406`)
- Updated README.md (`ecc3a34`)

## [0.2.7] - 2024-10-18
### Changed
- Thermodynamics with cells (`9ca9689`)
- First attemp GPU intracluster energy (`ee20e48`)
- Potgl cor (`9794685`)
- First partly working version (`0ed2e26`)
- V 0.2.9 almost ready (`3e067f7`)
- V corrected name of rdf_ij (`609646b`)
- Remove manged memory (`b4324ea`)
- Test problems in ladon (`45d6e21`)
- Output improvements, no thermodynamics computed when not all species selected. Program documentation improved (`72cb933`)

### Performance
- V 0.2.10 first fully parallel working version (2D to be tested) (`de0fc63`)

### Refactored
- Clean up and further comment code (`bec05d5`)

## [0.2.6] - 2024-10-16
### Added
- First attempt to implement cells in energy (`012255d`)

### Changed
- Potencl implemented (`8fa0896`)
- Potencl implemented, refactoring deltpots <-> potnbins (`4ab9188`)
- Changes in tmax (`6f27409`)
- Starting calculation of u using cells and GPU (`c73a708`)
- Corrected tmax (`37fd06d`)
- Corrected cellsp (`ccd9f9b`)
- Thermo_init corrected for reading potential tables (`d754bd4`)
- Working version after reorganization. Cluster potential energies working (`67194d7`)

## [0.2.5] - 2024-10-15
### Added
- Add tools (`77e64ef`)
- Add various controls for dynamics and input files (`f3b7ef1`)

### Changed
- Poteng2d developed, improving (`ba2ef21`)
- Potengcl implementing (`4e01f79`)
- Adding read_nc tool (`f03db42`)
- Time units (`6f728e0`)
- Remove initialization of profiles from sq_init (`dec0e9f`)
- Modified printout of large number of g(r) and ntypes/nsp (`1c70210`)
- Intracluster poteng testing (`870aaea`)
- Potmargin fixed! (`7ddcdfa`)
- Potenghistomixcl refactoring (`01a3c8f`)
- Major reorganization of main program (`e5b6390`)
- First proper compilation V0.2.5 (`091cf48`)

### Fixed
- Remove refrences to JSON and correct boundaries of profile (`35cdad8`)
- Correct for rdf to be calculated when density profiles are required (`6e3db50`)
- Poteng old gpus bug fixed (`d1452b0`)

## [0.2.4] - 2024-10-15
### Changed
- Testing poteng and various changes (`1f3bebe`)
- Poteng developed, improving (`ff87fed`)
- Removed log.cuf (`7e03294`)
- Revert potengcl core (`00c0cf0`)
- Prepare for potengcl (`d06bc05`)
- Implementing potengcl (`9859d5d`)
- Thermo_init_cluster implemented (`a1cec54`)
- Potengcl implemented (`b598737`)
- PrintPotEngCl implemented (`655e439`)
- Cells for potential calculation removed (`0bdb7cf`)

## [0.2.3] - 2024-10-15
### Changed
- Pot energy implemented, remains to be tested (`15c5e28`)
- First version with thermodynamics computed from lammps sored computes (`9bfaef8`)
- Changed input (`f9acebd`)
- Corrected epotperatom (`b745db4`)

### Documentation
- Corrected README.md (`c7b2072`)
- Update README.md (`6276270`)

## [0.2.2] - 2024-03-11
### Changed
- Thermo execution control (`f6e396d`)

### Fixed
- Correct x1 and x1 in log.f90 for 1 component (`20913fc`)
- Correct format in dynamics output (`ef05bd0`)

## [0.2.1] - 2024-02-14
### Changed
- Improve printout of dynamics and S(Q)'s (`da4bcea`)
- Some comments added (`e3ce6c4`)
- Corrected output of velocities in centers.lammpstrj (`9a1de63`)

### Documentation
- Update file README.md (`a8342f4`)
- Update README.md (`b59e1d9`)
- Update readme (`37302a3`)

## [0.2.0] - 2024-02-10
### Changed
- Dynamics finally working (`5f97d59`)
- Potential rc & tables ready in CPU & GPU memory (`2821e86`)
- Added Fourier transforms (`85913d4`)

## [0.1.0] - 2024-02-07
### Added
- Add velocities (`dce8894`)

### Changed
- Initial commit (`13d0e42`)
- Init repo (`a09a5d7`)
- First commit with dynamics (`4e8283f`)
- First trial of dynamics (`5f6714b`)
- First working version of dynamics (`dfdde25`)
- Added velocities (`3c9fdc2`)
- Added dyn_clear (`3c46308`)
- Working version with velocities ? (`33a97b4`)
- Corrected initialization (`9e5dbdb`)
- Still searching (`23c3237`)

### Fixed
- Correct vel assignments (`956cf81`)
- Seeking the error (`f0a8a62`)
- Correct transfers (`4359d01`)

### Refactored
- Reorganize clean up and initialization (`5907dab`)
