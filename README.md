# PnetCDF I/O Benchmark using GCRM Application I/O Kernel

This software benchmarks the performance of
[PnetCDF](https://parallel-netcdf.github.io/) method implementing the I/O
kernel of Global Cloud Resolving Model (GCRM) simulation codes. GCRM was
developed at the Colorado State University
(http://kiwi.atmos.colostate.edu/gcrm). GCRM's I/O module, called GIO library,
was developed at Pacific Northwest National Laboratory
(https://svn.pnl.gov/gcrm). GCRM and GIO are written in Fortran90.

This software package extracts the I/O kernel from GCRM and converts it into C
for more flexible parameter settings (e.g. number of MPI processes, max level)
and dynamic memory allocation for buffer management.

## Note on memory requirement
GCRM is a memory-intensive parallel application which allocates a significant
amount of memory per process when running. This software preserves such
characteristics to reflect its high demands on memory space.  For detailed
information about the variable sizes and the amount of memory required, please
refer to documents in directory `./doc`.

## Software requirements
 1. MPI compiler
 2. PnetCDF version 1.4.0 or later

## Build instructions
Run commands below.
```
./configure --with-pnetcdf=/path/to/PnetCDF MPICC=your_MPI_C_compiler
make
```
It creates an executable file named `gcrm_io` under directory `./run`.

## Run Instructions
Several input parameter files are available under directory `./inputs` which
should be copied over to directory `.run`. The parameter file can be customized
to run different problem sizes.

The executable `gcrm_io` takes one command-line argument which is the input
parameter file name that contains all the parameters settings. If no
command-line argument is given, it reads the file named `zgrd.in`. The
following parameter files are provided, each specifying a particular problem
size.
```
zgrd.in.r6     - refine level ==  6, corresponding to grid size 112    Km
zgrd.in.r7     - refine level ==  7, corresponding to grid size  55.9  Km
zgrd.in.r8     - refine level ==  8, corresponding to grid size  27.9  Km
zgrd.in.r9     - refine level ==  9, corresponding to grid size  14.0  Km
zgrd.in.r10    - refine level == 10, corresponding to grid size   6.98 Km
zgrd.in.r11    - refine level == 11, corresponding to grid size   3.49 Km
```

### Parameters determine the number of solution variables
Parameter `io_desc_file` points to file `ZGrd.desc`, which describes all GCRM
variables, including all grid and field variables with information about their
size, type, and dimension.

Parameter `io_config_file` points to a file describing the I/O configuration.

Parameter `physics_mode` indicates if additional physics variables will be
added to the benchmark. Enabling this mode will demand more memory space on
each MPI process.

Four pre-defined configurations are provided in `./inputs`.
```
ZGrd.LZ.fcfg.nophys      writes 11 solution variables, one variable per file,
                         no physics variables. Use it for testing only.
ZGrd.CP.lyr.fcfg         writes 38 solution variables, one per netCDF file.
                         It can only be used when physics_mode is enabled.
ZGrd.CP.lyr.fcfg.groups  writes 38 solution variables in 12 netCDF files,
                         some variables that are often accessed together are
                         saved in the same file. It can only be used when
                         physics_mode is enabled.
ZGrd.CP.lyr.fcfg.one*    writes 38 solution variables in one netCDF file, all
                         variables are saved in the same file. It can only be
                         used when physics_mode is enabled.
```

### Parameter for setting refine level
```
level_max
```
### Parameter for setting frequency of writes
The following 4 parameters determine/mimic the length of simulation:
 * end_days
 * end_hours
 * end_minutes
 * end_seconds

Together with parameter `cdf_output_frequency`, the number of writes is
determined. For example, if the length of simulation is set to 1 minutes
(`end_minutes` == 1) and `cdf_output_frequency` is set to 60 (seconds), then
there will be 2 writes (including one at the initial dump).

Please note that `cdf_output_frequency` is the amount of time advanced in
each simulation iteration. Depending on the value of `level_max`, set
`cdf_output_frequency` correspondingly. Below is the possible value pairs.
```
for level_max        =   3,   4,  5,  6,  7, 8, 9, 10, 11, 12
cdf_output_frequency = 900, 120, 60, 30, 15, 8, 4,  2,  1,  1
```

### Parameters specify the directories to store the output files
Users should edit the input parameter file and change the following two
parameters `output_path` and `cdf_output_path` to the locations where the
output files are to be saved.
```
parameter "output_path"     points to the path for miscellaneous outputs
parameter "cdf_output_path" points to the path for netcdf output files
```

### Note on restart options
Although restart options are provided (see the bottom of file `zgrd.in`),
not all of them are supported in this release, in particular, reading part
is not completed. The only options taking effect are
 * restart_interval
 * restart_output
 * restart_overwrite

### I/O methods
Currently there are four I/O options (parameter `iotype`):
 * nonblocking_collective
 * nonblocking_independent
 * blocking_collective
 * blocking_independent

which correspond to using PnetCDF nonblocking APIs with collective flush,
nonblocking APIs with independent flush, blocking collective APIs, and blocking
independent APIs, respectively. The original GCRM has two additional options
namely direct and interleaved, but they are not implemented in this benchmark.

### Number of processes to run and command-line arguments
GCRM runs only on certain numbers of MPI processes, depending on the value of
sbdmn_iota. The number of processes must be able to divide evenly
```
M = 10 x 4 ^ sbdmn_iota
```
`M` is the total number of subdomains. This is to ensure the numbers of
subdomains assigned to individual MPI processes are the same.

`sbdmn_iota` is the first index i of array `nblocks` below such that
`nblocks[i]` is divisible by the number of MPI processes.
```
int nblocks[10] = {10, 40, 160, 640, 2560, 10240, 40960, 163840, 655360, 2621440};
```
For example, if the number of MPI processes, `nprocs`, is 80, then `sbdmn_iota`
will be set to 2, because `nblocks[2]` will be picked 160 which is the minimum
number in `nblocks[]` divisible by `nprocs`.

Note that because GCRM is memory intensive, the number of processes to run for
a refine level, set in `level_max`, must be sufficiently large to avoid running
out of memory. In `doc/data\ summary.pdf`, the 3rd page shows a table of the
memory usage for various combinations of number of processes and `level_max`
when writing all physics variables is enabled and no restart write is
performed.  The memory usage will be reported for different settings.

Again, the default input parameter file is `zgrd.in` (when no command-line
argument is given.) Users can specify a different parameters file, such as
```
mpiexec -n 8 ./gcrm_io zgrd.in.r6
```

## Example output from stdout
Below is the output on screen of an example run on 20 MPI processes using the
input files `ZGrd.CP.lyr.fcfg`, `ZGrd.desc`, and `zgrd.in` copied from folder
`inputs`.
```
% mpiexec -n 20 ./gcrm_io

while loop num_dumps=1 time_ZGrd=0.000000 gio_default_frequency=60
while loop num_dumps=2 time_ZGrd=60.000000 gio_default_frequency=60
GIO: ----------------------------------------------------------
GIO statistics: (cumulative across all MPI processes)
GIO: total_time_in_API                           1566.9720
GIO: time_in_nf_put_var                             0.0000
GIO: time_in_nf_put_var_grid                        6.6997
GIO: time_in_nf_put_att                             0.0336
GIO: time_in_nf_def_dim                             0.0063
GIO: time_in_nf_def_var                             0.0059
GIO: time_in_update_time                            0.0036
GIO: time_in_nf_create                            226.8757
GIO: time_in_nf_open                               30.7878
GIO: time_in_nf_close                               5.6107
GIO: time_in_nf_enddef                              2.1618
GIO: time_in_nf_inq_varid                           0.0612
GIO: time_in_nf_inq_dimlen                          0.0000
GIO: time_in_nf_iput                                4.8521
GIO: time_in_nf_iget                                0.0000
GIO: time_in_nf_wait                             1246.1515
GIO: time_in_API_copy                               1.3175
GIO: time_in_avgs                                   0.1935
GIO: ----------------------------------------------------------
GIO: bytes_API_write (Bytes):                   3654442380
GIO: bytes_API_write (MiB):                      3485.1478
GIO: bytes_API_write (GiB):                         3.4035
GIO: bandwidth for writes (MiB/sec):               44.4826
GIO: bandwidth for writes (GiB/sec):                0.0434
GIO: ----------------------------------------------------------
---- MPI file info used ----
MPI File Info: [ 0] key =     romio_pvfs2_debugmask, value = 0
MPI File Info: [ 1] key =           striping_factor, value = 0
MPI File Info: [ 2] key =             striping_unit, value = 0
MPI File Info: [ 3] key =    romio_pvfs2_posix_read, value = disable
MPI File Info: [ 4] key =   romio_pvfs2_posix_write, value = disable
MPI File Info: [ 5] key =    romio_pvfs2_dtype_read, value = disable
MPI File Info: [ 6] key =   romio_pvfs2_dtype_write, value = disable
MPI File Info: [ 7] key =   romio_pvfs2_listio_read, value = disable
MPI File Info: [ 8] key =  romio_pvfs2_listio_write, value = disable
MPI File Info: [ 9] key =            cb_buffer_size, value = 16777216
MPI File Info: [10] key =             romio_cb_read, value = enable
MPI File Info: [11] key =            romio_cb_write, value = enable
MPI File Info: [12] key =                  cb_nodes, value = 8
MPI File Info: [13] key =         romio_no_indep_rw, value = true
MPI File Info: [14] key =              romio_cb_pfr, value = disable
MPI File Info: [15] key =         romio_cb_fr_types, value = aar
MPI File Info: [16] key =     romio_cb_fr_alignment, value = 1
MPI File Info: [17] key =     romio_cb_ds_threshold, value = 0
MPI File Info: [18] key =         romio_cb_alltoall, value = automatic
MPI File Info: [19] key =        ind_rd_buffer_size, value = 4194304
MPI File Info: [20] key =             romio_ds_read, value = automatic
MPI File Info: [21] key =            romio_ds_write, value = disable
MPI File Info: [22] key =            cb_config_list, value = *:1
MPI File Info: [23] key =      nc_header_align_size, value = 2048
MPI File Info: [24] key =         nc_var_align_size, value = 1
MPI File Info: [25] key = nc_header_read_chunk_size, value = 0
====  gcrm-io-pnetcdf 1.0.0   released on  22 Nov 2013 ====
----  Run-time parameters ---------------------------------
----  Number of processes                            = 20
----  level_max  (global horizontal grid resolution) = 5
----  sbdmn_iota (see grid_params.h)                 = 1
----  level_glbl (see grid_params.h)                 = 2
----  km         (number of vertical layers)         = 256
----  cell_max   (global number of cells)            = 10242
----  im         (local  number of cells along i)    = 18
----  jm         (local  number of cells along j)    = 18
----  nsdm_glbl  (global number of blocks)           = 40
----  nsdm       (local  number of blocks)           = 2
----  using physics variables                        = enabled
----  number of files written                        = 40
----  number of grid  variables written              = 16
----  number of field variables written              = 74
----  number of snapshot dumps                       = 2
---------------------------------------------------------------
Timing results (max among all processes)
init     time=             1.43 sec
Max comp time=             2.23 sec
Max I/O  time=            77.58 sec
finalize time=             0.06 sec
---------------------------------------------------------------
Write  amount=          3654.51 MB      =        3485.22 MiB
Read   amount=             9.08 MB      =           8.66 MiB
I/O    amount=          3663.60 MB      =        3493.88 MiB
             =             3.66 GB      =           3.41 GiB
I/O bandwidth=            47.22 MB/sec  =          45.04 MiB/sec
             =             0.05 GB/sec  =           0.04 GiB/sec
---------------------------------------------------------------
memory MAX usage (among 20 procs) =        315.78 MiB
memory MIN usage (among 20 procs) =        315.48 MiB
memory AVG usage (among 20 procs) =        315.62 MiB
```

## Output files created in the output folder.
```
% ls -lhgG
total 3.5G
-rw------- 1  21M Jul 18 15:24 cloud_ice_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 cloud_water_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 divergence_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 exner_lyr_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 geopotential_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 graupel_mmr_19010101_000000.nc
-rw------- 1 2.0M Jul 18 15:23 grid.nc
-rw------- 1  21M Jul 18 15:24 heat_flux_vdiff_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 heating_latent_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 heating_lw_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 heating_lw_cs_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 heating_sw_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 heating_sw_cs_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 ke_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 mass_19010101_000000.nc
-rw------- 1  83K Jul 18 15:24 olr_19010101_000000.nc
-rw------- 1  83K Jul 18 15:24 prec_frz_19010101_000000.nc
-rw------- 1  83K Jul 18 15:24 prec_tot_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 pressure_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 qci_tend_micro_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 qcw_tend_micro_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 qgr_tend_micro_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 qrw_tend_micro_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 qsn_tend_micro_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 qwv_tend_micro_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 rain_mmr_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 rel_vorticity_19010101_000000.nc
-rw------- 1 1.4G Jul 18 15:24 restart.nc19010101_000000.nc
-rw------- 1 1.4G Jul 18 15:25 restart.nc19010101_000100.nc
-rw------- 1  21M Jul 18 15:24 snow_mmr_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 strm_func_19010101_000000.nc
-rw------- 1  83K Jul 18 15:24 swinc_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 temperature_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 vel_pot_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 vorticity_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 water_vapor_19010101_000000.nc
-rw------- 1  61M Jul 18 15:24 wind_19010101_000000.nc
-rw------- 1  41M Jul 18 15:24 wind_crn_ew_19010101_000000.nc
-rw------- 1  41M Jul 18 15:24 wind_crn_ns_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 wtr_flux_vdiff_19010101_000000.nc
-rw------- 1  21M Jul 18 15:24 w_vert_19010101_000000.nc
```

## Questions/Comments:
email: wkliao@eecs.northwestern.edu

Copyright (C) 2013, Northwestern University

See COPYRIGHT notice in top-level directory.

