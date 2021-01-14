# chips_ozstar_tools
Repository for scripts to run CHIPS related tasks on ozstar

## run_CHIPS.py
This is a wrapper to submit jobs to run CHIPS with the correct order of dependencies of running the jobs. Requires that you put an environment variable `$CHIPS_ENV_VARIABLES` in your `~/.bashrc` to points to you CHIPS environment variables. An explicit example is adding this line:

```
export CHIPS_ENV_VARIABLES=/fred/oz048/MWA/CODE/CHIPS/JL_CHIPS_OUT/ozstar_env_variables.sh
```
to the `~/.bashrc`, where the contents of `/fred/oz048/MWA/CODE/CHIPS/JL_CHIPS_OUT/ozstar_env_variables.sh` are:

```
#!/bin/bash

CODEDIR=/fred/oz048/MWA/CODE/CHIPS/chips_2018/bin/
OUTPUTDIR=/fred/oz048/MWA/CODE/CHIPS/JL_CHIPS_OUT/
OBSDIR=/fred/oz048/MWA/CODE/CHIPS/chips/obsinfo_files/
INPUTDIR=/fred/oz048/MWA/data/
BEAMDIR=/fred/oz048/MWA/CODE/CHIPS/beam_files/
PBSDIR=/fred/oz048/MWA/CODE/CHIPS/chips/scripts/
PLOTSDIR=/fred/oz048/MWA/CODE/CHIPS/chips_plots/

export CODEDIR
export OUTPUTDIR
export FINALDIR
export INPUTDIR
export BEAMDIR
export OBSDIR
export PBSDIR
```
will guarantee the outputs are written to the correct `OUTPUTDIR` - if you haven't already, create your own output folder in `/fred/oz048/MWA/CODE/CHIPS/`, create your own `ozstar_env_variables.sh` pointed to the `OUTPUTDIR` and then point `$CHIPS_ENV_VARIABLES` to it inside you `~/.bashrc`.

`python run_CHIPS.py --help` details each arguement. As a bare minimum, you need to tell `run_CHIPS.py` where the uvfits data live, how many observations you want to run, and what band they are. There is a combination of 3 arguements: `--data_dir, --uvfits_dir, --uvfits_tag` which manifest into looking for data in `/data_dir/obsID/uvfits_dir/uvfits_tag%02d.uvfits`.  Note it always takes the obsID from the list in `--obs_list`. There are many defaults so you may not have to use every arguement. By default, if `run_CHIPS.py` cannot find the uvfits, it'll throw an error telling you where it looked for them and exit. You can switch this behaviour off with `--no_uvfits_check`.

An example command using the wrapper is

```sh
python /fred/oz048/MWA/CODE/chips_ozstar_tools/run_CHIPS.py \
   --obs_list=/fred/oz048/jline/test_woden/list_test.txt \
   --output_tag=new_CHIPS_wrapper \
   --band=high --obs_range=0,3 \
   --data_dir=/fred/oz048/MWA/data/2013/RTS_residuals \
   --no_clean --no_run
```
where `list_test.txt` contains obs ID numbers. `--no_clean` means the clean script will not be automatically run. `--no_run` also means the jobs will not be submitted to the queue, so you can check the scripts make sense. Remove this option and the jobs will be submitted to the queue automagically.



## plot_chipsout_general.pro  

This is an IDL script to perform the final binning, masking, and normalizations for CHIPS outputs given either FHD or RTS inputs. Diagnostic plots are also made.

Various libraries are required. Please see the the [dependencies](https://github.com/EoRImaging/eppsilon#dependencies) of eppsilon if a local installation is desired. Otherwise, add these paths to your `IDL_PATH` to your environment:
```
/fred/oz048/MWA/CODE/CHIPS/PS_libs/coyote
/fred/oz048/MWA/CODE/CHIPS/PS_libs/fhdps_utils
/fred/oz048/MWA/CODE/CHIPS/PS_libs/Healpix_3.11
/fred/oz048/MWA/CODE/CHIPS/PS_libs/pro/astro
/fred/oz048/MWA/CODE/CHIPS/PS
```

To run, open up an interactive IDL session   
```
idl
```

Then, compile the program if it not in your environment path or if you are not currently in its directory  
```
.compile <path>/<to>/<chips_ozstar_tools>/plot_chipsout_general.pro
```

You can now run the script with your specific keywords, for example   
```
plot_chipsout_general, 'my_output_tag', initials='NB', RTS=1, twoD=1, n_freq=384, band='high', lssa_num=0, beam_point_weight = [0.,0.,0.,0.,1.,0.,0.,0.,0.]
```

### keywords

General options:   
* output_tag: REQUIRED. Enter the output_tag that was used to run CHIPS to generate N-S and E-W power spectra. Alternatively, add a string array of two output_tags: `['xx_0.iter.output_tag1', 'xx_0.iter.output_tag2']` to specify pols directly or `['output_tag1', 'output_tag2']` which defaults to the `xx` pols only.
* initials: Enter in your initials that were used to create your personal CHIPS_OUT directory (Default: empty).   
* FHD, RTS: Specify FHD or RTS inputs with those keywords, set to 1 in case of True.   
* oneD, twoD: Specify 1D and 2D plotting mode, set to 1 in case of True.   
* n_freq: Enter the number of frequency channels, defaults to 384.   
* band, base_freq: Enter the band name (high, low, ultra-low, defaults to high) or the lowest frequency in the band in Hz. Setting band will default the base_freq to the expected value. Do not use if analyzing subbands.   
* chan_width: Enter the frequency resolution in Hz, defaults to 80e3.   
* lssa_num: Specify whether kriging has been performed on the input (0: kriging applied, 1: kriging not applied, defaults to 0).   
* band_point_weight: An array of 9 integers to indicate what pointings, from -4 to 4, were used. Default is just zenith. For example, if pointings -2 to 2 were used, then band_point_weight = [0,0,1,1,1,1,1,0,0].  
* obs_volume_file: Optional observational volume csv file. Defaults to corresponding obs volume file in observation_volumes
* set_tsys: Either set the tsys directly (i.e. set_tsys=200.) or set to 1 to scale the tsys to an updated theoretical number. For devel purposes and pipeline comparison only. Default: 0

1D cutting options:    
* wedge_cut: 0 = no cut, 3.0 = horizon. Default is 0.   
* wedge_angle: Alternatevely, provide sky angle in degrees to cut. Will overwrite wedge_cut keyword. 103.7 is the horizon, 120 is buffer in Beardsley et al. 2016 and Barry et al. 2019b 
* kperp_min_1D_lambda: Minimum k perpendicular in lambda. Default 10 wavelengths
* kperp_max_1D_lambda: Maximum k perpendicular in lambda. Default 50 wavelengths
* kperp_min_1D_Mpc: Minimum k perpendicular in Mpc^-1. Default is unset. Supercedes kperp_min_1D_lambda if set.
* kperp_max_1D_Mpc: Maximum k perpendicular in Mpc^-1. Default is unset. Supercedes kperp_max_1D_lambda if set.
* kpar_min_1D: Minimum k parallel in inverse Mpc. Default 0 Mpc^-1. Value of 0.105 Mpc^-1 used in Beardsley et al. 2016 and Barry et al. 2019b   
* kpar_max_1D: Maximum k parallel in inverse Mpc. Default is 10 Mpc^-1    

Plotting options:
* plot_max_1D: 1D plot maximum range. Default is 10^11 (units are either mK^2 (mK_units=1) or mK^2 h^-3 Mpc^3 (mK_units=0))
* plot_min_1D: 1D plot minimum range. Default is 10^3 (units are either mK^2 (mK_units=1) or mK^2 h^-3 Mpc^3 (mK_units=0))
* plot_max_2D: 2D plot maximum range. Default is 10^15 mK^2 h^-3 Mpc^3
* plot_min_2D: 2D plot minimum range. Default is 10^3 mK^2 h^-3 Mpc^3
* mK_units: 1D plots in mK^2 units. Default is 1.

Input and output options:    
* output_dir: Output directory for plots. Defaults to your home directory.
* input_dir: Input directory. Default is built off of initials and lssa_num inputs  
* pdf: Set to one to make pdfs instead of pngs   

If anything else is required, please inquire, and we will add an interface.
