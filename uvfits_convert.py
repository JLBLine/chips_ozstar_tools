import numpy
import os
import copy
import subprocess
from astropy.io import fits
import argparse



#********************************
def main():

    # Parse the command line inputs. 
    parser=argparse.ArgumentParser(description=
        "Convert outputs from RTS (per coarse band uvfits) and FHD (visibility IDL save files) into CHIPS-compatible " +
        "(per coarse band uvfits) or FHD-image->eppsilon-compatible (calibrated uvfits over full band) outputs. " +
        "Example call to turn FHD output into CHIPS uvfits: "
        "python uvfits_convert.py --chips_format --specific_dir=/path/to/dir/fhd_my_test_name --obsfile_name=~/obs.txt." +
        "Example call to turn RTS outputs into FHD-image->eppsilon uvfits: " +
        "python uvfits_convert.py --fhd_epp_format --specific_dir=/path/to/RTS/dir --sub_dir=my_time_stamp --obsfile_name=~/obs.txt")    
    parser.add_argument("-c", "--chips_format", default=None, action='store_true', 
        help="Convert to CHIPS format: uvfits per coarse band")
    parser.add_argument("-e", "--fhd_epp_format", default=None, action='store_true', 
        help="Convert to FHD/eppsilon format: one uvfits")
    parser.add_argument("-s", "--specific_dir", default=None, required=True,
        help="Full path to the FHD or RTS common directory. For example, an FHD common directory would be " +
        "/path/to/dir/fhd_my_test_name/ and an RTS common directory would be /path/to/data/dir/. " + 
        "An additional subdirectory for nested RTS analysis runs can be specified with --sub_dir")
    parser.add_argument("--sub_dir", default=None, 
        help="Subdirectory to specify a specific analysis run within the common directory, common with RTS " +
        "processing.")
    parser.add_argument("-o", "--obsfile_name", required=True,
        help="Name of the file containing the observation list")
    parser.add_argument("--obsinfo_dir", default=None,
        help="Specify path to the common directory where CHIPS observation info files are kept to create obsinfo " +
        "files if they do not already exist")
    parser.add_argument("--coarse_bands", default=24, type=int,
        help="Specify number of coarse bands in the observations. Default is 24.")
    parser.add_argument("-d", "--make_dirty", default=None, action='store_true',
        help="Make dirty uvfits.")
    parser.add_argument("-m", "--make_model", default=None, action='store_true',
        help="Make model uvfits")
    parser.add_argument("-r", "--make_residual", default=True, action='store_true',
        help="Make residual uvfits. Default.")
    parser.add_argument("--dry_run", default=None, action='store_true',
        help="Dry run. Run the code without submitting anything.")
    args = parser.parse_args()

    if not args.chips_format and not args.fhd_epp_format:
        exit("Please specify either --chips_format or --fhd_epp_format")
        
    if args.chips_format:
        out = fhd_vis_to_chips_format(specific_dir=args.specific_dir, obsfile_name=args.obsfile_name, 
                                      obsinfo_dir=args.obsinfo_dir, coarse_bands=args.coarse_bands, 
                                      make_dirty=args.make_dirty, make_model=args.make_model, 
                                      make_residual=args.make_residual, dry_run=args.dry_run)
    
    if args.fhd_epp_format:
        out = rts_vis_to_fhd_format(specific_dir=args.specific_dir, obsfile_name=args.obsfile_name,
                                    sub_dir=args.sub_dir, coarse_bands=args.coarse_bands, 
                                    dry_run=args.dry_run)

#********************************

#********************************
# Turn FHD sav directories into CHIPS compatible uvfits per coarse band
# Example call:
# python uvfits_convert.py --chips_format --specific_dir=/path/to/dir/fhd_my_test_name --obsfile_name=~/obs.txt
def fhd_vis_to_chips_format(specific_dir=None, obsfile_name=None, obsinfo_dir=None, coarse_bands=None,
                            make_dirty=None, make_model=None, make_residual=None, dry_run=None):
    from pyuvdata import UVData

    make_dirty_init=make_dirty
    make_model_init=make_model
    make_residual_init=make_residual

    # Get obsids to download
    obsfile = open(obsfile_name, "r")
    obsids = [line.split( ) for line in obsfile.readlines()]
    obsids = [obs[0] for obs in obsids]
    obsfile.close()

    print('Output directory for uvfits is %s/CHIPS_input' %(specific_dir))
    if not os.path.exists(specific_dir + '/CHIPS_input'):
           os.mkdir(specific_dir + '/CHIPS_input')

    for obs_id in obsids:
        # Initialize
        model_exists = 0
        dirty_exists = 0
        res_exists = 0
        make_dirty=make_dirty_init
        make_model=make_model_init
        make_residual=make_residual_init
        script_flag=0

        # Make the dir in the FHD output folder
        if not os.path.exists(specific_dir + '/CHIPS_input/' + obs_id):
            os.mkdir(specific_dir + '/CHIPS_input/' + obs_id)

        # Set CHIPS shell script flag to make obs info if not present already
        if obsinfo_dir:
            if (not os.path.exists(obsinfo_dir + 'obsid' + obs_id + '.txt')) or \
               (not os.path.exists(obsinfo_dir + 'obsidra' + obs_id + '.txt')) or \
               (not os.path.exists(obsinfo_dir + 'obsidband' + obs_id + '.txt')):
                script_flag = 1

        # Check for CHIPS compatible uvfits for specified obs, skip if they already exist
        if (os.path.exists(specific_dir + '/CHIPS_input/' + obs_id + '/uv_model_' + str(1).zfill(2) + '.uvfits')) and \
           (os.path.exists(specific_dir + '/CHIPS_input/' + obs_id + '/uv_model_' + str(coarse_bands).zfill(2) + '.uvfits')):
            model_exists = 1

        if (os.path.exists(specific_dir + '/CHIPS_input/' + obs_id + '/uv_dirty_' + str(1).zfill(2) + '.uvfits')) and \
           (os.path.exists(specific_dir + '/CHIPS_input/' + obs_id + '/uv_dirty_' + str(coarse_bands).zfill(2) + '.uvfits')):
            dirty_exists = 1

        if (os.path.exists(specific_dir + '/CHIPS_input/' + obs_id + '/uv_res_' + str(1).zfill(2) + '.uvfits')) and \
           (os.path.exists(specific_dir + '/CHIPS_input/' + obs_id + '/uv_res_' + str(coarse_bands).zfill(2) + '.uvfits')):
            res_exists = 1

        if make_residual and res_exists:
            make_residual = 0
            print("Residual already exists for " + obs_id)

        if make_dirty and dirty_exists:
            make_dirty = 0
            print("Dirty already exists for " + obs_id)

        if make_model and model_exists:
            make_model = 0
            print("Model already exists for " + obs_id)

        # Construct the list of files
        fhd_prefix_vis = specific_dir + '/vis_data/' + obs_id + '_'
        fhd_prefix_set = specific_dir + '/metadata/' + obs_id + '_'

        if not(os.path.exists(fhd_prefix_vis + 'vis_XX.sav')) or not(os.path.exists(fhd_prefix_vis + 'vis_YY.sav')):
            print(obs_id + ' FHD visibilities are missing')
            continue

        fhd_files = [fhd_prefix_vis + 'flags.sav', fhd_prefix_vis + 'vis_XX.sav', fhd_prefix_set + 'params.sav',
                    fhd_prefix_vis + 'vis_YY.sav', fhd_prefix_vis + 'vis_model_XX.sav',
                    fhd_prefix_vis + 'vis_model_YY.sav', fhd_prefix_set + 'settings.txt',
                    fhd_prefix_set + 'layout.sav']
        if make_model and not(make_dirty) and not(make_residual):
            fhd_files = [fhd_prefix_vis + 'flags.sav', fhd_prefix_set + 'params.sav',
                    fhd_prefix_vis + 'vis_model_XX.sav', fhd_prefix_vis + 'vis_model_YY.sav', 
                    fhd_prefix_set + 'settings.txt', fhd_prefix_set + 'layout.sav']

        if make_dirty:
            print("Creating dirty vis for " + obs_id)
        if make_model:
            print("Creating model vis for " + obs_id)
        if make_residual:
            print("Creating residual vis for " + obs_id)

        if dry_run:
            continue

        # Read in all fhd files into UV class
        if make_dirty or make_residual:
            UV_dirty = UVData()
            UV_dirty.read_fhd(fhd_files,run_check=False,check_extra=False,run_check_acceptability=False)
        if make_model or make_residual:
            UV_model = UVData()
            UV_model.read_fhd(fhd_files,use_model=1)

        # Find the field center from the settings file in order to rephase later
        # Find the number of frequencies and the frequency resolution for parsing into per coarse channel
        # Settings file may contain multiple instances of string, so grab only the first occurance 
        with open(fhd_prefix_set + 'settings.txt') as settings_file:
            for line in settings_file:
                if ("ORIG_PHASERA" in line) and ('field_center_ra' not in locals()):
                    field_center_ra = float(line.split()[1])
                if ("ORIG_PHASEDEC" in line) and ('field_center_dec' not in locals()):
                    field_center_dec = float(line.split()[1])
                if ("N_FREQ" in line) and ('n_freq' not in locals()):
                    n_freq = int(line.split()[1])

        # Find the upper and lower bound for each coarse band
        num_per_coarse = n_freq / coarse_bands
        mod_list = numpy.array([(x % num_per_coarse) for x in range(n_freq)])
        top_coarse_ind = numpy.where(mod_list == (num_per_coarse - 1)) 
        bottom_coarse_ind = numpy.where(mod_list == 0) 

        # Split up UV class into coarse bands and save uvfits
        # Visibilities need to be unphased from observation center and rephased to field center for RTS
        for coarse_i in range(len(bottom_coarse_ind[0])):
            if make_dirty:
                UV1 = copy.deepcopy(UV_dirty)
                UV1.unphase_to_drift()
                UV1.phase(field_center_ra,field_center_dec*numpy.pi/180,epoch="J2000")
                UV1.select(frequencies=UV_dirty.freq_array[0][bottom_coarse_ind[0][coarse_i]:top_coarse_ind[0][coarse_i]+1])
                UV1.write_uvfits(specific_dir + '/CHIPS_input/' + obs_id + '/uv_dirty_'+str(coarse_i+1).zfill(2)+'.uvfits', spoof_nonessential=True)
                out = update_header(specific_dir=specific_dir, obs_id=obs_id, coarse_i=course_i, field_center_ra=field_center_ra, 
                                    field_center_dec=field_center_dec, type_str='dirty')
            if make_model:
                UV1 = copy.deepcopy(UV_model)
                UV1.unphase_to_drift()
                UV1.phase(field_center_ra,field_center_dec*numpy.pi/180,epoch="J2000")
                UV1.select(frequencies=UV_model.freq_array[0][bottom_coarse_ind[0][coarse_i]:top_coarse_ind[0][coarse_i]+1])
                UV1.write_uvfits(specific_dir + '/CHIPS_input/' + obs_id + '/uv_model_'+str(coarse_i+1).zfill(2)+'.uvfits', spoof_nonessential=True)
                out = update_header(specific_dir=specific_dir, obs_id=obs_id, coarse_i=coarse_i, field_center_ra=field_center_ra, 
                                    field_center_dec=field_center_dec, type_str='model')
            if make_residual:
                res_data = UV_dirty.data_array - UV_model.data_array
                UV1 = copy.deepcopy(UV_dirty)
                UV1.data_array = res_data
                UV1.unphase_to_drift()
                UV1.phase(field_center_ra,field_center_dec*numpy.pi/180,epoch="J2000")
                UV1.select(frequencies=UV_dirty.freq_array[0][bottom_coarse_ind[0][coarse_i]:top_coarse_ind[0][coarse_i]+1])
                UV1.write_uvfits(specific_dir + '/CHIPS_input/' + obs_id + '/uv_res_'+str(coarse_i+1).zfill(2)+'.uvfits', spoof_nonessential=True)
                out = update_header(specific_dir=specific_dir, obs_id=obs_id, coarse_i=coarse_i, field_center_ra=field_center_ra, 
                                    field_center_dec=field_center_dec, type_str='res')

        # Run obs info generator for CHIPS input if it doesn't already exist
        # Create temporary text file of the obsid for obs_info input
        if script_flag:
            obsfile_new = open(specific_dir + '/CHIPS_input/' + obs_id + '.txt','w')
            obsfile_new.write(obs_id)
            obsfile_new.close()
            print('Making observation info for ' + obs_id)
            # OZSTAR-SPECIFIC FILE PATH
            subprocess.call(['/fred/oz048/MWA/CODE/CHIPS/chips/scripts/make_obs_info_ozstar.sh',specific_dir + '/CHIPS_input/' + obs_id + '.txt'])
            os.remove(specific_dir + '/CHIPS_input/' + obs_id + '.txt')

#********************************

#********************************
# Updates the uvfits header during generation to be CHIPS compatible
def update_header(specific_dir=None, obs_id=None, coarse_i=None, field_center_ra=None, field_center_dec=None, type_str=None):

    with fits.open(specific_dir + '/CHIPS_input/' + obs_id + '/' + "uv_" + type_str + "_" + 
                   "%02d" % (coarse_i+1,) +".uvfits", mode='update') as hdu:
        header = hdu[0].header
        crval7 = header['CRVAL7'] #dec
        crval7 = field_center_dec
        crval6 = header['CRVAL6'] #ra
        crval6 = field_center_ra
        ctype5 = header['CTYPE5'] #if
        ctype7 = header['CTYPE7']
        ctype6 = header['CTYPE6']
        crval5 = header['CRVAL5']
        cdelt5 = header['CDELT5']
        header['CTYPE7']=ctype5
        header['CTYPE5']=ctype6
        header['CTYPE6']=ctype7
        header['CRVAL7']=crval5
        header['CRVAL5']=crval6
        header['CRVAL6']=crval7
        header.set('OBSRA',crval6)
        header.set('OBSDEC',crval7)
        hdu[0].header = header
        hdu.flush()

#********************************

#********************************
# Turn RTS uvfits per coarse band outputs into FHD compatible full uvfits
# Example call:
# python uvfits_covert.py --fhd_epp_format --specific_dir=/path/to/RTS/dir --sub_dir=my_time_stamp --obsfile_name=~/obs.txt
def rts_vis_to_fhd_format(specific_dir=None, obsfile_name=None, sub_dir=None, coarse_bands=None, dry_run=None):
    from pyuvdata import UVData

    # Get obsids to download
    obsfile = open(obsfile_name, "r")
    obsids = [line.split( ) for line in obsfile.readlines()]
    obsids = [obs[0] for obs in obsids]
    obsfile.close()

    # Remove uvdump files afterwards
    #cleanup = 1

    for obs_id in obsids:
        if not os.path.exists(specific_dir + '/' + obs_id + '/' + sub_dir):
            print('Directory ' + specific_dir + '/' + obs_id + '/' + sub_dir + ' does not exist')
            print('Download/create uvdump files for ' + str(obs_id) + ' first')
            continue
        else:
            dir_obs = specific_dir + '/' + obs_id + '/' + sub_dir + '/'

        if os.path.exists(dir_obs + obs_id + '.uvfits'):
            print('uvfits already exists for ' + obs_id)
            continue

        if not os.path.exists(dir_obs + 'uvdump_01.uvfits'):
            print('Download/create uvdump files for ' + str(obs_id) + ' first using uvdump naming convention.')
            continue

        if dry_run:
            continue

        UV_total = UVData()

        for coarse_i in range(coarse_bands):
            UV = UVData()
            UV.read(dir_obs + 'uvdump_' + str(coarse_i+1).zfill(2)+'.uvfits')
            if coarse_i == 0:
                UV_total = UV
            else:
                UV_total = UV_total + UV

        # Get pointing centre ra/dec from metafits file, assumed to be in same location as uvdump files
        with fits.open(dir_obs + str(obs_id) + '.metafits') as hdu:
            header = hdu[0].header
            ra_obs = header['RA']
            dec_obs = header['DEC']
            ra_obs = ra_obs * numpy.pi/180 #pyuvdata expects radians
            dec_obs = dec_obs * numpy.pi/180 #pyuvdata expects radians

        # FHD expects uvfits to be phased to the observation ra/dec.
        # RTS uvfits are phased to the field centre by default.
        UV_total.unphase_to_drift()
        UV_total.phase(ra_obs,dec_obs,epoch="J2000")

        print("Writing " + dir_obs + obs_id + ".uvfits")
        UV_total.write_uvfits(dir_obs + obs_id + '.uvfits', spoof_nonessential=True)

        #if cleanup:
        #    files = os.listdir(dir_obs)
        #    for file_i in files:
        #        if file_i.startswith("uvdump"):
        #            os.remove(dir_obs + file_i)

#********************************

#********************************

if __name__ == '__main__':
    main()
