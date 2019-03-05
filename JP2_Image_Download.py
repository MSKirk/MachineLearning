import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from sunpy.net import hek, helioviewer
from sunpy.time import parse_time
from sunpy.coordinates import frames
from sunpy.io import read_file_header
from sunpy.map import Map

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import numpy as np

from imageio import imwrite

import datetime
from datetime import timedelta
import os, glob, shutil
from mahotas.polygon import fill_polygon
import csv




class Jp2ImageDownload:

    def __init__(self, save_dir, full_image_set=False, dt=timedelta(minutes=30),
                 tstart='2012/05/30 23:59:59', tend='2012/05/31 23:59:59'):
        """
        Set up initial start and ending times and save directory.

        :param save_dir: Directory to contain the full subdirectory tree of images and masks
        :param full_image_set: boolean if true will run through the entire SDO catalogue
        :param tstart: starting time string of any format read by parse_time
        :param tend: ending time string of any format read by parse_time
        :param dt: limit the time gap between downloaded image and requested hek time. Set to None for no limit.
        """

        # if full_image_set:
        #     tstart = '2010/05/31 23:59:59'
        #     tend = '2018/12/31 23:59:59'

        self.save_dir = os.path.abspath(save_dir)

        # Subdirectory for discarding the files from incomplete sets
        self.reject_dir = os.path.join(self.save_dir, 'incomplete_set')
        # Create the "rejection" subdirectory to move the images from incomplete groups
        os.makedirs(self.reject_dir, exist_ok=True)
        # Directory to save the label masks
        self.label_save_dir = os.path.join(self.save_dir, 'label_masks')
        os.makedirs(self.label_save_dir, exist_ok=True)

        self.tstart = parse_time(tstart)
        self.tend = parse_time(tend)
        self.hek_times = []

        self.date_list, self.date_string_list  = self.gen_date_list()
        # time gap tolerated between requested time and image time s.t. time_in - dt < actual image time < time_in + dt
        self.dt = dt
        # Build lists of measurements that will be requested from the SDO data archive
        self.aia_wav = ['1600', '1700', '94', '131', '171', '193', '211', '304', '335']
        self.hmi_segments = ['continuum', 'magnetogram']
        # Actual list of requested measurement
        self.measurements_req = self.aia_wav + self.hmi_segments
        # Build list of strings inferred from the jp2000 filenames showing what measurement they come from
        self.inst_file_map = ['SDO_AIA_AIA_{:s}'.format(wav) for wav in self.aia_wav]
        self.inst_file_map.append('SDO_HMI_HMI_continuum')
        self.inst_file_map.append('SDO_HMI_HMI_magnetogram')
        # store the list of files in the save directory
        self.jp2f = []
        self.mask_hek_time_map_csv = os.path.join(self.save_dir, 'mask_hek_time_map.csv')
        self.hek_time_jp2_map_csv = os.path.join(self.save_dir, 'hek_time_jp2_map.csv')
        self.rejected_hek_csv = os.path.join(self.save_dir, 'rejected_hek.csv')
        self.missed_downloads_csv = os.path.join(self.save_dir, 'missed_downloads.csv')
        self.rejected_hek_events = []
        self.missed_downloads = []
        self.download_flag = True
        self.do_plot = True


    def gen_date_list(self):
        """
        Create list of datetimes by day

        :return: list of datetimes broken into day segments and their string representation
        """
        time_between = self.tend - self.tstart
        ndays = max(1, time_between.days)
        date_list = [self.tstart + timedelta(days=dd) for dd in range(0, ndays+1)]
        date_string_list = [tt.strftime('%Y/%m/%d %H:%M:%S') for tt in date_list]
        return date_list, date_string_list


    def list_missing_measurements(self, hek_time, jp2_files, jp2_datetimes):
        """
        In case of an incomplete or earlier run that got interrupted, check what data are missing using the
        times and instrument descriptors parsed from the jpeg2000 file names.
        If no jpeg2000 files are present, returns the entire list measurements needed.

        :param hek_time: datetime time of the hek event
        :param jp2_files: list of jpeg2000 file paths already downloaded
        :param jp2_datetimes: list of datetime times of the jpeg2000 files
        :return: list of missing measurements that needs to be downloaded if available.
        """

        if not jp2_files:
            return self.measurements_req
        # Get the index of files whose time fall within a 2hr window centered on the hek event time
        tmatches = [idx for idx, jp2_time in enumerate(jp2_datetimes) if hek_time - self.dt < jp2_time < hek_time + self.dt]
        # Get the series of files at these indices and get which instruments & measurements they are
        inst_names = [os.path.basename(jp2_files[idx])[25:-4] for idx in tmatches]
        # Check what measurements are missing. Return missing ones in a list.
        measurements = [self.measurements_req[i] for i, x in enumerate(self.inst_file_map) if x not in inst_names]
        # if that measurements_ list is not empty, we need to download those.
        return measurements


    def download_images(self):
        """
        Run through the complete set of dates as determined in the self.__init__ call by day
        Sets up a subdirectory tree to organize the images

        """
        # IF any, move previously discarded files to the save directory (this to be consistent with the cleanup pass)
        discared_files = sorted(glob.glob(os.path.join(self.reject_dir, '*.jp2')))
        for f in discared_files:
            print('restoring discarded file: ' + os.path.basename(f))
            shutil.move(f, self.save_dir)

        # for ii, download_date in enumerate(self.date_list):
        # First download the images. Feature extraction will be done separately

        # Get a list of existing files (if any)
        jp2f = sorted(glob.glob(os.path.join(self.save_dir, '*.jp2')))
        # Parse filenames to get the actual image time
        jp2_datetimes = [datetime_from_filename(filename) for filename in jp2f]
        # Get hek event times
        _, self.hek_times = get_hek_result(self.tstart, self.tend)
        # Initialize the list of missed_downloads list that will contain any hek time entry that failed to download
        self.missed_downloads = []
        for i, time_in in enumerate(self.hek_times):
            # Check that we haven't already downloaded some data for that hek time.
            # Skip them if we did and download only the missing ones.
            hek_time = parse_time(time_in)
            measurements = self.list_missing_measurements(hek_time, jp2f, jp2_datetimes)

            if measurements:
                print('Checking available data for hek time {:s} at index {:d}'.format(
                    hek_time.strftime('%Y/%m/%d %H:%M:%S'), i))
                try:
                    image_files, image_times = download_sdo_images(hek_time, measurements, dt=self.dt, save_path=self.save_dir)
                except ValueError:
                    self.missed_downloads.append(time_in)
                    continue
                for measure_idx, fpath in enumerate(image_files):
                    if fpath is None:
                        print('...Skipped measurement {:s} at time {:s} (too far from hek time) '.format(measurements[measure_idx],
                                                                                 image_times[measure_idx].strftime(
                                                                                     '%Y_%m_%d %H:%M:%S')))
                    else:
                        print('...downloaded file(s) {:s}'.format(fpath))
            else:
                print('skipping already downloaded data for hek time {:s}'.format(
                    hek_time.strftime('%Y/%m/%d %H:%M:%S')))

        if self.missed_downloads:
            with open(self.missed_downloads_csv, 'w') as csvFile:
                writer = csv.writer(csvFile)
                writer.writerows(self.missed_downloads)
            csvFile.close()
            self.download_flag = True
        else:
            self.download_flag = False

        #self.data_cleanup()
        print('Finished download')


    def data_cleanup(self):
        """
        Cleanup the downloaded image to have only complete groups in the training set
        """
        print('data cleanup...')
        hek_time_jp2_map = []
        # List the downloaded images
        downloaded_files = sorted(glob.glob(os.path.join(self.save_dir, '*.jp2')))
        # Parse filenames to get the actual image time
        jp2_datetimes = [datetime_from_filename(filename) for filename in downloaded_files]
        n_incomplete_groups = 0
        rejected_hek_events = []
        # Loop over all the hek times and test if we have all measurements after the download process
        # Reject the whole group if that's not the case
        for i, time_in in enumerate(self.hek_times):
            hek_time = parse_time(time_in)
            # Get the index of files whose time fall within a 2hr window centered on the hek event time
            tmatches = [i for i, file_time in enumerate(jp2_datetimes) if hek_time - self.dt < file_time < hek_time + self.dt]
            if len(tmatches) != len(self.inst_file_map):
                n_incomplete_groups += 1
                # Move files in that incomplete group for rejection in seperate directory
                # TODO: We must also reject the corresponding hek event entry
                rejected_hek_events.append([self.hek_times.index(time_in), time_in])
                for f in tmatches:
                    # The same file may have already been move if the hek time matched again to that file
                    # so check first it exist
                    if os.path.isfile(downloaded_files[f]):
                        shutil.move(downloaded_files[f], os.path.join(self.reject_dir, os.path.basename(downloaded_files[f])))
            else:
                # Append all files of the group to the hek_time <-> jp2 map
                jp2_basenames = [os.path.basename(downloaded_files[t]) for t in tmatches]
                hek_time_jp2_map_entry = [hek_time] + jp2_basenames
                hek_time_jp2_map.append(hek_time_jp2_map_entry)

        # Write hek_time_jp2_map to a csv file
        with open(self.hek_time_jp2_map_csv, 'w') as csvFile:
            writer = csv.writer(csvFile)
            writer.writerows(hek_time_jp2_map)
        csvFile.close()
        # Write the csv of rejected events
        with open(self.rejected_hek_csv, 'w') as csvFile:
            writer = csv.writer(csvFile)
            writer.writerows(rejected_hek_events)
        csvFile.close()
        print('data cleanup finished.')



    def make_labels(self):
        """
        Create feature masks from the downloaded images for the start and end time defined for self. Save them to disk
        in the subdirectory 'label_masks' under the image directory.
        """
        if self.do_plot:
            matplotlib.rcParams.update({'font.size': 18})
        # Get a list of existing files (if any)
        jp2f = sorted(glob.glob(os.path.join(self.save_dir, '*.jp2')))
        # Parse filenames to get the actual image time
        jp2_datetimes = [datetime_from_filename(filename) for filename in jp2f]

        # Use the curated hek results from the cleanup pass instead of querying the hek again.
        # Otherwise this will use an uncurated list of hek times, and inconsistent map of jp2 <-> hek times
        result, times = get_hek_result(self.tstart, self.tend)
        # Read the csv for rejected events
        self.rejected_hek_events = []
        with open(self.rejected_hek_csv) as csvfile:
            readcsv = csv.reader(csvfile, delimiter=',')
            for row in readcsv:
                self.rejected_hek_events.append(row[1])
        csvfile.close()
        # Filter out the rejected events
        for time in self.rejected_hek_events:
            idx = times.index(time)
            del result[idx]
            del times[idx]


        ch = [elem for elem in result if elem['event_type'] == 'CH']
        ar = [elem for elem in result if elem['event_type'] == 'AR']
        ss = [elem for elem in result if elem['event_type'] == 'SS']

        mask_time_map = []

        for i, time_in in enumerate(times):
            hek_time = parse_time(time_in)
            # Get closest image
            nearest_datetime = nearest(jp2_datetimes, hek_time)
            nearest_file = jp2f[jp2_datetimes.index(nearest_datetime)]
            print('...processing hek time: {:s} at index {:d} '.format(hek_time.strftime('%Y/%m/%d %H:%M:%S'), i))
            print('......using nearest image at time: {:s}'.format(nearest_datetime.strftime('%Y/%m/%d %H:%M:%S')))

            ch_list = [elem for elem in ch if elem['event_starttime'] == time_in]
            ar_list = [elem for elem in ar if elem['event_starttime'] == time_in]
            ss_list = [elem for elem in ss if elem['event_starttime'] == time_in]
            # The above 3 lists have typically only 1 that's not empty. Let's explicitly tell to not process any empty label list.
            if ch_list:
                ch_mask, ch_file_path = gen_label_mask(ch_list, nearest_file, hek_time, 'CH', save_path=self.label_save_dir, do_plot=self.do_plot)
                mask_time_map.append([os.path.basename(ch_file_path), time_in])
            if ar_list:
                ar_mask, ar_file_path = gen_label_mask(ar_list, nearest_file, hek_time, 'AR', save_path=self.label_save_dir, do_plot=self.do_plot)
                mask_time_map.append([os.path.basename(ar_file_path), time_in])
            if ss_list:
                ss_mask, ss_file_path = gen_label_mask(ss_list, nearest_file, hek_time, 'SS', save_path=self.label_save_dir, do_plot=self.do_plot)
                mask_time_map.append([os.path.basename(ss_file_path), time_in])

        # Write mask_time_map to a csv file
        with open(self.mask_hek_time_map_csv, 'w') as csvFile:
            writer = csv.writer(csvFile)
            writer.writerows(mask_time_map)
        csvFile.close()


def gen_label_mask(label_list, image_filepath, hek_time, label, save_path=None, do_plot=False):
    """
    Create a binary mask of feature locations. If there is no feature, an empty mask is returned

    :param label_list: List of HEK dictionaries
    :param image_filepath: file path of image closest to label. Its header is used for converting hek event coordinates.
    :param hek_time: Datetime used to write the filenames of the mask files
    :param label: The event type being processed: either 'CH', 'AR' or 'SS'.
    :param save_path: directory of the mask files and optionally the figures.
    :param do_plot: True will plot figures and save them to save_path
    :return: binary mask of feature and list of polygon vertices tuples.
    """

    print('.........processing ' + label)
    jp2_header = get_header(image_filepath)[0]
    mask_shape = (4096, 4096)
    jp2_shape = (jp2_header['NAXIS1'], jp2_header['NAXIS2'])
    mask = np.zeros(mask_shape, dtype=np.int)
    dummy_array = np.empty(mask_shape, dtype=np.int)
    verts_yx_list = []

    if jp2_shape != mask.shape:
        raise ValueError('Mask and WCS array shapes do not agree.')

    aia_map = Map(dummy_array, jp2_header)

    fig, ax = None, None
    if do_plot:
        fig = plt.figure(figsize=(10, 10))
        ax = plt.subplot(projection=aia_map)
        aia_map.plot(axes=ax)
    # parsing boundary coordinate string for each feature
    for labels in label_list:

        p1 = labels["hpc_boundcc"][9:-2]
        p2 = p1.split(',')
        p3 = [v.split(" ") for v in p2]
        print(labels["hpc_boundcc"])
        print(p3)

        boundary_coords = SkyCoord([(float(v[0]), float(v[1])) * u.arcsec for v in p3], frame=aia_map.coordinate_frame)

        # Vertices  of feature
        pixel_verts = aia_map.world_to_pixel(boundary_coords)
        print(pixel_verts)
        verts_x = np.array([x for x in pixel_verts[0].value if not np.isnan(x)])
        verts_y = np.array([y for y in pixel_verts[1].value if not np.isnan(y)])
        verts_yx = np.round(np.array((verts_y, verts_x)).T).astype(np.int)
        verts_yx_list.append(verts_yx)
        # The mask is populated in-place -> accross different instances of the same hek event, the mask builds itself.
        fill_polygon(verts_yx, mask)

        if do_plot:
            ax.plot_coord(boundary_coords, color='r')

    mask_file_path = write_mask(mask, hek_time, label, save_path=save_path)


    if do_plot:
        label_map = Map(mask, jp2_header)
        label_map.plot(axes=ax)
        label_map.draw_limb()
        ax.set_title(label + ' @ hek_time: ' + hek_time.strftime('%Y/%m/%d %H:%M:%S'))
        # TODO: fix the tight_layout issue
        #plt.tight_layout()
        plt.savefig(os.path.join(save_path, os.path.basename(image_filepath)[0:-4] + '_plot_' + label + '.png'))
        plt.close(fig)


    return mask, mask_file_path





def datetime_from_filename(filepath):
    basename = os.path.basename(filepath)[0:20]
    file_time_str= basename[:4] + '-' + basename[5:7] + '-' + basename[8:10] + 'T' + basename[12:14] + ':' + basename[15:17] \
               + ':' + basename[18:20]
    file_datetime = parse_time(file_time_str)
    return file_datetime


def get_hek_result(time_start, time_end):
    client = hek.HEKClient()
    result = client.search(hek.attrs.Time(time_start, time_end), hek.attrs.FRM.Name == 'SPoCA')  # CH and AR
    result += client.search(hek.attrs.Time(time_start, time_end), hek.attrs.FRM.Name == 'EGSO_SFC')  # SS
    times = list(set([elem["event_starttime"] for elem in result]))
    times.sort()
    return result, times


def download_sdo_images(time_in, measurements, dt, save_path=''):
    """
    Download a complete set of the SDO images in AIA and HMI for a given time, with optional tolerance on time difference

    :param time_in: requested datetime for JP2 image download.
    :param measurements: list of string of measurement names: AIA wavelength: '193', '94',... or HMI segment name: 'continuum' or 'magnetogram'
    :param dt: time difference tolerated between requested time and available image time so that time_in - dt < actual image time < time_in + dt
    :param save_path: save path for downloaded images
    :return: full file path of example AIA image downloaded (335 channel). Set to None if time is off limits
    """

    hv = helioviewer.HelioviewerClient()

    filepaths = []
    image_times = []
    for measure in measurements:
        if measure is not 'continuum' and measure is not 'magnetogram':
            kwargs = {'observatory': 'SDO', 'instrument': 'AIA', 'detector': 'AIA', 'measurement': measure}
        else:
            kwargs = {'observatory': 'SDO', 'instrument': 'HMI', 'detector': 'HMI', 'measurement': measure}

        if dt is not None:
            # Check how far requested time in metadata is from requested hek time
            metadata = hv.get_closest_image(time_in, **kwargs)
            image_time = metadata['date']
            if time_in - dt < image_time < time_in + dt:
                filepath = hv.download_jp2(time_in, directory=save_path, **kwargs)
            else:
                # Do not download if actual image time is too far from requested time
                filepath = None
        else:
            image_time = None
            filepath = hv.download_jp2(time_in, directory=save_path, **kwargs)

        filepaths.append(filepath)
        image_times.append(image_time)

    return filepaths, image_times


def get_header(filepath):
    """
    Reads the header from the file and sanitizes it to Fits standards.

    Parameters
    ----------
    filepath : `str`
        The file to be read

    Returns
    -------
    headers : list
        A list of headers read from the file
    """

    sunpyhead = read_file_header(filepath)

    for subhead in sunpyhead:
        # Remove newlines from comment and history
        if 'comment' in subhead:
            subhead['comment'] = subhead['comment'].replace("\n", "")
        if 'history' in subhead:
            subhead['history'] = subhead['history'].replace("\n", "")

        badkeys = []

        # dumps header keywords that are NaN
        for key, value in subhead.items():
            if type(value) in (int, float):
                if np.isnan(value):
                    badkeys += [key]

        for key in badkeys:
            subhead.pop(key, None)

    return sunpyhead


def mahotas_fill_polygon(verts_yx, mask_shape):
    mahotas_poly_mask = np.zeros(mask_shape, dtype=np.int)
    fill_polygon(verts_yx, mahotas_poly_mask)
    return mahotas_poly_mask


def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))


def write_mask(mask_in, label_date, label='Mask', save_path=''):
    """
    Save the feature mask as an image

    :param mask_in: Mask of labels as numpy array
    :param label_date: time of the label detection as Datetime
    :param label: String of feature type e.g. "SS", "CH", or "AR"
    :param save_path: save path for saving mask as image
    :return: path to the saved file
    """

    label_mask_name = label_date.strftime("%Y_%m_%d__%H_%M_%S")+'__'+label+'_mask.npz'
    save_mask_name = os.path.join(save_path, label_mask_name)

    np.savez_compressed(save_mask_name, mask_in.astype(np.bool))

    return save_mask_name



