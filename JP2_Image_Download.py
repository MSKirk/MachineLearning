from sunpy.net import hek, helioviewer
from sunpy.time import parse_time
from sunpy.coordinates import frames
from sunpy.io import read_file_header

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import numpy as np

from imageio import imwrite

import datetime
from datetime import timedelta
import os, glob, shutil
from mahotas.polygon import fill_polygon





class Jp2ImageDownload:

    def __init__(self, parent_dir='', full_image_set=False, dt=timedelta(minutes=30),
                 tstart='2012/05/30 23:59:59', tend='2012/05/31 23:59:59'):
        """
        Set up initial start and ending times and save directory.

        :param parent_dir: Directory to contain the full subdirectory tree of images and masks
        :param full_image_set: boolean if true will run through the entire SDO catalogue
        :param tstart: starting time string of any format read by parse_time
        :param tend: ending time string of any format read by parse_time
        :param dt: limit the time gap between downloaded image and requested hek time. Set to None for no limit.
        """

        if full_image_set:
            tstart = '2010/05/31 23:59:59'
            tend = '2018/12/31 23:59:59'

        if parent_dir:
            self.parent_dir = os.path.abspath(parent_dir)
        else:
            self.parent_dir = os.path.abspath(os.curdir)
            print('Using current directory as save directory: {:s}'.format(self.parent_dir))

        # Format string for the date-based subdirectories created under self.save_dir
        self.data_subdir_format = '%Y_%m_%d'
        # Subdirectory for discarding the files from incomplete sets
        self.reject_subdir = 'incomplet_set'

        self.tstart = parse_time(tstart)
        self.tend = parse_time(tend)

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

        # Let's separate the image download from the initialization
        # self.download_images()


    def gen_date_list(self):
        """
        Create list of datetimes by day

        :return: list of datetimes broken into day segments and their string representation
        """
        time_between = self.tend - self.tstart
        ndays = max(1, time_between.days)
        date_list = [self.tstart + timedelta(days=dd) for dd in range(0, ndays)]
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

        for ii, download_date in enumerate(self.date_list):
            subdir = download_date.strftime(self.data_subdir_format)
            save_path = os.path.join(self.parent_dir, subdir)
            os.makedirs(save_path, exist_ok=True)
            # First download the images. Feature extraction will be done separately
            tstart = download_date
            tend = tstart + timedelta(days=1)

            # Get a list of existing files (if any)
            jp2f = sorted(glob.glob(os.path.join(save_path, '*.jp2')))
            # Parse filenames to get the actual image time
            jp2_datetimes = [datetime_from_filename(filename) for filename in jp2f]
            # Get hek event times
            _, times = get_hek_result(tstart, tend)

            for i, time_in in enumerate(times):
                # Check that we haven't already downloaded some data for that hek time.
                # Skip them if we did and download only the missing ones.
                hek_time = parse_time(time_in)
                measurements = self.list_missing_measurements(hek_time, jp2f, jp2_datetimes)

                if measurements:
                    print('Checking available data for hek time {:s} at index {:d}'.format(
                        hek_time.strftime('%Y/%m/%d %H:%M:%S'), i))
                    image_files, image_times = download_sdo_images(hek_time, measurements, dt=self.dt, save_path=save_path)
                    for measure_idx, fpath in enumerate(image_files):
                        if fpath is None:
                            print('...Skipped measurement {:s} at time {:s} '.format(measurements[measure_idx],
                                                                                     image_times[measure_idx].strftime(
                                                                                         '%Y_%m_%d %H:%M:%S')))
                        else:
                            print('...downloaded file(s) {:s}'.format(fpath))
                else:
                    print('skipping already downloaded data for hek time {:s}'.format(
                        hek_time.strftime('%Y/%m/%d %H:%M:%S')))

            self.data_cleanup(times, save_path)
            print('Finished download of group {:d}/{:d}'.format(ii+1, len(self.date_list)))


    def data_cleanup(self, times, save_path):
        """
        Cleanup the downloaded image to have only complete groups in the training set

        :param times: list of all hek times
        :param save_path: directory where the jpeg2000 files were downloaded
        :return:
        """

        # List the downloaded images
        downloaded_files = sorted(glob.glob(os.path.join(save_path, '*.jp2')))
        # Parse filenames to get the actual image time
        jp2_datetimes = [datetime_from_filename(filename) for filename in downloaded_files]
        n_incomplete_groups = 0
        # Create the "rejection" subdirectory to move the image.
        reject_dir = os.path.join(save_path, self.reject_subdir)
        os.makedirs(reject_dir, exist_ok=True)
        # Loop over all the hek times and test if we have all measurements after the download process
        # Reject the whole group if that's not the case
        for i, time_in in enumerate(times):
            hek_time = parse_time(time_in)
            # Get the index of files whose time fall within a 2hr window centered on the hek event time
            tmatches = [i for i, file_time in enumerate(jp2_datetimes) if hek_time - self.dt < file_time < hek_time + self.dt]
            if len(tmatches) != len(self.inst_file_map):
                n_incomplete_groups += 1
                # Move files in that incomplete group for rejection in seperate directory
                for f in tmatches:
                    shutil.move(downloaded_files[f], reject_dir)


    def make_labels(self, save_dir):
        """
        Download jp2 SDO images and create feature masks for a given time period

        :param save_dir: directory where label-masks are saved
        :return:
        """

        for ii, download_date in enumerate(self.date_list):
            subdir = download_date.strftime(self.data_subdir_format)
            save_path = os.path.join(self.parent_dir, subdir)
            os.makedirs(save_path, exist_ok=True)
            # First download the images. Feature extraction will be done separately
            tstart = download_date
            tend = tstart + timedelta(days=1)

            # Get a list of existing files (if any)
            jp2f = sorted(glob.glob(os.path.join(save_path, '*.jp2')))
            # Parse filenames to get the actual image time
            jp2_datetimes = [datetime_from_filename(filename) for filename in jp2f]

            result, times = get_hek_result(tstart, tend)

            ch = [elem for elem in result if elem['event_type'] == 'CH']
            ar = [elem for elem in result if elem['event_type'] == 'AR']
            ss = [elem for elem in result if elem['event_type'] == 'SS']

            for i, time_in in enumerate(times):
                #image_file = self.get_all_sdo_images(time_in, save_path=save_dir)
                # Get closest image
                nearest_datetime = nearest(jp2_datetimes, parse_time(time_in))
                nearest_file = jp2f[jp2_datetimes.index(nearest_datetime)]

                ch_list = [elem for elem in ch if elem['event_starttime'] == time_in]
                ar_list = [elem for elem in ar if elem['event_starttime'] == time_in]
                ss_list = [elem for elem in ss if elem['event_starttime'] == time_in]

                ch_mask, _ = self.gen_label_mask(time_in, ch_list, nearest_file)
                ar_mask, _ = self.gen_label_mask(time_in, ar_list, nearest_file)
                ss_mask, _ = self.gen_label_mask(time_in, ss_list, nearest_file)



                # self.write_mask(ch_mask, time_in, 'CH', save_path=save_dir)
                # self.write_mask(ar_mask, time_in, 'AR', save_path=save_dir)
                # self.write_mask(ss_mask, time_in, 'SS', save_path=save_dir)



    def gen_label_mask(self, label_time, label_list, image_filepath):
        """
        Create a binary mask of feature locations. If there is no feature, an empty mask is returned

        :param label_time: Time string of the pattern detection
        :param label_list: List of HEK dictionaries
        :param image_filepath: file path of image in which the feature is detected
        :return: binary mask of feature
        """

        aia_wcs = WCS(self.get_header(image_filepath)[0])

        mask = np.zeros([4096,4096])
        verts_yx_list = []

        if aia_wcs.array_shape != mask.shape:
            raise ValueError('Mask and WCS array shapes do not agree.')

        # parsing boundary coordinate string for each feature
        # Even if there is no feature, an empty mask is returned
        # This is ok but need to be classified as
        for labels in label_list:

            p1 = labels["hpc_boundcc"][9:-2]
            p2 = p1.split(',')
            p3 = [v.split(" ") for v in p2]

            label_date = parse_time(label_time)

            label_boundary = SkyCoord([(float(v[0]), float(v[1])) * u.arcsec for v in p3], obstime=label_date,
                                        frame=frames.Helioprojective)

            # Vertices  of feature
            # R: This can contain NaN. The output of SkyCoord needs to be sanitized
            pixel_vertex = label_boundary.to_pixel(aia_wcs)
            verts_x = np.array([x for x in pixel_vertex[0] if not np.isnan(x)])
            verts_y = np.array([y for y in pixel_vertex[1] if not np.isnan(y)])
            verts_yx = np.round(np.array((verts_y, verts_x)).T).astype(np.int)
            verts_yx_list.append(verts_yx)
            fill_polygon(verts_yx, mask)

            # Need to add in contours between vertices
            # Need to add in feature filling in mask

        return mask, verts_yx_list





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


def write_mask(mask_in, label_time, label='Mask', save_path=''):
    """
    Save the feature mask as an image

    :param mask_in: Mask of feature as np array
    :param label_time: time of the label detection as DateTime
    :param label: String of feature type e.g. "SS", "CH", or "AR"
    :param save_path: save path for saving mask as image
    :return: none
    """

    label_date = parse_time(label_time)

    label_mask_name = label_date.strftime("%Y_%m_%d__%H_%M_%S")+'__'+label+'_mask.jp2'
    save_mask_name = os.path.join(save_path, label_mask_name)

    imwrite(save_mask_name, mask_in.astype('uint8'))


