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

import numpy as np

from imageio import imwrite

import datetime
from datetime import timedelta
import os, glob, shutil
from mahotas.polygon import fill_polygon
import csv

TIME_FORMAT = '%Y/%m/%d %H:%M:%S'
FILE_TIME_FORMAT = '%Y_%m_%dT%H_%M_%S'

class Jp2ImageDownload:

    def __init__(self, save_dir, dt=timedelta(minutes=30),
                 tstart='2012/05/30 23:59:59', tend='2012/05/31 23:59:59'):
        """
        Set up initial start and ending times and save directory.

        :param save_dir: parent directory to contain the year/month-based-name subdirectories of images, masks, and csv files
        :param full_image_set: boolean if true will run through the entire SDO catalogue
        :param tstart: starting time string of any format read by parse_time
        :param tend: ending time string of any format read by parse_time
        :param dt: limit the time gap between downloaded image and requested hek time. Set to None for no limit.
        """

        # if full_image_set:
        #     tstart = '2010/05/31 23:59:59'
        #     tend = '2018/12/31 23:59:59'

        self.tstart = parse_time(tstart)
        self.tend = parse_time(tend)
        # Build subdirectory name from queried year and month
        subdir = parse_time(tstart).strftime('%Y_%m')
        self.save_dir = os.path.join(save_dir, subdir)
        self.save_dir_jp2 = os.path.join(save_dir, subdir, 'jp2')
        os.makedirs(self.save_dir_jp2, exist_ok=True)
        # Subdirectory for discarding the files from incomplete sets
        self.reject_dir = os.path.join(self.save_dir, 'incomplete_set')
        # Create the "rejection" subdirectory to move the images from incomplete groups
        os.makedirs(self.reject_dir, exist_ok=True)
        # Directory to save the label masks
        self.label_save_dir = os.path.join(self.save_dir, 'label_masks')
        os.makedirs(self.label_save_dir, exist_ok=True)

        self.hek_times = None
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
        self.inst_file_map.sort()
        # store the list of files in the save directory
        self.jp2f = []
        self.jp2_datetimes = []
        self.jp2_measurements = []
        self.generic_filenames = []
        self.mask_hek_time_map_csv = os.path.join(self.save_dir, 'label_jp2_map.csv')
        self.hek_time_jp2_map_csv = os.path.join(self.save_dir, 'hek_time_jp2_map.csv')
        self.rejected_hek_csv = os.path.join(self.save_dir, 'rejected_hek.csv')
        self.missed_downloads_csv = os.path.join(self.save_dir, 'missed_downloads.csv')
        self.blank_hek_events_csv = os.path.join(self.save_dir, 'blank_hek_events.csv')
        self.hek_time_jp2_map = []
        self.rejected_hek_events = []
        self.missed_downloads = []
        self.blank_hek_events = []
        self.download_flag = True
        self.do_plot = True

        matplotlib.rcParams.update({'font.size': 18})


    def download_images(self):
        """
        Run through the complete set of dates as determined in the self.__init__ call by day
        Sets up a subdirectory tree to organize the images

        """
        # IF any, move previously discarded files to the save directory (this to be consistent with the cleanup pass)
        discared_files = sorted(glob.glob(os.path.join(self.reject_dir, '*.jp2')))
        for f in discared_files:
            print('restoring discarded file: ' + os.path.basename(f))
            shutil.move(f, self.save_dir_jp2)

        # for ii, download_date in enumerate(self.date_list):
        # First download the images. Feature extraction will be done separately

        # Get a list of existing files (if any)
        self.jp2f = sorted(glob.glob(os.path.join(self.save_dir_jp2, '*.jp2')))
        # Parse filenames to get the actual image time
        self.jp2_datetimes = [datetime_from_filename(filename).strftime(FILE_TIME_FORMAT) for filename in self.jp2f]
        self.jp2_measurements = [inst_from_filename(filename) for filename in self.jp2f]
        self.generic_filenames = [date + '__' + measurement for date, measurement in zip(self.jp2_datetimes, self.jp2_measurements)]
        # Get hek event times
        if self.missed_downloads:
            self.hek_times = self.missed_downloads.copy()
        else:
            print('Querying times to the HEK...')
            _, self.hek_times = get_hek_result(self.tstart, self.tend)
        # Initialize the list of missed_downloads list that will contain any hek time entry that failed to download
        self.missed_downloads = []


        downloaded_events = []
        if os.path.isfile(self.hek_time_jp2_map_csv):
            with open(self.hek_time_jp2_map_csv) as csvfile:
                readcsv = csv.reader(csvfile, delimiter=',')
                for row in readcsv:
                    downloaded_events.append(row)
            csvfile.close()

        downloaded_event_times = [row[0] for row in downloaded_events]

        for i, time_in in enumerate(self.hek_times):
            # Check that we haven't already downloaded some data for that hek time.
            # Skip them if we did and download only the missing ones.
            hek_time = parse_time(time_in)

            if time_in not in downloaded_event_times: # TODO: Check instead a previous csv file to see if this hek time has been populated
                print('Checking available data for hek time {:s} at index {:d}'.format(
                    hek_time.strftime('%Y/%m/%d %H:%M:%S'), i))
                try:
                    image_files, image_times = self.download_sdo_images(hek_time, self.measurements_req, dt=self.dt, save_path=self.save_dir_jp2)
                    if len(image_files) < len(self.measurements_req):
                        # That event is incomplete, must be rejected.
                        print('Rejecting event (incomplete)')
                        self.rejected_hek_events.append(hek_time)
                    else:
                        jp2_basenames = [os.path.basename(f) for f in image_files]
                        downloaded_events.append([time_in] + jp2_basenames)
                except ValueError:
                    # This includes JSONDecodeError, occurs when something between the client and the server goes wrong.
                    # This should be added to the missed download, which will be subject to new download attempts
                    print('Exception raised by helioviewer client. Appending to missed_downloads.')
                    self.missed_downloads.append(time_in)
                    continue
                # TODO: Catch also what's thrown by the helioviewer client when the json response does not contain a valid key
                except KeyError:
                    print('Helioviewer KeyError. Skipping.')
                    # And these ones are just rubbish data, we must NOT attempt a new download
                    continue
            else:
                print('skipping already downloaded data for hek time {:s}'.format(
                    hek_time.strftime('%Y/%m/%d %H:%M:%S')))

        # Chronological ordering of downloaded_events
        downloaded_events.sort()


        with open(self.hek_time_jp2_map_csv, 'w') as csvFile:
            writer = csv.writer(csvFile)
            writer.writerows(downloaded_events)
        csvFile.close()


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


    def download_sdo_images(self, time_in, measurements, dt, save_path=''):
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

            requested_file_measurement = '{:s}_{:s}_{:s}_{:s}'.\
                format(kwargs['observatory'], kwargs['instrument'], kwargs['detector'], measure)

            # Check how far requested time in metadata is from requested hek time
            metadata = hv.get_closest_image(time_in, **kwargs)
            image_time = metadata['date']
            # Build the generic filename that match this
            generic_fname = image_time.strftime(FILE_TIME_FORMAT) + '__' + requested_file_measurement

            if generic_fname not in self.generic_filenames:
                if time_in - dt < image_time < time_in + dt:
                    filepath = hv.download_jp2(time_in, directory=save_path, overwrite=True, **kwargs)
                    print('...downloaded file(s) {:s}'.format(filepath))
                    filepaths.append(filepath)
                else:
                    # Do not download if actual image time is too far from requested time
                    print('...Skipped measurement {:s} at time {:s} (too far from hek time) '.format(
                        measure, image_time.strftime(TIME_FORMAT)))
            else:
                fidx = self.generic_filenames.index(generic_fname)
                filepaths.append(self.jp2f[fidx])
                print('...skipping already downloaded file with generic name: ' + generic_fname)

            image_times.append(image_time)

        return filepaths, image_times


    def data_cleanup(self):
        """
        Cleanup the downloaded image to have only complete groups in the training set
        """


        print('data cleanup...')

        if not self.hek_time_jp2_map:
            with open(self.hek_time_jp2_map_csv) as csvfile:
                readcsv = csv.reader(csvfile, delimiter=',')
                for row in readcsv:
                    self.hek_time_jp2_map.append(row)
            csvfile.close()

        # List the downloaded images
        downloaded_files = sorted(glob.glob(os.path.join(self.save_dir_jp2, '*.jp2')))
        ## move the file present in the directory but that are not listed in hek_time_jp2_map.
        # list the valid files that have an entry in the hek time as a complete set
        valid_sets = [row[1:] for row in self.hek_time_jp2_map]
        # Need to unroll all of them in one single flat list, convert them back with a full path.
        flat_valid_files = [os.path.join(self.save_dir_jp2, file) for file_set in valid_sets for file in file_set]
        # Now for each file in the directory that does not match any file in flat_valid_files, discard them
        for file in downloaded_files:
            if file not in flat_valid_files:
                shutil.move(file, os.path.join(self.reject_dir, os.path.basename(file)))


        # Write hek_time_jp2_map to a csv file
        # with open(self.hek_time_jp2_map_csv, 'w') as csvFile:
        #     writer = csv.writer(csvFile)
        #     writer.writerows(hek_time_jp2_map)
        # csvFile.close()
        # # Write the csv of rejected events
        # with open(self.rejected_hek_csv, 'w') as csvFile:
        #     writer = csv.writer(csvFile)
        #     writer.writerows(rejected_hek_events)
        # csvFile.close()
        print('data cleanup finished.')


    def make_labels(self):
        """
        Create feature masks from the downloaded images for the start and end time defined for self. Save them to disk
        in the subdirectory 'label_masks' under the image directory.
        """


        # Get a list of existing files (if any)
        self.jp2f = sorted(glob.glob(os.path.join(self.save_dir_jp2, '*.jp2')))
        # Parse filenames to get the actual image time
        jp2_datetimes = [datetime_from_filename(filename) for filename in self.jp2f]
        # Use the curated hek results from the cleanup pass instead of querying the hek again.
        # Otherwise this will use an uncurated list of hek times, and inconsistent map of jp2 <-> hek times
        results, _ = get_hek_result(self.tstart, self.tend)
        # Read the csv for rejected events
        # self.rejected_hek_events = []
        # with open(self.rejected_hek_csv) as csvfile:
        #     readcsv = csv.reader(csvfile, delimiter=',')
        #     for row in readcsv:
        #         self.rejected_hek_events.append(row[1])
        # csvfile.close()
        # # Filter out the rejected events
        # for time in self.rejected_hek_events:
        #     idx = times.index(time)
        #     # TODO: fix the inconsistent indexing between times and results
        #     del results[idx]
        #     del times[idx]

        # Read the mapping of hek times to jp2 files
        hek_time_jp2_map = []
        with open(self.hek_time_jp2_map_csv) as csvfile:
            readcsv = csv.reader(csvfile, delimiter=',')
            for row in readcsv:
                hek_time_jp2_map.append(row)
        csvfile.close()

        times = [row[0] for row in hek_time_jp2_map]

        ch = [elem for elem in results if elem['event_type'] == 'CH']
        ar = [elem for elem in results if elem['event_type'] == 'AR']
        ss = [elem for elem in results if elem['event_type'] == 'SS']

        mask_time_map = []

        for i, time_in in enumerate(times):
            hek_time = parse_time(time_in)
            # Get all jp2 files that map to that hek time
            jp2f_at_hek_time = hek_time_jp2_map[i][1:]

            # Get closest image
            nearest_datetime = nearest(jp2_datetimes, hek_time)
            nearest_file = self.jp2f[jp2_datetimes.index(nearest_datetime)]
            print('...processing hek time: {:s} at index {:d} '.format(hek_time.strftime('%Y/%m/%d %H:%M:%S'), i))
            print('......using nearest image at time: {:s}'.format(nearest_datetime.strftime('%Y/%m/%d %H:%M:%S')))
            # Extract metadata for each class and at the specific date time_in
            ch_list = [elem for elem in ch if elem['event_starttime'] == time_in]
            ar_list = [elem for elem in ar if elem['event_starttime'] == time_in]
            ss_list = [elem for elem in ss if elem['event_starttime'] == time_in]
            # The above 3 lists have typically only 1 that's not empty. Let's explicitly tell to not process any empty label list.
            if ch_list:
                ch_mask, ch_file_path, ch_blanks = gen_label_mask(ch_list, nearest_file, hek_time, 'CH', save_path=self.label_save_dir, do_plot=self.do_plot)
                mask_time_map.append([os.path.basename(ch_file_path), time_in] + jp2f_at_hek_time)
                self.blank_hek_events += ch_blanks
            if ar_list:
                ar_mask, ar_file_path, ar_blanks = gen_label_mask(ar_list, nearest_file, hek_time, 'AR', save_path=self.label_save_dir, do_plot=self.do_plot)
                mask_time_map.append([os.path.basename(ar_file_path), time_in] + jp2f_at_hek_time)
                self.blank_hek_events += ar_blanks
            if ss_list:
                ss_mask, ss_file_path, ss_blanks = gen_label_mask(ss_list, nearest_file, hek_time, 'SS', save_path=self.label_save_dir, do_plot=self.do_plot)
                mask_time_map.append([os.path.basename(ss_file_path), time_in] + jp2f_at_hek_time)
                self.blank_hek_events += ss_blanks

        # Write mask_time_map to a csv file
        with open(self.mask_hek_time_map_csv, 'w') as csvFile:
            writer = csv.writer(csvFile)
            writer.writerows(mask_time_map)
        csvFile.close()

        # Create the csv file that will contain the "blank" hek events, i.e, event that have a hek entry but no
        # hpc_boundcc coordinates
        with open(self.blank_hek_events_csv, 'w') as outcsv:
            writer = csv.writer(outcsv)
            writer.writerow(['frm_specificid', 'event_starttime'])
            writer.writerows(self.blank_hek_events)
        outcsv.close()




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
    blank_hek_elems = []

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
        if len(p3) > 2:
            # Convert coordinates of polygon vertices from helioprojective cartesian (HPC) to pixel "image" coordinates.
            boundary_coords = SkyCoord([(float(v[0]), float(v[1])) * u.arcsec for v in p3], frame=aia_map.coordinate_frame)
            pixel_verts = aia_map.world_to_pixel(boundary_coords)
            #print(pixel_verts)
            verts_x = np.array([x for x in pixel_verts[0].value if not np.isnan(x)])
            verts_y = np.array([y for y in pixel_verts[1].value if not np.isnan(y)])
            verts_yx = np.round(np.array((verts_y, verts_x)).T).astype(np.int)
            verts_yx_list.append(verts_yx)
            # The mask is populated in-place -> accross different instances of the same hek event, the mask builds itself.
            fill_polygon(verts_yx, mask)

            if do_plot:
                ax.plot_coord(boundary_coords, color='r')
        else:
            blank_hek_elems.append([labels['frm_specificid'], labels['event_starttime']])
    # In the very rare case where the HPC coordinates of all elements in the hek result are blank, raise that specifically to let us know.
    if len(blank_hek_elems) == len(label_list):
        raise Exception('All elements at hek_time {:s} are blank for label {:s}'.format(hek_time.strftime(TIME_FORMAT), label))
        # TODO: handle such error in a non-interruptive manner.

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


    return mask, mask_file_path, blank_hek_elems


def datetime_from_filename(filepath):
    basename = os.path.basename(filepath)
    file_time_str= basename[:4] + '-' + basename[5:7] + '-' + basename[8:10] + 'T' + basename[12:14] + ':' + basename[15:17] \
               + ':' + basename[18:20]
    file_datetime = parse_time(file_time_str)
    return file_datetime


def inst_from_filename(filepath):
    basename = os.path.basename(filepath)
    idx = basename.index('SDO')
    inst_str = basename[idx:-4]
    return inst_str



def get_hek_result(time_start, time_end):
    client = hek.HEKClient()
    results = client.search(hek.attrs.Time(time_start, time_end), hek.attrs.FRM.Name == 'SPoCA')  # CH and AR
    results += client.search(hek.attrs.Time(time_start, time_end), hek.attrs.FRM.Name == 'EGSO_SFC')  # SS
    times = list(set([elem["event_starttime"] for elem in results]))
    times.sort()
    return results, times



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



