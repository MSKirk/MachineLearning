from sunpy.net import hek, helioviewer
from sunpy.time import parse_time
from sunpy.coordinates import frames
from sunpy.io import read_file_header
from sunpy.map import Map

import astropy.units as u
from astropy.coordinates import SkyCoord
#from astropy.wcs import WCS

import numpy as np

from imageio import imwrite

import datetime
import os


class Jp2ImageDownload:

    def __init__(self, save_dir='', full_image_set=False, tstart='2012/05/30 23:59:59', tend='2012/05/31 23:59:59'):
        """
        Set up initial start and ending times and save directory.

        :param save_dir: Directory to contain the full subdirectory tree of images and masks
        :param full_image_set: boolean if true will run through the entire SDO catalogue
        :param tstart: starting time string of any format read by parse_time
        :param tend: ending time string of any format read by parse_time
        """

        if full_image_set:
            tstart = '2010/05/31 23:59:59'
            tend = '2018/12/31 23:59:59'

        if save_dir:
            self.save_dir = os.path.abspath(save_dir)
        else:
            self.save_dir = os.path.abspath(os.curdir)

        self.tstart = parse_time(tstart)
        self.tend = parse_time(tend)

        self.date_list = self.gen_date_list()
        self.date_string_list = [tt.strftime('%Y/%m/%d %H:%M:%S') for tt in self.date_list]
        self.download_images()

    def download_images(self):
        """
        Run through the complete set of dates as determined in the self.__init__ call by day
        Sets up a subdirectory tree to organize the images
        """

        for ii, download_date in enumerate(self.date_list):
            directories = download_date.strftime("%Y/%m/%d")
            save_path = os.path.join(self.save_dir, directories)
            os.makedirs(save_path, exist_ok=True)

            self.get_feature_images(self.date_string_list[ii + 1], self.date_string_list[ii], save_path)


    def get_all_sdo_images(self, time_in, save_path=''):
        """
        Download a complete set of the SDO images in AIA and HMI for a given time

        :param time_in: string time for JP2 image download
        :param save_path: save path for downloaded images
        :return: full file path of example AIA image downloaded (335 channel)
        """

        hv = helioviewer.HelioviewerClient()

        wavelnths = ['1600', '1700', '094', '131', '171', '211', '304', '335', '193']
        measurement = ['continuum', 'magnetogram']

        for wav in wavelnths:
            aia_filepath = hv.download_jp2(time_in, observatory='SDO', instrument='AIA', detector='AIA', measurement=wav,
                                       directory=save_path, overwrite=True)

        for measure in measurement:
            hmi_filepath = hv.download_jp2(time_in, observatory='SDO', instrument='HMI', detector='HMI',
                                       measurement=measure, directory=save_path, overwrite=True)

        return aia_filepath

    def get_feature_images(self, time_start, time_end, save_dir):
        """
        Download jp2 SDO images and create feature masks for a given time period

        :param time_start: String start time for HEK query
        :param time_end: string end time for HEK query
        :param save_dir: directory where images are saved
        :return:
        """

        client = hek.HEKClient()
        result = client.search(hek.attrs.Time(time_start, time_end), hek.attrs.FRM.Name == 'SPoCA')  # CH and AR
        result += client.search(hek.attrs.Time(time_start, time_end), hek.attrs.FRM.Name == 'EGSO_SFC')  # SS

        times = list(set([elem["event_starttime"] for elem in result]))
        times.sort()

        ch = [elem for elem in result if elem['event_type'] == 'CH']
        ar = [elem for elem in result if elem['event_type'] == 'AR']
        ss = [elem for elem in result if elem['event_type'] == 'SS']

        for time_in in times:
            image_file = self.get_all_sdo_images(time_in, save_path=save_dir)
            ch_mask = self.gen_feature_mask(time_in, [elem for elem in ch if elem['event_starttime'] == time_in], image_file)
            ar_mask = self.gen_feature_mask(time_in, [elem for elem in ar if elem['event_starttime'] == time_in], image_file)
            ss_mask = self.gen_feature_mask(time_in, [elem for elem in ss if elem['event_starttime'] == time_in], image_file)

            self.write_feature_mask(ch_mask, time_in, 'CH', save_path=save_dir)
            self.write_feature_mask(ar_mask, time_in, 'AR', save_path=save_dir)
            self.write_feature_mask(ss_mask, time_in, 'SS', save_path=save_dir)

    def gen_date_list(self):
        """
        Create list of datetimes by day

        :return: list of datetimes broken into day segments
        """
        time_between = self.tend - self.tstart
        return [self.tend - datetime.timedelta(days=dd) for dd in range(0, time_between.days + 1)]

    def gen_feature_mask(self, feature_time, feature_list, image_filepath):
        """
        Create a binary mask of feature locations

        :param feature_time: Time string of the feature detection
        :param feature_list: List of HEK dictionaries
        :param image_filepath: file path of image in which the feature is detected
        :return: binary mask of feature
        """

        #aia_wcs = WCS(self.get_header(image_filepath)[0])
        aia_map = Map(image_filepath)

        mask = np.zeros([4096,4096])

        #if aia_wcs.array_shape != mask.shape:
        #    raise ValueError('Mask and WCS array shapes do not agree.')

        # parsing boundary coordinate string for each feature
        for feature in feature_list:
            p1 = feature["hpc_boundcc"][9:-2]
            p2 = p1.split(',')
            p3 = [v.split(" ") for v in p2]
            #feature_date = parse_time(feature_time)

            feature_boundary = SkyCoord([(float(v[0]), float(v[1])) * u.arcsec for v in p3], frame=aia_map.coordinate_frame)

            # Vertices  of feature
            pixel_vertex = aia_map.world_to_pixel(feature_boundary)
            mask[np.round(pixel_vertex[0]).astype('int'), np.round(pixel_vertex[1]).astype('int')] = True

            # Need to add in contours between vertices
            # Need to add in feature filling in mask

        return mask

    def write_feature_mask(self, mask_in, feature_time, feature_type='Mask', save_path=''):
        """
        Save the feature mask as an image

        :param mask_in: Mask of feature as np array
        :param feature_time: time of the feature detection as DateTime
        :param feature_type: String of feature type e.g. "SS", "CH", or "AR"
        :param save_path: save path for saving mask as image
        :return: none
        """

        feature_date = parse_time(feature_time)

        feature_mask_name = feature_date.strftime("%Y_%m_%d__%H_%M_%S")+'__'+feature_type+'_mask.jp2'
        save_mask_name = os.path.join(save_path, feature_mask_name)

        imwrite(save_mask_name, mask_in.astype('uint8'))

    def get_header(self, filepath):
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


