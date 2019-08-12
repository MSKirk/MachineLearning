import os
import JP2_Image_Download as Jpd
from dateutil.rrule import rrule, MONTHLY
import datetime
from sunpy.time import parse_time
import logging  


if __name__ == '__main__':

    # SET THESE PARAMETERS:
    save_dir = os.path.abspath('/Users/mskirk/Desktop/MLDataTest')
    #save_dir = os.path.abspath('/Volumes/RAPH_1TB/Data/Michael/Hek_project')

    # DO NOT DOWNLOAD ANYTHING BEFORE 2010/12 (no HMI data in helioviewer)
    start_date = '2017/01/01 00:00:00'  # inclusive
    end_date = '2017/06/01 00:00:00'  # not inclusive

    # SHOULDN'T NEED TO CHANGE ANYTHING BELOW THIS --------->>

    begin_list = [dt for dt in rrule(MONTHLY, dtstart=parse_time(start_date).to_datetime(), until=parse_time(end_date).to_datetime())]
    end_list = [elem - datetime.timedelta(minutes=30) for elem in begin_list[1:]]
    del begin_list[-1]


    for tstart, tend in zip(begin_list, end_list):

        j = Jpd.Jp2ImageDownload(save_dir, tstart=tstart, tend=tend)

        logging.basicConfig(format='%(asctime)s %(message)s', filename=os.path.join(j.save_dir, 'logger.log'), level=logging.INFO)

        while j.missed_downloads_flag:
            try:
                j.download_images()
                if j.missed_downloads:
                    print('Missed downloads. Making new attempt(s)...')
            except ConnectionResetError:
                print('HEK server error. Trying again...')
                logging.warning('HEK server error raised ConnectionResetError')
                continue

        j.data_cleanup()

        # Test the creation of label-masks

        while j.missed_labels_flag:
            try:
                j.make_labels()
            except ConnectionResetError:
                print('HEK server error during make_labels(). Trying again...')
                logging.warning('HEK server error raised ConnectionResetError during make_labels()')
                continue
