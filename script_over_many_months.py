import os
import JP2_Image_Download as Jpd
from dateutil.rrule import rrule, MONTHLY
import datetime
from sunpy.time import parse_time
import logging  


# SET THESE PARAMETERS

#save_dir = os.path.join(os.path.expanduser('~'), 'Data/michael/MachineLearning/Hek_project/')
save_dir = os.path.abspath('/Volumes/SolarData/LabledImages')

start_date = '2010/06/01 00:00:00'  # inclusive
end_date = '2013/06/01 00:00:00'  # not inclusive

# SHOULDN'T NEED TO CHANGE BELOW THIS
begin_list = [dt for dt in rrule(MONTHLY, dtstart=parse_time(start_date), until=parse_time(end_date))]
end_list = [ii - datetime.timedelta(minutes=30) for ii in begin_list]

for tstart, tend in zip(begin_list[:-1], end_list[1:]):

    j = Jpd.Jp2ImageDownload(save_dir, tstart=tstart, tend=tend)
    logging.basicConfig(format='%(asctime)s %(message)s', filename=os.path.join(j.save_dir, 'logger.log'), level=logging.DEBUG)

    while j.download_flag:
        try:
            j.download_images()
            if j.missed_downloads:
                print('Missed downloads. Making new attempt(s)...')
        except ConnectionResetError:
            print('HEK server error. Trying again...')
            logging.debug('HEK server error raised ConnectionResetError')
            continue

    j.data_cleanup()

    # Test the creation of label-masks
    j.make_labels()

