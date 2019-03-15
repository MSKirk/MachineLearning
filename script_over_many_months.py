import os
import JP2_Image_Download as Jpd
from dateutil.rrule import rrule, MONTHLY 
from sunpy.time import parse_time
import logging  


#save_dir = os.path.join(os.path.expanduser('~'), 'Data/michael/MachineLearning/Hek_project/')
save_dir = os.path.abspath('/Volumes/SolarData/LabledImages')

begin_list = [dt for dt in rrule(MONTHLY, dtstart=parse_time('2010/07/01 00:00:00') , until=parse_time('2013/05/31 23:30:00'))] 
end_list = [dt for dt in rrule(MONTHLY, dtstart=parse_time('2010/07/31 23:30:00') , until=parse_time('2013/05/31 23:30:00'))] 

for tstart, tend in zip(begin_list, end_list):

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

