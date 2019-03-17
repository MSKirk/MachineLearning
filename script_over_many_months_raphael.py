import os
import JP2_Image_Download as Jpd
from dateutil.rrule import rrule, MONTHLY
from sunpy.time import parse_time
import logging


save_dir = os.path.abspath('/Volumes/RAPH_1TB/Data/Michael/Hek_project')

begin_list = [dt for dt in rrule(MONTHLY, dtstart=parse_time('2016/06/01 00:00:00') , until=parse_time('2019/02/01 00:00:00'))]
end_list = [dt for dt in rrule(MONTHLY, dtstart=parse_time('2016/06/30 23:30:00') , until=parse_time('2019/02/28 23:30:00'))]

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
