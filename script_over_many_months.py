import os
import JP2_Image_Download as Jpd
from dateutil.rrule import rrule, MONTHLY 
from sunpy.time import parse_time  


#save_dir = os.path.join(os.path.expanduser('~'), 'Data/michael/MachineLearning/Hek_project/')
save_dir = os.path.abspath('/Volumes/SolarData/LabledImages')

begin_list = [dt for dt in rrule(MONTHLY, dtstart=parse_time('2010/07/01 00:00:00') , until=parse_time('2012/05/31 23:30:00'))] 
end_list = [dt for dt in rrule(MONTHLY, dtstart=parse_time('2010/07/31 23:30:00') , until=parse_time('2012/05/31 23:30:00'))] 

for tstart, tend in zip(begin_list, end_list):

	j = Jpd.Jp2ImageDownload(save_dir, tstart=tstart.strftime('%Y/%m/%d %H:%M:%S'), tend=tend.strftime('%Y/%m/%d %H:%M:%S'))

	while j.download_flag:
		j.download_images()
		if j.missed_downloads:
			print('Missed downloads. Making new attempt(s)...')

	j.data_cleanup()

	# Test the creation of label-masks
	j.make_labels()
