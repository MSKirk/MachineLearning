import os
import JP2_Image_Download as Jpd

save_dir = os.path.join(os.path.expanduser('~'), 'Dev/michael/JP2000')

tstart = '2012/06/01 00:00:00'
tend = '2012/06/02 23:59:59'
j = Jpd.Jp2ImageDownload(save_dir, tstart=tstart, tend=tend)
j.download_images()
# Test the creation of label-masks
j.make_labels()