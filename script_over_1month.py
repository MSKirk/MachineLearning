import os
import JP2_Image_Download as Jpd


save_dir = os.path.join(os.path.expanduser('~'), 'Data/michael/MachineLearning/Hek_project/')

tstart = '2012/07/01 00:00:00'
tend = '2012/07/31 23:00:00'
j = Jpd.Jp2ImageDownload(save_dir, tstart=tstart, tend=tend)

while j.download_flag:
    j.download_images()
    if j.missed_downloads:
        print('Missed downloads. Making new attempt(s)...')

j.data_cleanup()

# Test the creation of label-masks
j.make_labels()


