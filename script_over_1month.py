import os
import JP2_Image_Download as Jpd
import logging

if __name__ == '__main__':

    save_dir = os.path.abspath('/Volumes/RAPH_1TB/Data/Michael/Hek_project')

    tstart = '2017/04/01 00:00:00'
    tend = '2017/04/30 23:30:00'
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
    j.make_labels()


