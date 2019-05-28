import os
import JP2_Image_Download as Jpd
import logging

if __name__ == '__main__':

    save_dir = os.path.abspath('/Users/rattie/Data/Michael/MachineLearning/Hek_project/test_few_days')

    tstart = '2013/02/01 00:00:00'
    tend = '2013/02/01 23:30:00'
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


