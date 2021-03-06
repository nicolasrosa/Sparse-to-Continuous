# ========
#  README
# ========
# Apollo Scape
# Uses Depth Maps: measures distances [close - LOW values, far - HIGH values]
# Image: (2710, 3384, 3) uint8
# Depth: (2710, 3384)    uint16

# -----
# Dataset Guidelines
# -----
# Raw Depth image to Depth (meters):
# In the depth image, the depth value is save as unsigned short int format. It can be easily read in OpenCV as:
#   cv::Mat depth_u16 = cv::imread ( depth_path, CV_LOAD_IMAGE_ANYDEPTH);
#
# The absolute depth value in meter can be obtained as:
#   double depth_value = depth_u16.at(row, col) / 200.00;
# -----


# ===========
#  Libraries
# ===========
import glob
import os

from modules.args import args
from .dataset import Dataset


# ===================
#  Class Declaration
# ===================
class Apolloscape(Dataset):
    def __init__(self, **kwargs):
        super(Apolloscape, self).__init__(**kwargs)

    def get_filenames_lists(self, mode, test_split='', test_file_path=''):
        file_path = self.get_file_path(mode, test_split, test_file_path)

        if os.path.exists(file_path):
            image_filenames, depth_filenames = self.read_text_file(file_path, self.dataset_path)
        else:
            print("[Dataloader] '%s' doesn't exist..." % file_path)
            print("[Dataloader] Searching files using glob (This may take a while)...")

            # Finds input images and labels inside the list of folders.
            image_filenames_tmp = glob.glob(self.dataset_path + "ColorImage/*/*/*.jpg")
            depth_filenames_tmp = glob.glob(self.dataset_path + "Depth/*/*/*.png")

            image_filenames_aux = [os.path.splitext(os.path.split(image)[1])[0] for image in image_filenames_tmp]
            depth_filenames_aux = [os.path.splitext(os.path.split(depth)[1])[0] for depth in depth_filenames_tmp]

            # TODO: Add Comment
            image_filenames, depth_filenames, n2, m2 = self.search_pairs(image_filenames_tmp, depth_filenames_tmp,
                                                                         image_filenames_aux, depth_filenames_aux)

            # Splits Train/Test Subsets
            divider = int(n2 * self.ratio)

            if mode == 'train':
                image_filenames = image_filenames[:divider]
                depth_filenames = depth_filenames[:divider]
            elif mode == 'test':
                # Defines Testing Subset
                image_filenames = image_filenames[divider:]
                depth_filenames = depth_filenames[divider:]

            n3, m3 = len(image_filenames), len(depth_filenames)

            print('%s_image_set: %d/%d' % (mode, n3, n2))
            print('%s_depth_set: %d/%d' % (mode, m3, m2))

            if args.debug:
                filenames = list(zip(image_filenames[:10], depth_filenames[:10]))
                for i in filenames:
                    print(i)
                input("Continue...")

            # TODO: Acredito que dê pra mover a chamada dessa função para fora
            self.save_list(image_filenames, depth_filenames, self.name, mode, self.dataset_path)

        return image_filenames, depth_filenames
