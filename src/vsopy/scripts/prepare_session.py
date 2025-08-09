import argparse
import astropy.units as u
import sys
from astropy.table import vstack
from vsopy import util, data

import faulthandler
faulthandler.enable()

def parse_args():
    parser = argparse.ArgumentParser(
        description='Prepare session image directory for measurement'
    )
    parser.add_argument('-O', '--object', type=str, required=True, help='Object name')
    parser.add_argument('-t', '--tag', type=str, required=True, help='Tag (date)')
    parser.add_argument('-w', '--work-dir', type=str, required=True, help='Work directory')
    parser.add_argument('-i', '--image-dir', type=str, required=True,
                        default=None, help='Image directory')
    parser.add_argument('-a', '--aperture', type=float, nargs=3,
                        default=[5.0, 10.0, 15.0],
                        help='Aperture radii in arcsec (default: 5.0, 10.0, 15.0)')
    parser.add_argument('-F', '--fov', type=float, default=60.0,
                        help='Field of view in arcmin (default: 60.0)')
    parser.add_argument('--overwrite', action='store_true',
                        default=False, help='Overwrite output files')

    return parser.parse_args()

def main():
    args = parse_args()

    session = util.Session(tag=args.tag, name=args.object)
    work_layout = util.WorkLayout(args.work_dir)
    session_layout = work_layout.get_session(session)
    img_layout = util.ImageLayout(args.image_dir)

    settings = util.Settings(session_layout.settings_file_path)
    if (args.overwrite):
        bands = [str(d.name) for d in img_layout.get_images(session).lights_dir.iterdir() if d.is_dir()]
        settings.set_bands(bands)
        settings.set_aperture(util.Aperture(args.aperture[0],
                                            args.aperture[1],
                                            args.aperture[2]))
        settings.save()

    images = util.session_image_list(img_layout.get_images(session))
    images.meta.update({'object': args.object})
    images.write(session_layout.images_file_path,
                format='ascii.ecsv', overwrite=args.overwrite)

    batches, batch_images = util.batch_session_images(session_layout.images_file_path)
    batches.meta.update({'object': args.object})
    batches.write(session_layout.batches_file_path,
                  format='ascii.ecsv', overwrite=args.overwrite)
    batch_images.write(session_layout.batch_images_file_path,
                       format='ascii.ecsv', overwrite=args.overwrite)

    star_data = data.StarData(work_layout.charts_dir)
    centroids, sequence = star_data.collect_stars(args.object, args.fov * u.arcmin)
    centroids.write(session_layout.centroid_file_path,
                   format='ascii.ecsv', overwrite=args.overwrite)
    sequence.write(session_layout.sequence_file_path,
                  format='ascii.ecsv', overwrite=args.overwrite)

    return 0

# Example: python3 list_images.py -O RR_Lyr -t 20230704 -w /home/user/work -i /home/user/img

if __name__ == '__main__':
    sys.exit(main())