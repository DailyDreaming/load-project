import os
import glob
import shutil


def link_project_metadata(project_dir, file_pattern='*.json'):
    assert os.path.isdir(project_dir)
    assert len(os.path.basename(project_dir)) == len('4a95101c-9ffc-4f30-a809-f04518a23803')
    bundle_dir = os.path.join(project_dir, 'bundle')
    hca_dir = os.path.join(project_dir, 'hca')

    if os.path.exists(hca_dir):
        if os.path.exists(bundle_dir):
            shutil.rmtree(bundle_dir)

        os.mkdir(bundle_dir)

        for filename in glob.glob(os.path.join(hca_dir, file_pattern)):
            new_filename = os.path.join(bundle_dir, os.path.basename(filename))
            os.link(filename, new_filename)
    else:
        raise RuntimeError(f'{hca_dir} does not exist.  Nothing to link over.')
