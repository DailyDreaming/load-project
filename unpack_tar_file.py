from pathlib import Path
import shutil
from tarfile import TarFile


def main():
    projects_dir = Path.cwd() / 'projects'
    for project_dir in projects_dir.iterdir():
        if project_dir.is_dir() and not project_dir.is_symlink():
            geo_dir = project_dir / 'geo'
            for geo_file in geo_dir.iterdir():
                if geo_file.name.endswith('_RAW.tar'):
                    dest_dir = geo_dir / geo_file.stem
                    completion_file = dest_dir / '.complete'
                    if completion_file.exists():
                        print('Expansion of', dest_dir, 'already complete')
                    else:
                        with TarFile.open(geo_file.as_posix()) as tar_file:
                            assert completion_file.name not in tar_file.getnames()
                            if dest_dir.exists():
                                print(f'Removing partially expanded', dest_dir)
                                shutil.rmtree(dest_dir.as_posix())
                            print('Expanding', dest_dir)
                            dest_dir.mkdir()
                            tar_file.extractall(dest_dir.as_posix())
                            completion_file.touch()
                            print('Expansion of', dest_dir, 'is complete')



if __name__ == '__main__':
    main()
