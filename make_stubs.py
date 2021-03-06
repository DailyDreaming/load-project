import os
from _pathlib import Path


cwd = Path(os.getcwd())
child_links = [x for x in cwd.iterdir() if x.is_symlink()]
ids = [(os.readlink(str(link)), link) for link in child_links]
ids = sorted(ids)

for uuid, geo in ids:
    print(f'''
class {geo.name}(Converter):
    """
    {uuid}
    """
    
    def _convert(self):
        raise NotImplementedError()

''')
