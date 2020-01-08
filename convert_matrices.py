from abc import (
    ABCMeta,
    abstractmethod,
)
import logging
import os
from pathlib import Path
import shutil
import sys
from typing import (
    Iterable,
    Optional,
)

from dataclasses import (
    astuple,
    dataclass,
)

from csv2mtx import convert_csv_to_mtx

log = logging.getLogger(__file__)


@dataclass(frozen=True)
class Matrix:
    mtx: str
    genes: str
    barcodes: str


@dataclass(frozen=True)
class CSV:
    name: str
    sep: str = ','
    rows_are_genes: bool = True


class Converter(metaclass=ABCMeta):

    def __init__(self, project_dir: Path):
        self.project_dir = project_dir

    @property
    def matrices_dir(self) -> Path:
        return self.project_dir / 'matrices'

    @property
    def geo_dir(self) -> Path:
        return self.project_dir / 'geo'

    @property
    def bundle_dir(self) -> Path:
        return self.project_dir / 'bundle'

    @property
    def zip_file(self) -> Path:
        return self.bundle_dir / 'matrix.mtx.zip'

    def matrix_dir(self, input_: str) -> Path:
        return self.matrices_dir / input_.replace('/', '__')

    def convert(self):
        if self.zip_file.exists():
            log.info('Final matrix already exists for project %s; moving on.', self.project_dir)
        else:
            self._convert()
            self._create_zip()

    def _create_zip(self):
        os.makedirs(self.zip_file.parent.as_posix(), exist_ok=True)
        atomic_make_archive(self.zip_file, root_dir=self.matrices_dir)

    @abstractmethod
    def _convert(self):
        raise NotImplementedError()

    std_matrix = Matrix(
        mtx='matrix.mtx.gz',
        genes='genes.tsv.gz',
        barcodes='barcodes.tsv.gz'
    )

    def _link_matrix(self, src: 'Matrix'):
        assert all(name.endswith('.gz') for name in astuple(src))
        dst_dir = self.matrix_dir(src.mtx)
        dst_dir.mkdir(parents=True, exist_ok=True)
        dst = self.std_matrix
        for src_name, dst_name in zip(astuple(src), astuple(dst)):
            idempotent_link(self.geo_dir / src_name, dst_dir / dst_name)

    def _link_matrices(self, matrices):
        for matrix in matrices:
            self._link_matrix(matrix)

    def _convert_csvs(self, csvs: Iterable[CSV]):
        for csv in csvs:
            convert_csv_to_mtx(input_file=self.geo_dir / csv.name,
                               output_dir=self.matrix_dir(csv.name),
                               delimiter=csv.sep,
                               rows_are_genes=csv.rows_are_genes)


def atomic_make_archive(dst: Path, root_dir: Path):
    tmp_stem = dst.parent / (dst.stem + '.tmp')
    tmp = tmp_stem.parent / (tmp_stem.name + '.zip')
    try:
        # make_archive adds it's own .zip at the end
        shutil.make_archive(tmp_stem, 'zip', root_dir=root_dir.as_posix())
    except BaseException as e:
        tmp.unlink()
        raise e
    else:
        tmp.rename(dst)


def inode(file: Path, missing_ok=False) -> Optional[int]:
    try:
        return file.stat().st_ino
    except FileNotFoundError:
        if missing_ok:
            return None
        else:
            raise


def idempotent_link(src: Path, dst: Path):
    if inode(src) != inode(dst, missing_ok=True):
        try:
            dst.unlink()
        except FileNotFoundError:
            pass
        os.link(src.as_posix(), dst.as_posix())


class GSE107909(Converter):
    """
    04ba7269-1301-5758-8f13-025565326f66
    """

    def _convert(self):
        self._convert_csvs([
            CSV("GSE107909_RAW/GSM2883183_PLNr9c.csv.gz"),
            CSV("GSE107909_RAW/GSM2883182_PLN++.csv.gz")
        ])


class GSE117089(Converter):
    """
    061ec9d5-9acf-54db-9eee-555136d5ce41
    """

    def _convert(self):
        self._link_matrices([
            Matrix(
                barcodes='GSE117089_RAW/GSM3271042_RNA_only_A549_cell.txt.gz',
                genes='GSE117089_RAW/GSM3271042_RNA_only_A549_gene.txt.gz',
                mtx='GSE117089_RAW/GSM3271042_RNA_only_A549_gene_count.txt.gz',
            ),
            Matrix(
                barcodes='GSE117089_RAW/GSM3271043_ATAC_only_A549_cell.txt.gz',
                genes='GSE117089_RAW/GSM3271043_ATAC_only_A549_peak.txt.gz',
                mtx='GSE117089_RAW/GSM3271043_ATAC_only_A549_peak_count.txt.gz',
            ),
            Matrix(
                barcodes='GSE117089_RAW/GSM3271044_RNA_mouse_kidney_cell.txt.gz',
                genes='GSE117089_RAW/GSM3271044_RNA_mouse_kidney_gene.txt.gz',
                mtx='GSE117089_RAW/GSM3271044_RNA_mouse_kidney_gene_count.txt.gz',
            ),
            Matrix(
                barcodes='GSE117089_RAW/GSM3271045_ATAC_mouse_kidney_cell.txt.gz',
                genes='GSE117089_RAW/GSM3271045_ATAC_mouse_kidney_peak.txt.gz',
                mtx='GSE117089_RAW/GSM3271045_ATAC_mouse_kidney_peak_count.txt.gz',
            )
        ])


class GSE114557(Converter):
    """
    06917f50-92aa-5e58-8376-aae1d888e8b7
    """

    def _convert(self):
        raise NotImplementedError()


class GSE131736(Converter):
    """
    069198f7-c2b5-5d39-988c-1cb12db4f28a
    """

    def _convert(self):
        raise NotImplementedError()


class GSE67835(Converter):
    """
    06a318d9-54d8-5e41-aab5-f2d682fba690
    """

    def _convert(self):
        raise NotImplementedError()


class GSE102580(Converter):
    """
    06f8848d-9c54-5829-92d3-d334809ad1e2
    """

    def _convert(self):
        raise NotImplementedError()


class GSE107585(Converter):
    """
    096b7311-2bf7-5e61-9afb-d65c24a71243
    """

    def _convert(self):
        raise NotImplementedError()


class GSE106273(Converter):
    """
    099c02da-23b2-5748-8618-92bc6770dc51
    """

    def _convert(self):
        self._link_matrices([
            Matrix(
                mtx='GSE106273_combined_matrix.tsv.gz',
                genes='GSE106273_combined_genes.tsv.gz',
                barcodes='GSE106273_combined_barcodes.tsv.gz',
            ),
        ])


class GSE130430(Converter):
    """
    0a8f2289-5862-5bf0-8c27-0885453de788
    """

    def _convert(self):
        raise NotImplementedError()


class GSE86469(Converter):
    """
    0f4d7e06-5f77-5614-8cd6-123f555dc9b1
    """

    def _convert(self):
        raise NotImplementedError()


class GSE129798(Converter):
    """
    11d48902-5824-5520-836b-7fc78ed02a61
    """

    def _convert(self):
        raise NotImplementedError()


class GSE126836(Converter):
    """
    1a72144b-f4e8-5fd5-b46e-6a40eee8b6e6
    """

    def _convert(self):
        raise NotImplementedError()


class GSE81608(Converter):
    """
    1a7ccf4f-a500-5aa6-b31a-2fff80cf8f08
    """

    def _convert(self):
        raise NotImplementedError()


class GSE97104(Converter):
    """
    1a85f2da-5aaa-5bc2-a6ea-9584742118e6
    """

    def _convert(self):
        raise NotImplementedError()


class GSE113197(Converter):
    """
    1cafb09c-e0dc-536b-b166-5cb6debfc3cf
    """

    def _convert(self):
        raise NotImplementedError()


class GSE110499(Converter):
    """
    227f5c51-389c-576c-b4d3-e4da53b89f79
    """

    def _convert(self):
        raise NotImplementedError()


class GSE36552(Converter):
    """
    2b4c411d-35b1-5f97-9d6a-c9331c7f679a
    """

    def _convert(self):
        raise NotImplementedError()


class GSE132044(Converter):
    """
    2cdd0744-6422-57fd-8fdd-9ac2bb8bf257
    """

    def _convert(self):
        raise NotImplementedError()


class GSE130636(Converter):
    """
    2ff83170-ec52-5b24-9962-161c558f52ba
    """

    def _convert(self):
        raise NotImplementedError()


class GSE81383(Converter):
    """
    30fb622d-6629-527e-a681-6e2ba143af3d
    """

    def _convert(self):
        raise NotImplementedError()


class GSE116237(Converter):
    """
    31d48835-7d9f-52ae-8cdc-ae227b63dd2c
    """

    def _convert(self):
        raise NotImplementedError()


class GSE114802(Converter):
    """
    36a7d62a-57ae-59af-ae8d-7dd2a8f1422e
    """

    def _convert(self):
        raise NotImplementedError()


class GSE124472(Converter):
    """
    389ad9f9-4a14-5a3d-b971-45dc3baf95f1
    """

    def _convert(self):
        raise NotImplementedError()


class GSE84465(Converter):
    """
    39bfc05a-44ca-507a-bbf5-156bd35c5c74
    """

    def _convert(self):
        raise NotImplementedError()


class GSE134881(Converter):
    """
    3fe16b18-e782-542b-b308-de9b26e7f69c
    """

    def _convert(self):
        raise NotImplementedError()


class GSE128639(Converter):
    """
    427157c0-993f-56ae-9e40-8cfe40ef81c5
    """

    def _convert(self):
        raise NotImplementedError()


class GSE118127(Converter):
    """
    458cbaeb-c4b7-5537-b1f3-a5d537478112
    """

    def _convert(self):
        raise NotImplementedError()


class GSE81905(Converter):
    """
    4f5c0011-416d-5e8e-8eb6-f7cb5b0140a5
    """

    def _convert(self):
        raise NotImplementedError()


class GSE94820(Converter):
    """
    5016f45e-b86a-57ce-984e-a50605641d08
    """

    def _convert(self):
        raise NotImplementedError()


class GSE81904(Converter):
    """
    51a21599-a014-5c5a-9760-d5bdeb80f741
    """

    def _convert(self):
        raise NotImplementedError()


class GSE116470(Converter):
    """
    53fe8e51-1646-5f68-96ac-e4d4fde67b93
    """

    def _convert(self):
        raise NotImplementedError()


class GSE124494(Converter):
    """
    56483fc6-ab20-5495-bb93-8cd2ce8a322a
    """

    def _convert(self):
        raise NotImplementedError()


class GSE135889(Converter):
    """
    56d9146d-bc73-5327-9615-05931f1863f6
    """

    def _convert(self):
        raise NotImplementedError()


class GSE111727(Converter):
    """
    57c9a7c8-5dc5-551b-b4af-7d37d8a87f64
    """

    def _convert(self):
        raise NotImplementedError()


class GSE84147(Converter):
    """
    5cd871a3-96ab-52ed-a7c9-77e91278c13d
    """

    def _convert(self):
        raise NotImplementedError()


class GSE93374(Converter):
    """
    60ec348b-ff28-5d47-b0d6-b787f1885c9c
    """

    def _convert(self):
        raise NotImplementedError()


class GSE127969(Converter):
    """
    62437ea1-3d06-5f22-b9de-a7d934138dd5
    """

    def _convert(self):
        raise NotImplementedError()


class GSE75478(Converter):
    """
    682f2474-f875-5e0a-bf99-2b102c8c6193
    """

    def _convert(self):
        raise NotImplementedError()


class GSE100618(Converter):
    """
    6ae3cbfe-200e-5c03-a74a-edd266b8182b
    """

    def _convert(self):
        raise NotImplementedError()


class GSE75367(Converter):
    """
    6b892786-989c-5844-b191-6d420e328fdf
    """

    def _convert(self):
        raise NotImplementedError()


class GSE70580(Converter):
    """
    6fb6d88c-7023-53fb-967b-ef95b2f6f5a0
    """

    def _convert(self):
        raise NotImplementedError()


class GSE130606(Converter):
    """
    73fd591e-2310-5983-ba8a-8079c0d0b758
    """

    def _convert(self):
        raise NotImplementedError()


class GSE75688(Converter):
    """
    789850ec-3540-5023-9767-fb8a4d2a21fc
    """

    def _convert(self):
        raise NotImplementedError()


class GSE89232(Converter):
    """
    7acdc227-c543-5b0c-8bd8-c6fa4e30310a
    """

    def _convert(self):
        raise NotImplementedError()


class GSE107618(Converter):
    """
    7dd4e06c-c889-511c-a4d4-45b74088caa8
    """

    def _convert(self):
        raise NotImplementedError()


class GSE132802(Converter):
    """
    7eedae3a-b350-5a3a-ab2d-8dcebc4a37b2
    """

    def _convert(self):
        raise NotImplementedError()


class GSE75140(Converter):
    """
    80ad934f-66ed-5c21-8f9a-7d3b0f58bcab
    """

    def _convert(self):
        raise NotImplementedError()


class GSE130473(Converter):
    """
    86963b4f-1e8e-5691-9ba3-465f3a789428
    """

    def _convert(self):
        raise NotImplementedError()


class GSE96583(Converter):
    """
    88564eae-cceb-5eee-957d-5f0a251fb177
    """

    def _convert(self):
        raise NotImplementedError()


class GSE90806(Converter):
    """
    8d1bf054-faad-5ee8-a67e-f9b8f379e6c3
    """

    def _convert(self):
        raise NotImplementedError()


class GSE76312(Converter):
    """
    932ae148-c3c2-5b07-91c0-2083cafe0dc1
    """

    def _convert(self):
        raise NotImplementedError()


class GSE93593(Converter):
    """
    99ab18ff-ca15-5d7d-9c2d-5c0b537fb1c2
    """

    def _convert(self):
        raise NotImplementedError()


class GSE92280(Converter):
    """
    9d65c4d0-c048-5c4f-8278-85dac99ea2ae
    """

    def _convert(self):
        raise NotImplementedError()


class GSE103354(Converter):
    """
    9fc2d285-804a-5989-956f-1843a0f11673
    """

    def _convert(self):
        raise NotImplementedError()


class GSE102596(Converter):
    """
    a0b6322d-1da3-5481-8768-84227ad4dd1e
    """

    def _convert(self):
        raise NotImplementedError()


class GSE44183(Converter):
    """
    aa372dea-8469-5b80-9007-18c16a21655d
    """

    def _convert(self):
        raise NotImplementedError()


class GSE103275(Converter):
    """
    b0f40b69-943f-5959-9457-c8e53c2d480e
    """

    def _convert(self):
        raise NotImplementedError()


class GSE110154(Converter):
    """
    b26137d3-a709-5492-aa74-0d783e6b628b
    """

    def _convert(self):
        raise NotImplementedError()


class GSE86473(Converter):
    """
    b48f6f16-1b5a-5055-9e14-a8920e1bcaad
    """

    def _convert(self):
        raise NotImplementedError()


class GSE114374(Converter):
    """
    b4b128d5-61e5-510e-9d91-a151b94fbb99
    """

    def _convert(self):
        raise NotImplementedError()


class GSE89322(Converter):
    """
    b55c0638-d86b-5665-9ad1-0d45b937770a
    """

    def _convert(self):
        raise NotImplementedError()


class GSE86146(Converter):
    """
    b5a0936b-a351-54ac-8d7d-0af6926e0bdc
    """

    def _convert(self):
        raise NotImplementedError()


class GSE99795(Converter):
    """
    ba75d697-8712-5e75-b35b-fd3f2b66cae5
    """

    def _convert(self):
        raise NotImplementedError()


class GSE81547(Converter):
    """
    bd4ebaac-7bcb-5069-b3ea-3a13887092e8
    """

    def _convert(self):
        raise NotImplementedError()


class GSE115469(Converter):
    """
    bdfe2399-b8a6-5b6a-9f0a-a5fd81d08ff4
    """

    def _convert(self):
        raise NotImplementedError()


class GSE111586(Converter):
    """
    bfddbefc-f2fd-5815-89a9-a94ed667be82
    """

    def _convert(self):
        raise NotImplementedError()


class GSE132040(Converter):
    """
    c2e2302f-4077-5394-9ee2-78a0ec94cbb7
    """

    def _convert(self):
        raise NotImplementedError()


class GSE84133(Converter):
    """
    c366f1f5-27aa-5157-a142-110e492a3e52
    """

    def _convert(self):
        raise NotImplementedError()


class GSE75659(Converter):
    """
    cac7f9f2-0592-5617-9530-f63803c49f8b
    """

    def _convert(self):
        raise NotImplementedError()


class GSE109822(Converter):
    """
    cc2112b7-9df1-5910-a7c6-6e41203130fa
    """

    def _convert(self):
        raise NotImplementedError()


class GSE109979(Converter):
    """
    d136eea6-03f9-5f02-86c7-c677b4c80164
    """

    def _convert(self):
        raise NotImplementedError()


class GSE131181(Converter):
    """
    d26d2ae7-4355-5ac1-8476-2e514973097e
    """

    def _convert(self):
        raise NotImplementedError()


class GSE107746(Converter):
    """
    d36952f4-cfa7-5d03-b4b6-db2c31dd41c6
    """

    def _convert(self):
        raise NotImplementedError()


class GSE131685(Converter):
    """
    d9117a4f-36e0-5912-b8cd-744a0c5306c7
    """

    def _convert(self):
        raise NotImplementedError()


class GSE108041(Converter):
    """
    d9258dc7-985e-533a-9fc2-3ad9cc7e32ca
    """

    def _convert(self):
        raise NotImplementedError()


class GSE106540(Converter):
    """
    dd761426-cc9c-5fbd-8bec-93a4cc4eb999
    """

    def _convert(self):
        raise NotImplementedError()


class GSE132566(Converter):
    """
    df1875a9-1a6a-58e0-8fd2-dac0ebb3b1b2
    """

    def _convert(self):
        raise NotImplementedError()


class GSE83139(Converter):
    """
    e8579e71-7472-5671-85b0-9841a4d06d5a
    """

    def _convert(self):
        raise NotImplementedError()


class GSE76381(Converter):
    """
    ebb8c1be-6739-57b8-9ce3-aa67caa900b4
    """

    def _convert(self):
        raise NotImplementedError()


class GSE117498(Converter):
    """
    ed008b9b-0039-5ec2-a557-f082d4ba1810
    """

    def _convert(self):
        raise NotImplementedError()


class GSE114396(Converter):
    """
    ef6f570f-a991-5528-9649-fbf06e6eb896
    """

    def _convert(self):
        raise NotImplementedError()


class GSE109488(Converter):
    """
    efddbf1e-a5cf-5e61-a38f-6a9f228e07c2
    """

    def _convert(self):
        raise NotImplementedError()


class GSE57872(Converter):
    """
    f0caef8c-f839-539d-aa17-61fe04e6d3dd
    """

    def _convert(self):
        raise NotImplementedError()


class GSE108291(Converter):
    """
    f10b6bef-febb-58dd-83ee-1180d076e53f
    """

    def _convert(self):
        raise NotImplementedError()


class GSE73727(Converter):
    """
    ff0a4d85-a1c7-571c-97ae-d964eee7ecad
    """

    def _convert(self):
        raise NotImplementedError()


def main(projects: Path):
    failed_projects = {}
    succeeded_projects = set()
    for project_dir in projects.iterdir():
        if project_dir.is_symlink():
            try:
                converter_class = globals()[project_dir.name]
                converter = converter_class(project_dir)
                converter.convert()
            except Exception as e:
                failed_projects[project_dir] = e
                log.exception('Failed to process project', exc_info=True)
            else:
                succeeded_projects.add(project_dir)

    print('Failed projects', file=sys.stderr)
    for p in failed_projects:
        print(p, file=sys.stderr)

    print('\nSucceeded projects')
    for project_dir in succeeded_projects:
        print(project_dir)


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',
                        level=logging.DEBUG)
    main(Path.cwd() / 'projects')
