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
    List,
    Optional,
)

from dataclasses import (
    astuple,
    dataclass,
)

from csv2mtx import (
    RowFilter,
    convert_csv_to_mtx,
)

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
    row_filter: Optional[RowFilter] = None


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

    def _link_matrices(self, *matrices: Matrix):
        for matrix in matrices:
            assert all(name.endswith('.gz') for name in astuple(matrix))
            dst_dir = self.matrix_dir(matrix.mtx)
            dst_dir.mkdir(parents=True, exist_ok=True)
            dst = self.std_matrix
            for src_name, dst_name in zip(astuple(matrix), astuple(dst)):
                idempotent_link(self.geo_dir / src_name, dst_dir / dst_name)

    def _convert_csvs(self, *csvs: CSV):
        for csv in csvs:
            convert_csv_to_mtx(input_file=self.geo_dir / csv.name,
                               output_dir=self.matrix_dir(csv.name),
                               delimiter=csv.sep,
                               rows_are_genes=csv.rows_are_genes,
                               row_filter=csv.row_filter)


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


class PostponedImplementationError(NotImplementedError):
    """This project has been examined, but postponed for some reason"""


class HDF5ConversionError(PostponedImplementationError):
    """A placeholder until we have HDF5 conversion integrated"""


class GSE107909(Converter):
    """
    04ba7269-1301-5758-8f13-025565326f66
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE107909_RAW/GSM2883183_PLNr9c.csv.gz'),
            CSV('GSE107909_RAW/GSM2883182_PLN++.csv.gz')
        )


class GSE117089(Converter):
    """
    061ec9d5-9acf-54db-9eee-555136d5ce41
    """

    def _convert(self):
        self._link_matrices(
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
        )


class GSE114557(Converter):
    """
    06917f50-92aa-5e58-8376-aae1d888e8b7
    """

    def _convert(self):
        raise PostponedImplementationError('No recognizable CSVs / MTXs.')


class GSE131736(Converter):
    """
    069198f7-c2b5-5d39-988c-1cb12db4f28a
    """

    def _convert(self):
        raise HDF5ConversionError()


class GSE67835(Converter):
    """
    06a318d9-54d8-5e41-aab5-f2d682fba690
    """

    def _convert(self):
        raise PostponedImplementationError('https://github.com/DailyDreaming/load-project/issues/43')


class GSE102580(Converter):
    """
    06f8848d-9c54-5829-92d3-d334809ad1e2
    """

    def _convert(self):
        self._convert_csvs(*[
            CSV(name=f, sep='\t', row_filter=self._filter) for f in (
                'GSE102580_filtered_normalized_counts_human.tsv.gz',
                'GSE102580_filtered_normalized_counts_human_viral_transduction.tsv.gz',
                'GSE102580_filtered_normalized_counts_mouse.tsv.gz',
            )])

    def _filter(self, row: List[str]):
        # Skip lines that start with #
        if row[0].startswith('#'):
            return True


class GSE107585(Converter):
    """
    096b7311-2bf7-5e61-9afb-d65c24a71243
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE107585_Mouse_kidney_single_cell_datamatrix.txt.gz', sep='\t'),
        )


class GSE106273(Converter):
    """
    099c02da-23b2-5748-8618-92bc6770dc51
    """

    def _convert(self):
        self._link_matrices(
            Matrix(
                mtx='GSE106273_combined_matrix.tsv.gz',
                genes='GSE106273_combined_genes.tsv.gz',
                barcodes='GSE106273_combined_barcodes.tsv.gz',
            ),
        )


class GSE130430(Converter):
    """
    0a8f2289-5862-5bf0-8c27-0885453de788
    """

    def _convert(self):
        self._link_matrices(
            Matrix(
                barcodes='GSE130430_RAW/GSM3738536_24y_F_BM_barcodes.tsv.gz',
                genes='GSE130430_RAW/GSM3738536_24y_F_BM_genes.tsv.gz',
                mtx='GSE130430_RAW/GSM3738536_24y_F_BM_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE130430_RAW/GSM3738537_25y_F_1_BM_barcodes.tsv.gz',
                genes='GSE130430_RAW/GSM3738537_25y_F_1_BM_genes.tsv.gz',
                mtx='GSE130430_RAW/GSM3738537_25y_F_1_BM_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE130430_RAW/GSM3738538_25y_F_2_BM_barcodes.tsv.gz',
                genes='GSE130430_RAW/GSM3738538_25y_F_2_BM_genes.tsv.gz',
                mtx='GSE130430_RAW/GSM3738538_25y_F_2_BM_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE130430_RAW/GSM3738539_26y_F_BM_barcodes.tsv.gz',
                genes='GSE130430_RAW/GSM3738539_26y_F_BM_genes.tsv.gz',
                mtx='GSE130430_RAW/GSM3738539_26y_F_BM_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE130430_RAW/GSM3738540_25y_M_BM_barcodes.tsv.gz',
                genes='GSE130430_RAW/GSM3738540_25y_M_BM_genes.tsv.gz',
                mtx='GSE130430_RAW/GSM3738540_25y_M_BM_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE130430_RAW/GSM3738541_35y_M_BM_barcodes.tsv.gz',
                genes='GSE130430_RAW/GSM3738541_35y_M_BM_genes.tsv.gz',
                mtx='GSE130430_RAW/GSM3738541_35y_M_BM_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE130430_RAW/GSM3738542_25y_F_2_PB_barcodes.tsv.gz',
                genes='GSE130430_RAW/GSM3738542_25y_F_2_PB_genes.tsv.gz',
                mtx='GSE130430_RAW/GSM3738542_25y_F_2_PB_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE130430_RAW/GSM3738543_26y_F_PB_barcodes.tsv.gz',
                genes='GSE130430_RAW/GSM3738543_26y_F_PB_genes.tsv.gz',
                mtx='GSE130430_RAW/GSM3738543_26y_F_PB_matrix.mtx.gz',
            ),
        )


class GSE86469(Converter):
    """
    0f4d7e06-5f77-5614-8cd6-123f555dc9b1
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE86469_GEO.islet.single.cell.processed.data.RSEM.raw.expected.counts.csv.gz'),
        )


class GSE129798(Converter):
    """
    11d48902-5824-5520-836b-7fc78ed02a61
    """

    def _convert(self):
        raise PostponedImplementationError('https://github.com/DailyDreaming/load-project/issues/66')
        # There are two matrices present but the smaller one has 100% barcode
        # overlap with the larger one. The larger one is also named "final" so
        # I'm only including the larger one.
        # noinspection PyUnreachableCode
        prefix = 'GSE129798_Mouse_Adult_DGE_final/Mouse_Adult_DGE_final'
        self._link_matrices(
            Matrix(
                mtx=prefix + '/matrix.mtx',
                genes=prefix + '/genes.tsv',
                barcodes=prefix + '/barcodes.tsv'
            )
        )


class GSE126836(Converter):
    """
    1a72144b-f4e8-5fd5-b46e-6a40eee8b6e6
    """

    def _convert(self):
        # Excluding only_Mal and only_Fem since they have 100% barcode overlap
        # with region_neur_Mal and region_neur_Fem respectively
        self._link_matrices(*[
            Matrix(
                mtx=prefix + '_matrix.mtx.gz',
                genes=prefix + '_genes.csv.gz',
                barcodes=prefix + '_barcodes.csv.gz'
            )
            for prefix in (
                'GSE126836_BNST_region_neur_Fem',
                'GSE126836_BNST_region_neur_Mal',
                'GSE126836_SN_MD5534',
                'GSE126836_SN_MD5828',
                'GSE126836_SN_MD5840',
                'GSE126836_SN_MD5862',
                'GSE126836_SN_MD5893',
                'GSE126836_SN_MD6060',
                'GSE126836_SN_MD6063'
            )
        ])


class GSE81608(Converter):
    """
    1a7ccf4f-a500-5aa6-b31a-2fff80cf8f08
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE81608_human_islets_rpkm.txt.gz', sep='\t')
        )


class GSE97104(Converter):
    """
    1a85f2da-5aaa-5bc2-a6ea-9584742118e6
    """

    def _convert(self):
        raise PostponedImplementationError('Invalid CSV syntax')
        # NOTE: this file contains a comment (#) and multiple blank lines at the beginning,
        # not sure if Daniel's script handles this

        # There are other files in RAW but the comments claims that this one is
        # the result of their concatenation
        # noinspection PyUnreachableCode
        self._convert_csvs(
            CSV('GSE97104_all_umi.mtx.txt.gz', sep='\t')
        )


class GSE113197(Converter):
    """
    1cafb09c-e0dc-536b-b166-5cb6debfc3cf
    """

    def _convert(self):
        raise PostponedImplementationError('https://github.com/DailyDreaming/load-project/issues/43')


class GSE110499(Converter):
    """
    227f5c51-389c-576c-b4d3-e4da53b89f79
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE110499_GEO_processed_MM_10X_raw_UMI_count_martix.txt.gz', sep='\t'),
            CSV('GSE110499_GEO_processed_MM_raw_TPM_matrix.txt.gz', sep='\t'),
        )


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
        self._link_matrices(
            Matrix(
                mtx='GSE132044_cortex_mm10_count_matrix.mtx.gz',
                genes='GSE132044_cortex_mm10_gene.tsv.gz',
                barcodes='GSE132044_cortex_mm10_cell.tsv.gz',
            ),
            Matrix(
                mtx='GSE132044_mixture_hg19_mm10_count_matrix.mtx.gz',
                genes='GSE132044_mixture_hg19_mm10_gene.tsv.gz',
                barcodes='GSE132044_mixture_hg19_mm10_cell.tsv.gz',
            ),
            Matrix(
                mtx='GSE132044_pbmc_hg38_count_matrix.mtx.gz',
                genes='GSE132044_pbmc_hg38_gene.tsv.gz',
                barcodes='GSE132044_pbmc_hg38_cell.tsv.gz',
            ),
        )


class GSE130636(Converter):
    """
    2ff83170-ec52-5b24-9962-161c558f52ba
    """

    def _convert(self):
        self._convert_csvs(
            self._csv('GSE130636_RAW/GSM3745992_fovea_donor_1_expression.tsv.gz'),
            self._csv('GSE130636_RAW/GSM3745993_fovea_donor_2_expression.tsv.gz'),
            self._csv('GSE130636_RAW/GSM3745994_fovea_donor_3_expression.tsv.gz'),
            self._csv('GSE130636_RAW/GSM3745995_peripheral_donor_1_expression.tsv.gz'),
            self._csv('GSE130636_RAW/GSM3745996_peripheral_donor_2_expression.tsv.gz'),
            self._csv('GSE130636_RAW/GSM3745997_peripheral_donor_3_expression.tsv.gz'),
        )

    def _csv(self, name):
        return CSV(name, sep='\t', rows_are_genes=False, row_filter=self._filter)

    def _filter(self, row: List[str]):
        del row[1]


class GSE81383(Converter):
    """
    30fb622d-6629-527e-a681-6e2ba143af3d
    """

    def _convert(self):
        # There are 2 csvs with identical data, one just has quotes around the
        # values.

        self._convert_csvs(
            CSV('GSE81383_data_melanoma_scRNAseq_BT_2015-07-02.txt.gz', sep='\t')
        )


class GSE116237(Converter):
    """
    31d48835-7d9f-52ae-8cdc-ae227b63dd2c
    """

    def _convert(self):
        # This one has \r for line breaks so this will probably fail
        self._convert_csvs(
            CSV('GSE116237_scRNAseq_expressionMatrix.txt.gz')
        )


class GSE114802(Converter):
    """
    36a7d62a-57ae-59af-ae8d-7dd2a8f1422e
    """

    def _convert(self):
        self._link_matrices(*[
            Matrix(
                mtx=prefix + 'matrix.mtx.gz',
                barcodes=prefix + 'barcodes.tsv.gz',
                genes=prefix + 'genes.tsv.gz'
            )
            for prefix in (
                'GSE114802_org4_',
                'GSE114802_org_'
            )
        ])


class GSE124472(Converter):
    """
    389ad9f9-4a14-5a3d-b971-45dc3baf95f1
    """

    def _convert(self):
        raise PostponedImplementationError('https://github.com/DailyDreaming/load-project/issues/66')
        # noinspection PyUnreachableCode
        self._link_matrices([
            Matrix(
                mtx=str(prefix / 'matrix.mtx'),
                genes=str(prefix / 'genes.tsv'),
                barcodes=str(prefix / 'barcodes.tsv')
            )
            for prefix
            in [
                Path('GSE124472_RAW') / p
                for p
                in [
                    'GSM3534656_H17w_Z1_raw_counts',
                    'GSM3534657_H17w_Z2_raw_counts',
                    'GSM3534658_H15w_Z1_raw_counts',
                    'GSM3534659_H15w_Z2_raw_counts',
                    'GSM3534660_HuOrg_D16_1_raw_counts',
                    'GSM3534661_HuOrg_D16_2_raw_counts',
                    'GSM3534662_HuOrg_D28_1_raw_counts',
                    'GSM3534663_HuOrg_D28_2_raw_counts'
                ]
            ]
        ])


class GSE84465(Converter):
    """
    39bfc05a-44ca-507a-bbf5-156bd35c5c74
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE84465_GBM_All_data.csv.gz', sep=' ')
        )


class GSE134881(Converter):
    """
    3fe16b18-e782-542b-b308-de9b26e7f69c
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE134881_RISKfpkmAll.csv.gz')
        )


class GSE128639(Converter):
    """
    427157c0-993f-56ae-9e40-8cfe40ef81c5
    """

    def _convert(self):
        # Seems that the same samples are spread across 3 files with different genes in each?
        raise PostponedImplementationError('Confusing organisation of files')


class GSE118127(Converter):
    """
    458cbaeb-c4b7-5537-b1f3-a5d537478112
    """

    def _convert(self):
        raise HDF5ConversionError()


class GSE81905(Converter):
    """
    4f5c0011-416d-5e8e-8eb6-f7cb5b0140a5
    """

    def _convert(self):
        raise PostponedImplementationError('BAM files, (probably not expression data)')


class GSE94820(Converter):
    """
    5016f45e-b86a-57ce-984e-a50605641d08
    """

    def _convert(self):
        self._convert_csvs(*[
            CSV(csv, sep='\t')
            for csv in (
                'GSE94820_raw.expMatrix_DCnMono.discovery.set.submission.txt.gz',
                'GSE94820_raw.expMatrix_deeper.characterization.set.submission.txt.gz'
            )
        ])


class GSE81904(Converter):
    """
    51a21599-a014-5c5a-9760-d5bdeb80f741
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE81904_BipolarUMICounts_Cell2016.txt.gz', sep='\t')
        )


class GSE116470(Converter):
    """
    53fe8e51-1646-5f68-96ac-e4d4fde67b93
    """

    def _convert(self):
        # Might be like. An mtx with the genes and barcodes defined in comments????
        raise PostponedImplementationError('Confusing data')


class GSE124494(Converter):
    """
    56483fc6-ab20-5495-bb93-8cd2ce8a322a
    """

    def _convert(self):
        dir_ = 'GSE124494_RAW/'
        self._link_matrices(*[
            Matrix(
                mtx=dir_ + prefix + 'matrix.mtx.gz',
                genes=dir_ + prefix + 'genes.tsv.gz',
                barcodes=dir_ + prefix + 'barcodes.tsv.gz'
            )
            for prefix in [
                'GSM3535276_AXLN1_',
                'GSM3535277_AXLN2_',
                'GSM3535278_AXLN3_',
                'GSM3535279_AXLN4_',
                'GSM3535280_HNLN1_',
                'GSM3535281_HNLN2_'
            ]
        ])


class GSE135889(Converter):
    """
    56d9146d-bc73-5327-9615-05931f1863f6
    """

    def _convert(self):
        # Doesn't seem to be any data in this one
        raise NotImplementedError()


class GSE111727(Converter):
    """
    57c9a7c8-5dc5-551b-b4af-7d37d8a87f64
    """

    def _convert(self):
        # SKipped this one because I'm tired and couldn;t figure out file relations
        raise NotImplementedError()


class GSE84147(Converter):
    """
    5cd871a3-96ab-52ed-a7c9-77e91278c13d
    """

    def _convert(self):
        # Couldn;t understand file formats
        raise NotImplementedError()


class GSE93374(Converter):
    """
    60ec348b-ff28-5d47-b0d6-b787f1885c9c
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE93374_Merged_all_020816_DGE.txt.gz', sep='\t')
        )


class GSE127969(Converter):
    """
    62437ea1-3d06-5f22-b9de-a7d934138dd5
    """

    def _convert(self):
        # Fails with a UnicodeError due to fancy quotation marks (0x93 0x94)
        raise PostponedImplementationError('Not UTF-8')
        # noinspection PyUnreachableCode
        self._convert_csvs(
            CSV('GSE127969_counts_TPM_ALL.csv.gz', sep='\t')
        )


class GSE75478(Converter):
    """
    682f2474-f875-5e0a-bf99-2b102c8c6193
    """

    def _convert(self):
        self._convert_csvs(*[
            CSV(csv)
            for csv in (
                'GSE75478_transcriptomics_raw_filtered_I1.csv.gz',
                'GSE75478_transcriptomics_raw_filtered_I2.csv.gz'
            )
        ])


class GSE100618(Converter):
    """
    6ae3cbfe-200e-5c03-a74a-edd266b8182b
    """

    def _convert(self):
        self._convert_csvs(CSV('GSE100618_HTSeq_counts.txt.gz', sep=' '))


class GSE75367(Converter):
    """
    6b892786-989c-5844-b191-6d420e328fdf
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE75367_readCounts.txt.gz', sep='\t')
        )


class GSE70580(Converter):
    """
    6fb6d88c-7023-53fb-967b-ef95b2f6f5a0
    """

    def _convert(self):
        raise PostponedImplementationError('https://github.com/DailyDreaming/load-project/issues/43')


class GSE130606(Converter):
    """
    73fd591e-2310-5983-ba8a-8079c0d0b758
    """

    def _convert(self):
        self._link_matrices(
            Matrix(
                mtx='GSE130606_matrix.mtx.gz',
                genes='GSE130606_genes.tsv.gz',
                barcodes='GSE130606_barcodes.tsv.gz'
            )
        )


class GSE75688(Converter):
    """
    789850ec-3540-5023-9767-fb8a4d2a21fc
    """

    def _convert(self):
        self._convert_csvs(
            CSV(
                'GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt.gz',
                sep='\t',
                row_filter=self._filter,
            )
        )

    def _filter(self, row: List[str]):
        del row[2]  # gene_type
        del row[0]  # gene_id, we'll use gene name (row[1])


class GSE89232(Converter):
    """
    7acdc227-c543-5b0c-8bd8-c6fa4e30310a
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE89232_expMatrix.txt.gz', sep='\t', row_filter=self._filter)
        )

    def _filter(self, row: List[str]):
        if len(row) == 957:
            row.insert(0, '')  # header row is one short
        elif len(row) == 958:
            pass
        else:
            assert False, len(row)


class GSE107618(Converter):
    """
    7dd4e06c-c889-511c-a4d4-45b74088caa8
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE107618_Merge.TPM.csv.gz', row_filter=self._filter)
        )

    def _filter(self, row: List[str]):
        for i, v in enumerate(row):
            if v == '':
                row[i] = '0'


class GSE132802(Converter):
    """
    7eedae3a-b350-5a3a-ab2d-8dcebc4a37b2
    """

    def _convert(self):
        raise HDF5ConversionError()


class GSE75140(Converter):
    """
    80ad934f-66ed-5c21-8f9a-7d3b0f58bcab
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE75140_hOrg.fetal.master.data.frame.txt', rows_are_genes=False, sep='\t', row_filter=self._filter)
        )

    def _filter(self, row: List[str]):
        del row[-1]


class GSE130473(Converter):
    """
    86963b4f-1e8e-5691-9ba3-465f3a789428
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE130473_Series_count_matrix.csv.gz')
        )


class GSE96583(Converter):
    """
    88564eae-cceb-5eee-957d-5f0a251fb177
    """

    def _convert(self):
        self._link_matrices(
            Matrix(
                mtx='GSE96583_RAW/GSM2560245_A.mat.gz',
                genes='GSE96583_batch1.genes.tsv.gz',
                barcodes='GSE96583_RAW/GSM2560245_barcodes.tsv.gz'),
            Matrix(
                mtx='GSE96583_RAW/GSM2560246_B.mat.gz',
                genes='GSE96583_batch1.genes.tsv.gz',
                barcodes='GSE96583_RAW/GSM2560246_barcodes.tsv.gz'),
            Matrix(
                mtx='GSE96583_RAW/GSM2560247_C.mat.gz',
                genes='GSE96583_batch1.genes.tsv.gz',
                barcodes='GSE96583_RAW/GSM2560247_barcodes.tsv.gz'),
            Matrix(
                mtx='GSE96583_RAW/GSM2560248_2.1.mtx.gz',
                genes='GSE96583_batch2.genes.tsv.gz',
                barcodes='GSE96583_RAW/GSM2560248_barcodes.tsv.gz'),
            Matrix(
                mtx='GSE96583_RAW/GSM2560249_2.2.mtx.gz',
                genes='GSE96583_batch2.genes.tsv.gz',
                barcodes='GSE96583_RAW/GSM2560249_barcodes.tsv.gz'),
        )


class GSE90806(Converter):
    """
    8d1bf054-faad-5ee8-a67e-f9b8f379e6c3
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE90806_RIP-Cre_ARC_GeneCounts.csv.gz')
        )


class GSE76312(Converter):
    """
    932ae148-c3c2-5b07-91c0-2083cafe0dc1
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE76312_Giustacchini_Thongjuea_et.al_Nat.Med.RPKM.txt.gz', sep='\t', row_filter=self._filter)
        )

    def _filter(self, row: List[str]):
        if len(row) == 2287:
            row.insert(0, '')  # header row is one short
        elif len(row) == 2288:
            pass
        else:
            assert False, len(row)


class GSE93593(Converter):
    """
    99ab18ff-ca15-5d7d-9c2d-5c0b537fb1c2
    """

    def _convert(self):
        # ignoring GSE93593_tpm.csv.gz
        self._convert_csvs(
            CSV('GSE93593_counts.csv.gz')
        )


class GSE92280(Converter):
    """
    9d65c4d0-c048-5c4f-8278-85dac99ea2ae
    """

    def _convert(self):
        raise PostponedImplementationError('No supplementary files')


class GSE103354(Converter):
    """
    9fc2d285-804a-5989-956f-1843a0f11673
    """

    def _convert(self):
        self._link_matrices(
            Matrix(
                barcodes='GSE103354_RAW/GSM3314330_Tp0_Rep1_GFP_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314330_Tp0_Rep1_GFP_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314330_Tp0_Rep1_GFP_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314331_Tp0_Rep1_Tom_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314331_Tp0_Rep1_Tom_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314331_Tp0_Rep1_Tom_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314332_Tp0_Rep2_GFP_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314332_Tp0_Rep2_GFP_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314332_Tp0_Rep2_GFP_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314333_Tp0_Rep2_Tom_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314333_Tp0_Rep2_Tom_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314333_Tp0_Rep2_Tom_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314334_Tp0_Rep3_GFP_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314334_Tp0_Rep3_GFP_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314334_Tp0_Rep3_GFP_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314335_Tp0_Rep3_Tom_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314335_Tp0_Rep3_Tom_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314335_Tp0_Rep3_Tom_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314336_Tp30_Rep1_GFP_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314336_Tp30_Rep1_GFP_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314336_Tp30_Rep1_GFP_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314337_Tp30_Rep1_Tom_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314337_Tp30_Rep1_Tom_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314337_Tp30_Rep1_Tom_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314338_Tp30_Rep2_GFP_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314338_Tp30_Rep2_GFP_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314338_Tp30_Rep2_GFP_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314339_Tp30_Rep2_Tom_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314339_Tp30_Rep2_Tom_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314339_Tp30_Rep2_Tom_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314340_Tp30_Rep3_GFP_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314340_Tp30_Rep3_GFP_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314340_Tp30_Rep3_GFP_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314341_Tp30_Rep3_Tom_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314341_Tp30_Rep3_Tom_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314341_Tp30_Rep3_Tom_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314342_Tp60_Rep1_GFP_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314342_Tp60_Rep1_GFP_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314342_Tp60_Rep1_GFP_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314343_Tp60_Rep1_Tom_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314343_Tp60_Rep1_Tom_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314343_Tp60_Rep1_Tom_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314344_Tp60_Rep2_GFP_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314344_Tp60_Rep2_GFP_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314344_Tp60_Rep2_GFP_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314345_Tp60_Rep2_Tom_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314345_Tp60_Rep2_Tom_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314345_Tp60_Rep2_Tom_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314346_Tp60_Rep3_GFP_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314346_Tp60_Rep3_GFP_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314346_Tp60_Rep3_GFP_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314347_Tp60_Rep3_Tom_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314347_Tp60_Rep3_Tom_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314347_Tp60_Rep3_Tom_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314348_M1_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314348_M1_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314348_M1_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314349_M2_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314349_M2_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314349_M2_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314350_M3_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314350_M3_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314350_M3_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314351_M4_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314351_M4_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314351_M4_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314352_M5_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314352_M5_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314352_M5_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103354_RAW/GSM3314353_M6_barcodes.tsv.gz',
                genes='GSE103354_RAW/GSM3314353_M6_genes.tsv.gz',
                mtx='GSE103354_RAW/GSM3314353_M6_matrix.mtx.gz',
            )
        )


class GSE102596(Converter):
    """
    a0b6322d-1da3-5481-8768-84227ad4dd1e
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE102596_RAW/GSM2741551_count-table-human16w.tsv.gz', sep='\t', row_filter=self._filter)
        )

    def _filter(self, row: List[str]):
        if len(row) == 3745:
            row.insert(0, '')  # header row is one short
        elif len(row) == 3746:
            pass
        else:
            assert False, len(row)


class GSE44183(Converter):
    """
    aa372dea-8469-5b80-9007-18c16a21655d
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE44183_human_expression_mat.txt.gz', sep='\t'),
            CSV('GSE44183_mouse_expression_mat.txt.gz', sep='\t', row_filter=self._filter),
        )

    def _filter(self, row: List[str]):
        for i, v in enumerate(row):
            if v == '':
                row[i] = '0'


class GSE103275(Converter):
    """
    b0f40b69-943f-5959-9457-c8e53c2d480e
    """

    def _convert(self):
        self._link_matrices(
            Matrix(
                barcodes='GSE103275_RAW/GSM2510616_P4-barcodes.tsv.gz',
                genes='GSE103275_RAW/GSM2510616_P4-genes.tsv.gz',
                mtx='GSE103275_RAW/GSM2510616_P4-matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103275_RAW/GSM2510617_P7-barcodes.tsv.gz',
                genes='GSE103275_RAW/GSM2510617_P7-genes.tsv.gz',
                mtx='GSE103275_RAW/GSM2510617_P7-matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103275_RAW/GSM2759554_5wk-1-barcodes.tsv.gz',
                genes='GSE103275_RAW/GSM2759554_5wk-1-genes.tsv.gz',
                mtx='GSE103275_RAW/GSM2759554_5wk-1-matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103275_RAW/GSM2759555_5wk-2-barcodes.tsv.gz',
                genes='GSE103275_RAW/GSM2759555_5wk-2-genes.tsv.gz',
                mtx='GSE103275_RAW/GSM2759555_5wk-2-matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103275_RAW/GSM2759556_P7D-barcodes.tsv.gz',
                genes='GSE103275_RAW/GSM2759556_P7D-genes.tsv.gz',
                mtx='GSE103275_RAW/GSM2759556_P7D-matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE103275_RAW/GSM2759557_P7E-barcodes.tsv.gz',
                genes='GSE103275_RAW/GSM2759557_P7E-genes.tsv.gz',
                mtx='GSE103275_RAW/GSM2759557_P7E-matrix.mtx.gz',
            ),
        )


class GSE110154(Converter):
    """
    b26137d3-a709-5492-aa74-0d783e6b628b
    """

    def _convert(self):
        raise PostponedImplementationError('https://github.com/DailyDreaming/load-project/issues/43')


class GSE86473(Converter):
    """
    b48f6f16-1b5a-5055-9e14-a8920e1bcaad
    """

    def _convert(self):
        raise PostponedImplementationError('No supplementary files')


class GSE114374(Converter):
    """
    b4b128d5-61e5-510e-9d91-a151b94fbb99
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE114374_Human_HC_expression_matrix.txt.gz', sep='\t'),
            CSV('GSE114374_Human_UC_expression_matrix.txt.gz', sep='\t'),
            CSV('GSE114374_Mouse_DSS_expression_matrix.txt.gz', sep='\t'),
            CSV('GSE114374_Mouse_HC_expression_matrix.txt.gz', sep='\t'),
        )


class GSE89322(Converter):
    """
    b55c0638-d86b-5665-9ad1-0d45b937770a
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE89322_bulk_counts.txt.gz', sep='\t'),
            CSV('GSE89322_single_cell_counts.txt.gz', sep='\t'),
        )


class GSE86146(Converter):
    """
    b5a0936b-a351-54ac-8d7d-0af6926e0bdc
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE86146_RAW/GSM2295850_F_10W_embryo1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2295851_F_14W_embryo1_1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2295852_F_14W_embryo1_2_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2295853_F_14W_embryo1_3_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2295854_F_18W_embryo1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306011_F_5W_embryo1_and_2_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306012_F_7W_embryo1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306013_F_8W_embryo1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306014_F_12W_embryo1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306015_F_18W_embryo2_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306016_F_20W_embryo1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306017_F_20W_embryo2_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306018_F_23W_embryo1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306019_F_23W_embryo2_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306020_F_24W_embryo1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306021_F_24W_embryo2_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306022_F_26W_embryo1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306023_M_4W_embryo1_and_F_11W_embryo1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306024_M_9W_embryo1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306025_M_10W_embryo1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306026_M_10W_embryo2_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306027_M_12W_embryo1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306028_M_19W_embryo1_101_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306029_M_19W_embryo1_24_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306030_M_19W_embryo1_26_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306031_M_19W_embryo2_102_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306032_M_19W_embryo2_103_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306033_M_19W_embryo2_104_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306034_M_20W_embryo1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306035_M_21W_embryo1_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306036_M_21W_embryo2_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306037_M_21W_embryo3_10_2_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306038_M_21W_embryo3_10_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306039_M_21W_embryo3_17_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306040_M_25W_embryo1_101_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306041_M_25W_embryo1_102_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306042_M_25W_embryo1_103_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306043_M_25W_embryo1_104_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306044_M_25W_embryo1_105_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306045_M_25W_embryo1_106_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306046_M_25W_embryo1_107_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306047_M_25W_embryo1_24_gene_expression.txt.gz', sep='\t'),
            CSV('GSE86146_RAW/GSM2306048_M_25W_embryo1_26_gene_expression.txt.gz', sep='\t'),
        )


class GSE99795(Converter):
    """
    ba75d697-8712-5e75-b35b-fd3f2b66cae5
    """

    def _convert(self):
        raise PostponedImplementationError('No gene by cell matrices.')


class GSE81547(Converter):
    """
    bd4ebaac-7bcb-5069-b3ea-3a13887092e8
    """

    def _convert(self):
        raise PostponedImplementationError('https://github.com/DailyDreaming/load-project/issues/43')


class GSE115469(Converter):
    """
    bdfe2399-b8a6-5b6a-9f0a-a5fd81d08ff4
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE115469_Data.csv.gz')
        )


class GSE111586(Converter):
    """
    bfddbefc-f2fd-5815-89a9-a94ed667be82
    """

    def _convert(self):
        raise PostponedImplementationError('Similar, but not gene by cell matrices.')


class GSE132040(Converter):
    """
    c2e2302f-4077-5394-9ee2-78a0ec94cbb7
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE132040_190214_A00111_0269_AHH3J3DSXX_190214_A00111_0270_BHHMFWDSXX.csv.gz')
        )


class GSE84133(Converter):
    """
    c366f1f5-27aa-5157-a142-110e492a3e52
    """

    def _convert(self):
        raise PostponedImplementationError('non-numeric descriptive columns')
        # noinspection PyUnreachableCode
        self._convert_csvs(
            CSV('GSE84133_RAW/' + csv, rows_are_genes=False)
            for csv in (
                'GSM2230757_human1_umifm_counts.csv.gz',
                'GSM2230758_human2_umifm_counts.csv.gz',
                'GSM2230759_human3_umifm_counts.csv.gz',
                'GSM2230760_human4_umifm_counts.csv.gz',
                'GSM2230761_mouse1_umifm_counts.csv.gz',
                'GSM2230762_mouse2_umifm_counts.csv.gz'
            ))


class GSE75659(Converter):
    """
    cac7f9f2-0592-5617-9530-f63803c49f8b
    """

    def _convert(self):
        raise PostponedImplementationError('https://github.com/DailyDreaming/load-project/issues/43')


class GSE109822(Converter):
    """
    cc2112b7-9df1-5910-a7c6-6e41203130fa
    """

    def _convert(self):
        self._convert_csvs(*[
            CSV(csv)
            for csv in (
                'GSE109822_CD3145.csv.gz',
                'GSE109822_CD90.csv.gz',
                'GSE109822_dermis.csv.gz'
            )
        ])


class GSE109979(Converter):
    """
    d136eea6-03f9-5f02-86c7-c677b4c80164
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE109979_329Cell_RPKM.txt.gz', sep='\t', row_filter=self._filter),
        )

    def _filter(self, row: List[str]):
        # this csv has 2 header rows, ignore the header row with no first element
        if row[0] == '':
            return True


class GSE131181(Converter):
    """
    d26d2ae7-4355-5ac1-8476-2e514973097e
    """

    def _convert(self):
        # ignoring GSE131181_e10.5.meta.data.csv.gz
        # ignoring GSE131181_e13.5.meta.data.csv.gz
        self._convert_csvs(
            CSV('GSE131181_e10.5.raw.data.csv.gz', row_filter=self._filter1),
            # CSV('GSE131181_e10.5.scale.data.csv.gz', row_filter=self._filter2),  # 4gb too big?
            CSV('GSE131181_e13.5.raw.data.csv.gz', row_filter=self._filter3),
        )

    def _filter1(self, row: List[str]):
        if len(row) == 31479:
            row.insert(0, '')  # header row is one short
        elif len(row) == 31480:
            pass
        else:
            assert False, len(row)

    def _filter2(self, row: List[str]):
        if len(row) == 22685:
            row.insert(0, '')  # header row is one short
        elif len(row) == 22686:
            pass
        else:
            assert False, len(row)

    def _filter3(self, row: List[str]):
        if len(row) == 63102:
            row.insert(0, '')  # header row is one short
        elif len(row) == 63103:
            pass
        else:
            assert False, len(row)


class GSE107746(Converter):
    """
    d36952f4-cfa7-5d03-b4b6-db2c31dd41c6
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE107746_Folliculogenesis_FPKM.log2.txt.gz', sep='\t'),
        )


class GSE131685(Converter):
    """
    d9117a4f-36e0-5912-b8cd-744a0c5306c7
    """

    def _convert(self):
        self._link_matrices(
            Matrix(
                barcodes='GSE131685_RAW/GSM4145204_kidney1_barcodes.tsv.gz',
                genes='GSE131685_RAW/GSM4145204_kidney1_features.tsv.gz',
                mtx='GSE131685_RAW/GSM4145204_kidney1_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE131685_RAW/GSM4145205_kidney2_barcodes.tsv.gz',
                genes='GSE131685_RAW/GSM4145205_kidney2_features.tsv.gz',
                mtx='GSE131685_RAW/GSM4145205_kidney2_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE131685_RAW/GSM4145206_kidney3_barcodes.tsv.gz',
                genes='GSE131685_RAW/GSM4145206_kidney3_features.tsv.gz',
                mtx='GSE131685_RAW/GSM4145206_kidney3_matrix.mtx.gz',
            ),
        )


class GSE108041(Converter):
    """
    d9258dc7-985e-533a-9fc2-3ad9cc7e32ca
    """

    def _convert(self):
        # ignoring GSE108041_RAW/GSM2888370_Uninfected.h5
        # ignoring GSE108041_RAW/GSM2888371_6h.h5
        # ignoring GSE108041_RAW/GSM2888372_8h.h5
        # ignoring GSE108041_RAW/GSM2888373_8h-2.h5
        # ignoring GSE108041_RAW/GSM2888374_10h.h5
        # ignoring GSE108041_README.txt
        # ignoring GSE108041_reference_fasta.fasta.gz
        # ignoring GSE108041_reference_genes.gtf.gz
        raise HDF5ConversionError()


class GSE106540(Converter):
    """
    dd761426-cc9c-5fbd-8bec-93a4cc4eb999
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE106540_SC_raw_counts.txt.gz', sep='\t'),
        )


class GSE132566(Converter):
    """
    df1875a9-1a6a-58e0-8fd2-dac0ebb3b1b2
    """

    def _convert(self):
        raise PostponedImplementationError('https://github.com/DailyDreaming/load-project/issues/43')


class GSE83139(Converter):
    """
    e8579e71-7472-5671-85b0-9841a4d06d5a
    """

    def _convert(self):
        raise PostponedImplementationError('looks like a SAM file??')


class GSE76381(Converter):
    """
    ebb8c1be-6739-57b8-9ce3-aa67caa900b4
    """

    def _convert(self):
        raise PostponedImplementationError('.cef files with multi-line column headers')


class GSE117498(Converter):
    """
    ed008b9b-0039-5ec2-a557-f082d4ba1810
    """

    def _convert(self):
        self._convert_csvs(*[
            CSV('GSE117498_RAW/' + csv, sep='\t')
            for csv in (
                'GSM3305359_HSC.raw_counts.tsv.gz',
                'GSM3305360_MPP.raw_counts.tsv.gz',
                'GSM3305361_MLP.raw_counts.tsv.gz',
                'GSM3305362_PreBNK.raw_counts.tsv.gz',
                'GSM3305363_MEP.raw_counts.tsv.gz',
                'GSM3305364_CMP.raw_counts.tsv.gz',
                'GSM3305365_GMP.raw_counts.tsv.gz',
                'GSM3305366_LinNegCD34PosCD164Pos.raw_counts.tsv.gz',
                'GSM3305367_LinNegCD34NegCD164high.raw_counts.tsv.gz',
                'GSM3305368_LinNegCD34lowCD164high.raw_counts.tsv.gz',
                'GSM3305369_LinNegCD34NegCD164low.raw_counts.tsv.gz'
            )
        ])


class GSE114396(Converter):
    """
    ef6f570f-a991-5528-9649-fbf06e6eb896
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE114396_RAW/GSM3141011_ILC_count.csv.gz')
        )


class GSE109488(Converter):
    """
    efddbf1e-a5cf-5e61-a38f-6a9f228e07c2
    """

    def _convert(self):
        self._convert_csvs(*[
            CSV(f'GSE109488_RAW/{csv}_gene_expression_TPM.txt.gz', sep='\t')
            for csv in (
                'GSM2944263_HEK7W2C1',
                'GSM2944264_HEK7W2C2',
                'GSM2944265_HEK7W2D1',
                'GSM2944266_HEK7W2D2',
                'GSM2944267_HEK7W2P1',
                'GSM2944268_HEK7W2P2',
                'GSM2944269_HEK9WC1',
                'GSM2944270_HEK9WC2',
                'GSM2944271_HEK9WD1',
                'GSM2944272_HEK9WD2',
                'GSM2944273_HEK9WP1',
                'GSM2944274_HEK9WP2',
                'GSM2944275_HEK13WA',
                'GSM2944276_HEK13WB',
                'GSM2944277_HEK19W1C1',
                'GSM2944278_HEK19W1C2',
                'GSM2944279_HEK19W1D1',
                'GSM2944280_HEK19W1D2',
                'GSM2944281_HEK19W1P1',
                'GSM2944282_HEK19W1P2',
                'GSM2944283_HEK19W2C1',
                'GSM2944284_HEK19W2C2',
                'GSM2944285_HEK19W2D1',
                'GSM2944286_HEK19W2D2',
                'GSM2944287_HEK19W2P1',
                'GSM2944288_HEK19W2P2',
                'GSM2944289_HEK25WC1',
                'GSM2944290_HEK25WC2',
                'GSM2944291_HEK25WD1',
                'GSM2944292_HEK25WD2',
                'GSM2944293_HEK25WP1',
                'GSM2944294_HEK25WP2',
                'GSM2944295_K7W1C1',
                'GSM2944296_K7W1C2',
                'GSM2944297_K7W1D1',
                'GSM2944298_K7W1D2',
                'GSM2944299_K7W1P1',
                'GSM2944300_K7W1P2',
                'GSM2944301_K8W1A1',
                'GSM2944302_K8W1B1',
                'GSM2944303_K8W1B2',
                'GSM2944304_K10W2C1',
                'GSM2944305_K10W2C2',
                'GSM2944306_K10W2D1',
                'GSM2944307_K10W2D2',
                'GSM2944308_K10W2P1',
                'GSM2944309_K10W2P2',
                'GSM2944310_K22W2C11',
                'GSM2944311_K22W2C12',
                'GSM2944312_K22W2C21',
                'GSM2944313_K22W2C22',
                'GSM2944314_K22W2C31',
                'GSM2944315_K22W2C32',
                'GSM2944316_K22W2D11',
                'GSM2944317_K22W2D12',
                'GSM2944318_K22W2D21',
                'GSM2944319_K22W2D22',
                'GSM2944320_K22W2D31',
                'GSM2944321_K22W2D32',
                'GSM2944322_K22W2P11',
                'GSM2944323_K22W2P12',
                'GSM2944324_K22W2P21',
                'GSM2944325_K22W2P22',
                'GSM2944326_K22W2P31',
                'GSM2944327_K22W2P32',
                'GSM2944328_K24W1C',
                'GSM2944329_K24W1D',
                'GSM2944330_K24W1P',
                'GSM2944331_K801',
                'GSM2944332_K802'
            )
        ])


class GSE57872(Converter):
    """
    f0caef8c-f839-539d-aa17-61fe04e6d3dd
    """

    def _convert(self):
        self._convert_csvs(
            CSV('GSE57872_GBM_data_matrix.txt.gz', sep='\t')
        )


class GSE108291(Converter):
    """
    f10b6bef-febb-58dd-83ee-1180d076e53f
    """

    def _convert(self):
        self._link_matrices(
            Matrix(
                barcodes='GSE108291_kid_barcodes.tsv.gz',
                genes='GSE108291_kid_genes.tsv.gz',
                mtx='GSE108291_kid_matrix.mtx.gz',
            ),
            Matrix(
                barcodes='GSE108291_org4_barcodes.tsv.gz',
                genes='GSE108291_org4_genes.tsv.gz',
                mtx='GSE108291_org4_matrix.mtx.gz',
            ),
        )
        self._convert_csvs(
            CSV('GSE108291_org_counts.csv.gz')
        )


class GSE73727(Converter):
    """
    ff0a4d85-a1c7-571c-97ae-d964eee7ecad
    """

    def _convert(self):
        raise PostponedImplementationError('No recognizable matrices.')


def main(projects: Path):
    not_implemented_projects = []
    failed_projects = []
    succeeded_projects = []
    for project_dir in sorted(projects.iterdir()):
        if project_dir.is_symlink():
            # noinspection PyBroadException
            try:
                converter_class = globals()[project_dir.name]
                converter = converter_class(project_dir)
                converter.convert()
            except NotImplementedError:
                not_implemented_projects.append(project_dir)
            except Exception:
                failed_projects.append(project_dir)
                log.exception('Failed to process project', exc_info=True)
            else:
                succeeded_projects.append(project_dir)

    print_projects('not implemented', not_implemented_projects, file=sys.stderr)
    print_projects('failed', failed_projects, file=sys.stderr)
    print_projects('succeeded', succeeded_projects)


def print_projects(title, projects, file=None):
    if len(projects) > 0:
        print('Projects', title, file=file)
        for p in projects:
            print(p, file=file)
        print()


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',
                        level=logging.DEBUG)
    main(Path.cwd() / 'projects')
