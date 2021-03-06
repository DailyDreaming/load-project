import os

"""
Although this parser can be rerun to reproduce the cell counts given,
I don't see much need to?  This should probably be only done once?

It's a decent record though if anyone finds fault with the cell number and wants to see how it was determined.
"""

already_seen = ['GSE81547',
                'GSE81904',
                'GSE84133',
                'GSE93593',
                'GSE95435',
                'GSE106540',
                'GSE107585',
                'GSE115469',
                'GSE124472',
                'GSE110154',
                'GSE44183',
                'GSE75367',
                'GSE75659',
                'GSE75688',
                'GSE81383',
                'GSE81608',
                'GSE97104',
                'GSE99795',
                'GSE110154',
                'GSE116237',
                'GSE124494',
                'GSE130606',
                'GSE118127',
                'GSE135889']

skip_headers = ['gene_id', 'gene_name', 'gene_type', '""', 'A-6', 'gene.id']


def count_cells_in_file_headers(paths, delimiter, line_number=0):
    """Some files have thousands of delimited fields on the first line for each cell."""
    cell_count = 0
    current_line_number = 0
    for path in paths:
        with open(path, 'r') as file:
            for line in file:
                if current_line_number == line_number:
                    # only return the relevant header line for the file, which may not be the first line
                    cell_count += len([l for l in line.split(delimiter) if l.strip() and l.strip() not in skip_headers])
                    break
                current_line_number += 1
    return cell_count


def lines_in_files(paths, header_prefix):
    """Counts lines in files that have one line per cell."""
    cell_count = 0
    for path in paths:
        with open(path, 'r') as file:
            for line in file:
                if line.strip() and not line.strip().startswith(header_prefix):
                    cell_count += 1
    return cell_count


def count_files(dir_path, prefix='', suffix=''):
    """Use if there is one file per cell."""
    return len([f for f in os.listdir(dir_path) if f.startswith(prefix) and f.endswith(suffix)])


def cell_counts(accession_code):
    try:
        if accession_code == 'GSE81547':
            count = count_files(dir_path="accessions/GSE81547", suffix='.csv.gz')
        elif accession_code == 'GSE81904':
            count = count_cells_in_file_headers(paths=["accessions/GSE81904/GSE81904_BipolarUMICounts_Cell2016.txt"],
                                                delimiter='\t')
        elif accession_code == 'GSE84133':
            count = lines_in_files(paths=['accessions/GSE84133/GSM2230757_human1_umifm_counts.csv',
                                          'accessions/GSE84133/GSM2230758_human2_umifm_counts.csv',
                                          'accessions/GSE84133/GSM2230759_human3_umifm_counts.csv',
                                          'accessions/GSE84133/GSM2230760_human4_umifm_counts.csv',
                                          'accessions/GSE84133/GSM2230761_mouse1_umifm_counts.csv',
                                          'accessions/GSE84133/GSM2230762_mouse2_umifm_counts.csv'],
                                   header_prefix=',barcode')
        elif accession_code == 'GSE93593':
            count = count_cells_in_file_headers(paths=['accessions/GSE93593/GSE93593_counts.csv'],
                                                delimiter=',')
        elif accession_code == 'GSE95435':
            count = count_cells_in_file_headers(paths=['accessions/GSE95435/GSE95435_P150057_full_gene_count_table.txt'],
                                                delimiter='\t')
        elif accession_code == 'GSE106540':
            count = count_cells_in_file_headers(paths=['accessions/GSE106540/GSE106540_SC_raw_counts.txt'],
                                                delimiter='\t')
            # GSE107585
        elif accession_code == 'GSE107585':
            count = count_cells_in_file_headers(paths=['accessions/GSE107585/GSE107585_Mouse_kidney_single_cell_datamatrix.txt'],
                                                delimiter='\t')
        elif accession_code == 'GSE115469':
            count = count_cells_in_file_headers(paths=['accessions/GSE115469/GSE115469_Data.csv'],
                                                delimiter=',')
            # GSE124472
        elif accession_code == 'GSE124472':
            count = lines_in_files(paths=['accessions/GSE124472/GSM3534656_H17w_Z1_raw_counts/barcodes.tsv',
                                          'accessions/GSE124472/GSM3534657_H17w_Z2_raw_counts/barcodes.tsv',
                                          'accessions/GSE124472/GSM3534658_H15w_Z1_raw_counts/barcodes.tsv',
                                          'accessions/GSE124472/GSM3534659_H15w_Z2_raw_counts/barcodes.tsv',
                                          'accessions/GSE124472/GSM3534660_HuOrg_D16_1_raw_counts/barcodes.tsv',
                                          'accessions/GSE124472/GSM3534661_HuOrg_D16_2_raw_counts/barcodes.tsv',
                                          'accessions/GSE124472/GSM3534662_HuOrg_D28_1_raw_counts/barcodes.tsv',
                                          'accessions/GSE124472/GSM3534663_HuOrg_D28_2_raw_counts/barcodes.tsv'],
                                   header_prefix='thereisnoheaderprefixforthis >.>')
        elif accession_code == 'GSE110154':
            count = count_files(dir_path="accessions/GSE110154", suffix='.csv.gz')
        elif accession_code == 'GSE44183':
            count = count_cells_in_file_headers(paths=['accessions/GSE44183/GSE44183_human_expression_mat.txt',
                                                       'accessions/GSE44183/GSE44183_mouse_expression_mat.txt'],
                                                delimiter='\t')
        elif accession_code == 'GSE75367':
            count = count_cells_in_file_headers(paths=['accessions/GSE75367/GSE75367_readCounts.txt'],
                                                delimiter='\t')
        elif accession_code == 'GSE75659':
            count = count_files(dir_path="accessions/GSE75659", suffix='.txt.gz')
        elif accession_code == 'GSE75688':
            count = count_cells_in_file_headers(paths=['accessions/GSE75688/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt'],
                                                delimiter='\t')
        elif accession_code == 'GSE81383':
            count = count_cells_in_file_headers(paths=['accessions/GSE81383/GSE81383_data_melanoma_scRNAseq_BT_Mel.txt'],
                                                delimiter='\t')
        elif accession_code == 'GSE81608':
            count = count_cells_in_file_headers(paths=['accessions/GSE81608/GSE81608_human_islets_rpkm.txt'],
                                                delimiter='\t')
        elif accession_code == 'GSE97104':
            count = count_cells_in_file_headers(paths=['accessions/GSE97104/GSE97104_all_umi.mtx.txt'],
                                                delimiter='\t',
                                                line_number=3)
        elif accession_code == 'GSE99795':
            count = count_cells_in_file_headers(paths=['accessions/GSE99795/GSE99795_raw.txt'],
                                                delimiter='\t') - 8  # header has one massive first category
        elif accession_code == 'GSE110154':
            count = count_files(dir_path="accessions/GSE110154", suffix='.csv.gz')
        elif accession_code == 'GSE116237':
            count = count_cells_in_file_headers(paths=['accessions/GSE116237/GSE116237_scRNAseq_expressionMatrix.txt'],
                                                delimiter=',')
        elif accession_code == 'GSE124494':
            count = lines_in_files(paths=['accessions/GSE124494/GSM3535276_AXLN1_barcodes.tsv',
                                          'accessions/GSE124494/GSM3535277_AXLN2_barcodes.tsv',
                                          'accessions/GSE124494/GSM3535278_AXLN3_barcodes.tsv',
                                          'accessions/GSE124494/GSM3535279_AXLN4_barcodes.tsv',
                                          'accessions/GSE124494/GSM3535280_HNLN1_barcodes.tsv',
                                          'accessions/GSE124494/GSM3535281_HNLN2_barcodes.tsv'],
                                   header_prefix='thereisnoheaderprefixforthis >.>')
        elif accession_code == 'GSE130606':
            count = lines_in_files(paths=['accessions/GSE130606/GSE130606_barcodes.tsv'],
                                   header_prefix='thereisnoheaderprefixforthis >.>')
        elif accession_code in ('GSE118127', 'GSE135889'):
            count = 'Accession Has Special Problems.'
        else:
            count = 'Accession Not Found In Parser List.'
        print(f'{accession_code}: {count}')
    except FileNotFoundError as e:
        print(f'{accession_code}: {e}')


for i in already_seen:
    cell_counts(accession_code=i)

print("""
Cell counts expected:

     GSE81547: 2544
     GSE81904: 44994
     GSE84133: 10455
     GSE93593: 1733
     GSE95435: 96
    GSE106540: 2244
    GSE107585: 43745
    GSE115469: 8444
    GSE124472: 22687
    GSE110154: 1860  # 10
     GSE44183: 48
     GSE75367: 74
     GSE75659: 1318
     GSE75688: 563
     GSE81383: 307
     GSE81608: 1600
     GSE97104: 35016
     GSE99795: 147
    GSE116237: 674
    GSE124494: 33257  # 20
    GSE130606: 7924

    GSE118127: h5 files that need parsing
    GSE135889: No Supplementary Files
    """)
