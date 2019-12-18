import os


def count_cells_in_file_headers(paths, delimiter):
    """"""
    cell_count = 0
    for path in paths:
        with open(path, 'r') as file:
            for line in file:
                # only return the first line for each file, which should be the header
                cell_count += len([l for l in line.split(delimiter) if l.strip() and l.strip() not in ['""', 'A-6']])
                break
    return cell_count


def lines_in_files(paths, header_prefix):
    """"""
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
        elif accession_code == 'GSE107585':
            count = count_cells_in_file_headers(paths=['accessions/GSE107585/GSE107585_Mouse_kidney_single_cell_datamatrix.txt'],
                                                delimiter='\t')
        elif accession_code == 'GSE115469':
            count = count_cells_in_file_headers(paths=['accessions/GSE115469/GSE115469_Data.csv'],
                                                delimiter=',')
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
        else:
            count = 'Accession Not Found In Parser List.'
        print(f'{accession_code}: {count}')
    except FileNotFoundError as e:
        print(f'{accession_code}: {e}')


for i in ['GSE81547', 'GSE81904', 'GSE84133', 'GSE93593', 'GSE95435', 'GSE106540', 'GSE107585', 'GSE115469', 'GSE124472']:
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
    """)

