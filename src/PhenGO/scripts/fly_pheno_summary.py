
import os as _os, sys as _sys
if __name__ == "__main__" and not __package__:
    import importlib.util as _ilu
    _here   = _os.path.dirname(_os.path.abspath(__file__))
    _src    = _os.path.normpath(_os.path.join(_here, '..', '..'))
    _phengo = _os.path.dirname(_here)
    if _src not in _sys.path:
        _sys.path.insert(0, _src)
    for _n, _f, _d in [
        ('PhenGO',         _os.path.join(_phengo, '__init__.py'), _phengo),
        ('PhenGO.scripts', _os.path.join(_here,   '__init__.py'), _here),
    ]:
        if _n not in _sys.modules:
            _sp = _ilu.spec_from_file_location(_n, _f, submodule_search_locations=[_d])
            _mo = _ilu.module_from_spec(_sp)
            _mo.__path__ = [_d]
            _sys.modules[_n] = _mo
            _sp.loader.exec_module(_mo)
    __package__ = 'PhenGO.scripts'
    del _ilu, _here, _src, _phengo, _n, _f, _d, _sp, _mo

import argparse
import csv
import os
import logging
from ..constants import configure_logger

def parse_input_file(input_file):
    """
    Reads the tab-delimited input file and returns a list of records (dicts).
    """
    records = []
    with open(input_file, 'r', encoding='utf-8') as f:
        for line in f:
            # Skip empty or comment lines
            if not line.strip() or line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue  # skip incomplete lines

            allele_full = parts[0]
            fb_id = parts[1]
            phenotype = parts[2]
            reference = parts[3]

            # Extract gene name: before first '['
            if '[' in allele_full:
                gene = allele_full.split('[')[0]
            else:
                gene = allele_full  # fallback

            # Extract context: after 'with' if present
            if 'with' in phenotype:
                context = phenotype.split('with', 1)[1].strip()
            else:
                context = ''

            # Simplify phenotype
            phen_lower = phenotype.lower()
            if 'lethal' in phen_lower or 'die' in phen_lower:
                simple_pheno = 'lethal'
            elif 'viable' in phen_lower:
                simple_pheno = 'viable'
            else:
                simple_pheno = 'other'

            record = {
                'gene': gene,
                'allele': allele_full,
                'fb_id': fb_id,
                'phenotype': phenotype,
                'context': context,
                'simple_pheno': simple_pheno,
                'reference': reference
            }
            records.append(record)

    return records

def write_raw_csv(records, output_file):
    """
    Write raw observations CSV.
    """
    fieldnames = ['gene', 'allele', 'context', 'phenotype', 'simple_pheno']
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for rec in records:
            writer.writerow({
                'gene': rec['gene'],
                'allele': rec['allele'],
                'context': rec['context'],
                'phenotype': rec['phenotype'],
                'simple_pheno': rec['simple_pheno']
            })
    logger = logging.getLogger('PhenGO.fly_pheno_summary')
    logger.info(f"✅ Raw observations saved to: {output_file}")

def summarise_by_gene(records):
    """
    Returns a list of summary dicts per gene.
    """
    summary = {}
    for rec in records:
        gene = rec['gene']
        if gene not in summary:
            summary[gene] = {'ever_lethal': False, 'ever_viable': False}

        if rec['simple_pheno'] == 'lethal':
            summary[gene]['ever_lethal'] = True
        elif rec['simple_pheno'] == 'viable':
            summary[gene]['ever_viable'] = True

    # Build list of summary rows
    summary_rows = []
    for gene, stats in summary.items():
        if stats['ever_lethal'] and stats['ever_viable']:
            text = 'sometimes lethal, sometimes viable (context-dependent)'
        elif stats['ever_lethal']:
            text = 'always lethal in data'
        elif stats['ever_viable']:
            text = 'always viable in data'
        else:
            text = 'no lethal/viable phenotypes found'

        summary_rows.append({
            'gene': gene,
            'ever_lethal': str(stats['ever_lethal']),
            'ever_viable': str(stats['ever_viable']),
            'summary': text
        })
    return summary_rows

def write_summary_csv(summary_rows, output_file):
    """
    Write the summary CSV.
    """
    fieldnames = ['gene', 'ever_lethal', 'ever_viable', 'summary']
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)
    logger = logging.getLogger('PhenGO.fly_pheno_summary')
    logger.info(f"✅ Summary table saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Summarise FlyBase allele phenotype data by gene (no pandas)")
    parser.add_argument('--input', required=True, help='Input tab-delimited file like allele_phenotypic_data_fb_2017_05')
    parser.add_argument('--raw_csv', required=True, help='Output raw observations CSV')
    parser.add_argument('--summary_csv', required=True, help='Output summary CSV per gene')
    args = parser.parse_args()

    # Configure logger (console only)
    configure_logger('PhenGO.fly_pheno_summary', enable_file=False)
    logger = logging.getLogger('PhenGO.fly_pheno_summary')

    if not os.path.exists(args.input):
        logger.error(f"❌ Input file not found: {args.input}")
        return

    logger.info(f"📦 Reading input: {args.input}")
    records = parse_input_file(args.input)

    logger.info(f"📝 Writing raw observations CSV...")
    write_raw_csv(records, args.raw_csv)

    logger.info(f"📊 Summarising per gene...")
    summary_rows = summarise_by_gene(records)

    logger.info(f"📝 Writing summary CSV...")
    write_summary_csv(summary_rows, args.summary_csv)

    logger.info("✅ Done.")

if __name__ == '__main__':
    main()
