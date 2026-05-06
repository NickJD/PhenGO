#!/usr/bin/env python3
"""PhenGO toolkit entrypoint

Provides a single CLI to run small utility tools bundled with PhenGO.

Usage examples:
  python -m PhenGO.tools arff-validate -i data.arff -o outdir
  python -m PhenGO.tools report -i results/ -o report.html
  python -m PhenGO.tools go-map -i go_ids.txt --obo go.obo -o mapped.tsv
"""
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

import sys
import argparse
import subprocess
import logging
from ..constants import configure_logger


def run_module(module_name, argv):
    cmd = [sys.executable, '-m', module_name] + argv
    return subprocess.call(cmd)


def main():
    configure_logger('PhenGO.tools', enable_file=False)
    logger = logging.getLogger('PhenGO.tools')

    parser = argparse.ArgumentParser(prog='phengo-tools', description='PhenGO small utilities')
    sub = parser.add_subparsers(dest='command')

    p_av = sub.add_parser('arff-validate', help='Validate an ARFF file')
    p_av.add_argument('-i', '--input', required=True)
    p_av.add_argument('-o', '--output-dir', default='.')
    p_av.add_argument('--fail-on-error', action='store_true')

    p_rg = sub.add_parser('report', help='Generate HTML report from results dir')
    p_rg.add_argument('-i', '--input-dir', required=True)
    p_rg.add_argument('-o', '--output', default='report.html')

    p_gm = sub.add_parser('go-map', help='Map GO IDs using OBO file')
    p_gm.add_argument('-i', '--input', required=True)
    p_gm.add_argument('--obo', required=True)
    p_gm.add_argument('-o', '--output', default='go_mapped.tsv')

    args, rest = parser.parse_known_args()
    if not args.command:
        parser.print_help()
        sys.exit(1)

    if args.command == 'arff-validate':
        argv = ['-i', args.input, '-o', args.output_dir]
        if args.fail_on_error:
            argv.append('--fail-on-error')
        rc = run_module('PhenGO.scripts.arff_validator', argv)
        sys.exit(rc)

    if args.command == 'report':
        argv = ['-i', args.input_dir, '-o', args.output]
        rc = run_module('PhenGO.scripts.report_generator', argv)
        sys.exit(rc)

    if args.command == 'go-map':
        argv = ['-i', args.input, '--obo', args.obo, '-o', args.output]
        rc = run_module('PhenGO.scripts.go_mapper', argv)
        sys.exit(rc)


if __name__ == '__main__':
    main()

