"""Execute one bounded hosted research job and retain every command result."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('kind', choices=('prediction', 'sensitivity'))
    parser.add_argument('--configuration')
    args = parser.parse_args()
    protocol = json.loads(Path('data/manifests/sensitivity-protocol.json').read_text())
    names = {row['name'] for row in protocol['configurations']}
    if args.kind == 'sensitivity' and args.configuration not in names:
        raise ValueError('sensitivity jobs require one declared configuration')
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    if subprocess.check_output(['git', 'status', '--porcelain'], text=True).strip():
        raise ValueError('hosted research requires a clean source checkout')
    name = args.configuration if args.kind == 'sensitivity' else 'prediction'
    root = Path('runs') / name
    root.mkdir(parents=True, exist_ok=False)
    inputs = Path('runs/inputs')
    write(root / 'ENVIRONMENT.json', {
        'source_commit': commit, 'python': sys.version, 'platform': platform.platform(),
        'versions': {package: importlib.metadata.version(package) for package in ('numpy', 'scipy', 'pyarrow')},
        'hosted_run': os.environ.get('GITHUB_RUN_ID'), 'attempt': os.environ.get('GITHUB_RUN_ATTEMPT'),
        'job': os.environ.get('GITHUB_JOB'), 'workers': 4,
    })
    from diachron import _core
    binary = Path(_core.__file__)
    write(root / 'SOURCE_INPUTS.json', {
        'source_commit': commit, 'native_binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
        'tracked_file_sha256': {name: hashlib.sha256(Path(name).read_bytes()).hexdigest()
                                for name in subprocess.check_output(['git', 'ls-files', '-z'], text=True).split('\0')
                                if name and Path(name).is_file()},
    })
    commands = []

    def run(label, arguments):
        start = time.time()
        command = [sys.executable, '-m', 'diachron.cli.main', *arguments]
        with (root / f'{label}.stdout.json').open('w') as stdout, (root / f'{label}.stderr.log').open('w') as stderr:
            process = subprocess.Popen(command, stdout=stdout, stderr=stderr, start_new_session=True)
            try:
                code = process.wait(timeout=18000)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                code = 124
        record = {'label': label, 'arguments': command, 'exit': code,
                  'seconds': time.time() - start, 'source_commit': commit}
        commands.append(record)
        with (root / 'COMMANDS.jsonl').open('a') as stream:
            stream.write(json.dumps(record) + '\n')
        return code

    status = {'source_commit': commit, 'kind': args.kind, 'configuration': args.configuration,
              'status': 'failed', 'full_project_acceptance': False}
    try:
        if run('acquire', ['acquire', '--out', str(inputs)]):
            return 1
        source_args = ['--source', str(inputs / 'iecor'), '--catalog', str(inputs / 'glottolog53.zip'),
                       '--edh-source', str(inputs / 'edh'), '--edh-geography', str(inputs / 'edh-geography.json')]
        command = 'predict' if args.kind == 'prediction' else 'sensitivity'
        fit = root / 'fit'
        arguments = [command, '--study', 'iecor', '--out', str(fit), '--chains', '4',
                     '--warmup', '3000', '--iterations', '12000', '--jobs', '4', *source_args]
        if args.kind == 'sensitivity':
            arguments += ['--configuration', args.configuration]
        if run('fit', arguments):
            return 1
        if args.kind == 'prediction':
            target = fit
        else:
            result = json.loads((fit / 'sensitivity.json').read_text())
            target = fit / args.configuration / result['configs'][0]['model']['name']
        if run('verify', ['verify', '--evidence', str(target), '--expected-commit', commit, *source_args]):
            return 1
        status['status'] = 'executed_and_verified'
        return 0
    finally:
        status['commands'] = commands
        write(root / 'JOB.json', status)


if __name__ == '__main__':
    raise SystemExit(main())
