"""Package the built workbench and a bounded saved-data export in one HTML file."""

import argparse
import hashlib
import json
from pathlib import Path
import re


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--build', type=Path, default=Path('ui/dist'))
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    build = args.build.resolve()
    template = (build / 'index.html').read_text(encoding='utf-8')

    def asset(reference, closing):
        path = (build / reference.lstrip('/')).resolve()
        if not path.is_relative_to(build) or not path.is_file() or path.stat().st_size > 2_000_000:
            raise ValueError('invalid workbench build asset')
        content = path.read_text(encoding='utf-8')
        if closing in content.lower():
            raise ValueError('inline asset contains an HTML closing tag')
        return content

    template, scripts = re.subn(
        r'<script type="module" crossorigin src="([^"]+)"></script>',
        lambda match: '<script type="module">' + asset(match[1], '</script') + '</script>', template,
    )
    template, styles = re.subn(
        r'<link rel="stylesheet" crossorigin href="([^"]+)">',
        lambda match: '<style>' + asset(match[1], '</style') + '</style>', template,
    )
    if scripts != 1 or styles != 1 or template.count('</body>') != 1:
        raise ValueError('expected a single built workbench entry point')
    if args.data.stat().st_size > 16_000_000:
        raise ValueError('workbench data exceeds standalone size budget')
    raw = args.data.read_bytes()
    payload = json.loads(raw)
    if not payload.get('bundles') or not payload.get('matrix'):
        raise ValueError('workbench data requires saved model bundles and lexical context')
    identity = hashlib.sha256(raw).hexdigest()
    # JSON is inert data. Escape markup that could terminate the script element.
    inline = raw.decode('utf-8').replace('<', '\\u003c')
    template = template.replace('</body>',
        f'<script id="workbench-data" type="application/json">{inline}</script>\n</body>')
    template = template.replace('</head>', f'<meta name="diachron-data-sha256" content="{identity}">\n</head>')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(template, encoding='utf-8')
    print(json.dumps({'out': str(args.out), 'data_sha256': identity, 'models': list(payload['bundles']),
                      'html_bytes': args.out.stat().st_size}))


if __name__ == '__main__':
    main()
